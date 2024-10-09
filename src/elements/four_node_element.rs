use crate::simulation::{Simulation};
use super::base_element::{BaseElement, Material, ElementFields, ElementType};
use nalgebra as na;
use na::{DMatrix, DVector};
use crate::utilities::{compute_von_mises};
use serde::{Serialize, Deserialize};
use log::{debug, trace};

const TOLERANCE: f64 = 1e-10;

#[derive(Serialize, Deserialize, Debug)] 
pub struct FourNodeElement {
    id: usize,
    connectivity: Vec<usize>,
    material: Material,
    #[serde(skip, default = "default_deformation_gradient")]
    deformation_gradient: DMatrix<f64>, // 2x2 for 2D elements
    #[serde(skip, default = "default_zero_matrix")]
    stiffness: DMatrix<f64>,
    #[serde(skip, default = "default_zero_matrix")]
    mass: DMatrix<f64>,
    #[serde(skip, default = "Vec::new")]//lumped mass matrix
    lumped_mass: Vec<f64>,
    active: bool,
}

fn default_deformation_gradient() -> DMatrix<f64> {
    DMatrix::<f64>::identity(2, 2)
}

fn default_zero_matrix() -> DMatrix<f64> {
    DMatrix::zeros(8, 8) // 4 nodes * 2 DOF per node
}


fn get_area_of_triangle(a: &na::Vector3<f64>, b: &na::Vector3<f64>, c: &na::Vector3<f64>) -> f64 {
    let ab = b - a;
    let ac = c - a;
    ab.cross(&ac).magnitude() / 2.0
}

impl FourNodeElement {
    pub fn new(id: usize, connectivity: Vec<usize>, material: Material) -> Self {
        assert_eq!(
            connectivity.len(),
            4,
            "4 nodes required for a four-node plane element"
        );
        FourNodeElement {
            id,
            connectivity,
            material,
            deformation_gradient: DMatrix::<f64>::identity(2, 2),
            stiffness: default_zero_matrix(),
            mass: default_zero_matrix(),
            lumped_mass: Vec::new(),
            active: false,
        }
    }

    fn center(&self, simulation: &Simulation) -> na::Vector3<f64> {
        let x_pos = self.get_x_vector(simulation);
        let center = (x_pos[0] + x_pos[1] + x_pos[2] + x_pos[3]) / 4.0;
        center  
    }

    fn get_gauss_points() -> Vec<(f64, f64, f64)> { // xi, eta, weight
        trace!("Generating Gauss points for four-node element");
        let mut gauss_points: Vec<(f64, f64, f64)> = Vec::new();
        let a = 1.0 / (3.0 as f64).sqrt();
        // Gauss points for 2x2 quadrature in a plane
        gauss_points.push((-a, -a, 1.0));
        gauss_points.push((a, -a, 1.0));
        gauss_points.push((a, a, 1.0));
        gauss_points.push((-a, a, 1.0));
        gauss_points
    }

    fn get_plane_normal(&self, simulation: &Simulation) -> na::Vector3<f64> {
        let node1 = simulation.get_node(self.connectivity[0]).unwrap();
        let node2 = simulation.get_node(self.connectivity[1]).unwrap();
        let node3 = simulation.get_node(self.connectivity[2]).unwrap();
        let vector1 = node2.position - node1.position;
        let vector2 = node3.position - node1.position;
        let normal = vector1.cross(&vector2);
        normal
    }



    fn is_point_inside(&self, point: na::Vector3<f64>, simulation: &Simulation) -> bool {
        let nodes = self.get_x_vector(simulation);
        
        // Split the quadrilateral into two triangles
        let triangle1 = [nodes[0], nodes[1], nodes[2]];
        let triangle2 = [nodes[0], nodes[2], nodes[3]];
        
        // Check if the point is inside either triangle
        self.point_in_triangle(&point, &triangle1) || self.point_in_triangle(&point, &triangle2)
    }

    fn point_in_triangle(&self, p: &na::Vector3<f64>, triangle: &[na::Vector3<f64>; 3]) -> bool {
        let v0 = triangle[2] - triangle[0];
        let v1 = triangle[1] - triangle[0];
        let v2 = p - triangle[0];

        let dot00 = v0.dot(&v0);
        let dot01 = v0.dot(&v1);
        let dot02 = v0.dot(&v2);
        let dot11 = v1.dot(&v1);
        let dot12 = v1.dot(&v2);

        let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
        let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
        let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

        (u >= -TOLERANCE) && (v >= -TOLERANCE) && (u + v <= 1.0 + TOLERANCE)
    }

    fn get_area(&self, simulation: &Simulation) -> f64 {
        let x_pos = self.get_x_vector(simulation);
        let area = get_area_of_triangle(&x_pos[0], &x_pos[1], &x_pos[2]) + get_area_of_triangle(&x_pos[0], &x_pos[2], &x_pos[3]);
        area
    }

    fn get_x_vector(&self, simulation: &Simulation) -> Vec<na::Vector3<f64>> {
        let x = self.get_x(simulation);
        let n1 = x.column(0);
        let n2 = x.column(1);
        let n3 = x.column(2);
        let n4 = x.column(3);
        vec![na::Vector3::new(n1[0], n1[1], n1[2]), na::Vector3::new(n2[0], n2[1], n2[2]), na::Vector3::new(n3[0], n3[1], n3[2]), na::Vector3::new(n4[0], n4[1], n4[2])]    
    }

    


}

#[typetag::serde]
impl BaseElement for FourNodeElement {

    fn get_signed_distance_vector(&self, point: na::Vector3<f64>, simulation: &Simulation) -> Option<na::Vector3<f64>> {
        // Get the plane normal
        let normal = self.get_plane_normal(simulation).normalize();
        let nodes = self.get_x_vector(simulation);
        let center = self.center(simulation);
        let signed_distance = normal.dot(&(center - point));

        if signed_distance > 0.0 {
            return None;
        }

        let projected_point = point - normal * signed_distance;
        
        // Check if the projected point is within the element
        if self.is_point_inside(projected_point, simulation) {
            Some(signed_distance * normal)
        } else {
            None
        }
    }


    fn get_id(&self) -> usize {
        self.id
    }

    fn get_connectivity(&self) -> &Vec<usize> {
        &self.connectivity
    }

    fn get_material(&self) -> &Material {
        &self.material
    }

    fn get_deformation_gradient(&self) -> &DMatrix<f64> {
        &self.deformation_gradient
    }

    fn get_shape_functions(&self, xi: f64, eta: f64, _zeta: f64) -> DVector<f64> {
        // Shape functions for 2D bilinear quadrilateral element
        let mut shape_functions: DVector<f64> = DVector::<f64>::zeros(4);
        shape_functions[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
        shape_functions[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
        shape_functions[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
        shape_functions[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
        shape_functions
    }

    fn get_x(&self, simulation: &Simulation) -> DMatrix<f64> {
        // Compute X matrix in 3D coordinates for each node
        let mut X = DMatrix::<f64>::zeros(3, 4);
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            X[(0, i)] = node.position[0];
            X[(1, i)] = node.position[1];
            X[(2, i)] = node.position[2];
        }
        X
    }

    fn get_u(&self, simulation: &Simulation) -> DVector<f64> {
        // Compute displacement vector u in 3D (though we're using 2D calculations)
        let mut u = DVector::<f64>::zeros(12); // 4 nodes * 3 DOF per node (in 3D)
        for (i, node_id) in self.connectivity.iter().enumerate() {
            let node = simulation.get_node(*node_id).unwrap();
            u[3 * i] = node.displacement[0];
            u[3 * i + 1] = node.displacement[1];
            u[3 * i + 2] = node.displacement[2];
        }
        u
    }

    fn get_shape_derivatives(&self, xi: f64, eta: f64, _zeta: f64) -> DMatrix<f64> {
        // Return 2x4 matrix of shape function derivatives (for plane elements)
        let matrix_data = vec![
            [-0.25 * (1.0 - eta), -0.25 * (1.0 - xi)],
            [0.25 * (1.0 - eta), -0.25 * (1.0 + xi)],
            [0.25 * (1.0 + eta), 0.25 * (1.0 + xi)],
            [-0.25 * (1.0 + eta), 0.25 * (1.0 - xi)],
        ];
        DMatrix::from_row_slice(4, 2, &matrix_data.concat())
    }

    fn compute_jacobian_matrix(&self, xi: f64, eta: f64, _zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        todo!()
    }

    fn compute_b(&self, xi: f64, eta: f64, _zeta: f64, simulation: &Simulation) -> DMatrix<f64> {
        todo!()
    }
    fn compute_stiffness(&self, simulation: &Simulation) -> DMatrix<f64> {
        todo!()
    }
    fn get_stiffness(&self) -> &DMatrix<f64> {
        &self.stiffness
    }
    fn compute_mass(&self, simulation: &Simulation) -> DMatrix<f64> {
        todo!()
    }   
    fn get_mass(&self) -> &DMatrix<f64> {
        &self.mass
    }
    fn set_lumped_mass(&mut self, lumped_mass: &DMatrix<f64>) -> f64 {
        todo!()
    }   
    fn compute_stress(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DVector<f64> {
        todo!()
    }
    fn compute_strain(&self, xi: f64, eta: f64, zeta: f64, simulation: &Simulation) -> DVector<f64> {
        todo!()
    }

    fn get_lumped_mass(&self) -> &Vec<f64> {
        &self.lumped_mass
    }

    fn type_name(&self) -> ElementType {
        ElementType::Quad
    }           
}