//! Section cut rendering - compute and render cross-section polygons where the 
//! clipping plane intersects elements, with field value interpolation.

use crate::state::{MeshState, SimulationResults, ColorMode, ClipAxis};
use crate::render_cache::ViewTransform;
use eframe::egui;

/// A point on the section cut surface with interpolated field value
#[derive(Clone, Debug)]
pub struct CutPoint {
    /// 3D position in world space
    pub position: [f32; 3],
    /// Parametric coordinate within the element (for interpolation)
    pub xi: f32,
    pub eta: f32,
    pub zeta: f32,
    /// Interpolated field value at this point
    pub field_value: f64,
}

/// Element type enumeration for shape function selection
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum CutElementType {
    /// 8-node linear hexahedral (C3D8)
    Hex8,
    /// 20-node quadratic hexahedral (C3D20)  
    Hex20,
    /// 4-node linear tetrahedral (C3D4)
    Tet4,
    /// Unknown/unsupported
    Unknown,
}

impl CutElementType {
    pub fn from_str(s: &str) -> Self {
        match s {
            "C3D8" => CutElementType::Hex8,
            "C3D20" => CutElementType::Hex20,
            "C3D4" => CutElementType::Tet4,
            _ => CutElementType::Unknown,
        }
    }
    
    pub fn num_nodes(&self) -> usize {
        match self {
            CutElementType::Hex8 => 8,
            CutElementType::Hex20 => 20,
            CutElementType::Tet4 => 4,
            CutElementType::Unknown => 0,
        }
    }
}

/// A polygon formed by the intersection of a plane with an element
#[derive(Clone, Debug)]
pub struct CutPolygon {
    /// Vertices of the polygon (in order, forming a convex polygon)
    pub vertices: Vec<CutPoint>,
    /// Element ID this polygon came from
    pub element_id: usize,
    /// Center of the polygon (for depth sorting)
    pub center: [f32; 3],
    /// Normal vector of the polygon (same as clip plane normal)
    pub normal: [f32; 3],
    /// Element type for shape function selection
    pub element_type: CutElementType,
    /// Nodal field values for shape function interpolation (all nodes of the element)
    pub nodal_field_values: Vec<f64>,
    /// Node positions in world space (for inverse mapping if needed)
    pub node_positions: Vec<[f32; 3]>,
}

/// Edge of a hex element (node indices)
const HEX_EDGES: [[usize; 2]; 12] = [
    // Bottom face edges
    [0, 1], [1, 2], [2, 3], [3, 0],
    // Top face edges
    [4, 5], [5, 6], [6, 7], [7, 4],
    // Vertical edges connecting top and bottom
    [0, 4], [1, 5], [2, 6], [3, 7],
];

/// Natural coordinates of hex element corners (standard isoparametric ordering)
const HEX_NATURAL_COORDS: [[f32; 3]; 8] = [
    [-1.0, -1.0, -1.0], // 0
    [ 1.0, -1.0, -1.0], // 1
    [ 1.0,  1.0, -1.0], // 2
    [-1.0,  1.0, -1.0], // 3
    [-1.0, -1.0,  1.0], // 4
    [ 1.0, -1.0,  1.0], // 5
    [ 1.0,  1.0,  1.0], // 6
    [-1.0,  1.0,  1.0], // 7
];

/// Pre-compute nodal field values for element-based fields
/// Uses actual nodal values from solver when available, falls back to 
/// element-averaged values for backwards compatibility.
fn precompute_nodal_field_values(
    _mesh_state: &MeshState,
    results: &Option<SimulationResults>,
    color_mode: ColorMode,
    stress_component: usize,
    strain_component: usize,
) -> std::collections::HashMap<usize, f64> {
    use std::collections::HashMap;
    
    let mut nodal_values: HashMap<usize, f64> = HashMap::new();
    
    // For displacement mode, we don't need pre-computation (already nodal)
    // For solid mode, we don't need field values
    let results = match results {
        Some(r) => r,
        None => return nodal_values,
    };
    
    match color_mode {
        ColorMode::Solid | ColorMode::Displacement => {
            // No pre-computation needed
            return nodal_values;
        }
        ColorMode::VonMises => {
            // Use actual nodal von Mises values from solver
            for (node_id, &val) in results.nodal_von_mises.iter().enumerate() {
                nodal_values.insert(node_id, val);
            }
        }
        ColorMode::Stress => {
            // Use actual nodal stress values from solver
            for (node_id, stress) in results.nodal_stress.iter().enumerate() {
                if stress_component < 6 {
                    nodal_values.insert(node_id, stress[stress_component]);
                }
            }
        }
        ColorMode::Strain => {
            // Use actual nodal strain values from solver
            for (node_id, strain) in results.nodal_strain.iter().enumerate() {
                if strain_component < 6 {
                    nodal_values.insert(node_id, strain[strain_component]);
                }
            }
        }
    }
    
    nodal_values
}

/// Compute section cut polygons for all elements intersected by the clipping plane
pub fn compute_section_cuts(
    mesh_state: &MeshState,
    results: &Option<SimulationResults>,
    _clip_axis: ClipAxis,
    clip_normal: [f32; 3],
    clip_position: f32, // Normalized -1 to 1
    displacement_scale: f32,
    color_mode: ColorMode,
    stress_component: usize,
    strain_component: usize,
) -> Vec<CutPolygon> {
    let mut cut_polygons = Vec::new();
    
    let mesh = &mesh_state.mesh;
    let bounds = &mesh_state.bounds;
    
    // Pre-compute averaged nodal values for element-based fields (Von Mises, Stress, Strain)
    // This enables smooth interpolation across element boundaries
    let nodal_field_map = precompute_nodal_field_values(
        mesh_state, results, color_mode, stress_component, strain_component
    );
    
    // Calculate actual clip plane position in world coordinates
    let center = bounds.center();
    let half_diag = bounds.diagonal() * 0.5;
    let clip_offset = clip_position * half_diag;
    
    // Plane equation: normal · (point - plane_point) = 0
    // Or: normal · point = d, where d = normal · plane_point
    let plane_d = clip_normal[0] * center[0] + 
                  clip_normal[1] * center[1] + 
                  clip_normal[2] * center[2] + clip_offset;
    
    // Process each element
    for (&el_id, element) in &mesh.elements {
        let conn = &element.connectivity;
        let el_type = CutElementType::from_str(&element.el_type);
        
        // Currently only support hex elements (8 or 20 node)
        if conn.len() < 8 {
            continue;
        }
        
        // For now, use first 8 nodes for plane intersection (corner nodes)
        // but store ALL nodal values for proper shape function interpolation
        let num_nodes = el_type.num_nodes().max(conn.len());
        
        // Get ALL node positions and field values for the element
        let mut all_node_positions: Vec<[f32; 3]> = Vec::with_capacity(num_nodes);
        let mut all_node_field_values: Vec<f64> = Vec::with_capacity(num_nodes);
        
        for (_local_idx, &node_id) in conn.iter().enumerate() {
            if let Some(node) = mesh.nodes.get(&node_id) {
                let mut pos = [
                    node.coordinates[0] as f32,
                    node.coordinates[1] as f32,
                    node.coordinates[2] as f32,
                ];
                
                // Apply displacement if available
                let mut disp_mag = 0.0;
                if let Some(res) = results {
                    if node_id * 3 + 2 < res.displacements.len() {
                        let dx = res.displacements[node_id * 3] as f32;
                        let dy = res.displacements[node_id * 3 + 1] as f32;
                        let dz = res.displacements[node_id * 3 + 2] as f32;
                        pos[0] += dx * displacement_scale;
                        pos[1] += dy * displacement_scale;
                        pos[2] += dz * displacement_scale;
                        disp_mag = ((dx * dx + dy * dy + dz * dz) as f64).sqrt();
                    }
                }
                all_node_positions.push(pos);
                
                // Get field value for this node - use pre-computed averaged values
                // for element-based fields (Von Mises, Stress, Strain)
                let field_val = match color_mode {
                    ColorMode::Solid => 0.5,
                    ColorMode::Displacement => disp_mag,
                    ColorMode::VonMises | ColorMode::Stress | ColorMode::Strain => {
                        // Use pre-computed averaged nodal value
                        nodal_field_map.get(&node_id).copied().unwrap_or(0.0)
                    }
                };
                all_node_field_values.push(field_val);
            }
        }
        
        if all_node_positions.len() < 8 {
            continue;
        }
        
        // Use corner nodes (first 8) for plane intersection test
        let corner_positions: Vec<[f32; 3]> = all_node_positions.iter().take(8).cloned().collect();
        let corner_field_values: Vec<f64> = all_node_field_values.iter().take(8).cloned().collect();
        
        // Compute signed distances from each corner node to the plane
        let mut distances: [f32; 8] = [0.0; 8];
        for (i, pos) in corner_positions.iter().enumerate() {
            distances[i] = pos[0] * clip_normal[0] + 
                          pos[1] * clip_normal[1] + 
                          pos[2] * clip_normal[2] - plane_d;
        }
        
        // Check if element is intersected by the plane
        let mut has_positive = false;
        let mut has_negative = false;
        for d in distances.iter() {
            if *d > 1e-6 { has_positive = true; }
            if *d < -1e-6 { has_negative = true; }
        }
        
        if !has_positive || !has_negative {
            continue; // Element doesn't cross the plane
        }
        
        // Find intersection points on each edge
        let mut cut_points: Vec<CutPoint> = Vec::new();
        
        for edge in HEX_EDGES.iter() {
            let i0 = edge[0];
            let i1 = edge[1];
            let d0 = distances[i0];
            let d1 = distances[i1];
            
            // Check if this edge crosses the plane
            if (d0 > 1e-6 && d1 < -1e-6) || (d0 < -1e-6 && d1 > 1e-6) {
                // Compute intersection parameter t along edge
                let t = d0 / (d0 - d1);
                
                // Interpolate position using corner nodes
                let pos = [
                    corner_positions[i0][0] + t * (corner_positions[i1][0] - corner_positions[i0][0]),
                    corner_positions[i0][1] + t * (corner_positions[i1][1] - corner_positions[i0][1]),
                    corner_positions[i0][2] + t * (corner_positions[i1][2] - corner_positions[i0][2]),
                ];
                
                // Interpolate natural coordinates
                let xi = HEX_NATURAL_COORDS[i0][0] + t * (HEX_NATURAL_COORDS[i1][0] - HEX_NATURAL_COORDS[i0][0]);
                let eta = HEX_NATURAL_COORDS[i0][1] + t * (HEX_NATURAL_COORDS[i1][1] - HEX_NATURAL_COORDS[i0][1]);
                let zeta = HEX_NATURAL_COORDS[i0][2] + t * (HEX_NATURAL_COORDS[i1][2] - HEX_NATURAL_COORDS[i0][2]);
                
                // For the cut point field value, use shape function interpolation
                // For now, use linear interpolation along edge (will be overridden in rendering)
                let field_value = corner_field_values[i0] + (t as f64) * (corner_field_values[i1] - corner_field_values[i0]);
                
                cut_points.push(CutPoint {
                    position: pos,
                    xi,
                    eta,
                    zeta,
                    field_value,
                });
            }
        }
        
        if cut_points.len() < 3 {
            continue; // Need at least 3 points to form a polygon
        }
        
        // Sort points to form a proper convex polygon
        // Use angle from centroid method
        let centroid = compute_centroid(&cut_points);
        sort_points_by_angle(&mut cut_points, &centroid, &clip_normal);
        
        // Compute polygon center for depth sorting
        let polygon_center = [
            cut_points.iter().map(|p| p.position[0]).sum::<f32>() / cut_points.len() as f32,
            cut_points.iter().map(|p| p.position[1]).sum::<f32>() / cut_points.len() as f32,
            cut_points.iter().map(|p| p.position[2]).sum::<f32>() / cut_points.len() as f32,
        ];
        
        cut_polygons.push(CutPolygon {
            vertices: cut_points,
            element_id: el_id,
            center: polygon_center,
            normal: clip_normal,
            element_type: el_type,
            nodal_field_values: all_node_field_values,
            node_positions: all_node_positions,
        });
    }
    
    cut_polygons
}

/// Get field value for a node/element based on color mode
fn get_field_value(
    results: &Option<SimulationResults>,
    element_id: usize,
    _node_id: usize,
    _local_idx: usize,
    disp_mag: f64,
    color_mode: ColorMode,
    stress_component: usize,
    strain_component: usize,
) -> f64 {
    match results {
        None => 0.0,
        Some(res) => {
            match color_mode {
                ColorMode::Solid => 0.5, // Neutral value
                ColorMode::Displacement => disp_mag,
                ColorMode::VonMises => {
                    res.von_mises.get(&element_id).copied().unwrap_or(0.0)
                }
                ColorMode::Stress => {
                    res.stresses.get(&element_id)
                        .and_then(|s| s.get(stress_component))
                        .copied()
                        .unwrap_or(0.0)
                }
                ColorMode::Strain => {
                    res.strains.get(&element_id)
                        .and_then(|s| s.get(strain_component))
                        .copied()
                        .unwrap_or(0.0)
                }
            }
        }
    }
}

/// Compute centroid of cut points
fn compute_centroid(points: &[CutPoint]) -> [f32; 3] {
    let n = points.len() as f32;
    [
        points.iter().map(|p| p.position[0]).sum::<f32>() / n,
        points.iter().map(|p| p.position[1]).sum::<f32>() / n,
        points.iter().map(|p| p.position[2]).sum::<f32>() / n,
    ]
}

/// Sort points by angle around centroid to form a proper polygon
fn sort_points_by_angle(points: &mut [CutPoint], centroid: &[f32; 3], normal: &[f32; 3]) {
    if points.len() < 3 {
        return;
    }
    
    // Create a local coordinate system on the plane
    // Find a vector not parallel to normal for cross product
    let ref_vec = if normal[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    
    // u = normal x ref_vec (perpendicular to normal)
    let mut u = [
        normal[1] * ref_vec[2] - normal[2] * ref_vec[1],
        normal[2] * ref_vec[0] - normal[0] * ref_vec[2],
        normal[0] * ref_vec[1] - normal[1] * ref_vec[0],
    ];
    
    // Normalize u
    let u_len = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    if u_len < 1e-6 {
        return; // Can't create coordinate system
    }
    u[0] /= u_len;
    u[1] /= u_len;
    u[2] /= u_len;
    
    // v = normal x u (perpendicular to both)
    let v = [
        normal[1] * u[2] - normal[2] * u[1],
        normal[2] * u[0] - normal[0] * u[2],
        normal[0] * u[1] - normal[1] * u[0],
    ];
    
    // Compute angle for each point
    let mut angles: Vec<(usize, f32)> = points.iter().enumerate().map(|(i, p)| {
        let dx = p.position[0] - centroid[0];
        let dy = p.position[1] - centroid[1];
        let dz = p.position[2] - centroid[2];
        
        // Project onto plane coordinates
        let px = dx * u[0] + dy * u[1] + dz * u[2];
        let py = dx * v[0] + dy * v[1] + dz * v[2];
        
        let angle = py.atan2(px);
        (i, angle)
    }).collect();
    
    // Sort by angle
    angles.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    
    // Reorder points
    let old_points: Vec<_> = points.to_vec();
    for (new_idx, (old_idx, _)) in angles.iter().enumerate() {
        points[new_idx] = old_points[*old_idx].clone();
    }
}

/// A projected section cut polygon ready for depth-sorted rendering
pub struct ProjectedCutPolygon {
    pub vertices: Vec<egui::Pos2>,
    pub depth: f32,
    pub color: egui::Color32,
}

/// A triangle with per-vertex colors for smooth shading
pub struct ProjectedTriangle {
    pub v0: egui::Pos2,
    pub v1: egui::Pos2,
    pub v2: egui::Pos2,
    pub c0: egui::Color32,
    pub c1: egui::Color32,
    pub c2: egui::Color32,
    pub depth: f32,
}

/// Generate projected section cut polygons with depth values for integration
/// with the main face rendering pipeline.
/// Uses fan triangulation from centroid with subdivision for smooth gradients
/// without missing edge coverage.
pub fn generate_section_cut_shapes(
    rect: egui::Rect,
    view_transform: &ViewTransform,
    cut_polygons: &[CutPolygon],
    results: &Option<SimulationResults>,
    color_mode: ColorMode,
    camera_z_fn: impl Fn([f32; 3]) -> f32,
) -> Vec<ProjectedCutPolygon> {
    if cut_polygons.is_empty() {
        return Vec::new();
    }
    
    let rect_center = rect.center();
    let half_width = rect.width() * 0.5;
    let half_height = rect.height() * 0.5;
    let aspect = rect.width() / rect.height();
    
    // Get field value range for normalization
    let (min_val, max_val) = compute_field_range(cut_polygons, results, color_mode);
    let val_range = (max_val - min_val).max(1e-10);
    
    let mut result = Vec::new();
    
    // Subdivision level for smooth gradients within each fan triangle
    // Higher = smoother but more triangles (6 gives very smooth gradients)
    let subdiv_level = 6;
    
    for polygon in cut_polygons {
        let n_verts = polygon.vertices.len();
        if n_verts < 3 {
            continue;
        }
        
        // Compute depth from polygon center
        let depth = camera_z_fn(polygon.center);
        
        // Project centroid
        let centroid_screen = match view_transform.project_close(
            polygon.center,
            rect_center,
            half_width,
            half_height,
            aspect,
        ) {
            Some(p) if p.x.is_finite() && p.y.is_finite() => p,
            _ => continue,
        };
        
        // Compute centroid parametric coords (average of vertex coords)
        let centroid_xi: f32 = polygon.vertices.iter().map(|v| v.xi).sum::<f32>() / n_verts as f32;
        let centroid_eta: f32 = polygon.vertices.iter().map(|v| v.eta).sum::<f32>() / n_verts as f32;
        let centroid_zeta: f32 = polygon.vertices.iter().map(|v| v.zeta).sum::<f32>() / n_verts as f32;
        
        // Compute centroid color
        let centroid_field = interpolate_field_value(polygon, centroid_xi, centroid_eta, centroid_zeta);
        let centroid_t = ((centroid_field - min_val) / val_range).clamp(0.0, 1.0) as f32;
        let _centroid_color = value_to_color(centroid_t);
        
        // Project all vertices
        let mut vertex_data: Vec<(egui::Pos2, f32, f32, f32, egui::Color32)> = Vec::with_capacity(n_verts);
        let mut all_valid = true;
        
        for v in &polygon.vertices {
            if let Some(screen_pos) = view_transform.project_close(
                v.position,
                rect_center,
                half_width,
                half_height,
                aspect,
            ) {
                if screen_pos.x.is_finite() && screen_pos.y.is_finite() 
                   && rect.expand(200.0).contains(screen_pos) {
                    let field_val = interpolate_field_value(polygon, v.xi, v.eta, v.zeta);
                    let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
                    let color = value_to_color(t);
                    vertex_data.push((screen_pos, v.xi, v.eta, v.zeta, color));
                } else {
                    all_valid = false;
                    break;
                }
            } else {
                all_valid = false;
                break;
            }
        }
        
        if !all_valid || vertex_data.len() < 3 {
            continue;
        }
        
        // Fan triangulation from centroid to each edge, with subdivision
        for i in 0..vertex_data.len() {
            let j = (i + 1) % vertex_data.len();
            
            let (p0, xi0, eta0, zeta0, _c0) = vertex_data[i];
            let (p1, xi1, eta1, zeta1, _c1) = vertex_data[j];
            
            // Subdivide this fan triangle into smaller triangles
            // Using barycentric coordinates (t_c, t_0, t_1) where t_c + t_0 + t_1 = 1
            
            // Helper to get screen position and color for a barycentric point
            let get_point = |tc: f32, t0: f32, t1: f32| -> (egui::Pos2, egui::Color32) {
                let px = centroid_screen.x * tc + p0.x * t0 + p1.x * t1;
                let py = centroid_screen.y * tc + p0.y * t0 + p1.y * t1;
                
                // Interpolate parametric coords
                let xi = centroid_xi * tc + xi0 * t0 + xi1 * t1;
                let eta = centroid_eta * tc + eta0 * t0 + eta1 * t1;
                let zeta = centroid_zeta * tc + zeta0 * t0 + zeta1 * t1;
                
                // Compute color via shape function interpolation
                let field_val = interpolate_field_value(polygon, xi, eta, zeta);
                let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
                let color = value_to_color(t);
                
                (egui::pos2(px, py), color)
            };
            
            let step = 1.0 / subdiv_level as f32;
            
            // Create subdivision triangles within the fan triangle
            for si in 0..subdiv_level {
                for sj in 0..(subdiv_level - si) {
                    // Barycentric coords for corners
                    let t0_base = si as f32 * step;
                    let t1_base = sj as f32 * step;
                    
                    // First triangle
                    let (pa, ca) = get_point(1.0 - t0_base - t1_base, t0_base, t1_base);
                    let (pb, cb) = get_point(1.0 - (t0_base + step) - t1_base, t0_base + step, t1_base);
                    let (pc, cc) = get_point(1.0 - t0_base - (t1_base + step), t0_base, t1_base + step);
                    
                    // Slightly expand triangle from centroid to eliminate sub-pixel gaps
                    let expand_tri = |p0: egui::Pos2, p1: egui::Pos2, p2: egui::Pos2| -> Vec<egui::Pos2> {
                        let cx = (p0.x + p1.x + p2.x) / 3.0;
                        let cy = (p0.y + p1.y + p2.y) / 3.0;
                        let expand = 0.02;
                        vec![
                            egui::pos2(p0.x + (p0.x - cx) * expand, p0.y + (p0.y - cy) * expand),
                            egui::pos2(p1.x + (p1.x - cx) * expand, p1.y + (p1.y - cy) * expand),
                            egui::pos2(p2.x + (p2.x - cx) * expand, p2.y + (p2.y - cy) * expand),
                        ]
                    };
                    
                    let avg_color = average_color3(ca, cb, cc);
                    result.push(ProjectedCutPolygon {
                        vertices: expand_tri(pa, pb, pc),
                        depth,
                        color: avg_color,
                    });
                    
                    // Second triangle (if not on the hypotenuse edge)
                    if si + sj + 1 < subdiv_level {
                        let (pd, cd) = get_point(
                            1.0 - (t0_base + step) - (t1_base + step), 
                            t0_base + step, 
                            t1_base + step
                        );
                        let avg_color2 = average_color3(cb, cd, cc);
                        result.push(ProjectedCutPolygon {
                            vertices: expand_tri(pb, pd, pc),
                            depth,
                            color: avg_color2,
                        });
                    }
                }
            }
        }
    }
    
    result
}

/// Average three colors
fn average_color3(c0: egui::Color32, c1: egui::Color32, c2: egui::Color32) -> egui::Color32 {
    let r = ((c0.r() as u32 + c1.r() as u32 + c2.r() as u32) / 3) as u8;
    let g = ((c0.g() as u32 + c1.g() as u32 + c2.g() as u32) / 3) as u8;
    let b = ((c0.b() as u32 + c1.b() as u32 + c2.b() as u32) / 3) as u8;
    let a = ((c0.a() as u32 + c1.a() as u32 + c2.a() as u32) / 3) as u8;
    egui::Color32::from_rgba_unmultiplied(r, g, b, a)
}

/// Compute signed area of 2D polygon (shoelace formula)
fn compute_2d_polygon_area(vertices: &[egui::Pos2]) -> f32 {
    let mut area = 0.0f32;
    let n = vertices.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += vertices[i].x * vertices[j].y;
        area -= vertices[j].x * vertices[i].y;
    }
    area * 0.5
}

/// Compute field value range across all cut polygons
fn compute_field_range(
    cut_polygons: &[CutPolygon],
    results: &Option<SimulationResults>,
    color_mode: ColorMode,
) -> (f64, f64) {
    // For solid color, return neutral range
    if color_mode == ColorMode::Solid {
        return (0.0, 1.0);
    }
    
    // Get range from results for global consistency
    if let Some(res) = results {
        match color_mode {
            ColorMode::Solid => (0.0, 1.0),
            ColorMode::Displacement => (0.0, res.stats.max_displacement.max(1e-10)),
            ColorMode::VonMises => (0.0, res.stats.max_von_mises.max(1e-10)),
            ColorMode::Stress => {
                // Use polygon local range or compute from stresses
                let max_val = res.stresses.values()
                    .filter_map(|s| s.first())
                    .map(|v| v.abs())
                    .fold(1e-10f64, |a: f64, b| a.max(b));
                (-max_val, max_val)
            }
            ColorMode::Strain => {
                let max_val = res.strains.values()
                    .filter_map(|s| s.first())
                    .map(|v| v.abs())
                    .fold(1e-10f64, |a: f64, b| a.max(b));
                (-max_val, max_val)
            }
        }
    } else {
        // Compute from polygon field values
        let mut min_val = f64::MAX;
        let mut max_val = f64::MIN;
        for poly in cut_polygons {
            for v in &poly.vertices {
                min_val = min_val.min(v.field_value);
                max_val = max_val.max(v.field_value);
            }
        }
        if min_val > max_val {
            (0.0, 1.0)
        } else {
            (min_val, max_val)
        }
    }
}

/// Convert normalized value (0-1) to color using jet colormap
fn value_to_color(t: f32) -> egui::Color32 {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    
    egui::Color32::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8
    )
}

/// Cache for section cut data to avoid recomputing every frame
#[derive(Default)]
pub struct SectionCutCache {
    /// Cached cut polygons
    pub polygons: Vec<CutPolygon>,
    /// Hash of parameters used to compute the cache
    pub param_hash: u64,
    /// Whether the cache is valid
    pub valid: bool,
}

impl SectionCutCache {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Check if cache is valid for current parameters
    pub fn is_valid(&self, new_hash: u64) -> bool {
        self.valid && self.param_hash == new_hash
    }
    
    /// Compute parameter hash
    pub fn compute_hash(
        clip_position: f32,
        clip_axis: ClipAxis,
        flip: bool,
        displacement_scale: f32,
        color_mode: ColorMode,
        stress_component: usize,
        strain_component: usize,
        mesh_version: u64,
        results_version: u64,
    ) -> u64 {
        let axis_val = match clip_axis {
            ClipAxis::X => 0u64,
            ClipAxis::Y => 1u64,
            ClipAxis::Z => 2u64,
            ClipAxis::Custom => 3u64,
        };
        let color_val = match color_mode {
            ColorMode::Solid => 0u64,
            ColorMode::Displacement => 1u64,
            ColorMode::VonMises => 2u64,
            ColorMode::Stress => 3u64,
            ColorMode::Strain => 4u64,
        };
        
        let mut hash = 0u64;
        hash ^= (clip_position.to_bits() as u64).wrapping_mul(31);
        hash ^= axis_val.wrapping_mul(37);
        hash ^= (flip as u64).wrapping_mul(41);
        hash ^= (displacement_scale.to_bits() as u64).wrapping_mul(43);
        hash ^= color_val.wrapping_mul(47);
        hash ^= (stress_component as u64).wrapping_mul(53);
        hash ^= (strain_component as u64).wrapping_mul(59);
        hash ^= mesh_version.wrapping_mul(61);
        hash ^= results_version.wrapping_mul(67);
        hash
    }
    
    /// Invalidate the cache
    pub fn invalidate(&mut self) {
        self.valid = false;
        self.polygons.clear();
    }
    
    /// Invalidate just the hash to force recomputation on next check
    pub fn invalidate_hash(&mut self) {
        self.param_hash = u64::MAX;  // Set to impossible value
    }
    
    /// Update cache with new polygons
    pub fn update(&mut self, polygons: Vec<CutPolygon>, param_hash: u64) {
        self.polygons = polygons;
        self.param_hash = param_hash;
        self.valid = true;
    }
}


// ============================================================================
// SHAPE FUNCTION INTERPOLATION
// ============================================================================

/// Natural coordinates for 20-node serendipity element
pub const HEX20_NATURAL_COORDS: [[f32; 3]; 20] = [
    // Corner nodes (0-7)
    [-1.0, -1.0, -1.0], // 0
    [ 1.0, -1.0, -1.0], // 1
    [ 1.0,  1.0, -1.0], // 2
    [-1.0,  1.0, -1.0], // 3
    [-1.0, -1.0,  1.0], // 4
    [ 1.0, -1.0,  1.0], // 5
    [ 1.0,  1.0,  1.0], // 6
    [-1.0,  1.0,  1.0], // 7
    // Mid-edge nodes on bottom face (z=-1)
    [ 0.0, -1.0, -1.0], // 8
    [ 1.0,  0.0, -1.0], // 9
    [ 0.0,  1.0, -1.0], // 10
    [-1.0,  0.0, -1.0], // 11
    // Mid-edge nodes on top face (z=+1)
    [ 0.0, -1.0,  1.0], // 12
    [ 1.0,  0.0,  1.0], // 13
    [ 0.0,  1.0,  1.0], // 14
    [-1.0,  0.0,  1.0], // 15
    // Mid-edge nodes on vertical edges
    [-1.0, -1.0,  0.0], // 16
    [ 1.0, -1.0,  0.0], // 17
    [ 1.0,  1.0,  0.0], // 18
    [-1.0,  1.0,  0.0], // 19
];

/// Evaluate 8-node linear hex shape functions at (xi, eta, zeta)
pub fn hex8_shape_functions(xi: f32, eta: f32, zeta: f32) -> [f64; 8] {
    let xi = xi as f64;
    let eta = eta as f64;
    let zeta = zeta as f64;
    
    let xi_m = 1.0 - xi;
    let xi_p = 1.0 + xi;
    let eta_m = 1.0 - eta;
    let eta_p = 1.0 + eta;
    let zeta_m = 1.0 - zeta;
    let zeta_p = 1.0 + zeta;
    
    [
        0.125 * xi_m * eta_m * zeta_m,
        0.125 * xi_p * eta_m * zeta_m,
        0.125 * xi_p * eta_p * zeta_m,
        0.125 * xi_m * eta_p * zeta_m,
        0.125 * xi_m * eta_m * zeta_p,
        0.125 * xi_p * eta_m * zeta_p,
        0.125 * xi_p * eta_p * zeta_p,
        0.125 * xi_m * eta_p * zeta_p,
    ]
}

/// Evaluate 20-node quadratic serendipity hex shape functions at (xi, eta, zeta)
pub fn hex20_shape_functions(xi: f32, eta: f32, zeta: f32) -> [f64; 20] {
    let xi = xi as f64;
    let eta = eta as f64;
    let zeta = zeta as f64;
    
    let xi2 = xi * xi;
    let eta2 = eta * eta;
    let zeta2 = zeta * zeta;
    
    let mut n = [0.0f64; 20];
    
    // Corner nodes (0-7) - serendipity formulation
    let corners = [
        (-1.0, -1.0, -1.0),
        ( 1.0, -1.0, -1.0),
        ( 1.0,  1.0, -1.0),
        (-1.0,  1.0, -1.0),
        (-1.0, -1.0,  1.0),
        ( 1.0, -1.0,  1.0),
        ( 1.0,  1.0,  1.0),
        (-1.0,  1.0,  1.0),
    ];
    
    for i in 0..8 {
        let (xi_i, eta_i, zeta_i) = corners[i];
        let xi_term = 1.0 + xi_i * xi;
        let eta_term = 1.0 + eta_i * eta;
        let zeta_term = 1.0 + zeta_i * zeta;
        n[i] = 0.125 * xi_term * eta_term * zeta_term * (xi_i * xi + eta_i * eta + zeta_i * zeta - 2.0);
    }
    
    // Mid-edge nodes on bottom face (z = -1)
    n[8] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 - zeta);
    n[9] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 - zeta);
    n[10] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 - zeta);
    n[11] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 - zeta);
    
    // Mid-edge nodes on top face (z = +1)
    n[12] = 0.25 * (1.0 - xi2) * (1.0 - eta) * (1.0 + zeta);
    n[13] = 0.25 * (1.0 + xi) * (1.0 - eta2) * (1.0 + zeta);
    n[14] = 0.25 * (1.0 - xi2) * (1.0 + eta) * (1.0 + zeta);
    n[15] = 0.25 * (1.0 - xi) * (1.0 - eta2) * (1.0 + zeta);
    
    // Mid-edge nodes on vertical edges
    n[16] = 0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 - zeta2);
    n[17] = 0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - zeta2);
    n[18] = 0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - zeta2);
    n[19] = 0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 - zeta2);
    
    n
}

/// Interpolate field value using shape functions
pub fn interpolate_field_value(
    polygon: &CutPolygon,
    xi: f32, 
    eta: f32, 
    zeta: f32,
) -> f64 {
    match polygon.element_type {
        CutElementType::Hex8 => {
            let n = hex8_shape_functions(xi, eta, zeta);
            let mut value = 0.0;
            for i in 0..8.min(polygon.nodal_field_values.len()) {
                value += n[i] * polygon.nodal_field_values[i];
            }
            value
        }
        CutElementType::Hex20 => {
            let n = hex20_shape_functions(xi, eta, zeta);
            let mut value = 0.0;
            for i in 0..20.min(polygon.nodal_field_values.len()) {
                value += n[i] * polygon.nodal_field_values[i];
            }
            value
        }
        CutElementType::Tet4 | CutElementType::Unknown => {
            // Fallback to linear interpolation from corner values
            let n = hex8_shape_functions(xi, eta, zeta);
            let mut value = 0.0;
            for i in 0..8.min(polygon.nodal_field_values.len()) {
                value += n[i] * polygon.nodal_field_values[i];
            }
            value
        }
    }
}

/// Generate a rasterized texture for accurate field visualization using shape functions
/// Works in a local 2D coordinate system on the cut polygon
/// Returns pixel colors for a rectangular region covering the cut polygon
pub fn rasterize_cut_polygon(
    polygon: &CutPolygon,
    resolution: usize,  // Pixels per element side
    min_val: f64,
    max_val: f64,
) -> (Vec<[u8; 4]>, usize, usize, [f32; 2], [f32; 2]) {
    if polygon.vertices.len() < 3 || polygon.node_positions.len() < 8 {
        return (vec![], 0, 0, [0.0, 0.0], [0.0, 0.0]);
    }
    
    // Build a 2D local coordinate system on the polygon plane
    // Use polygon normal and compute tangent vectors
    let n = &polygon.normal;
    
    // Choose a reference direction not parallel to normal
    let ref_dir = if n[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    
    // u_axis = ref_dir x normal (tangent 1)
    let u_axis = [
        ref_dir[1] * n[2] - ref_dir[2] * n[1],
        ref_dir[2] * n[0] - ref_dir[0] * n[2],
        ref_dir[0] * n[1] - ref_dir[1] * n[0],
    ];
    let u_len = (u_axis[0]*u_axis[0] + u_axis[1]*u_axis[1] + u_axis[2]*u_axis[2]).sqrt();
    let u_axis = [u_axis[0]/u_len, u_axis[1]/u_len, u_axis[2]/u_len];
    
    // v_axis = normal x u_axis (tangent 2)
    let v_axis = [
        n[1] * u_axis[2] - n[2] * u_axis[1],
        n[2] * u_axis[0] - n[0] * u_axis[2],
        n[0] * u_axis[1] - n[1] * u_axis[0],
    ];
    
    // Project polygon vertices to 2D local coords
    let center = polygon.center;
    let verts_2d: Vec<[f32; 2]> = polygon.vertices.iter().map(|v| {
        let dx = v.position[0] - center[0];
        let dy = v.position[1] - center[1];
        let dz = v.position[2] - center[2];
        let u = dx * u_axis[0] + dy * u_axis[1] + dz * u_axis[2];
        let v = dx * v_axis[0] + dy * v_axis[1] + dz * v_axis[2];
        [u, v]
    }).collect();
    
    // Compute bounding box in local 2D
    let mut min_u = f32::MAX;
    let mut max_u = f32::MIN;
    let mut min_v = f32::MAX;
    let mut max_v = f32::MIN;
    for p in &verts_2d {
        min_u = min_u.min(p[0]);
        max_u = max_u.max(p[0]);
        min_v = min_v.min(p[1]);
        max_v = max_v.max(p[1]);
    }
    
    let range_u = (max_u - min_u).max(1e-6);
    let range_v = (max_v - min_v).max(1e-6);
    let aspect = range_u / range_v;
    
    let (tex_width, tex_height) = if aspect > 1.0 {
        (resolution, ((resolution as f32) / aspect).max(1.0) as usize)
    } else {
        (((resolution as f32) * aspect).max(1.0) as usize, resolution)
    };
    
    let val_range = (max_val - min_val).max(1e-10);
    
    // Create pixel buffer
    let mut pixels = Vec::with_capacity(tex_width * tex_height);
    
    for j in 0..tex_height {
        for i in 0..tex_width {
            let frac_u = i as f32 / (tex_width - 1).max(1) as f32;
            let frac_v = j as f32 / (tex_height - 1).max(1) as f32;
            
            let local_u = min_u + frac_u * range_u;
            let local_v = min_v + frac_v * range_v;
            
            // Check if point is inside the polygon (in 2D local space)
            if !point_in_polygon_local([local_u, local_v], &verts_2d) {
                pixels.push([0, 0, 0, 0]);
                continue;
            }
            
            // Use generalized barycentric coordinates (mean value coordinates) to interpolate
            // parametric coords from the polygon vertices
            let (xi, eta, zeta) = interpolate_parametric_coords_barycentric(
                [local_u, local_v],
                &verts_2d,
                &polygon.vertices,
            );
            
            // Interpolate field value using shape functions at the parametric coords
            let field_val = interpolate_field_value(polygon, xi, eta, zeta);
            
            // Map to color
            let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
            let color = jet_colormap(t);
            pixels.push(color);
        }
    }
    
    (pixels, tex_width, tex_height, [min_u, min_v], [max_u, max_v])
}

/// Public wrapper for parametric coordinate interpolation
pub fn interpolate_parametric_coords_barycentric_pub(
    point: [f32; 2],
    verts_2d: &[[f32; 2]],
    cut_vertices: &[CutPoint],
) -> (f32, f32, f32) {
    interpolate_parametric_coords_barycentric(point, verts_2d, cut_vertices)
}

/// Interpolate parametric coordinates using generalized barycentric coordinates
/// (Mean Value Coordinates for arbitrary convex polygons)
fn interpolate_parametric_coords_barycentric(
    point: [f32; 2],
    verts_2d: &[[f32; 2]],
    cut_vertices: &[CutPoint],
) -> (f32, f32, f32) {
    let n = verts_2d.len();
    if n < 3 || n != cut_vertices.len() {
        // Fallback
        let xi = cut_vertices.iter().map(|v| v.xi).sum::<f32>() / n as f32;
        let eta = cut_vertices.iter().map(|v| v.eta).sum::<f32>() / n as f32;
        let zeta = cut_vertices.iter().map(|v| v.zeta).sum::<f32>() / n as f32;
        return (xi, eta, zeta);
    }
    
    // Compute Mean Value Coordinates (Floater 2003)
    // w_i = (tan(α_{i-1}/2) + tan(α_i/2)) / r_i
    // where α_i is the angle at point between vertex i and vertex i+1
    // and r_i is the distance from point to vertex i
    
    let mut weights = vec![0.0f32; n];
    let mut total_weight = 0.0f32;
    
    // Compute vectors from point to each vertex
    let mut r = Vec::with_capacity(n);
    let mut d = Vec::with_capacity(n);
    
    for v in verts_2d {
        let dx = v[0] - point[0];
        let dy = v[1] - point[1];
        let dist = (dx*dx + dy*dy).sqrt();
        r.push(dist);
        d.push([dx, dy]);
    }
    
    // Check if point is very close to a vertex
    for i in 0..n {
        if r[i] < 1e-6 {
            return (cut_vertices[i].xi, cut_vertices[i].eta, cut_vertices[i].zeta);
        }
    }
    
    // Compute angles and weights
    for i in 0..n {
        let i_prev = (i + n - 1) % n;
        let i_next = (i + 1) % n;
        
        // Angle at point between vertex i-1 and vertex i
        let cos_alpha_prev = (d[i_prev][0] * d[i][0] + d[i_prev][1] * d[i][1]) / (r[i_prev] * r[i]);
        let alpha_prev = cos_alpha_prev.clamp(-1.0, 1.0).acos();
        
        // Angle at point between vertex i and vertex i+1
        let cos_alpha_next = (d[i][0] * d[i_next][0] + d[i][1] * d[i_next][1]) / (r[i] * r[i_next]);
        let alpha_next = cos_alpha_next.clamp(-1.0, 1.0).acos();
        
        // Mean value weight
        let tan_half_prev = (alpha_prev * 0.5).tan();
        let tan_half_next = (alpha_next * 0.5).tan();
        
        weights[i] = (tan_half_prev + tan_half_next) / r[i];
        total_weight += weights[i];
    }
    
    // Normalize weights and interpolate parametric coordinates
    if total_weight.abs() < 1e-10 {
        // Fallback to simple average
        let xi = cut_vertices.iter().map(|v| v.xi).sum::<f32>() / n as f32;
        let eta = cut_vertices.iter().map(|v| v.eta).sum::<f32>() / n as f32;
        let zeta = cut_vertices.iter().map(|v| v.zeta).sum::<f32>() / n as f32;
        return (xi, eta, zeta);
    }
    
    let mut xi = 0.0f32;
    let mut eta = 0.0f32;
    let mut zeta = 0.0f32;
    
    for i in 0..n {
        let w = weights[i] / total_weight;
        xi += w * cut_vertices[i].xi;
        eta += w * cut_vertices[i].eta;
        zeta += w * cut_vertices[i].zeta;
    }
    
    (xi, eta, zeta)
}

/// Point-in-polygon test using 2D local coordinates
fn point_in_polygon_local(point: [f32; 2], vertices: &[[f32; 2]]) -> bool {
    let n = vertices.len();
    if n < 3 {
        return false;
    }
    
    let mut inside = false;
    let mut j = n - 1;
    
    for i in 0..n {
        let vi = vertices[i];
        let vj = vertices[j];
        
        if ((vi[1] > point[1]) != (vj[1] > point[1])) &&
           (point[0] < (vj[0] - vi[0]) * (point[1] - vi[1]) / (vj[1] - vi[1]) + vi[0])
        {
            inside = !inside;
        }
        j = i;
    }
    
    inside
}

/// Jet colormap: maps 0-1 to blue-cyan-green-yellow-red
fn jet_colormap(t: f32) -> [u8; 4] {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    
    [
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
        255, // Full alpha
    ]
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_hex8_shape_functions_at_corners() {
        // At each corner, only one shape function should be 1.0, others 0.0
        let corners = [
            (-1.0, -1.0, -1.0),
            ( 1.0, -1.0, -1.0),
            ( 1.0,  1.0, -1.0),
            (-1.0,  1.0, -1.0),
            (-1.0, -1.0,  1.0),
            ( 1.0, -1.0,  1.0),
            ( 1.0,  1.0,  1.0),
            (-1.0,  1.0,  1.0),
        ];
        
        for (corner_idx, (xi, eta, zeta)) in corners.iter().enumerate() {
            let n = hex8_shape_functions(*xi as f32, *eta as f32, *zeta as f32);
            
            for (i, &ni) in n.iter().enumerate() {
                if i == corner_idx {
                    assert!((ni - 1.0).abs() < 1e-10, 
                        "At corner {}, N[{}] should be 1.0, got {}", corner_idx, i, ni);
                } else {
                    assert!(ni.abs() < 1e-10, 
                        "At corner {}, N[{}] should be 0.0, got {}", corner_idx, i, ni);
                }
            }
        }
    }
    
    #[test]
    fn test_hex8_shape_functions_at_center() {
        // At center (0,0,0), all shape functions should be 1/8 = 0.125
        let n = hex8_shape_functions(0.0, 0.0, 0.0);
        
        for (i, &ni) in n.iter().enumerate() {
            assert!((ni - 0.125).abs() < 1e-10, 
                "At center, N[{}] should be 0.125, got {}", i, ni);
        }
    }
    
    #[test]
    fn test_hex8_linear_field_interpolation() {
        // If nodal values vary linearly with xi, interpolation should reproduce this
        // Nodal field values for f = xi (ranges from -1 to 1)
        // Corner 0,3,4,7 have xi=-1, corners 1,2,5,6 have xi=1
        let nodal_values = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0f64];
        
        // Test at several interior points
        let test_points = [
            (0.0_f32, 0.0, 0.0, 0.0),   // center: xi=0
            (0.5, 0.0, 0.0, 0.5),        // xi=0.5
            (-0.5, 0.0, 0.0, -0.5),      // xi=-0.5
            (0.25, 0.0, 0.0, 0.25),      // xi=0.25
            (1.0, 0.0, 0.0, 1.0),        // xi=1.0 (edge)
        ];
        
        for (xi, eta, zeta, expected) in test_points {
            let n = hex8_shape_functions(xi, eta, zeta);
            let mut interpolated = 0.0;
            for i in 0..8 {
                interpolated += n[i] * nodal_values[i];
            }
            assert!((interpolated - expected as f64).abs() < 1e-10, 
                "At ({},{},{}): expected {}, got {}", xi, eta, zeta, expected, interpolated);
        }
    }
    
    #[test]
    fn test_shape_functions_partition_of_unity() {
        // Shape functions should sum to 1.0 at any point
        let test_points = [
            (0.0_f32, 0.0, 0.0),
            (0.5, 0.3, -0.2),
            (-0.7, 0.8, 0.1),
            (0.99, -0.99, 0.0),
        ];
        
        for (xi, eta, zeta) in test_points {
            let n = hex8_shape_functions(xi, eta, zeta);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, 
                "At ({},{},{}): shape functions sum to {}, should be 1.0", xi, eta, zeta, sum);
        }
    }
    
    #[test]
    fn test_mean_value_coords_triangle() {
        // For a triangle, mean value coordinates should reduce to barycentric
        let vertices = vec![
            [0.0f32, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
        ];
        
        let cut_vertices = vec![
            CutPoint { position: [0.0, 0.0, 0.0], xi: -1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [1.0, 0.0, 0.0], xi:  1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [0.5, 1.0, 0.0], xi:  0.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
        ];
        
        // At centroid, all barycentric coords should be 1/3
        let centroid = [0.5, 1.0/3.0];
        let (xi, eta, zeta) = interpolate_parametric_coords_barycentric(centroid, &vertices, &cut_vertices);
        
        // Expected: xi = (-1 + 1 + 0)/3 = 0, eta = (-1 + -1 + 1)/3 = -1/3
        let expected_xi = 0.0;
        let expected_eta = -1.0/3.0;
        
        assert!((xi - expected_xi).abs() < 0.1, "xi: expected ~{}, got {}", expected_xi, xi);
        assert!((eta - expected_eta).abs() < 0.1, "eta: expected ~{}, got {}", expected_eta, eta);
    }
    
    #[test]
    fn test_mean_value_coords_at_vertex() {
        // At a vertex, the parametric coords should exactly match that vertex
        let vertices = vec![
            [0.0f32, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ];
        
        let cut_vertices = vec![
            CutPoint { position: [0.0, 0.0, 0.0], xi: -1.0, eta: -1.0, zeta: 0.5, field_value: 0.0 },
            CutPoint { position: [1.0, 0.0, 0.0], xi:  1.0, eta: -1.0, zeta: 0.5, field_value: 0.0 },
            CutPoint { position: [1.0, 1.0, 0.0], xi:  1.0, eta:  1.0, zeta: 0.5, field_value: 0.0 },
            CutPoint { position: [0.0, 1.0, 0.0], xi: -1.0, eta:  1.0, zeta: 0.5, field_value: 0.0 },
        ];
        
        // Test at first vertex
        let (xi, eta, zeta) = interpolate_parametric_coords_barycentric([0.0, 0.0], &vertices, &cut_vertices);
        assert!((xi - (-1.0)).abs() < 1e-5, "xi at vertex 0: expected -1, got {}", xi);
        assert!((eta - (-1.0)).abs() < 1e-5, "eta at vertex 0: expected -1, got {}", eta);
        
        // Test at third vertex
        let (xi, eta, zeta) = interpolate_parametric_coords_barycentric([1.0, 1.0], &vertices, &cut_vertices);
        assert!((xi - 1.0).abs() < 1e-5, "xi at vertex 2: expected 1, got {}", xi);
        assert!((eta - 1.0).abs() < 1e-5, "eta at vertex 2: expected 1, got {}", eta);
    }
    
    #[test]
    fn test_mean_value_coords_center_of_square() {
        // At center of a square, should be average of all corners
        let vertices = vec![
            [0.0f32, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ];
        
        let cut_vertices = vec![
            CutPoint { position: [0.0, 0.0, 0.0], xi: -1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [1.0, 0.0, 0.0], xi:  1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [1.0, 1.0, 0.0], xi:  1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [0.0, 1.0, 0.0], xi: -1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
        ];
        
        let (xi, eta, _zeta) = interpolate_parametric_coords_barycentric([0.5, 0.5], &vertices, &cut_vertices);
        
        // At center of square with these corners, xi and eta should both be 0
        assert!(xi.abs() < 0.1, "xi at center: expected ~0, got {}", xi);
        assert!(eta.abs() < 0.1, "eta at center: expected ~0, got {}", eta);
    }
    
    #[test]
    fn test_rasterize_linear_field_shows_gradient() {
        // Create a cut polygon representing a square section through an element
        // with a linear field varying in x direction
        let cut_vertices = vec![
            CutPoint { position: [0.0, 0.0, 0.0], xi: -1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [1.0, 0.0, 0.0], xi:  1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [1.0, 1.0, 0.0], xi:  1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
            CutPoint { position: [0.0, 1.0, 0.0], xi: -1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
        ];
        
        // Create a polygon with the cut
        // Nodal field values: f = xi varies from -1 to 1
        // corners 0,3 have xi=-1, corners 1,2 have xi=1
        let nodal_field_values = vec![-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0f64];
        
        // Node positions for a unit cube centered at origin (shifted for our polygon)
        let node_positions = vec![
            [-0.5_f32, -0.5, -0.5],  // corner 0 at xi=-1
            [ 0.5, -0.5, -0.5],      // corner 1 at xi=1
            [ 0.5,  0.5, -0.5],      // corner 2 at xi=1
            [-0.5,  0.5, -0.5],      // corner 3 at xi=-1
            [-0.5, -0.5,  0.5],      // corner 4 at xi=-1
            [ 0.5, -0.5,  0.5],      // corner 5 at xi=1
            [ 0.5,  0.5,  0.5],      // corner 6 at xi=1
            [-0.5,  0.5,  0.5],      // corner 7 at xi=-1
        ];
        
        let polygon = CutPolygon {
            vertices: cut_vertices,
            element_id: 0,
            center: [0.5, 0.5, 0.0],
            normal: [0.0, 0.0, 1.0],  // Cut plane normal is Z
            element_type: CutElementType::Hex8,
            nodal_field_values,
            node_positions,
        };
        
        // Rasterize with a decent resolution
        let (pixels, tex_w, tex_h, _, _) = rasterize_cut_polygon(&polygon, 16, -1.0, 1.0);
        
        assert!(tex_w > 0 && tex_h > 0, "Should produce non-empty texture");
        assert!(!pixels.is_empty(), "Should produce pixels");
        
        // Count non-transparent pixels
        let non_transparent: Vec<_> = pixels.iter().filter(|p| p[3] > 0).collect();
        assert!(non_transparent.len() > 0, "Should have some non-transparent pixels");
        
        // Check that there's color variation (not all the same)
        let first_color = non_transparent[0];
        let has_variation = non_transparent.iter().any(|p| p[0] != first_color[0] || p[1] != first_color[1] || p[2] != first_color[2]);
        assert!(has_variation, "Should have color gradient, not flat color");
        
        println!("Rasterization test: {}x{} texture, {} non-transparent pixels, color variation: {}", 
                 tex_w, tex_h, non_transparent.len(), has_variation);
    }
    
    #[test]
    fn test_field_varies_within_element() {
        // This test specifically verifies that field values at different points
        // within the same cut polygon give different interpolated values
        // (proving that we're not just returning a constant)
        
        let polygon = CutPolygon {
            vertices: vec![
                CutPoint { position: [0.0, 0.0, 0.0], xi: -1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
                CutPoint { position: [1.0, 0.0, 0.0], xi:  1.0, eta: -1.0, zeta: 0.0, field_value: 0.0 },
                CutPoint { position: [1.0, 1.0, 0.0], xi:  1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
                CutPoint { position: [0.0, 1.0, 0.0], xi: -1.0, eta:  1.0, zeta: 0.0, field_value: 0.0 },
            ],
            element_id: 0,
            center: [0.5, 0.5, 0.0],
            normal: [0.0, 0.0, 1.0],
            element_type: CutElementType::Hex8,
            // Nodal values representing f = xi (linear in x)
            nodal_field_values: vec![-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0],
            node_positions: vec![
                [-0.5_f32, -0.5, -0.5], [ 0.5, -0.5, -0.5],
                [ 0.5,  0.5, -0.5], [-0.5,  0.5, -0.5],
                [-0.5, -0.5,  0.5], [ 0.5, -0.5,  0.5],
                [ 0.5,  0.5,  0.5], [-0.5,  0.5,  0.5],
            ],
        };
        
        // Test at the center: should get xi=0, field=0
        let field_center = interpolate_field_value(&polygon, 0.0, 0.0, 0.0);
        assert!(field_center.abs() < 0.01, "At center, field should be ~0, got {}", field_center);
        
        // Test near the right edge: should get positive field
        let field_right = interpolate_field_value(&polygon, 0.8, 0.0, 0.0);
        assert!(field_right > 0.5, "At xi=0.8, field should be ~0.8, got {}", field_right);
        
        // Test near the left edge: should get negative field
        let field_left = interpolate_field_value(&polygon, -0.8, 0.0, 0.0);
        assert!(field_left < -0.5, "At xi=-0.8, field should be ~-0.8, got {}", field_left);
        
        // Verify monotonic variation
        assert!(field_left < field_center && field_center < field_right,
                "Field should increase from left to right: {} < {} < {}", 
                field_left, field_center, field_right);
                
        println!("Field interpolation test passed: left={:.3}, center={:.3}, right={:.3}", 
                 field_left, field_center, field_right);
    }
    
    #[test]
    fn test_cut_polygon_with_y_varying_field() {
        // Simulate a cut through a hex element with an X-normal plane
        // For an X-normal cut, xi should be approximately constant for all cut vertices
        // but eta and zeta should vary - so a field varying with eta should show gradient
        
        // A cut through a unit cube at x=0.3 (xi ≈ -0.4)
        let t = 0.3_f32;
        let xi_cut = -1.0 + t * 2.0;  // = -0.4
        
        // Cut vertex parametric coords (all have same xi from the cut plane)
        let cut_vertices = vec![
            CutPoint { position: [0.3, 0.0, 0.0], xi: xi_cut, eta: -1.0, zeta: -1.0, field_value: 0.0 },
            CutPoint { position: [0.3, 1.0, 0.0], xi: xi_cut, eta:  1.0, zeta: -1.0, field_value: 0.0 },
            CutPoint { position: [0.3, 1.0, 1.0], xi: xi_cut, eta:  1.0, zeta:  1.0, field_value: 0.0 },
            CutPoint { position: [0.3, 0.0, 1.0], xi: xi_cut, eta: -1.0, zeta:  1.0, field_value: 0.0 },
        ];
        
        let polygon = CutPolygon {
            vertices: cut_vertices.clone(),
            element_id: 0,
            center: [0.3, 0.5, 0.5],
            normal: [1.0, 0.0, 0.0],  // X-normal
            element_type: CutElementType::Hex8,
            // Field varies with eta (y direction): corners at y=0 have eta=-1, at y=1 have eta=1
            nodal_field_values: vec![
                -1.0, -1.0, 1.0, 1.0,  // z=0 face: corners 0,1 have eta=-1, corners 2,3 have eta=1
                -1.0, -1.0, 1.0, 1.0,  // z=1 face: corners 4,5 have eta=-1, corners 6,7 have eta=1
            ],
            node_positions: vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],  // y=0, z=0
                [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],  // y=1, z=0
                [0.0, 0.0, 1.0], [1.0, 0.0, 1.0],  // y=0, z=1
                [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],  // y=1, z=1
            ],
        };
        
        // Test interpolation at different eta values on the cut plane
        let field_bottom = interpolate_field_value(&polygon, xi_cut, -1.0, 0.0);
        let field_mid = interpolate_field_value(&polygon, xi_cut, 0.0, 0.0);
        let field_top = interpolate_field_value(&polygon, xi_cut, 1.0, 0.0);
        
        println!("Y-varying field test:");
        println!("  At eta=-1: {:.3}", field_bottom);
        println!("  At eta=0: {:.3}", field_mid);
        println!("  At eta=1: {:.3}", field_top);
        
        assert!(field_bottom < -0.8, "Field at eta=-1 should be ~-1: {}", field_bottom);
        assert!(field_mid.abs() < 0.2, "Field at eta=0 should be ~0: {}", field_mid);
        assert!(field_top > 0.8, "Field at eta=1 should be ~1: {}", field_top);
        
        // Now test the full rasterization
        let (pixels, tex_w, tex_h, _, _) = rasterize_cut_polygon(&polygon, 16, -1.0, 1.0);
        
        // Count unique colors to see if there's actual variation
        let non_transparent: Vec<_> = pixels.iter().filter(|p| p[3] > 0).collect();
        let unique_colors: std::collections::HashSet<_> = non_transparent.iter().map(|p| (p[0], p[1], p[2])).collect();
        
        println!("  Texture: {}x{}, {} non-transparent, {} unique colors", 
                 tex_w, tex_h, non_transparent.len(), unique_colors.len());
        
        assert!(unique_colors.len() > 1, "Should have multiple colors in the gradient, got {}", unique_colors.len());
    }
}
