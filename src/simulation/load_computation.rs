use log::{debug, warn};
use nalgebra::DVector;
#[cfg(feature = "native")]
use rayon::prelude::*;
#[cfg(feature = "native")]
use std::sync::Arc;

use super::Simulation;

impl Simulation {
    pub fn compute_force_vector(&mut self, displacement: &DVector<f64>) -> DVector<f64> {
        #[cfg(feature = "native")]
        if self.worker_count > 1 {
            return self.compute_force_vector_threaded(displacement);
        }
        self.compute_force_vector_single(displacement)
    }

    pub fn compute_force_vector_single(&mut self, displacement: &DVector<f64>) -> DVector<f64> {
        let active_ids = self.active_elements();
        let mut force_vector = DVector::zeros(self.nodes.len() * self.dofs);
        let elements = std::mem::take(&mut self.elements);
        for active_id in active_ids {
            elements[&active_id].add_force(self, &mut force_vector);
        }
        self.elements = elements;
        force_vector
    }

    #[cfg(feature = "native")]
    pub fn compute_force_vector_threaded(&mut self, displacement: &DVector<f64>) -> DVector<f64> {
        let active_ids = self.active_elements();
        let elements = Arc::new(std::mem::take(&mut self.elements));
        let dofs = self.dofs;
        let n_dofs = self.nodes.len() * dofs;

        // Use rayon parallel iterator with map-reduce pattern
        let force_vector = active_ids
            .par_iter()
            .map(|&id| {
                let element = &elements[&id];
                let connectivity = element.get_connectivity();
                let f = element.compute_force(displacement);
                
                // Create a sparse contribution for this element
                let mut local_force = DVector::zeros(n_dofs);
                for (i, &node_id) in connectivity.iter().enumerate() {
                    for dof in 0..dofs {
                        local_force[node_id * dofs + dof] += f[i * dofs + dof];
                    }
                }
                local_force
            })
            .reduce(
                || DVector::zeros(n_dofs),
                |mut acc, local| {
                    acc += local;
                    acc
                },
            );

        self.elements = Arc::try_unwrap(elements).ok().unwrap();
        force_vector
    }
}
