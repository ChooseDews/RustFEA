use log::{debug, warn};
use nalgebra::DVector;
use rayon::iter::IntoParallelRefIterator;
use std::sync::{mpsc, Arc};
use std::thread;

use super::Simulation;

impl Simulation {
    pub fn compute_force_vector(&mut self, displacement: &DVector<f64>) -> DVector<f64> {
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

    pub fn compute_force_vector_threaded(&mut self, displacement: &DVector<f64>) -> DVector<f64> {
        let active_ids = self.active_elements();
        let elements_count = active_ids.len();
        let elements = Arc::new(std::mem::take(&mut self.elements));
        let mut handles = Vec::new();
        let chunk_size = (active_ids.len() + self.worker_count - 1) / self.worker_count;
        let chunks = self
            .active_elements()
            .chunks(chunk_size)
            .map(|c| c.to_vec())
            .collect::<Vec<_>>();
        let dofs = self.dofs;
        let displacement = Arc::new(displacement.clone());
        let n_dofs = self.nodes.len() * dofs;

        for chunk in chunks {
            let elements = Arc::clone(&elements);
            let mut force_vector = DVector::zeros(n_dofs);
            let displacement = Arc::clone(&displacement);
            let handle = thread::spawn(move || {
                for &id in &chunk {
                    let element = &elements[&id];
                    let connectivity = element.get_connectivity();
                    let f = element.compute_force(&displacement);
                    for (i, &node_id) in connectivity.iter().enumerate() {
                        for dof in 0..dofs {
                            force_vector[node_id * dofs + dof] += f[i * dofs + dof];
                        }
                    }
                }
                force_vector
            });
            handles.push(handle);
        }
        let force_vector = handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .fold(DVector::zeros(n_dofs), |mut acc, result| {
                acc += result;
                acc
            });

        self.elements = Arc::try_unwrap(elements).ok().unwrap();
        force_vector
    }
}
