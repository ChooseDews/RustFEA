//! GPU-accelerated mesh renderer using wgpu
//! 
//! This module provides high-performance 3D rendering using wgpu.
//! For the initial implementation, we use egui's software renderer in viewport.rs.
//! This module can be expanded for GPU acceleration when needed.

use crate::state::{MeshState, Vertex, MeshRenderData, CameraState, ColorMode};

/// GPU mesh renderer (placeholder for future wgpu implementation)
pub struct MeshRenderer {
    // TODO: wgpu resources
    // device: wgpu::Device,
    // queue: wgpu::Queue,
    // pipeline: wgpu::RenderPipeline,
    // vertex_buffer: wgpu::Buffer,
    // index_buffer: wgpu::Buffer,
}

impl MeshRenderer {
    pub fn new() -> Self {
        Self {
            // TODO: Initialize wgpu
        }
    }
    
    /// Generate render data from mesh
    pub fn generate_render_data(mesh_state: &MeshState, results: Option<&crate::state::SimulationResults>, color_mode: ColorMode, scale: f32) -> MeshRenderData {
        let mesh = &mesh_state.mesh;
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        let mut wireframe_indices = Vec::new();
        
        // Process each element
        for (_el_id, element) in &mesh.elements {
            let conn = &element.connectivity;
            
            if conn.len() < 8 {
                continue; // Skip non-brick elements for now
            }
            
            // Get base vertex index
            let base_idx = vertices.len() as u32;
            
            // Add vertices for this element
            for (_local_idx, &node_id) in conn.iter().enumerate() {
                if let Some(node) = mesh.nodes.get(&node_id) {
                    let mut pos = [
                        node.coordinates[0] as f32,
                        node.coordinates[1] as f32,
                        node.coordinates[2] as f32,
                    ];
                    
                    // Apply displacement
                    if let Some(res) = results {
                        if node_id * 3 + 2 < res.displacements.len() {
                            pos[0] += res.displacements[node_id * 3] as f32 * scale;
                            pos[1] += res.displacements[node_id * 3 + 1] as f32 * scale;
                            pos[2] += res.displacements[node_id * 3 + 2] as f32 * scale;
                        }
                    }
                    
                    // Compute color
                    let color = match color_mode {
                        ColorMode::Solid => [0.4, 0.6, 0.9, 1.0], // Cornflower blue
                        ColorMode::Displacement => {
                            if let Some(res) = results {
                                if node_id * 3 + 2 < res.displacements.len() {
                                    let dx = res.displacements[node_id * 3];
                                    let dy = res.displacements[node_id * 3 + 1];
                                    let dz = res.displacements[node_id * 3 + 2];
                                    let mag = (dx * dx + dy * dy + dz * dz).sqrt();
                                    let t = (mag / res.stats.max_displacement.max(1e-10)) as f32;
                                    value_to_color_array(t.clamp(0.0, 1.0))
                                } else {
                                    [0.5, 0.5, 0.5, 1.0]
                                }
                            } else {
                                [0.5, 0.5, 0.5, 1.0]
                            }
                        }
                        _ => [0.4, 0.6, 0.9, 1.0],
                    };
                    
                    // Placeholder normal (should compute from geometry)
                    let normal = [0.0, 1.0, 0.0];
                    
                    vertices.push(Vertex::new(pos, normal, color));
                }
            }
            
            // Face indices for 8-node brick
            let faces = [
                [0, 1, 2, 3], // Bottom
                [4, 7, 6, 5], // Top
                [0, 4, 5, 1], // Front
                [2, 6, 7, 3], // Back
                [0, 3, 7, 4], // Left
                [1, 5, 6, 2], // Right
            ];
            
            for face in &faces {
                // Two triangles per quad
                indices.push(base_idx + face[0] as u32);
                indices.push(base_idx + face[1] as u32);
                indices.push(base_idx + face[2] as u32);
                
                indices.push(base_idx + face[0] as u32);
                indices.push(base_idx + face[2] as u32);
                indices.push(base_idx + face[3] as u32);
                
                // Wireframe edges
                wireframe_indices.push(base_idx + face[0] as u32);
                wireframe_indices.push(base_idx + face[1] as u32);
                wireframe_indices.push(base_idx + face[1] as u32);
                wireframe_indices.push(base_idx + face[2] as u32);
                wireframe_indices.push(base_idx + face[2] as u32);
                wireframe_indices.push(base_idx + face[3] as u32);
                wireframe_indices.push(base_idx + face[3] as u32);
                wireframe_indices.push(base_idx + face[0] as u32);
            }
        }
        
        MeshRenderData {
            vertices,
            indices,
            wireframe_indices,
        }
    }
}

fn value_to_color_array(t: f32) -> [f32; 4] {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    [r, g, b, 1.0]
}

/// View-projection matrix computation
pub fn compute_view_matrix(camera: &CameraState) -> [[f32; 4]; 4] {
    let eye = camera.eye_position();
    let target = camera.target;
    let up = [0.0, 1.0, 0.0];
    
    // Compute look-at matrix
    let f = [
        target[0] - eye[0],
        target[1] - eye[1],
        target[2] - eye[2],
    ];
    let f_len = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
    let f = [f[0] / f_len, f[1] / f_len, f[2] / f_len];
    
    let s = [
        f[1] * up[2] - f[2] * up[1],
        f[2] * up[0] - f[0] * up[2],
        f[0] * up[1] - f[1] * up[0],
    ];
    let s_len = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    let s = [s[0] / s_len, s[1] / s_len, s[2] / s_len];
    
    let u = [
        s[1] * f[2] - s[2] * f[1],
        s[2] * f[0] - s[0] * f[2],
        s[0] * f[1] - s[1] * f[0],
    ];
    
    [
        [s[0], u[0], -f[0], 0.0],
        [s[1], u[1], -f[1], 0.0],
        [s[2], u[2], -f[2], 0.0],
        [
            -s[0] * eye[0] - s[1] * eye[1] - s[2] * eye[2],
            -u[0] * eye[0] - u[1] * eye[1] - u[2] * eye[2],
            f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2],
            1.0,
        ],
    ]
}

pub fn compute_projection_matrix(fov: f32, aspect: f32, near: f32, far: f32) -> [[f32; 4]; 4] {
    let f = 1.0 / (fov / 2.0).tan();
    
    [
        [f / aspect, 0.0, 0.0, 0.0],
        [0.0, f, 0.0, 0.0],
        [0.0, 0.0, (far + near) / (near - far), -1.0],
        [0.0, 0.0, (2.0 * far * near) / (near - far), 0.0],
    ]
}
