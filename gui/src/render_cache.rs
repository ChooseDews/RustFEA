//! Render cache for optimized viewport rendering
//!
//! This module provides caching and spatial acceleration to avoid
//! recomputing geometry and colors every frame.

use crate::state::{CameraState, ColorMode, SimulationResults, MeshState};
use eframe::egui;
use std::collections::HashMap;

/// Cached face data for rendering
#[derive(Clone)]
#[allow(dead_code)]
pub struct CachedFace {
    /// World-space corner positions (4 corners for quad face)
    pub corners: [[f32; 3]; 4],
    /// Precomputed face normal
    pub normal: [f32; 3],
    /// Face center in world space
    pub center: [f32; 3],
    /// Element ID this face belongs to
    pub element_id: usize,
    /// Face index within the element (0-5 for hex element)
    pub face_index: usize,
    /// Local node indices within the element
    pub node_indices: [usize; 4],
    /// Node IDs for color lookup
    pub node_ids: [usize; 4],
}

/// Cached colors for nodes/elements - avoids recomputation every frame
#[derive(Clone)]
#[allow(dead_code)]
pub struct ColorCache {
    /// Node colors by node_id (for displacement mode)
    pub node_colors: HashMap<usize, egui::Color32>,
    /// Element colors by element_id (for stress/strain/vm modes)
    pub element_colors: HashMap<usize, egui::Color32>,
    /// Color mode these were computed for
    pub color_mode: ColorMode,
    /// Stress/strain component these were computed for
    pub stress_component: usize,
    pub strain_component: usize,
    /// Results version (to detect changes)
    pub results_version: u64,
}

impl Default for ColorCache {
    fn default() -> Self {
        Self {
            node_colors: HashMap::new(),
            element_colors: HashMap::new(),
            color_mode: ColorMode::Solid,
            stress_component: 0,
            strain_component: 0,
            results_version: 0,
        }
    }
}

/// Axis-aligned bounding box for frustum culling
#[derive(Clone, Copy, Debug)]
pub struct AABB {
    pub min: [f32; 3],
    pub max: [f32; 3],
}

impl AABB {
    pub fn new(min: [f32; 3], max: [f32; 3]) -> Self {
        Self { min, max }
    }
    
    pub fn center(&self) -> [f32; 3] {
        [
            (self.min[0] + self.max[0]) * 0.5,
            (self.min[1] + self.max[1]) * 0.5,
            (self.min[2] + self.max[2]) * 0.5,
        ]
    }
    
    pub fn radius(&self) -> f32 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        (dx * dx + dy * dy + dz * dz).sqrt() * 0.5
    }
    
    /// Check if AABB is potentially visible from camera
    /// Uses conservative sphere-frustum test
    pub fn is_potentially_visible(&self, camera: &CameraState, _viewport_aspect: f32) -> bool {
        let eye = camera.eye_position();
        let center = self.center();
        let radius = self.radius();
        
        // Vector from eye to AABB center
        let to_center = [
            center[0] - eye[0],
            center[1] - eye[1],
            center[2] - eye[2],
        ];
        
        // Distance to center
        let dist_sq = to_center[0] * to_center[0] + 
                      to_center[1] * to_center[1] + 
                      to_center[2] * to_center[2];
        let dist = dist_sq.sqrt();
        
        // If too close to camera, always render
        if dist < radius * 2.0 {
            return true;
        }
        
        // Get forward direction
        let forward = [
            camera.target[0] - eye[0],
            camera.target[1] - eye[1],
            camera.target[2] - eye[2],
        ];
        let forward_len = (forward[0] * forward[0] + forward[1] * forward[1] + forward[2] * forward[2]).sqrt();
        if forward_len < 1e-6 {
            return true;
        }
        let forward = [forward[0] / forward_len, forward[1] / forward_len, forward[2] / forward_len];
        
        // Dot product to check if in front of camera
        let dot = to_center[0] * forward[0] + to_center[1] * forward[1] + to_center[2] * forward[2];
        
        // Behind camera check (with margin for bounding sphere)
        // Be conservative - only reject if clearly behind
        if dot < -radius * 2.0 {
            return false;
        }
        
        // For objects in front of camera, be very conservative
        // Only reject if the object is way outside the FOV
        // This avoids false negatives from complex angle calculations
        let cos_angle = (dot / dist).clamp(-1.0, 1.0);
        
        // If center is roughly in front of camera (within ~120 degrees), render it
        // cos(120°) = -0.5, so if cos_angle > -0.5 we're within 120° of view direction
        cos_angle > -0.5
    }
}

/// Simple octree node for spatial acceleration
pub struct OctreeNode {
    pub bounds: AABB,
    pub faces: Vec<usize>,  // Indices into face array
    pub children: Option<Box<[OctreeNode; 8]>>,
}

impl OctreeNode {
    pub fn new(bounds: AABB) -> Self {
        Self {
            bounds,
            faces: Vec::new(),
            children: None,
        }
    }
    
    /// Get octant index for a point
    fn octant_index(center: [f32; 3], point: [f32; 3]) -> usize {
        let mut idx = 0;
        if point[0] >= center[0] { idx |= 1; }
        if point[1] >= center[1] { idx |= 2; }
        if point[2] >= center[2] { idx |= 4; }
        idx
    }
    
    /// Get child bounds for an octant
    fn child_bounds(parent: &AABB, octant: usize) -> AABB {
        let center = parent.center();
        let mut min = parent.min;
        let mut max = parent.max;
        
        if octant & 1 != 0 { min[0] = center[0]; } else { max[0] = center[0]; }
        if octant & 2 != 0 { min[1] = center[1]; } else { max[1] = center[1]; }
        if octant & 4 != 0 { min[2] = center[2]; } else { max[2] = center[2]; }
        
        AABB::new(min, max)
    }
    
    /// Insert a face into the octree
    pub fn insert(&mut self, face_idx: usize, face_center: [f32; 3], max_depth: u32, current_depth: u32) {
        // If at max depth or few faces, store here
        if current_depth >= max_depth || (self.children.is_none() && self.faces.len() < 16) {
            self.faces.push(face_idx);
            return;
        }
        
        // Create children if needed
        if self.children.is_none() {
            let center = self.bounds.center();
            let children = [
                OctreeNode::new(Self::child_bounds(&self.bounds, 0)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 1)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 2)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 3)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 4)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 5)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 6)),
                OctreeNode::new(Self::child_bounds(&self.bounds, 7)),
            ];
            self.children = Some(Box::new(children));
            
            // Move existing faces to children
            let old_faces = std::mem::take(&mut self.faces);
            // Note: We'd need face centers for these too - simplified: keep in parent
            self.faces = old_faces;
        }
        
        // Insert into appropriate child
        if let Some(ref mut children) = self.children {
            let octant = Self::octant_index(self.bounds.center(), face_center);
            children[octant].insert(face_idx, face_center, max_depth, current_depth + 1);
        }
    }
    
    /// Collect all visible face indices
    pub fn collect_visible(&self, camera: &CameraState, aspect: f32, out: &mut Vec<usize>) {
        if !self.bounds.is_potentially_visible(camera, aspect) {
            return;
        }
        
        // Add all faces in this node
        out.extend(self.faces.iter().copied());
        
        // Recurse into children
        if let Some(ref children) = self.children {
            for child in children.iter() {
                child.collect_visible(camera, aspect, out);
            }
        }
    }
}

/// Main render cache for the viewport
#[allow(dead_code)]
pub struct RenderCache {
    /// Cached face geometry
    pub faces: Vec<CachedFace>,
    /// Spatial acceleration structure
    pub octree: Option<OctreeNode>,
    /// Color cache
    pub colors: ColorCache,
    /// Displacement scale used for cached positions
    pub cached_displacement_scale: f32,
    /// Mesh version (to detect mesh changes)
    pub mesh_version: u64,
    /// Whether face normals point outward (for consistent backface culling)
    pub normals_validated: bool,
    
    // === Persistent buffers to reduce allocations ===
    /// Reusable buffer for visible face indices
    pub visible_faces_buffer: Vec<usize>,
    /// Unique edges for wireframe (to avoid drawing shared edges multiple times)
    pub unique_edges: Vec<([f32; 3], [f32; 3])>,
    /// Edge to face mapping for quick lookup
    pub edge_face_map: std::collections::HashSet<(u64, u64)>,
    
    // === Frame timing for performance monitoring ===
    /// Number of faces rendered last frame (for stats display)
    pub last_rendered_faces: usize,
    /// Number of triangles rendered last frame
    pub last_rendered_triangles: usize,
    
    // === Cached view transform for projection ===
    /// Cached view transform (avoids recomputing per-vertex)
    pub view_transform: Option<ViewTransform>,
}

/// Cached view transform for fast projection
#[derive(Clone)]
#[allow(dead_code)]
pub struct ViewTransform {
    pub eye: [f32; 3],
    pub forward: [f32; 3],
    pub right: [f32; 3],
    pub up: [f32; 3],
    pub fov_factor: f32,
    pub near_plane: f32,
    /// Camera state hash for change detection
    pub camera_hash: u64,
    /// Use orthographic projection
    pub orthographic: bool,
    /// Camera distance (for orthographic scale)
    pub distance: f32,
}

impl ViewTransform {
    /// Create a view transform from camera state
    pub fn from_camera(camera: &CameraState) -> Option<Self> {
        let eye = camera.eye_position();
        
        // View direction
        let view_x = camera.target[0] - eye[0];
        let view_y = camera.target[1] - eye[1];
        let view_z = camera.target[2] - eye[2];
        let view_len = (view_x * view_x + view_y * view_y + view_z * view_z).sqrt();
        
        if view_len < 1e-6 {
            return None; // Eye at target
        }
        
        // Normalize view direction
        let forward = [view_x / view_len, view_y / view_len, view_z / view_len];
        
        // Choose up vector that's not parallel to forward
        let world_up = if forward[1].abs() > 0.99 {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        
        // Right vector = forward × up
        let right = [
            forward[1] * world_up[2] - forward[2] * world_up[1],
            forward[2] * world_up[0] - forward[0] * world_up[2],
            forward[0] * world_up[1] - forward[1] * world_up[0],
        ];
        let right_len = (right[0] * right[0] + right[1] * right[1] + right[2] * right[2]).sqrt();
        
        if right_len < 1e-6 {
            return None;
        }
        
        let right = [right[0] / right_len, right[1] / right_len, right[2] / right_len];
        
        // True up = right × forward
        let up = [
            right[1] * forward[2] - right[2] * forward[1],
            right[2] * forward[0] - right[0] * forward[2],
            right[0] * forward[1] - right[1] * forward[0],
        ];
        
        let fov_factor = (camera.fov / 2.0).tan();
        let near_plane = (camera.distance * 0.01).max(0.001);
        
        // Simple hash of camera state for change detection
        let camera_hash = {
            let bits = |f: f32| -> u64 { f.to_bits() as u64 };
            bits(camera.yaw) ^ bits(camera.pitch).rotate_left(16) 
                ^ bits(camera.distance).rotate_left(32) 
                ^ bits(camera.target[0]).rotate_left(48)
                ^ if camera.orthographic { 1 } else { 0 }
        };
        
        Some(Self {
            eye,
            forward,
            right,
            up,
            fov_factor,
            near_plane,
            camera_hash,
            orthographic: camera.orthographic,
            distance: camera.distance,
        })
    }
    
    /// Project a 3D point using the cached transform
    #[inline]
    pub fn project(&self, point: [f32; 3], rect_center: egui::Pos2, half_width: f32, half_height: f32, aspect: f32) -> Option<egui::Pos2> {
        // Vector from eye to point
        let rel = [
            point[0] - self.eye[0],
            point[1] - self.eye[1],
            point[2] - self.eye[2],
        ];
        
        // Camera space coordinates
        let cam_x = rel[0] * self.right[0] + rel[1] * self.right[1] + rel[2] * self.right[2];
        let cam_y = rel[0] * self.up[0] + rel[1] * self.up[1] + rel[2] * self.up[2];
        let cam_z = rel[0] * self.forward[0] + rel[1] * self.forward[1] + rel[2] * self.forward[2];
        
        // Behind camera or near plane (for perspective) - orthographic shows everything
        if !self.orthographic && cam_z < self.near_plane {
            return None;
        }
        
        let (ndc_x, ndc_y) = if self.orthographic {
            // Orthographic projection: scale by camera distance to maintain size
            // The view size is distance * tan(fov/2) * 2 for perspective, use same scale
            let ortho_scale = 1.0 / (self.distance * self.fov_factor);
            (cam_x * ortho_scale / aspect, cam_y * ortho_scale)
        } else {
            // Perspective projection
            let inv_z = 1.0 / (cam_z * self.fov_factor);
            (cam_x * inv_z / aspect, cam_y * inv_z)
        };
        
        // Convert to screen coordinates
        let screen_x = rect_center.x + ndc_x * half_width;
        let screen_y = rect_center.y - ndc_y * half_height;
        
        // Reject extreme values
        let max_extent = half_width.max(half_height) * 4.0;
        if screen_x.abs() > max_extent + rect_center.x.abs() || 
           screen_y.abs() > max_extent + rect_center.y.abs() ||
           !screen_x.is_finite() || !screen_y.is_finite() {
            return None;
        }
        
        Some(egui::pos2(screen_x, screen_y))
    }
    
    /// Project a 3D point with a very small near plane - for section cuts and clip planes
    /// that need to be visible even when very close to the camera.
    #[inline]
    pub fn project_close(&self, point: [f32; 3], rect_center: egui::Pos2, half_width: f32, half_height: f32, aspect: f32) -> Option<egui::Pos2> {
        // Vector from eye to point
        let rel = [
            point[0] - self.eye[0],
            point[1] - self.eye[1],
            point[2] - self.eye[2],
        ];
        
        // Camera space coordinates
        let cam_x = rel[0] * self.right[0] + rel[1] * self.right[1] + rel[2] * self.right[2];
        let cam_y = rel[0] * self.up[0] + rel[1] * self.up[1] + rel[2] * self.up[2];
        let cam_z = rel[0] * self.forward[0] + rel[1] * self.forward[1] + rel[2] * self.forward[2];
        
        // Use a much smaller near plane for section cuts (0.0001 instead of distance*0.01)
        // Orthographic mode has no near plane culling
        let min_near = 0.0001;
        if !self.orthographic && cam_z < min_near {
            return None;
        }
        
        let (ndc_x, ndc_y) = if self.orthographic {
            let ortho_scale = 1.0 / (self.distance * self.fov_factor);
            (cam_x * ortho_scale / aspect, cam_y * ortho_scale)
        } else {
            // For very close points, clamp cam_z to avoid extreme perspective distortion
            let safe_z = cam_z.max(min_near);
            let inv_z = 1.0 / (safe_z * self.fov_factor);
            (cam_x * inv_z / aspect, cam_y * inv_z)
        };
        
        // Convert to screen coordinates
        let screen_x = rect_center.x + ndc_x * half_width;
        let screen_y = rect_center.y - ndc_y * half_height;
        
        // Allow a larger extent for close objects that may appear large on screen
        let max_extent = half_width.max(half_height) * 10.0;
        if screen_x.abs() > max_extent + rect_center.x.abs() || 
           screen_y.abs() > max_extent + rect_center.y.abs() ||
           !screen_x.is_finite() || !screen_y.is_finite() {
            return None;
        }
        
        Some(egui::pos2(screen_x, screen_y))
    }
}

impl Default for RenderCache {
    fn default() -> Self {
        Self {
            faces: Vec::new(),
            octree: None,
            colors: ColorCache::default(),
            cached_displacement_scale: 1.0,
            mesh_version: 0,
            normals_validated: false,
            visible_faces_buffer: Vec::new(),
            unique_edges: Vec::new(),
            edge_face_map: std::collections::HashSet::new(),
            last_rendered_faces: 0,
            last_rendered_triangles: 0,
            view_transform: None,
        }
    }
}

impl RenderCache {
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Invalidate the cache, forcing a rebuild on next render
    pub fn invalidate(&mut self) {
        self.faces.clear();
        self.octree = None;
        self.mesh_version = self.mesh_version.wrapping_add(1);
    }
    
    /// Check if cache is valid for current state
    pub fn is_valid(
        &self,
        mesh_state: &MeshState,
        results: &Option<SimulationResults>,
        displacement_scale: f32,
        color_mode: ColorMode,
        stress_component: usize,
        strain_component: usize,
    ) -> bool {
        // Check mesh hasn't changed
        if self.faces.is_empty() {
            return false;
        }
        
        // Check displacement scale (affects geometry)
        if results.is_some() && (self.cached_displacement_scale - displacement_scale).abs() > 1e-6 {
            return false;
        }
        
        // Check color settings
        if self.colors.color_mode != color_mode {
            return false;
        }
        if color_mode == ColorMode::Stress && self.colors.stress_component != stress_component {
            return false;
        }
        if color_mode == ColorMode::Strain && self.colors.strain_component != strain_component {
            return false;
        }
        
        true
    }
    
    /// Rebuild cache from mesh
    pub fn rebuild(
        &mut self,
        mesh_state: &MeshState,
        results: &Option<SimulationResults>,
        displacement_scale: f32,
        color_mode: ColorMode,
        stress_component: usize,
        strain_component: usize,
    ) {
        self.faces.clear();
        self.octree = None;
        self.cached_displacement_scale = displacement_scale;
        
        let mesh = &mesh_state.mesh;
        
        // Face ordering for 8-node brick element
        let face_indices = [
            [0, 1, 2, 3], // Bottom
            [4, 7, 6, 5], // Top
            [0, 4, 5, 1], // Front
            [2, 6, 7, 3], // Back
            [0, 3, 7, 4], // Left
            [1, 5, 6, 2], // Right
        ];
        
        // Extract all faces from all elements
        for (&el_id, element) in &mesh.elements {
            let conn = &element.connectivity;
            if conn.len() < 8 {
                continue;
            }
            
            // Get node positions
            let mut positions: Vec<[f32; 3]> = Vec::with_capacity(8);
            for &node_id in conn.iter().take(8) {
                if let Some(node) = mesh.nodes.get(&node_id) {
                    let mut pos = [
                        node.coordinates[0] as f32,
                        node.coordinates[1] as f32,
                        node.coordinates[2] as f32,
                    ];
                    
                    // Apply displacement if available
                    if let Some(res) = results {
                        if node_id * 3 + 2 < res.displacements.len() {
                            pos[0] += res.displacements[node_id * 3] as f32 * displacement_scale;
                            pos[1] += res.displacements[node_id * 3 + 1] as f32 * displacement_scale;
                            pos[2] += res.displacements[node_id * 3 + 2] as f32 * displacement_scale;
                        }
                    }
                    positions.push(pos);
                }
            }
            
            if positions.len() < 8 {
                continue;
            }
            
            // Create faces
            for (face_idx, face) in face_indices.iter().enumerate() {
                let corners = [
                    positions[face[0]],
                    positions[face[1]],
                    positions[face[2]],
                    positions[face[3]],
                ];
                
                // Compute face normal
                let v1 = [
                    corners[1][0] - corners[0][0],
                    corners[1][1] - corners[0][1],
                    corners[1][2] - corners[0][2],
                ];
                let v2 = [
                    corners[2][0] - corners[0][0],
                    corners[2][1] - corners[0][1],
                    corners[2][2] - corners[0][2],
                ];
                let normal = [
                    v1[1] * v2[2] - v1[2] * v2[1],
                    v1[2] * v2[0] - v1[0] * v2[2],
                    v1[0] * v2[1] - v1[1] * v2[0],
                ];
                
                // Normalize
                let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
                let normal = if len > 1e-6 {
                    [normal[0] / len, normal[1] / len, normal[2] / len]
                } else {
                    [0.0, 1.0, 0.0]
                };
                
                let center = [
                    (corners[0][0] + corners[1][0] + corners[2][0] + corners[3][0]) * 0.25,
                    (corners[0][1] + corners[1][1] + corners[2][1] + corners[3][1]) * 0.25,
                    (corners[0][2] + corners[1][2] + corners[2][2] + corners[3][2]) * 0.25,
                ];
                
                let node_ids = [
                    conn[face[0]],
                    conn[face[1]],
                    conn[face[2]],
                    conn[face[3]],
                ];
                
                self.faces.push(CachedFace {
                    corners,
                    normal,
                    center,
                    element_id: el_id,
                    face_index: face_idx,
                    node_indices: *face,
                    node_ids,
                });
            }
        }
        
        // Build octree if mesh is large enough
        if self.faces.len() > 100 {
            let bounds = AABB::new(mesh_state.bounds.min, mesh_state.bounds.max);
            let mut octree = OctreeNode::new(bounds);
            
            for (i, face) in self.faces.iter().enumerate() {
                octree.insert(i, face.center, 6, 0);
            }
            
            self.octree = Some(octree);
        }
        
        // Build unique edges for wireframe optimization
        self.build_unique_edges();
        
        // Rebuild colors
        self.rebuild_colors(mesh_state, results, color_mode, stress_component, strain_component);
    }
    
    /// Build unique edges to avoid drawing shared edges multiple times
    fn build_unique_edges(&mut self) {
        self.unique_edges.clear();
        self.edge_face_map.clear();
        
        // Hash function for edge endpoints (order-independent)
        fn edge_key(p1: [f32; 3], p2: [f32; 3]) -> (u64, u64) {
            // Quantize to avoid floating point issues
            let quantize = |v: f32| -> u64 { (v * 10000.0) as i64 as u64 };
            
            let k1 = (quantize(p1[0]) << 42) | (quantize(p1[1]) << 21) | quantize(p1[2]);
            let k2 = (quantize(p2[0]) << 42) | (quantize(p2[1]) << 21) | quantize(p2[2]);
            
            if k1 < k2 { (k1, k2) } else { (k2, k1) }
        }
        
        for face in &self.faces {
            let edges = [
                (face.corners[0], face.corners[1]),
                (face.corners[1], face.corners[2]),
                (face.corners[2], face.corners[3]),
                (face.corners[3], face.corners[0]),
            ];
            
            for (p1, p2) in edges {
                let key = edge_key(p1, p2);
                if self.edge_face_map.insert(key) {
                    self.unique_edges.push((p1, p2));
                }
            }
        }
    }
    
    /// Rebuild just the colors (cheaper than full rebuild)
    pub fn rebuild_colors(
        &mut self,
        mesh_state: &MeshState,
        results: &Option<SimulationResults>,
        color_mode: ColorMode,
        stress_component: usize,
        strain_component: usize,
    ) {
        self.colors.node_colors.clear();
        self.colors.element_colors.clear();
        self.colors.color_mode = color_mode;
        self.colors.stress_component = stress_component;
        self.colors.strain_component = strain_component;
        
        let mesh = &mesh_state.mesh;
        
        match color_mode {
            ColorMode::Solid => {
                // No need to cache - solid color is trivial
            }
            ColorMode::Displacement => {
                if let Some(res) = results {
                    let max_disp = res.stats.max_displacement.max(1e-10);
                    for (&node_id, _node) in &mesh.nodes {
                        if node_id * 3 + 2 < res.displacements.len() {
                            let dx = res.displacements[node_id * 3];
                            let dy = res.displacements[node_id * 3 + 1];
                            let dz = res.displacements[node_id * 3 + 2];
                            let mag = (dx * dx + dy * dy + dz * dz).sqrt();
                            let t = (mag / max_disp) as f32;
                            self.colors.node_colors.insert(node_id, value_to_color(t.clamp(0.0, 1.0)));
                        }
                    }
                }
            }
            ColorMode::VonMises => {
                if let Some(res) = results {
                    let max_vm = res.stats.max_von_mises.max(1e-10);
                    for (&el_id, &vm) in &res.von_mises {
                        let t = (vm / max_vm) as f32;
                        self.colors.element_colors.insert(el_id, value_to_color(t.clamp(0.0, 1.0)));
                    }
                }
            }
            ColorMode::Stress => {
                if let Some(res) = results {
                    let max_val = res.stresses.values()
                        .filter_map(|s| s.get(stress_component))
                        .map(|v| v.abs())
                        .fold(1e-10f64, |a, b| a.max(b));
                    for (&el_id, stress) in &res.stresses {
                        if let Some(&val) = stress.get(stress_component) {
                            let t = ((val / max_val) * 0.5 + 0.5) as f32;
                            self.colors.element_colors.insert(el_id, value_to_color(t.clamp(0.0, 1.0)));
                        }
                    }
                }
            }
            ColorMode::Strain => {
                if let Some(res) = results {
                    let max_val = res.strains.values()
                        .filter_map(|s| s.get(strain_component))
                        .map(|v| v.abs())
                        .fold(1e-10f64, |a, b| a.max(b));
                    for (&el_id, strain) in &res.strains {
                        if let Some(&val) = strain.get(strain_component) {
                            let t = ((val / max_val) * 0.5 + 0.5) as f32;
                            self.colors.element_colors.insert(el_id, value_to_color(t.clamp(0.0, 1.0)));
                        }
                    }
                }
            }
        }
    }
    
    /// Get visible face indices using frustum culling
    pub fn get_visible_faces(&self, camera: &CameraState, aspect: f32) -> Vec<usize> {
        if let Some(ref octree) = self.octree {
            let mut visible = Vec::with_capacity(self.faces.len() / 2);
            octree.collect_visible(camera, aspect, &mut visible);
            visible
        } else {
            // No octree - return all faces
            (0..self.faces.len()).collect()
        }
    }
    
    /// Get color for a face (uses cache)
    pub fn get_face_color(&self, face: &CachedFace, color_mode: ColorMode) -> egui::Color32 {
        match color_mode {
            ColorMode::Solid => egui::Color32::from_rgb(100, 149, 237),
            ColorMode::Displacement => {
                // Average node colors
                let colors: Vec<_> = face.node_ids.iter()
                    .filter_map(|id| self.colors.node_colors.get(id))
                    .copied()
                    .collect();
                if colors.is_empty() {
                    egui::Color32::GRAY
                } else {
                    average_colors(&colors)
                }
            }
            ColorMode::VonMises | ColorMode::Stress | ColorMode::Strain => {
                self.colors.element_colors.get(&face.element_id)
                    .copied()
                    .unwrap_or(egui::Color32::GRAY)
            }
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

fn average_colors(colors: &[egui::Color32]) -> egui::Color32 {
    if colors.is_empty() {
        return egui::Color32::GRAY;
    }
    let mut r = 0u32;
    let mut g = 0u32;
    let mut b = 0u32;
    let mut a = 0u32;
    
    for c in colors {
        r += c.r() as u32;
        g += c.g() as u32;
        b += c.b() as u32;
        a += c.a() as u32;
    }
    
    let n = colors.len() as u32;
    egui::Color32::from_rgba_unmultiplied(
        (r / n) as u8,
        (g / n) as u8,
        (b / n) as u8,
        (a / n) as u8
    )
}

/// Darken a color by a factor (0.0 = black, 1.0 = unchanged)
pub fn darken_color(color: egui::Color32, factor: f32) -> egui::Color32 {
    egui::Color32::from_rgba_unmultiplied(
        (color.r() as f32 * factor) as u8,
        (color.g() as f32 * factor) as u8,
        (color.b() as f32 * factor) as u8,
        color.a()
    )
}
