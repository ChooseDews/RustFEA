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
        if conn.len() < 8 {
            continue; // Only support hex elements for now
        }
        
        // Get node positions (with displacement if available)
        let mut node_positions: Vec<[f32; 3]> = Vec::with_capacity(8);
        let mut node_field_values: Vec<f64> = Vec::with_capacity(8);
        
        for (local_idx, &node_id) in conn.iter().take(8).enumerate() {
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
                node_positions.push(pos);
                
                // Get field value for this node/element based on color mode
                let field_val = get_field_value(
                    results, el_id, node_id, local_idx, disp_mag,
                    color_mode, stress_component, strain_component
                );
                node_field_values.push(field_val);
            }
        }
        
        if node_positions.len() < 8 {
            continue;
        }
        
        // Compute signed distances from each node to the plane
        let mut distances: [f32; 8] = [0.0; 8];
        for (i, pos) in node_positions.iter().enumerate() {
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
                
                // Interpolate position
                let pos = [
                    node_positions[i0][0] + t * (node_positions[i1][0] - node_positions[i0][0]),
                    node_positions[i0][1] + t * (node_positions[i1][1] - node_positions[i0][1]),
                    node_positions[i0][2] + t * (node_positions[i1][2] - node_positions[i0][2]),
                ];
                
                // Interpolate natural coordinates
                let xi = HEX_NATURAL_COORDS[i0][0] + t * (HEX_NATURAL_COORDS[i1][0] - HEX_NATURAL_COORDS[i0][0]);
                let eta = HEX_NATURAL_COORDS[i0][1] + t * (HEX_NATURAL_COORDS[i1][1] - HEX_NATURAL_COORDS[i0][1]);
                let zeta = HEX_NATURAL_COORDS[i0][2] + t * (HEX_NATURAL_COORDS[i1][2] - HEX_NATURAL_COORDS[i0][2]);
                
                // Interpolate field value
                let field_value = node_field_values[i0] + (t as f64) * (node_field_values[i1] - node_field_values[i0]);
                
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

/// Generate projected section cut polygons with depth values for integration
/// with the main face rendering pipeline
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
    
    let mut result = Vec::with_capacity(cut_polygons.len());
    
    for polygon in cut_polygons {
        if polygon.vertices.len() < 3 {
            continue;
        }
        
        // Project vertices using project_close for better near-plane handling
        let mut projected: Vec<egui::Pos2> = Vec::new();
        let mut all_valid = true;
        let mut avg_field = 0.0f64;
        
        for vertex in &polygon.vertices {
            if let Some(screen_pos) = view_transform.project_close(
                vertex.position,
                rect_center,
                half_width,
                half_height,
                aspect,
            ) {
                // Check for valid screen position (allow larger bounds for section cuts)
                if screen_pos.x.is_finite() && screen_pos.y.is_finite() 
                   && rect.expand(200.0).contains(screen_pos) {
                    projected.push(screen_pos);
                    avg_field += vertex.field_value;
                } else {
                    all_valid = false;
                    break;
                }
            } else {
                all_valid = false;
                break;
            }
        }
        
        if !all_valid || projected.len() < 3 {
            continue;
        }
        
        // Check that polygon is valid (non-zero area)
        let area = compute_2d_polygon_area(&projected);
        if area.abs() < 1.0 {
            continue; // Skip degenerate polygons
        }
        
        // Compute average field value for color
        avg_field /= projected.len() as f64;
        let t = ((avg_field - min_val) / val_range).clamp(0.0, 1.0) as f32;
        let fill_color = value_to_color(t);
        
        // Compute depth from polygon center
        let depth = camera_z_fn(polygon.center);
        
        result.push(ProjectedCutPolygon {
            vertices: projected,
            depth,
            color: fill_color,
        });
    }
    
    result
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
    
    /// Update cache with new polygons
    pub fn update(&mut self, polygons: Vec<CutPolygon>, param_hash: u64) {
        self.polygons = polygons;
        self.param_hash = param_hash;
        self.valid = true;
    }
}
