//! Application state management for RustFEA GUI

use rust_fea::elements::Material;
use rust_fea::mesh::MeshAssembly;
use rust_fea::mesh::NodeGroup;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

/// Test mode configuration passed to the app
#[derive(Clone, Default)]
pub struct TestMode {
    pub enabled: bool,
    pub screenshot_path: Option<String>,
}

/// User settings that persist across sessions
#[derive(Clone, Serialize, Deserialize)]
pub struct UserSettings {
    /// Display settings
    pub display: DisplaySettingsSerializable,
    /// Default solver type
    pub default_solver: String,
    /// Recent files (paths as strings for serialization)
    pub recent_files: Vec<String>,
    /// Max recent files to track
    pub max_recent_files: usize,
    /// Auto-fit camera when loading mesh
    pub auto_fit_on_load: bool,
    /// Show welcome screen on startup
    pub show_welcome: bool,
    /// Animation speed multiplier
    pub default_animation_speed: f32,
    /// Default displacement scale factor
    pub default_displacement_scale: f32,
    /// Window state (if we want to restore position/size)
    pub window_maximized: bool,
}

impl Default for UserSettings {
    fn default() -> Self {
        Self {
            display: DisplaySettingsSerializable::default(),
            default_solver: "Direct".to_string(),
            recent_files: Vec::new(),
            max_recent_files: 10,
            auto_fit_on_load: true,
            show_welcome: true,
            default_animation_speed: 1.0,
            default_displacement_scale: 1.0,
            window_maximized: false,
        }
    }
}

impl UserSettings {
    /// Load settings from disk (native only)
    #[cfg(feature = "native")]
    pub fn load() -> Self {
        if let Some(config_dir) = dirs::config_dir() {
            let settings_path = config_dir.join("RustFEA").join("settings.json");
            if settings_path.exists() {
                if let Ok(contents) = std::fs::read_to_string(&settings_path) {
                    if let Ok(settings) = serde_json::from_str(&contents) {
                        return settings;
                    }
                }
            }
        }
        Self::default()
    }

    #[cfg(not(feature = "native"))]
    pub fn load() -> Self {
        // In WASM, try to load from localStorage
        #[cfg(target_arch = "wasm32")]
        {
            if let Some(window) = web_sys::window() {
                if let Ok(Some(storage)) = window.local_storage() {
                    if let Ok(Some(data)) = storage.get_item("rustfea_settings") {
                        if let Ok(settings) = serde_json::from_str(&data) {
                            return settings;
                        }
                    }
                }
            }
        }
        Self::default()
    }

    /// Save settings to disk (native only)
    #[cfg(feature = "native")]
    pub fn save(&self) {
        if let Some(config_dir) = dirs::config_dir() {
            let app_dir = config_dir.join("RustFEA");
            let _ = std::fs::create_dir_all(&app_dir);
            let settings_path = app_dir.join("settings.json");
            if let Ok(json) = serde_json::to_string_pretty(self) {
                let _ = std::fs::write(&settings_path, json);
            }
        }
    }

    #[cfg(not(feature = "native"))]
    pub fn save(&self) {
        // In WASM, save to localStorage
        #[cfg(target_arch = "wasm32")]
        {
            if let Some(window) = web_sys::window() {
                if let Ok(Some(storage)) = window.local_storage() {
                    if let Ok(json) = serde_json::to_string(self) {
                        let _ = storage.set_item("rustfea_settings", &json);
                    }
                }
            }
        }
    }

    /// Add a file to recent files list
    pub fn add_recent_file(&mut self, path: &std::path::Path) {
        let path_str = path.to_string_lossy().to_string();
        // Remove if already exists
        self.recent_files.retain(|p| p != &path_str);
        // Add to front
        self.recent_files.insert(0, path_str);
        // Trim to max
        self.recent_files.truncate(self.max_recent_files);
    }
}

/// Serializable version of DisplaySettings for persistence
#[derive(Clone, Serialize, Deserialize)]
pub struct DisplaySettingsSerializable {
    pub show_grid: bool,
    pub grid_spacing: f32,
    pub grid_size: i32,
    pub background_color: [u8; 3],
    pub wireframe_color: [u8; 3],
    pub face_color: [u8; 3],
    pub show_axis: bool,
}

impl Default for DisplaySettingsSerializable {
    fn default() -> Self {
        Self {
            show_grid: true,
            grid_spacing: 1.0,
            grid_size: 10,
            background_color: [30, 30, 35],
            wireframe_color: [40, 40, 40],
            face_color: [100, 149, 237],
            show_axis: true,
        }
    }
}

impl From<&DisplaySettings> for DisplaySettingsSerializable {
    fn from(d: &DisplaySettings) -> Self {
        Self {
            show_grid: d.show_grid,
            grid_spacing: d.grid_spacing,
            grid_size: d.grid_size,
            background_color: d.background_color,
            wireframe_color: d.wireframe_color,
            face_color: d.face_color,
            show_axis: d.show_axis,
        }
    }
}

impl DisplaySettingsSerializable {
    pub fn apply_to(&self, d: &mut DisplaySettings) {
        d.show_grid = self.show_grid;
        d.grid_spacing = self.grid_spacing;
        d.grid_size = self.grid_size;
        d.background_color = self.background_color;
        d.wireframe_color = self.wireframe_color;
        d.face_color = self.face_color;
        d.show_axis = self.show_axis;
    }
}

/// Undo/Redo action types
#[derive(Clone)]
pub enum UndoAction {
    /// Boundary condition added
    AddBoundaryCondition(BoundaryConditionConfig),
    /// Boundary condition removed (index, config)
    RemoveBoundaryCondition(usize, BoundaryConditionConfig),
    /// Boundary condition modified (index, old config)
    ModifyBoundaryCondition(usize, BoundaryConditionConfig),
    /// Material added
    AddMaterial(MaterialConfig),
    /// Material removed (index, config)
    RemoveMaterial(usize, MaterialConfig),
    /// Material modified (index, old config)
    ModifyMaterial(usize, MaterialConfig),
    /// Node group created (mesh index, group name)
    CreateNodeGroup(usize, String),
    /// Node group deleted (mesh index, group name, node group)
    DeleteNodeGroup(usize, String, NodeGroup),
    /// Mesh transform applied (mesh index, old transform values)
    MeshTransform(usize, [f64; 3], [f64; 3], [f64; 3]), // translation, scale, rotation
}

/// Undo/Redo stack for tracking changes
#[derive(Default)]
pub struct UndoStack {
    /// Actions that can be undone
    undo_stack: Vec<UndoAction>,
    /// Actions that can be redone
    redo_stack: Vec<UndoAction>,
    /// Maximum stack size
    max_size: usize,
}

impl UndoStack {
    pub fn new() -> Self {
        Self {
            undo_stack: Vec::new(),
            redo_stack: Vec::new(),
            max_size: 50,
        }
    }

    /// Push an action onto the undo stack
    pub fn push(&mut self, action: UndoAction) {
        self.undo_stack.push(action);
        self.redo_stack.clear(); // Clear redo when new action is performed

        // Trim to max size
        if self.undo_stack.len() > self.max_size {
            self.undo_stack.remove(0);
        }
    }

    /// Pop an action from the undo stack and push its inverse to redo
    pub fn pop_undo(&mut self) -> Option<UndoAction> {
        self.undo_stack.pop()
    }

    /// Push an action to redo stack (called after undoing)
    pub fn push_redo(&mut self, action: UndoAction) {
        self.redo_stack.push(action);
    }

    /// Pop an action from the redo stack
    pub fn pop_redo(&mut self) -> Option<UndoAction> {
        self.redo_stack.pop()
    }

    /// Check if undo is available
    pub fn can_undo(&self) -> bool {
        !self.undo_stack.is_empty()
    }

    /// Check if redo is available
    pub fn can_redo(&self) -> bool {
        !self.redo_stack.is_empty()
    }

    /// Get the name of the next undo action
    pub fn undo_description(&self) -> Option<&'static str> {
        self.undo_stack.last().map(|a| match a {
            UndoAction::AddBoundaryCondition(_) => "Add Boundary Condition",
            UndoAction::RemoveBoundaryCondition(_, _) => "Remove Boundary Condition",
            UndoAction::ModifyBoundaryCondition(_, _) => "Modify Boundary Condition",
            UndoAction::AddMaterial(_) => "Add Material",
            UndoAction::RemoveMaterial(_, _) => "Remove Material",
            UndoAction::ModifyMaterial(_, _) => "Modify Material",
            UndoAction::CreateNodeGroup(_, _) => "Create Node Group",
            UndoAction::DeleteNodeGroup(_, _, _) => "Delete Node Group",
            UndoAction::MeshTransform(_, _, _, _) => "Transform Mesh",
        })
    }

    /// Get the name of the next redo action
    pub fn redo_description(&self) -> Option<&'static str> {
        self.redo_stack.last().map(|a| match a {
            UndoAction::AddBoundaryCondition(_) => "Add Boundary Condition",
            UndoAction::RemoveBoundaryCondition(_, _) => "Remove Boundary Condition",
            UndoAction::ModifyBoundaryCondition(_, _) => "Modify Boundary Condition",
            UndoAction::AddMaterial(_) => "Add Material",
            UndoAction::RemoveMaterial(_, _) => "Remove Material",
            UndoAction::ModifyMaterial(_, _) => "Modify Material",
            UndoAction::CreateNodeGroup(_, _) => "Create Node Group",
            UndoAction::DeleteNodeGroup(_, _, _) => "Delete Node Group",
            UndoAction::MeshTransform(_, _, _, _) => "Transform Mesh",
        })
    }

    /// Clear both stacks
    pub fn clear(&mut self) {
        self.undo_stack.clear();
        self.redo_stack.clear();
    }
}

/// Statistics overlay configuration
#[derive(Clone)]
pub struct StatsOverlay {
    /// Show the stats overlay
    pub visible: bool,
    /// Show mesh statistics
    pub show_mesh_stats: bool,
    /// Show performance stats (FPS, render time)
    pub show_performance: bool,
    /// Show camera info
    pub show_camera_info: bool,
    /// Show result statistics when available
    pub show_result_stats: bool,
    /// Position (0 = top-left, 1 = top-right, 2 = bottom-left, 3 = bottom-right)
    pub position: u8,
}

impl Default for StatsOverlay {
    fn default() -> Self {
        Self {
            visible: false,
            show_mesh_stats: true,
            show_performance: true,
            show_camera_info: false,
            show_result_stats: true,
            position: 0, // top-left
        }
    }
}

/// Application-wide state
#[derive(Default)]
pub struct AppState {
    /// Current project path
    pub project_path: Option<PathBuf>,

    /// Loaded meshes
    pub meshes: Vec<MeshState>,

    /// Current mesh index
    pub current_mesh_idx: Option<usize>,

    /// Simulation setup
    pub simulation_config: SimulationConfig,

    /// Simulation results
    pub results: Option<SimulationResults>,

    /// UI state
    pub ui_state: UiState,

    /// Is simulation currently running
    pub is_running: bool,

    /// Simulation progress (0.0 - 1.0)
    pub progress: f32,

    /// Status message
    pub status_message: String,

    /// Detailed solve progress tracking
    pub solve_progress: SolveProgress,

    /// Undo/Redo stack
    pub undo_stack: UndoStack,

    /// User settings (persisted)
    pub user_settings: UserSettings,

    /// Statistics overlay
    pub stats_overlay: StatsOverlay,

    /// Frame timing for performance stats
    pub frame_times: Vec<f32>,
}

impl AppState {
    pub fn new() -> Self {
        // Load user settings from disk
        let user_settings = UserSettings::load();

        let mut state = Self {
            status_message: "Ready".to_string(),
            user_settings,
            undo_stack: UndoStack::new(),
            stats_overlay: StatsOverlay::default(),
            frame_times: Vec::with_capacity(120), // Store last ~2 seconds at 60fps
            ..Default::default()
        };

        // Apply saved display settings
        state
            .user_settings
            .display
            .apply_to(&mut state.ui_state.display_settings);
        state.ui_state.displacement_scale = state.user_settings.default_displacement_scale;
        state.ui_state.animation.speed = state.user_settings.default_animation_speed;

        state
    }

    /// Save current settings to disk
    pub fn save_settings(&mut self) {
        // Update serializable settings from current state
        self.user_settings.display =
            DisplaySettingsSerializable::from(&self.ui_state.display_settings);
        self.user_settings.default_displacement_scale = self.ui_state.displacement_scale;
        self.user_settings.default_animation_speed = self.ui_state.animation.speed;
        self.user_settings.save();
    }

    /// Record frame time for performance stats
    pub fn record_frame_time(&mut self, dt: f32) {
        self.frame_times.push(dt);
        if self.frame_times.len() > 120 {
            self.frame_times.remove(0);
        }
    }

    /// Get average FPS from recent frames
    pub fn average_fps(&self) -> f32 {
        if self.frame_times.is_empty() {
            return 0.0;
        }
        let avg_dt = self.frame_times.iter().sum::<f32>() / self.frame_times.len() as f32;
        if avg_dt > 0.0 {
            1.0 / avg_dt
        } else {
            0.0
        }
    }

    pub fn current_mesh(&self) -> Option<&MeshState> {
        self.current_mesh_idx.and_then(|idx| self.meshes.get(idx))
    }

    pub fn current_mesh_mut(&mut self) -> Option<&mut MeshState> {
        self.current_mesh_idx
            .and_then(|idx| self.meshes.get_mut(idx))
    }

    pub fn add_mesh(&mut self, mesh: MeshAssembly, name: String, path: Option<PathBuf>) {
        let mesh_state = MeshState::from_mesh(mesh, name, path);
        self.meshes.push(mesh_state);
        self.current_mesh_idx = Some(self.meshes.len() - 1);
    }
}

/// State for a single mesh
pub struct MeshState {
    /// The actual mesh data
    pub mesh: MeshAssembly,

    /// Display name
    pub name: String,

    /// Source file path
    pub path: Option<PathBuf>,

    /// Computed bounding box
    pub bounds: BoundingBox,

    /// GPU-ready vertex data
    pub render_data: Option<MeshRenderData>,

    /// Selected node groups
    pub selected_node_groups: Vec<String>,

    /// Selected element groups
    pub selected_element_groups: Vec<String>,
}

impl MeshState {
    pub fn from_mesh(mesh: MeshAssembly, name: String, path: Option<PathBuf>) -> Self {
        let bounds = compute_bounding_box(&mesh);
        Self {
            mesh,
            name,
            path,
            bounds,
            render_data: None,
            selected_node_groups: Vec::new(),
            selected_element_groups: Vec::new(),
        }
    }
}

/// Bounding box for mesh
#[derive(Clone, Copy, Debug, Default)]
pub struct BoundingBox {
    pub min: [f32; 3],
    pub max: [f32; 3],
}

impl BoundingBox {
    pub fn center(&self) -> [f32; 3] {
        [
            (self.min[0] + self.max[0]) / 2.0,
            (self.min[1] + self.max[1]) / 2.0,
            (self.min[2] + self.max[2]) / 2.0,
        ]
    }

    pub fn diagonal(&self) -> f32 {
        let dx = self.max[0] - self.min[0];
        let dy = self.max[1] - self.min[1];
        let dz = self.max[2] - self.min[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

pub fn compute_bounding_box_from_mesh(mesh: &MeshAssembly) -> BoundingBox {
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];

    for node in mesh.nodes.values() {
        for i in 0..3 {
            let coord = node.coordinates[i] as f32;
            min[i] = min[i].min(coord);
            max[i] = max[i].max(coord);
        }
    }

    // Handle empty mesh
    if min[0].is_infinite() {
        return BoundingBox::default();
    }

    BoundingBox { min, max }
}

fn compute_bounding_box(mesh: &MeshAssembly) -> BoundingBox {
    compute_bounding_box_from_mesh(mesh)
}

/// GPU-ready mesh data
pub struct MeshRenderData {
    pub vertices: Vec<Vertex>,
    pub indices: Vec<u32>,
    pub wireframe_indices: Vec<u32>,
}

/// Vertex for rendering
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub color: [f32; 4],
}

impl Vertex {
    pub fn new(position: [f32; 3], normal: [f32; 3], color: [f32; 4]) -> Self {
        Self {
            position,
            normal,
            color,
        }
    }
}

/// Simulation configuration
#[derive(Clone)]
pub struct SimulationConfig {
    /// Solver type
    pub solver: SolverType,

    /// Degrees of freedom per node (default: 3 for 3D solid mechanics)
    pub dofs: usize,

    /// Materials defined
    pub materials: Vec<MaterialConfig>,

    /// Boundary conditions
    pub boundary_conditions: Vec<BoundaryConditionConfig>,

    /// Output settings
    pub output: OutputConfig,

    /// Explicit solver settings
    pub explicit_settings: ExplicitSettings,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            solver: SolverType::Direct,
            dofs: 3, // 3 DOFs per node for 3D solid mechanics
            materials: vec![MaterialConfig::default()],
            boundary_conditions: Vec::new(),
            output: OutputConfig::default(),
            explicit_settings: ExplicitSettings::default(),
        }
    }
}

impl SimulationConfig {
    pub fn new() -> Self {
        Self::default()
    }
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum SolverType {
    #[default]
    Direct,
    Explicit,
}

/// Material configuration
#[derive(Clone)]
pub struct MaterialConfig {
    pub id: usize,
    pub name: String,
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub density: f64,
}

impl Default for MaterialConfig {
    fn default() -> Self {
        Self {
            id: 1,
            name: "Steel".to_string(),
            youngs_modulus: 200e9,
            poissons_ratio: 0.3,
            density: 7850.0,
        }
    }
}

impl MaterialConfig {
    pub fn to_material(&self) -> Material {
        Material::new(self.youngs_modulus, self.poissons_ratio, self.density)
    }
}

/// Boundary condition configuration
#[derive(Clone)]
pub enum BoundaryConditionConfig {
    Fixed(FixedBcConfig),
    Load(LoadBcConfig),
    Torque(TorqueBcConfig),
    Contact(ContactBcConfig),
    Pressure(PressureBcConfig),
    Traction(TractionBcConfig),
    BodyForce(BodyForceBcConfig),
}

#[derive(Clone)]
pub struct FixedBcConfig {
    pub name: String,
    pub node_group: String,
    pub constrain_x: Option<f64>,
    pub constrain_y: Option<f64>,
    pub constrain_z: Option<f64>,
}

impl Default for FixedBcConfig {
    fn default() -> Self {
        Self {
            name: "Fixed".to_string(),
            node_group: String::new(),
            constrain_x: Some(0.0),
            constrain_y: Some(0.0),
            constrain_z: Some(0.0),
        }
    }
}

#[derive(Clone)]
pub struct LoadBcConfig {
    pub name: String,
    pub node_group: String,
    pub force_x: f64,
    pub force_y: f64,
    pub force_z: f64,
}

impl Default for LoadBcConfig {
    fn default() -> Self {
        Self {
            name: "Load".to_string(),
            node_group: String::new(),
            force_x: 0.0,
            force_y: 0.0,
            force_z: 0.0,
        }
    }
}

#[derive(Clone)]
pub struct TorqueBcConfig {
    pub name: String,
    pub node_group: String,
    pub axis_point: [f64; 3],
    pub axis_direction: [f64; 3],
    pub magnitude: f64,
}

impl Default for TorqueBcConfig {
    fn default() -> Self {
        Self {
            name: "Torque".to_string(),
            node_group: String::new(),
            axis_point: [0.0, 0.0, 0.0],
            axis_direction: [0.0, 0.0, 1.0],
            magnitude: 1000.0,
        }
    }
}

#[derive(Clone)]
pub struct ContactBcConfig {
    pub name: String,
    pub primary_surface: String,
    pub secondary_surface: String,
}

impl Default for ContactBcConfig {
    fn default() -> Self {
        Self {
            name: "Contact".to_string(),
            primary_surface: String::new(),
            secondary_surface: String::new(),
        }
    }
}

/// Pressure BC - applied normal to surface elements
#[derive(Clone)]
pub struct PressureBcConfig {
    pub name: String,
    /// Surface element group (elements defining the surface)
    pub element_group: String,
    /// Pressure magnitude (positive = compression into surface)
    pub pressure: f64,
}

impl Default for PressureBcConfig {
    fn default() -> Self {
        Self {
            name: "Pressure".to_string(),
            element_group: String::new(),
            pressure: 1e6, // 1 MPa default
        }
    }
}

/// Traction BC - arbitrary surface force per unit area
#[derive(Clone)]
pub struct TractionBcConfig {
    pub name: String,
    /// Surface element group
    pub element_group: String,
    /// Traction type
    pub traction_type: TractionTypeConfig,
}

#[derive(Clone)]
pub enum TractionTypeConfig {
    /// Uniform traction in global coordinates [N/m²]
    Uniform { fx: f64, fy: f64, fz: f64 },
    /// Normal traction (positive = tension outward)
    Normal { magnitude: f64 },
    /// Tangential (shear) along surface
    Shear { magnitude: f64 },
}

impl Default for TractionTypeConfig {
    fn default() -> Self {
        TractionTypeConfig::Uniform {
            fx: 0.0,
            fy: 0.0,
            fz: -1e6,
        }
    }
}

impl Default for TractionBcConfig {
    fn default() -> Self {
        Self {
            name: "Traction".to_string(),
            element_group: String::new(),
            traction_type: TractionTypeConfig::Uniform {
                fx: 0.0,
                fy: 0.0,
                fz: -1e6,
            },
        }
    }
}

/// Body Force BC - distributed force per unit volume (gravity, centrifugal)
#[derive(Clone)]
pub struct BodyForceBcConfig {
    pub name: String,
    /// Element group to apply to (empty = all elements)
    pub element_group: String,
    /// Body force type
    pub force_type: BodyForceTypeConfig,
}

#[derive(Clone)]
pub enum BodyForceTypeConfig {
    /// Gravity acceleration vector [m/s²]
    Gravity { gx: f64, gy: f64, gz: f64 },
    /// Centrifugal force about an axis
    Centrifugal {
        axis_point: [f64; 3],
        axis_direction: [f64; 3],
        angular_velocity: f64, // rad/s
    },
    /// Uniform body force per unit volume [N/m³]
    Uniform { fx: f64, fy: f64, fz: f64 },
}

impl Default for BodyForceTypeConfig {
    fn default() -> Self {
        // Default to Earth gravity in -Y direction
        BodyForceTypeConfig::Gravity {
            gx: 0.0,
            gy: -9.81,
            gz: 0.0,
        }
    }
}

impl Default for BodyForceBcConfig {
    fn default() -> Self {
        Self {
            name: "Body Force".to_string(),
            element_group: String::new(),
            force_type: BodyForceTypeConfig::default(),
        }
    }
}

/// Output configuration
#[derive(Clone)]
pub struct OutputConfig {
    pub vtk_path: Option<PathBuf>,
    pub save_stiffness_matrix: bool,
    pub stiffness_matrix_path: Option<PathBuf>,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            vtk_path: None,
            save_stiffness_matrix: false,
            stiffness_matrix_path: None,
        }
    }
}

/// Explicit solver settings
#[derive(Clone)]
pub struct ExplicitSettings {
    pub time_steps: usize,
    pub time_step_override: Option<f64>,
    pub vtk_save_steps: usize,
    pub state_save_steps: usize,
}

impl Default for ExplicitSettings {
    fn default() -> Self {
        Self {
            time_steps: 1000,
            time_step_override: None,
            vtk_save_steps: 100,
            state_save_steps: 100,
        }
    }
}

/// Simulation results
pub struct SimulationResults {
    /// Displacement field (per node, 3 components)
    pub displacements: Vec<f64>,

    /// Stress field (per element) - averaged for backwards compatibility
    pub stresses: HashMap<usize, Vec<f64>>,

    /// Strain field (per element) - averaged for backwards compatibility
    pub strains: HashMap<usize, Vec<f64>>,

    /// Von Mises stress (per element) - averaged for backwards compatibility
    pub von_mises: HashMap<usize, f64>,

    /// Nodal Von Mises stress (per node) - for smooth interpolation
    pub nodal_von_mises: Vec<f64>,

    /// Nodal stress components (per node, 6 components: xx, yy, zz, xy, yz, xz)
    pub nodal_stress: Vec<[f64; 6]>,

    /// Nodal strain components (per node, 6 components: xx, yy, zz, xy, yz, xz)
    pub nodal_strain: Vec<[f64; 6]>,

    /// Result statistics
    pub stats: ResultStats,

    /// Time steps (for explicit solver)
    pub time_steps: Vec<TimeStepResult>,
}

#[derive(Clone, Default)]
pub struct ResultStats {
    pub max_displacement: f64,
    pub min_displacement: f64,
    pub max_von_mises: f64,
    pub solver_time_ms: u64,
    /// Detailed solve phase timing breakdown
    pub phase_timing: Option<SolvePhaseTimings>,
}

#[derive(Clone)]
pub struct TimeStepResult {
    pub time: f64,
    pub iteration: usize,
    pub max_displacement: f64,
    pub kinetic_energy: f64,
}

/// Detailed timing breakdown for solve phases
#[derive(Clone, Default)]
pub struct SolvePhaseTimings {
    /// Total wall-clock time (ms)
    pub total_ms: u64,
    /// Individual phase timings
    pub phases: Vec<SolvePhaseEntry>,
}

/// A single solve phase entry for timing/progress tracking
#[derive(Clone)]
pub struct SolvePhaseEntry {
    /// Phase name (e.g., "Assembling stiffness matrix")
    pub name: String,
    /// Phase category for grouping
    pub category: SolvePhaseCategory,
    /// Duration in milliseconds
    pub duration_ms: u64,
    /// Additional details (e.g., matrix size, NNZ count)
    pub details: Option<String>,
    /// Start time relative to solve start (ms)
    pub start_offset_ms: u64,
}

/// Categories for solve phases
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub enum SolvePhaseCategory {
    /// Initialization and setup
    Init,
    /// Matrix assembly
    Assembly,
    /// Linear solve / time integration
    Solve,
    /// Post-processing (stress computation, etc.)
    PostProcess,
    /// Other/misc
    Other,
}

impl SolvePhaseCategory {
    /// Get display color for the category (RGBA)
    pub fn color(&self) -> [u8; 4] {
        match self {
            SolvePhaseCategory::Init => [100, 149, 237, 255], // Cornflower blue
            SolvePhaseCategory::Assembly => [255, 165, 0, 255], // Orange
            SolvePhaseCategory::Solve => [50, 205, 50, 255],  // Lime green
            SolvePhaseCategory::PostProcess => [147, 112, 219, 255], // Medium purple
            SolvePhaseCategory::Other => [128, 128, 128, 255], // Gray
        }
    }

    /// Get display name
    pub fn name(&self) -> &'static str {
        match self {
            SolvePhaseCategory::Init => "Initialization",
            SolvePhaseCategory::Assembly => "Assembly",
            SolvePhaseCategory::Solve => "Solve",
            SolvePhaseCategory::PostProcess => "Post-Processing",
            SolvePhaseCategory::Other => "Other",
        }
    }
}

/// Live solve progress for displaying during solve
#[derive(Clone, Default)]
pub struct SolveProgress {
    /// Currently active phase
    pub current_phase: Option<String>,
    /// Current phase category
    pub current_category: Option<SolvePhaseCategory>,
    /// Phases completed so far
    pub completed_phases: Vec<SolvePhaseEntry>,
    /// Elapsed time (ms)
    pub elapsed_ms: u64,
    /// For explicit: current step / total steps
    pub step_progress: Option<(usize, usize)>,
    /// Estimated time remaining (ms), if available
    pub estimated_remaining_ms: Option<u64>,
}

/// UI state
#[derive(Default)]
pub struct UiState {
    /// Active panel
    pub active_panel: ActivePanel,

    /// Show wireframe
    pub show_wireframe: bool,

    /// Show mesh faces
    pub show_faces: bool,

    /// Show nodes
    pub show_nodes: bool,

    /// Show node groups
    pub show_node_groups: bool,

    /// Show boundary conditions visualization
    pub show_boundary_conditions: bool,

    /// Color mode for results
    pub color_mode: ColorMode,

    /// Selected stress component (0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz)
    pub stress_component: usize,

    /// Selected strain component (0=xx, 1=yy, 2=zz, 3=xy, 4=yz, 5=xz)
    pub strain_component: usize,

    /// Result scale factor (for displacement visualization)
    pub displacement_scale: f32,

    /// Camera state
    pub camera: CameraState,

    /// BC editor state
    pub bc_editor: BcEditorState,

    /// Material editor state
    pub material_editor: MaterialEditorState,

    /// Primitive creation dialog
    pub primitive_dialog_open: bool,
    pub primitive_config: PrimitiveConfig,

    /// Node group creator
    pub node_group_creator: NodeGroupCreator,

    /// Rename dialog
    pub rename_dialog_open: bool,
    pub rename_buffer: String,

    /// Time step playback
    pub current_time_step: usize,
    pub playback_active: bool,
    pub playback_speed: f32, // steps per second

    /// Example configuration dialog
    pub example_dialog_open: bool,
    pub example_config: crate::examples::ExampleConfig,

    /// Clipping plane for section view
    pub clipping_plane: ClippingPlane,

    /// Recent files list
    pub recent_files: RecentFiles,

    /// Screenshot settings
    pub screenshot_settings: ScreenshotSettings,

    /// Screenshot dialog open
    pub screenshot_dialog_open: bool,

    /// Request to export 2D section view as PNG
    pub section_export_requested: bool,

    /// Animation playback state  
    pub animation: AnimationState,

    /// About dialog open
    pub about_dialog_open: bool,

    /// Preferences dialog open
    pub preferences_dialog_open: bool,

    /// Show keyboard shortcuts help
    pub show_shortcuts_help: bool,

    /// Display settings (colors, grid, etc.)
    pub display_settings: DisplaySettings,

    /// Measurement tool state
    pub measurement_tool: MeasurementTool,

    /// Mesh editing state (transform, selection, context menus)
    pub mesh_edit: MeshEditState,

    /// Statistics overlay
    pub stats_overlay: StatsOverlay,

    /// Simulation progress panel (floating window)
    pub sim_progress_panel_open: bool,
}

impl UiState {
    pub fn new() -> Self {
        Self {
            show_faces: true,
            show_wireframe: true,
            show_nodes: false,
            show_node_groups: false,
            show_boundary_conditions: true,
            displacement_scale: 1.0,
            camera: CameraState::new(),
            primitive_config: PrimitiveConfig::default(),
            node_group_creator: NodeGroupCreator::default(),
            playback_speed: 10.0, // 10 steps per second
            example_config: crate::examples::ExampleConfig::default(),
            ..Default::default()
        }
    }
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum ActivePanel {
    #[default]
    Mesh,
    Setup,
    Run,
    Results,
}

#[derive(Clone, Copy, PartialEq, Default, Debug)]
pub enum ColorMode {
    #[default]
    Solid,
    Displacement,
    VonMises,
    Stress,
    Strain,
}

/// Camera state for 3D view
#[derive(Clone)]
pub struct CameraState {
    pub yaw: f32,
    pub pitch: f32,
    pub distance: f32,
    pub target: [f32; 3],
    pub fov: f32,
    /// Use orthographic projection instead of perspective
    pub orthographic: bool,
}

impl Default for CameraState {
    fn default() -> Self {
        Self::new()
    }
}

impl CameraState {
    pub fn new() -> Self {
        Self {
            yaw: 45.0_f32.to_radians(),
            pitch: 30.0_f32.to_radians(),
            distance: 5.0,
            target: [0.0, 0.0, 0.0],
            fov: 45.0_f32.to_radians(),
            orthographic: false,
        }
    }

    pub fn eye_position(&self) -> [f32; 3] {
        let x = self.target[0] + self.distance * self.pitch.cos() * self.yaw.sin();
        let y = self.target[1] + self.distance * self.pitch.sin();
        let z = self.target[2] + self.distance * self.pitch.cos() * self.yaw.cos();
        [x, y, z]
    }

    pub fn fit_to_bounds(&mut self, bounds: &BoundingBox) {
        self.target = bounds.center();
        self.distance = bounds.diagonal() * 1.5;
    }

    /// Set camera to top-down view (looking along -Y axis)
    pub fn set_top_view(&mut self) {
        self.yaw = 0.0;
        self.pitch = 89.9_f32.to_radians(); // Nearly 90 degrees to avoid gimbal lock issues
    }

    /// Set camera to front view (looking along -Z axis)
    pub fn set_front_view(&mut self) {
        self.yaw = 0.0;
        self.pitch = 0.0;
    }

    /// Set camera to side view (looking along -X axis)
    pub fn set_side_view(&mut self) {
        self.yaw = 90.0_f32.to_radians();
        self.pitch = 0.0;
    }

    /// Set camera to isometric view
    pub fn set_iso_view(&mut self) {
        self.yaw = 45.0_f32.to_radians();
        self.pitch = 30.0_f32.to_radians();
    }
}

/// BC editor state
#[derive(Default)]
pub struct BcEditorState {
    pub selected_bc_idx: Option<usize>,
    pub editing_bc: Option<BoundaryConditionConfig>,
    pub new_bc_type: NewBcType,
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum NewBcType {
    #[default]
    Fixed,
    Load,
    Torque,
    Contact,
    Pressure,
    Traction,
    BodyForce,
}

/// Material editor state
#[derive(Default)]
pub struct MaterialEditorState {
    pub selected_material_idx: Option<usize>,
    pub editing_material: Option<MaterialConfig>,
}

/// Primitive type for mesh generation
#[derive(Clone, Copy, PartialEq, Default)]
pub enum PrimitiveType {
    #[default]
    Block,
    Cylinder,
}

/// Configuration for primitive mesh generation
#[derive(Clone)]
pub struct PrimitiveConfig {
    pub primitive_type: PrimitiveType,
    pub name: String,
    // Block parameters
    pub block_size: [f64; 3],        // width, height, depth
    pub block_divisions: [usize; 3], // divisions along each axis
    pub block_origin: [f64; 3],      // origin position
    // Cylinder parameters
    pub cyl_radius: f64,
    pub cyl_height: f64,
    pub cyl_radial_divisions: usize,
    pub cyl_height_divisions: usize,
    pub cyl_axis: [f64; 3],   // axis direction
    pub cyl_origin: [f64; 3], // base center position
}

impl Default for PrimitiveConfig {
    fn default() -> Self {
        Self {
            primitive_type: PrimitiveType::Block,
            name: "Primitive".to_string(),
            block_size: [1.0, 1.0, 1.0],
            block_divisions: [2, 2, 2],
            block_origin: [0.0, 0.0, 0.0],
            cyl_radius: 0.5,
            cyl_height: 1.0,
            cyl_radial_divisions: 8,
            cyl_height_divisions: 4,
            cyl_axis: [0.0, 1.0, 0.0],
            cyl_origin: [0.0, 0.0, 0.0],
        }
    }
}

/// State for creating new node groups
#[derive(Default)]
pub struct NodeGroupCreator {
    pub active: bool,
    pub name: String,
    pub selection_mode: NodeSelectionMode,
    pub selected_nodes: Vec<usize>,
    // Box selection bounds
    pub box_min: [f64; 3],
    pub box_max: [f64; 3],
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum NodeSelectionMode {
    #[default]
    BoxSelect,
    Manual,
    Plane,
}

/// Clipping plane state for section view
#[derive(Clone)]
pub struct ClippingPlane {
    /// Whether clipping is enabled
    pub enabled: bool,
    /// Plane normal direction
    pub normal: [f32; 3],
    /// Plane position (offset from mesh center)
    pub position: f32,
    /// Which axis the plane is aligned to
    pub axis: ClipAxis,
    /// Show the plane outline
    pub show_plane: bool,
    /// Flip the clipping direction
    pub flip: bool,
    /// Show section cut surface with interpolated field values
    pub show_section_surface: bool,
    /// Show 2D cross-section view window
    pub show_2d_view: bool,
}

impl Default for ClippingPlane {
    fn default() -> Self {
        Self {
            enabled: false,
            normal: [1.0, 0.0, 0.0],
            position: 0.0,
            axis: ClipAxis::X,
            show_plane: true,
            flip: false,
            show_section_surface: true,
            show_2d_view: false,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum ClipAxis {
    #[default]
    X,
    Y,
    Z,
    Custom,
}

impl ClipAxis {
    pub fn to_normal(&self) -> [f32; 3] {
        match self {
            ClipAxis::X => [1.0, 0.0, 0.0],
            ClipAxis::Y => [0.0, 1.0, 0.0],
            ClipAxis::Z => [0.0, 0.0, 1.0],
            ClipAxis::Custom => [1.0, 0.0, 0.0], // Custom requires manual setting
        }
    }
}

/// Recent files tracking
#[derive(Clone, Default)]
pub struct RecentFiles {
    pub files: Vec<RecentFile>,
    pub max_files: usize,
}

impl RecentFiles {
    pub fn new() -> Self {
        Self {
            files: Vec::new(),
            max_files: 10,
        }
    }

    pub fn add(&mut self, path: PathBuf, name: String) {
        // Remove if already exists
        self.files.retain(|f| f.path != path);

        // Add to front
        self.files.insert(
            0,
            RecentFile {
                path,
                name,
                timestamp: std::time::SystemTime::now(),
            },
        );

        // Trim to max
        self.files.truncate(self.max_files);
    }
}

#[derive(Clone)]
pub struct RecentFile {
    pub path: PathBuf,
    pub name: String,
    pub timestamp: std::time::SystemTime,
}

/// Screenshot settings
#[derive(Clone)]
pub struct ScreenshotSettings {
    pub width: u32,
    pub height: u32,
    pub transparent_background: bool,
    pub include_ui: bool,
}

impl Default for ScreenshotSettings {
    fn default() -> Self {
        Self {
            width: 1920,
            height: 1080,
            transparent_background: false,
            include_ui: false,
        }
    }
}

/// Animation playback state
#[derive(Clone)]
pub struct AnimationState {
    pub playing: bool,
    pub loop_playback: bool,
    pub speed: f32, // Multiplier (1.0 = realtime if available, otherwise steps/sec)
    pub current_frame: usize,
    pub last_frame_time: Option<web_time::Instant>,
}

impl Default for AnimationState {
    fn default() -> Self {
        Self {
            playing: false,
            loop_playback: true,
            speed: 1.0,
            current_frame: 0,
            last_frame_time: None,
        }
    }
}

/// Display settings for viewport
#[derive(Clone)]
pub struct DisplaySettings {
    /// Show grid
    pub show_grid: bool,
    /// Grid spacing (world units)
    pub grid_spacing: f32,
    /// Grid size (number of lines)
    pub grid_size: i32,
    /// Background color
    pub background_color: [u8; 3],
    /// Wireframe color
    pub wireframe_color: [u8; 3],
    /// Default face color
    pub face_color: [u8; 3],
    /// Axis indicator visible
    pub show_axis: bool,
}

impl Default for DisplaySettings {
    fn default() -> Self {
        Self {
            show_grid: true,
            grid_spacing: 1.0,
            grid_size: 10,
            background_color: [30, 30, 35],
            wireframe_color: [40, 40, 40],
            face_color: [100, 149, 237], // Cornflower blue
            show_axis: true,
        }
    }
}

/// Measurement tool state
#[derive(Clone, Default)]
pub struct MeasurementTool {
    /// Tool is active
    pub active: bool,
    /// First point selected
    pub point1: Option<[f64; 3]>,
    /// Second point selected  
    pub point2: Option<[f64; 3]>,
    /// Measured distance
    pub distance: Option<f64>,
    /// History of measurements
    pub history: Vec<Measurement>,
}

#[derive(Clone)]
pub struct Measurement {
    pub point1: [f64; 3],
    pub point2: [f64; 3],
    pub distance: f64,
    pub label: String,
}

/// Mesh transform state for editing meshes after creation
#[derive(Clone)]
pub struct MeshTransform {
    /// Translation offset
    pub translation: [f64; 3],
    /// Scale factors
    pub scale: [f64; 3],
    /// Rotation angles (degrees) around each axis
    pub rotation: [f64; 3],
    /// Transform mode currently active
    pub mode: TransformMode,
    /// Whether transform gizmo is visible
    pub show_gizmo: bool,
}

impl Default for MeshTransform {
    fn default() -> Self {
        Self {
            translation: [0.0, 0.0, 0.0],
            scale: [1.0, 1.0, 1.0],
            rotation: [0.0, 0.0, 0.0],
            mode: TransformMode::None,
            show_gizmo: false,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum TransformMode {
    #[default]
    None,
    Translate,
    Scale,
    Rotate,
}

/// Face selection state for selecting mesh faces
#[derive(Clone, Default)]
pub struct FaceSelection {
    /// Face selection mode active
    pub active: bool,
    /// Selected face indices (element_id, face_index)
    pub selected_faces: Vec<(usize, usize)>,
    /// Hovered face (for highlighting)
    pub hovered_face: Option<(usize, usize)>,
    /// Selection mode
    pub mode: FaceSelectionMode,
    /// What to do with selected faces
    pub purpose: SelectionPurpose,
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum FaceSelectionMode {
    #[default]
    Single,
    Add,    // Shift-click to add
    Remove, // Ctrl-click to remove
    Box,    // Box select
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum SelectionPurpose {
    #[default]
    General,
    CreateNodeGroup,
    ApplyBC,
    CreateSurface,
}

/// Context menu state
#[derive(Clone, Default)]
pub struct ContextMenuState {
    /// Menu is open
    pub open: bool,
    /// Menu position
    pub position: [f32; 2],
    /// What was right-clicked
    pub target: ContextMenuTarget,
    /// Associated data
    pub target_id: Option<usize>,
}

#[derive(Clone, Copy, PartialEq, Default)]
pub enum ContextMenuTarget {
    #[default]
    Viewport,
    Mesh,
    Face,
    NodeGroup,
    ElementGroup,
    BoundaryCondition,
}

/// Mesh editing state - combines transform, selection, and editing tools
#[derive(Clone, Default)]
pub struct MeshEditState {
    /// Current mesh transform
    pub transform: MeshTransform,
    /// Face selection
    pub face_selection: FaceSelection,
    /// Context menu
    pub context_menu: ContextMenuState,
    /// Mesh density editing (for re-meshing primitives)
    pub density_edit: Option<DensityEdit>,
    /// Currently editing primitive (for live preview)
    pub editing_primitive: bool,
}

/// Density editing for re-meshing
#[derive(Clone)]
pub struct DensityEdit {
    /// Target element size
    pub target_size: f64,
    /// Divisions override
    pub divisions: [usize; 3],
    /// Refinement region (box bounds)
    pub refine_region: Option<([f64; 3], [f64; 3])>,
}
