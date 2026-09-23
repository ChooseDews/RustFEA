//! Main application module for RustFEA GUI

use eframe::egui;
#[cfg(not(target_arch = "wasm32"))]
use std::sync::mpsc::{self, Receiver, Sender};
#[cfg(not(target_arch = "wasm32"))]
use std::thread;

use crate::render_cache::RenderCache;
use crate::renderer::MeshRenderer;
use crate::section_cut::SectionCutCache;
use crate::state::{AppState, ActivePanel, ColorMode, SimulationResults, SolvePhaseEntry, SolvePhaseCategory, SolveProgress};
use crate::ui;

/// Messages from simulation thread to UI (native only)
#[cfg(not(target_arch = "wasm32"))]
pub enum SimMessage {
    /// Simple progress update (progress fraction, message)
    Progress(f32, String),
    /// Detailed phase update with timing
    PhaseUpdate(SolveProgress),
    /// Phase completed - add to completed list
    PhaseCompleted(SolvePhaseEntry),
    /// Simulation completed successfully
    Completed(Box<SimulationResults>),
    /// Simulation error
    Error(String),
}

/// Main application struct
pub struct FeaApp {
    /// Application state
    pub state: AppState,
    
    /// 3D mesh renderer
    pub renderer: Option<MeshRenderer>,
    
    /// Render cache for optimized viewport rendering
    pub render_cache: RenderCache,
    
    /// Section cut cache for clipping plane cross-sections
    pub section_cut_cache: SectionCutCache,
    
    #[cfg(not(target_arch = "wasm32"))]
    sim_receiver: Option<Receiver<SimMessage>>,
    #[cfg(not(target_arch = "wasm32"))]
    sim_sender: Option<Sender<bool>>,
    
    /// Pending simulation to run on next frame (WASM only - allows UI to update first)
    #[cfg(target_arch = "wasm32")]
    pending_simulation: Option<PendingSimulation>,
    
    /// Test mode for automated screenshots
    #[cfg(not(target_arch = "wasm32"))]
    test_mode: crate::state::TestMode,
    #[cfg(not(target_arch = "wasm32"))]
    test_frame_count: u32,
}

/// Pending simulation data for deferred execution (WASM)
#[cfg(target_arch = "wasm32")]
struct PendingSimulation {
    simulation: rust_fea::simulation::Simulation,
    solver_type: crate::state::SolverType,
    explicit_settings: crate::state::ExplicitSettings,
    /// Frame delay counter - wait this many frames before running (lets UI render)
    frames_to_wait: u32,
}

impl FeaApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        Self::new_with_test_mode(cc, crate::state::TestMode::default())
    }
    
    #[cfg(not(target_arch = "wasm32"))]
    pub fn new_with_test_mode(cc: &eframe::CreationContext<'_>, test_mode: crate::state::TestMode) -> Self {
        setup_custom_style(&cc.egui_ctx);
        
        let mut state = AppState::new();
        state.ui_state = crate::state::UiState::new();
        
        Self {
            state,
            renderer: None,
            render_cache: RenderCache::new(),
            section_cut_cache: SectionCutCache::new(),
            sim_receiver: None,
            sim_sender: None,
            test_mode,
            test_frame_count: 0,
        }
    }
    
    #[cfg(target_arch = "wasm32")]
    pub fn new_with_test_mode(cc: &eframe::CreationContext<'_>, _test_mode: crate::state::TestMode) -> Self {
        setup_custom_style(&cc.egui_ctx);
        
        let mut state = AppState::new();
        state.ui_state = crate::state::UiState::new();
        
        Self {
            state,
            renderer: None,
            render_cache: RenderCache::new(),
            section_cut_cache: SectionCutCache::new(),
            pending_simulation: None,
        }
    }
    
    /// Load mesh from file
    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_mesh_file(&mut self, path: std::path::PathBuf) {
        use rust_fea::io::mesh_reader::read_file_single_body;
        
        self.state.status_message = format!("Loading mesh from {:?}...", path);
        
        let path_str = path.to_string_lossy().to_string();
        match std::panic::catch_unwind(|| read_file_single_body(&path_str)) {
            Ok(mesh) => {
                let name = path.file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_else(|| "Unnamed Mesh".to_string());
                
                self.state.add_mesh(mesh, name.clone(), Some(path));
                self.state.status_message = format!("Loaded mesh: {}", name);
                
                if let Some(mesh_state) = self.state.current_mesh() {
                    let bounds = mesh_state.bounds;
                    self.state.ui_state.camera.fit_to_bounds(&bounds);
                }
                
                self.renderer = None;
            }
            Err(e) => {
                self.state.status_message = format!("Failed to load mesh: {:?}", e);
            }
        }
    }
    
    /// Load mesh from bytes (for WASM file uploads)
    /// Supports .inp (Abaqus), .json (JSON mesh), and .bin (binary) formats
    pub fn load_mesh_from_bytes(&mut self, name: String, data: &[u8]) {
        self.state.status_message = format!("Loading mesh: {}...", name);
        
        // Determine format from extension
        let extension = name.rsplit('.').next().unwrap_or("").to_lowercase();
        
        let result: Result<rust_fea::mesh::MeshAssembly, String> = match extension.as_str() {
            "inp" => {
                // Parse as Abaqus INP using internal parser
                let content = match std::str::from_utf8(data) {
                    Ok(s) => s,
                    Err(e) => {
                        self.state.status_message = format!("Invalid UTF-8 in file: {}", e);
                        return;
                    }
                };
                parse_inp_from_string(content)
            }
            "json" => {
                // Parse as JSON mesh
                match serde_json::from_slice::<rust_fea::mesh::MeshAssembly>(data) {
                    Ok(mesh) => Ok(mesh),
                    Err(e) => Err(format!("JSON parse error: {}", e)),
                }
            }
            "bin" => {
                // Parse as binary mesh (bincode)
                match bincode::deserialize::<rust_fea::mesh::MeshAssembly>(data) {
                    Ok(mesh) => Ok(mesh),
                    Err(e) => Err(format!("Binary parse error: {}", e)),
                }
            }
            _ => {
                Err(format!("Unsupported file format: .{}", extension))
            }
        };
        
        match result {
            Ok(mesh) => {
                self.state.add_mesh(mesh, name.clone(), None);
                self.state.status_message = format!("Loaded mesh: {}", name);
                
                if let Some(mesh_state) = self.state.current_mesh() {
                    let bounds = mesh_state.bounds;
                    self.state.ui_state.camera.fit_to_bounds(&bounds);
                }
                
                self.renderer = None;
            }
            Err(e) => {
                self.state.status_message = format!("Failed to load mesh: {}", e);
            }
        }
    }
    
    /// Load project from bytes (for WASM file uploads)
    pub fn load_project_from_bytes(&mut self, name: String, data: &[u8]) {
        self.state.status_message = format!("Loading project: {}...", name);
        
        let content = match std::str::from_utf8(data) {
            Ok(s) => s,
            Err(e) => {
                self.state.status_message = format!("Invalid UTF-8 in project file: {}", e);
                return;
            }
        };
        
        match crate::project_io::parse_project_toml(content) {
            Ok((project, config)) => {
                self.state.simulation_config = config;
                self.state.status_message = format!("Loaded project: {}", project.name);
                
                // Note: In WASM we can't load the referenced mesh file automatically
                // The user will need to import the mesh separately
                if project.mesh.is_some() {
                    self.state.status_message = format!(
                        "Loaded project: {}. Please import the mesh file separately.",
                        project.name
                    );
                }
            }
            Err(e) => {
                self.state.status_message = format!("Failed to load project: {}", e);
            }
        }
    }
    
    /// Load project bundle (ZIP) from bytes (for WASM file uploads)
    /// A bundle contains project.toml + mesh.json
    pub fn load_bundle_from_bytes(&mut self, name: String, data: &[u8]) {
        self.state.status_message = format!("Loading project bundle: {}...", name);
        
        // Try to extract the bundle
        let bundle = match crate::web_file_io::extract_project_bundle(data) {
            Ok(b) => b,
            Err(e) => {
                // Maybe it's just a TOML file with wrong extension
                if let Ok(content) = std::str::from_utf8(data) {
                    if content.trim_start().starts_with('[') || content.contains("name =") {
                        // Looks like TOML, try loading as project
                        self.load_project_from_bytes(name, data);
                        return;
                    }
                }
                self.state.status_message = format!("Failed to extract bundle: {}", e);
                return;
            }
        };
        
        // Load the project config first
        if let Some(toml_content) = &bundle.project_toml {
            match crate::project_io::parse_project_toml(toml_content) {
                Ok((project, config)) => {
                    self.state.simulation_config = config;
                    self.state.status_message = format!("Loaded project: {}", project.name);
                }
                Err(e) => {
                    self.state.status_message = format!("Failed to parse project.toml: {}", e);
                    return;
                }
            }
        }
        
        // Load the mesh if present
        if let Some(mesh_json) = &bundle.mesh_json {
            match serde_json::from_str::<rust_fea::mesh::MeshAssembly>(mesh_json) {
                Ok(assembly) => {
                    let mesh_name = name.trim_end_matches(".rfea").trim_end_matches(".zip").to_string();
                    self.state.add_mesh(assembly, mesh_name.clone(), None);
                    
                    // Fit camera to mesh
                    if let Some(mesh_state) = self.state.current_mesh() {
                        let bounds = mesh_state.bounds;
                        self.state.ui_state.camera.fit_to_bounds(&bounds);
                    }
                    
                    self.state.status_message = format!(
                        "Loaded project bundle with mesh: {}",
                        mesh_name
                    );
                }
                Err(e) => {
                    self.state.status_message = format!(
                        "Project loaded, but failed to parse mesh: {}",
                        e
                    );
                }
            }
        } else if bundle.project_toml.is_some() {
            self.state.status_message = "Project loaded (no mesh in bundle)".to_string();
        } else {
            self.state.status_message = "Bundle contains no project or mesh data".to_string();
        }
    }
    
    /// Create simulation from current state
    fn build_simulation(&self) -> Option<rust_fea::simulation::Simulation> {
        let mesh_state = self.state.current_mesh()?;
        let config = &self.state.simulation_config;
        
        let mut mesh = mesh_state.mesh.clone();
        mesh.single_body();
        
        let mut simulation = rust_fea::simulation::Simulation::from_mesh(mesh, config.dofs);
        
        // Add boundary conditions
        for bc_config in &config.boundary_conditions {
            match bc_config {
                crate::state::BoundaryConditionConfig::Fixed(cfg) => {
                    let node_ids: Vec<usize> = simulation.mesh.get_nodes_in_group(&cfg.node_group);
                    let values = vec![cfg.constrain_x, cfg.constrain_y, cfg.constrain_z];
                    let bc = rust_fea::bc::FixedCondition::new(node_ids, values);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::Load(cfg) => {
                    let node_ids: Vec<usize> = simulation.mesh.get_nodes_in_group(&cfg.node_group);
                    let force = nalgebra::DVector::from_vec(vec![cfg.force_x, cfg.force_y, cfg.force_z]);
                    let bc = rust_fea::bc::LoadCondition::new(node_ids, force);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::Torque(cfg) => {
                    let node_ids: Vec<usize> = simulation.mesh.get_nodes_in_group(&cfg.node_group);
                    let bc = rust_fea::bc::TorqueCondition::new_from_vec(
                        node_ids,
                        cfg.axis_point.to_vec(),
                        cfg.axis_direction.to_vec(),
                        cfg.magnitude,
                    );
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::Contact(cfg) => {
                    let mut bc = rust_fea::bc::NormalContact::new(
                        cfg.primary_surface.clone(),
                        cfg.secondary_surface.clone(),
                    );
                    let primary_nodes = simulation.mesh.get_nodes_in_group(&cfg.primary_surface);
                    let secondary_nodes = simulation.mesh.get_nodes_in_group(&cfg.secondary_surface);
                    bc.set_contact_surfaces_nodes(primary_nodes, secondary_nodes);
                    
                    let primary_elements = simulation.mesh.get_elements_in_group(&cfg.primary_surface);
                    let secondary_elements = simulation.mesh.get_elements_in_group(&cfg.secondary_surface);
                    bc.set_contact_surfaces_elements(primary_elements, secondary_elements);
                    
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::Pressure(cfg) => {
                    // Get element IDs for the surface
                    let element_ids = simulation.mesh.get_elements_in_group(&cfg.element_group);
                    let bc = rust_fea::bc::PressureCondition::new(element_ids, cfg.pressure);
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::Traction(cfg) => {
                    let element_ids = simulation.mesh.get_elements_in_group(&cfg.element_group);
                    let bc = match &cfg.traction_type {
                        crate::state::TractionTypeConfig::Uniform { fx, fy, fz } => {
                            rust_fea::bc::Traction::new(
                                element_ids,
                                nalgebra::Vector3::new(*fx, *fy, *fz)
                            )
                        }
                        crate::state::TractionTypeConfig::Normal { magnitude } => {
                            rust_fea::bc::Traction::new_normal(element_ids, *magnitude)
                        }
                        crate::state::TractionTypeConfig::Shear { magnitude } => {
                            rust_fea::bc::Traction::new_shear_xi(element_ids, *magnitude)
                        }
                    };
                    simulation.add_boundary_condition(Box::new(bc));
                }
                crate::state::BoundaryConditionConfig::BodyForce(cfg) => {
                    let bc = match &cfg.force_type {
                        crate::state::BodyForceTypeConfig::Gravity { gx, gy, gz } => {
                            rust_fea::bc::BodyForce::gravity_from_vec(*gx, *gy, *gz)
                        }
                        crate::state::BodyForceTypeConfig::Centrifugal { axis_point, axis_direction, angular_velocity } => {
                            rust_fea::bc::BodyForce::centrifugal(
                                nalgebra::Vector3::new(axis_point[0], axis_point[1], axis_point[2]),
                                nalgebra::Vector3::new(axis_direction[0], axis_direction[1], axis_direction[2]),
                                *angular_velocity
                            )
                        }
                        crate::state::BodyForceTypeConfig::Uniform { fx, fy, fz } => {
                            rust_fea::bc::BodyForce::uniform(
                                nalgebra::Vector3::new(*fx, *fy, *fz)
                            )
                        }
                    };
                    simulation.add_boundary_condition(Box::new(bc));
                }
            }
        }
        
        Some(simulation)
    }
    
    /// Start simulation
    pub fn start_simulation(&mut self) {
        if self.state.is_running {
            return;
        }
        
        let simulation = match self.build_simulation() {
            Some(sim) => sim,
            None => {
                self.state.status_message = "Cannot start: no mesh loaded".to_string();
                return;
            }
        };
        
        self.state.is_running = true;
        self.state.progress = 0.0;
        self.state.status_message = "Starting simulation...".to_string();
        self.state.solve_progress = SolveProgress::default(); // Reset solve progress
        self.state.ui_state.sim_progress_panel_open = true; // Open progress panel
        
        let solver_type = self.state.simulation_config.solver;
        let explicit_settings = self.state.simulation_config.explicit_settings.clone();
        
        #[cfg(not(target_arch = "wasm32"))]
        {
            let (tx, rx) = mpsc::channel();
            let (cmd_tx, _cmd_rx) = mpsc::channel();
            
            self.sim_receiver = Some(rx);
            self.sim_sender = Some(cmd_tx);
            
            thread::spawn(move || {
                run_simulation_threaded(simulation, solver_type, explicit_settings, tx);
            });
        }
        
        #[cfg(target_arch = "wasm32")]
        {
            // Defer simulation - wait for UI to fully render the progress panel
            self.pending_simulation = Some(PendingSimulation {
                simulation,
                solver_type,
                explicit_settings,
                frames_to_wait: 30,  // Wait ~500ms at 60fps for UI to settle
            });
        }
    }
    
    /// Stop running simulation
    pub fn stop_simulation(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(sender) = &self.sim_sender {
            let _ = sender.send(false);
        }
        self.state.is_running = false;
        // Clear solve progress when stopping
        self.state.solve_progress = SolveProgress::default();
    }
    
    /// Check for simulation updates
    fn poll_simulation(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(receiver) = &self.sim_receiver {
            while let Ok(msg) = receiver.try_recv() {
                match msg {
                    SimMessage::Progress(progress, message) => {
                        self.state.progress = progress;
                        self.state.status_message = message.clone();
                        // Update current phase name for legacy messages
                        self.state.solve_progress.current_phase = Some(message);
                    }
                    SimMessage::PhaseUpdate(progress) => {
                        self.state.solve_progress = progress;
                        if let Some(ref phase) = self.state.solve_progress.current_phase {
                            self.state.status_message = phase.clone();
                        }
                    }
                    SimMessage::PhaseCompleted(entry) => {
                        self.state.solve_progress.completed_phases.push(entry);
                    }
                    SimMessage::Completed(results) => {
                        self.state.is_running = false;
                        self.state.progress = 1.0;
                        self.state.results = Some(*results);
                        self.state.status_message = "Simulation completed!".to_string();
                        self.state.ui_state.active_panel = ActivePanel::Results;
                        // Auto-switch to displacement coloring if currently solid
                        if self.state.ui_state.color_mode == ColorMode::Solid {
                            self.state.ui_state.color_mode = ColorMode::Displacement;
                        }
                        self.renderer = None;
                    }
                    SimMessage::Error(error) => {
                        self.state.is_running = false;
                        self.state.status_message = format!("Error: {}", error);
                    }
                }
            }
        }
    }
}

impl eframe::App for FeaApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_simulation();
        
        // Run deferred simulation (WASM) - wait a few frames for UI to render first
        #[cfg(target_arch = "wasm32")]
        if let Some(mut pending) = self.pending_simulation.take() {
            if pending.frames_to_wait > 0 {
                // Still waiting - decrement counter and put it back
                pending.frames_to_wait -= 1;
                self.pending_simulation = Some(pending);
                ctx.request_repaint(); // Keep repainting during countdown
            } else {
                // Ready to run
                let result = run_simulation_sync(pending.simulation, pending.solver_type, pending.explicit_settings);
                match result {
                    Ok(results) => {
                        self.state.is_running = false;
                        self.state.progress = 1.0;
                        self.state.results = Some(results);
                        self.state.status_message = "Simulation completed!".to_string();
                        // Auto-switch to displacement coloring if currently solid
                        if self.state.ui_state.color_mode == ColorMode::Solid {
                            self.state.ui_state.color_mode = ColorMode::Displacement;
                        }
                        self.renderer = None;
                    }
                    Err(e) => {
                        self.state.is_running = false;
                        self.state.status_message = format!("Error: {}", e);
                    }
                }
            }
        }
        
        // Check for pending file uploads (WASM)
        if let Some(pending) = crate::web_file_io::take_pending_file() {
            match pending.file_type {
                crate::web_file_io::FileType::Mesh => {
                    self.load_mesh_from_bytes(pending.name, &pending.data);
                }
                crate::web_file_io::FileType::Project => {
                    self.load_project_from_bytes(pending.name, &pending.data);
                }
                crate::web_file_io::FileType::Bundle => {
                    self.load_bundle_from_bytes(pending.name, &pending.data);
                }
            }
        }
        
        // Handle test mode automation
        #[cfg(not(target_arch = "wasm32"))]
        if self.test_mode.enabled {
            self.test_frame_count += 1;
            
            // Frame 2: Load the torque shaft example with Fine mesh
            if self.test_frame_count == 2 {
                println!("[Test] Loading torque shaft example (Fine mesh)...");
                let config = crate::examples::ExampleConfig {
                    example_type: crate::examples::ExampleType::TorqueShaft,
                    resolution: crate::examples::MeshResolution::Fine,
                    ..Default::default()
                };
                let example = crate::examples::load_example_with_config(&config);
                
                let name = example.name.clone();
                self.state.add_mesh(example.mesh, name.clone(), None);
                
                // Set up the simulation config from the example
                self.state.simulation_config.boundary_conditions = example.boundary_conditions;
                self.state.simulation_config.solver = example.solver_type;
                
                if let Some(mesh_state) = self.state.current_mesh() {
                    let bounds = mesh_state.bounds;
                    self.state.ui_state.camera.fit_to_bounds(&bounds);
                }
                println!("[Test] Loaded mesh: {}", name);
            }
            
            // Frame 5: Switch to top view
            if self.test_frame_count == 5 {
                println!("[Test] Switching to top view...");
                self.state.ui_state.camera.set_top_view();
                if let Some(mesh_state) = self.state.current_mesh() {
                    let bounds = mesh_state.bounds;
                    self.state.ui_state.camera.fit_to_bounds(&bounds);
                }
                println!("[Test] Top view set");
            }
            
            // Frame 10: Take screenshot and exit
            if self.test_frame_count == 10 {
                if let Some(ref screenshot_path) = self.test_mode.screenshot_path {
                    println!("[Test] Taking screenshot to: {}", screenshot_path);
                    // Request the screenshot - egui will handle it
                    ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(Default::default()));
                }
            }
            
            // Frame 15: Exit
            if self.test_frame_count >= 15 {
                println!("[Test] Test complete, exiting...");
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
            
            // Keep requesting repaints in test mode
            ctx.request_repaint();
        }
        
        // Handle keyboard shortcuts
        ctx.input(|i| {
            // Handle screenshot callback in test mode
            #[cfg(not(target_arch = "wasm32"))]
            if self.test_mode.enabled {
                for event in &i.raw.events {
                    if let egui::Event::Screenshot { image, .. } = event {
                        if let Some(ref path) = self.test_mode.screenshot_path {
                            println!("[Test] Saving screenshot ({} x {})...", image.width(), image.height());
                            // Convert to image and save
                            let pixels: Vec<u8> = image.pixels.iter()
                                .flat_map(|c| [c.r(), c.g(), c.b(), c.a()])
                                .collect();
                            if let Some(img) = image::RgbaImage::from_raw(
                                image.width() as u32,
                                image.height() as u32,
                                pixels
                            ) {
                                if let Err(e) = img.save(path) {
                                    eprintln!("[Test] Failed to save screenshot: {}", e);
                                } else {
                                    println!("[Test] Screenshot saved to: {}", path);
                                }
                            }
                        }
                    }
                }
            }
            
            // F - Fit view
            if i.key_pressed(egui::Key::F) && !i.modifiers.any() {
                if let Some(mesh) = self.state.current_mesh() {
                    let bounds = mesh.bounds;
                    self.state.ui_state.camera.fit_to_bounds(&bounds);
                }
            }
            // W - Toggle wireframe
            if i.key_pressed(egui::Key::W) && !i.modifiers.any() {
                self.state.ui_state.show_wireframe = !self.state.ui_state.show_wireframe;
                self.renderer = None;
            }
            // S - Toggle solid faces
            if i.key_pressed(egui::Key::S) && !i.modifiers.any() {
                self.state.ui_state.show_faces = !self.state.ui_state.show_faces;
                self.renderer = None;
            }
            // N - Toggle nodes
            if i.key_pressed(egui::Key::N) && !i.modifiers.any() {
                self.state.ui_state.show_nodes = !self.state.ui_state.show_nodes;
                self.renderer = None;
            }
            // B - Toggle boundary conditions
            if i.key_pressed(egui::Key::B) && !i.modifiers.any() {
                self.state.ui_state.show_boundary_conditions = !self.state.ui_state.show_boundary_conditions;
            }
            // Space - Play/pause (when results available)
            if i.key_pressed(egui::Key::Space) {
                if let Some(results) = &self.state.results {
                    if !results.time_steps.is_empty() {
                        self.state.ui_state.playback_active = !self.state.ui_state.playback_active;
                    }
                }
            }
            // 1-4 - Switch panels
            if i.key_pressed(egui::Key::Num1) {
                self.state.ui_state.active_panel = crate::state::ActivePanel::Mesh;
            }
            if i.key_pressed(egui::Key::Num2) {
                self.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
            }
            if i.key_pressed(egui::Key::Num3) {
                self.state.ui_state.active_panel = crate::state::ActivePanel::Run;
            }
            if i.key_pressed(egui::Key::Num4) {
                self.state.ui_state.active_panel = crate::state::ActivePanel::Results;
            }
        });
        
        // Handle time step playback
        if self.state.ui_state.playback_active {
            if let Some(results) = &self.state.results {
                let num_steps = results.time_steps.len();
                if num_steps > 0 {
                    // Advance frame based on playback speed
                    self.state.ui_state.current_time_step += 1;
                    if self.state.ui_state.current_time_step >= num_steps {
                        self.state.ui_state.current_time_step = 0; // Loop
                    }
                    self.renderer = None;
                    
                    // Request repaint after delay based on playback speed
                    let delay_ms = (1000.0 / self.state.ui_state.playback_speed) as u64;
                    ctx.request_repaint_after(std::time::Duration::from_millis(delay_ms));
                } else {
                    self.state.ui_state.playback_active = false;
                }
            }
        }
        
        if self.state.is_running {
            ctx.request_repaint();
        }
        
        ui::menu_bar::show(ctx, self);
        ui::menu_bar::show_example_dialog(ctx, self);
        ui::side_panel::show(ctx, self);
        ui::status_bar::show(ctx, self);
        ui::viewport::show(ctx, self);
        
        // Show About dialog
        if self.state.ui_state.about_dialog_open {
            show_about_dialog(ctx, self);
        }
        
        // Show Screenshot dialog
        #[cfg(not(target_arch = "wasm32"))]
        if self.state.ui_state.screenshot_dialog_open {
            show_screenshot_dialog(ctx, self);
        }
        
        // Show Preferences dialog
        if self.state.ui_state.preferences_dialog_open {
            show_preferences_dialog(ctx, self);
        }
        
        // Show Simulation Progress panel (floating window)
        if self.state.ui_state.sim_progress_panel_open {
            show_simulation_progress_panel(ctx, self);
        }
        
        // Show 2D Section View window when enabled
        if self.state.ui_state.clipping_plane.enabled && self.state.ui_state.clipping_plane.show_2d_view {
            show_2d_section_view(ctx, self);
        }
    }
}

/// Show the About dialog
fn show_about_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    egui::Window::new("About RustFEA")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.vertical_centered(|ui| {
                ui.heading("RustFEA");
                ui.label("Finite Element Analysis in Rust");
                ui.add_space(8.0);
                
                ui.label(format!("Version: {}", env!("CARGO_PKG_VERSION")));
                ui.add_space(8.0);
                
                ui.horizontal(|ui| {
                    ui.label("A project by");
                    if ui.link("John Dews-Flick").clicked() {
                        let url = "https://johndews.com";
                        #[cfg(not(target_arch = "wasm32"))]
                        {
                            let _ = open::that(url);
                        }
                        #[cfg(target_arch = "wasm32")]
                        {
                            if let Some(window) = web_sys::window() {
                                let _ = window.open_with_url_and_target(url, "_blank");
                            }
                        }
                    }
                });
                
                ui.add_space(4.0);
                
                ui.horizontal(|ui| {
                    ui.label("Source code on");
                    if ui.link("GitHub").clicked() {
                        let url = "https://github.com/ChooseDews/RustFEA";
                        #[cfg(not(target_arch = "wasm32"))]
                        {
                            let _ = open::that(url);
                        }
                        #[cfg(target_arch = "wasm32")]
                        {
                            if let Some(window) = web_sys::window() {
                                let _ = window.open_with_url_and_target(url, "_blank");
                            }
                        }
                    }
                });
                
                ui.add_space(12.0);
                
                ui.label("A modular FEA library with support for:");
                ui.label("• 3D solid mechanics");
                ui.label("• Direct and explicit solvers");
                ui.label("• Contact analysis");
                ui.label("• Multiple boundary conditions");
                
                ui.add_space(16.0);
                ui.label("Built with:");
                ui.horizontal(|ui| {
                    ui.label("• egui/eframe for GUI");
                });
                ui.horizontal(|ui| {
                    ui.label("• nalgebra for linear algebra");
                });
                ui.horizontal(|ui| {
                    ui.label("• russell_sparse for solvers");
                });
                
                ui.add_space(16.0);
                
                if ui.button("Close").clicked() {
                    app.state.ui_state.about_dialog_open = false;
                }
            });
        });
}

/// Show the Screenshot dialog
#[cfg(not(target_arch = "wasm32"))]
fn show_screenshot_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    egui::Window::new("📷 Save Screenshot")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.label("Save the current viewport as an image file.");
            ui.add_space(8.0);
            
            egui::Grid::new("screenshot_options")
                .num_columns(2)
                .spacing([20.0, 8.0])
                .show(ui, |ui| {
                    ui.label("Width:");
                    ui.add(egui::DragValue::new(&mut app.state.ui_state.screenshot_settings.width)
                        .range(640..=7680)
                        .suffix(" px"));
                    ui.end_row();
                    
                    ui.label("Height:");
                    ui.add(egui::DragValue::new(&mut app.state.ui_state.screenshot_settings.height)
                        .range(480..=4320)
                        .suffix(" px"));
                    ui.end_row();
                });
            
            ui.add_space(8.0);
            
            // Preset buttons
            ui.horizontal(|ui| {
                ui.label("Presets:");
                if ui.button("1080p").clicked() {
                    app.state.ui_state.screenshot_settings.width = 1920;
                    app.state.ui_state.screenshot_settings.height = 1080;
                }
                if ui.button("4K").clicked() {
                    app.state.ui_state.screenshot_settings.width = 3840;
                    app.state.ui_state.screenshot_settings.height = 2160;
                }
                if ui.button("Square").clicked() {
                    app.state.ui_state.screenshot_settings.width = 1024;
                    app.state.ui_state.screenshot_settings.height = 1024;
                }
            });
            
            ui.add_space(16.0);
            
            ui.horizontal(|ui| {
                if ui.button("Save...").clicked() {
                    // Use file dialog to choose save location
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("PNG Image", &["png"])
                        .set_file_name("screenshot.png")
                        .save_file()
                    {
                        // Request screenshot
                        ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(Default::default()));
                        app.state.status_message = format!("Screenshot saved to: {}", path.display());
                        app.state.ui_state.screenshot_dialog_open = false;
                    }
                }
                
                if ui.button("Cancel").clicked() {
                    app.state.ui_state.screenshot_dialog_open = false;
                }
            });
        });
}

/// Show the Preferences dialog
fn show_preferences_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    egui::Window::new("Preferences")
        .collapsible(false)
        .resizable(true)
        .default_width(450.0)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("Display Settings");
                ui.add_space(8.0);
                
                // Grid settings
                ui.group(|ui| {
                    ui.label("Grid");
                    ui.checkbox(&mut app.state.ui_state.display_settings.show_grid, "Show Grid");
                    
                    ui.horizontal(|ui| {
                        ui.label("Grid Spacing:");
                        ui.add(egui::DragValue::new(&mut app.state.ui_state.display_settings.grid_spacing)
                            .range(0.1..=100.0)
                            .speed(0.1));
                    });
                    
                    ui.horizontal(|ui| {
                        ui.label("Grid Size:");
                        ui.add(egui::DragValue::new(&mut app.state.ui_state.display_settings.grid_size)
                            .range(1..=50));
                    });
                });
                
                ui.add_space(8.0);
                
                // Background color
                ui.group(|ui| {
                    ui.label("Colors");
                    
                    ui.horizontal(|ui| {
                        ui.label("Background:");
                        let color = app.state.ui_state.display_settings.background_color;
                        let mut color32 = egui::Color32::from_rgb(color[0], color[1], color[2]);
                        if ui.color_edit_button_srgba(&mut color32).changed() {
                            app.state.ui_state.display_settings.background_color = [color32.r(), color32.g(), color32.b()];
                        }
                        
                        // Presets
                        if ui.small_button("Dark").clicked() {
                            app.state.ui_state.display_settings.background_color = [30, 30, 35];
                        }
                        if ui.small_button("Light").clicked() {
                            app.state.ui_state.display_settings.background_color = [200, 200, 200];
                        }
                        if ui.small_button("Blue").clicked() {
                            app.state.ui_state.display_settings.background_color = [20, 30, 50];
                        }
                    });
                    
                    ui.horizontal(|ui| {
                        ui.label("Wireframe:");
                        let mut color32 = egui::Color32::from_rgb(
                            app.state.ui_state.display_settings.wireframe_color[0],
                            app.state.ui_state.display_settings.wireframe_color[1],
                            app.state.ui_state.display_settings.wireframe_color[2],
                        );
                        if ui.color_edit_button_srgba(&mut color32).changed() {
                            app.state.ui_state.display_settings.wireframe_color = [color32.r(), color32.g(), color32.b()];
                        }
                    });
                    
                    ui.horizontal(|ui| {
                        ui.label("Default Mesh:");
                        let mut color32 = egui::Color32::from_rgb(
                            app.state.ui_state.display_settings.face_color[0],
                            app.state.ui_state.display_settings.face_color[1],
                            app.state.ui_state.display_settings.face_color[2],
                        );
                        if ui.color_edit_button_srgba(&mut color32).changed() {
                            app.state.ui_state.display_settings.face_color = [color32.r(), color32.g(), color32.b()];
                        }
                    });
                });
                
                ui.add_space(8.0);
                
                // Viewport settings
                ui.group(|ui| {
                    ui.label("Viewport");
                    ui.checkbox(&mut app.state.ui_state.display_settings.show_axis, "Show Axis Indicator");
                    ui.checkbox(&mut app.state.user_settings.auto_fit_on_load, "Auto-fit camera when loading mesh");
                });
                
                ui.add_space(8.0);
                
                // Stats overlay settings
                ui.group(|ui| {
                    ui.label("Statistics Overlay");
                    ui.checkbox(&mut app.state.ui_state.stats_overlay.visible, "Show Stats Overlay (I)");
                    
                    ui.add_enabled_ui(app.state.ui_state.stats_overlay.visible, |ui| {
                        ui.indent("stats_opts", |ui| {
                            ui.checkbox(&mut app.state.ui_state.stats_overlay.show_mesh_stats, "Mesh Statistics");
                            ui.checkbox(&mut app.state.ui_state.stats_overlay.show_result_stats, "Result Statistics");
                            ui.checkbox(&mut app.state.ui_state.stats_overlay.show_performance, "Performance (FPS)");
                            ui.checkbox(&mut app.state.ui_state.stats_overlay.show_camera_info, "Camera Info");
                            
                            ui.horizontal(|ui| {
                                ui.label("Position:");
                                egui::ComboBox::from_id_salt("stats_pos")
                                    .selected_text(match app.state.ui_state.stats_overlay.position {
                                        0 => "Top Left",
                                        1 => "Top Right",
                                        2 => "Bottom Left",
                                        _ => "Bottom Right",
                                    })
                                    .show_ui(ui, |ui| {
                                        ui.selectable_value(&mut app.state.ui_state.stats_overlay.position, 0, "Top Left");
                                        ui.selectable_value(&mut app.state.ui_state.stats_overlay.position, 1, "Top Right");
                                        ui.selectable_value(&mut app.state.ui_state.stats_overlay.position, 2, "Bottom Left");
                                        ui.selectable_value(&mut app.state.ui_state.stats_overlay.position, 3, "Bottom Right");
                                    });
                            });
                        });
                    });
                });
                
                ui.add_space(8.0);
                
                // Animation defaults
                ui.group(|ui| {
                    ui.label("Animation Defaults");
                    ui.horizontal(|ui| {
                        ui.label("Speed:");
                        ui.add(egui::Slider::new(&mut app.state.user_settings.default_animation_speed, 0.1..=10.0)
                            .suffix("x"));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Displacement Scale:");
                        ui.add(egui::Slider::new(&mut app.state.user_settings.default_displacement_scale, 0.1..=100.0)
                            .logarithmic(true));
                    });
                });
                
                ui.add_space(16.0);
                
                ui.horizontal(|ui| {
                    if ui.button("Reset to Defaults").clicked() {
                        app.state.ui_state.display_settings = crate::state::DisplaySettings::default();
                        app.state.user_settings = crate::state::UserSettings::default();
                        app.state.ui_state.stats_overlay = crate::state::StatsOverlay::default();
                    }
                    
                    if ui.button("Save Settings").clicked() {
                        app.state.save_settings();
                        app.state.status_message = "Settings saved".to_string();
                    }
                    
                    if ui.button("Close").clicked() {
                        app.state.ui_state.preferences_dialog_open = false;
                    }
                });
            });
        });
}

/// Show the Simulation Progress floating panel
fn show_simulation_progress_panel(ctx: &egui::Context, app: &mut FeaApp) {
    use crate::state::ActivePanel;
    
    let is_running = app.state.is_running;
    let has_results = app.state.results.is_some();
    
    let title = if is_running {
        "⏳ Simulation Running..."
    } else if has_results {
        "✅ Simulation Complete"
    } else {
        "Simulation"
    };
    
    let mut open = app.state.ui_state.sim_progress_panel_open;
    
    // Get screen size for centering
    let screen_rect = ctx.screen_rect();
    let center_x = screen_rect.center().x;
    let center_y = screen_rect.center().y;
    
    egui::Window::new(title)
        .open(&mut open)
        .collapsible(true)
        .resizable(true)
        .default_width(320.0)
        .min_width(280.0)
        .pivot(egui::Align2::CENTER_CENTER)
        .default_pos([center_x, center_y])  // Centered, but movable
        .show(ctx, |ui| {
            if is_running {
                // Running state - show progress
                ui.add(
                    egui::ProgressBar::new(app.state.progress)
                        .show_percentage()
                        .animate(true)
                );
                
                ui.add_space(8.0);
                
                // Status / Elapsed time
                ui.horizontal(|ui| {
                    if app.state.solve_progress.elapsed_ms == 0 {
                        // Show "preparing" message when we haven't started tracking yet
                        ui.label(egui::RichText::new("⏳ Preparing simulation...").color(egui::Color32::from_rgb(150, 150, 200)));
                    } else {
                        ui.label("Elapsed:");
                        ui.label(format_duration(app.state.solve_progress.elapsed_ms));
                        
                        if let Some(remaining) = app.state.solve_progress.estimated_remaining_ms {
                            ui.separator();
                            ui.label("ETA:");
                            ui.label(format_duration(remaining));
                        }
                    }
                });
                
                // Current phase
                if let Some(ref phase) = app.state.solve_progress.current_phase {
                    ui.add_space(4.0);
                    ui.horizontal(|ui| {
                        if let Some(cat) = app.state.solve_progress.current_category {
                            let color = category_color(cat);
                            ui.label(egui::RichText::new(category_icon(cat)).color(color));
                        }
                        ui.label(egui::RichText::new(phase).strong());
                    });
                }
                
                // Step progress for explicit solver
                if let Some((current, total)) = app.state.solve_progress.step_progress {
                    ui.add_space(4.0);
                    let step_progress = current as f32 / total as f32;
                    ui.add(
                        egui::ProgressBar::new(step_progress)
                            .text(format!("Step {}/{}", current, total))
                    );
                }
                
                ui.add_space(8.0);
                
                // Completed phases
                if !app.state.solve_progress.completed_phases.is_empty() {
                    ui.collapsing("Completed Phases", |ui| {
                        for entry in &app.state.solve_progress.completed_phases {
                            ui.horizontal(|ui| {
                                let color = category_color(entry.category);
                                ui.label(egui::RichText::new(category_icon(entry.category)).color(color));
                                ui.label(&entry.name);
                                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                    ui.label(egui::RichText::new(format_duration(entry.duration_ms)).monospace());
                                });
                            });
                        }
                    });
                }
                
                ui.add_space(12.0);
                
                // Cancel button
                ui.horizontal(|ui| {
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("Cancel").clicked() {
                            app.stop_simulation();
                        }
                    });
                });
                
            } else if let Some(ref results) = app.state.results {
                // Completed state - show results summary
                ui.label(egui::RichText::new("Simulation completed successfully!").strong().color(egui::Color32::from_rgb(100, 200, 100)));
                
                ui.add_space(8.0);
                
                // Time stats
                egui::Grid::new("result_stats")
                    .num_columns(2)
                    .spacing([20.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Total Time:");
                        ui.label(egui::RichText::new(format_duration(results.stats.solver_time_ms)).strong());
                        ui.end_row();
                        
                        ui.label("Max Displacement:");
                        ui.label(format!("{:.4e}", results.stats.max_displacement));
                        ui.end_row();
                        
                        ui.label("Max Von Mises:");
                        ui.label(format!("{:.4e}", results.stats.max_von_mises));
                        ui.end_row();
                    });
                
                // Phase timing breakdown
                if let Some(ref timing) = results.stats.phase_timing {
                    ui.add_space(8.0);
                    egui::CollapsingHeader::new("Phase Timing Breakdown")
                        .default_open(true)
                        .show(ui, |ui| {
                            show_phase_timing_in_panel(ui, timing);
                        });
                }
                
                ui.add_space(12.0);
                
                // Action buttons
                ui.horizontal(|ui| {
                    if ui.button("View Results →").clicked() {
                        app.state.ui_state.active_panel = ActivePanel::Results;
                        app.state.ui_state.sim_progress_panel_open = false;
                    }
                    
                    if ui.button("Close").clicked() {
                        app.state.ui_state.sim_progress_panel_open = false;
                    }
                });
                
            } else {
                // Error or cancelled state
                ui.label(&app.state.status_message);
                
                ui.add_space(12.0);
                
                if ui.button("Close").clicked() {
                    app.state.ui_state.sim_progress_panel_open = false;
                }
            }
        });
    
    // Only update if the window's X button was clicked (open becomes false)
    // Don't overwrite changes from the Close buttons inside the window
    if !open {
        app.state.ui_state.sim_progress_panel_open = false;
    }
}

/// Show phase timing breakdown in the progress panel
fn show_phase_timing_in_panel(ui: &mut egui::Ui, timing: &crate::state::SolvePhaseTimings) {
    use crate::state::SolvePhaseCategory;
    
    let total_ms = timing.total_ms.max(1) as f32;
    let bar_width = 240.0_f32.min(ui.available_width() - 20.0);
    
    // Group by category
    let mut category_totals: std::collections::HashMap<SolvePhaseCategory, u64> = std::collections::HashMap::new();
    for phase in &timing.phases {
        *category_totals.entry(phase.category).or_insert(0) += phase.duration_ms;
    }
    
    // Category bar
    ui.label("Time by Category:");
    let (rect, _response) = ui.allocate_exact_size(
        egui::vec2(bar_width, 20.0),
        egui::Sense::hover()
    );
    
    let painter = ui.painter();
    let mut x_offset = rect.left();
    
    for cat in &[SolvePhaseCategory::Init, SolvePhaseCategory::Assembly, 
                 SolvePhaseCategory::Solve, SolvePhaseCategory::PostProcess] {
        if let Some(&cat_ms) = category_totals.get(cat) {
            let width = (cat_ms as f32 / total_ms) * bar_width;
            if width > 1.0 {
                let color = category_color(*cat);
                painter.rect_filled(
                    egui::Rect::from_min_size(
                        egui::pos2(x_offset, rect.top()),
                        egui::vec2(width, rect.height())
                    ),
                    2.0,
                    color
                );
                x_offset += width;
            }
        }
    }
    
    ui.add_space(4.0);
    
    // Legend - each on its own line
    for cat in &[SolvePhaseCategory::Init, SolvePhaseCategory::Assembly, 
                 SolvePhaseCategory::Solve, SolvePhaseCategory::PostProcess] {
        if let Some(&cat_ms) = category_totals.get(cat) {
            let pct = (cat_ms as f32 / total_ms) * 100.0;
            let color = category_color(*cat);
            ui.horizontal(|ui| {
                let (swatch_rect, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
                ui.painter().rect_filled(swatch_rect, 2.0, color);
                ui.label(format!("{}: {} ({:.1}%)", cat.name(), format_duration(cat_ms), pct));
            });
        }
    }
    
    ui.add_space(8.0);
    
    // Detailed phase list
    ui.collapsing("Phase Details", |ui| {
        for phase in &timing.phases {
            ui.horizontal(|ui| {
                let color = category_color(phase.category);
                ui.label(egui::RichText::new(category_icon(phase.category)).color(color));
                ui.label(&phase.name);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let pct = (phase.duration_ms as f32 / total_ms) * 100.0;
                    ui.label(format!("{:.1}%", pct));
                    ui.label(format_duration(phase.duration_ms));
                });
            });
            if let Some(ref details) = phase.details {
                ui.indent("detail", |ui| {
                    ui.label(egui::RichText::new(details).small().weak());
                });
            }
        }
    });
}

/// Format duration in ms to human-readable string
fn format_duration(ms: u64) -> String {
    if ms < 1000 {
        format!("{}ms", ms)
    } else if ms < 60000 {
        format!("{:.2}s", ms as f64 / 1000.0)
    } else {
        let mins = ms / 60000;
        let secs = (ms % 60000) / 1000;
        format!("{}m {}s", mins, secs)
    }
}

/// Show the 2D cross-section view window
fn show_2d_section_view(ctx: &egui::Context, app: &mut FeaApp) {
    use crate::section_cut;
    use crate::state::ClipAxis;
    
    let mut open = app.state.ui_state.clipping_plane.show_2d_view;
    
    // Determine axis labels based on clip plane orientation
    let (h_axis_label, v_axis_label) = match app.state.ui_state.clipping_plane.axis {
        ClipAxis::X => ("Y", "Z"),
        ClipAxis::Y => ("X", "Z"),
        ClipAxis::Z => ("X", "Y"),
        ClipAxis::Custom => ("U", "V"),
    };
    
    let title = format!("2D Section View ({}-{} Plane)", h_axis_label, v_axis_label);
    
    egui::Window::new(title)
        .open(&mut open)
        .collapsible(true)
        .resizable(true)
        .default_width(400.0)
        .default_height(400.0)
        .min_width(200.0)
        .min_height(200.0)
        .default_pos([100.0, 100.0])  // Ensure window appears on screen
        .show(ctx, |ui| {
            // Get mesh and results if available
            let mesh_state = app.state.current_mesh();
            let results = &app.state.results;
            
            if mesh_state.is_none() {
                ui.centered_and_justified(|ui| {
                    ui.label("No mesh loaded");
                });
                return;
            }
            
            let mesh_state = mesh_state.unwrap();
            let ui_state = &app.state.ui_state;
            
            // Get the clipping plane parameters
            let clip_normal = if ui_state.clipping_plane.flip {
                [-ui_state.clipping_plane.normal[0], 
                 -ui_state.clipping_plane.normal[1], 
                 -ui_state.clipping_plane.normal[2]]
            } else {
                ui_state.clipping_plane.normal
            };
            
            // Compute section cuts if cache is invalid
            let mesh_version = app.render_cache.mesh_version;
            let results_version = if results.is_some() { 1u64 } else { 0u64 };
            
            let param_hash = section_cut::SectionCutCache::compute_hash(
                ui_state.clipping_plane.position,
                ui_state.clipping_plane.axis,
                ui_state.clipping_plane.flip,
                ui_state.displacement_scale,
                ui_state.color_mode,
                ui_state.stress_component,
                ui_state.strain_component,
                mesh_version,
                results_version,
            );
            
            // Check if we need to recompute
            if !app.section_cut_cache.is_valid(param_hash) {
                let polygons = section_cut::compute_section_cuts(
                    mesh_state,
                    results,
                    ui_state.clipping_plane.axis,
                    clip_normal,
                    ui_state.clipping_plane.position,
                    ui_state.displacement_scale,
                    ui_state.color_mode,
                    ui_state.stress_component,
                    ui_state.strain_component,
                );
                app.section_cut_cache.update(polygons, param_hash);
            }
            
            if app.section_cut_cache.polygons.is_empty() {
                ui.centered_and_justified(|ui| {
                    ui.label("No section cut at current position");
                });
                return;
            }
            
            // Extract values we need after the ui.horizontal closure
            // to avoid holding the borrow of app.state.ui_state across it
            let color_mode = app.state.ui_state.color_mode;
            let clipping_axis = app.state.ui_state.clipping_plane.axis;
            
            // Show field selector and export button
            ui.horizontal(|ui| {
                ui.label("Field:");
                
                // Field selector dropdown
                let field_options = [
                    (crate::state::ColorMode::Solid, "Solid"),
                    (crate::state::ColorMode::Displacement, "Displacement"),
                    (crate::state::ColorMode::VonMises, "Von Mises"),
                    (crate::state::ColorMode::Stress, "Stress"),
                    (crate::state::ColorMode::Strain, "Strain"),
                ];
                
                let current_name = field_options.iter()
                    .find(|(mode, _)| *mode == app.state.ui_state.color_mode)
                    .map(|(_, name)| *name)
                    .unwrap_or("Solid");
                
                egui::ComboBox::from_id_salt("section_field_selector")
                    .selected_text(current_name)
                    .show_ui(ui, |ui| {
                        for (mode, name) in field_options {
                            if ui.selectable_value(&mut app.state.ui_state.color_mode, mode, name).changed() {
                                // Invalidate cache when color mode changes
                                app.section_cut_cache.invalidate_hash();
                                app.render_cache.invalidate();
                            }
                        }
                    });
                
                #[cfg(feature = "native")]
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("📷 Export PNG").clicked() {
                        app.state.ui_state.section_export_requested = true;
                    }
                });
                
                #[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    if ui.button("📷 Export PNG").clicked() {
                        app.state.ui_state.section_export_requested = true;
                    }
                });
            });
            ui.separator();
            
            // Get viewport rect
            let available = ui.available_size();
            let (rect, _response) = ui.allocate_exact_size(available, egui::Sense::hover());
            
            let painter = ui.painter_at(rect);
            
            // Fill background
            painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(30, 30, 35));
            
            // Compute bounds of the section in the plane's local coordinates
            let (min_u, max_u, min_v, max_v) = compute_section_bounds_2d(
                &app.section_cut_cache.polygons,
                clipping_axis,
            );
            
            let range_u = (max_u - min_u).max(0.001);
            let range_v = (max_v - min_v).max(0.001);
            
            // Add margin
            let margin = 0.05;
            let padded_range_u = range_u * (1.0 + 2.0 * margin);
            let padded_range_v = range_v * (1.0 + 2.0 * margin);
            let center_u = (min_u + max_u) * 0.5;
            let center_v = (min_v + max_v) * 0.5;
            
            // Compute scale to fit within rect while maintaining aspect ratio
            let scale_u = rect.width() / padded_range_u;
            let scale_v = rect.height() / padded_range_v;
            let scale = scale_u.min(scale_v);
            
            // Transform function: 2D section coords -> screen coords
            let transform = |u: f32, v: f32| -> egui::Pos2 {
                let screen_x = rect.center().x + (u - center_u) * scale;
                let screen_y = rect.center().y - (v - center_v) * scale;  // Flip Y
                egui::pos2(screen_x, screen_y)
            };
            
            // Get field range for coloring
            let (min_val, max_val) = compute_field_range_from_nodal_values(&app.section_cut_cache.polygons);
            let val_range = (max_val - min_val).max(1e-10);
            
            // Draw section cut polygons using fan triangulation from centroid
            // Each subdivision point uses shape function interpolation for accurate field values
            for polygon in &app.section_cut_cache.polygons {
                let n_verts = polygon.vertices.len();
                if n_verts < 3 {
                    continue;
                }
                
                // Store vertex screen positions and parametric coords
                let vertex_data: Vec<(egui::Pos2, f32, f32, f32)> = polygon.vertices.iter()
                    .map(|v| {
                        let (u, v_coord) = project_to_plane_coords(v.position, clipping_axis);
                        let screen_pos = transform(u, v_coord);
                        (screen_pos, v.xi, v.eta, v.zeta)
                    })
                    .collect();
                
                // Compute centroid screen position and parametric coords
                let centroid_pos = egui::pos2(
                    vertex_data.iter().map(|(p, _, _, _)| p.x).sum::<f32>() / n_verts as f32,
                    vertex_data.iter().map(|(p, _, _, _)| p.y).sum::<f32>() / n_verts as f32,
                );
                let centroid_xi = polygon.vertices.iter().map(|v| v.xi).sum::<f32>() / n_verts as f32;
                let centroid_eta = polygon.vertices.iter().map(|v| v.eta).sum::<f32>() / n_verts as f32;
                let centroid_zeta = polygon.vertices.iter().map(|v| v.zeta).sum::<f32>() / n_verts as f32;
                
                // Draw fan triangles from centroid to each edge
                for i in 0..n_verts {
                    let j = (i + 1) % n_verts;
                    let (p0, xi0, eta0, zeta0) = vertex_data[i];
                    let (p1, xi1, eta1, zeta1) = vertex_data[j];
                    
                    // Subdivide this triangle for smoother gradients
                    let subdiv = 6; // Higher subdivision for smoother appearance
                    
                    // Helper to get screen position and parametric coords for a barycentric point
                    let get_point_data = |ti: f32, tj: f32| -> (egui::Pos2, f32, f32, f32) {
                        let tk = 1.0 - ti - tj;
                        let px = centroid_pos.x * tk + p0.x * ti + p1.x * tj;
                        let py = centroid_pos.y * tk + p0.y * ti + p1.y * tj;
                        
                        // Interpolate parametric coords
                        let xi = centroid_xi * tk + xi0 * ti + xi1 * tj;
                        let eta = centroid_eta * tk + eta0 * ti + eta1 * tj;
                        let zeta = centroid_zeta * tk + zeta0 * ti + zeta1 * tj;
                        
                        (egui::pos2(px, py), xi, eta, zeta)
                    };
                    
                    // Helper to compute color using shape function interpolation
                    let get_color = |xi: f32, eta: f32, zeta: f32| -> egui::Color32 {
                        let field_val = section_cut::interpolate_field_value(polygon, xi, eta, zeta);
                        let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
                        value_to_color_egui(t)
                    };
                    
                    let step = 1.0 / subdiv as f32;
                    
                    // Create subdivision triangles
                    for si in 0..subdiv {
                        for sj in 0..(subdiv - si) {
                            let t0 = si as f32 / subdiv as f32;
                            let t1 = sj as f32 / subdiv as f32;
                            
                            // First triangle
                            let (pa, xi_a, eta_a, zeta_a) = get_point_data(t0, t1);
                            let (pb, xi_b, eta_b, zeta_b) = get_point_data(t0 + step, t1);
                            let (pc, xi_c, eta_c, zeta_c) = get_point_data(t0, t1 + step);
                            
                            let ca = get_color(xi_a, eta_a, zeta_a);
                            let cb = get_color(xi_b, eta_b, zeta_b);
                            let cc = get_color(xi_c, eta_c, zeta_c);
                            
                            // Slightly expand triangle from centroid to eliminate sub-pixel gaps
                            let expand_triangle = |p0: egui::Pos2, p1: egui::Pos2, p2: egui::Pos2, expand: f32| -> Vec<egui::Pos2> {
                                let cx = (p0.x + p1.x + p2.x) / 3.0;
                                let cy = (p0.y + p1.y + p2.y) / 3.0;
                                vec![
                                    egui::pos2(p0.x + (p0.x - cx) * expand, p0.y + (p0.y - cy) * expand),
                                    egui::pos2(p1.x + (p1.x - cx) * expand, p1.y + (p1.y - cy) * expand),
                                    egui::pos2(p2.x + (p2.x - cx) * expand, p2.y + (p2.y - cy) * expand),
                                ]
                            };
                            
                            let avg_color = average_colors_3(ca, cb, cc);
                            painter.add(egui::Shape::convex_polygon(
                                expand_triangle(pa, pb, pc, 0.02),
                                avg_color,
                                egui::Stroke::NONE,
                            ));
                            
                            // Second triangle (if not on the hypotenuse)
                            if si + sj + 1 < subdiv {
                                let (pd, xi_d, eta_d, zeta_d) = get_point_data(t0 + step, t1 + step);
                                let cd = get_color(xi_d, eta_d, zeta_d);
                                let avg_color2 = average_colors_3(cb, cd, cc);
                                painter.add(egui::Shape::convex_polygon(
                                    expand_triangle(pb, pd, pc, 0.02),
                                    avg_color2,
                                    egui::Stroke::NONE,
                                ));
                            }
                        }
                    }
                }
                
                // No element boundary outlines - for seamless appearance
            }
            
            // Draw axis labels
            let label_color = egui::Color32::from_gray(180);
            painter.text(
                egui::pos2(rect.right() - 20.0, rect.center().y),
                egui::Align2::CENTER_CENTER,
                h_axis_label,
                egui::FontId::proportional(14.0),
                label_color,
            );
            painter.text(
                egui::pos2(rect.center().x, rect.top() + 15.0),
                egui::Align2::CENTER_CENTER,
                v_axis_label,
                egui::FontId::proportional(14.0),
                label_color,
            );
            
            // Draw color legend
            draw_color_legend(&painter, rect, min_val, max_val, &color_mode);
        });
    
    app.state.ui_state.clipping_plane.show_2d_view = open;
    
    // Handle section export request (outside the Window closure to avoid borrow conflicts)
    #[cfg(feature = "native")]
    if app.state.ui_state.section_export_requested {
        app.state.ui_state.section_export_requested = false;
        export_2d_section_png(app, 1920, 1080);
    }
    
    #[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
    if app.state.ui_state.section_export_requested {
        app.state.ui_state.section_export_requested = false;
        export_2d_section_png_wasm(app, 1920, 1080);
    }
}

/// Export the 2D section view as a PNG image
#[cfg(feature = "native")]
fn export_2d_section_png(app: &mut FeaApp, width: u32, height: u32) {
    use crate::section_cut;
    
    let polygons = &app.section_cut_cache.polygons;
    if polygons.is_empty() {
        app.state.status_message = "No section cut to export".to_string();
        return;
    }
    
    let axis = app.state.ui_state.clipping_plane.axis;
    let color_mode = app.state.ui_state.color_mode;
    
    // Compute bounds
    let (min_u, max_u, min_v, max_v) = compute_section_bounds_2d(polygons, axis);
    let range_u = (max_u - min_u).max(0.001);
    let range_v = (max_v - min_v).max(0.001);
    
    // Add margin
    let margin = 0.05;
    let padded_range_u = range_u * (1.0 + 2.0 * margin);
    let padded_range_v = range_v * (1.0 + 2.0 * margin);
    let center_u = (min_u + max_u) * 0.5;
    let center_v = (min_v + max_v) * 0.5;
    
    // Compute scale to fit while maintaining aspect ratio
    let scale_u = width as f32 / padded_range_u;
    let scale_v = height as f32 / padded_range_v;
    let scale = scale_u.min(scale_v);
    
    // Get field range for coloring
    let (min_val, max_val) = compute_field_range_from_nodal_values(polygons);
    let val_range = (max_val - min_val).max(1e-10);
    
    // Create image buffer
    let mut img_buffer = image::RgbaImage::from_pixel(width, height, image::Rgba([30, 30, 35, 255]));
    
    // Transform function: 2D section coords -> pixel coords
    let transform = |u: f32, v: f32| -> (f32, f32) {
        let px = (width as f32 / 2.0) + (u - center_u) * scale;
        let py = (height as f32 / 2.0) - (v - center_v) * scale;  // Flip Y
        (px, py)
    };
    
    // Helper to get color for a parametric point
    let get_color = |polygon: &section_cut::CutPolygon, xi: f32, eta: f32, zeta: f32| -> [u8; 4] {
        let field_val = section_cut::interpolate_field_value(polygon, xi, eta, zeta);
        let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
        let (r, g, b) = value_to_color_rgb(t);
        [r, g, b, 255]
    };
    
    // Draw all polygons using fan triangulation
    for polygon in polygons {
        let n_verts = polygon.vertices.len();
        if n_verts < 3 {
            continue;
        }
        
        // Compute centroid
        let centroid_pos: (f32, f32) = {
            let sum_u: f32 = polygon.vertices.iter()
                .map(|v| project_to_plane_coords(v.position, axis).0)
                .sum();
            let sum_v: f32 = polygon.vertices.iter()
                .map(|v| project_to_plane_coords(v.position, axis).1)
                .sum();
            (sum_u / n_verts as f32, sum_v / n_verts as f32)
        };
        let centroid_xi = polygon.vertices.iter().map(|v| v.xi).sum::<f32>() / n_verts as f32;
        let centroid_eta = polygon.vertices.iter().map(|v| v.eta).sum::<f32>() / n_verts as f32;
        let centroid_zeta = polygon.vertices.iter().map(|v| v.zeta).sum::<f32>() / n_verts as f32;
        
        // Draw fan triangles
        for i in 0..n_verts {
            let j = (i + 1) % n_verts;
            
            let v0 = &polygon.vertices[i];
            let v1 = &polygon.vertices[j];
            let (u0, uv0) = project_to_plane_coords(v0.position, axis);
            let (u1, uv1) = project_to_plane_coords(v1.position, axis);
            
            // Subdivide triangle for smooth gradients
            let subdiv = 8;
            let step = 1.0 / subdiv as f32;
            
            for si in 0..subdiv {
                for sj in 0..(subdiv - si) {
                    let t0 = si as f32 / subdiv as f32;
                    let t1 = sj as f32 / subdiv as f32;
                    
                    // First sub-triangle
                    let points_a = [
                        (t0, t1),
                        (t0 + step, t1),
                        (t0, t1 + step),
                    ];
                    
                    draw_sub_triangle(
                        &mut img_buffer,
                        &points_a,
                        centroid_pos, (u0, uv0), (u1, uv1),
                        centroid_xi, centroid_eta, centroid_zeta,
                        v0.xi, v0.eta, v0.zeta,
                        v1.xi, v1.eta, v1.zeta,
                        &transform,
                        polygon,
                        &get_color,
                    );
                    
                    // Second sub-triangle (if not on hypotenuse)
                    if si + sj + 1 < subdiv {
                        let points_b = [
                            (t0 + step, t1),
                            (t0 + step, t1 + step),
                            (t0, t1 + step),
                        ];
                        draw_sub_triangle(
                            &mut img_buffer,
                            &points_b,
                            centroid_pos, (u0, uv0), (u1, uv1),
                            centroid_xi, centroid_eta, centroid_zeta,
                            v0.xi, v0.eta, v0.zeta,
                            v1.xi, v1.eta, v1.zeta,
                            &transform,
                            polygon,
                            &get_color,
                        );
                    }
                }
            }
        }
    }
    
    // Draw color legend
    draw_color_legend_to_image(&mut img_buffer, min_val, max_val, &color_mode);
    
    // Draw axis labels
    draw_axis_labels_to_image(&mut img_buffer, axis, width, height);
    
    // Get PNG data
    let mut png_data = Vec::new();
    {
        use image::codecs::png::PngEncoder;
        use image::ImageEncoder;
        let encoder = PngEncoder::new(&mut png_data);
        encoder.write_image(
            img_buffer.as_raw(),
            width,
            height,
            image::ExtendedColorType::Rgba8,
        ).ok();
    }
    
    // Save/download the image
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("PNG Image", &["png"])
            .set_file_name("section_cut.png")
            .save_file()
        {
            if let Err(e) = std::fs::write(&path, &png_data) {
                app.state.status_message = format!("Failed to save image: {}", e);
            } else {
                app.state.status_message = format!("Section view exported to {:?}", path);
            }
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        crate::web_file_io::download_file("section_cut.png", &png_data, "image/png");
        app.state.status_message = "Section view exported as section_cut.png".to_string();
    }
}

/// Draw a subdivided triangle to the image buffer
#[cfg(feature = "native")]
fn draw_sub_triangle<F, G>(
    img: &mut image::RgbaImage,
    bary_points: &[(f32, f32); 3],
    centroid: (f32, f32),
    v0: (f32, f32),
    v1: (f32, f32),
    c_xi: f32, c_eta: f32, c_zeta: f32,
    xi0: f32, eta0: f32, zeta0: f32,
    xi1: f32, eta1: f32, zeta1: f32,
    transform: &F,
    polygon: &crate::section_cut::CutPolygon,
    get_color: &G,
)
where
    F: Fn(f32, f32) -> (f32, f32),
    G: Fn(&crate::section_cut::CutPolygon, f32, f32, f32) -> [u8; 4],
{
    // Convert barycentric to world then screen coords
    let get_point = |ti: f32, tj: f32| -> ((f32, f32), f32, f32, f32) {
        let tk = 1.0 - ti - tj;
        let u = centroid.0 * tk + v0.0 * ti + v1.0 * tj;
        let v = centroid.1 * tk + v0.1 * ti + v1.1 * tj;
        let xi = c_xi * tk + xi0 * ti + xi1 * tj;
        let eta = c_eta * tk + eta0 * ti + eta1 * tj;
        let zeta = c_zeta * tk + zeta0 * ti + zeta1 * tj;
        (transform(u, v), xi, eta, zeta)
    };
    
    let (p0, xi_a, eta_a, zeta_a) = get_point(bary_points[0].0, bary_points[0].1);
    let (p1, xi_b, eta_b, zeta_b) = get_point(bary_points[1].0, bary_points[1].1);
    let (p2, xi_c, eta_c, zeta_c) = get_point(bary_points[2].0, bary_points[2].1);
    
    // Average color for the triangle
    let c0 = get_color(polygon, xi_a, eta_a, zeta_a);
    let c1 = get_color(polygon, xi_b, eta_b, zeta_b);
    let c2 = get_color(polygon, xi_c, eta_c, zeta_c);
    let avg_color = [
        ((c0[0] as u32 + c1[0] as u32 + c2[0] as u32) / 3) as u8,
        ((c0[1] as u32 + c1[1] as u32 + c2[1] as u32) / 3) as u8,
        ((c0[2] as u32 + c1[2] as u32 + c2[2] as u32) / 3) as u8,
        255,
    ];
    
    // Rasterize the triangle
    fill_triangle_to_image(img, p0, p1, p2, avg_color);
}

/// Fill a triangle with a solid color using scanline rasterization
#[cfg(feature = "native")]
fn fill_triangle_to_image(
    img: &mut image::RgbaImage,
    p0: (f32, f32),
    p1: (f32, f32),
    p2: (f32, f32),
    color: [u8; 4],
) {
    let width = img.width() as i32;
    let height = img.height() as i32;
    
    // Bounding box
    let min_x = (p0.0.min(p1.0).min(p2.0).floor() as i32).max(0);
    let max_x = (p0.0.max(p1.0).max(p2.0).ceil() as i32).min(width - 1);
    let min_y = (p0.1.min(p1.1).min(p2.1).floor() as i32).max(0);
    let max_y = (p0.1.max(p1.1).max(p2.1).ceil() as i32).min(height - 1);
    
    // Edge function for point-in-triangle test
    let edge = |a: (f32, f32), b: (f32, f32), p: (f32, f32)| -> f32 {
        (p.0 - a.0) * (b.1 - a.1) - (p.1 - a.1) * (b.0 - a.0)
    };
    
    let area = edge(p0, p1, p2);
    if area.abs() < 0.001 {
        return; // Degenerate triangle
    }
    
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let p = (x as f32 + 0.5, y as f32 + 0.5);
            let w0 = edge(p1, p2, p);
            let w1 = edge(p2, p0, p);
            let w2 = edge(p0, p1, p);
            
            // Check if point is inside triangle (same sign for all edges)
            if (w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0) || (w0 <= 0.0 && w1 <= 0.0 && w2 <= 0.0) {
                img.put_pixel(x as u32, y as u32, image::Rgba(color));
            }
        }
    }
}

/// Convert normalized value (0-1) to RGB tuple
fn value_to_color_rgb(t: f32) -> (u8, u8, u8) {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    ((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

/// Draw color legend to image
#[cfg(feature = "native")]
fn draw_color_legend_to_image(
    img: &mut image::RgbaImage,
    _min_val: f64,
    _max_val: f64,
    _color_mode: &ColorMode,
) {
    let _width = img.width();
    let height = img.height();
    
    let legend_width = 25;
    let legend_height = (height as f32 * 0.5) as u32;
    let legend_x = 15;
    let legend_top = (height - legend_height) / 2;
    
    // Draw gradient bar
    for y in 0..legend_height {
        let t = 1.0 - (y as f32 / legend_height as f32);
        let (r, g, b) = value_to_color_rgb(t);
        for x in 0..legend_width {
            img.put_pixel(legend_x + x, legend_top + y, image::Rgba([r, g, b, 255]));
        }
    }
    
    // Draw border
    let border_color = image::Rgba([100, 100, 100, 255]);
    for x in 0..legend_width {
        img.put_pixel(legend_x + x, legend_top, border_color);
        img.put_pixel(legend_x + x, legend_top + legend_height - 1, border_color);
    }
    for y in 0..legend_height {
        img.put_pixel(legend_x, legend_top + y, border_color);
        img.put_pixel(legend_x + legend_width - 1, legend_top + y, border_color);
    }
    
    // Note: text rendering would require a font library, keeping it simple
    // The values are communicated by the gradient itself
}

/// Draw axis labels to image (simplified - just indicators in corners)
#[cfg(feature = "native")]
fn draw_axis_labels_to_image(
    img: &mut image::RgbaImage,
    axis: crate::state::ClipAxis,
    _width: u32,
    height: u32,
) {
    use crate::state::ClipAxis;
    
    let (_h_label, _v_label) = match axis {
        ClipAxis::X => ("Y", "Z"),
        ClipAxis::Y => ("X", "Z"),
        ClipAxis::Z => ("X", "Y"),
        ClipAxis::Custom => ("U", "V"),
    };
    
    // Draw small axis indicator arrows in the bottom-left corner
    let arrow_color = image::Rgba([180, 180, 180, 255]);
    let origin_x = 50;
    let origin_y = height - 50;
    let arrow_len = 30;
    
    // Horizontal arrow (right)
    for x in 0..arrow_len {
        img.put_pixel(origin_x + x, origin_y, arrow_color);
    }
    // Arrow head
    img.put_pixel(origin_x + arrow_len - 2, origin_y - 1, arrow_color);
    img.put_pixel(origin_x + arrow_len - 2, origin_y + 1, arrow_color);
    img.put_pixel(origin_x + arrow_len - 3, origin_y - 2, arrow_color);
    img.put_pixel(origin_x + arrow_len - 3, origin_y + 2, arrow_color);
    
    // Vertical arrow (up)
    for y in 0..arrow_len {
        img.put_pixel(origin_x, origin_y - y, arrow_color);
    }
    // Arrow head
    img.put_pixel(origin_x - 1, origin_y - arrow_len + 2, arrow_color);
    img.put_pixel(origin_x + 1, origin_y - arrow_len + 2, arrow_color);
    img.put_pixel(origin_x - 2, origin_y - arrow_len + 3, arrow_color);
    img.put_pixel(origin_x + 2, origin_y - arrow_len + 3, arrow_color);
}

/// Project a 3D point to 2D coordinates on the section plane
fn project_to_plane_coords(pos: [f32; 3], axis: crate::state::ClipAxis) -> (f32, f32) {
    use crate::state::ClipAxis;
    match axis {
        ClipAxis::X => (pos[1], pos[2]),  // Y-Z plane
        ClipAxis::Y => (pos[0], pos[2]),  // X-Z plane
        ClipAxis::Z => (pos[0], pos[1]),  // X-Y plane
        ClipAxis::Custom => (pos[0], pos[1]),  // Default to X-Y
    }
}

/// Export the 2D section view as a PNG image (WASM version)
#[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
fn export_2d_section_png_wasm(app: &mut FeaApp, width: u32, height: u32) {
    use crate::section_cut;
    
    let polygons = &app.section_cut_cache.polygons;
    if polygons.is_empty() {
        app.state.status_message = "No section cut to export".to_string();
        return;
    }
    
    let axis = app.state.ui_state.clipping_plane.axis;
    
    // Compute bounds
    let (min_u, max_u, min_v, max_v) = compute_section_bounds_2d(polygons, axis);
    let range_u = (max_u - min_u).max(0.001);
    let range_v = (max_v - min_v).max(0.001);
    
    // Add margin
    let margin = 0.05;
    let padded_range_u = range_u * (1.0 + 2.0 * margin);
    let padded_range_v = range_v * (1.0 + 2.0 * margin);
    let center_u = (min_u + max_u) * 0.5;
    let center_v = (min_v + max_v) * 0.5;
    
    // Compute scale to fit while maintaining aspect ratio
    let scale_u = width as f32 / padded_range_u;
    let scale_v = height as f32 / padded_range_v;
    let scale = scale_u.min(scale_v);
    
    // Get field range for coloring
    let (min_val, max_val) = compute_field_range_from_nodal_values(polygons);
    let val_range = (max_val - min_val).max(1e-10);
    
    // Create RGBA buffer (4 bytes per pixel)
    let mut buffer: Vec<u8> = vec![0; (width * height * 4) as usize];
    
    // Fill with background color
    for pixel in buffer.chunks_exact_mut(4) {
        pixel[0] = 30;  // R
        pixel[1] = 30;  // G
        pixel[2] = 35;  // B
        pixel[3] = 255; // A
    }
    
    // Transform function: 2D section coords -> pixel coords
    let transform = |u: f32, v: f32| -> (f32, f32) {
        let px = (width as f32 / 2.0) + (u - center_u) * scale;
        let py = (height as f32 / 2.0) - (v - center_v) * scale;  // Flip Y
        (px, py)
    };
    
    // Helper to get color for a parametric point
    let get_color = |polygon: &section_cut::CutPolygon, xi: f32, eta: f32, zeta: f32| -> [u8; 4] {
        let field_val = section_cut::interpolate_field_value(polygon, xi, eta, zeta);
        let t = ((field_val - min_val) / val_range).clamp(0.0, 1.0) as f32;
        let (r, g, b) = value_to_color_rgb(t);
        [r, g, b, 255]
    };
    
    // Draw all polygons using fan triangulation
    for polygon in polygons {
        let n_verts = polygon.vertices.len();
        if n_verts < 3 {
            continue;
        }
        
        // Compute centroid
        let centroid_pos: (f32, f32) = {
            let sum_u: f32 = polygon.vertices.iter()
                .map(|v| project_to_plane_coords(v.position, axis).0)
                .sum();
            let sum_v: f32 = polygon.vertices.iter()
                .map(|v| project_to_plane_coords(v.position, axis).1)
                .sum();
            (sum_u / n_verts as f32, sum_v / n_verts as f32)
        };
        let centroid_xi = polygon.vertices.iter().map(|v| v.xi).sum::<f32>() / n_verts as f32;
        let centroid_eta = polygon.vertices.iter().map(|v| v.eta).sum::<f32>() / n_verts as f32;
        let centroid_zeta = polygon.vertices.iter().map(|v| v.zeta).sum::<f32>() / n_verts as f32;
        
        // Draw fan triangles
        for i in 0..n_verts {
            let j = (i + 1) % n_verts;
            
            let v0 = &polygon.vertices[i];
            let v1 = &polygon.vertices[j];
            let (u0, uv0) = project_to_plane_coords(v0.position, axis);
            let (u1, uv1) = project_to_plane_coords(v1.position, axis);
            
            // Subdivide triangle for smooth gradients
            let subdiv = 6;  // Slightly less than native for performance
            let step = 1.0 / subdiv as f32;
            
            for si in 0..subdiv {
                for sj in 0..(subdiv - si) {
                    let t0 = si as f32 / subdiv as f32;
                    let t1 = sj as f32 / subdiv as f32;
                    
                    // First sub-triangle
                    draw_sub_triangle_wasm(
                        &mut buffer, width, height,
                        [(t0, t1), (t0 + step, t1), (t0, t1 + step)],
                        centroid_pos, (u0, uv0), (u1, uv1),
                        centroid_xi, centroid_eta, centroid_zeta,
                        v0.xi, v0.eta, v0.zeta,
                        v1.xi, v1.eta, v1.zeta,
                        &transform, polygon, &get_color,
                    );
                    
                    // Second sub-triangle (if not on hypotenuse)
                    if si + sj + 1 < subdiv {
                        draw_sub_triangle_wasm(
                            &mut buffer, width, height,
                            [(t0 + step, t1), (t0 + step, t1 + step), (t0, t1 + step)],
                            centroid_pos, (u0, uv0), (u1, uv1),
                            centroid_xi, centroid_eta, centroid_zeta,
                            v0.xi, v0.eta, v0.zeta,
                            v1.xi, v1.eta, v1.zeta,
                            &transform, polygon, &get_color,
                        );
                    }
                }
            }
        }
    }
    
    // Draw color legend
    draw_color_legend_wasm(&mut buffer, width, height);
    
    // Encode as PNG
    let mut png_data = Vec::new();
    {
        let mut encoder = png::Encoder::new(&mut png_data, width, height);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        if let Ok(mut writer) = encoder.write_header() {
            if writer.write_image_data(&buffer).is_err() {
                app.state.status_message = "Failed to encode PNG".to_string();
                return;
            }
        } else {
            app.state.status_message = "Failed to create PNG encoder".to_string();
            return;
        }
    }
    
    // Trigger browser download
    crate::web_file_io::download_file("section_cut.png", &png_data, "image/png");
    app.state.status_message = "Section view exported as section_cut.png".to_string();
}

/// Draw a subdivided triangle to buffer (WASM version)
#[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
fn draw_sub_triangle_wasm<F, G>(
    buffer: &mut [u8],
    width: u32,
    height: u32,
    bary_points: [(f32, f32); 3],
    centroid: (f32, f32),
    v0: (f32, f32),
    v1: (f32, f32),
    c_xi: f32, c_eta: f32, c_zeta: f32,
    xi0: f32, eta0: f32, zeta0: f32,
    xi1: f32, eta1: f32, zeta1: f32,
    transform: &F,
    polygon: &crate::section_cut::CutPolygon,
    get_color: &G,
)
where
    F: Fn(f32, f32) -> (f32, f32),
    G: Fn(&crate::section_cut::CutPolygon, f32, f32, f32) -> [u8; 4],
{
    let get_point = |ti: f32, tj: f32| -> ((f32, f32), f32, f32, f32) {
        let tk = 1.0 - ti - tj;
        let u = centroid.0 * tk + v0.0 * ti + v1.0 * tj;
        let v = centroid.1 * tk + v0.1 * ti + v1.1 * tj;
        let xi = c_xi * tk + xi0 * ti + xi1 * tj;
        let eta = c_eta * tk + eta0 * ti + eta1 * tj;
        let zeta = c_zeta * tk + zeta0 * ti + zeta1 * tj;
        (transform(u, v), xi, eta, zeta)
    };
    
    let (p0, xi_a, eta_a, zeta_a) = get_point(bary_points[0].0, bary_points[0].1);
    let (p1, xi_b, eta_b, zeta_b) = get_point(bary_points[1].0, bary_points[1].1);
    let (p2, xi_c, eta_c, zeta_c) = get_point(bary_points[2].0, bary_points[2].1);
    
    // Average color for the triangle
    let c0 = get_color(polygon, xi_a, eta_a, zeta_a);
    let c1 = get_color(polygon, xi_b, eta_b, zeta_b);
    let c2 = get_color(polygon, xi_c, eta_c, zeta_c);
    let avg_color = [
        ((c0[0] as u16 + c1[0] as u16 + c2[0] as u16) / 3) as u8,
        ((c0[1] as u16 + c1[1] as u16 + c2[1] as u16) / 3) as u8,
        ((c0[2] as u16 + c1[2] as u16 + c2[2] as u16) / 3) as u8,
        255,
    ];
    
    // Simple triangle rasterization
    fill_triangle_wasm(buffer, width, height, [p0, p1, p2], avg_color);
}

/// Fill a triangle in the buffer (WASM version)
#[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
fn fill_triangle_wasm(
    buffer: &mut [u8],
    width: u32,
    height: u32,
    points: [(f32, f32); 3],
    color: [u8; 4],
) {
    let min_x = points.iter().map(|p| p.0).fold(f32::MAX, f32::min).max(0.0) as u32;
    let max_x = points.iter().map(|p| p.0).fold(f32::MIN, f32::max).min(width as f32 - 1.0) as u32;
    let min_y = points.iter().map(|p| p.1).fold(f32::MAX, f32::min).max(0.0) as u32;
    let max_y = points.iter().map(|p| p.1).fold(f32::MIN, f32::max).min(height as f32 - 1.0) as u32;
    
    let (x0, y0) = points[0];
    let (x1, y1) = points[1];
    let (x2, y2) = points[2];
    
    let area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    if area.abs() < 0.001 {
        return;
    }
    
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let px = x as f32 + 0.5;
            let py = y as f32 + 0.5;
            
            let w0 = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0);
            let w1 = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1);
            let w2 = (x0 - x2) * (py - y2) - (y0 - y2) * (px - x2);
            
            let inside = if area > 0.0 {
                w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0
            } else {
                w0 <= 0.0 && w1 <= 0.0 && w2 <= 0.0
            };
            
            if inside {
                let idx = ((y * width + x) * 4) as usize;
                if idx + 3 < buffer.len() {
                    buffer[idx] = color[0];
                    buffer[idx + 1] = color[1];
                    buffer[idx + 2] = color[2];
                    buffer[idx + 3] = color[3];
                }
            }
        }
    }
}

/// Draw color legend to buffer (WASM version)
#[cfg(all(feature = "wasm-bindgen", not(feature = "native")))]
fn draw_color_legend_wasm(buffer: &mut [u8], width: u32, height: u32) {
    let legend_width = 25u32;
    let legend_height = (height as f32 * 0.5) as u32;
    let legend_x = 15u32;
    let legend_top = (height - legend_height) / 2;
    
    // Draw gradient bar
    for y in 0..legend_height {
        let t = 1.0 - (y as f32 / legend_height as f32);
        let (r, g, b) = value_to_color_rgb(t);
        for x in 0..legend_width {
            let px = legend_x + x;
            let py = legend_top + y;
            let idx = ((py * width + px) * 4) as usize;
            if idx + 3 < buffer.len() {
                buffer[idx] = r;
                buffer[idx + 1] = g;
                buffer[idx + 2] = b;
                buffer[idx + 3] = 255;
            }
        }
    }
    
    // Draw border
    let border = [100u8, 100, 100, 255];
    for x in 0..legend_width {
        // Top border
        let idx = ((legend_top * width + legend_x + x) * 4) as usize;
        if idx + 3 < buffer.len() {
            buffer[idx..idx+4].copy_from_slice(&border);
        }
        // Bottom border
        let idx = (((legend_top + legend_height - 1) * width + legend_x + x) * 4) as usize;
        if idx + 3 < buffer.len() {
            buffer[idx..idx+4].copy_from_slice(&border);
        }
    }
    for y in 0..legend_height {
        // Left border
        let idx = (((legend_top + y) * width + legend_x) * 4) as usize;
        if idx + 3 < buffer.len() {
            buffer[idx..idx+4].copy_from_slice(&border);
        }
        // Right border
        let idx = (((legend_top + y) * width + legend_x + legend_width - 1) * 4) as usize;
        if idx + 3 < buffer.len() {
            buffer[idx..idx+4].copy_from_slice(&border);
        }
    }
}

/// Compute bounds of section cut in 2D plane coordinates
fn compute_section_bounds_2d(
    polygons: &[crate::section_cut::CutPolygon],
    axis: crate::state::ClipAxis,
) -> (f32, f32, f32, f32) {
    let mut min_u = f32::MAX;
    let mut max_u = f32::MIN;
    let mut min_v = f32::MAX;
    let mut max_v = f32::MIN;
    
    for polygon in polygons {
        for vertex in &polygon.vertices {
            let (u, v) = project_to_plane_coords(vertex.position, axis);
            min_u = min_u.min(u);
            max_u = max_u.max(u);
            min_v = min_v.min(v);
            max_v = max_v.max(v);
        }
    }
    
    if min_u > max_u {
        (0.0, 1.0, 0.0, 1.0)
    } else {
        (min_u, max_u, min_v, max_v)
    }
}

/// Compute field value range from polygons (using nodal values for accurate range)
fn compute_field_range_from_nodal_values(polygons: &[crate::section_cut::CutPolygon]) -> (f64, f64) {
    let mut min_val = f64::MAX;
    let mut max_val = f64::MIN;
    
    for polygon in polygons {
        // Use all nodal field values for accurate range
        for val in &polygon.nodal_field_values {
            min_val = min_val.min(*val);
            max_val = max_val.max(*val);
        }
    }
    
    if min_val > max_val {
        (0.0, 1.0)
    } else {
        (min_val, max_val)
    }
}

/// Compute field value range from polygon vertices (fallback)
fn compute_field_range_from_polygons(polygons: &[crate::section_cut::CutPolygon]) -> (f64, f64) {
    let mut min_val = f64::MAX;
    let mut max_val = f64::MIN;
    
    for polygon in polygons {
        for vertex in &polygon.vertices {
            min_val = min_val.min(vertex.field_value);
            max_val = max_val.max(vertex.field_value);
        }
    }
    
    if min_val > max_val {
        (0.0, 1.0)
    } else {
        (min_val, max_val)
    }
}

/// Point-in-polygon test for 2D coordinates
fn point_in_polygon_2d(point: [f32; 2], vertices: &[[f32; 2]]) -> bool {
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

/// Convert normalized value (0-1) to jet colormap color
fn value_to_color_egui(t: f32) -> egui::Color32 {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    
    egui::Color32::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}

/// Average three colors
fn average_colors_3(c0: egui::Color32, c1: egui::Color32, c2: egui::Color32) -> egui::Color32 {
    let r = ((c0.r() as u32 + c1.r() as u32 + c2.r() as u32) / 3) as u8;
    let g = ((c0.g() as u32 + c1.g() as u32 + c2.g() as u32) / 3) as u8;
    let b = ((c0.b() as u32 + c1.b() as u32 + c2.b() as u32) / 3) as u8;
    let a = ((c0.a() as u32 + c1.a() as u32 + c2.a() as u32) / 3) as u8;
    egui::Color32::from_rgba_unmultiplied(r, g, b, a)
}

/// Convert normalized value (0-1) to jet colormap color
fn value_to_jet_color(t: f32) -> egui::Color32 {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    
    egui::Color32::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}

/// Draw a color legend on the 2D section view
fn draw_color_legend(
    painter: &egui::Painter,
    rect: egui::Rect,
    min_val: f64,
    max_val: f64,
    color_mode: &ColorMode,
) {
    let legend_width = 20.0;
    let legend_height = rect.height() * 0.6;
    let legend_x = rect.left() + 10.0;
    let legend_top = rect.center().y - legend_height * 0.5;
    
    // Draw gradient bar
    let steps = 50;
    let step_height = legend_height / steps as f32;
    
    for i in 0..steps {
        let t = 1.0 - (i as f32 / steps as f32);  // Top is high, bottom is low
        let color = value_to_jet_color(t);
        let y = legend_top + i as f32 * step_height;
        painter.rect_filled(
            egui::Rect::from_min_size(
                egui::pos2(legend_x, y),
                egui::vec2(legend_width, step_height + 0.5),
            ),
            0.0,
            color,
        );
    }
    
    // Draw border
    painter.rect_stroke(
        egui::Rect::from_min_size(
            egui::pos2(legend_x, legend_top),
            egui::vec2(legend_width, legend_height),
        ),
        0.0,
        egui::Stroke::new(1.0_f32, egui::Color32::from_gray(100)),
        egui::StrokeKind::Outside,
    );
    
    // Draw labels
    let label_x = legend_x + legend_width + 5.0;
    let label_color = egui::Color32::from_gray(200);
    let font = egui::FontId::proportional(10.0);
    
    // Max value (top)
    painter.text(
        egui::pos2(label_x, legend_top),
        egui::Align2::LEFT_TOP,
        format_value(max_val),
        font.clone(),
        label_color,
    );
    
    // Min value (bottom)
    painter.text(
        egui::pos2(label_x, legend_top + legend_height),
        egui::Align2::LEFT_BOTTOM,
        format_value(min_val),
        font.clone(),
        label_color,
    );
    
    // Color mode label (rotated title would be nice but just put at top for now)
    painter.text(
        egui::pos2(legend_x + legend_width * 0.5, legend_top - 5.0),
        egui::Align2::CENTER_BOTTOM,
        format!("{:?}", color_mode),
        egui::FontId::proportional(9.0),
        egui::Color32::from_gray(150),
    );
}

/// Format a value for display in the legend
fn format_value(val: f64) -> String {
    if val.abs() < 0.001 || val.abs() >= 10000.0 {
        format!("{:.2e}", val)
    } else if val.abs() < 1.0 {
        format!("{:.4}", val)
    } else {
        format!("{:.2}", val)
    }
}

/// Get color for a phase category
fn category_color(cat: crate::state::SolvePhaseCategory) -> egui::Color32 {
    use crate::state::SolvePhaseCategory;
    match cat {
        SolvePhaseCategory::Init => egui::Color32::from_rgb(100, 149, 237),      // Cornflower blue
        SolvePhaseCategory::Assembly => egui::Color32::from_rgb(255, 165, 0),    // Orange
        SolvePhaseCategory::Solve => egui::Color32::from_rgb(50, 205, 50),       // Lime green
        SolvePhaseCategory::PostProcess => egui::Color32::from_rgb(186, 85, 211), // Medium orchid
        SolvePhaseCategory::Other => egui::Color32::GRAY,
    }
}

/// Get icon for a phase category
fn category_icon(cat: crate::state::SolvePhaseCategory) -> &'static str {
    use crate::state::SolvePhaseCategory;
    match cat {
        SolvePhaseCategory::Init => "⚙",
        SolvePhaseCategory::Assembly => "🔧",
        SolvePhaseCategory::Solve => "📊",
        SolvePhaseCategory::PostProcess => "📈",
        SolvePhaseCategory::Other => "•",
    }
}

fn setup_custom_style(ctx: &egui::Context) {
    use egui::{Color32, FontId, FontFamily, CornerRadius, Stroke, Shadow, Vec2, Margin};
    use egui::style::{Widgets, WidgetVisuals, Selection, HandleShape};
    
    // =========================================================================
    // Custom Fonts - Inter (UI) + JetBrains Mono (code) + Noto Sans Symbols 2 (icons) + Remix Icons
    // =========================================================================
    let mut fonts = egui::FontDefinitions::default();
    
    // Load Remix Icons font for UI icons (2800+ icons)
    fonts.font_data.insert(
        "remix_icons".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/remixicon.ttf")).into(),
    );
    
    // Load Inter font family
    fonts.font_data.insert(
        "inter_regular".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/Inter-Regular.ttf")).into(),
    );
    fonts.font_data.insert(
        "inter_medium".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/Inter-Medium.ttf")).into(),
    );
    fonts.font_data.insert(
        "inter_bold".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/Inter-Bold.ttf")).into(),
    );
    
    // Load JetBrains Mono for monospace
    fonts.font_data.insert(
        "jetbrains_mono".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/JetBrainsMono-Regular.ttf")).into(),
    );
    
    // Load Noto Sans Symbols 2 for geometric shapes and symbols
    fonts.font_data.insert(
        "noto_symbols".to_owned(),
        egui::FontData::from_static(include_bytes!("../assets/NotoSansSymbols2-Regular.ttf")).into(),
    );
    
    // Set Inter as the primary proportional font with Remix Icons and Noto Symbols as fallback
    fonts.families.entry(FontFamily::Proportional).or_default()
        .insert(0, "inter_regular".to_owned());
    fonts.families.entry(FontFamily::Proportional).or_default()
        .push("remix_icons".to_owned());
    fonts.families.entry(FontFamily::Proportional).or_default()
        .push("noto_symbols".to_owned());
    
    // Set JetBrains Mono as the primary monospace font with Remix Icons and Noto Symbols fallback
    fonts.families.entry(FontFamily::Monospace).or_default()
        .insert(0, "jetbrains_mono".to_owned());
    fonts.families.entry(FontFamily::Monospace).or_default()
        .push("remix_icons".to_owned());
    fonts.families.entry(FontFamily::Monospace).or_default()
        .push("noto_symbols".to_owned());
    
    ctx.set_fonts(fonts);
    
    let mut style = (*ctx.style()).clone();
    
    // =========================================================================
    // Typography - Clean, readable fonts
    // =========================================================================
    style.text_styles.insert(
        egui::TextStyle::Heading,
        FontId::new(20.0, FontFamily::Proportional),
    );
    style.text_styles.insert(
        egui::TextStyle::Body,
        FontId::new(14.0, FontFamily::Proportional),
    );
    style.text_styles.insert(
        egui::TextStyle::Button,
        FontId::new(14.0, FontFamily::Proportional),
    );
    style.text_styles.insert(
        egui::TextStyle::Small,
        FontId::new(12.0, FontFamily::Proportional),
    );
    style.text_styles.insert(
        egui::TextStyle::Monospace,
        FontId::new(13.0, FontFamily::Monospace),
    );
    
    // =========================================================================
    // Spacing - More breathing room for modern look
    // =========================================================================
    style.spacing.item_spacing = Vec2::new(8.0, 6.0);
    style.spacing.window_margin = Margin::same(12);
    style.spacing.button_padding = Vec2::new(10.0, 5.0);
    style.spacing.menu_margin = Margin::same(8);
    style.spacing.indent = 20.0;
    style.spacing.interact_size = Vec2::new(44.0, 22.0);  // Slightly taller touch targets
    style.spacing.slider_width = 140.0;
    style.spacing.combo_width = 120.0;
    style.spacing.text_edit_width = 200.0;
    style.spacing.icon_width = 16.0;
    style.spacing.icon_width_inner = 10.0;
    style.spacing.icon_spacing = 6.0;
    style.spacing.tooltip_width = 300.0;
    style.spacing.combo_height = 240.0;
    style.spacing.indent_ends_with_horizontal_line = false;
    
    // Scroll bar styling
    style.spacing.scroll.bar_width = 10.0;
    style.spacing.scroll.handle_min_length = 24.0;
    style.spacing.scroll.bar_inner_margin = 3.0;
    style.spacing.scroll.bar_outer_margin = 2.0;
    
    // =========================================================================
    // Color Palette - Modern dark theme with blue accent
    // =========================================================================
    // Base colors
    let bg_dark = Color32::from_rgb(24, 26, 32);           // Main background
    let bg_medium = Color32::from_rgb(32, 35, 42);         // Panel background
    let bg_light = Color32::from_rgb(42, 46, 56);          // Elevated surfaces
    let bg_hover = Color32::from_rgb(52, 58, 70);          // Hover state
    let bg_active = Color32::from_rgb(62, 68, 82);         // Active/pressed state
    
    // Text colors
    let text_primary = Color32::from_rgb(230, 233, 240);   // Primary text
    let text_secondary = Color32::from_rgb(160, 168, 180); // Secondary/muted text
    
    // Accent colors - Modern blue
    let accent = Color32::from_rgb(66, 133, 244);          // Primary accent (Google blue-ish)
    let accent_hover = Color32::from_rgb(90, 152, 255);    // Lighter on hover
    let accent_muted = Color32::from_rgb(45, 95, 170);     // Subtle accent
    
    // Stroke colors
    let stroke_subtle = Color32::from_rgb(55, 60, 72);     // Subtle borders
    
    // Status colors
    let warn_color = Color32::from_rgb(255, 180, 70);      // Warning orange
    let error_color = Color32::from_rgb(255, 100, 100);    // Error red
    let hyperlink = Color32::from_rgb(100, 170, 255);      // Links
    
    // =========================================================================
    // Visuals - Core visual settings
    // =========================================================================
    let mut visuals = style.visuals.clone();
    visuals.dark_mode = true;
    
    // Window styling
    visuals.window_corner_radius = CornerRadius::same(8);
    visuals.window_shadow = Shadow {
        offset: [0, 4],
        blur: 16,
        spread: 0,
        color: Color32::from_black_alpha(100),
    };
    visuals.window_fill = bg_medium;
    visuals.window_stroke = Stroke::new(1.0_f32, stroke_subtle);
    visuals.window_highlight_topmost = true;
    
    // Panel and popup styling
    visuals.panel_fill = bg_medium;
    visuals.popup_shadow = Shadow {
        offset: [0, 3],
        blur: 12,
        spread: 0,
        color: Color32::from_black_alpha(80),
    };
    visuals.menu_corner_radius = CornerRadius::same(6);
    
    // Background colors
    visuals.extreme_bg_color = bg_dark;
    visuals.faint_bg_color = Color32::from_rgb(28, 30, 38);
    visuals.code_bg_color = Color32::from_rgb(35, 38, 48);
    
    // Text colors
    visuals.override_text_color = None;
    visuals.warn_fg_color = warn_color;
    visuals.error_fg_color = error_color;
    visuals.hyperlink_color = hyperlink;
    
    // UI elements
    visuals.resize_corner_size = 12.0;
    visuals.clip_rect_margin = 3.0;
    visuals.button_frame = true;
    visuals.collapsing_header_frame = false;
    visuals.indent_has_left_vline = true;
    visuals.striped = true;
    visuals.slider_trailing_fill = true;
    visuals.handle_shape = HandleShape::Circle;
    visuals.interact_cursor = Some(egui::CursorIcon::PointingHand);
    visuals.image_loading_spinners = true;
    
    // Selection styling
    visuals.selection = Selection {
        bg_fill: accent_muted,
        stroke: Stroke::new(1.0_f32, accent),
    };
    
    // =========================================================================
    // Widget Visuals - State-specific styling
    // =========================================================================
    visuals.widgets = Widgets {
        // Non-interactive widgets (labels, panel backgrounds)
        noninteractive: WidgetVisuals {
            bg_fill: bg_medium,
            weak_bg_fill: bg_light,
            bg_stroke: Stroke::new(1.0_f32, stroke_subtle),
            corner_radius: CornerRadius::same(4),
            fg_stroke: Stroke::new(1.0_f32, text_secondary),
            expansion: 0.0,
        },
        // Interactive widgets at rest
        inactive: WidgetVisuals {
            bg_fill: bg_light,
            weak_bg_fill: bg_light,
            bg_stroke: Stroke::new(1.0_f32, stroke_subtle),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.0_f32, text_primary),
            expansion: 0.0,
        },
        // Hovered interactive widgets
        hovered: WidgetVisuals {
            bg_fill: bg_hover,
            weak_bg_fill: bg_hover,
            bg_stroke: Stroke::new(1.5_f32, accent),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.5_f32, text_primary),
            expansion: 1.0,
        },
        // Active (clicked/focused) widgets
        active: WidgetVisuals {
            bg_fill: bg_active,
            weak_bg_fill: bg_active,
            bg_stroke: Stroke::new(2.0_f32, accent_hover),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(2.0_f32, text_primary),
            expansion: 1.0,
        },
        // Open (e.g., combo box with menu open)
        open: WidgetVisuals {
            bg_fill: bg_hover,
            weak_bg_fill: bg_hover,
            bg_stroke: Stroke::new(1.5_f32, accent),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.5_f32, text_primary),
            expansion: 1.0,
        },
    };
    
    // Text cursor styling
    visuals.text_cursor.stroke = Stroke::new(2.0_f32, accent);
    visuals.text_cursor.preview = false;
    visuals.text_cursor.blink = true;
    visuals.text_cursor.on_duration = 0.5;
    visuals.text_cursor.off_duration = 0.5;
    
    style.visuals = visuals;
    
    // =========================================================================
    // Interaction - Responsive feel
    // =========================================================================
    style.interaction.tooltip_delay = 0.3;
    style.interaction.show_tooltips_only_when_still = false;
    style.interaction.selectable_labels = true;
    style.interaction.multi_widget_text_select = true;
    
    // Animation settings
    style.animation_time = 0.12;  // Snappy animations
    
    ctx.set_style(style);
}

// ============================================================================
// Simulation runners
// ============================================================================

/// Phase tracker for detailed timing
#[cfg(not(target_arch = "wasm32"))]
struct PhaseTracker {
    start: web_time::Instant,
    phase_start: web_time::Instant,
    completed_phases: Vec<SolvePhaseEntry>,
    current_phase: Option<String>,
    current_category: Option<SolvePhaseCategory>,
    tx: Sender<SimMessage>,
}

#[cfg(not(target_arch = "wasm32"))]
impl PhaseTracker {
    fn new(tx: Sender<SimMessage>) -> Self {
        let now = web_time::Instant::now();
        Self {
            start: now,
            phase_start: now,
            completed_phases: Vec::new(),
            current_phase: None,
            current_category: None,
            tx,
        }
    }
    
    fn begin_phase(&mut self, name: &str, category: SolvePhaseCategory) {
        // End previous phase if any
        if let Some(prev_name) = self.current_phase.take() {
            let duration = self.phase_start.elapsed().as_millis() as u64;
            let start_offset = (self.phase_start - self.start).as_millis() as u64;
            self.completed_phases.push(SolvePhaseEntry {
                name: prev_name,
                category: self.current_category.unwrap_or(SolvePhaseCategory::Other),
                duration_ms: duration,
                details: None,
                start_offset_ms: start_offset,
            });
        }
        
        self.current_phase = Some(name.to_string());
        self.current_category = Some(category);
        self.phase_start = web_time::Instant::now();
        
        self.send_update();
    }
    
    fn set_phase_detail(&mut self, detail: &str) {
        // Send update with detail info
        let progress = SolveProgress {
            current_phase: self.current_phase.clone(),
            current_category: self.current_category,
            completed_phases: self.completed_phases.clone(),
            elapsed_ms: self.start.elapsed().as_millis() as u64,
            step_progress: None,
            estimated_remaining_ms: None,
        };
        let _ = self.tx.send(SimMessage::PhaseUpdate(progress));
        
        // Update status message with detail
        if let Some(ref phase) = self.current_phase {
            let _ = self.tx.send(SimMessage::Progress(0.0, format!("{}: {}", phase, detail)));
        }
    }
    
    fn update_step_progress(&mut self, current: usize, total: usize, extra_info: Option<&str>) {
        let progress = SolveProgress {
            current_phase: self.current_phase.clone(),
            current_category: self.current_category,
            completed_phases: self.completed_phases.clone(),
            elapsed_ms: self.start.elapsed().as_millis() as u64,
            step_progress: Some((current, total)),
            estimated_remaining_ms: if current > 0 {
                let elapsed = self.phase_start.elapsed().as_millis() as u64;
                let per_step = elapsed / current as u64;
                Some(per_step * (total - current) as u64)
            } else {
                None
            },
        };
        let _ = self.tx.send(SimMessage::PhaseUpdate(progress));
        
        // Also send progress message
        let pct = current as f32 / total as f32;
        let msg = match extra_info {
            Some(info) => format!("Step {}/{}: {}", current, total, info),
            None => format!("Step {}/{}", current, total),
        };
        let _ = self.tx.send(SimMessage::Progress(pct, msg));
    }
    
    fn send_update(&self) {
        let progress = SolveProgress {
            current_phase: self.current_phase.clone(),
            current_category: self.current_category,
            completed_phases: self.completed_phases.clone(),
            elapsed_ms: self.start.elapsed().as_millis() as u64,
            step_progress: None,
            estimated_remaining_ms: None,
        };
        let _ = self.tx.send(SimMessage::PhaseUpdate(progress));
    }
    
    fn finish(mut self) -> Vec<SolvePhaseEntry> {
        // End final phase
        if let Some(prev_name) = self.current_phase.take() {
            let duration = self.phase_start.elapsed().as_millis() as u64;
            let start_offset = (self.phase_start - self.start).as_millis() as u64;
            self.completed_phases.push(SolvePhaseEntry {
                name: prev_name,
                category: self.current_category.unwrap_or(SolvePhaseCategory::Other),
                duration_ms: duration,
                details: None,
                start_offset_ms: start_offset,
            });
        }
        self.completed_phases
    }
    
    fn total_elapsed_ms(&self) -> u64 {
        self.start.elapsed().as_millis() as u64
    }
}

/// Run simulation in background thread (native only)
#[cfg(not(target_arch = "wasm32"))]
fn run_simulation_threaded(
    mut simulation: rust_fea::simulation::Simulation,
    solver_type: crate::state::SolverType,
    explicit_settings: crate::state::ExplicitSettings,
    tx: Sender<SimMessage>,
) {
    use crate::state::SolvePhaseTimings;
    
    let mut tracker = PhaseTracker::new(tx.clone());
    
    tracker.begin_phase("Initializing simulation", SolvePhaseCategory::Init);
    simulation.initialize();
    simulation.one_time_init();
    
    let result = match solver_type {
        crate::state::SolverType::Direct => {
            run_direct_solver_with_tracking(&mut simulation, &mut tracker)
        }
        crate::state::SolverType::Explicit => {
            run_explicit_solver_with_tracking(&mut simulation, &explicit_settings, &mut tracker)
        }
    };
    
    let phases = tracker.finish();
    let total_ms = phases.iter().map(|p| p.duration_ms).sum();
    
    match result {
        Ok(mut results) => {
            results.stats.solver_time_ms = total_ms;
            results.stats.phase_timing = Some(SolvePhaseTimings {
                total_ms,
                phases,
            });
            let _ = tx.send(SimMessage::Completed(Box::new(results)));
        }
        Err(e) => {
            let _ = tx.send(SimMessage::Error(e));
        }
    }
}

/// Direct solver with detailed phase tracking
#[cfg(not(target_arch = "wasm32"))]
fn run_direct_solver_with_tracking(
    simulation: &mut rust_fea::simulation::Simulation,
    tracker: &mut PhaseTracker,
) -> Result<SimulationResults, String> {
    tracker.begin_phase("Assembling stiffness matrix", SolvePhaseCategory::Assembly);
    
    let (stiffness, load_vector) = simulation.assemble(
        rust_fea::simulation::AssemblyOutputType::SymmetricUpper
    );
    
    let nnz = stiffness.len();
    tracker.set_phase_detail(&format!("{} non-zeros", nnz));
    
    tracker.begin_phase("Converting to sparse format", SolvePhaseCategory::Assembly);
    
    let n = simulation.nodes.len() * simulation.dofs;
    let mut rows = Vec::with_capacity(nnz * 2);
    let mut cols = Vec::with_capacity(nnz * 2);
    let mut vals = Vec::with_capacity(nnz * 2);
    
    for ((i, j), v) in &stiffness {
        rows.push(*i);
        cols.push(*j);
        vals.push(*v);
        if i != j {
            rows.push(*j);
            cols.push(*i);
            vals.push(*v);
        }
    }
    
    tracker.set_phase_detail(&format!("{}x{} matrix, {} entries", n, n, vals.len()));
    
    tracker.begin_phase("Solving linear system (Ax=b)", SolvePhaseCategory::Solve);
    
    let displacement = rust_fea::solver::direct_solve_triplet(
        n, &rows, &cols, &vals, &load_vector,
    ).map_err(|e| format!("Solver failed: {:?}", e))?;
    
    tracker.begin_phase("Updating node displacements", SolvePhaseCategory::PostProcess);
    
    // Update node displacements
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        let dx = displacement[node_id * 3];
        let dy = displacement[node_id * 3 + 1];
        let dz = displacement[node_id * 3 + 2];
        node.set_displacement(dx, dy, dz);
    }
    
    tracker.begin_phase("Computing stress/strain fields", SolvePhaseCategory::PostProcess);
    
    simulation.compute_result_fields();
    
    let (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain) = extract_stress_results(simulation);
    
    tracker.begin_phase("Processing results", SolvePhaseCategory::PostProcess);
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
        nodal_von_mises,
        nodal_stress,
        nodal_strain,
        stats: crate::state::ResultStats {
            max_displacement: max_disp,
            min_displacement: min_disp,
            max_von_mises: max_vm,
            solver_time_ms: 0, // Will be filled in by caller
            phase_timing: None, // Will be filled in by caller
        },
        time_steps: Vec::new(),
    })
}

/// Explicit solver with detailed phase tracking
#[cfg(not(target_arch = "wasm32"))]
fn run_explicit_solver_with_tracking(
    simulation: &mut rust_fea::simulation::Simulation,
    settings: &crate::state::ExplicitSettings,
    tracker: &mut PhaseTracker,
) -> Result<SimulationResults, String> {
    use nalgebra::DVector;
    
    tracker.begin_phase("Computing element matrices", SolvePhaseCategory::Assembly);
    
    simulation.compute_all_element_stiffness();
    simulation.compute_all_element_mass();
    
    tracker.begin_phase("Assembling mass matrix", SolvePhaseCategory::Assembly);
    
    let mass_diag = simulation.compute_global_mass_matrix_diagonal();
    
    // Get material from first element
    let active_el_ids = simulation.active_elements();
    let first_el_id = *active_el_ids.first()
        .ok_or("No active elements")?;
    let material = simulation.get_element(first_el_id)
        .ok_or("Cannot get element")?
        .get_material();
    
    tracker.begin_phase("Computing critical time step", SolvePhaseCategory::Init);
    
    let wave_speed = ((material.youngs_modulus / material.density) 
        * (1.0 - material.poisson_ratio) 
        / ((1.0 + material.poisson_ratio) * (1.0 - 2.0 * material.poisson_ratio))).sqrt();
    
    let dt_crit = simulation.mesh.compute_dt(wave_speed);
    let dt = settings.time_step_override.unwrap_or(dt_crit * 0.5);
    
    tracker.set_phase_detail(&format!("dt={:.2e}s (crit={:.2e}s)", dt, dt_crit));
    
    let n = simulation.nodes.len() * simulation.dofs;
    let dofs = simulation.dofs;
    let total_steps = settings.time_steps;
    let save_interval = settings.vtk_save_steps.max(1);
    
    // Initialize state vectors
    let mut u = DVector::zeros(n);
    let mut u_dot = DVector::zeros(n);
    let mut u_half_dot;
    
    // IMPORTANT: assemble_global_force populates fixed_global_nodal_values via handle_bc()
    // This MUST be called BEFORE get_specified_bc()
    simulation.assemble_global_force();
    let f_ext = simulation.load_vector.clone();
    
    // Get fixed BC DOFs and values (now populated)
    let bc_values = simulation.get_specified_bc();
    
    // Apply initial BCs
    for (dof, val) in &bc_values {
        if *dof < n {
            u[*dof] = *val;
        }
    }
    
    // Update node displacements from initial u vector
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        if node_id * dofs + 2 < n {
            node.set_displacement(u[node_id * dofs], u[node_id * dofs + 1], u[node_id * dofs + 2]);
        }
    }
    
    tracker.begin_phase("Computing initial forces", SolvePhaseCategory::Assembly);
    
    let mut f_int = simulation.compute_force_vector(&u);
    
    let mut time_steps_data = Vec::new();
    let mut current_time = 0.0;
    
    tracker.begin_phase("Time integration", SolvePhaseCategory::Solve);
    
    for step in 0..total_steps {
        // Compute acceleration
        let mut u_ddot = DVector::zeros(n);
        for i in 0..n {
            if mass_diag[i].abs() > 1e-30 {
                u_ddot[i] = (f_ext[i] - f_int[i]) / mass_diag[i];
            }
        }
        
        // Half-step velocity
        u_half_dot = &u_dot + 0.5 * dt * &u_ddot;
        
        // Update displacement
        u = &u + dt * &u_half_dot;
        
        // Apply boundary conditions
        for (dof, val) in &bc_values {
            if *dof < n {
                u[*dof] = *val;
                u_dot[*dof] = 0.0;
                u_half_dot[*dof] = 0.0;
            }
        }
        
        // Update node displacements
        for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
            if node_id * dofs + 2 < n {
                node.set_displacement(u[node_id * dofs], u[node_id * dofs + 1], u[node_id * dofs + 2]);
            }
        }
        
        // Compute new internal forces
        f_int = simulation.compute_force_vector(&u);
        
        // New acceleration
        let mut new_u_ddot = DVector::zeros(n);
        for i in 0..n {
            if mass_diag[i].abs() > 1e-30 {
                new_u_ddot[i] = (f_ext[i] - f_int[i]) / mass_diag[i];
            }
        }
        
        // Complete velocity update
        u_dot = &u_half_dot + 0.5 * dt * &new_u_ddot;
        u_dot *= 0.9995; // Damping
        
        // Apply BC to velocity
        for (dof, _) in &bc_values {
            if *dof < n {
                u_dot[*dof] = 0.0;
            }
        }
        
        current_time += dt;
        
        // Check for instability
        let max_u = u.iter().map(|x| x.abs()).fold(0.0f64, |a, b| a.max(b));
        if max_u.is_nan() || max_u.is_infinite() || max_u > 1e10 {
            return Err(format!(
                "Numerical instability at step {} (t={:.4e}s). max_u={:.2e}", 
                step, current_time, max_u
            ));
        }
        
        // Record at intervals
        if step % save_interval == 0 || step == total_steps - 1 {
            let max_disp = (0..n/dofs).map(|i| {
                let dx = u[i * dofs];
                let dy = u[i * dofs + 1];
                let dz = u[i * dofs + 2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            }).fold(0.0f64, |a, b| a.max(b));
            
            let ke: f64 = (0..n).map(|i| 0.5 * mass_diag[i] * u_dot[i] * u_dot[i]).sum();
            
            time_steps_data.push(crate::state::TimeStepResult {
                time: current_time,
                iteration: step,
                max_displacement: max_disp,
                kinetic_energy: ke,
            });
        }
        
        // Progress updates
        if step % (total_steps / 100).max(1) == 0 {
            let extra = format!("t={:.4e}s, max_u={:.2e}", current_time, max_u);
            tracker.update_step_progress(step, total_steps, Some(&extra));
        }
    }
    
    tracker.begin_phase("Computing final stress fields", SolvePhaseCategory::PostProcess);
    
    let displacement: Vec<f64> = u.iter().cloned().collect();
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        if node_id * 3 + 2 < displacement.len() {
            node.set_displacement(
                displacement[node_id * 3],
                displacement[node_id * 3 + 1],
                displacement[node_id * 3 + 2]
            );
        }
    }
    
    simulation.compute_result_fields();
    
    let (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain) = extract_stress_results(simulation);
    
    tracker.begin_phase("Processing results", SolvePhaseCategory::PostProcess);
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
        nodal_von_mises,
        nodal_stress,
        nodal_strain,
        stats: crate::state::ResultStats {
            max_displacement: max_disp,
            min_displacement: min_disp,
            max_von_mises: max_vm,
            solver_time_ms: 0,
            phase_timing: None,
        },
        time_steps: time_steps_data,
    })
}

/// Run simulation synchronously (WASM) with timing
#[cfg(target_arch = "wasm32")]
fn run_simulation_sync(
    mut simulation: rust_fea::simulation::Simulation,
    solver_type: crate::state::SolverType,
    explicit_settings: crate::state::ExplicitSettings,
) -> Result<SimulationResults, String> {
    use crate::state::{SolvePhaseTimings, SolvePhaseEntry, SolvePhaseCategory};
    
    let total_start = web_time::Instant::now();
    let mut phases: Vec<SolvePhaseEntry> = Vec::new();
    
    // Phase: Initialize
    let phase_start = web_time::Instant::now();
    simulation.initialize();
    simulation.one_time_init();
    phases.push(SolvePhaseEntry {
        name: "Initializing simulation".to_string(),
        category: SolvePhaseCategory::Init,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: None,
        start_offset_ms: 0,
    });
    
    let mut result = match solver_type {
        crate::state::SolverType::Direct => run_direct_solver_impl_timed(&mut simulation, &mut phases, total_start),
        crate::state::SolverType::Explicit => run_explicit_solver_impl_timed(&mut simulation, &explicit_settings, &mut phases, total_start),
    }?;
    
    let total_ms = total_start.elapsed().as_millis() as u64;
    result.stats.solver_time_ms = total_ms;
    result.stats.phase_timing = Some(SolvePhaseTimings { total_ms, phases });
    
    Ok(result)
}

/// Direct solver with timing (WASM)
#[cfg(target_arch = "wasm32")]
fn run_direct_solver_impl_timed(
    simulation: &mut rust_fea::simulation::Simulation,
    phases: &mut Vec<crate::state::SolvePhaseEntry>,
    total_start: web_time::Instant,
) -> Result<SimulationResults, String> {
    use crate::state::{SolvePhaseEntry, SolvePhaseCategory};
    
    // Phase: Assembly
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    let (stiffness, load_vector) = simulation.assemble(
        rust_fea::simulation::AssemblyOutputType::SymmetricUpper
    );
    
    phases.push(SolvePhaseEntry {
        name: "Assembling stiffness matrix".to_string(),
        category: SolvePhaseCategory::Assembly,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: Some(format!("{} entries", stiffness.len())),
        start_offset_ms: start_offset,
    });
    
    // Phase: Matrix conversion
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    let n = simulation.nodes.len() * simulation.dofs;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    
    for ((i, j), v) in &stiffness {
        rows.push(*i);
        cols.push(*j);
        vals.push(*v);
        if i != j {
            rows.push(*j);
            cols.push(*i);
            vals.push(*v);
        }
    }
    
    phases.push(SolvePhaseEntry {
        name: "Converting to sparse format".to_string(),
        category: SolvePhaseCategory::Assembly,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: Some(format!("{}x{} matrix, {} entries", n, n, vals.len())),
        start_offset_ms: start_offset,
    });
    
    // Phase: Solve
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    let displacement = rust_fea::solver::direct_solve_triplet(
        n, &rows, &cols, &vals, &load_vector,
    ).map_err(|e| format!("Solver failed: {:?}", e))?;
    
    phases.push(SolvePhaseEntry {
        name: "Solving linear system (Ax=b)".to_string(),
        category: SolvePhaseCategory::Solve,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: Some(format!("{} DOFs", n)),
        start_offset_ms: start_offset,
    });
    
    // Phase: Update displacements
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        let dx = displacement[node_id * 3];
        let dy = displacement[node_id * 3 + 1];
        let dz = displacement[node_id * 3 + 2];
        node.set_displacement(dx, dy, dz);
    }
    
    phases.push(SolvePhaseEntry {
        name: "Updating node displacements".to_string(),
        category: SolvePhaseCategory::PostProcess,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: None,
        start_offset_ms: start_offset,
    });
    
    // Phase: Compute results
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    simulation.compute_result_fields();
    let (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain) = extract_stress_results(simulation);
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    phases.push(SolvePhaseEntry {
        name: "Computing stress/strain fields".to_string(),
        category: SolvePhaseCategory::PostProcess,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: None,
        start_offset_ms: start_offset,
    });
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
        nodal_von_mises,
        nodal_stress,
        nodal_strain,
        stats: crate::state::ResultStats {
            max_displacement: max_disp,
            min_displacement: min_disp,
            max_von_mises: max_vm,
            solver_time_ms: 0, // Will be set by caller
            phase_timing: None, // Will be set by caller
        },
        time_steps: Vec::new(),
    })
}

/// Explicit solver with timing (WASM)
#[cfg(target_arch = "wasm32")]
fn run_explicit_solver_impl_timed(
    simulation: &mut rust_fea::simulation::Simulation,
    settings: &crate::state::ExplicitSettings,
    phases: &mut Vec<crate::state::SolvePhaseEntry>,
    total_start: web_time::Instant,
) -> Result<SimulationResults, String> {
    use crate::state::{SolvePhaseEntry, SolvePhaseCategory};
    
    // Just wrap the non-timed impl and capture overall time
    let phase_start = web_time::Instant::now();
    let start_offset = (phase_start - total_start).as_millis() as u64;
    
    let result = run_explicit_solver_impl(simulation, settings, |_, _| {});
    
    phases.push(SolvePhaseEntry {
        name: format!("Explicit time integration ({} steps)", settings.time_steps),
        category: SolvePhaseCategory::Solve,
        duration_ms: phase_start.elapsed().as_millis() as u64,
        details: None,
        start_offset_ms: start_offset,
    });
    
    result
}

/// Core direct solver implementation
fn run_direct_solver_impl<F>(
    simulation: &mut rust_fea::simulation::Simulation,
    mut progress: F,
) -> Result<SimulationResults, String>
where
    F: FnMut(f32, String),
{
    progress(0.3, "Assembling stiffness matrix...".to_string());
    
    let (stiffness, load_vector) = simulation.assemble(
        rust_fea::simulation::AssemblyOutputType::SymmetricUpper
    );
    
    progress(0.5, "Converting to sparse matrix...".to_string());
    
    let n = simulation.nodes.len() * simulation.dofs;
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    
    for ((i, j), v) in &stiffness {
        rows.push(*i);
        cols.push(*j);
        vals.push(*v);
        if i != j {
            rows.push(*j);
            cols.push(*i);
            vals.push(*v);
        }
    }
    
    progress(0.7, "Solving linear system...".to_string());
    
    let displacement = rust_fea::solver::direct_solve_triplet(
        n, &rows, &cols, &vals, &load_vector,
    ).map_err(|e| format!("Solver failed: {:?}", e))?;
    
    progress(0.85, "Computing stress fields...".to_string());
    
    // Update node displacements
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        let dx = displacement[node_id * 3];
        let dy = displacement[node_id * 3 + 1];
        let dz = displacement[node_id * 3 + 2];
        node.set_displacement(dx, dy, dz);
    }
    
    simulation.compute_result_fields();
    
    // Extract stress, strain, and von Mises results - keyed by ELEMENT ID for visualization
    let (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain) = extract_stress_results(simulation);
    
    progress(0.95, "Processing results...".to_string());
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
        nodal_von_mises,
        nodal_stress,
        nodal_strain,
        stats: crate::state::ResultStats {
            max_displacement: max_disp,
            min_displacement: min_disp,
            max_von_mises: max_vm,
            solver_time_ms: 0,
            phase_timing: None,
        },
        time_steps: Vec::new(),
    })
}

/// Core explicit (time-stepping) solver implementation
fn run_explicit_solver_impl<F>(
    simulation: &mut rust_fea::simulation::Simulation,
    settings: &crate::state::ExplicitSettings,
    mut progress: F,
) -> Result<SimulationResults, String>
where
    F: FnMut(f32, String),
{
    use nalgebra::DVector;
    
    progress(0.1, "Computing element matrices...".to_string());
    
    // Compute stiffness and mass for all elements
    simulation.compute_all_element_stiffness();
    simulation.compute_all_element_mass();
    
    // Assemble global mass matrix (diagonal)
    let mass_diag = simulation.compute_global_mass_matrix_diagonal();
    
    // Get material from first element for wave speed calculation
    let active_el_ids = simulation.active_elements();
    let first_el_id = *active_el_ids.first()
        .ok_or("No active elements")?;
    let material = simulation.get_element(first_el_id)
        .ok_or("Cannot get element")?
        .get_material();
    
    let wave_speed = ((material.youngs_modulus / material.density) 
        * (1.0 - material.poisson_ratio) 
        / ((1.0 + material.poisson_ratio) * (1.0 - 2.0 * material.poisson_ratio))).sqrt();
    
    // Compute critical time step
    let dt_crit = simulation.mesh.compute_dt(wave_speed);
    // Use more conservative safety factor (0.5 instead of 0.9) for stability
    let dt = settings.time_step_override.unwrap_or(dt_crit * 0.5);
    
    progress(0.15, format!("dt={:.2e} s (crit={:.2e} s, c={:.0} m/s)", dt, dt_crit, wave_speed));
    
    let n = simulation.nodes.len() * simulation.dofs;
    let dofs = simulation.dofs;
    let total_steps = settings.time_steps;
    let save_interval = settings.vtk_save_steps.max(1);
    
    // Initialize state vectors (following library's approach)
    let mut u = DVector::zeros(n);
    let mut u_dot = DVector::zeros(n);
    let mut u_half_dot;
    
    // IMPORTANT: assemble_global_force populates fixed_global_nodal_values via handle_bc()
    // This MUST be called BEFORE get_specified_bc()
    simulation.assemble_global_force();
    let f_ext = simulation.load_vector.clone();
    
    // Get fixed BC DOFs and values (now populated)
    let bc_values = simulation.get_specified_bc();
    
    // Apply initial BCs to displacement
    for (dof, val) in &bc_values {
        if *dof < n {
            u[*dof] = *val;
        }
    }
    
    // Update node displacements from initial u vector
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        if node_id * dofs + 2 < n {
            node.set_displacement(u[node_id * dofs], u[node_id * dofs + 1], u[node_id * dofs + 2]);
        }
    }
    
    // Compute initial internal forces (reads from node displacements)
    let mut f_int = simulation.compute_force_vector(&u);
    
    let mut time_steps_data = Vec::new();
    let mut current_time = 0.0;
    
    progress(0.2, "Starting time integration...".to_string());
    
    // Explicit time integration (Velocity Verlet / Leapfrog)
    for step in 0..total_steps {
        // Compute residual and acceleration
        // a = (F_ext - F_int) / M
        let mut u_ddot = DVector::zeros(n);
        for i in 0..n {
            if mass_diag[i].abs() > 1e-30 {
                u_ddot[i] = (f_ext[i] - f_int[i]) / mass_diag[i];
            }
        }
        
        // Half-step velocity: v(t+dt/2) = v(t) + 0.5*dt*a(t)
        u_half_dot = &u_dot + 0.5 * dt * &u_ddot;
        
        // Update displacement: u(t+dt) = u(t) + dt*v(t+dt/2)
        u = &u + dt * &u_half_dot;
        
        // Apply boundary conditions to displacement
        for (dof, val) in &bc_values {
            if *dof < n {
                u[*dof] = *val;
                u_dot[*dof] = 0.0;
                u_half_dot[*dof] = 0.0;
            }
        }
        
        // CRITICAL: Update node displacements BEFORE computing forces
        for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
            if node_id * dofs + 2 < n {
                node.set_displacement(u[node_id * dofs], u[node_id * dofs + 1], u[node_id * dofs + 2]);
            }
        }
        
        // Compute new internal forces (now reads updated node displacements)
        f_int = simulation.compute_force_vector(&u);
        
        // Compute new acceleration
        let mut new_u_ddot = DVector::zeros(n);
        for i in 0..n {
            if mass_diag[i].abs() > 1e-30 {
                new_u_ddot[i] = (f_ext[i] - f_int[i]) / mass_diag[i];
            }
        }
        
        // Complete velocity update: v(t+dt) = v(t+dt/2) + 0.5*dt*a(t+dt)
        u_dot = &u_half_dot + 0.5 * dt * &new_u_ddot;
        
        // Apply velocity damping for stability
        u_dot *= 0.9995;
        
        // Apply BC to velocity
        for (dof, _) in &bc_values {
            if *dof < n {
                u_dot[*dof] = 0.0;
            }
        }
        
        current_time += dt;
        
        // Check for numerical instability (NaN or Inf)
        let max_u = u.iter().map(|x| x.abs()).fold(0.0f64, |a, b| a.max(b));
        if max_u.is_nan() || max_u.is_infinite() || max_u > 1e10 {
            return Err(format!(
                "Numerical instability detected at step {} (t={:.4e}s). \
                Try reducing time step or check boundary conditions. max_u={:.2e}", 
                step, current_time, max_u
            ));
        }
        
        // Record time step data at intervals
        if step % save_interval == 0 || step == total_steps - 1 {
            let max_disp = (0..n/dofs).map(|i| {
                let dx = u[i * dofs];
                let dy = u[i * dofs + 1];
                let dz = u[i * dofs + 2];
                (dx * dx + dy * dy + dz * dz).sqrt()
            }).fold(0.0f64, |a, b| a.max(b));
            
            // Kinetic energy: 0.5 * m * v^2
            let ke: f64 = (0..n).map(|i| 0.5 * mass_diag[i] * u_dot[i] * u_dot[i]).sum();
            
            time_steps_data.push(crate::state::TimeStepResult {
                time: current_time,
                iteration: step,
                max_displacement: max_disp,
                kinetic_energy: ke,
            });
        }
        
        // Progress update
        if step % (total_steps / 20).max(1) == 0 {
            let pct = 0.2 + 0.7 * (step as f32 / total_steps as f32);
            progress(pct, format!("Step {}/{} (t={:.4e} s, max_u={:.2e})", 
                step, total_steps, current_time, 
                u.iter().map(|x| x.abs()).fold(0.0f64, |a,b| a.max(b))));
        }
    }
    
    progress(0.9, "Computing final stress fields...".to_string());
    
    // Update node displacements for stress computation
    let displacement: Vec<f64> = u.iter().cloned().collect();
    for (node_id, node) in simulation.nodes.iter_mut().enumerate() {
        if node_id * 3 + 2 < displacement.len() {
            let dx = displacement[node_id * 3];
            let dy = displacement[node_id * 3 + 1];
            let dz = displacement[node_id * 3 + 2];
            node.set_displacement(dx, dy, dz);
        }
    }
    
    simulation.compute_result_fields();
    
    // Extract stress, strain, and von Mises results
    let (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain) = extract_stress_results(simulation);
    
    progress(0.95, "Processing results...".to_string());
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
        nodal_von_mises,
        nodal_stress,
        nodal_strain,
        stats: crate::state::ResultStats {
            max_displacement: max_disp,
            min_displacement: min_disp,
            max_von_mises: max_vm,
            solver_time_ms: 0,
            phase_timing: None,
        },
        time_steps: time_steps_data,
    })
}

/// Extract stress, strain, and von Mises results from simulation, keyed by element ID
/// Also extracts nodal values for smooth field interpolation
fn extract_stress_results(
    simulation: &rust_fea::simulation::Simulation,
) -> (
    std::collections::HashMap<usize, Vec<f64>>,  // stresses (per element)
    std::collections::HashMap<usize, Vec<f64>>,  // strains (per element)
    std::collections::HashMap<usize, f64>,       // von_mises (per element)
    f64,                                         // max_vm
    Vec<f64>,                                    // nodal_von_mises
    Vec<[f64; 6]>,                               // nodal_stress
    Vec<[f64; 6]>,                               // nodal_strain
) {
    let mut stresses = std::collections::HashMap::new();
    let mut strains = std::collections::HashMap::new();
    let mut von_mises = std::collections::HashMap::new();
    let mut max_vm = 0.0f64;
    
    // Library uses s_xx, s_yy, etc. for stress fields and e_xx, e_yy, etc. for strain
    let stress_fields = ["s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"];
    let strain_fields = ["e_xx", "e_yy", "e_zz", "e_xy", "e_yz", "e_xz"];
    
    // Extract nodal field values
    let num_nodes = simulation.mesh.nodes.len();
    let mut nodal_von_mises = vec![0.0; num_nodes];
    let mut nodal_stress = vec![[0.0; 6]; num_nodes];
    let mut nodal_strain = vec![[0.0; 6]; num_nodes];
    
    // Extract nodal von Mises
    if let Some(vm_field) = simulation.node_fields.get("vm") {
        for (i, &val) in vm_field.iter().enumerate() {
            if i < nodal_von_mises.len() {
                nodal_von_mises[i] = val;
                max_vm = max_vm.max(val);
            }
        }
    }
    
    // Extract nodal stress components
    for (comp_idx, field_name) in stress_fields.iter().enumerate() {
        if let Some(field) = simulation.node_fields.get(*field_name) {
            for (nid, &val) in field.iter().enumerate() {
                if nid < nodal_stress.len() {
                    nodal_stress[nid][comp_idx] = val;
                }
            }
        }
    }
    
    // Extract nodal strain components
    for (comp_idx, field_name) in strain_fields.iter().enumerate() {
        if let Some(field) = simulation.node_fields.get(*field_name) {
            for (nid, &val) in field.iter().enumerate() {
                if nid < nodal_strain.len() {
                    nodal_strain[nid][comp_idx] = val;
                }
            }
        }
    }
    
    // Also compute per-element averages for backwards compatibility
    for el_id in simulation.active_elements() {
        if let Some(element) = simulation.get_element(el_id) {
            let conn = element.get_connectivity();
            let mut elem_stress = vec![0.0; 6];
            let mut elem_strain = vec![0.0; 6];
            let mut elem_vm = 0.0;
            let mut count = 0;
            
            for &nid in conn {
                // Accumulate stress components
                for (i, field_name) in stress_fields.iter().enumerate() {
                    if let Some(field) = simulation.node_fields.get(*field_name) {
                        if nid < field.len() {
                            elem_stress[i] += field[nid];
                        }
                    }
                }
                
                // Accumulate strain components
                for (i, field_name) in strain_fields.iter().enumerate() {
                    if let Some(field) = simulation.node_fields.get(*field_name) {
                        if nid < field.len() {
                            elem_strain[i] += field[nid];
                        }
                    }
                }
                
                // Accumulate von Mises
                if let Some(vm_field) = simulation.node_fields.get("vm") {
                    if nid < vm_field.len() {
                        elem_vm += vm_field[nid];
                    }
                }
                count += 1;
            }
            
            // Average over element nodes
            if count > 0 {
                for v in elem_stress.iter_mut() {
                    *v /= count as f64;
                }
                for v in elem_strain.iter_mut() {
                    *v /= count as f64;
                }
                elem_vm /= count as f64;
            }
            
            stresses.insert(el_id, elem_stress);
            strains.insert(el_id, elem_strain);
            von_mises.insert(el_id, elem_vm);
        }
    }
    
    (stresses, strains, von_mises, max_vm, nodal_von_mises, nodal_stress, nodal_strain)
}

/// Compute displacement statistics
fn compute_displacement_stats(displacement: &[f64]) -> (f64, f64) {
    let mut max_disp = 0.0f64;
    let mut min_disp = f64::MAX;
    
    for i in 0..displacement.len() / 3 {
        let dx = displacement[i * 3];
        let dy = displacement[i * 3 + 1];
        let dz = displacement[i * 3 + 2];
        let mag = (dx * dx + dy * dy + dz * dz).sqrt();
        max_disp = max_disp.max(mag);
        min_disp = min_disp.min(mag);
    }
    
    if min_disp == f64::MAX {
        min_disp = 0.0;
    }
    
    (max_disp, min_disp)
}


/// Parse Abaqus INP format from a string (for WASM file uploads)
/// This is a simplified parser that handles basic INP files exported from Gmsh
fn parse_inp_from_string(content: &str) -> Result<rust_fea::mesh::MeshAssembly, String> {
    use std::collections::HashMap;
    use rust_fea::mesh::{MeshAssembly, MeshElement, MeshNode, NodeGroup, ElementGroup};
    
    let mut mesh = MeshAssembly::empty();
    let lines: Vec<&str> = content.lines().collect();
    
    #[derive(PartialEq)]
    enum Block { None, Node, Element, Elset, Nset }
    let mut current_block = Block::None;
    let mut current_params: HashMap<String, String> = HashMap::new();
    let mut temp_nset_nodes: Vec<usize> = Vec::new();
    let mut temp_elset_elements: Vec<usize> = Vec::new();
    
    for line in lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with("**") {
            continue;
        }
        
        // Check for section headers
        if line.starts_with('*') {
            // Parse previous block's data if needed
            if current_block == Block::Nset {
                if let Some(name) = current_params.get("NSET") {
                    if !temp_nset_nodes.is_empty() {
                        mesh.node_groups.insert(name.clone(), NodeGroup { 
                            name: name.clone(), 
                            nodes: temp_nset_nodes.clone() 
                        });
                    }
                }
                temp_nset_nodes.clear();
            }
            if current_block == Block::Elset {
                if let Some(name) = current_params.get("ELSET") {
                    if !temp_elset_elements.is_empty() {
                        let el_type = current_params.get("TYPE")
                            .cloned()
                            .unwrap_or_else(|| "C3D8".to_string());
                        mesh.element_groups.insert(name.clone(), ElementGroup {
                            name: name.clone(),
                            elements: temp_elset_elements.clone(),
                            el_type,
                        });
                    }
                }
                temp_elset_elements.clear();
            }
            
            // Parse new header
            current_params.clear();
            let upper = line.to_uppercase();
            
            if upper.starts_with("*NODE") {
                current_block = Block::Node;
            } else if upper.starts_with("*ELEMENT") {
                current_block = Block::Element;
            } else if upper.starts_with("*ELSET") {
                current_block = Block::Elset;
            } else if upper.starts_with("*NSET") {
                current_block = Block::Nset;
            } else {
                current_block = Block::None;
            }
            
            // Parse parameters from header line
            for part in line.split(',').skip(1) {
                let part = part.trim();
                if let Some(idx) = part.find('=') {
                    let key = part[..idx].trim().to_uppercase();
                    let value = part[idx + 1..].trim().to_string();
                    current_params.insert(key, value);
                }
            }
            continue;
        }
        
        // Parse data lines
        match current_block {
            Block::Node => {
                // Format: id, x, y, z
                let parts: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
                if parts.len() >= 4 {
                    if let (Ok(id), Ok(x), Ok(y), Ok(z)) = (
                        parts[0].parse::<usize>(),
                        parts[1].parse::<f64>(),
                        parts[2].parse::<f64>(),
                        parts[3].parse::<f64>(),
                    ) {
                        mesh.nodes.insert(id, MeshNode {
                            id,
                            coordinates: vec![x, y, z],
                        });
                    }
                }
            }
            Block::Element => {
                // Format: id, n1, n2, n3, n4, ... (may span multiple lines)
                let parts: Vec<&str> = line.split(',').map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
                if parts.len() >= 2 {
                    if let Ok(id) = parts[0].parse::<usize>() {
                        let nodes: Vec<usize> = parts[1..]
                            .iter()
                            .filter_map(|s| s.parse::<usize>().ok())
                            .collect();
                        
                        if !nodes.is_empty() {
                            let el_type = current_params.get("TYPE")
                                .map(|s| s.to_uppercase())
                                .unwrap_or_else(|| "C3D8".to_string());
                            mesh.elements.insert(id, MeshElement {
                                id,
                                connectivity: nodes,
                                el_type,
                                name: String::new(),
                            });
                        }
                    }
                }
            }
            Block::Nset => {
                // Node set: comma-separated node IDs
                for part in line.split(',') {
                    let part = part.trim();
                    if !part.is_empty() {
                        if let Ok(id) = part.parse::<usize>() {
                            temp_nset_nodes.push(id);
                        }
                    }
                }
            }
            Block::Elset => {
                // Element set: comma-separated element IDs
                for part in line.split(',') {
                    let part = part.trim();
                    if !part.is_empty() {
                        if let Ok(id) = part.parse::<usize>() {
                            temp_elset_elements.push(id);
                        }
                    }
                }
            }
            Block::None => {}
        }
    }
    
    // Handle any remaining sets
    if current_block == Block::Nset {
        if let Some(name) = current_params.get("NSET") {
            if !temp_nset_nodes.is_empty() {
                mesh.node_groups.insert(name.clone(), NodeGroup { 
                    name: name.clone(), 
                    nodes: temp_nset_nodes 
                });
            }
        }
    }
    if current_block == Block::Elset {
        if let Some(name) = current_params.get("ELSET") {
            if !temp_elset_elements.is_empty() {
                let el_type = current_params.get("TYPE")
                    .cloned()
                    .unwrap_or_else(|| "C3D8".to_string());
                mesh.element_groups.insert(name.clone(), ElementGroup {
                    name: name.clone(),
                    elements: temp_elset_elements,
                    el_type,
                });
            }
        }
    }
    
    if mesh.nodes.is_empty() {
        return Err("No nodes found in INP file".to_string());
    }
    if mesh.elements.is_empty() {
        return Err("No elements found in INP file".to_string());
    }
    
    Ok(mesh)
}
