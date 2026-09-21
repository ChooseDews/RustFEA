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
    
    /// Test mode for automated screenshots
    #[cfg(not(target_arch = "wasm32"))]
    test_mode: crate::state::TestMode,
    #[cfg(not(target_arch = "wasm32"))]
    test_frame_count: u32,
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
            // Run synchronously on WASM
            let result = run_simulation_sync(simulation, solver_type, explicit_settings);
            match result {
                Ok(results) => {
                    self.state.is_running = false;
                    self.state.progress = 1.0;
                    self.state.results = Some(results);
                    self.state.status_message = "Simulation completed!".to_string();
                    self.state.ui_state.active_panel = ActivePanel::Results;
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
        .default_width(400.0)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
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
                    let mut color = app.state.ui_state.display_settings.background_color;
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
            });
            
            ui.add_space(8.0);
            
            // Other settings
            ui.group(|ui| {
                ui.label("Viewport");
                ui.checkbox(&mut app.state.ui_state.display_settings.show_axis, "Show Axis Indicator");
            });
            
            ui.add_space(16.0);
            
            ui.horizontal(|ui| {
                if ui.button("Reset to Defaults").clicked() {
                    app.state.ui_state.display_settings = crate::state::DisplaySettings::default();
                }
                
                if ui.button("Close").clicked() {
                    app.state.ui_state.preferences_dialog_open = false;
                }
            });
        });
}

fn setup_custom_style(ctx: &egui::Context) {
    use egui::{Color32, FontId, FontFamily, CornerRadius, Stroke, Shadow, Vec2, Margin};
    use egui::style::{Widgets, WidgetVisuals, Selection, HandleShape};
    
    // =========================================================================
    // Custom Fonts - Inter (UI) + JetBrains Mono (code) + Noto Sans Symbols 2 (icons)
    // =========================================================================
    let mut fonts = egui::FontDefinitions::default();
    
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
    
    // Set Inter as the primary proportional font with Noto Symbols as fallback
    fonts.families.entry(FontFamily::Proportional).or_default()
        .insert(0, "inter_regular".to_owned());
    fonts.families.entry(FontFamily::Proportional).or_default()
        .push("noto_symbols".to_owned());
    
    // Set JetBrains Mono as the primary monospace font with Noto Symbols fallback
    fonts.families.entry(FontFamily::Monospace).or_default()
        .insert(0, "jetbrains_mono".to_owned());
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
    visuals.window_stroke = Stroke::new(1.0, stroke_subtle);
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
        stroke: Stroke::new(1.0, accent),
    };
    
    // =========================================================================
    // Widget Visuals - State-specific styling
    // =========================================================================
    visuals.widgets = Widgets {
        // Non-interactive widgets (labels, panel backgrounds)
        noninteractive: WidgetVisuals {
            bg_fill: bg_medium,
            weak_bg_fill: bg_light,
            bg_stroke: Stroke::new(1.0, stroke_subtle),
            corner_radius: CornerRadius::same(4),
            fg_stroke: Stroke::new(1.0, text_secondary),
            expansion: 0.0,
        },
        // Interactive widgets at rest
        inactive: WidgetVisuals {
            bg_fill: bg_light,
            weak_bg_fill: bg_light,
            bg_stroke: Stroke::new(1.0, stroke_subtle),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.0, text_primary),
            expansion: 0.0,
        },
        // Hovered interactive widgets
        hovered: WidgetVisuals {
            bg_fill: bg_hover,
            weak_bg_fill: bg_hover,
            bg_stroke: Stroke::new(1.5, accent),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.5, text_primary),
            expansion: 1.0,
        },
        // Active (clicked/focused) widgets
        active: WidgetVisuals {
            bg_fill: bg_active,
            weak_bg_fill: bg_active,
            bg_stroke: Stroke::new(2.0, accent_hover),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(2.0, text_primary),
            expansion: 1.0,
        },
        // Open (e.g., combo box with menu open)
        open: WidgetVisuals {
            bg_fill: bg_hover,
            weak_bg_fill: bg_hover,
            bg_stroke: Stroke::new(1.5, accent),
            corner_radius: CornerRadius::same(6),
            fg_stroke: Stroke::new(1.5, text_primary),
            expansion: 1.0,
        },
    };
    
    // Text cursor styling
    visuals.text_cursor.stroke = Stroke::new(2.0, accent);
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
    start: std::time::Instant,
    phase_start: std::time::Instant,
    completed_phases: Vec<SolvePhaseEntry>,
    current_phase: Option<String>,
    current_category: Option<SolvePhaseCategory>,
    tx: Sender<SimMessage>,
}

#[cfg(not(target_arch = "wasm32"))]
impl PhaseTracker {
    fn new(tx: Sender<SimMessage>) -> Self {
        let now = std::time::Instant::now();
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
        self.phase_start = std::time::Instant::now();
        
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
    
    let (stresses, strains, von_mises, max_vm) = extract_stress_results(simulation);
    
    tracker.begin_phase("Processing results", SolvePhaseCategory::PostProcess);
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
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
    let mut u_half_dot = DVector::zeros(n);
    
    // Get fixed BC DOFs and values
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
    simulation.assemble_global_force();
    let f_ext = simulation.load_vector.clone();
    
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
    
    let (stresses, strains, von_mises, max_vm) = extract_stress_results(simulation);
    
    tracker.begin_phase("Processing results", SolvePhaseCategory::PostProcess);
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
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

/// Run simulation synchronously (WASM)
#[cfg(target_arch = "wasm32")]
fn run_simulation_sync(
    mut simulation: rust_fea::simulation::Simulation,
    solver_type: crate::state::SolverType,
    explicit_settings: crate::state::ExplicitSettings,
) -> Result<SimulationResults, String> {
    simulation.initialize();
    simulation.one_time_init();
    
    match solver_type {
        crate::state::SolverType::Direct => run_direct_solver_impl(&mut simulation, |_, _| {}),
        crate::state::SolverType::Explicit => run_explicit_solver_impl(&mut simulation, &explicit_settings, |_, _| {}),
    }
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
    let (stresses, strains, von_mises, max_vm) = extract_stress_results(simulation);
    
    progress(0.95, "Processing results...".to_string());
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
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
    let mut u_half_dot = DVector::zeros(n);
    
    // Get fixed BC DOFs and values
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
    
    // Assemble external force vector
    simulation.assemble_global_force();
    let f_ext = simulation.load_vector.clone();
    
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
    let (stresses, strains, von_mises, max_vm) = extract_stress_results(simulation);
    
    progress(0.95, "Processing results...".to_string());
    
    let (max_disp, min_disp) = compute_displacement_stats(&displacement);
    
    Ok(SimulationResults {
        displacements: displacement,
        stresses,
        strains,
        von_mises,
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
fn extract_stress_results(
    simulation: &rust_fea::simulation::Simulation,
) -> (
    std::collections::HashMap<usize, Vec<f64>>,  // stresses
    std::collections::HashMap<usize, Vec<f64>>,  // strains
    std::collections::HashMap<usize, f64>,       // von_mises
    f64                                          // max_vm
) {
    let mut stresses = std::collections::HashMap::new();
    let mut strains = std::collections::HashMap::new();
    let mut von_mises = std::collections::HashMap::new();
    let mut max_vm = 0.0f64;
    
    // Library uses s_xx, s_yy, etc. for stress fields and e_xx, e_yy, etc. for strain
    let stress_fields = ["s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"];
    let strain_fields = ["e_xx", "e_yy", "e_zz", "e_xy", "e_yz", "e_xz"];
    
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
                
                // Accumulate von Mises - library uses "vm" field name
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
            max_vm = max_vm.max(elem_vm);
        }
    }
    
    (stresses, strains, von_mises, max_vm)
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
