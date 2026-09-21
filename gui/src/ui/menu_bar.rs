//! Top menu bar for RustFEA GUI

use eframe::egui;
use crate::app::FeaApp;
use crate::examples::{ExampleType, load_example, load_example_with_config, MeshResolution};
use crate::project_io;

pub fn show(ctx: &egui::Context, app: &mut FeaApp) {
    egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
        egui::menu::bar(ui, |ui| {
            // File menu
            ui.menu_button("File", |ui| {
                if ui.button("📂 Open Project...").clicked() {
                    open_project_dialog(app);
                    ui.close_menu();
                }
                
                if ui.button("💾 Save Project").clicked() {
                    save_project(app);
                    ui.close_menu();
                }
                
                if ui.button("💾 Save Project As...").clicked() {
                    save_project_as(app);
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("📥 Import Mesh...").clicked() {
                    import_mesh_dialog(app);
                    ui.close_menu();
                }
                
                // Recent files submenu
                #[cfg(not(target_arch = "wasm32"))]
                if !app.state.ui_state.recent_files.files.is_empty() {
                    ui.menu_button("📋 Recent Files", |ui| {
                        let mut file_to_load: Option<std::path::PathBuf> = None;
                        
                        for recent in &app.state.ui_state.recent_files.files {
                            let label = format!("{}", recent.name);
                            if ui.button(&label).on_hover_text(recent.path.to_string_lossy().to_string()).clicked() {
                                file_to_load = Some(recent.path.clone());
                                ui.close_menu();
                            }
                        }
                        
                        if let Some(path) = file_to_load {
                            app.load_mesh_file(path);
                        }
                    });
                }
                
                // Load Example submenu
                ui.menu_button("📚 Load Example", |ui| {
                    // Cantilever Beam
                    ui.horizontal(|ui| {
                        if ui.button("🔩 Cantilever Beam").clicked() {
                            app.state.ui_state.example_config.example_type = crate::examples::ExampleType::CantileverBeam;
                            app.state.ui_state.example_dialog_open = true;
                            ui.close_menu();
                        }
                    });
                    ui.label("  Fixed beam with tip load (bending)");
                    
                    ui.add_space(4.0);
                    
                    // Torque Shaft
                    ui.horizontal(|ui| {
                        if ui.button("🔧 Torque Shaft").clicked() {
                            app.state.ui_state.example_config.example_type = crate::examples::ExampleType::TorqueShaft;
                            app.state.ui_state.example_dialog_open = true;
                            ui.close_menu();
                        }
                    });
                    ui.label("  Cylinder with applied torque (torsion)");
                    
                    ui.add_space(4.0);
                    
                    // Contact Blocks
                    ui.horizontal(|ui| {
                        if ui.button("📦 Contact Blocks").clicked() {
                            app.state.ui_state.example_config.example_type = crate::examples::ExampleType::ContactBlocks;
                            app.state.ui_state.example_dialog_open = true;
                            ui.close_menu();
                        }
                    });
                    ui.label("  Two blocks with contact (explicit solver)");
                    
                    ui.separator();
                    ui.label("Click to customize mesh & loads");
                });
                
                ui.separator();
                
                if ui.button("📤 Export Results...").clicked() {
                    export_results_dialog(app);
                    ui.close_menu();
                }
                
                #[cfg(not(target_arch = "wasm32"))]
                if ui.button("📷 Save Screenshot...").clicked() {
                    app.state.ui_state.screenshot_dialog_open = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                #[cfg(not(target_arch = "wasm32"))]
                if ui.button("🚪 Exit").clicked() {
                    std::process::exit(0);
                }
            });
            
            // Edit menu
            ui.menu_button("Edit", |ui| {
                if ui.button("⟲ Undo").on_hover_text("Coming soon").clicked() {
                    app.state.status_message = "Undo not yet implemented".to_string();
                    ui.close_menu();
                }
                
                if ui.button("⟳ Redo").on_hover_text("Coming soon").clicked() {
                    app.state.status_message = "Redo not yet implemented".to_string();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("⚙ Preferences...").clicked() {
                    app.state.ui_state.preferences_dialog_open = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("🗑 Clear Workspace").on_hover_text("Remove all meshes and reset to default state").clicked() {
                    // Reset state to default while preserving UI preferences
                    let display_settings = app.state.ui_state.display_settings.clone();
                    app.state = crate::state::AppState::new();
                    app.state.ui_state.display_settings = display_settings;
                    app.renderer = None;
                    app.render_cache.invalidate();
                    app.state.status_message = "Workspace cleared".to_string();
                    ui.close_menu();
                }
            });
            
            // View menu
            ui.menu_button("View", |ui| {
                if ui.checkbox(&mut app.state.ui_state.show_faces, "Show Faces").on_hover_text("(F)").changed() {
                    app.renderer = None;
                }
                if ui.checkbox(&mut app.state.ui_state.show_wireframe, "Show Wireframe").on_hover_text("(W)").changed() {
                    app.renderer = None;
                }
                if ui.checkbox(&mut app.state.ui_state.show_nodes, "Show Nodes").on_hover_text("(N)").changed() {
                    app.renderer = None;
                }
                ui.checkbox(&mut app.state.ui_state.show_node_groups, "Show Node Groups");
                ui.checkbox(&mut app.state.ui_state.show_boundary_conditions, "Show Boundary Conditions").on_hover_text("(B)");
                
                ui.separator();
                
                // Clipping plane
                if ui.checkbox(&mut app.state.ui_state.clipping_plane.enabled, "✂ Clipping Plane")
                    .on_hover_text("Section view (C)")
                    .changed() 
                {
                    app.render_cache.invalidate();
                }
                
                ui.separator();
                
                if ui.button("🎯 Reset Camera").on_hover_text("(Home)").clicked() {
                    if let Some(mesh) = app.state.current_mesh() {
                        let bounds = mesh.bounds;
                        app.state.ui_state.camera.fit_to_bounds(&bounds);
                    }
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("📐 Front View").on_hover_text("(1)").clicked() {
                    app.state.ui_state.camera.set_front_view();
                    ui.close_menu();
                }
                
                if ui.button("📐 Top View").on_hover_text("(3)").clicked() {
                    app.state.ui_state.camera.set_top_view();
                    ui.close_menu();
                }
                
                if ui.button("📐 Side View").on_hover_text("(2)").clicked() {
                    app.state.ui_state.camera.set_side_view();
                    ui.close_menu();
                }
                
                if ui.button("📐 Isometric").on_hover_text("(0)").clicked() {
                    app.state.ui_state.camera.set_iso_view();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("⌨ Keyboard Shortcuts...").on_hover_text("(?)").clicked() {
                    app.state.ui_state.show_shortcuts_help = true;
                    ui.close_menu();
                }
            });
            
            // Simulation menu
            ui.menu_button("Simulation", |ui| {
                let can_run = app.state.current_mesh().is_some() && !app.state.is_running;
                
                ui.add_enabled_ui(can_run, |ui| {
                    if ui.button("▶ Run Simulation").clicked() {
                        app.start_simulation();
                        ui.close_menu();
                    }
                });
                
                ui.add_enabled_ui(app.state.is_running, |ui| {
                    if ui.button("⏹ Stop Simulation").clicked() {
                        app.stop_simulation();
                        ui.close_menu();
                    }
                });
                
                ui.separator();
                
                if ui.button("🔧 Solver Settings...").clicked() {
                    // TODO: Open solver settings dialog
                    ui.close_menu();
                }
            });
            
            // Help menu
            ui.menu_button("Help", |ui| {
                if ui.button("⌨ Keyboard Shortcuts").clicked() {
                    app.state.ui_state.show_shortcuts_help = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("📖 Documentation").clicked() {
                    // Open browser to documentation (can add URL open later)
                    app.state.status_message = "Documentation: https://github.com/your-repo/rust_fea".to_string();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("ℹ About RustFEA").clicked() {
                    app.state.ui_state.about_dialog_open = true;
                    ui.close_menu();
                }
            });
            
            // Spacer to push Run button to the right
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                // Run Simulation button - prominent in top bar
                let can_run = app.state.current_mesh().is_some() && !app.state.is_running;
                
                if app.state.is_running {
                    // Show stop button when running
                    if ui.button("⏹ Stop").clicked() {
                        app.stop_simulation();
                    }
                    // Show progress
                    ui.spinner();
                    ui.label(&app.state.status_message);
                } else {
                    ui.add_enabled_ui(can_run, |ui| {
                        let run_btn = egui::Button::new("▶ Run Simulation")
                            .fill(egui::Color32::from_rgb(46, 125, 50)); // Green
                        if ui.add(run_btn).clicked() {
                            app.start_simulation();
                        }
                    });
                    
                    if !can_run && app.state.current_mesh().is_none() {
                        ui.label("Load a mesh to run");
                    }
                }
            });
        });
    });
}

fn open_project_dialog(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("RustFEA Files", &["rfea", "toml"])
            .add_filter("TOML Project", &["toml"])
            .add_filter("Project Bundle", &["rfea"])
            .pick_file()
        {
            load_project(app, path);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        // Use browser file picker for bundles (includes both .rfea and .toml)
        crate::web_file_io::open_bundle_file_picker();
        app.state.status_message = "Select a project file...".to_string();
    }
}

#[cfg(feature = "native")]
fn load_project(app: &mut FeaApp, path: std::path::PathBuf) {
    app.state.status_message = format!("Loading project from {:?}...", path);
    
    match project_io::load_project(&path) {
        Ok((project, config)) => {
            // Load the simulation config
            app.state.simulation_config = config;
            
            // Try to load the mesh if specified
            if let Some(mesh_path) = &project.mesh {
                // Resolve relative paths relative to project file
                let mesh_path = if std::path::Path::new(mesh_path).is_relative() {
                    path.parent()
                        .map(|p| p.join(mesh_path))
                        .unwrap_or_else(|| std::path::PathBuf::from(mesh_path))
                } else {
                    std::path::PathBuf::from(mesh_path)
                };
                
                if mesh_path.exists() {
                    app.load_mesh_file(mesh_path);
                } else {
                    app.state.status_message = format!(
                        "Project loaded, but mesh not found: {}",
                        mesh_path.display()
                    );
                }
            }
            
            app.state.project_path = Some(path);
            app.renderer = None;
            
            // Fit camera if mesh was loaded
            if let Some(mesh_state) = app.state.current_mesh() {
                let bounds = mesh_state.bounds;
                app.state.ui_state.camera.fit_to_bounds(&bounds);
                app.state.status_message = format!("Project loaded: {}", project.name);
            } else if project.mesh.is_none() {
                app.state.status_message = format!(
                    "Project loaded: {} (no mesh specified)",
                    project.name
                );
            }
            
            // Switch to setup panel
            app.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
        }
        Err(e) => {
            app.state.status_message = format!("Failed to load project: {}", e);
        }
    }
}

fn save_project(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = &app.state.project_path.clone() {
            do_save_project(app, &path);
        } else {
            save_project_as(app);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        // In browser, always download with a generated name
        download_project_as(app, "project.toml");
    }
}

fn save_project_as(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("TOML Project", &["toml"])
            .save_file()
        {
            do_save_project(app, &path);
            app.state.project_path = Some(path);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        // Generate filename from mesh name or use default
        let filename = app.state.current_mesh()
            .map(|m| format!("{}.toml", m.name.replace(" ", "_")))
            .unwrap_or_else(|| "project.toml".to_string());
        download_project_as(app, &filename);
    }
}

#[cfg(not(feature = "native"))]
fn download_project_as(app: &mut FeaApp, filename: &str) {
    let mesh_name = app.state.current_mesh()
        .map(|m| m.name.clone())
        .unwrap_or_else(|| "Untitled".to_string());
    
    // Serialize the project TOML
    let toml_str = match project_io::serialize_project_toml(&mesh_name, None, &app.state.simulation_config) {
        Ok(s) => s,
        Err(e) => {
            app.state.status_message = format!("Failed to serialize project: {}", e);
            return;
        }
    };
    
    // Serialize mesh to JSON if available
    let mesh_json = app.state.current_mesh().map(|m| {
        serde_json::to_string_pretty(&m.mesh).ok()
    }).flatten();
    
    // Create and download ZIP bundle
    let project_name = filename.trim_end_matches(".toml").trim_end_matches(".rfea");
    crate::web_file_io::download_project_bundle(
        project_name,
        &toml_str,
        mesh_json.as_deref(),
    );
    app.state.status_message = format!("Project bundle downloaded as {}.rfea", project_name);
}

#[cfg(feature = "native")]
fn do_save_project(app: &mut FeaApp, path: &std::path::Path) {
    let mesh_path = app.state.current_mesh().and_then(|m| m.path.as_ref().map(|p| p.as_path()));
    let name = app.state.current_mesh()
        .map(|m| m.name.as_str())
        .unwrap_or("Untitled");
    
    match project_io::save_project(path, name, mesh_path, &app.state.simulation_config) {
        Ok(()) => {
            app.state.status_message = format!("Project saved to {}", path.display());
        }
        Err(e) => {
            app.state.status_message = format!("Failed to save: {}", e);
        }
    }
}

fn import_mesh_dialog(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("Mesh Files", &["inp", "bin", "json", "xz"])
            .add_filter("Abaqus INP", &["inp"])
            .add_filter("Binary Mesh", &["bin", "bin.xz"])
            .add_filter("JSON Mesh", &["json", "json.xz"])
            .pick_file()
        {
            app.load_mesh_file(path);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        // Use browser file picker
        crate::web_file_io::open_mesh_file_picker();
        app.state.status_message = "Select a mesh file...".to_string();
    }
}

fn export_results_dialog(app: &mut FeaApp) {
    if app.state.results.is_none() {
        app.state.status_message = "No results to export".to_string();
        return;
    }
    
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("VTK File", &["vtk"])
            .add_filter("JSON Results", &["json"])
            .save_file()
        {
            // TODO: Export results
            app.state.status_message = format!("Results exported to {:?}", path);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        // Export results as a simple JSON summary
        if let Some(results) = &app.state.results {
            let summary = format!(
                r#"{{
  "max_displacement": {},
  "min_displacement": {},
  "max_von_mises": {},
  "solver_time_ms": {},
  "num_nodes": {},
  "num_time_steps": {}
}}"#,
                results.stats.max_displacement,
                results.stats.min_displacement,
                results.stats.max_von_mises,
                results.stats.solver_time_ms,
                results.displacements.len() / 3,
                results.time_steps.len()
            );
            crate::web_file_io::download_file("results_summary.json", summary.as_bytes(), "application/json");
            app.state.status_message = "Results summary downloaded".to_string();
        }
    }
}

/// Load a built-in example into the application
fn load_example_into_app(app: &mut FeaApp, example_type: ExampleType) {
    let example = load_example(example_type);
    
    // Clear previous results
    app.state.results = None;
    
    // Add the mesh
    app.state.add_mesh(example.mesh, example.name.clone(), None);
    
    // Set up boundary conditions
    app.state.simulation_config.boundary_conditions = example.boundary_conditions;
    
    // Set solver type
    app.state.simulation_config.solver = example.solver_type;
    
    // Fit camera to new mesh
    if let Some(mesh) = app.state.current_mesh() {
        let bounds = mesh.bounds;
        app.state.ui_state.camera.fit_to_bounds(&bounds);
    }
    
    // Invalidate renderer
    app.renderer = None;
    
    // Switch to Setup panel so user can see the BCs
    app.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
    
    // Status message with description
    app.state.status_message = format!("Loaded: {} - {}", example.name, example.description);
}

/// Load example with custom configuration
fn load_example_with_custom_config(app: &mut FeaApp) {
    let config = &app.state.ui_state.example_config;
    let example = load_example_with_config(config);
    
    // Clear previous results
    app.state.results = None;
    
    // Add the mesh
    app.state.add_mesh(example.mesh, example.name.clone(), None);
    
    // Set up boundary conditions
    app.state.simulation_config.boundary_conditions = example.boundary_conditions;
    
    // Set solver type
    app.state.simulation_config.solver = example.solver_type;
    
    // Fit camera to new mesh
    if let Some(mesh) = app.state.current_mesh() {
        let bounds = mesh.bounds;
        app.state.ui_state.camera.fit_to_bounds(&bounds);
    }
    
    // Invalidate renderer
    app.renderer = None;
    
    // Switch to Setup panel
    app.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
    
    // Status message
    app.state.status_message = format!("Loaded: {}", example.description);
}

/// Show the example configuration dialog
pub fn show_example_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    if !app.state.ui_state.example_dialog_open {
        return;
    }
    
    let mut open = app.state.ui_state.example_dialog_open;
    
    egui::Window::new("📚 Load Example")
        .open(&mut open)
        .resizable(false)
        .collapsible(false)
        .default_width(400.0)
        .show(ctx, |ui| {
            let config = &mut app.state.ui_state.example_config;
            
            // Example type selection
            ui.heading("Select Example");
            ui.horizontal(|ui| {
                let prev_type = config.example_type;
                if ui.selectable_value(&mut config.example_type, ExampleType::CantileverBeam, "🔧 Cantilever Beam").clicked() 
                    && prev_type != ExampleType::CantileverBeam {
                    // Reset to default params for this example
                    config.load_magnitude = 10000.0;
                }
                if ui.selectable_value(&mut config.example_type, ExampleType::TorqueShaft, "⚙ Torque Shaft").clicked()
                    && prev_type != ExampleType::TorqueShaft {
                    config.load_magnitude = 5000.0;
                }
                if ui.selectable_value(&mut config.example_type, ExampleType::ContactBlocks, "📦 Contact Blocks").clicked()
                    && prev_type != ExampleType::ContactBlocks {
                    config.load_magnitude = 50000.0;
                }
            });
            ui.label(config.example_type.description());
            
            ui.add_space(12.0);
            ui.separator();
            ui.add_space(8.0);
            
            // Mesh resolution
            ui.heading("Mesh Resolution");
            ui.horizontal(|ui| {
                ui.selectable_value(&mut config.resolution, MeshResolution::Coarse, "Coarse");
                ui.selectable_value(&mut config.resolution, MeshResolution::Medium, "Medium");
                ui.selectable_value(&mut config.resolution, MeshResolution::Fine, "Fine");
                ui.selectable_value(&mut config.resolution, MeshResolution::VeryFine, "Very Fine");
                ui.selectable_value(&mut config.resolution, MeshResolution::Custom, "Custom");
            });
            
            // Show estimated mesh size
            let params = if config.resolution == MeshResolution::Custom {
                config.custom_params.clone()
            } else {
                crate::examples::ExampleMeshParams::default().with_resolution(config.resolution)
            };
            
            let (est_nodes, est_elements) = estimate_mesh_size(config.example_type, &params);
            ui.label(format!("Estimated: ~{} nodes, ~{} elements", est_nodes, est_elements));
            
            ui.add_space(8.0);
            
            // Custom parameters (when Custom is selected)
            if config.resolution == MeshResolution::Custom {
                ui.collapsing("Custom Mesh Parameters", |ui| {
                    match config.example_type {
                        ExampleType::CantileverBeam => {
                            egui::Grid::new("beam_params")
                                .num_columns(2)
                                .spacing([10.0, 4.0])
                                .show(ui, |ui| {
                                    ui.label("Length (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_length).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Height (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_height).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Width (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_width).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Divisions (X):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_nx).range(2..=100));
                                    ui.end_row();
                                    
                                    ui.label("Divisions (Y):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_ny).range(1..=20));
                                    ui.end_row();
                                    
                                    ui.label("Divisions (Z):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.beam_nz).range(1..=20));
                                    ui.end_row();
                                });
                        }
                        ExampleType::TorqueShaft => {
                            egui::Grid::new("shaft_params")
                                .num_columns(2)
                                .spacing([10.0, 4.0])
                                .show(ui, |ui| {
                                    ui.label("Radius (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.shaft_radius).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Length (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.shaft_length).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Radial divisions:");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.shaft_n_radial).range(4..=32));
                                    ui.end_row();
                                    
                                    ui.label("Axial divisions:");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.shaft_n_height).range(2..=50));
                                    ui.end_row();
                                    
                                    ui.label("Radial layers:");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.shaft_n_layers).range(1..=8));
                                    ui.end_row();
                                });
                        }
                        ExampleType::ContactBlocks => {
                            egui::Grid::new("contact_params")
                                .num_columns(2)
                                .spacing([10.0, 4.0])
                                .show(ui, |ui| {
                                    ui.label("Block size (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.block_size).speed(0.1));
                                    ui.end_row();
                                    
                                    ui.label("Gap (m):");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.block_gap).speed(0.001));
                                    ui.end_row();
                                    
                                    ui.label("Divisions:");
                                    ui.add(egui::DragValue::new(&mut config.custom_params.block_divisions).range(1..=20));
                                    ui.end_row();
                                });
                        }
                    }
                });
            }
            
            ui.add_space(8.0);
            ui.separator();
            ui.add_space(8.0);
            
            // Load parameters
            ui.heading("Load Parameters");
            ui.horizontal(|ui| {
                ui.label("Load magnitude:");
                ui.add(egui::DragValue::new(&mut config.load_magnitude).speed(100.0));
                ui.label(match config.example_type {
                    ExampleType::CantileverBeam => "N",
                    ExampleType::TorqueShaft => "N·m",
                    ExampleType::ContactBlocks => "N",
                });
            });
            
            ui.horizontal(|ui| {
                ui.label("Scale factor:");
                ui.add(egui::DragValue::new(&mut config.scale).speed(0.1).range(0.1..=10.0));
            });
            
            ui.add_space(16.0);
            
            // Action buttons
            ui.horizontal(|ui| {
                if ui.button("✓ Load Example").clicked() {
                    app.state.ui_state.example_dialog_open = false;
                    load_example_with_custom_config(app);
                }
                
                if ui.button("Cancel").clicked() {
                    app.state.ui_state.example_dialog_open = false;
                }
            });
        });
    
    app.state.ui_state.example_dialog_open = open;
}

/// Estimate mesh size for given example and parameters
fn estimate_mesh_size(example_type: ExampleType, params: &crate::examples::ExampleMeshParams) -> (usize, usize) {
    match example_type {
        ExampleType::CantileverBeam => {
            let nodes = (params.beam_nx + 1) * (params.beam_ny + 1) * (params.beam_nz + 1);
            let elements = params.beam_nx * params.beam_ny * params.beam_nz;
            (nodes, elements)
        }
        ExampleType::TorqueShaft => {
            let center_nodes = params.shaft_n_height + 1;
            let layer_nodes = params.shaft_n_layers * params.shaft_n_radial * (params.shaft_n_height + 1);
            let nodes = center_nodes + layer_nodes;
            let elements = params.shaft_n_height * params.shaft_n_radial * params.shaft_n_layers;
            (nodes, elements)
        }
        ExampleType::ContactBlocks => {
            let nodes_per_block = (params.block_divisions + 1).pow(3);
            let elements_per_block = params.block_divisions.pow(3);
            (nodes_per_block * 2, elements_per_block * 2)
        }
    }
}
