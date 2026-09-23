//! Top menu bar for RustFEA GUI

use eframe::egui;
use crate::app::FeaApp;
use crate::examples::{ExampleType, ExampleElementType, load_example, load_example_with_config, MeshResolution};
use crate::icons;
use crate::project_io;
use crate::state::ActivePanel;

pub fn show(ctx: &egui::Context, app: &mut FeaApp) {
    egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
        egui::menu::bar(ui, |ui| {
            // File menu
            ui.menu_button("File", |ui| {
                if ui.button(format!("{} Open Project...", icons::FOLDER_OPEN)).clicked() {
                    open_project_dialog(app);
                    ui.close_menu();
                }
                
                if ui.button(format!("{} Save Project", icons::SAVE)).clicked() {
                    save_project(app);
                    ui.close_menu();
                }
                
                if ui.button(format!("{} Save Project As...", icons::SAVE)).clicked() {
                    save_project_as(app);
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Import Mesh...", icons::FILE_ADD)).clicked() {
                    import_mesh_dialog(app);
                    ui.close_menu();
                }
                
                // Recent files submenu
                #[cfg(not(target_arch = "wasm32"))]
                if !app.state.ui_state.recent_files.files.is_empty() {
                    ui.menu_button("Recent Files", |ui| {
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
                ui.menu_button("Load Example", |ui| {
                    // Cantilever Beam
                    ui.horizontal(|ui| {
                        if ui.button("Cantilever Beam").clicked() {
                            app.state.ui_state.example_config.example_type = crate::examples::ExampleType::CantileverBeam;
                            app.state.ui_state.example_dialog_open = true;
                            ui.close_menu();
                        }
                    });
                    ui.label("  Fixed beam with tip load (bending)");
                    
                    ui.add_space(4.0);
                    
                    // Torque Shaft
                    ui.horizontal(|ui| {
                        if ui.button("Torque Shaft").clicked() {
                            app.state.ui_state.example_config.example_type = crate::examples::ExampleType::TorqueShaft;
                            app.state.ui_state.example_dialog_open = true;
                            ui.close_menu();
                        }
                    });
                    ui.label("  Cylinder with applied torque (torsion)");
                    
                    ui.add_space(4.0);
                    
                    // Contact Blocks
                    ui.horizontal(|ui| {
                        if ui.button("Contact Blocks").clicked() {
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
                
                if ui.button(format!("{} Export Results...", icons::DOWNLOAD)).clicked() {
                    export_results_dialog(app);
                    ui.close_menu();
                }
                
                #[cfg(not(target_arch = "wasm32"))]
                if ui.button(format!("{} Save Screenshot...", icons::SCREENSHOT)).clicked() {
                    app.state.ui_state.screenshot_dialog_open = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                #[cfg(not(target_arch = "wasm32"))]
                if ui.button(format!("{} Exit", icons::CLOSE)).clicked() {
                    std::process::exit(0);
                }
            });
            
            // Edit menu
            ui.menu_button("Edit", |ui| {
                let can_undo = app.state.undo_stack.can_undo();
                let can_redo = app.state.undo_stack.can_redo();
                
                let undo_text = if let Some(desc) = app.state.undo_stack.undo_description() {
                    format!("{} Undo {}", icons::ARROW_LEFT, desc)
                } else {
                    format!("{} Undo", icons::ARROW_LEFT)
                };
                
                let redo_text = if let Some(desc) = app.state.undo_stack.redo_description() {
                    format!("{} Redo {}", icons::ARROW_RIGHT, desc)
                } else {
                    format!("{} Redo", icons::ARROW_RIGHT)
                };
                
                ui.add_enabled_ui(can_undo, |ui| {
                    if ui.button(&undo_text).on_hover_text("Ctrl+Z").clicked() {
                        perform_undo(app);
                        ui.close_menu();
                    }
                });
                
                ui.add_enabled_ui(can_redo, |ui| {
                    if ui.button(&redo_text).on_hover_text("Ctrl+Shift+Z").clicked() {
                        perform_redo(app);
                        ui.close_menu();
                    }
                });
                
                ui.separator();
                
                if ui.button(format!("{} Preferences...", icons::SETTINGS)).clicked() {
                    app.state.ui_state.preferences_dialog_open = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Clear Workspace", icons::DELETE)).on_hover_text("Remove all meshes and reset to default state").clicked() {
                    // Reset state to default while preserving UI preferences
                    let display_settings = app.state.ui_state.display_settings.clone();
                    let user_settings = app.state.user_settings.clone();
                    app.state = crate::state::AppState::new();
                    app.state.ui_state.display_settings = display_settings;
                    app.state.user_settings = user_settings;
                    app.renderer = None;
                    app.render_cache.invalidate();
                    app.state.status_message = "Workspace cleared".to_string();
                    ui.close_menu();
                }
            });
            
            // View menu
            ui.menu_button("View", |ui| {
                if ui.checkbox(&mut app.state.ui_state.show_faces, format!("{} Show Faces", icons::FACES)).on_hover_text("(F)").changed() {
                    app.renderer = None;
                }
                if ui.checkbox(&mut app.state.ui_state.show_wireframe, format!("{} Show Wireframe", icons::WIREFRAME)).on_hover_text("(W)").changed() {
                    app.renderer = None;
                }
                if ui.checkbox(&mut app.state.ui_state.show_nodes, format!("{} Show Nodes", icons::POINTS)).on_hover_text("(N)").changed() {
                    app.renderer = None;
                }
                ui.checkbox(&mut app.state.ui_state.show_node_groups, format!("{} Show Node Groups", icons::NODE_TREE));
                ui.checkbox(&mut app.state.ui_state.show_boundary_conditions, format!("{} Show Boundary Conditions", icons::MARKUP)).on_hover_text("(B)");
                
                ui.separator();
                
                // Clipping plane
                if ui.checkbox(&mut app.state.ui_state.clipping_plane.enabled, format!("{} Clipping Plane", icons::SCISSORS))
                    .on_hover_text("Section view (C)")
                    .changed() 
                {
                    app.render_cache.invalidate();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Reset Camera", icons::FOCUS)).on_hover_text("(Home)").clicked() {
                    if let Some(mesh) = app.state.current_mesh() {
                        let bounds = mesh.bounds;
                        app.state.ui_state.camera.fit_to_bounds(&bounds);
                    }
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Front View", icons::VIEW_FRONT)).on_hover_text("(1)").clicked() {
                    app.state.ui_state.camera.set_front_view();
                    ui.close_menu();
                }
                
                if ui.button(format!("{} Top View", icons::VIEW_TOP)).on_hover_text("(3)").clicked() {
                    app.state.ui_state.camera.set_top_view();
                    ui.close_menu();
                }
                
                if ui.button(format!("{} Side View", icons::VIEW_LEFT)).on_hover_text("(2)").clicked() {
                    app.state.ui_state.camera.set_side_view();
                    ui.close_menu();
                }
                
                if ui.button(format!("{} Isometric", icons::VIEW_ISO)).on_hover_text("(0)").clicked() {
                    app.state.ui_state.camera.set_iso_view();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Keyboard Shortcuts...", icons::QUESTION)).on_hover_text("(?)").clicked() {
                    app.state.ui_state.show_shortcuts_help = true;
                    ui.close_menu();
                }
            });
            
            // Simulation menu
            ui.menu_button("Simulation", |ui| {
                let can_run = app.state.current_mesh().is_some() && !app.state.is_running;
                
                ui.add_enabled_ui(can_run, |ui| {
                    if ui.button(format!("{} Run Simulation", icons::PLAY)).clicked() {
                        app.start_simulation();
                        ui.close_menu();
                    }
                });
                
                ui.add_enabled_ui(app.state.is_running, |ui| {
                    if ui.button(format!("{} Stop Simulation", icons::STOP)).clicked() {
                        app.stop_simulation();
                        ui.close_menu();
                    }
                });
                
                ui.separator();
                
                if ui.button(format!("{} Solver Settings...", icons::SETTINGS)).clicked() {
                    app.state.ui_state.active_panel = ActivePanel::Setup;
                    ui.close_menu();
                }
            });
            
            // Help menu
            ui.menu_button("Help", |ui| {
                if ui.button(format!("{} Keyboard Shortcuts", icons::QUESTION)).clicked() {
                    app.state.ui_state.show_shortcuts_help = true;
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button(format!("{} Documentation", icons::EXTERNAL_LINK)).clicked() {
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
                    app.state.status_message = "Opening GitHub documentation...".to_string();
                    ui.close_menu();
                }
                
                ui.separator();
                
                if ui.button("About RustFEA").clicked() {
                    app.state.ui_state.about_dialog_open = true;
                    ui.close_menu();
                }
            });
            
            // Center attribution - use separator and centered layout
            ui.separator();
            ui.centered_and_justified(|ui| {
                if ui.link("Made by John Dews-Flick").clicked() {
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
            ui.separator();
            
            // Spacer to push Run button to the right
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                // Run Simulation button - prominent in top bar
                let can_run = app.state.current_mesh().is_some() && !app.state.is_running;
                
                if app.state.is_running {
                    // Show stop button when running
                    if ui.button("■ Stop").clicked() {
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
    
    // Clear existing meshes and add the new one
    let mesh_state = crate::state::MeshState::from_mesh(
        example.mesh,
        example.name.clone(),
        None
    );
    app.state.meshes.clear();
    app.state.meshes.push(mesh_state);
    app.state.current_mesh_idx = Some(0);
    
    // Set up boundary conditions
    app.state.simulation_config.boundary_conditions = example.boundary_conditions;
    
    // Set solver type
    app.state.simulation_config.solver = example.solver_type;
    
    // Fit camera to new mesh
    if let Some(mesh) = app.state.current_mesh() {
        let bounds = mesh.bounds;
        app.state.ui_state.camera.fit_to_bounds(&bounds);
    }
    
    // Disable wireframe for cleaner view on examples
    app.state.ui_state.show_wireframe = false;
    
    // Ensure faces and boundary conditions are visible
    app.state.ui_state.show_faces = true;
    app.state.ui_state.show_boundary_conditions = true;
    
    // Invalidate renderer and cache
    app.renderer = None;
    app.render_cache.invalidate();
    
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
    
    // Clear existing meshes and add the new one
    let mesh_state = crate::state::MeshState::from_mesh(
        example.mesh,
        example.name.clone(),
        None
    );
    app.state.meshes.clear();
    app.state.meshes.push(mesh_state);
    app.state.current_mesh_idx = Some(0);
    
    // Set up boundary conditions
    app.state.simulation_config.boundary_conditions = example.boundary_conditions;
    
    // Set solver type
    app.state.simulation_config.solver = example.solver_type;
    
    // Fit camera to new mesh
    if let Some(mesh) = app.state.current_mesh() {
        let bounds = mesh.bounds;
        app.state.ui_state.camera.fit_to_bounds(&bounds);
    }
    
    // Disable wireframe for cleaner view on examples
    app.state.ui_state.show_wireframe = false;
    
    // Ensure faces and boundary conditions are visible
    app.state.ui_state.show_faces = true;
    app.state.ui_state.show_boundary_conditions = true;
    
    // Invalidate renderer and cache
    app.renderer = None;
    app.render_cache.invalidate();
    
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
    
    egui::Window::new("Load Example")
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
                if ui.selectable_value(&mut config.example_type, ExampleType::CantileverBeam, "Cantilever Beam").clicked() 
                    && prev_type != ExampleType::CantileverBeam {
                    // Reset to default params for this example
                    config.load_magnitude = 10000.0;
                }
                if ui.selectable_value(&mut config.example_type, ExampleType::TorqueShaft, "Torque Shaft").clicked()
                    && prev_type != ExampleType::TorqueShaft {
                    config.load_magnitude = 5000.0;
                }
                if ui.selectable_value(&mut config.example_type, ExampleType::ContactBlocks, "Contact Blocks").clicked()
                    && prev_type != ExampleType::ContactBlocks {
                    config.load_magnitude = 50000.0;
                }
            });
            ui.label(config.example_type.description());
            
            ui.add_space(12.0);
            ui.separator();
            ui.add_space(8.0);
            
            // Element type selection (only for cantilever beam)
            if config.example_type == ExampleType::CantileverBeam {
                ui.heading("Element Type");
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut config.element_type, ExampleElementType::C3D8, "C3D8 (8-node)");
                    ui.selectable_value(&mut config.element_type, ExampleElementType::C3D20, "C3D20 (20-node)");
                });
                ui.label(config.element_type.description());
                
                ui.add_space(12.0);
                ui.separator();
                ui.add_space(8.0);
            }
            
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
            
            let (est_nodes, est_elements) = estimate_mesh_size(config.example_type, config.element_type, &params);
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
                if ui.button("Load Example").clicked() {
                    app.state.ui_state.example_dialog_open = false;
                    load_example_with_custom_config(app);
                }
                
                if ui.button("Cancel").clicked() {
                    app.state.ui_state.example_dialog_open = false;
                }
            });
        });
    
    // Only update from window close button if we didn't explicitly close via buttons
    if open == false {
        app.state.ui_state.example_dialog_open = false;
    }
}

/// Estimate mesh size for given example and parameters
fn estimate_mesh_size(example_type: ExampleType, element_type: ExampleElementType, params: &crate::examples::ExampleMeshParams) -> (usize, usize) {
    match example_type {
        ExampleType::CantileverBeam => {
            let elements = params.beam_nx * params.beam_ny * params.beam_nz;
            let nodes = match element_type {
                ExampleElementType::C3D8 => {
                    (params.beam_nx + 1) * (params.beam_ny + 1) * (params.beam_nz + 1)
                }
                ExampleElementType::C3D20 => {
                    // C3D20 serendipity: corner nodes + edge midpoint nodes (NO face/body centers)
                    // Corner nodes: (nx+1) * (ny+1) * (nz+1)
                    // Edge midpoints (x-direction): nx * (ny+1) * (nz+1)
                    // Edge midpoints (y-direction): (nx+1) * ny * (nz+1)
                    // Edge midpoints (z-direction): (nx+1) * (ny+1) * nz
                    let corners = (params.beam_nx + 1) * (params.beam_ny + 1) * (params.beam_nz + 1);
                    let edge_x = params.beam_nx * (params.beam_ny + 1) * (params.beam_nz + 1);
                    let edge_y = (params.beam_nx + 1) * params.beam_ny * (params.beam_nz + 1);
                    let edge_z = (params.beam_nx + 1) * (params.beam_ny + 1) * params.beam_nz;
                    corners + edge_x + edge_y + edge_z
                }
            };
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

/// Perform undo operation (public for keyboard shortcuts)
pub fn perform_undo_public(app: &mut FeaApp) {
    perform_undo(app);
}

/// Perform redo operation (public for keyboard shortcuts)
pub fn perform_redo_public(app: &mut FeaApp) {
    perform_redo(app);
}

/// Perform undo operation
fn perform_undo(app: &mut FeaApp) {
    use crate::state::UndoAction;
    
    if let Some(action) = app.state.undo_stack.pop_undo() {
        let redo_action = match &action {
            UndoAction::AddBoundaryCondition(bc) => {
                // Undo add = remove the last BC
                if let Some(idx) = app.state.simulation_config.boundary_conditions.iter().position(|b| {
                    // Match by name since we don't have exact equality
                    match (b, bc) {
                        (crate::state::BoundaryConditionConfig::Fixed(a), crate::state::BoundaryConditionConfig::Fixed(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::Load(a), crate::state::BoundaryConditionConfig::Load(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::Torque(a), crate::state::BoundaryConditionConfig::Torque(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::Contact(a), crate::state::BoundaryConditionConfig::Contact(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::Pressure(a), crate::state::BoundaryConditionConfig::Pressure(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::Traction(a), crate::state::BoundaryConditionConfig::Traction(b)) => a.name == b.name,
                        (crate::state::BoundaryConditionConfig::BodyForce(a), crate::state::BoundaryConditionConfig::BodyForce(b)) => a.name == b.name,
                        _ => false,
                    }
                }) {
                    let removed = app.state.simulation_config.boundary_conditions.remove(idx);
                    UndoAction::RemoveBoundaryCondition(idx, removed)
                } else {
                    return;
                }
            }
            UndoAction::RemoveBoundaryCondition(idx, bc) => {
                // Undo remove = add it back at the same position
                let idx = (*idx).min(app.state.simulation_config.boundary_conditions.len());
                app.state.simulation_config.boundary_conditions.insert(idx, bc.clone());
                UndoAction::AddBoundaryCondition(bc.clone())
            }
            UndoAction::ModifyBoundaryCondition(idx, old_bc) => {
                // Undo modify = restore old value
                if let Some(current) = app.state.simulation_config.boundary_conditions.get(*idx) {
                    let redo_bc = current.clone();
                    app.state.simulation_config.boundary_conditions[*idx] = old_bc.clone();
                    UndoAction::ModifyBoundaryCondition(*idx, redo_bc)
                } else {
                    return;
                }
            }
            UndoAction::AddMaterial(mat) => {
                // Undo add = remove the material
                if let Some(idx) = app.state.simulation_config.materials.iter().position(|m| m.name == mat.name) {
                    let removed = app.state.simulation_config.materials.remove(idx);
                    UndoAction::RemoveMaterial(idx, removed)
                } else {
                    return;
                }
            }
            UndoAction::RemoveMaterial(idx, mat) => {
                // Undo remove = add it back
                let idx = (*idx).min(app.state.simulation_config.materials.len());
                app.state.simulation_config.materials.insert(idx, mat.clone());
                UndoAction::AddMaterial(mat.clone())
            }
            UndoAction::ModifyMaterial(idx, old_mat) => {
                if let Some(current) = app.state.simulation_config.materials.get(*idx) {
                    let redo_mat = current.clone();
                    app.state.simulation_config.materials[*idx] = old_mat.clone();
                    UndoAction::ModifyMaterial(*idx, redo_mat)
                } else {
                    return;
                }
            }
            UndoAction::CreateNodeGroup(mesh_idx, name) => {
                // Undo create = delete the group
                if let Some(mesh) = app.state.meshes.get_mut(*mesh_idx) {
                    if let Some(nodes) = mesh.mesh.node_groups.remove(name) {
                        UndoAction::DeleteNodeGroup(*mesh_idx, name.clone(), nodes)
                    } else {
                        return;
                    }
                } else {
                    return;
                }
            }
            UndoAction::DeleteNodeGroup(mesh_idx, name, nodes) => {
                // Undo delete = recreate the group
                if let Some(mesh) = app.state.meshes.get_mut(*mesh_idx) {
                    mesh.mesh.node_groups.insert(name.clone(), nodes.clone());
                    UndoAction::CreateNodeGroup(*mesh_idx, name.clone())
                } else {
                    return;
                }
            }
            UndoAction::MeshTransform(mesh_idx, trans, scale, rot) => {
                // For now just record the current transform (would need inverse transform logic)
                UndoAction::MeshTransform(*mesh_idx, *trans, *scale, *rot)
            }
        };
        
        app.state.undo_stack.push_redo(redo_action);
        app.state.status_message = "Undone".to_string();
    }
}

/// Perform redo operation
fn perform_redo(app: &mut FeaApp) {
    use crate::state::UndoAction;
    
    if let Some(action) = app.state.undo_stack.pop_redo() {
        let undo_action = match &action {
            UndoAction::AddBoundaryCondition(bc) => {
                // Redo add = add it back
                app.state.simulation_config.boundary_conditions.push(bc.clone());
                UndoAction::AddBoundaryCondition(bc.clone())
            }
            UndoAction::RemoveBoundaryCondition(idx, _bc) => {
                // Redo remove = remove it again
                if *idx < app.state.simulation_config.boundary_conditions.len() {
                    let removed = app.state.simulation_config.boundary_conditions.remove(*idx);
                    UndoAction::RemoveBoundaryCondition(*idx, removed)
                } else {
                    return;
                }
            }
            UndoAction::ModifyBoundaryCondition(idx, new_bc) => {
                if let Some(current) = app.state.simulation_config.boundary_conditions.get(*idx) {
                    let undo_bc = current.clone();
                    app.state.simulation_config.boundary_conditions[*idx] = new_bc.clone();
                    UndoAction::ModifyBoundaryCondition(*idx, undo_bc)
                } else {
                    return;
                }
            }
            UndoAction::AddMaterial(mat) => {
                app.state.simulation_config.materials.push(mat.clone());
                UndoAction::AddMaterial(mat.clone())
            }
            UndoAction::RemoveMaterial(idx, _mat) => {
                if *idx < app.state.simulation_config.materials.len() {
                    let removed = app.state.simulation_config.materials.remove(*idx);
                    UndoAction::RemoveMaterial(*idx, removed)
                } else {
                    return;
                }
            }
            UndoAction::ModifyMaterial(idx, new_mat) => {
                if let Some(current) = app.state.simulation_config.materials.get(*idx) {
                    let undo_mat = current.clone();
                    app.state.simulation_config.materials[*idx] = new_mat.clone();
                    UndoAction::ModifyMaterial(*idx, undo_mat)
                } else {
                    return;
                }
            }
            UndoAction::CreateNodeGroup(mesh_idx, name) => {
                // Can't redo create without the nodes - would need better tracking
                UndoAction::CreateNodeGroup(*mesh_idx, name.clone())
            }
            UndoAction::DeleteNodeGroup(mesh_idx, name, nodes) => {
                if let Some(mesh) = app.state.meshes.get_mut(*mesh_idx) {
                    mesh.mesh.node_groups.remove(name);
                    UndoAction::DeleteNodeGroup(*mesh_idx, name.clone(), nodes.clone())
                } else {
                    return;
                }
            }
            UndoAction::MeshTransform(mesh_idx, trans, scale, rot) => {
                UndoAction::MeshTransform(*mesh_idx, *trans, *scale, *rot)
            }
        };
        
        app.state.undo_stack.push(undo_action);
        app.state.status_message = "Redone".to_string();
    }
}
