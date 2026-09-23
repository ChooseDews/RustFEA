//! Mesh panel for importing and managing meshes

use eframe::egui;
use crate::app::FeaApp;
use crate::state::{PrimitiveType, NodeSelectionMode};

pub fn show(ui: &mut egui::Ui, app: &mut FeaApp) {
    ui.heading("Mesh");
    ui.add_space(8.0);
    
    // Import buttons with tooltips
    ui.horizontal(|ui| {
        if ui.button("Import Mesh")
            .on_hover_text("Import mesh from file\nSupported: .inp (Abaqus), .bin, .json")
            .clicked() 
        {
            import_mesh(app);
        }
        
        if ui.button("+ New Primitive")
            .on_hover_text("Create a new primitive mesh\n(Block or Cylinder)")
            .clicked() 
        {
            app.state.ui_state.primitive_dialog_open = true;
        }
    });
    
    ui.add_space(8.0);
    ui.separator();
    ui.add_space(8.0);
    
    // Mesh list
    if app.state.meshes.is_empty() {
        ui.colored_label(egui::Color32::GRAY, "No meshes loaded");
        ui.label("Import a mesh file (.inp, .bin, .json) or create a primitive");
        
        // Quick start tips
        ui.add_space(8.0);
        ui.group(|ui| {
            ui.label("Quick Start:");
            ui.label("Use File > Load Example for sample meshes");
            ui.label("• Click New Primitive for simple shapes");
            ui.label("• Use Gmsh to create complex meshes");
        });
    } else {
        ui.heading("Loaded Meshes");
        
        let mut selected_changed = false;
        let mut new_selection = app.state.current_mesh_idx;
        let mut remove_idx = None;
        let mut rename_request = None;
        
        for (idx, mesh) in app.state.meshes.iter().enumerate() {
            let selected = app.state.current_mesh_idx == Some(idx);
            
            ui.horizontal(|ui| {
                let label = ui.selectable_label(selected, format!("● {}", mesh.name));
                if label.clicked() {
                    new_selection = Some(idx);
                    selected_changed = true;
                }
                label.on_hover_ui(|ui| {
                    ui.label(&mesh.name);
                    ui.separator();
                    ui.label(format!("Nodes: {}", mesh.mesh.nodes.len()));
                    ui.label(format!("Elements: {}", mesh.mesh.elements.len()));
                    if let Some(path) = &mesh.path {
                        ui.label(format!("Path: {}", path.display()));
                    }
                });
                
                // Only show edit controls for selected mesh
                if selected {
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        // Delete button
                        if ui.small_button("x").on_hover_text("Remove mesh").clicked() {
                            remove_idx = Some(idx);
                        }
                        // Rename button  
                        if ui.small_button("✏").on_hover_text("Rename mesh").clicked() {
                            rename_request = Some(idx);
                        }
                    });
                }
            });
        }
        
        // Handle mesh removal
        if let Some(idx) = remove_idx {
            app.state.meshes.remove(idx);
            // Update selection
            if app.state.meshes.is_empty() {
                app.state.current_mesh_idx = None;
            } else if app.state.current_mesh_idx == Some(idx) {
                // Select previous mesh or first one
                app.state.current_mesh_idx = Some(idx.saturating_sub(1).min(app.state.meshes.len() - 1));
            } else if let Some(current) = app.state.current_mesh_idx {
                if current > idx {
                    app.state.current_mesh_idx = Some(current - 1);
                }
            }
            app.renderer = None;
            app.state.status_message = "Mesh removed".to_string();
        }
        
        // Handle rename request
        if rename_request.is_some() {
            app.state.ui_state.rename_dialog_open = true;
        }
        
        if selected_changed {
            app.state.current_mesh_idx = new_selection;
            if let Some(mesh) = app.state.current_mesh() {
                let bounds = mesh.bounds;
                app.state.ui_state.camera.fit_to_bounds(&bounds);
            }
            app.renderer = None; // Invalidate renderer
        }
    }
    
    // Rename dialog
    show_rename_dialog(ui.ctx(), app);
    
    // Current mesh details
    if app.state.current_mesh().is_some() {
        ui.add_space(16.0);
        ui.separator();
        ui.add_space(8.0);
        
        ui.heading("Mesh Details");
        
        // Extract data we need before mutable borrow
        let (mesh_name, mesh_path, node_count, elem_count, body_count, bounds, element_groups) = {
            let mesh = app.state.current_mesh().unwrap();
            (
                mesh.name.clone(),
                mesh.path.clone(),
                mesh.mesh.nodes.len(),
                mesh.mesh.elements.len(),
                mesh.mesh.bodies.len(),
                mesh.bounds,
                mesh.mesh.element_groups.iter()
                    .map(|(name, group)| (name.clone(), group.elements.len(), group.el_type.clone()))
                    .collect::<Vec<_>>(),
            )
        };
        
        ui.horizontal(|ui| {
            ui.label("Name:");
            ui.strong(&mesh_name);
        });
        
        if let Some(path) = &mesh_path {
            ui.horizontal(|ui| {
                ui.label("Path:");
                ui.label(path.to_string_lossy().to_string());
            });
        }
        
        ui.add_space(8.0);
        
        egui::Grid::new("mesh_stats")
            .num_columns(2)
            .spacing([20.0, 4.0])
            .show(ui, |ui| {
                ui.label("Nodes:");
                ui.label(format!("{}", node_count));
                ui.end_row();
                
                ui.label("Elements:");
                ui.label(format!("{}", elem_count));
                ui.end_row();
                
                ui.label("Bodies:");
                ui.label(format!("{}", body_count));
                ui.end_row();
            });
        
        ui.add_space(8.0);
        
        // Bounding box
        ui.collapsing("Bounding Box", |ui| {
            egui::Grid::new("bbox_grid")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Min:");
                    ui.label(format!("({:.3}, {:.3}, {:.3})", 
                        bounds.min[0], bounds.min[1], bounds.min[2]));
                    ui.end_row();
                    
                    ui.label("Max:");
                    ui.label(format!("({:.3}, {:.3}, {:.3})", 
                        bounds.max[0], bounds.max[1], bounds.max[2]));
                    ui.end_row();
                    
                    ui.label("Diagonal:");
                    ui.label(format!("{:.3}", bounds.diagonal()));
                    ui.end_row();
                });
        });
        
        // Node groups section
        ui.add_space(8.0);
        show_node_groups_section(ui, app);
        
        // Element groups
        if !element_groups.is_empty() {
            ui.add_space(8.0);
            ui.collapsing("Element Groups", |ui| {
                for (name, count, el_type) in &element_groups {
                    ui.horizontal(|ui| {
                        ui.label(format!("🧊 {} ({} elements, {})", 
                            name, count, el_type));
                    });
                }
            });
        }
    }
    
    // View options
    ui.add_space(16.0);
    ui.separator();
    ui.add_space(8.0);
    
    ui.heading("Display Options");
    
    ui.checkbox(&mut app.state.ui_state.show_faces, "Show Faces");
    ui.checkbox(&mut app.state.ui_state.show_wireframe, "Show Wireframe");
    ui.checkbox(&mut app.state.ui_state.show_nodes, "Show Nodes");
    ui.checkbox(&mut app.state.ui_state.show_node_groups, "Highlight Node Groups");
    ui.checkbox(&mut app.state.ui_state.show_boundary_conditions, "Show Boundary Conditions");
    
    // Mesh transform controls
    if app.state.current_mesh().is_some() {
        ui.add_space(16.0);
        ui.separator();
        ui.add_space(8.0);
        show_transform_controls(ui, app);
        
        ui.add_space(16.0);
        ui.separator();
        ui.add_space(8.0);
        show_face_selection_controls(ui, app);
    }
    
    // Primitive creation dialog
    show_primitive_dialog(ui.ctx(), app);
}

/// Show mesh transform controls (scale, translate, rotate)
fn show_transform_controls(ui: &mut egui::Ui, app: &mut FeaApp) {
    use crate::state::TransformMode;
    
    let mut should_apply = false;
    let mut should_reset = false;
    
    ui.collapsing("Transform", |ui| {
        let transform = &mut app.state.ui_state.mesh_edit.transform;
        
        // Mode selection
        ui.horizontal(|ui| {
            ui.label("Mode:");
            ui.selectable_value(&mut transform.mode, TransformMode::None, "None");
            ui.selectable_value(&mut transform.mode, TransformMode::Translate, "Move");
            ui.selectable_value(&mut transform.mode, TransformMode::Scale, "Scale");
            ui.selectable_value(&mut transform.mode, TransformMode::Rotate, "Rotate");
        });
        
        ui.add_space(4.0);
        
        match transform.mode {
            TransformMode::Translate => {
                ui.label("Translation:");
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut transform.translation[0]).prefix("X: ").speed(0.01));
                    ui.add(egui::DragValue::new(&mut transform.translation[1]).prefix("Y: ").speed(0.01));
                    ui.add(egui::DragValue::new(&mut transform.translation[2]).prefix("Z: ").speed(0.01));
                });
            }
            TransformMode::Scale => {
                ui.label("Scale:");
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut transform.scale[0]).prefix("X: ").speed(0.01).range(0.01..=100.0));
                    ui.add(egui::DragValue::new(&mut transform.scale[1]).prefix("Y: ").speed(0.01).range(0.01..=100.0));
                    ui.add(egui::DragValue::new(&mut transform.scale[2]).prefix("Z: ").speed(0.01).range(0.01..=100.0));
                });
                
                // Uniform scale toggle
                ui.horizontal(|ui| {
                    if ui.small_button("Uniform").clicked() {
                        let avg = (transform.scale[0] + transform.scale[1] + transform.scale[2]) / 3.0;
                        transform.scale = [avg, avg, avg];
                    }
                    if ui.small_button("Reset").clicked() {
                        transform.scale = [1.0, 1.0, 1.0];
                    }
                });
            }
            TransformMode::Rotate => {
                ui.label("Rotation (degrees):");
                ui.horizontal(|ui| {
                    ui.add(egui::DragValue::new(&mut transform.rotation[0]).prefix("X: ").speed(1.0).suffix("°"));
                    ui.add(egui::DragValue::new(&mut transform.rotation[1]).prefix("Y: ").speed(1.0).suffix("°"));
                    ui.add(egui::DragValue::new(&mut transform.rotation[2]).prefix("Z: ").speed(1.0).suffix("°"));
                });
                
                // Quick rotation buttons
                ui.horizontal(|ui| {
                    if ui.small_button("+90° X").clicked() { transform.rotation[0] += 90.0; }
                    if ui.small_button("+90° Y").clicked() { transform.rotation[1] += 90.0; }
                    if ui.small_button("+90° Z").clicked() { transform.rotation[2] += 90.0; }
                });
            }
            TransformMode::None => {
                ui.label("Select a transform mode to edit the mesh");
            }
        }
        
        if transform.mode != TransformMode::None {
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Apply Transform").clicked() {
                    should_apply = true;
                }
                if ui.button("↺ Reset").clicked() {
                    should_reset = true;
                }
            });
        }
    });
    
    // Execute deferred actions
    if should_apply {
        apply_mesh_transform(app);
    }
    if should_reset {
        app.state.ui_state.mesh_edit.transform = crate::state::MeshTransform::default();
    }
}

/// Apply the current transform to the mesh
fn apply_mesh_transform(app: &mut FeaApp) {
    let transform = &app.state.ui_state.mesh_edit.transform;
    
    if let Some(mesh_idx) = app.state.current_mesh_idx {
        if let Some(mesh_state) = app.state.meshes.get_mut(mesh_idx) {
            let mesh = &mut mesh_state.mesh;
            
            // Apply transformations to all nodes
            for node in mesh.nodes.values_mut() {
                let mut coords = node.coordinates.clone();
                
                // Scale (around origin)
                coords[0] *= transform.scale[0];
                coords[1] *= transform.scale[1];
                coords[2] *= transform.scale[2];
                
                // Rotate (simple Euler angles around origin)
                let (rx, ry, rz) = (
                    transform.rotation[0].to_radians(),
                    transform.rotation[1].to_radians(),
                    transform.rotation[2].to_radians()
                );
                
                // Rotation around X
                if rx.abs() > 1e-10 {
                    let (y, z) = (coords[1], coords[2]);
                    coords[1] = y * rx.cos() - z * rx.sin();
                    coords[2] = y * rx.sin() + z * rx.cos();
                }
                
                // Rotation around Y
                if ry.abs() > 1e-10 {
                    let (x, z) = (coords[0], coords[2]);
                    coords[0] = x * ry.cos() + z * ry.sin();
                    coords[2] = -x * ry.sin() + z * ry.cos();
                }
                
                // Rotation around Z
                if rz.abs() > 1e-10 {
                    let (x, y) = (coords[0], coords[1]);
                    coords[0] = x * rz.cos() - y * rz.sin();
                    coords[1] = x * rz.sin() + y * rz.cos();
                }
                
                // Translate
                coords[0] += transform.translation[0];
                coords[1] += transform.translation[1];
                coords[2] += transform.translation[2];
                
                node.coordinates = coords;
            }
            
            // Recompute bounding box
            mesh_state.bounds = crate::state::compute_bounding_box_from_mesh(&mesh_state.mesh);
            
            app.state.status_message = "Transform applied".to_string();
            app.render_cache.invalidate();
            app.renderer = None;
        }
    }
    
    // Reset transform state
    app.state.ui_state.mesh_edit.transform = crate::state::MeshTransform::default();
}

/// Show face selection controls
fn show_face_selection_controls(ui: &mut egui::Ui, app: &mut FeaApp) {
    use crate::state::{FaceSelectionMode, SelectionPurpose};
    
    // Track pending actions to execute after UI drawing
    let mut pending_action: Option<&str> = None;
    let mut should_clear = false;
    
    ui.collapsing("🎯 Face Selection", |ui| {
        let selection = &mut app.state.ui_state.mesh_edit.face_selection;
        
        // Toggle face selection mode
        let was_active = selection.active;
        ui.checkbox(&mut selection.active, "Enable Face Selection");
        
        if selection.active && !was_active {
            // Will set status message after collapsing scope
        }
        
        if selection.active {
            ui.add_space(4.0);
            
            // Selection mode
            ui.horizontal(|ui| {
                ui.label("Mode:");
                ui.selectable_value(&mut selection.mode, FaceSelectionMode::Single, "Single");
                ui.selectable_value(&mut selection.mode, FaceSelectionMode::Add, "Add (+)");
                ui.selectable_value(&mut selection.mode, FaceSelectionMode::Remove, "Remove (-)");
            });
            
            // Purpose
            ui.horizontal(|ui| {
                ui.label("For:");
                egui::ComboBox::from_id_salt("selection_purpose")
                    .selected_text(match selection.purpose {
                        SelectionPurpose::General => "General",
                        SelectionPurpose::CreateNodeGroup => "Create Node Group",
                        SelectionPurpose::ApplyBC => "Apply Boundary Condition",
                        SelectionPurpose::CreateSurface => "Create Surface Group",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut selection.purpose, SelectionPurpose::General, "General");
                        ui.selectable_value(&mut selection.purpose, SelectionPurpose::CreateNodeGroup, "Create Node Group");
                        ui.selectable_value(&mut selection.purpose, SelectionPurpose::ApplyBC, "Apply Boundary Condition");
                        ui.selectable_value(&mut selection.purpose, SelectionPurpose::CreateSurface, "Create Surface Group");
                    });
            });
            
            // Selected count
            let count = selection.selected_faces.len();
            ui.label(format!("Selected: {} faces", count));
            
            // Actions based on selection - defer execution
            if count > 0 {
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    match selection.purpose {
                        SelectionPurpose::CreateNodeGroup => {
                            if ui.button("📍 Create Node Group").clicked() {
                                pending_action = Some("create_node_group");
                            }
                        }
                        SelectionPurpose::ApplyBC => {
                            if ui.button("Apply Fixed BC").clicked() {
                                pending_action = Some("apply_fixed");
                            }
                            if ui.button("Apply Load").clicked() {
                                pending_action = Some("apply_load");
                            }
                        }
                        SelectionPurpose::CreateSurface => {
                            if ui.button("Create Surface").clicked() {
                                pending_action = Some("create_surface");
                            }
                        }
                        SelectionPurpose::General => {
                            ui.label("Choose a purpose above");
                        }
                    }
                    
                    if ui.button("✗ Clear").clicked() {
                        should_clear = true;
                    }
                });
            }
            
            // Help text
            ui.add_space(4.0);
            ui.colored_label(egui::Color32::GRAY, "Shift+click to add, Ctrl+click to remove");
        }
    });
    
    // Execute deferred actions after UI drawing is complete
    if let Some(action) = pending_action {
        match action {
            "create_node_group" => create_node_group_from_faces(app),
            "apply_fixed" => create_bc_from_faces(app, "fixed"),
            "apply_load" => create_bc_from_faces(app, "load"),
            "create_surface" => create_surface_from_faces(app),
            _ => {}
        }
    }
    
    if should_clear {
        app.state.ui_state.mesh_edit.face_selection.selected_faces.clear();
    }
}

/// Create a node group from selected faces
fn create_node_group_from_faces(app: &mut FeaApp) {
    let faces = app.state.ui_state.mesh_edit.face_selection.selected_faces.clone();
    if faces.is_empty() {
        return;
    }
    
    if let Some(mesh_idx) = app.state.current_mesh_idx {
        if let Some(mesh_state) = app.state.meshes.get_mut(mesh_idx) {
            // Collect all nodes from selected faces
            let mut nodes = std::collections::HashSet::new();
            
            // Face ordering for 8-node brick element
            let face_indices = [
                [0, 1, 2, 3], // Bottom
                [4, 7, 6, 5], // Top
                [0, 4, 5, 1], // Front
                [2, 6, 7, 3], // Back
                [0, 3, 7, 4], // Left
                [1, 5, 6, 2], // Right
            ];
            
            for (elem_id, face_idx) in &faces {
                if let Some(element) = mesh_state.mesh.elements.get(elem_id) {
                    if face_idx < &6 && element.connectivity.len() >= 8 {
                        let face = &face_indices[*face_idx];
                        for &corner in face {
                            nodes.insert(element.connectivity[corner]);
                        }
                    }
                }
            }
            
            // Create node group
            let group_num = mesh_state.mesh.node_groups.len() + 1;
            let name = format!("FaceGroup_{}", group_num);
            let node_group = rust_fea::mesh::NodeGroup {
                nodes: nodes.into_iter().collect(),
                name: name.clone(),
            };
            mesh_state.mesh.node_groups.insert(name.clone(), node_group);
            
            app.state.status_message = format!("Created node group '{}' from {} faces", name, faces.len());
        }
    }
    
    // Clear selection
    app.state.ui_state.mesh_edit.face_selection.selected_faces.clear();
}

/// Create BC from selected faces
fn create_bc_from_faces(app: &mut FeaApp, bc_type: &str) {
    // First create a node group
    create_node_group_from_faces(app);
    
    // Get the name of the group we just created
    let group_name = app.state.current_mesh()
        .and_then(|m| m.mesh.node_groups.keys().last().cloned());
    
    if let Some(name) = group_name {
        match bc_type {
            "fixed" => {
                let bc = crate::state::BoundaryConditionConfig::Fixed(crate::state::FixedBcConfig {
                    name: format!("Fixed_{}", name),
                    node_group: name,
                    constrain_x: Some(0.0),
                    constrain_y: Some(0.0),
                    constrain_z: Some(0.0),
                });
                app.state.simulation_config.boundary_conditions.push(bc);
                app.state.status_message = "Created Fixed BC from face selection".to_string();
            }
            "load" => {
                let bc = crate::state::BoundaryConditionConfig::Load(crate::state::LoadBcConfig {
                    name: format!("Load_{}", name),
                    node_group: name,
                    force_x: 0.0,
                    force_y: -1000.0,
                    force_z: 0.0,
                });
                app.state.simulation_config.boundary_conditions.push(bc);
                app.state.status_message = "Created Load BC from face selection (edit in Setup)".to_string();
            }
            _ => {}
        }
        
        // Switch to Setup panel to edit the BC
        app.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
    }
}

/// Create a surface element group from selected faces
fn create_surface_from_faces(app: &mut FeaApp) {
    // For now, just create a node group - full surface element creation would need mesh operations
    create_node_group_from_faces(app);
    app.state.status_message = "Created surface group (as node group)".to_string();
}

fn show_node_groups_section(ui: &mut egui::Ui, app: &mut FeaApp) {
    let has_groups = app.state.current_mesh()
        .map(|m| !m.mesh.node_groups.is_empty())
        .unwrap_or(false);
    
    ui.collapsing("Node Groups", |ui| {
        // Show existing groups with context menu
        if has_groups {
            // Collect group info first to avoid borrow issues
            let groups_info: Vec<(String, usize)> = app.state.current_mesh()
                .map(|m| m.mesh.node_groups.iter()
                    .map(|(name, group)| (name.clone(), group.nodes.len()))
                    .collect())
                .unwrap_or_default();
            
            let mut delete_group: Option<String> = None;
            let mut action: Option<(String, &str)> = None;
            
            for (name, node_count) in &groups_info {
                let response = ui.horizontal(|ui| {
                    ui.label(format!("📍 {} ({} nodes)", name, node_count))
                }).response;
                
                // Right-click context menu
                response.context_menu(|ui| {
                    if ui.button("🎯 Highlight").clicked() {
                        action = Some((name.clone(), "highlight"));
                        ui.close_menu();
                    }
                    if ui.button("Apply Fixed BC").clicked() {
                        action = Some((name.clone(), "fixed"));
                        ui.close_menu();
                    }
                    if ui.button("Apply Load").clicked() {
                        action = Some((name.clone(), "load"));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Delete").clicked() {
                        delete_group = Some(name.clone());
                        ui.close_menu();
                    }
                });
            }
            
            // Handle actions after iteration
            if let Some((name, act)) = action {
                match act {
                    "highlight" => {
                        app.state.ui_state.show_node_groups = true;
                    }
                    "fixed" => {
                        let bc = crate::state::BoundaryConditionConfig::Fixed(crate::state::FixedBcConfig {
                            name: format!("Fixed_{}", name),
                            node_group: name.clone(),
                            constrain_x: Some(0.0),
                            constrain_y: Some(0.0),
                            constrain_z: Some(0.0),
                        });
                        app.state.simulation_config.boundary_conditions.push(bc);
                        app.state.status_message = format!("Applied Fixed BC to '{}'", name);
                    }
                    "load" => {
                        let bc = crate::state::BoundaryConditionConfig::Load(crate::state::LoadBcConfig {
                            name: format!("Load_{}", name),
                            node_group: name.clone(),
                            force_x: 0.0,
                            force_y: -1000.0,
                            force_z: 0.0,
                        });
                        app.state.simulation_config.boundary_conditions.push(bc);
                        app.state.status_message = format!("Applied Load to '{}' - edit in Setup", name);
                    }
                    _ => {}
                }
            }
            
            // Handle group deletion
            if let Some(name) = delete_group {
                if let Some(mesh_idx) = app.state.current_mesh_idx {
                    if let Some(mesh_state) = app.state.meshes.get_mut(mesh_idx) {
                        mesh_state.mesh.node_groups.remove(&name);
                        app.state.status_message = format!("Deleted node group '{}'", name);
                    }
                }
            }
        } else {
            ui.label("No node groups defined");
        }
        
        ui.add_space(8.0);
        
        // Create new group button
        if ui.button("+ Create Node Group").clicked() {
            app.state.ui_state.node_group_creator.active = true;
            app.state.ui_state.node_group_creator.name = format!("Group_{}", 
                app.state.current_mesh().map(|m| m.mesh.node_groups.len()).unwrap_or(0) + 1);
        }
        
        // Alternative: face selection
        if ui.button("🎯 Select from Faces").clicked() {
            app.state.ui_state.mesh_edit.face_selection.active = true;
            app.state.ui_state.mesh_edit.face_selection.purpose = crate::state::SelectionPurpose::CreateNodeGroup;
            app.state.status_message = "Click faces in viewport to select".to_string();
        }
        
        // Node group creator UI
        if app.state.ui_state.node_group_creator.active {
            ui.add_space(8.0);
            ui.group(|ui| {
                ui.heading("New Node Group");
                
                ui.horizontal(|ui| {
                    ui.label("Name:");
                    ui.text_edit_singleline(&mut app.state.ui_state.node_group_creator.name);
                });
                
                ui.horizontal(|ui| {
                    ui.label("Selection:");
                    egui::ComboBox::from_id_salt("selection_mode")
                        .selected_text(match app.state.ui_state.node_group_creator.selection_mode {
                            NodeSelectionMode::BoxSelect => "Box Select",
                            NodeSelectionMode::Manual => "Manual",
                            NodeSelectionMode::Plane => "Plane",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut app.state.ui_state.node_group_creator.selection_mode,
                                NodeSelectionMode::BoxSelect,
                                "Box Select"
                            );
                            ui.selectable_value(
                                &mut app.state.ui_state.node_group_creator.selection_mode,
                                NodeSelectionMode::Plane,
                                "Plane"
                            );
                        });
                });
                
                match app.state.ui_state.node_group_creator.selection_mode {
                    NodeSelectionMode::BoxSelect => {
                        show_box_selection_ui(ui, app);
                    }
                    NodeSelectionMode::Plane => {
                        show_plane_selection_ui(ui, app);
                    }
                    NodeSelectionMode::Manual => {
                        ui.label("Click nodes in viewport to select");
                    }
                }
                
                // Preview count
                let count = app.state.ui_state.node_group_creator.selected_nodes.len();
                ui.label(format!("Selected: {} nodes", count));
                
                ui.horizontal(|ui| {
                    if ui.button("Create").clicked() && count > 0 {
                        create_node_group(app);
                    }
                    if ui.button("✗ Cancel").clicked() {
                        app.state.ui_state.node_group_creator.active = false;
                        app.state.ui_state.node_group_creator.selected_nodes.clear();
                    }
                });
            });
        }
    });
}

fn show_box_selection_ui(ui: &mut egui::Ui, app: &mut FeaApp) {
    // Get mesh bounds first (immutable borrow)
    let bounds = app.state.current_mesh()
        .map(|m| m.bounds)
        .unwrap_or_default();
    
    // Now get mutable borrow for creator
    let creator = &mut app.state.ui_state.node_group_creator;
    
    ui.label("Min corner:");
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(&mut creator.box_min[0]).prefix("X: ").speed(0.1));
        ui.add(egui::DragValue::new(&mut creator.box_min[1]).prefix("Y: ").speed(0.1));
        ui.add(egui::DragValue::new(&mut creator.box_min[2]).prefix("Z: ").speed(0.1));
    });
    
    ui.label("Max corner:");
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(&mut creator.box_max[0]).prefix("X: ").speed(0.1));
        ui.add(egui::DragValue::new(&mut creator.box_max[1]).prefix("Y: ").speed(0.1));
        ui.add(egui::DragValue::new(&mut creator.box_max[2]).prefix("Z: ").speed(0.1));
    });
    
    // Quick presets - use the bounds we captured above
    ui.horizontal(|ui| {
        if ui.small_button("Min X Face").clicked() {
            creator.box_min = [bounds.min[0] as f64 - 0.001, bounds.min[1] as f64, bounds.min[2] as f64];
            creator.box_max = [bounds.min[0] as f64 + 0.001, bounds.max[1] as f64, bounds.max[2] as f64];
        }
        if ui.small_button("Max X Face").clicked() {
            creator.box_min = [bounds.max[0] as f64 - 0.001, bounds.min[1] as f64, bounds.min[2] as f64];
            creator.box_max = [bounds.max[0] as f64 + 0.001, bounds.max[1] as f64, bounds.max[2] as f64];
        }
    });
    ui.horizontal(|ui| {
        if ui.small_button("Min Y Face").clicked() {
            creator.box_min = [bounds.min[0] as f64, bounds.min[1] as f64 - 0.001, bounds.min[2] as f64];
            creator.box_max = [bounds.max[0] as f64, bounds.min[1] as f64 + 0.001, bounds.max[2] as f64];
        }
        if ui.small_button("Max Y Face").clicked() {
            creator.box_min = [bounds.min[0] as f64, bounds.max[1] as f64 - 0.001, bounds.min[2] as f64];
            creator.box_max = [bounds.max[0] as f64, bounds.max[1] as f64 + 0.001, bounds.max[2] as f64];
        }
    });
    ui.horizontal(|ui| {
        if ui.small_button("Min Z Face").clicked() {
            creator.box_min = [bounds.min[0] as f64, bounds.min[1] as f64, bounds.min[2] as f64 - 0.001];
            creator.box_max = [bounds.max[0] as f64, bounds.max[1] as f64, bounds.min[2] as f64 + 0.001];
        }
        if ui.small_button("Max Z Face").clicked() {
            creator.box_min = [bounds.min[0] as f64, bounds.min[1] as f64, bounds.max[2] as f64 - 0.001];
            creator.box_max = [bounds.max[0] as f64, bounds.max[1] as f64, bounds.max[2] as f64 + 0.001];
        }
    });
    
    // Drop mutable borrow before calling update function
    let _ = creator;
    
    if ui.button("Update Selection").clicked() {
        update_box_selection(app);
    }
}

fn show_plane_selection_ui(ui: &mut egui::Ui, app: &mut FeaApp) {
    // Get mesh bounds first (immutable borrow)
    let bounds = app.state.current_mesh()
        .map(|m| m.bounds)
        .unwrap_or_default();
    
    // Now get mutable borrow for creator
    let creator = &mut app.state.ui_state.node_group_creator;
    
    // Use box_min[0] for plane position, box_max[0] for tolerance
    ui.label("Plane selection (select nodes near a plane):");
    
    ui.horizontal(|ui| {
        ui.label("Axis:");
        if ui.small_button("X").clicked() {
            creator.box_min[1] = 0.0; // Use index 1 to store axis (0=X, 1=Y, 2=Z)
        }
        if ui.small_button("Y").clicked() {
            creator.box_min[1] = 1.0;
        }
        if ui.small_button("Z").clicked() {
            creator.box_min[1] = 2.0;
        }
    });
    
    ui.horizontal(|ui| {
        ui.label("Position:");
        ui.add(egui::DragValue::new(&mut creator.box_min[0]).speed(0.1));
    });
    
    ui.horizontal(|ui| {
        ui.label("Tolerance:");
        ui.add(egui::DragValue::new(&mut creator.box_max[0]).speed(0.01).range(0.0001..=1.0));
    });
    
    // Quick presets for plane selection - use bounds captured above
    ui.horizontal(|ui| {
        if ui.small_button("Min X").clicked() {
            creator.box_min[0] = bounds.min[0] as f64;
            creator.box_min[1] = 0.0;
            creator.box_max[0] = 0.001;
        }
        if ui.small_button("Max X").clicked() {
            creator.box_min[0] = bounds.max[0] as f64;
            creator.box_min[1] = 0.0;
            creator.box_max[0] = 0.001;
        }
        if ui.small_button("Min Y").clicked() {
            creator.box_min[0] = bounds.min[1] as f64;
            creator.box_min[1] = 1.0;
            creator.box_max[0] = 0.001;
        }
        if ui.small_button("Max Y").clicked() {
            creator.box_min[0] = bounds.max[1] as f64;
            creator.box_min[1] = 1.0;
            creator.box_max[0] = 0.001;
        }
    });
    
    // Drop mutable borrow before calling update function
    let _ = creator;
    
    if ui.button("Update Selection").clicked() {
        update_plane_selection(app);
    }
}

fn update_box_selection(app: &mut FeaApp) {
    // Extract box bounds first
    let min = app.state.ui_state.node_group_creator.box_min;
    let max = app.state.ui_state.node_group_creator.box_max;
    
    let mut selected = Vec::new();
    
    if let Some(mesh) = app.state.current_mesh() {
        for (&node_id, node) in &mesh.mesh.nodes {
            let pos = &node.coordinates;
            if pos[0] >= min[0] && pos[0] <= max[0] &&
               pos[1] >= min[1] && pos[1] <= max[1] &&
               pos[2] >= min[2] && pos[2] <= max[2] {
                selected.push(node_id);
            }
        }
    }
    
    app.state.ui_state.node_group_creator.selected_nodes = selected;
}

fn update_plane_selection(app: &mut FeaApp) {
    // Extract selection parameters first
    let position = app.state.ui_state.node_group_creator.box_min[0];
    let axis = app.state.ui_state.node_group_creator.box_min[1] as usize;
    let tolerance = app.state.ui_state.node_group_creator.box_max[0];
    
    let mut selected = Vec::new();
    
    if let Some(mesh) = app.state.current_mesh() {
        for (&node_id, node) in &mesh.mesh.nodes {
            let coord = node.coordinates[axis.min(2)];
            if (coord - position).abs() <= tolerance {
                selected.push(node_id);
            }
        }
    }
    
    app.state.ui_state.node_group_creator.selected_nodes = selected;
}

fn create_node_group(app: &mut FeaApp) {
    let name = app.state.ui_state.node_group_creator.name.clone();
    let nodes = app.state.ui_state.node_group_creator.selected_nodes.clone();
    
    if let Some(mesh_idx) = app.state.current_mesh_idx {
        if let Some(mesh_state) = app.state.meshes.get_mut(mesh_idx) {
            // Create node group
            let node_group = rust_fea::mesh::NodeGroup {
                nodes: nodes,
                name: name.clone(),
            };
            mesh_state.mesh.node_groups.insert(name.clone(), node_group);
            app.state.status_message = format!("Created node group '{}'", name);
        }
    }
    
    // Reset creator
    app.state.ui_state.node_group_creator.active = false;
    app.state.ui_state.node_group_creator.selected_nodes.clear();
}

fn show_rename_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    if !app.state.ui_state.rename_dialog_open {
        return;
    }
    
    // Get current mesh name for the dialog
    let current_name = app.state.current_mesh()
        .map(|m| m.name.clone())
        .unwrap_or_default();
    
    // Initialize rename buffer if empty
    if app.state.ui_state.rename_buffer.is_empty() {
        app.state.ui_state.rename_buffer = current_name.clone();
    }
    
    let mut close_dialog = false;
    let mut apply_rename = false;
    
    egui::Window::new("Rename Mesh")
        .collapsible(false)
        .resizable(false)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(&mut app.state.ui_state.rename_buffer);
            });
            
            ui.add_space(8.0);
            
            ui.horizontal(|ui| {
                if ui.button("Rename").clicked() {
                    apply_rename = true;
                    close_dialog = true;
                }
                if ui.button("✗ Cancel").clicked() {
                    close_dialog = true;
                }
            });
        });
    
    if apply_rename {
        let new_name = app.state.ui_state.rename_buffer.clone();
        if let Some(mesh_idx) = app.state.current_mesh_idx {
            if let Some(mesh_state) = app.state.meshes.get_mut(mesh_idx) {
                mesh_state.name = new_name.clone();
                app.state.status_message = format!("Renamed mesh to '{}'", new_name);
            }
        }
    }
    
    if close_dialog {
        app.state.ui_state.rename_dialog_open = false;
        app.state.ui_state.rename_buffer.clear();
    }
}

fn show_primitive_dialog(ctx: &egui::Context, app: &mut FeaApp) {
    if !app.state.ui_state.primitive_dialog_open {
        return;
    }
    
    egui::Window::new("Create Primitive Mesh")
        .collapsible(false)
        .resizable(true)
        .default_width(400.0)
        .show(ctx, |ui| {
            let config = &mut app.state.ui_state.primitive_config;
            
            ui.horizontal(|ui| {
                ui.label("Name:");
                ui.text_edit_singleline(&mut config.name);
            });
            
            ui.add_space(8.0);
            
            ui.horizontal(|ui| {
                ui.label("Type:");
                egui::ComboBox::from_id_salt("primitive_type")
                    .selected_text(match config.primitive_type {
                        PrimitiveType::Block => "🧊 Block",
                        PrimitiveType::Cylinder => "🔷 Cylinder",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut config.primitive_type, PrimitiveType::Block, "🧊 Block");
                        ui.selectable_value(&mut config.primitive_type, PrimitiveType::Cylinder, "🔷 Cylinder");
                    });
            });
            
            ui.add_space(8.0);
            ui.separator();
            ui.add_space(8.0);
            
            match config.primitive_type {
                PrimitiveType::Block => {
                    show_block_config_enhanced(ui, config);
                }
                PrimitiveType::Cylinder => {
                    show_cylinder_config_enhanced(ui, config);
                }
            }
            
            ui.add_space(16.0);
            
            ui.horizontal(|ui| {
                if ui.button("Create").clicked() {
                    create_primitive_mesh(app);
                    app.state.ui_state.primitive_dialog_open = false;
                }
                if ui.button("Create & Edit").on_hover_text("Create and keep dialog open for adjustments").clicked() {
                    create_primitive_mesh(app);
                    app.state.ui_state.mesh_edit.editing_primitive = true;
                }
                if ui.button("✗ Cancel").clicked() {
                    app.state.ui_state.primitive_dialog_open = false;
                    app.state.ui_state.mesh_edit.editing_primitive = false;
                }
            });
        });
}

fn show_block_config_enhanced(ui: &mut egui::Ui, config: &mut crate::state::PrimitiveConfig) {
    ui.heading("Block Parameters");
    
    // Presets
    ui.horizontal(|ui| {
        ui.label("Presets:");
        if ui.small_button("Unit Cube").clicked() {
            config.block_size = [1.0, 1.0, 1.0];
            config.block_divisions = [2, 2, 2];
            config.block_origin = [0.0, 0.0, 0.0];
        }
        if ui.small_button("Beam").clicked() {
            config.block_size = [10.0, 1.0, 1.0];
            config.block_divisions = [10, 2, 2];
            config.block_origin = [0.0, 0.0, 0.0];
        }
        if ui.small_button("Plate").clicked() {
            config.block_size = [5.0, 0.5, 5.0];
            config.block_divisions = [10, 1, 10];
            config.block_origin = [0.0, 0.0, 0.0];
        }
    });
    
    ui.add_space(8.0);
    
    // Size with sliders
    ui.group(|ui| {
        ui.label("Size:");
        ui.horizontal(|ui| {
            ui.label("Width (X):");
            ui.add(egui::DragValue::new(&mut config.block_size[0]).speed(0.1).range(0.01..=1000.0));
        });
        ui.horizontal(|ui| {
            ui.label("Height (Y):");
            ui.add(egui::DragValue::new(&mut config.block_size[1]).speed(0.1).range(0.01..=1000.0));
        });
        ui.horizontal(|ui| {
            ui.label("Depth (Z):");
            ui.add(egui::DragValue::new(&mut config.block_size[2]).speed(0.1).range(0.01..=1000.0));
        });
    });
    
    ui.add_space(4.0);
    
    // Mesh density controls
    ui.group(|ui| {
        ui.label("Mesh Density:");
        
        // Element count sliders
        let mut div_x = config.block_divisions[0] as i32;
        let mut div_y = config.block_divisions[1] as i32;
        let mut div_z = config.block_divisions[2] as i32;
        
        ui.horizontal(|ui| {
            ui.label("X divisions:");
            if ui.add(egui::Slider::new(&mut div_x, 1..=20).show_value(true)).changed() {
                config.block_divisions[0] = div_x.max(1) as usize;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Y divisions:");
            if ui.add(egui::Slider::new(&mut div_y, 1..=20).show_value(true)).changed() {
                config.block_divisions[1] = div_y.max(1) as usize;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Z divisions:");
            if ui.add(egui::Slider::new(&mut div_z, 1..=20).show_value(true)).changed() {
                config.block_divisions[2] = div_z.max(1) as usize;
            }
        });
        
        // Quick density presets
        ui.horizontal(|ui| {
            ui.label("Quick:");
            if ui.small_button("Coarse").clicked() {
                config.block_divisions = [1, 1, 1];
            }
            if ui.small_button("Medium").clicked() {
                config.block_divisions = [3, 3, 3];
            }
            if ui.small_button("Fine").clicked() {
                config.block_divisions = [6, 6, 6];
            }
            if ui.small_button("Very Fine").clicked() {
                config.block_divisions = [10, 10, 10];
            }
        });
        
        // Element size info
        let elem_size_x = config.block_size[0] / config.block_divisions[0] as f64;
        let elem_size_y = config.block_size[1] / config.block_divisions[1] as f64;
        let elem_size_z = config.block_size[2] / config.block_divisions[2] as f64;
        ui.colored_label(egui::Color32::GRAY, 
            format!("Element size: {:.3} × {:.3} × {:.3}", elem_size_x, elem_size_y, elem_size_z));
    });
    
    ui.add_space(4.0);
    
    // Origin
    ui.collapsing("Origin / Position", |ui| {
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut config.block_origin[0]).prefix("X: ").speed(0.1));
            ui.add(egui::DragValue::new(&mut config.block_origin[1]).prefix("Y: ").speed(0.1));
            ui.add(egui::DragValue::new(&mut config.block_origin[2]).prefix("Z: ").speed(0.1));
        });
        if ui.small_button("Center at Origin").clicked() {
            config.block_origin = [
                -config.block_size[0] / 2.0,
                -config.block_size[1] / 2.0,
                -config.block_size[2] / 2.0,
            ];
        }
    });
    
    // Preview info
    ui.add_space(8.0);
    let num_elements = config.block_divisions[0] * config.block_divisions[1] * config.block_divisions[2];
    let num_nodes = (config.block_divisions[0] + 1) * (config.block_divisions[1] + 1) * (config.block_divisions[2] + 1);
    ui.separator();
    ui.colored_label(egui::Color32::LIGHT_GREEN, 
        format!("{} elements, {} nodes", num_elements, num_nodes));
}

fn show_cylinder_config_enhanced(ui: &mut egui::Ui, config: &mut crate::state::PrimitiveConfig) {
    ui.heading("Cylinder Parameters");
    
    // Presets
    ui.horizontal(|ui| {
        ui.label("Presets:");
        if ui.small_button("Shaft").clicked() {
            config.cyl_radius = 0.1;
            config.cyl_height = 1.0;
            config.cyl_radial_divisions = 8;
            config.cyl_height_divisions = 8;
        }
        if ui.small_button("Pipe").clicked() {
            config.cyl_radius = 0.5;
            config.cyl_height = 2.0;
            config.cyl_radial_divisions = 12;
            config.cyl_height_divisions = 4;
        }
        if ui.small_button("Disk").clicked() {
            config.cyl_radius = 1.0;
            config.cyl_height = 0.1;
            config.cyl_radial_divisions = 16;
            config.cyl_height_divisions = 1;
        }
    });
    
    ui.add_space(8.0);
    
    // Dimensions
    ui.group(|ui| {
        ui.label("Dimensions:");
        ui.horizontal(|ui| {
            ui.label("Radius:");
            ui.add(egui::DragValue::new(&mut config.cyl_radius).speed(0.01).range(0.01..=1000.0));
        });
        ui.horizontal(|ui| {
            ui.label("Height:");
            ui.add(egui::DragValue::new(&mut config.cyl_height).speed(0.1).range(0.01..=1000.0));
        });
    });
    
    ui.add_space(4.0);
    
    // Mesh density
    ui.group(|ui| {
        ui.label("Mesh Density:");
        
        let mut radial = config.cyl_radial_divisions as i32;
        let mut height = config.cyl_height_divisions as i32;
        
        ui.horizontal(|ui| {
            ui.label("Circumferential:");
            if ui.add(egui::Slider::new(&mut radial, 4..=32).show_value(true)).changed() {
                config.cyl_radial_divisions = radial.max(4) as usize;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Axial:");
            if ui.add(egui::Slider::new(&mut height, 1..=20).show_value(true)).changed() {
                config.cyl_height_divisions = height.max(1) as usize;
            }
        });
        
        // Quick density presets
        ui.horizontal(|ui| {
            ui.label("Quick:");
            if ui.small_button("Coarse").clicked() {
                config.cyl_radial_divisions = 6;
                config.cyl_height_divisions = 2;
            }
            if ui.small_button("Medium").clicked() {
                config.cyl_radial_divisions = 12;
                config.cyl_height_divisions = 6;
            }
            if ui.small_button("Fine").clicked() {
                config.cyl_radial_divisions = 24;
                config.cyl_height_divisions = 12;
            }
        });
    });
    
    ui.add_space(4.0);
    
    // Axis and origin
    ui.collapsing("Axis & Position", |ui| {
        ui.label("Axis direction:");
        ui.horizontal(|ui| {
            if ui.small_button("X").clicked() { config.cyl_axis = [1.0, 0.0, 0.0]; }
            if ui.small_button("Y").clicked() { config.cyl_axis = [0.0, 1.0, 0.0]; }
            if ui.small_button("Z").clicked() { config.cyl_axis = [0.0, 0.0, 1.0]; }
        });
        
        ui.label("Base center:");
        ui.horizontal(|ui| {
            ui.add(egui::DragValue::new(&mut config.cyl_origin[0]).prefix("X: ").speed(0.1));
            ui.add(egui::DragValue::new(&mut config.cyl_origin[1]).prefix("Y: ").speed(0.1));
            ui.add(egui::DragValue::new(&mut config.cyl_origin[2]).prefix("Z: ").speed(0.1));
        });
    });
    
    // Preview
    ui.add_space(8.0);
    let num_elements = config.cyl_radial_divisions * config.cyl_height_divisions;
    let num_nodes = config.cyl_radial_divisions * (config.cyl_height_divisions + 1) + 2; // +2 for center nodes
    ui.separator();
    ui.colored_label(egui::Color32::LIGHT_GREEN, 
        format!("~{} elements, ~{} nodes", num_elements, num_nodes));
}

// Keep original functions for compatibility
fn show_block_config(ui: &mut egui::Ui, config: &mut crate::state::PrimitiveConfig) {
    show_block_config_enhanced(ui, config);
}

fn show_cylinder_config(ui: &mut egui::Ui, config: &mut crate::state::PrimitiveConfig) {
    show_cylinder_config_enhanced(ui, config);
}

fn create_primitive_mesh(app: &mut FeaApp) {
    // Clone the config data we need before mutating app
    let config = app.state.ui_state.primitive_config.clone();
    
    let mesh = match config.primitive_type {
        PrimitiveType::Block => create_block_mesh(&config),
        PrimitiveType::Cylinder => create_cylinder_mesh(&config),
    };
    
    let name = config.name.clone();
    app.state.add_mesh(mesh, name.clone(), None);
    app.state.status_message = format!("Created primitive mesh '{}'", name);
    
    // Fit camera to new mesh
    if let Some(mesh) = app.state.current_mesh() {
        let bounds = mesh.bounds;
        app.state.ui_state.camera.fit_to_bounds(&bounds);
    }
}

fn create_block_mesh(config: &crate::state::PrimitiveConfig) -> rust_fea::mesh::MeshAssembly {
    use rust_fea::mesh::{MeshAssembly, MeshElement, MeshNode, NodeGroup, ElementGroup};
    use std::collections::HashMap;
    
    let mut mesh = MeshAssembly::empty();
    mesh.name = config.name.clone();
    
    let [nx, ny, nz] = config.block_divisions;
    let [sx, sy, sz] = config.block_size;
    let [ox, oy, oz] = config.block_origin;
    
    // Create nodes
    let mut node_id = 0;
    let mut node_ids: Vec<Vec<Vec<usize>>> = Vec::new();
    
    for iz in 0..=nz {
        let mut plane = Vec::new();
        for iy in 0..=ny {
            let mut row = Vec::new();
            for ix in 0..=nx {
                let x = ox + (ix as f64 / nx as f64) * sx;
                let y = oy + (iy as f64 / ny as f64) * sy;
                let z = oz + (iz as f64 / nz as f64) * sz;
                
                mesh.nodes.insert(node_id, MeshNode { 
                    coordinates: vec![x, y, z],
                    id: node_id,
                });
                row.push(node_id);
                node_id += 1;
            }
            plane.push(row);
        }
        node_ids.push(plane);
    }
    
    // Create 8-node brick elements
    let mut elem_id = 0;
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                // Brick element node ordering (following standard convention)
                let connectivity = vec![
                    node_ids[iz][iy][ix],         // 0
                    node_ids[iz][iy][ix + 1],     // 1
                    node_ids[iz][iy + 1][ix + 1], // 2
                    node_ids[iz][iy + 1][ix],     // 3
                    node_ids[iz + 1][iy][ix],     // 4
                    node_ids[iz + 1][iy][ix + 1], // 5
                    node_ids[iz + 1][iy + 1][ix + 1], // 6
                    node_ids[iz + 1][iy + 1][ix],     // 7
                ];
                
                mesh.elements.insert(elem_id, MeshElement {
                    el_type: "C3D8".to_string(),
                    connectivity,
                    name: format!("Element_{}", elem_id),
                    id: elem_id,
                });
                elem_id += 1;
            }
        }
    }
    
    // Create boundary node groups for each face
    let mut node_groups = HashMap::new();
    
    // MinX face (x = ox)
    let mut min_x_nodes = Vec::new();
    for iz in 0..=nz {
        for iy in 0..=ny {
            min_x_nodes.push(node_ids[iz][iy][0]);
        }
    }
    node_groups.insert("MinX".to_string(), NodeGroup { 
        nodes: min_x_nodes,
        name: "MinX".to_string(),
    });
    
    // MaxX face
    let mut max_x_nodes = Vec::new();
    for iz in 0..=nz {
        for iy in 0..=ny {
            max_x_nodes.push(node_ids[iz][iy][nx]);
        }
    }
    node_groups.insert("MaxX".to_string(), NodeGroup { 
        nodes: max_x_nodes,
        name: "MaxX".to_string(),
    });
    
    // MinY face
    let mut min_y_nodes = Vec::new();
    for iz in 0..=nz {
        for ix in 0..=nx {
            min_y_nodes.push(node_ids[iz][0][ix]);
        }
    }
    node_groups.insert("MinY".to_string(), NodeGroup { 
        nodes: min_y_nodes,
        name: "MinY".to_string(),
    });
    
    // MaxY face
    let mut max_y_nodes = Vec::new();
    for iz in 0..=nz {
        for ix in 0..=nx {
            max_y_nodes.push(node_ids[iz][ny][ix]);
        }
    }
    node_groups.insert("MaxY".to_string(), NodeGroup { 
        nodes: max_y_nodes,
        name: "MaxY".to_string(),
    });
    
    // MinZ face
    let mut min_z_nodes = Vec::new();
    for iy in 0..=ny {
        for ix in 0..=nx {
            min_z_nodes.push(node_ids[0][iy][ix]);
        }
    }
    node_groups.insert("MinZ".to_string(), NodeGroup { 
        nodes: min_z_nodes,
        name: "MinZ".to_string(),
    });
    
    // MaxZ face
    let mut max_z_nodes = Vec::new();
    for iy in 0..=ny {
        for ix in 0..=nx {
            max_z_nodes.push(node_ids[nz][iy][ix]);
        }
    }
    node_groups.insert("MaxZ".to_string(), NodeGroup { 
        nodes: max_z_nodes,
        name: "MaxZ".to_string(),
    });
    
    mesh.node_groups = node_groups;
    
    // Create element group
    let all_elements: Vec<usize> = (0..elem_id).collect();
    mesh.element_groups.insert("All".to_string(), ElementGroup {
        el_type: "C3D8".to_string(),
        elements: all_elements,
        name: "All".to_string(),
    });
    
    // Create single body
    mesh.single_body();
    
    mesh
}

fn create_cylinder_mesh(config: &crate::state::PrimitiveConfig) -> rust_fea::mesh::MeshAssembly {
    use rust_fea::mesh::{MeshAssembly, MeshElement, MeshNode, NodeGroup, ElementGroup};
    use std::collections::HashMap;
    
    let mut mesh = MeshAssembly::empty();
    mesh.name = config.name.clone();
    
    let n_radial = config.cyl_radial_divisions;
    let n_height = config.cyl_height_divisions;
    let radius = config.cyl_radius;
    let height = config.cyl_height;
    let [ox, oy, oz] = config.cyl_origin;
    
    // Create nodes in cylindrical layers
    // For simplicity, create a hollow cylinder with inner radius = 0.5 * outer radius
    let inner_radius = radius * 0.5;
    
    let mut node_id = 0;
    // node_ids[height_layer][radial_position][inner/outer]
    let mut node_ids: Vec<Vec<[usize; 2]>> = Vec::new();
    
    for ih in 0..=n_height {
        let y = oy + (ih as f64 / n_height as f64) * height;
        let mut layer = Vec::new();
        
        for ir in 0..n_radial {
            let angle = 2.0 * std::f64::consts::PI * (ir as f64 / n_radial as f64);
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            
            // Inner node
            let inner_x = ox + inner_radius * cos_a;
            let inner_z = oz + inner_radius * sin_a;
            mesh.nodes.insert(node_id, MeshNode { 
                coordinates: vec![inner_x, y, inner_z],
                id: node_id,
            });
            let inner_id = node_id;
            node_id += 1;
            
            // Outer node
            let outer_x = ox + radius * cos_a;
            let outer_z = oz + radius * sin_a;
            mesh.nodes.insert(node_id, MeshNode { 
                coordinates: vec![outer_x, y, outer_z],
                id: node_id,
            });
            let outer_id = node_id;
            node_id += 1;
            
            layer.push([inner_id, outer_id]);
        }
        node_ids.push(layer);
    }
    
    // Create brick elements
    let mut elem_id = 0;
    for ih in 0..n_height {
        for ir in 0..n_radial {
            let ir_next = (ir + 1) % n_radial;
            
            // Bottom layer
            let b_inner = node_ids[ih][ir][0];
            let b_outer = node_ids[ih][ir][1];
            let b_inner_next = node_ids[ih][ir_next][0];
            let b_outer_next = node_ids[ih][ir_next][1];
            
            // Top layer
            let t_inner = node_ids[ih + 1][ir][0];
            let t_outer = node_ids[ih + 1][ir][1];
            let t_inner_next = node_ids[ih + 1][ir_next][0];
            let t_outer_next = node_ids[ih + 1][ir_next][1];
            
            // Create brick element
            let connectivity = vec![
                b_inner, b_outer, b_outer_next, b_inner_next,
                t_inner, t_outer, t_outer_next, t_inner_next,
            ];
            
            mesh.elements.insert(elem_id, MeshElement {
                el_type: "C3D8".to_string(),
                connectivity,
                name: format!("Element_{}", elem_id),
                id: elem_id,
            });
            elem_id += 1;
        }
    }
    
    // Create node groups
    let mut node_groups = HashMap::new();
    
    // Bottom face
    let bottom_nodes: Vec<usize> = node_ids[0].iter()
        .flat_map(|[inner, outer]| vec![*inner, *outer])
        .collect();
    node_groups.insert("Bottom".to_string(), NodeGroup { 
        nodes: bottom_nodes,
        name: "Bottom".to_string(),
    });
    
    // Top face
    let top_nodes: Vec<usize> = node_ids[n_height].iter()
        .flat_map(|[inner, outer]| vec![*inner, *outer])
        .collect();
    node_groups.insert("Top".to_string(), NodeGroup { 
        nodes: top_nodes,
        name: "Top".to_string(),
    });
    
    // Inner surface
    let inner_nodes: Vec<usize> = node_ids.iter()
        .flat_map(|layer| layer.iter().map(|[inner, _]| *inner))
        .collect();
    node_groups.insert("Inner".to_string(), NodeGroup { 
        nodes: inner_nodes,
        name: "Inner".to_string(),
    });
    
    // Outer surface
    let outer_nodes: Vec<usize> = node_ids.iter()
        .flat_map(|layer| layer.iter().map(|[_, outer]| *outer))
        .collect();
    node_groups.insert("Outer".to_string(), NodeGroup { 
        nodes: outer_nodes,
        name: "Outer".to_string(),
    });
    
    mesh.node_groups = node_groups;
    
    // Element group
    let all_elements: Vec<usize> = (0..elem_id).collect();
    mesh.element_groups.insert("All".to_string(), ElementGroup {
        el_type: "C3D8".to_string(),
        elements: all_elements,
        name: "All".to_string(),
    });
    
    // Create single body
    mesh.single_body();
    
    mesh
}

fn import_mesh(app: &mut FeaApp) {
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
