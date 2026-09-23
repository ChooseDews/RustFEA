//! 3D Viewport for mesh visualization using egui's painter
//!
//! This provides a software-rendered 3D view that works on both native and web.
//! For high-performance rendering with wgpu, see the renderer module.

use crate::app::FeaApp;
use crate::examples::{load_example_with_config, ExampleConfig, ExampleType, MeshResolution};
use crate::icons;
use crate::render_cache;
use crate::section_cut;
use crate::state::{BodyForceTypeConfig, BoundaryConditionConfig, ClipAxis, ColorMode};
use eframe::egui;

pub fn show(ctx: &egui::Context, app: &mut FeaApp) {
    // Handle keyboard shortcuts (global)
    handle_keyboard_shortcuts(ctx, app);

    egui::CentralPanel::default().show(ctx, |ui| {
        // Viewport toolbar at top
        ui.horizontal(|ui| {
            ui.label("View:");

            // View toggles with keyboard shortcut hints
            if ui
                .selectable_label(
                    app.state.ui_state.show_faces,
                    format!("{} Faces", icons::FACES),
                )
                .on_hover_text("Toggle face rendering (F)")
                .clicked()
            {
                app.state.ui_state.show_faces = !app.state.ui_state.show_faces;
                app.renderer = None;
            }
            if ui
                .selectable_label(
                    app.state.ui_state.show_wireframe,
                    format!("{} Wire", icons::WIREFRAME),
                )
                .on_hover_text("Toggle wireframe (W)")
                .clicked()
            {
                app.state.ui_state.show_wireframe = !app.state.ui_state.show_wireframe;
                app.renderer = None;
            }
            if ui
                .selectable_label(
                    app.state.ui_state.show_nodes,
                    format!("{} Nodes", icons::POINTS),
                )
                .on_hover_text("Toggle node display (N)")
                .clicked()
            {
                app.state.ui_state.show_nodes = !app.state.ui_state.show_nodes;
                app.renderer = None;
            }
            if ui
                .selectable_label(
                    app.state.ui_state.show_boundary_conditions,
                    format!("{} BCs", icons::MARKUP),
                )
                .on_hover_text("Toggle boundary conditions (B)")
                .clicked()
            {
                app.state.ui_state.show_boundary_conditions =
                    !app.state.ui_state.show_boundary_conditions;
            }

            // Clipping plane toggle
            if ui
                .selectable_label(
                    app.state.ui_state.clipping_plane.enabled,
                    format!("{} Clip", icons::SCISSORS),
                )
                .on_hover_text("Toggle section view / clipping plane (C)")
                .clicked()
            {
                app.state.ui_state.clipping_plane.enabled =
                    !app.state.ui_state.clipping_plane.enabled;
                app.render_cache.invalidate();
            }

            // Grid toggle
            if ui
                .selectable_label(
                    app.state.ui_state.display_settings.show_grid,
                    format!("{} Grid", icons::GRID),
                )
                .on_hover_text("Toggle grid visibility (G)")
                .clicked()
            {
                app.state.ui_state.display_settings.show_grid =
                    !app.state.ui_state.display_settings.show_grid;
            }

            // Stats overlay toggle
            if ui
                .selectable_label(
                    app.state.ui_state.stats_overlay.visible,
                    format!("{} Stats", icons::INFO),
                )
                .on_hover_text("Toggle statistics overlay (I)")
                .clicked()
            {
                app.state.ui_state.stats_overlay.visible =
                    !app.state.ui_state.stats_overlay.visible;
            }

            ui.separator();

            // Quick camera views with keyboard hints
            if ui.button("Front").on_hover_text("Front view (1)").clicked() {
                app.state.ui_state.camera.yaw = 0.0;
                app.state.ui_state.camera.pitch = 0.0;
            }
            if ui.button("Side").on_hover_text("Side view (2)").clicked() {
                app.state.ui_state.camera.yaw = std::f32::consts::FRAC_PI_2;
                app.state.ui_state.camera.pitch = 0.0;
            }
            if ui.button("Top").on_hover_text("Top view (3)").clicked() {
                app.state.ui_state.camera.yaw = 0.0;
                app.state.ui_state.camera.pitch = std::f32::consts::FRAC_PI_2 - 0.01;
            }
            if ui
                .button(format!("{} Iso", icons::VIEW_ISO))
                .on_hover_text("Isometric view (0)")
                .clicked()
            {
                app.state.ui_state.camera.yaw = 0.785; // 45°
                app.state.ui_state.camera.pitch = 0.524; // 30°
            }

            ui.separator();

            // Projection mode toggle
            let proj_label = if app.state.ui_state.camera.orthographic {
                "Ortho"
            } else {
                "Persp"
            };
            let proj_hover = if app.state.ui_state.camera.orthographic {
                "Switch to perspective projection (P)"
            } else {
                "Switch to orthographic projection (P)"
            };
            if ui.button(proj_label).on_hover_text(proj_hover).clicked() {
                app.state.ui_state.camera.orthographic = !app.state.ui_state.camera.orthographic;
            }

            ui.separator();

            // Fit view
            if ui
                .button(format!("{} Fit", icons::FOCUS))
                .on_hover_text("Fit view to mesh (Home)")
                .clicked()
            {
                if let Some(mesh) = app.state.current_mesh() {
                    let bounds = mesh.bounds;
                    app.state.ui_state.camera.fit_to_bounds(&bounds);
                }
            }

            // Screenshot button
            #[cfg(not(target_arch = "wasm32"))]
            if ui
                .button(icons::CAMERA)
                .on_hover_text("Save screenshot (Ctrl+S)")
                .clicked()
            {
                app.state.ui_state.screenshot_dialog_open = true;
            }

            // Keyboard shortcuts help
            if ui
                .button(icons::QUESTION)
                .on_hover_text("Show keyboard shortcuts (?)")
                .clicked()
            {
                app.state.ui_state.show_shortcuts_help = !app.state.ui_state.show_shortcuts_help;
            }
        });

        // Show clipping plane controls if enabled
        if app.state.ui_state.clipping_plane.enabled {
            show_clipping_controls(ui, app);
        }

        ui.separator();

        let available_size = ui.available_size();

        // Allocate the viewport area
        let (rect, response) =
            ui.allocate_exact_size(available_size, egui::Sense::click_and_drag());

        // Handle camera controls
        handle_camera_input(&response, app);

        // Draw background using display settings
        let painter = ui.painter_at(rect);
        let bg = &app.state.ui_state.display_settings.background_color;
        painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(bg[0], bg[1], bg[2]));

        // Draw grid
        draw_grid(&painter, rect, app);

        // Draw mesh if available
        if app.state.current_mesh().is_some() {
            // Extract needed data to avoid borrow conflicts
            let mesh_idx = app.state.current_mesh_idx.unwrap();
            draw_mesh_cached(&painter, rect, app, mesh_idx);

            // Draw selected faces highlight
            if app.state.ui_state.mesh_edit.face_selection.active {
                draw_selected_faces(&painter, rect, app);
            }

            // For boundary conditions and node selection, we can borrow immutably
            let mesh_state = app.state.current_mesh().unwrap();

            // Draw boundary conditions
            if app.state.ui_state.show_boundary_conditions {
                draw_boundary_conditions(&painter, rect, app, mesh_state);
            }

            // Draw node group selection preview
            if app.state.ui_state.node_group_creator.active {
                draw_selected_nodes(&painter, rect, app, mesh_state);
            }
        } else {
            // Show quick-load example buttons in empty state
            let center = rect.center();

            // Calculate total height of the UI block to center it properly
            let total_height = 24.0 + 20.0 + 16.0 + 15.0 + 32.0 + 15.0 + 20.0; // ~142px
            let top_y = center.y - total_height / 2.0;

            // Draw the "Import a mesh to begin" text
            painter.text(
                egui::pos2(center.x, top_y),
                egui::Align2::CENTER_TOP,
                "Import a mesh to begin",
                egui::FontId::proportional(24.0),
                egui::Color32::from_gray(100),
            );

            // "Open an example:" label
            painter.text(
                egui::pos2(center.x, top_y + 44.0),
                egui::Align2::CENTER_TOP,
                "Open an example:",
                egui::FontId::proportional(16.0),
                egui::Color32::from_gray(90),
            );

            // Three example buttons in a row (default to Medium resolution)
            let button_width = 140.0;
            let button_height = 32.0;
            let h_spacing = 12.0;
            let total_buttons_width = 3.0 * button_width + 2.0 * h_spacing;
            let start_x = center.x - total_buttons_width / 2.0;
            let buttons_y = top_y + 75.0;

            // Track which example was clicked
            let mut clicked_example: Option<ExampleType> = None;

            // Cantilever Beam button
            let btn_rect = egui::Rect::from_min_size(
                egui::pos2(start_x, buttons_y),
                egui::vec2(button_width, button_height),
            );
            if ui
                .put(btn_rect, egui::Button::new("Cantilever Beam"))
                .clicked()
            {
                clicked_example = Some(ExampleType::CantileverBeam);
            }

            // Torque Shaft button
            let btn_rect = egui::Rect::from_min_size(
                egui::pos2(start_x + button_width + h_spacing, buttons_y),
                egui::vec2(button_width, button_height),
            );
            if ui
                .put(btn_rect, egui::Button::new("Torque Shaft"))
                .clicked()
            {
                clicked_example = Some(ExampleType::TorqueShaft);
            }

            // Contact Blocks button
            let btn_rect = egui::Rect::from_min_size(
                egui::pos2(start_x + 2.0 * (button_width + h_spacing), buttons_y),
                egui::vec2(button_width, button_height),
            );
            if ui
                .put(btn_rect, egui::Button::new("Contact Blocks"))
                .clicked()
            {
                clicked_example = Some(ExampleType::ContactBlocks);
            }

            // Handle example loading (default to Medium resolution)
            if let Some(example_type) = clicked_example {
                let config = ExampleConfig {
                    example_type,
                    resolution: MeshResolution::Medium,
                    ..Default::default()
                };
                let example = load_example_with_config(&config);

                // Load the mesh using the MeshState constructor
                let mesh_state = crate::state::MeshState::from_mesh(
                    example.mesh,
                    example_type.name().to_string(),
                    None,
                );
                app.state.meshes.clear();
                app.state.meshes.push(mesh_state);
                app.state.current_mesh_idx = Some(0);

                // Load the boundary conditions into simulation config
                app.state.simulation_config.boundary_conditions = example.boundary_conditions;

                // Set the solver type in simulation config
                app.state.simulation_config.solver = example.solver_type;

                // Fit view to mesh
                if let Some(mesh_state) = app.state.current_mesh() {
                    let bounds = mesh_state.bounds.clone();
                    app.state.ui_state.camera.fit_to_bounds(&bounds);
                }

                // Clear renderer and cache
                app.renderer = None;
                app.render_cache.invalidate();

                app.state.status_message = format!("Loaded {} (Medium mesh)", example_type.name());
            }

            // "More Examples" link to open the examples panel
            let more_y = buttons_y + button_height + 18.0;
            let more_rect =
                egui::Rect::from_center_size(egui::pos2(center.x, more_y), egui::vec2(120.0, 24.0));
            if ui
                .put(more_rect, egui::Button::new("More Examples").small())
                .clicked()
            {
                app.state.ui_state.example_dialog_open = true;
            }
        }

        // Draw axis indicator
        if app.state.ui_state.display_settings.show_axis {
            draw_axis_indicator(&painter, rect, app);
        }

        // Draw view info
        draw_view_info(&painter, rect, app);

        // Handle face selection clicks
        if app.state.ui_state.mesh_edit.face_selection.active {
            handle_face_selection(&response, rect, app);
        }

        // Context menu (right-click)
        show_viewport_context_menu(&response, app);
    });
}

/// Handle face selection in viewport
fn handle_face_selection(response: &egui::Response, rect: egui::Rect, app: &mut FeaApp) {
    // Check for click (not drag)
    if response.clicked() {
        if let Some(pos) = response.interact_pointer_pos() {
            // Try to find a face at this position
            if let Some((elem_id, face_idx)) = find_face_at_position(pos, rect, app) {
                let selection = &mut app.state.ui_state.mesh_edit.face_selection;
                let face = (elem_id, face_idx);

                // Check modifiers for selection mode
                let add_mode = response.ctx.input(|i| i.modifiers.shift);
                let remove_mode = response
                    .ctx
                    .input(|i| i.modifiers.ctrl || i.modifiers.command);

                if remove_mode {
                    // Remove from selection
                    selection.selected_faces.retain(|f| f != &face);
                } else if add_mode {
                    // Add to selection if not already present
                    if !selection.selected_faces.contains(&face) {
                        selection.selected_faces.push(face);
                    }
                } else {
                    // Replace selection
                    selection.selected_faces.clear();
                    selection.selected_faces.push(face);
                }

                app.state.status_message =
                    format!("Selected {} faces", selection.selected_faces.len());
            }
        }
    }
}

/// Find which face (if any) is at the given screen position
fn find_face_at_position(
    pos: egui::Pos2,
    rect: egui::Rect,
    app: &FeaApp,
) -> Option<(usize, usize)> {
    let cache = &app.render_cache;
    let camera = &app.state.ui_state.camera;

    if cache.faces.is_empty() {
        return None;
    }

    let aspect = rect.width() / rect.height();
    let rect_center = rect.center();
    let half_width = rect.width() / 2.0;
    let half_height = rect.height() / 2.0;

    let view_transform = match crate::render_cache::ViewTransform::from_camera(camera) {
        Some(vt) => vt,
        None => return None,
    };

    // Find the closest face that contains this point
    let mut best_face: Option<(usize, usize, f32)> = None; // (elem_id, face_idx, depth)

    for face in &cache.faces {
        // Project face corners
        let proj: Vec<Option<egui::Pos2>> = face
            .corners
            .iter()
            .map(|c| view_transform.project(*c, rect_center, half_width, half_height, aspect))
            .collect();

        // Skip if any corner is behind camera
        if proj.iter().any(|p| p.is_none()) {
            continue;
        }

        let proj: Vec<egui::Pos2> = proj.into_iter().map(|p| p.unwrap()).collect();

        // Check if point is inside this quad (split into two triangles)
        let inside = point_in_triangle(pos, proj[0], proj[1], proj[2])
            || point_in_triangle(pos, proj[0], proj[2], proj[3]);

        if inside {
            // Calculate depth
            let depth = face.center[0] * camera.yaw.sin() + face.center[2] * camera.yaw.cos();

            if best_face.is_none() || depth < best_face.unwrap().2 {
                best_face = Some((face.element_id, face.face_index, depth));
            }
        }
    }

    best_face.map(|(e, f, _)| (e, f))
}

/// Check if a point is inside a triangle
fn point_in_triangle(p: egui::Pos2, a: egui::Pos2, b: egui::Pos2, c: egui::Pos2) -> bool {
    let v0 = egui::vec2(c.x - a.x, c.y - a.y);
    let v1 = egui::vec2(b.x - a.x, b.y - a.y);
    let v2 = egui::vec2(p.x - a.x, p.y - a.y);

    let dot00 = v0.x * v0.x + v0.y * v0.y;
    let dot01 = v0.x * v1.x + v0.y * v1.y;
    let dot02 = v0.x * v2.x + v0.y * v2.y;
    let dot11 = v1.x * v1.x + v1.y * v1.y;
    let dot12 = v1.x * v2.x + v1.y * v2.y;

    let inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    let u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
    let v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

    u >= 0.0 && v >= 0.0 && (u + v) <= 1.0
}

/// Draw selected faces with highlight
fn draw_selected_faces(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    let selection = &app.state.ui_state.mesh_edit.face_selection;
    if selection.selected_faces.is_empty() {
        return;
    }

    let cache = &app.render_cache;
    let camera = &app.state.ui_state.camera;

    let aspect = rect.width() / rect.height();
    let rect_center = rect.center();
    let half_width = rect.width() / 2.0;
    let half_height = rect.height() / 2.0;

    let view_transform = match crate::render_cache::ViewTransform::from_camera(camera) {
        Some(vt) => vt,
        None => return,
    };

    let highlight_color = egui::Color32::from_rgba_unmultiplied(255, 165, 0, 180); // Orange with transparency
    let outline_color = egui::Color32::from_rgb(255, 200, 0);

    for &(elem_id, face_idx) in &selection.selected_faces {
        // Find this face in the cache
        for face in &cache.faces {
            if face.element_id == elem_id && face.face_index == face_idx {
                // Project corners
                let proj: Vec<Option<egui::Pos2>> = face
                    .corners
                    .iter()
                    .map(|c| {
                        view_transform.project(*c, rect_center, half_width, half_height, aspect)
                    })
                    .collect();

                if proj.iter().any(|p| p.is_none()) {
                    continue;
                }

                let proj: Vec<egui::Pos2> = proj.into_iter().map(|p| p.unwrap()).collect();

                // Draw filled quad
                painter.add(egui::Shape::convex_polygon(
                    proj.clone(),
                    highlight_color,
                    egui::Stroke::new(2.0_f32, outline_color),
                ));

                break;
            }
        }
    }
}

/// Show context menu on right-click
fn show_viewport_context_menu(response: &egui::Response, app: &mut FeaApp) {
    response.context_menu(|ui| {
        ui.set_min_width(180.0);

        // View options
        ui.menu_button("📷 View", |ui| {
            if ui.button("Front").clicked() {
                app.state.ui_state.camera.set_front_view();
                ui.close_menu();
            }
            if ui.button("Side").clicked() {
                app.state.ui_state.camera.set_side_view();
                ui.close_menu();
            }
            if ui.button("Top").clicked() {
                app.state.ui_state.camera.set_top_view();
                ui.close_menu();
            }
            if ui.button("Isometric").clicked() {
                app.state.ui_state.camera.set_iso_view();
                ui.close_menu();
            }
            ui.separator();
            if ui.button("🎯 Fit to Mesh").clicked() {
                if let Some(mesh) = app.state.current_mesh() {
                    let bounds = mesh.bounds;
                    app.state.ui_state.camera.fit_to_bounds(&bounds);
                }
                ui.close_menu();
            }
        });

        ui.separator();

        // Display toggles
        ui.menu_button("👁 Display", |ui| {
            ui.checkbox(&mut app.state.ui_state.show_faces, "Faces");
            ui.checkbox(&mut app.state.ui_state.show_wireframe, "Wireframe");
            ui.checkbox(&mut app.state.ui_state.show_nodes, "Nodes");
            ui.checkbox(
                &mut app.state.ui_state.show_boundary_conditions,
                "Boundary Conditions",
            );
            ui.separator();
            ui.checkbox(
                &mut app.state.ui_state.clipping_plane.enabled,
                "Clipping Plane",
            );
            ui.checkbox(&mut app.state.ui_state.display_settings.show_grid, "Grid");
        });

        ui.separator();

        // Selection tools
        if app.state.current_mesh().is_some() {
            ui.menu_button("🎯 Selection", |ui| {
                let was_active = app.state.ui_state.mesh_edit.face_selection.active;
                if ui
                    .checkbox(
                        &mut app.state.ui_state.mesh_edit.face_selection.active,
                        "Face Selection Mode",
                    )
                    .changed()
                {
                    if app.state.ui_state.mesh_edit.face_selection.active && !was_active {
                        app.state.status_message = "Click faces to select".to_string();
                    }
                }

                if app.state.ui_state.mesh_edit.face_selection.active {
                    ui.separator();
                    if ui.button("Select All Faces").clicked() {
                        // Would need to implement
                        app.state.status_message = "Select all: not yet implemented".to_string();
                        ui.close_menu();
                    }
                    if ui.button("Clear Selection").clicked() {
                        app.state
                            .ui_state
                            .mesh_edit
                            .face_selection
                            .selected_faces
                            .clear();
                        ui.close_menu();
                    }
                }
            });

            // Quick actions if faces are selected
            let selected_count = app
                .state
                .ui_state
                .mesh_edit
                .face_selection
                .selected_faces
                .len();
            if selected_count > 0 {
                ui.separator();
                ui.label(format!("{} faces selected", selected_count));

                if ui.button("📍 Create Node Group").clicked() {
                    // Trigger node group creation from faces
                    app.state.ui_state.mesh_edit.face_selection.purpose =
                        crate::state::SelectionPurpose::CreateNodeGroup;
                    ui.close_menu();
                }
                if ui.button("Apply Fixed BC").clicked() {
                    app.state.ui_state.mesh_edit.face_selection.purpose =
                        crate::state::SelectionPurpose::ApplyBC;
                    ui.close_menu();
                }
            }

            ui.separator();

            // Mesh operations
            ui.menu_button("Mesh", |ui| {
                if ui.button("Transform...").clicked() {
                    app.state.ui_state.mesh_edit.transform.mode =
                        crate::state::TransformMode::Translate;
                    ui.close_menu();
                }
                if ui.button("+ Create Primitive...").clicked() {
                    app.state.ui_state.primitive_dialog_open = true;
                    ui.close_menu();
                }
            });
        } else {
            // No mesh loaded
            if ui.button("+ Create Primitive Mesh").clicked() {
                app.state.ui_state.primitive_dialog_open = true;
                ui.close_menu();
            }
            if ui.button("Import Mesh").clicked() {
                // Would need to trigger file dialog
                app.state.status_message = "Use File > Import Mesh".to_string();
                ui.close_menu();
            }
        }
    });
}

fn handle_camera_input(response: &egui::Response, app: &mut FeaApp) {
    let camera = &mut app.state.ui_state.camera;

    // Drag to rotate (only when not in face selection mode or no click)
    if response.dragged_by(egui::PointerButton::Primary) {
        let delta = response.drag_delta();
        camera.yaw += delta.x * 0.01;
        camera.pitch = (camera.pitch + delta.y * 0.01).clamp(
            -std::f32::consts::FRAC_PI_2 + 0.01,
            std::f32::consts::FRAC_PI_2 - 0.01,
        );
    }

    // Right-drag to pan (context menu handles right-click without drag)
    if response.dragged_by(egui::PointerButton::Secondary) {
        let delta = response.drag_delta();
        let pan_speed = camera.distance * 0.002;

        // Pan in screen space
        let cos_yaw = camera.yaw.cos();
        let sin_yaw = camera.yaw.sin();

        camera.target[0] -=
            (cos_yaw * delta.x - sin_yaw * delta.y * camera.pitch.sin()) * pan_speed;
        camera.target[1] += delta.y * camera.pitch.cos() * pan_speed;
        camera.target[2] -=
            (sin_yaw * delta.x + cos_yaw * delta.y * camera.pitch.sin()) * pan_speed;
    }

    // Scroll to zoom
    let scroll = response.ctx.input(|i| i.raw_scroll_delta.y);
    if scroll != 0.0 {
        let zoom_factor = 1.0 - scroll * 0.001;
        camera.distance = (camera.distance * zoom_factor).clamp(0.01, 1000.0);
    }
}

/// Clip a 3D line segment to the camera's near plane and project both endpoints.
/// Returns None if the entire segment is behind the camera.
fn clip_and_project_line(
    p1_3d: [f32; 3],
    p2_3d: [f32; 3],
    rect: egui::Rect,
    camera: &crate::state::CameraState,
) -> Option<(egui::Pos2, egui::Pos2)> {
    // For orthographic projection, no near-plane clipping needed
    if camera.orthographic {
        let proj1 = project_point(p1_3d, rect, camera)?;
        let proj2 = project_point(p2_3d, rect, camera)?;
        return Some((proj1, proj2));
    }

    let eye = camera.eye_position();

    // View direction (normalized)
    let view_x = camera.target[0] - eye[0];
    let view_y = camera.target[1] - eye[1];
    let view_z = camera.target[2] - eye[2];
    let view_len = (view_x * view_x + view_y * view_y + view_z * view_z).sqrt();
    if view_len < 1e-6 {
        return None;
    }
    let forward = [view_x / view_len, view_y / view_len, view_z / view_len];

    // Use a very small near plane for grid lines - they don't have depth-sorting issues
    // This allows the grid to remain visible even when zoomed in close
    let near_plane = 0.001;

    // Compute distance along view direction for each point
    let rel1 = [p1_3d[0] - eye[0], p1_3d[1] - eye[1], p1_3d[2] - eye[2]];
    let rel2 = [p2_3d[0] - eye[0], p2_3d[1] - eye[1], p2_3d[2] - eye[2]];

    let d1 = rel1[0] * forward[0] + rel1[1] * forward[1] + rel1[2] * forward[2];
    let d2 = rel2[0] * forward[0] + rel2[1] * forward[1] + rel2[2] * forward[2];

    // Both behind near plane - reject
    if d1 < near_plane && d2 < near_plane {
        return None;
    }

    // Clip the segment to the near plane if one point is behind
    let (clipped_p1, clipped_p2) = if d1 < near_plane {
        // p1 is behind, clip it
        let t = (near_plane - d1) / (d2 - d1);
        let new_p1 = [
            p1_3d[0] + t * (p2_3d[0] - p1_3d[0]),
            p1_3d[1] + t * (p2_3d[1] - p1_3d[1]),
            p1_3d[2] + t * (p2_3d[2] - p1_3d[2]),
        ];
        (new_p1, p2_3d)
    } else if d2 < near_plane {
        // p2 is behind, clip it
        let t = (near_plane - d2) / (d1 - d2);
        let new_p2 = [
            p2_3d[0] + t * (p1_3d[0] - p2_3d[0]),
            p2_3d[1] + t * (p1_3d[1] - p2_3d[1]),
            p2_3d[2] + t * (p1_3d[2] - p2_3d[2]),
        ];
        (p1_3d, new_p2)
    } else {
        (p1_3d, p2_3d)
    };

    // Project both points with a tiny near plane for grid rendering
    // We need a custom projection here that doesn't use project_point's larger near plane
    let fov_factor = (camera.fov / 2.0).tan();
    let aspect = rect.width() / rect.height();

    let project_grid_point = |point: [f32; 3]| -> Option<egui::Pos2> {
        let rel = [point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]];

        // Calculate camera basis vectors (same as forward calculation above)
        let world_up = if forward[1].abs() > 0.99 {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let right = [
            forward[1] * world_up[2] - forward[2] * world_up[1],
            forward[2] * world_up[0] - forward[0] * world_up[2],
            forward[0] * world_up[1] - forward[1] * world_up[0],
        ];
        let right_len = (right[0] * right[0] + right[1] * right[1] + right[2] * right[2]).sqrt();
        if right_len < 1e-6 {
            return None;
        }
        let right = [
            right[0] / right_len,
            right[1] / right_len,
            right[2] / right_len,
        ];
        let up = [
            right[1] * forward[2] - right[2] * forward[1],
            right[2] * forward[0] - right[0] * forward[2],
            right[0] * forward[1] - right[1] * forward[0],
        ];

        let cam_x = rel[0] * right[0] + rel[1] * right[1] + rel[2] * right[2];
        let cam_y = rel[0] * up[0] + rel[1] * up[1] + rel[2] * up[2];
        let cam_z = rel[0] * forward[0] + rel[1] * forward[1] + rel[2] * forward[2];

        // Use the tiny near plane we defined above
        if cam_z < near_plane {
            return None;
        }

        let ndc_x = cam_x / (cam_z * fov_factor * aspect);
        let ndc_y = cam_y / (cam_z * fov_factor);

        let screen_x = rect.center().x + ndc_x * rect.width() / 2.0;
        let screen_y = rect.center().y - ndc_y * rect.height() / 2.0;

        if !screen_x.is_finite() || !screen_y.is_finite() {
            return None;
        }

        Some(egui::pos2(screen_x, screen_y))
    };

    let proj1 = project_grid_point(clipped_p1)?;
    let proj2 = project_grid_point(clipped_p2)?;

    Some((proj1, proj2))
}

fn draw_grid(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    // Check if grid is enabled
    if !app.state.ui_state.display_settings.show_grid {
        return;
    }

    let camera = &app.state.ui_state.camera;

    // Disable grid when zoomed in too close (avoids rendering artifacts)
    if camera.distance < 0.5 {
        return;
    }

    let settings = &app.state.ui_state.display_settings;

    // Use display settings for grid
    let grid_size = settings.grid_size;
    let grid_spacing = settings.grid_spacing;
    let grid_color = egui::Color32::from_gray(50);

    // Expand rect for clipping
    let clip_rect = rect.expand(2.0);

    // Simple 2D line clipping (Cohen-Sutherland)
    let clip_line = |mut p1: egui::Pos2, mut p2: egui::Pos2| -> Option<(egui::Pos2, egui::Pos2)> {
        const INSIDE: u8 = 0;
        const LEFT: u8 = 1;
        const RIGHT: u8 = 2;
        const BOTTOM: u8 = 4;
        const TOP: u8 = 8;

        let outcode = |p: egui::Pos2| -> u8 {
            let mut code = INSIDE;
            if p.x < clip_rect.left() {
                code |= LEFT;
            } else if p.x > clip_rect.right() {
                code |= RIGHT;
            }
            if p.y < clip_rect.top() {
                code |= TOP;
            } else if p.y > clip_rect.bottom() {
                code |= BOTTOM;
            }
            code
        };

        let mut code1 = outcode(p1);
        let mut code2 = outcode(p2);

        for _ in 0..10 {
            // Max iterations to prevent infinite loop
            if (code1 | code2) == 0 {
                return Some((p1, p2));
            } else if (code1 & code2) != 0 {
                return None;
            } else {
                let code_out = if code1 != 0 { code1 } else { code2 };
                let dx = p2.x - p1.x;
                let dy = p2.y - p1.y;

                let (x, y) = if (code_out & TOP) != 0 {
                    (p1.x + dx * (clip_rect.top() - p1.y) / dy, clip_rect.top())
                } else if (code_out & BOTTOM) != 0 {
                    (
                        p1.x + dx * (clip_rect.bottom() - p1.y) / dy,
                        clip_rect.bottom(),
                    )
                } else if (code_out & RIGHT) != 0 {
                    (
                        clip_rect.right(),
                        p1.y + dy * (clip_rect.right() - p1.x) / dx,
                    )
                } else {
                    (clip_rect.left(), p1.y + dy * (clip_rect.left() - p1.x) / dx)
                };

                if !x.is_finite() || !y.is_finite() {
                    return None;
                }

                if code_out == code1 {
                    p1 = egui::pos2(x, y);
                    code1 = outcode(p1);
                } else {
                    p2 = egui::pos2(x, y);
                    code2 = outcode(p2);
                }
            }
        }
        None
    };

    // Direct projection for grid - no near plane culling at all
    let project_grid = |point: [f32; 3]| -> Option<egui::Pos2> {
        let eye = camera.eye_position();

        let view_x = camera.target[0] - eye[0];
        let view_y = camera.target[1] - eye[1];
        let view_z = camera.target[2] - eye[2];
        let view_len = (view_x * view_x + view_y * view_y + view_z * view_z).sqrt();
        if view_len < 1e-6 {
            return None;
        }

        let forward = [view_x / view_len, view_y / view_len, view_z / view_len];

        let world_up = if forward[1].abs() > 0.99 {
            [0.0, 0.0, 1.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let right = [
            forward[1] * world_up[2] - forward[2] * world_up[1],
            forward[2] * world_up[0] - forward[0] * world_up[2],
            forward[0] * world_up[1] - forward[1] * world_up[0],
        ];
        let right_len = (right[0] * right[0] + right[1] * right[1] + right[2] * right[2]).sqrt();
        if right_len < 1e-6 {
            return None;
        }
        let right = [
            right[0] / right_len,
            right[1] / right_len,
            right[2] / right_len,
        ];
        let up = [
            right[1] * forward[2] - right[2] * forward[1],
            right[2] * forward[0] - right[0] * forward[2],
            right[0] * forward[1] - right[1] * forward[0],
        ];

        let rel = [point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]];
        let cam_x = rel[0] * right[0] + rel[1] * right[1] + rel[2] * right[2];
        let cam_y = rel[0] * up[0] + rel[1] * up[1] + rel[2] * up[2];
        let cam_z = rel[0] * forward[0] + rel[1] * forward[1] + rel[2] * forward[2];

        // Only reject if actually behind camera (z <= 0), no near plane
        if cam_z <= 0.0 {
            return None;
        }

        let fov_factor = (camera.fov / 2.0).tan();
        let aspect = rect.width() / rect.height();

        let (ndc_x, ndc_y) = if camera.orthographic {
            let ortho_scale = 1.0 / (camera.distance * fov_factor);
            (cam_x * ortho_scale / aspect, cam_y * ortho_scale)
        } else {
            (
                cam_x / (cam_z * fov_factor * aspect),
                cam_y / (cam_z * fov_factor),
            )
        };

        let screen_x = rect.center().x + ndc_x * rect.width() / 2.0;
        let screen_y = rect.center().y - ndc_y * rect.height() / 2.0;

        if !screen_x.is_finite() || !screen_y.is_finite() {
            return None;
        }

        Some(egui::pos2(screen_x, screen_y))
    };

    for i in -grid_size..=grid_size {
        let x = i as f32 * grid_spacing;

        // Line along X axis (constant X, varies in Z)
        let p1_3d = [x, 0.0, -grid_size as f32 * grid_spacing];
        let p2_3d = [x, 0.0, grid_size as f32 * grid_spacing];

        if let (Some(p1), Some(p2)) = (project_grid(p1_3d), project_grid(p2_3d)) {
            if let Some((cp1, cp2)) = clip_line(p1, p2) {
                painter.line_segment([cp1, cp2], egui::Stroke::new(1.0_f32, grid_color));
            }
        }

        // Line along Z axis (constant Z, varies in X)
        let p1_3d = [-grid_size as f32 * grid_spacing, 0.0, x];
        let p2_3d = [grid_size as f32 * grid_spacing, 0.0, x];

        if let (Some(p1), Some(p2)) = (project_grid(p1_3d), project_grid(p2_3d)) {
            if let Some((cp1, cp2)) = clip_line(p1, p2) {
                painter.line_segment([cp1, cp2], egui::Stroke::new(1.0_f32, grid_color));
            }
        }
    }
}

fn draw_mesh_cached(painter: &egui::Painter, rect: egui::Rect, app: &mut FeaApp, mesh_idx: usize) {
    // Get all the data we need from app first
    let mesh_state = &app.state.meshes[mesh_idx];
    let camera = app.state.ui_state.camera.clone();
    let ui_state = &app.state.ui_state;

    // Get displacement scale and results
    let scale = ui_state.displacement_scale;
    let color_mode = ui_state.color_mode;
    let stress_component = ui_state.stress_component;
    let strain_component = ui_state.strain_component;
    let show_faces = ui_state.show_faces;
    let show_wireframe = ui_state.show_wireframe;
    let show_nodes = ui_state.show_nodes;

    // Get camera eye position for backface culling
    let eye = camera.eye_position();
    let aspect = rect.width() / rect.height();

    // === OPTIMIZATION: Cached view transform ===
    // Update view transform cache if camera changed
    let view_transform = match render_cache::ViewTransform::from_camera(&camera) {
        Some(vt) => vt,
        None => return, // Degenerate camera, skip rendering
    };

    // Pre-compute rect values for projection
    let rect_center = rect.center();
    let half_width = rect.width() / 2.0;
    let half_height = rect.height() / 2.0;

    // === OPTIMIZATION: Use render cache ===
    // Check if cache needs rebuild
    let results = &app.state.results;
    let cache = &mut app.render_cache;
    if !cache.is_valid(
        mesh_state,
        results,
        scale,
        color_mode,
        stress_component,
        strain_component,
    ) {
        cache.rebuild(
            mesh_state,
            results,
            scale,
            color_mode,
            stress_component,
            strain_component,
        );
    }
    cache.view_transform = Some(view_transform.clone());

    // Get visible faces using frustum culling
    let visible_face_indices = cache.get_visible_faces(&camera, aspect);

    // === OPTIMIZATION: LOD for very large meshes ===
    // For very large face counts, only render every Nth face (still looks good from distance)
    let total_faces = cache.faces.len();
    let face_stride = if total_faces > 50000 {
        4 // Show 25% of faces
    } else if total_faces > 20000 {
        2 // Show 50% of faces
    } else {
        1 // Show all faces
    };

    // Also reduce triangle validation overhead for large meshes
    let fast_mode = total_faces > 10000;

    // === CLIPPING PLANE SETUP ===
    let clip_enabled = ui_state.clipping_plane.enabled;
    let clip_normal = if ui_state.clipping_plane.flip {
        [
            -ui_state.clipping_plane.normal[0],
            -ui_state.clipping_plane.normal[1],
            -ui_state.clipping_plane.normal[2],
        ]
    } else {
        ui_state.clipping_plane.normal
    };

    // Calculate clip plane position based on mesh bounds
    let clip_pos = if clip_enabled {
        let bounds = &mesh_state.bounds;
        let center = bounds.center();
        let half_diag = bounds.diagonal() * 0.5;
        // Position ranges from -1 (one side of mesh) to +1 (other side)
        let offset = ui_state.clipping_plane.position * half_diag;

        // Clip plane distance from origin (for dot product test)
        center[0] * clip_normal[0]
            + center[1] * clip_normal[1]
            + center[2] * clip_normal[2]
            + offset
    } else {
        0.0
    };

    // === OPTIMIZATION: Pre-allocate shape vectors with estimated capacity ===
    // Use a simpler shape representation for better cache performance
    struct ProjectedFace {
        proj: [egui::Pos2; 4],
        depth: f32,
        color: egui::Color32,
        is_backface: bool,
        /// Per-corner colors for smooth shading (computed via shape functions)
        corner_colors: [egui::Color32; 4],
        /// Whether to use smooth shading (subdivide into triangles with interpolated colors)
        use_smooth_shading: bool,
    }

    let estimated_visible = visible_face_indices.len() / face_stride;
    let mut projected_faces: Vec<ProjectedFace> = Vec::with_capacity(estimated_visible);

    // === OPTIMIZATION: Batch process faces ===
    let max_screen_dist = rect.width().max(rect.height()) * 1.5;
    let expanded_rect = rect.expand(200.0);

    for (i, &face_idx) in visible_face_indices.iter().enumerate() {
        // Skip faces for LOD
        if i % face_stride != 0 {
            continue;
        }

        let face = &cache.faces[face_idx];

        // === CLIPPING PLANE CHECK ===
        // Skip faces that are entirely on the clipped side
        if clip_enabled {
            let face_dist = face.center[0] * clip_normal[0]
                + face.center[1] * clip_normal[1]
                + face.center[2] * clip_normal[2];
            if face_dist < clip_pos {
                continue; // Face center is on clipped side
            }
        }

        // === OPTIMIZATION: Early backface culling before projection ===
        let view_dir = [
            face.center[0] - eye[0],
            face.center[1] - eye[1],
            face.center[2] - eye[2],
        ];
        let dot = face.normal[0] * view_dir[0]
            + face.normal[1] * view_dir[1]
            + face.normal[2] * view_dir[2];
        let is_backface = dot > 0.0;

        // Skip backfaces early when wireframe is on (they'd be hidden anyway)
        if is_backface && show_wireframe && show_faces {
            continue;
        }

        // Project corners using cached view transform
        let proj = [
            view_transform.project(
                face.corners[0],
                rect_center,
                half_width,
                half_height,
                aspect,
            ),
            view_transform.project(
                face.corners[1],
                rect_center,
                half_width,
                half_height,
                aspect,
            ),
            view_transform.project(
                face.corners[2],
                rect_center,
                half_width,
                half_height,
                aspect,
            ),
            view_transform.project(
                face.corners[3],
                rect_center,
                half_width,
                half_height,
                aspect,
            ),
        ];

        // Skip if any point is behind camera
        if proj.iter().any(|p| p.is_none()) {
            continue;
        }

        let proj: [egui::Pos2; 4] = [
            proj[0].unwrap(),
            proj[1].unwrap(),
            proj[2].unwrap(),
            proj[3].unwrap(),
        ];

        // === OPTIMIZATION: Combined validation checks ===
        let mut valid = true;
        for p in &proj {
            if !expanded_rect.contains(*p) {
                valid = false;
                break;
            }
        }
        if !valid {
            continue;
        }

        // Check edge lengths (simplified - just check longest edges)
        let dx02 = proj[0].x - proj[2].x;
        let dy02 = proj[0].y - proj[2].y;
        let dx13 = proj[1].x - proj[3].x;
        let dy13 = proj[1].y - proj[3].y;
        let diag1_sq = dx02 * dx02 + dy02 * dy02;
        let diag2_sq = dx13 * dx13 + dy13 * dy13;

        if diag1_sq > max_screen_dist * max_screen_dist
            || diag2_sq > max_screen_dist * max_screen_dist
            || !diag1_sq.is_finite()
            || !diag2_sq.is_finite()
        {
            continue;
        }

        // Calculate depth
        let z0 = camera_space_z(face.corners[0], &camera);
        let z1 = camera_space_z(face.corners[1], &camera);
        let z2 = camera_space_z(face.corners[2], &camera);
        let z3 = camera_space_z(face.corners[3], &camera);
        let min_z = z0.min(z1).min(z2).min(z3);
        let center_z = (z0 + z1 + z2 + z3) * 0.25;
        let depth = center_z * 0.6 + min_z * 0.4;

        // Get color from cache (used for flat shading fallback in solid mode)
        let base_color = cache.get_face_color(face, color_mode);

        // Compute per-corner colors for smooth shading with lighting
        let (corner_colors, use_smooth_shading, color) = if color_mode != ColorMode::Solid {
            // Get colors for each corner using nodal values
            let c0 = cache.get_node_color(face.node_ids[0], face.element_id, color_mode);
            let c1 = cache.get_node_color(face.node_ids[1], face.element_id, color_mode);
            let c2 = cache.get_node_color(face.node_ids[2], face.element_id, color_mode);
            let c3 = cache.get_node_color(face.node_ids[3], face.element_id, color_mode);

            // Compute lighting based on face normal
            // Light comes from upper-front-right relative to view
            let light_dir = [0.3_f32, 0.5, 0.8]; // Normalized direction toward light
            let light_len = (light_dir[0] * light_dir[0]
                + light_dir[1] * light_dir[1]
                + light_dir[2] * light_dir[2])
                .sqrt();
            let light_dir = [
                light_dir[0] / light_len,
                light_dir[1] / light_len,
                light_dir[2] / light_len,
            ];

            // Normalize face normal
            let n_len = (face.normal[0] * face.normal[0]
                + face.normal[1] * face.normal[1]
                + face.normal[2] * face.normal[2])
                .sqrt();
            let n = if n_len > 0.0001 {
                [
                    face.normal[0] / n_len,
                    face.normal[1] / n_len,
                    face.normal[2] / n_len,
                ]
            } else {
                [0.0, 0.0, 1.0]
            };

            // Diffuse lighting (N dot L), clamped
            let n_dot_l = n[0] * light_dir[0] + n[1] * light_dir[1] + n[2] * light_dir[2];
            let diffuse = n_dot_l.abs().clamp(0.0, 1.0); // Use abs for both face orientations

            // Combine ambient + diffuse lighting
            let ambient = 0.4_f32;
            let light_factor = ambient + (1.0 - ambient) * diffuse;

            // Apply backface darkening
            let final_factor = if is_backface {
                light_factor * 0.5
            } else {
                light_factor
            };

            (
                [
                    crate::render_cache::darken_color(c0, final_factor),
                    crate::render_cache::darken_color(c1, final_factor),
                    crate::render_cache::darken_color(c2, final_factor),
                    crate::render_cache::darken_color(c3, final_factor),
                ],
                true,
                base_color, // fallback color for ProjectedFace
            )
        } else {
            // Solid color mode: apply same lighting to flat color
            // Light direction
            let light_dir = [0.3_f32, 0.5, 0.8];
            let light_len = (light_dir[0] * light_dir[0]
                + light_dir[1] * light_dir[1]
                + light_dir[2] * light_dir[2])
                .sqrt();
            let light_dir = [
                light_dir[0] / light_len,
                light_dir[1] / light_len,
                light_dir[2] / light_len,
            ];

            let n_len = (face.normal[0] * face.normal[0]
                + face.normal[1] * face.normal[1]
                + face.normal[2] * face.normal[2])
                .sqrt();
            let n = if n_len > 0.0001 {
                [
                    face.normal[0] / n_len,
                    face.normal[1] / n_len,
                    face.normal[2] / n_len,
                ]
            } else {
                [0.0, 0.0, 1.0]
            };

            let n_dot_l = n[0] * light_dir[0] + n[1] * light_dir[1] + n[2] * light_dir[2];
            let diffuse = n_dot_l.abs().clamp(0.0, 1.0);
            let ambient = 0.4_f32;
            let light_factor = ambient + (1.0 - ambient) * diffuse;
            let final_factor = if is_backface {
                light_factor * 0.5
            } else {
                light_factor
            };

            let lit_color = crate::render_cache::darken_color(base_color, final_factor);
            (
                [lit_color, lit_color, lit_color, lit_color],
                false,
                lit_color,
            )
        };

        projected_faces.push(ProjectedFace {
            proj,
            depth,
            color,
            is_backface,
            corner_colors,
            use_smooth_shading,
        });
    }

    // === OPTIMIZATION: Partial sort instead of full sort for large meshes ===
    // For very large face counts, use a simpler approach
    if projected_faces.len() > 5000 {
        // Use unstable sort (faster, O(n) stack space instead of O(n) heap)
        projected_faces.sort_unstable_by(|a, b| {
            b.depth
                .partial_cmp(&a.depth)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    } else {
        projected_faces.sort_by(|a, b| {
            match b.depth.partial_cmp(&a.depth) {
                Some(std::cmp::Ordering::Equal) => {
                    // Backfaces first at same depth
                    b.is_backface.cmp(&a.is_backface)
                }
                Some(ord) => ord,
                None => std::cmp::Ordering::Equal,
            }
        });
    }

    // === SECTION CUT: Generate projected polygons for unified depth sorting ===
    // Compute forward vector for depth calculation
    let forward = [
        camera.target[0] - eye[0],
        camera.target[1] - eye[1],
        camera.target[2] - eye[2],
    ];

    let section_cut_shapes =
        if clip_enabled && show_faces && ui_state.clipping_plane.show_section_surface {
            // Compute section cuts if needed (uses cache)
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

            // Generate projected shapes with depth
            let camera_z_fn = |point: [f32; 3]| -> f32 {
                let rel = [point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]];
                rel[0] * forward[0] + rel[1] * forward[1] + rel[2] * forward[2]
            };

            section_cut::generate_section_cut_shapes(
                rect,
                &view_transform,
                &app.section_cut_cache.polygons,
                results,
                ui_state.color_mode,
                camera_z_fn,
            )
        } else {
            Vec::new()
        };

    // === UNIFIED DEPTH SORTING ===
    // Create a combined list of (depth, shape_index, is_section_cut) for proper interleaving
    struct DepthSortedItem {
        depth: f32,
        face_idx: Option<usize>, // Index into projected_faces
        cut_idx: Option<usize>,  // Index into section_cut_shapes
    }

    let mut depth_items: Vec<DepthSortedItem> =
        Vec::with_capacity(projected_faces.len() + section_cut_shapes.len());

    for (i, pf) in projected_faces.iter().enumerate() {
        depth_items.push(DepthSortedItem {
            depth: pf.depth,
            face_idx: Some(i),
            cut_idx: None,
        });
    }

    for (i, cut) in section_cut_shapes.iter().enumerate() {
        depth_items.push(DepthSortedItem {
            depth: cut.depth,
            face_idx: None,
            cut_idx: Some(i),
        });
    }

    // Sort all items by depth (back to front)
    depth_items.sort_unstable_by(|a, b| {
        b.depth
            .partial_cmp(&a.depth)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // === OPTIMIZATION: Batch draw calls ===
    // Collect all triangles first, then all lines

    // Pre-collect shapes for batching
    let mut triangle_shapes: Vec<egui::Shape> =
        Vec::with_capacity(projected_faces.len() + section_cut_shapes.len());
    let mut line_shapes: Vec<egui::Shape> = Vec::with_capacity(if show_wireframe {
        projected_faces.len() * 4
    } else {
        0
    });

    let _min_edge = 2.0;
    let wire_color = egui::Color32::from_gray(40);
    let max_line_len = rect.width().min(rect.height()) * 0.5;
    let line_margin_rect = rect.expand(100.0);

    // Helper to interpolate between two colors (linear RGB for better gradients)
    let lerp_color = |c0: egui::Color32, c1: egui::Color32, t: f32| -> egui::Color32 {
        // Use linear interpolation in sRGB space (good enough for visualization)
        let r = (c0.r() as f32 * (1.0 - t) + c1.r() as f32 * t) as u8;
        let g = (c0.g() as f32 * (1.0 - t) + c1.g() as f32 * t) as u8;
        let b = (c0.b() as f32 * (1.0 - t) + c1.b() as f32 * t) as u8;
        egui::Color32::from_rgb(r, g, b)
    };

    // Helper to bilinearly interpolate color on quad
    let bilinear_color = |c: &[egui::Color32; 4], u: f32, v: f32| -> egui::Color32 {
        // c[0] = (0,0), c[1] = (1,0), c[2] = (1,1), c[3] = (0,1)
        let c_bottom = lerp_color(c[0], c[1], u);
        let c_top = lerp_color(c[3], c[2], u);
        lerp_color(c_bottom, c_top, v)
    };

    // Helper to bilinearly interpolate position on quad
    let bilinear_pos = |p: &[egui::Pos2; 4], u: f32, v: f32| -> egui::Pos2 {
        let p_bottom = egui::pos2(
            p[0].x * (1.0 - u) + p[1].x * u,
            p[0].y * (1.0 - u) + p[1].y * u,
        );
        let p_top = egui::pos2(
            p[3].x * (1.0 - u) + p[2].x * u,
            p[3].y * (1.0 - u) + p[2].y * u,
        );
        egui::pos2(
            p_bottom.x * (1.0 - v) + p_top.x * v,
            p_bottom.y * (1.0 - v) + p_top.y * v,
        )
    };

    // Calculate adaptive subdivision based on screen-space size and color variance
    let calc_subdivision = |proj: &[egui::Pos2; 4], colors: &[egui::Color32; 4]| -> usize {
        // Base subdivision on screen-space diagonal
        let diag1 = ((proj[2].x - proj[0].x).powi(2) + (proj[2].y - proj[0].y).powi(2)).sqrt();
        let diag2 = ((proj[3].x - proj[1].x).powi(2) + (proj[3].y - proj[1].y).powi(2)).sqrt();
        let max_diag = diag1.max(diag2);

        // Calculate color variance to determine if more subdivision helps
        let color_diff = |c0: egui::Color32, c1: egui::Color32| -> f32 {
            let dr = (c0.r() as f32 - c1.r() as f32).abs();
            let dg = (c0.g() as f32 - c1.g() as f32).abs();
            let db = (c0.b() as f32 - c1.b() as f32).abs();
            (dr + dg + db) / 3.0
        };

        // Max color difference across the quad
        let max_color_diff = color_diff(colors[0], colors[1])
            .max(color_diff(colors[1], colors[2]))
            .max(color_diff(colors[2], colors[3]))
            .max(color_diff(colors[3], colors[0]))
            .max(color_diff(colors[0], colors[2]))
            .max(color_diff(colors[1], colors[3]));

        // More subdivisions for larger faces with higher color variance
        if max_diag < 20.0 || max_color_diff < 10.0 {
            2 // Small face or uniform color: 2x2 = 8 triangles
        } else if max_diag < 50.0 || max_color_diff < 30.0 {
            3 // Medium: 3x3 = 18 triangles
        } else if max_diag < 100.0 {
            4 // Large with variance: 4x4 = 32 triangles
        } else {
            5 // Very large: 5x5 = 50 triangles for smooth gradients
        }
    };

    // Render in depth-sorted order
    for item in &depth_items {
        if let Some(face_idx) = item.face_idx {
            let pf = &projected_faces[face_idx];

            // Draw faces
            if show_faces && (!pf.is_backface || !show_wireframe) {
                let area = quad_area_2d(&pf.proj);
                if area > 1.0 && area.is_finite() {
                    if pf.use_smooth_shading && !fast_mode {
                        // Smooth shading: subdivide quad into triangles with interpolated colors
                        // Use adaptive subdivision based on screen size and color variance
                        let subdiv = calc_subdivision(&pf.proj, &pf.corner_colors);
                        let step = 1.0 / subdiv as f32;

                        for j in 0..subdiv {
                            for i in 0..subdiv {
                                let u0 = i as f32 * step;
                                let v0 = j as f32 * step;
                                let u1 = u0 + step;
                                let v1 = v0 + step;

                                // Four corners of this sub-quad
                                let p00 = bilinear_pos(&pf.proj, u0, v0);
                                let p10 = bilinear_pos(&pf.proj, u1, v0);
                                let p01 = bilinear_pos(&pf.proj, u0, v1);
                                let p11 = bilinear_pos(&pf.proj, u1, v1);

                                // Use centroid colors for each triangle to reduce banding
                                // Triangle 1: p00 - p10 - p11 (centroid at (u0+2*u1)/3, (v0+2*v1)/3)
                                let uc1 = (u0 + u1 + u1) / 3.0;
                                let vc1 = (v0 + v0 + v1) / 3.0;
                                let c1 = bilinear_color(&pf.corner_colors, uc1, vc1);

                                triangle_shapes.push(egui::Shape::convex_polygon(
                                    vec![p00, p10, p11],
                                    c1,
                                    egui::Stroke::NONE,
                                ));

                                // Triangle 2: p00 - p11 - p01 (centroid at (u0+u0+u1)/3, (v0+v1+v1)/3)
                                let uc2 = (u0 + u0 + u1) / 3.0;
                                let vc2 = (v0 + v1 + v1) / 3.0;
                                let c2 = bilinear_color(&pf.corner_colors, uc2, vc2);

                                triangle_shapes.push(egui::Shape::convex_polygon(
                                    vec![p00, p11, p01],
                                    c2,
                                    egui::Stroke::NONE,
                                ));
                            }
                        }
                    } else {
                        // Flat shading or fast mode: render as single quad
                        triangle_shapes.push(egui::Shape::convex_polygon(
                            vec![pf.proj[0], pf.proj[1], pf.proj[2], pf.proj[3]],
                            pf.color,
                            egui::Stroke::NONE,
                        ));
                    }
                }
            }

            // Draw wireframe for front faces only
            // Skip wireframe if section cut surface is shown (it would draw over the cut surface)
            let skip_wireframe_for_section =
                clip_enabled && ui_state.clipping_plane.show_section_surface;
            if show_wireframe && !pf.is_backface && !skip_wireframe_for_section {
                let edges = [
                    [pf.proj[0], pf.proj[1]],
                    [pf.proj[1], pf.proj[2]],
                    [pf.proj[2], pf.proj[3]],
                    [pf.proj[3], pf.proj[0]],
                ];

                if fast_mode {
                    // Fast path: skip edge validation
                    for edge in &edges {
                        line_shapes.push(egui::Shape::line_segment(
                            *edge,
                            egui::Stroke::new(1.0_f32, wire_color),
                        ));
                    }
                } else {
                    // Careful path: validate each edge
                    for edge in &edges {
                        let dx = edge[1].x - edge[0].x;
                        let dy = edge[1].y - edge[0].y;
                        let len_sq = dx * dx + dy * dy;
                        if len_sq < max_line_len * max_line_len && len_sq.is_finite() {
                            if line_margin_rect.contains(edge[0])
                                && line_margin_rect.contains(edge[1])
                            {
                                line_shapes.push(egui::Shape::line_segment(
                                    *edge,
                                    egui::Stroke::new(1.0_f32, wire_color),
                                ));
                            }
                        }
                    }
                }
            }
        } else if let Some(cut_idx) = item.cut_idx {
            // Draw section cut polygon (no stroke for seamless appearance)
            let cut = &section_cut_shapes[cut_idx];
            triangle_shapes.push(egui::Shape::convex_polygon(
                cut.vertices.clone(),
                cut.color,
                egui::Stroke::NONE,
            ));
        }
    }

    // === OPTIMIZATION: Batch add shapes ===
    // Update render stats before drawing
    app.render_cache.last_rendered_faces = projected_faces.len();
    app.render_cache.last_rendered_triangles = triangle_shapes.len();

    painter.extend(triangle_shapes);
    painter.extend(line_shapes);

    // Section cuts are now integrated into the main depth-sorted pipeline above

    // Draw clipping plane outline if enabled
    if clip_enabled && ui_state.clipping_plane.show_plane {
        draw_clip_plane_outline(&painter, rect, &view_transform, mesh_state, ui_state);
    }

    // Draw nodes (with LOD - skip some for very large meshes)
    if show_nodes {
        let mesh = &app.state.meshes[mesh_idx].mesh;
        let node_count = mesh.nodes.len();

        // LOD: for large meshes, only draw every Nth node
        let node_stride = if node_count > 10000 {
            (node_count / 5000).max(1)
        } else if node_count > 5000 {
            2
        } else {
            1
        };

        let mut node_shapes: Vec<egui::Shape> = Vec::with_capacity(node_count / node_stride + 1);

        for (i, (_node_id, node)) in mesh.nodes.iter().enumerate() {
            if i % node_stride != 0 {
                continue;
            }

            let pos = [
                node.coordinates[0] as f32,
                node.coordinates[1] as f32,
                node.coordinates[2] as f32,
            ];

            if let Some(proj) =
                view_transform.project(pos, rect_center, half_width, half_height, aspect)
            {
                if rect.expand(10.0).contains(proj) {
                    node_shapes.push(egui::Shape::circle_filled(proj, 2.0, egui::Color32::YELLOW));
                }
            }
        }

        painter.extend(node_shapes);
    }
}

/// Draw the clipping plane outline as a translucent rectangle with dashed edges
fn draw_clip_plane_outline(
    painter: &egui::Painter,
    rect: egui::Rect,
    view_transform: &crate::render_cache::ViewTransform,
    mesh_state: &crate::state::MeshState,
    ui_state: &crate::state::UiState,
) {
    let bounds = &mesh_state.bounds;
    let center = bounds.center();
    let half_diag = bounds.diagonal() * 0.5;

    // Calculate clip plane position
    let clip_normal = if ui_state.clipping_plane.flip {
        [
            -ui_state.clipping_plane.normal[0],
            -ui_state.clipping_plane.normal[1],
            -ui_state.clipping_plane.normal[2],
        ]
    } else {
        ui_state.clipping_plane.normal
    };

    let offset = ui_state.clipping_plane.position * half_diag;

    // Compute plane center point
    let plane_center = [
        center[0] + clip_normal[0] * offset,
        center[1] + clip_normal[1] * offset,
        center[2] + clip_normal[2] * offset,
    ];

    // Create orthogonal vectors on the plane (u and v perpendicular to normal)
    let up = if clip_normal[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };

    // u = normal x up
    let mut u = [
        clip_normal[1] * up[2] - clip_normal[2] * up[1],
        clip_normal[2] * up[0] - clip_normal[0] * up[2],
        clip_normal[0] * up[1] - clip_normal[1] * up[0],
    ];
    let u_len = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    if u_len < 1e-6 {
        return;
    }
    u[0] /= u_len;
    u[1] /= u_len;
    u[2] /= u_len;

    // v = normal x u
    let v = [
        clip_normal[1] * u[2] - clip_normal[2] * u[1],
        clip_normal[2] * u[0] - clip_normal[0] * u[2],
        clip_normal[0] * u[1] - clip_normal[1] * u[0],
    ];

    // Compute plane size to cover the mesh bounds (use largest extent)
    let plane_size = half_diag * 1.2;

    // Compute the four corners of the plane rectangle
    let corners: [[f32; 3]; 4] = [
        [
            plane_center[0] - u[0] * plane_size - v[0] * plane_size,
            plane_center[1] - u[1] * plane_size - v[1] * plane_size,
            plane_center[2] - u[2] * plane_size - v[2] * plane_size,
        ],
        [
            plane_center[0] + u[0] * plane_size - v[0] * plane_size,
            plane_center[1] + u[1] * plane_size - v[1] * plane_size,
            plane_center[2] + u[2] * plane_size - v[2] * plane_size,
        ],
        [
            plane_center[0] + u[0] * plane_size + v[0] * plane_size,
            plane_center[1] + u[1] * plane_size + v[1] * plane_size,
            plane_center[2] + u[2] * plane_size + v[2] * plane_size,
        ],
        [
            plane_center[0] - u[0] * plane_size + v[0] * plane_size,
            plane_center[1] - u[1] * plane_size + v[1] * plane_size,
            plane_center[2] - u[2] * plane_size + v[2] * plane_size,
        ],
    ];

    // Project corners using project_close for better near-plane handling
    let rect_center = rect.center();
    let half_width = rect.width() * 0.5;
    let half_height = rect.height() * 0.5;
    let aspect = rect.width() / rect.height();

    let projected: Vec<_> = corners
        .iter()
        .filter_map(|&c| {
            view_transform.project_close(c, rect_center, half_width, half_height, aspect)
        })
        .collect();

    if projected.len() != 4 {
        return; // Some corners behind camera
    }

    // Check all corners are within reasonable bounds
    let expanded_rect = rect.expand(500.0);
    if !projected.iter().all(|p| expanded_rect.contains(*p)) {
        return;
    }

    // Draw semi-transparent fill
    let fill_color = egui::Color32::from_rgba_unmultiplied(100, 150, 255, 30);
    painter.add(egui::Shape::convex_polygon(
        projected.clone(),
        fill_color,
        egui::Stroke::NONE,
    ));

    // Draw dashed outline
    let outline_color = egui::Color32::from_rgba_unmultiplied(100, 150, 255, 180);
    let dash_length = 8.0;
    let gap_length = 4.0;

    for i in 0..4 {
        let p0 = projected[i];
        let p1 = projected[(i + 1) % 4];
        draw_dashed_line(painter, p0, p1, outline_color, 1.5, dash_length, gap_length);
    }
}

/// Draw a dashed line between two points
fn draw_dashed_line(
    painter: &egui::Painter,
    start: egui::Pos2,
    end: egui::Pos2,
    color: egui::Color32,
    width: f32,
    dash_length: f32,
    gap_length: f32,
) {
    let dx = end.x - start.x;
    let dy = end.y - start.y;
    let length = (dx * dx + dy * dy).sqrt();

    if length < 1.0 {
        return;
    }

    let dir_x = dx / length;
    let dir_y = dy / length;

    let mut dist = 0.0;
    let cycle_length = dash_length + gap_length;

    while dist < length {
        let dash_start = dist;
        let dash_end = (dist + dash_length).min(length);

        let p0 = egui::pos2(start.x + dir_x * dash_start, start.y + dir_y * dash_start);
        let p1 = egui::pos2(start.x + dir_x * dash_end, start.y + dir_y * dash_end);

        painter.line_segment([p0, p1], egui::Stroke::new(width, color));

        dist += cycle_length;
    }
}

/// Compute color for a node based on the current color mode
/// Note: This function is kept for potential future use but currently unused
/// as colors are computed in the render_cache module.
#[allow(dead_code)]
fn compute_node_color(
    node_id: usize,
    element_id: usize,
    color_mode: ColorMode,
    stress_component: usize,
    strain_component: usize,
    results: &Option<crate::state::SimulationResults>,
) -> egui::Color32 {
    match color_mode {
        ColorMode::Solid => egui::Color32::from_rgb(100, 149, 237), // Cornflower blue

        ColorMode::Displacement => {
            if let Some(res) = results {
                if node_id * 3 + 2 < res.displacements.len() {
                    let dx = res.displacements[node_id * 3];
                    let dy = res.displacements[node_id * 3 + 1];
                    let dz = res.displacements[node_id * 3 + 2];
                    let mag = (dx * dx + dy * dy + dz * dz).sqrt();
                    let t = (mag / res.stats.max_displacement.max(1e-10)) as f32;
                    return value_to_color(t.clamp(0.0, 1.0));
                }
            }
            egui::Color32::GRAY
        }

        ColorMode::VonMises => {
            if let Some(res) = results {
                if let Some(&vm) = res.von_mises.get(&element_id) {
                    let max_vm = res.stats.max_von_mises.max(1e-10);
                    let t = (vm / max_vm) as f32;
                    return value_to_color(t.clamp(0.0, 1.0));
                }
            }
            egui::Color32::GRAY
        }

        ColorMode::Stress => {
            if let Some(res) = results {
                if let Some(stress) = res.stresses.get(&element_id) {
                    // Use selected stress component
                    if stress_component < stress.len() {
                        let val = stress[stress_component];
                        // Find max for this component across all elements
                        let max_val = res
                            .stresses
                            .values()
                            .filter_map(|s| s.get(stress_component))
                            .map(|v| v.abs())
                            .fold(1e-10f64, |a, b| a.max(b));
                        // Use signed normalization for stress (can be + or -)
                        let t = ((val / max_val) * 0.5 + 0.5) as f32;
                        return value_to_color(t.clamp(0.0, 1.0));
                    }
                }
            }
            egui::Color32::GRAY
        }

        ColorMode::Strain => {
            if let Some(res) = results {
                if let Some(strain) = res.strains.get(&element_id) {
                    // Use selected strain component
                    if strain_component < strain.len() {
                        let val = strain[strain_component];
                        // Find max for this component across all elements
                        let max_val = res
                            .strains
                            .values()
                            .filter_map(|s| s.get(strain_component))
                            .map(|v| v.abs())
                            .fold(1e-10f64, |a, b| a.max(b));
                        // Use signed normalization for strain
                        let t = ((val / max_val) * 0.5 + 0.5) as f32;
                        return value_to_color(t.clamp(0.0, 1.0));
                    }
                }
            }
            egui::Color32::GRAY
        }
    }
}

/// Draw boundary condition visualizations
fn draw_boundary_conditions(
    painter: &egui::Painter,
    rect: egui::Rect,
    app: &FeaApp,
    mesh_state: &crate::state::MeshState,
) {
    let camera = &app.state.ui_state.camera;
    let mesh = &mesh_state.mesh;

    for bc in &app.state.simulation_config.boundary_conditions {
        match bc {
            BoundaryConditionConfig::Fixed(cfg) => {
                if let Some(group) = mesh.node_groups.get(&cfg.node_group) {
                    // Draw fixed constraints as small triangles/cones
                    let color = egui::Color32::from_rgb(255, 100, 100); // Red for fixed

                    for &node_id in &group.nodes {
                        if let Some(node) = mesh.nodes.get(&node_id) {
                            let pos = [
                                node.coordinates[0] as f32,
                                node.coordinates[1] as f32,
                                node.coordinates[2] as f32,
                            ];

                            if let Some(proj) = project_point(pos, rect, camera) {
                                // Only draw if within viewport with margin
                                if rect.expand(30.0).contains(proj) {
                                    // Draw small triangle pointing to node
                                    draw_fixed_symbol(painter, proj, color, cfg);
                                }
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::Load(cfg) => {
                if let Some(group) = mesh.node_groups.get(&cfg.node_group) {
                    // Draw load arrows
                    let color = egui::Color32::from_rgb(100, 255, 100); // Green for loads

                    // Normalize force direction
                    let force_mag =
                        (cfg.force_x.powi(2) + cfg.force_y.powi(2) + cfg.force_z.powi(2)).sqrt();
                    if force_mag > 1e-10 {
                        let dir = [
                            (cfg.force_x / force_mag) as f32,
                            (cfg.force_y / force_mag) as f32,
                            (cfg.force_z / force_mag) as f32,
                        ];

                        // Arrow length based on mesh scale
                        let arrow_len = mesh_state.bounds.diagonal() * 0.1;

                        for &node_id in &group.nodes {
                            if let Some(node) = mesh.nodes.get(&node_id) {
                                let pos = [
                                    node.coordinates[0] as f32,
                                    node.coordinates[1] as f32,
                                    node.coordinates[2] as f32,
                                ];

                                // Arrow start (offset from node in opposite direction of force)
                                let start = [
                                    pos[0] - dir[0] * arrow_len,
                                    pos[1] - dir[1] * arrow_len,
                                    pos[2] - dir[2] * arrow_len,
                                ];

                                if let (Some(p1), Some(p2)) = (
                                    project_point(start, rect, camera),
                                    project_point(pos, rect, camera),
                                ) {
                                    // Only draw if both points are within viewport with margin
                                    let margin = 50.0;
                                    let expanded_rect = rect.expand(margin);
                                    if expanded_rect.contains(p1) && expanded_rect.contains(p2) {
                                        draw_arrow_clamped(painter, p1, p2, color, 2.0, 100.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::Torque(cfg) => {
                if let Some(group) = mesh.node_groups.get(&cfg.node_group) {
                    // Draw torque as circular arrows
                    let color = egui::Color32::from_rgb(255, 200, 100); // Orange for torque

                    // Draw axis
                    let axis_len = mesh_state.bounds.diagonal() * 0.2;
                    let axis_start = [
                        cfg.axis_point[0] as f32,
                        cfg.axis_point[1] as f32,
                        cfg.axis_point[2] as f32,
                    ];
                    let axis_end = [
                        axis_start[0] + cfg.axis_direction[0] as f32 * axis_len,
                        axis_start[1] + cfg.axis_direction[1] as f32 * axis_len,
                        axis_start[2] + cfg.axis_direction[2] as f32 * axis_len,
                    ];

                    if let (Some(p1), Some(p2)) = (
                        project_point(axis_start, rect, camera),
                        project_point(axis_end, rect, camera),
                    ) {
                        // Only draw if within bounds
                        let margin = 100.0;
                        let expanded_rect = rect.expand(margin);
                        if expanded_rect.contains(p1) && expanded_rect.contains(p2) {
                            painter.line_segment([p1, p2], egui::Stroke::new(3.0_f32, color));

                            // Draw circular arrow indicator at midpoint
                            let mid = egui::pos2((p1.x + p2.x) / 2.0, (p1.y + p2.y) / 2.0);
                            painter.circle_stroke(mid, 10.0, egui::Stroke::new(2.0_f32, color));
                        }
                    }

                    // Highlight affected nodes
                    for &node_id in &group.nodes {
                        if let Some(node) = mesh.nodes.get(&node_id) {
                            let pos = [
                                node.coordinates[0] as f32,
                                node.coordinates[1] as f32,
                                node.coordinates[2] as f32,
                            ];

                            if let Some(proj) = project_point(pos, rect, camera) {
                                if rect.expand(20.0).contains(proj) {
                                    painter.circle_filled(proj, 3.0, color);
                                }
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::Contact(cfg) => {
                // Draw contact surfaces with different colors
                let primary_color = egui::Color32::from_rgb(255, 100, 255); // Magenta
                let secondary_color = egui::Color32::from_rgb(100, 255, 255); // Cyan

                // Primary surface
                if let Some(group) = mesh.node_groups.get(&cfg.primary_surface) {
                    for &node_id in &group.nodes {
                        if let Some(node) = mesh.nodes.get(&node_id) {
                            let pos = [
                                node.coordinates[0] as f32,
                                node.coordinates[1] as f32,
                                node.coordinates[2] as f32,
                            ];

                            if let Some(proj) = project_point(pos, rect, camera) {
                                if rect.expand(20.0).contains(proj) {
                                    painter.circle_filled(proj, 4.0, primary_color);
                                }
                            }
                        }
                    }
                }

                // Secondary surface
                if let Some(group) = mesh.node_groups.get(&cfg.secondary_surface) {
                    for &node_id in &group.nodes {
                        if let Some(node) = mesh.nodes.get(&node_id) {
                            let pos = [
                                node.coordinates[0] as f32,
                                node.coordinates[1] as f32,
                                node.coordinates[2] as f32,
                            ];

                            if let Some(proj) = project_point(pos, rect, camera) {
                                if rect.expand(20.0).contains(proj) {
                                    painter.circle_stroke(
                                        proj,
                                        4.0,
                                        egui::Stroke::new(2.0_f32, secondary_color),
                                    );
                                }
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::Pressure(cfg) => {
                // Draw pressure as blue arrows pointing into surface
                let color = egui::Color32::from_rgb(100, 150, 255); // Light blue

                if let Some(group) = mesh.element_groups.get(&cfg.element_group) {
                    // Show pressure indicators on element group nodes
                    for &el_id in &group.elements {
                        if let Some(element) = mesh.elements.get(&el_id) {
                            // Draw at element center
                            let center = element
                                .connectivity
                                .iter()
                                .filter_map(|&nid| mesh.nodes.get(&nid))
                                .fold([0.0f32, 0.0, 0.0], |acc, n| {
                                    [
                                        acc[0] + n.coordinates[0] as f32,
                                        acc[1] + n.coordinates[1] as f32,
                                        acc[2] + n.coordinates[2] as f32,
                                    ]
                                });
                            let n = element.connectivity.len() as f32;
                            let center = [center[0] / n, center[1] / n, center[2] / n];

                            if let Some(proj) = project_point(center, rect, camera) {
                                // Draw pressure symbol (inward pointing arrows)
                                let size = 6.0;
                                painter.circle_stroke(
                                    proj,
                                    size,
                                    egui::Stroke::new(2.0_f32, color),
                                );
                                // Inner arrows pointing in
                                for i in 0..4 {
                                    let angle = i as f32 * std::f32::consts::FRAC_PI_2;
                                    let outer = egui::pos2(
                                        proj.x + angle.cos() * size,
                                        proj.y + angle.sin() * size,
                                    );
                                    let inner = egui::pos2(
                                        proj.x + angle.cos() * (size * 0.3),
                                        proj.y + angle.sin() * (size * 0.3),
                                    );
                                    painter.line_segment(
                                        [outer, inner],
                                        egui::Stroke::new(1.5_f32, color),
                                    );
                                }
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::Traction(cfg) => {
                // Draw traction as arrows
                let color = egui::Color32::from_rgb(255, 180, 100); // Orange

                if let Some(group) = mesh.element_groups.get(&cfg.element_group) {
                    for &el_id in &group.elements {
                        if let Some(element) = mesh.elements.get(&el_id) {
                            // Draw at element center
                            let center = element
                                .connectivity
                                .iter()
                                .filter_map(|&nid| mesh.nodes.get(&nid))
                                .fold([0.0f32, 0.0, 0.0], |acc, n| {
                                    [
                                        acc[0] + n.coordinates[0] as f32,
                                        acc[1] + n.coordinates[1] as f32,
                                        acc[2] + n.coordinates[2] as f32,
                                    ]
                                });
                            let n = element.connectivity.len() as f32;
                            let center = [center[0] / n, center[1] / n, center[2] / n];

                            if let Some(proj) = project_point(center, rect, camera) {
                                // Draw traction symbol
                                let size = 8.0;
                                painter.rect_stroke(
                                    egui::Rect::from_center_size(proj, egui::vec2(size, size)),
                                    2.0,
                                    egui::Stroke::new(2.0_f32, color),
                                    egui::StrokeKind::Outside,
                                );
                            }
                        }
                    }
                }
            }

            BoundaryConditionConfig::BodyForce(cfg) => {
                // Draw body force as arrows at mesh center
                let color = egui::Color32::from_rgb(200, 100, 200); // Purple

                // Get direction based on body force type
                let dir = match &cfg.force_type {
                    BodyForceTypeConfig::Gravity { gx, gy, gz } => {
                        let mag = (*gx * *gx + *gy * *gy + *gz * *gz).sqrt();
                        if mag > 1e-10 {
                            [(*gx / mag) as f32, (*gy / mag) as f32, (*gz / mag) as f32]
                        } else {
                            [0.0, -1.0, 0.0]
                        }
                    }
                    BodyForceTypeConfig::Uniform { fx, fy, fz } => {
                        let mag = (*fx * *fx + *fy * *fy + *fz * *fz).sqrt();
                        if mag > 1e-10 {
                            [(*fx / mag) as f32, (*fy / mag) as f32, (*fz / mag) as f32]
                        } else {
                            [0.0, 0.0, 0.0]
                        }
                    }
                    BodyForceTypeConfig::Centrifugal { axis_direction, .. } => {
                        // Show axis direction
                        let mag = (axis_direction[0].powi(2)
                            + axis_direction[1].powi(2)
                            + axis_direction[2].powi(2))
                        .sqrt();
                        if mag > 1e-10 {
                            [
                                (axis_direction[0] / mag) as f32,
                                (axis_direction[1] / mag) as f32,
                                (axis_direction[2] / mag) as f32,
                            ]
                        } else {
                            [0.0, 1.0, 0.0]
                        }
                    }
                };

                // Draw gravity arrows at multiple points
                let arrow_len = mesh_state.bounds.diagonal() * 0.15;
                let bounds = &mesh_state.bounds;

                // Draw at 4 corners of bounding box
                let positions = [
                    [
                        bounds.min[0] as f32,
                        bounds.max[1] as f32,
                        bounds.min[2] as f32,
                    ],
                    [
                        bounds.max[0] as f32,
                        bounds.max[1] as f32,
                        bounds.min[2] as f32,
                    ],
                    [
                        bounds.min[0] as f32,
                        bounds.max[1] as f32,
                        bounds.max[2] as f32,
                    ],
                    [
                        bounds.max[0] as f32,
                        bounds.max[1] as f32,
                        bounds.max[2] as f32,
                    ],
                ];

                for pos in &positions {
                    let end = [
                        pos[0] + dir[0] * arrow_len,
                        pos[1] + dir[1] * arrow_len,
                        pos[2] + dir[2] * arrow_len,
                    ];

                    if let (Some(p1), Some(p2)) = (
                        project_point(*pos, rect, camera),
                        project_point(end, rect, camera),
                    ) {
                        // Only draw if both points are within viewport with margin
                        let margin = 50.0;
                        let expanded_rect = rect.expand(margin);
                        if expanded_rect.contains(p1) && expanded_rect.contains(p2) {
                            draw_arrow_clamped(painter, p1, p2, color, 2.0, 100.0);
                        }
                    }
                }

                // Label
                let center = [
                    (bounds.min[0] + bounds.max[0]) as f32 / 2.0,
                    bounds.max[1] as f32 + arrow_len * 0.5,
                    (bounds.min[2] + bounds.max[2]) as f32 / 2.0,
                ];
                if let Some(proj) = project_point(center, rect, camera) {
                    // Only draw label if within viewport
                    if rect.expand(50.0).contains(proj) {
                        let label = match &cfg.force_type {
                            BodyForceTypeConfig::Gravity { .. } => "g",
                            BodyForceTypeConfig::Centrifugal { .. } => "ω",
                            BodyForceTypeConfig::Uniform { .. } => "b",
                        };
                        painter.text(
                            proj,
                            egui::Align2::CENTER_CENTER,
                            label,
                            egui::FontId::proportional(14.0),
                            color,
                        );
                    }
                }
            }
        }
    }
}

/// Draw fixed constraint symbol
fn draw_fixed_symbol(
    painter: &egui::Painter,
    pos: egui::Pos2,
    color: egui::Color32,
    cfg: &crate::state::FixedBcConfig,
) {
    let size = 8.0;

    // Draw based on which DOFs are constrained
    if cfg.constrain_x.is_some() {
        // Vertical line for X constraint
        painter.line_segment(
            [
                egui::pos2(pos.x - size, pos.y),
                egui::pos2(pos.x - size, pos.y + size),
            ],
            egui::Stroke::new(2.0_f32, color),
        );
        // Ground hatch
        for i in 0..3 {
            let y = pos.y + i as f32 * 3.0;
            painter.line_segment(
                [
                    egui::pos2(pos.x - size - 4.0, y + 4.0),
                    egui::pos2(pos.x - size, y),
                ],
                egui::Stroke::new(1.0_f32, color),
            );
        }
    }

    if cfg.constrain_y.is_some() {
        // Horizontal line for Y constraint (below node)
        painter.line_segment(
            [
                egui::pos2(pos.x - size / 2.0, pos.y + size),
                egui::pos2(pos.x + size / 2.0, pos.y + size),
            ],
            egui::Stroke::new(2.0_f32, color),
        );
        // Ground hatch
        for i in 0..3 {
            let x = pos.x - size / 2.0 + i as f32 * 4.0;
            painter.line_segment(
                [
                    egui::pos2(x, pos.y + size + 4.0),
                    egui::pos2(x + 4.0, pos.y + size),
                ],
                egui::Stroke::new(1.0_f32, color),
            );
        }
    }

    if cfg.constrain_z.is_some() {
        // Circle for Z constraint
        painter.circle_stroke(pos, 3.0, egui::Stroke::new(2.0_f32, color));
    }

    // If all constrained, draw filled triangle
    if cfg.constrain_x.is_some() && cfg.constrain_y.is_some() && cfg.constrain_z.is_some() {
        let points = vec![
            pos,
            egui::pos2(pos.x - size / 2.0, pos.y + size),
            egui::pos2(pos.x + size / 2.0, pos.y + size),
        ];
        painter.add(egui::Shape::convex_polygon(
            points,
            color.linear_multiply(0.5),
            egui::Stroke::new(1.0_f32, color),
        ));
    }
}

/// Draw selected nodes for node group creation
fn draw_selected_nodes(
    painter: &egui::Painter,
    rect: egui::Rect,
    app: &FeaApp,
    mesh_state: &crate::state::MeshState,
) {
    let camera = &app.state.ui_state.camera;
    let mesh = &mesh_state.mesh;
    let color = egui::Color32::from_rgb(0, 255, 128); // Bright green

    for &node_id in &app.state.ui_state.node_group_creator.selected_nodes {
        if let Some(node) = mesh.nodes.get(&node_id) {
            let pos = [
                node.coordinates[0] as f32,
                node.coordinates[1] as f32,
                node.coordinates[2] as f32,
            ];

            if let Some(proj) = project_point(pos, rect, camera) {
                painter.circle_filled(proj, 4.0, color);
                painter.circle_stroke(proj, 6.0, egui::Stroke::new(1.0_f32, egui::Color32::WHITE));
            }
        }
    }
}

fn draw_axis_indicator(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    let camera = &app.state.ui_state.camera;

    // Position in bottom-left corner
    let origin = egui::pos2(rect.left() + 50.0, rect.bottom() - 50.0);
    let axis_length = 30.0;

    // Rotation matrix based on camera angles
    let cos_yaw = camera.yaw.cos();
    let sin_yaw = camera.yaw.sin();
    let cos_pitch = camera.pitch.cos();
    let sin_pitch = camera.pitch.sin();

    // X axis (red)
    let x_dir = [cos_yaw, 0.0, sin_yaw];
    let x_end = egui::pos2(
        origin.x + x_dir[0] * axis_length,
        origin.y - x_dir[2] * axis_length * cos_pitch - sin_pitch * axis_length * 0.3,
    );
    painter.arrow(
        origin,
        x_end - origin,
        egui::Stroke::new(2.0_f32, egui::Color32::RED),
    );
    painter.text(
        x_end,
        egui::Align2::CENTER_CENTER,
        "X",
        egui::FontId::proportional(12.0),
        egui::Color32::RED,
    );

    // Y axis (green)
    let y_end = egui::pos2(origin.x, origin.y - axis_length * cos_pitch);
    painter.arrow(
        origin,
        y_end - origin,
        egui::Stroke::new(2.0_f32, egui::Color32::GREEN),
    );
    painter.text(
        y_end,
        egui::Align2::CENTER_CENTER,
        "Y",
        egui::FontId::proportional(12.0),
        egui::Color32::GREEN,
    );

    // Z axis (blue)
    let z_dir = [-sin_yaw, 0.0, cos_yaw];
    let z_end = egui::pos2(
        origin.x + z_dir[0] * axis_length,
        origin.y - z_dir[2] * axis_length * cos_pitch - sin_pitch * axis_length * 0.3,
    );
    painter.arrow(
        origin,
        z_end - origin,
        egui::Stroke::new(2.0_f32, egui::Color32::BLUE),
    );
    painter.text(
        z_end,
        egui::Align2::CENTER_CENTER,
        "Z",
        egui::FontId::proportional(12.0),
        egui::Color32::BLUE,
    );
}

fn draw_view_info(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    // Draw help text in top-right
    let help_text = "LMB: Rotate | RMB: Pan | Scroll: Zoom";
    painter.text(
        egui::pos2(rect.right() - 10.0, rect.top() + 10.0),
        egui::Align2::RIGHT_TOP,
        help_text,
        egui::FontId::proportional(12.0),
        egui::Color32::from_gray(150),
    );

    // Draw stats overlay if enabled
    if app.state.ui_state.stats_overlay.visible {
        draw_stats_overlay(painter, rect, app);
    }

    // Draw mesh statistics in bottom-left (simple version, always shown)
    if !app.state.ui_state.stats_overlay.visible {
        if let Some(mesh_state) = app.state.current_mesh() {
            let stats_text = format!(
                "Nodes: {} | Elements: {}",
                mesh_state.mesh.nodes.len(),
                mesh_state.mesh.elements.len()
            );
            painter.text(
                egui::pos2(rect.left() + 10.0, rect.bottom() - 10.0),
                egui::Align2::LEFT_BOTTOM,
                stats_text,
                egui::FontId::proportional(12.0),
                egui::Color32::from_gray(150),
            );
        }
    }

    // Draw color legend when showing results
    if app.state.results.is_some() && app.state.ui_state.color_mode != ColorMode::Solid {
        draw_color_legend(painter, rect, app);
    }
}

/// Draw detailed statistics overlay
fn draw_stats_overlay(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    let overlay = &app.state.ui_state.stats_overlay;

    // Determine position based on setting
    let (pos, _anchor) = match overlay.position {
        0 => (
            egui::pos2(rect.left() + 10.0, rect.top() + 30.0),
            egui::Align2::LEFT_TOP,
        ),
        1 => (
            egui::pos2(rect.right() - 10.0, rect.top() + 30.0),
            egui::Align2::RIGHT_TOP,
        ),
        2 => (
            egui::pos2(rect.left() + 10.0, rect.bottom() - 10.0),
            egui::Align2::LEFT_BOTTOM,
        ),
        _ => (
            egui::pos2(rect.right() - 10.0, rect.bottom() - 10.0),
            egui::Align2::RIGHT_BOTTOM,
        ),
    };

    let mut lines: Vec<String> = Vec::new();

    // Mesh statistics
    if overlay.show_mesh_stats {
        if let Some(mesh_state) = app.state.current_mesh() {
            lines.push(format!("═══ Mesh: {} ═══", mesh_state.name));
            lines.push(format!(
                "  Nodes: {}",
                format_number(mesh_state.mesh.nodes.len())
            ));
            lines.push(format!(
                "  Elements: {}",
                format_number(mesh_state.mesh.elements.len())
            ));
            lines.push(format!(
                "  DOFs: {}",
                format_number(mesh_state.mesh.nodes.len() * 3)
            ));

            // Bounding box size
            let dx = mesh_state.bounds.max[0] - mesh_state.bounds.min[0];
            let dy = mesh_state.bounds.max[1] - mesh_state.bounds.min[1];
            let dz = mesh_state.bounds.max[2] - mesh_state.bounds.min[2];
            lines.push(format!("  Size: {:.3}×{:.3}×{:.3}", dx, dy, dz));

            // Node/element groups
            if !mesh_state.mesh.node_groups.is_empty() {
                lines.push(format!(
                    "  Node Groups: {}",
                    mesh_state.mesh.node_groups.len()
                ));
            }
            if !mesh_state.mesh.element_groups.is_empty() {
                lines.push(format!(
                    "  Element Groups: {}",
                    mesh_state.mesh.element_groups.len()
                ));
            }
        } else {
            lines.push("No mesh loaded".to_string());
        }
    }

    // Result statistics
    if overlay.show_result_stats {
        if let Some(results) = &app.state.results {
            lines.push("═══ Results ═══".to_string());
            lines.push(format!(
                "  Max Disp: {:.4e}",
                results.stats.max_displacement
            ));
            lines.push(format!("  Max σ_VM: {:.4e}", results.stats.max_von_mises));
            lines.push(format!(
                "  Solve time: {:.1}s",
                results.stats.solver_time_ms as f64 / 1000.0
            ));

            if !results.time_steps.is_empty() {
                lines.push(format!("  Time steps: {}", results.time_steps.len()));
            }
        }
    }

    // Camera info
    if overlay.show_camera_info {
        lines.push("═══ Camera ═══".to_string());
        let camera = &app.state.ui_state.camera;
        lines.push(format!("  Yaw: {:.1}°", camera.yaw.to_degrees()));
        lines.push(format!("  Pitch: {:.1}°", camera.pitch.to_degrees()));
        lines.push(format!("  Dist: {:.2}", camera.distance));
        lines.push(format!(
            "  Mode: {}",
            if camera.orthographic {
                "Ortho"
            } else {
                "Persp"
            }
        ));
    }

    // Performance stats
    if overlay.show_performance {
        lines.push("═══ Performance ═══".to_string());
        let fps = app.state.average_fps();
        lines.push(format!("  FPS: {:.1}", fps));
        if let Some(dt) = app.state.frame_times.last() {
            lines.push(format!("  Frame: {:.1}ms", dt * 1000.0));
        }
    }

    // Draw background and text
    if !lines.is_empty() {
        let font = egui::FontId::monospace(11.0);
        let text_color = egui::Color32::from_gray(220);
        let bg_color = egui::Color32::from_rgba_unmultiplied(20, 20, 25, 200);

        // Calculate text bounds
        let line_height = 14.0;
        let max_width = lines
            .iter()
            .map(|l| l.len() as f32 * 7.0) // Approximate character width
            .fold(100.0f32, |a, b| a.max(b));
        let total_height = lines.len() as f32 * line_height + 10.0;

        // Compute actual rect position based on anchor
        let (bg_x, bg_y) = match overlay.position {
            0 => (pos.x, pos.y),
            1 => (pos.x - max_width - 10.0, pos.y),
            2 => (pos.x, pos.y - total_height),
            _ => (pos.x - max_width - 10.0, pos.y - total_height),
        };

        let bg_rect = egui::Rect::from_min_size(
            egui::pos2(bg_x, bg_y),
            egui::vec2(max_width + 10.0, total_height),
        );

        painter.rect_filled(bg_rect, 4.0, bg_color);
        painter.rect_stroke(
            bg_rect,
            4.0,
            egui::Stroke::new(1.0_f32, egui::Color32::from_gray(60)),
            egui::StrokeKind::Outside,
        );

        // Draw each line
        for (i, line) in lines.iter().enumerate() {
            painter.text(
                egui::pos2(bg_x + 5.0, bg_y + 5.0 + i as f32 * line_height),
                egui::Align2::LEFT_TOP,
                line,
                font.clone(),
                text_color,
            );
        }
    }
}

/// Format a number with thousands separators
fn format_number(n: usize) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

/// Draw a color legend/scale bar for results visualization
fn draw_color_legend(painter: &egui::Painter, rect: egui::Rect, app: &FeaApp) {
    let results = match &app.state.results {
        Some(r) => r,
        None => return,
    };

    // Legend position (right side, vertically centered)
    let legend_x = rect.right() - 60.0;
    let legend_top = rect.top() + 60.0;
    let legend_height = 200.0;
    let legend_width = 20.0;

    // Get the label and range based on color mode
    let (label, min_val, max_val) = match app.state.ui_state.color_mode {
        ColorMode::Solid => return,
        ColorMode::Displacement => (
            "Displacement (m)",
            results.stats.min_displacement,
            results.stats.max_displacement,
        ),
        ColorMode::VonMises => ("Von Mises (Pa)", 0.0, results.stats.max_von_mises),
        ColorMode::Stress => {
            let comp_names = ["σ_XX", "σ_YY", "σ_ZZ", "σ_XY", "σ_YZ", "σ_XZ"];
            let idx = app.state.ui_state.stress_component.min(5);
            let max_val = results
                .stresses
                .values()
                .filter_map(|s| s.get(idx))
                .map(|v| v.abs())
                .fold(1e-10f64, |a, b| a.max(b));
            (comp_names[idx], -max_val, max_val)
        }
        ColorMode::Strain => {
            let comp_names = ["ε_XX", "ε_YY", "ε_ZZ", "ε_XY", "ε_YZ", "ε_XZ"];
            let idx = app.state.ui_state.strain_component.min(5);
            let max_val = results
                .strains
                .values()
                .filter_map(|s| s.get(idx))
                .map(|v| v.abs())
                .fold(1e-10f64, |a, b| a.max(b));
            (comp_names[idx], -max_val, max_val)
        }
    };

    // Draw background
    let bg_rect = egui::Rect::from_min_size(
        egui::pos2(legend_x - 5.0, legend_top - 20.0),
        egui::vec2(legend_width + 55.0, legend_height + 40.0),
    );
    painter.rect_filled(bg_rect, 4.0, egui::Color32::from_black_alpha(180));

    // Draw label
    painter.text(
        egui::pos2(legend_x + legend_width / 2.0, legend_top - 8.0),
        egui::Align2::CENTER_BOTTOM,
        label,
        egui::FontId::proportional(11.0),
        egui::Color32::WHITE,
    );

    // Draw color bar gradient
    let num_bands = 50;
    let band_height = legend_height / num_bands as f32;

    for i in 0..num_bands {
        let t = 1.0 - (i as f32 / num_bands as f32);
        let color = value_to_color(t);
        let y = legend_top + i as f32 * band_height;

        painter.rect_filled(
            egui::Rect::from_min_size(
                egui::pos2(legend_x, y),
                egui::vec2(legend_width, band_height + 1.0),
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
        egui::Stroke::new(1.0_f32, egui::Color32::WHITE),
        egui::StrokeKind::Outside,
    );

    // Draw tick marks and labels
    let tick_values = [max_val, (max_val + min_val) / 2.0, min_val];
    let tick_positions = [0.0, 0.5, 1.0];

    for (val, pos) in tick_values.iter().zip(tick_positions.iter()) {
        let y = legend_top + pos * legend_height;

        // Tick mark
        painter.line_segment(
            [
                egui::pos2(legend_x + legend_width, y),
                egui::pos2(legend_x + legend_width + 4.0, y),
            ],
            egui::Stroke::new(1.0_f32, egui::Color32::WHITE),
        );

        // Label
        let label = format_scientific(*val);
        painter.text(
            egui::pos2(legend_x + legend_width + 6.0, y),
            egui::Align2::LEFT_CENTER,
            label,
            egui::FontId::proportional(10.0),
            egui::Color32::WHITE,
        );
    }
}

/// Format a number in compact scientific notation
fn format_scientific(val: f64) -> String {
    if val.abs() < 1e-10 {
        "0".to_string()
    } else if val.abs() >= 1e4 || val.abs() < 1e-2 {
        format!("{:.1e}", val)
    } else {
        format!("{:.2}", val)
    }
}

/// Project a 3D point to 2D screen coordinates
fn project_point(
    point: [f32; 3],
    rect: egui::Rect,
    camera: &crate::state::CameraState,
) -> Option<egui::Pos2> {
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
    // When looking straight down or up, world Y is parallel to view direction
    let world_up = if forward[1].abs() > 0.99 {
        // Looking nearly straight up or down - use Z as up reference
        [0.0, 0.0, 1.0]
    } else {
        // Normal case - use Y as up reference
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
        return None; // Degenerate case
    }

    let right = [
        right[0] / right_len,
        right[1] / right_len,
        right[2] / right_len,
    ];

    // True up = right × forward
    let up = [
        right[1] * forward[2] - right[2] * forward[1],
        right[2] * forward[0] - right[0] * forward[2],
        right[0] * forward[1] - right[1] * forward[0],
    ];

    // Vector from eye to point
    let rel = [point[0] - eye[0], point[1] - eye[1], point[2] - eye[2]];

    // Camera space coordinates
    let cam_x = rel[0] * right[0] + rel[1] * right[1] + rel[2] * right[2];
    let cam_y = rel[0] * up[0] + rel[1] * up[1] + rel[2] * up[2];
    let cam_z = rel[0] * forward[0] + rel[1] * forward[1] + rel[2] * forward[2];

    // Don't render points behind camera or too close to near plane (perspective only)
    // Points very close to the camera cause extreme projection artifacts
    // Use a proportional near plane that scales with camera distance for close-up viewing
    let near_plane = (camera.distance * 0.01).max(0.001);
    if !camera.orthographic && cam_z < near_plane {
        return None;
    }

    let fov_factor = (camera.fov / 2.0).tan();
    let aspect = rect.width() / rect.height();

    let (ndc_x, ndc_y) = if camera.orthographic {
        // Orthographic projection: scale by camera distance to maintain size
        let ortho_scale = 1.0 / (camera.distance * fov_factor);
        (cam_x * ortho_scale / aspect, cam_y * ortho_scale)
    } else {
        // Perspective projection
        (
            cam_x / (cam_z * fov_factor * aspect),
            cam_y / (cam_z * fov_factor),
        )
    };

    // Convert to screen coordinates
    let screen_x = rect.center().x + ndc_x * rect.width() / 2.0;
    let screen_y = rect.center().y - ndc_y * rect.height() / 2.0;

    // Reject points that project too far outside the viewport (prevents spike artifacts)
    // This happens when points are very close to the near plane
    let max_extent = rect.width().max(rect.height()) * 2.0; // Allow 2x viewport size
    if screen_x.abs() > max_extent + rect.center().x.abs()
        || screen_y.abs() > max_extent + rect.center().y.abs()
        || !screen_x.is_finite()
        || !screen_y.is_finite()
    {
        return None;
    }

    Some(egui::pos2(screen_x, screen_y))
}

/// Get camera-space Z coordinate for depth sorting
fn camera_space_z(point: [f32; 3], camera: &crate::state::CameraState) -> f32 {
    let eye = camera.eye_position();
    let dx = point[0] - eye[0];
    let dy = point[1] - eye[1];
    let dz = point[2] - eye[2];
    dx * dx + dy * dy + dz * dz // Distance squared
}

/// Calculate the area of a 2D triangle using the cross product formula
fn triangle_area_2d(p0: egui::Pos2, p1: egui::Pos2, p2: egui::Pos2) -> f32 {
    let v1x = p1.x - p0.x;
    let v1y = p1.y - p0.y;
    let v2x = p2.x - p0.x;
    let v2y = p2.y - p0.y;
    (v1x * v2y - v1y * v2x).abs() * 0.5
}

/// Calculate the area of a 2D quad (sum of two triangles)
fn quad_area_2d(proj: &[egui::Pos2; 4]) -> f32 {
    // Split into two triangles: [0,1,2] and [0,2,3]
    triangle_area_2d(proj[0], proj[1], proj[2]) + triangle_area_2d(proj[0], proj[2], proj[3])
}

/// Calculate the length of a 2D line segment
fn line_length(p0: egui::Pos2, p1: egui::Pos2) -> f32 {
    let dx = p1.x - p0.x;
    let dy = p1.y - p0.y;
    (dx * dx + dy * dy).sqrt()
}

/// Convert normalized value (0-1) to color using jet colormap
fn value_to_color(t: f32) -> egui::Color32 {
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);

    egui::Color32::from_rgb((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8)
}

#[allow(dead_code)]
fn average_colors(colors: &[egui::Color32]) -> egui::Color32 {
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
        (a / n) as u8,
    )
}

/// Average exactly 3 colors (optimized for triangle rendering)
fn average_colors_3(c0: egui::Color32, c1: egui::Color32, c2: egui::Color32) -> egui::Color32 {
    let r = ((c0.r() as u32 + c1.r() as u32 + c2.r() as u32) / 3) as u8;
    let g = ((c0.g() as u32 + c1.g() as u32 + c2.g() as u32) / 3) as u8;
    let b = ((c0.b() as u32 + c1.b() as u32 + c2.b() as u32) / 3) as u8;
    egui::Color32::from_rgb(r, g, b)
}

/// Darken a color by a factor (0.0 = black, 1.0 = unchanged)
#[allow(dead_code)]
fn darken_color(color: egui::Color32, factor: f32) -> egui::Color32 {
    egui::Color32::from_rgba_unmultiplied(
        (color.r() as f32 * factor) as u8,
        (color.g() as f32 * factor) as u8,
        (color.b() as f32 * factor) as u8,
        color.a(),
    )
}

/// Draw an arrow with clamped maximum length to prevent visual artifacts
fn draw_arrow_clamped(
    painter: &egui::Painter,
    from: egui::Pos2,
    to: egui::Pos2,
    color: egui::Color32,
    stroke_width: f32,
    max_length: f32,
) {
    let delta = to - from;
    let length = delta.length();

    // Skip if arrow is too small or has invalid length
    if length < 1.0 || !length.is_finite() {
        return;
    }

    // Clamp length if needed
    let (final_to, final_delta) = if length > max_length {
        let scale = max_length / length;
        let clamped_delta = delta * scale;
        (from + clamped_delta, clamped_delta)
    } else {
        (to, delta)
    };

    // Draw line
    painter.line_segment([from, final_to], egui::Stroke::new(stroke_width, color));

    // Draw arrowhead
    let arrow_size = (length.min(max_length) * 0.2).clamp(4.0, 12.0);
    let dir = final_delta.normalized();
    let perp = egui::vec2(-dir.y, dir.x);

    let tip = final_to;
    let left = tip - dir * arrow_size + perp * arrow_size * 0.5;
    let right = tip - dir * arrow_size - perp * arrow_size * 0.5;

    painter.add(egui::Shape::convex_polygon(
        vec![tip, left, right],
        color,
        egui::Stroke::NONE,
    ));
}

/// Handle keyboard shortcuts for the application
fn handle_keyboard_shortcuts(ctx: &egui::Context, app: &mut FeaApp) {
    // Check if any text widget has focus (skip shortcuts when typing)
    let has_text_focus = ctx.memory(|m| m.focused().is_some());
    if has_text_focus {
        return;
    }

    ctx.input(|i| {
        // View toggles
        if i.key_pressed(egui::Key::F) {
            app.state.ui_state.show_faces = !app.state.ui_state.show_faces;
            app.renderer = None;
        }
        if i.key_pressed(egui::Key::W) && !i.modifiers.ctrl && !i.modifiers.command {
            app.state.ui_state.show_wireframe = !app.state.ui_state.show_wireframe;
            app.renderer = None;
        }
        if i.key_pressed(egui::Key::N) {
            app.state.ui_state.show_nodes = !app.state.ui_state.show_nodes;
            app.renderer = None;
        }
        if i.key_pressed(egui::Key::B) {
            app.state.ui_state.show_boundary_conditions =
                !app.state.ui_state.show_boundary_conditions;
        }
        // Stats overlay toggle
        if i.key_pressed(egui::Key::I) && !i.modifiers.ctrl && !i.modifiers.command {
            app.state.ui_state.stats_overlay.visible = !app.state.ui_state.stats_overlay.visible;
        }
        // Grid toggle
        if i.key_pressed(egui::Key::G) && !i.modifiers.ctrl && !i.modifiers.command {
            app.state.ui_state.display_settings.show_grid =
                !app.state.ui_state.display_settings.show_grid;
        }
        if i.key_pressed(egui::Key::C) && !i.modifiers.ctrl && !i.modifiers.command {
            app.state.ui_state.clipping_plane.enabled = !app.state.ui_state.clipping_plane.enabled;
            app.render_cache.invalidate();
        }

        // Camera views (numpad/number keys)
        if i.key_pressed(egui::Key::Num1) {
            app.state.ui_state.camera.set_front_view();
        }
        if i.key_pressed(egui::Key::Num2) {
            app.state.ui_state.camera.set_side_view();
        }
        if i.key_pressed(egui::Key::Num3) {
            app.state.ui_state.camera.set_top_view();
        }
        if i.key_pressed(egui::Key::Num0) {
            app.state.ui_state.camera.set_iso_view();
        }

        // Projection mode toggle
        if i.key_pressed(egui::Key::P) && !i.modifiers.ctrl && !i.modifiers.command {
            app.state.ui_state.camera.orthographic = !app.state.ui_state.camera.orthographic;
        }

        // Fit to view (Home key)
        if i.key_pressed(egui::Key::Home) {
            if let Some(mesh) = app.state.current_mesh() {
                let bounds = mesh.bounds;
                app.state.ui_state.camera.fit_to_bounds(&bounds);
            }
        }

        // Help toggle
        if i.key_pressed(egui::Key::Slash) && i.modifiers.shift {
            // Shift+/ = ?
            app.state.ui_state.show_shortcuts_help = !app.state.ui_state.show_shortcuts_help;
        }

        // Animation controls (when results are available)
        if app.state.results.is_some() {
            if i.key_pressed(egui::Key::Space) {
                app.state.ui_state.animation.playing = !app.state.ui_state.animation.playing;
            }
            if i.key_pressed(egui::Key::ArrowLeft) {
                if app.state.ui_state.current_time_step > 0 {
                    app.state.ui_state.current_time_step -= 1;
                    app.renderer = None;
                }
            }
            if i.key_pressed(egui::Key::ArrowRight) {
                let max_steps = app
                    .state
                    .results
                    .as_ref()
                    .map(|r| r.time_steps.len().saturating_sub(1))
                    .unwrap_or(0);
                if app.state.ui_state.current_time_step < max_steps {
                    app.state.ui_state.current_time_step += 1;
                    app.renderer = None;
                }
            }
        }

        // Panel switching (Ctrl/Cmd + number)
        if i.modifiers.ctrl || i.modifiers.command {
            if i.key_pressed(egui::Key::Num1) {
                app.state.ui_state.active_panel = crate::state::ActivePanel::Mesh;
            }
            if i.key_pressed(egui::Key::Num2) {
                app.state.ui_state.active_panel = crate::state::ActivePanel::Setup;
            }
            if i.key_pressed(egui::Key::Num3) {
                app.state.ui_state.active_panel = crate::state::ActivePanel::Run;
            }
            if i.key_pressed(egui::Key::Num4) {
                app.state.ui_state.active_panel = crate::state::ActivePanel::Results;
            }

            // Undo (Ctrl+Z)
            if i.key_pressed(egui::Key::Z) && !i.modifiers.shift {
                crate::ui::menu_bar::perform_undo_public(app);
            }

            // Redo (Ctrl+Shift+Z)
            if i.key_pressed(egui::Key::Z) && i.modifiers.shift {
                crate::ui::menu_bar::perform_redo_public(app);
            }

            // Save settings (Ctrl+S)
            if i.key_pressed(egui::Key::S) {
                app.state.save_settings();
                app.state.status_message = "Settings saved".to_string();
            }
        }

        // Clipping plane axis (X, Y, Z when clipping is enabled)
        if app.state.ui_state.clipping_plane.enabled {
            if i.key_pressed(egui::Key::X) {
                app.state.ui_state.clipping_plane.axis = ClipAxis::X;
                app.state.ui_state.clipping_plane.normal = ClipAxis::X.to_normal();
                app.render_cache.invalidate();
            }
            if i.key_pressed(egui::Key::Y) {
                app.state.ui_state.clipping_plane.axis = ClipAxis::Y;
                app.state.ui_state.clipping_plane.normal = ClipAxis::Y.to_normal();
                app.render_cache.invalidate();
            }
            if i.key_pressed(egui::Key::Z) {
                app.state.ui_state.clipping_plane.axis = ClipAxis::Z;
                app.state.ui_state.clipping_plane.normal = ClipAxis::Z.to_normal();
                app.render_cache.invalidate();
            }
        }
    });

    // Show shortcuts help window if enabled
    if app.state.ui_state.show_shortcuts_help {
        show_shortcuts_window(ctx, app);
    }
}

/// Show clipping plane controls
fn show_clipping_controls(ui: &mut egui::Ui, app: &mut FeaApp) {
    ui.horizontal(|ui| {
        ui.label("Clip Axis:");

        let mut changed = false;

        if ui
            .selectable_label(app.state.ui_state.clipping_plane.axis == ClipAxis::X, "X")
            .on_hover_text("Clip along X axis")
            .clicked()
        {
            app.state.ui_state.clipping_plane.axis = ClipAxis::X;
            app.state.ui_state.clipping_plane.normal = ClipAxis::X.to_normal();
            changed = true;
        }
        if ui
            .selectable_label(app.state.ui_state.clipping_plane.axis == ClipAxis::Y, "Y")
            .on_hover_text("Clip along Y axis")
            .clicked()
        {
            app.state.ui_state.clipping_plane.axis = ClipAxis::Y;
            app.state.ui_state.clipping_plane.normal = ClipAxis::Y.to_normal();
            changed = true;
        }
        if ui
            .selectable_label(app.state.ui_state.clipping_plane.axis == ClipAxis::Z, "Z")
            .on_hover_text("Clip along Z axis")
            .clicked()
        {
            app.state.ui_state.clipping_plane.axis = ClipAxis::Z;
            app.state.ui_state.clipping_plane.normal = ClipAxis::Z.to_normal();
            changed = true;
        }

        ui.separator();

        // Position slider
        ui.label("Position:");
        let response = ui.add(
            egui::Slider::new(&mut app.state.ui_state.clipping_plane.position, -1.0..=1.0)
                .show_value(true)
                .suffix(""),
        );
        if response.changed() {
            changed = true;
        }

        ui.separator();

        // Flip direction
        if ui
            .checkbox(&mut app.state.ui_state.clipping_plane.flip, "Flip")
            .on_hover_text("Flip clipping direction")
            .changed()
        {
            changed = true;
        }

        // Section surface toggle
        if ui
            .checkbox(
                &mut app.state.ui_state.clipping_plane.show_section_surface,
                "Section Fill",
            )
            .on_hover_text("Show cross-section surface with field values")
            .changed()
        {
            changed = true;
        }

        // 2D cross-section view toggle
        ui.checkbox(
            &mut app.state.ui_state.clipping_plane.show_2d_view,
            "2D View",
        )
        .on_hover_text("Show 2D orthographic view of the cross-section");

        if changed {
            app.render_cache.invalidate();
        }
    });
}

/// Show keyboard shortcuts help window
fn show_shortcuts_window(ctx: &egui::Context, app: &mut FeaApp) {
    egui::Window::new("⌨ Keyboard Shortcuts")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.heading("View Controls");
            egui::Grid::new("shortcuts_view")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("F");
                    ui.label("Toggle faces");
                    ui.end_row();
                    ui.strong("W");
                    ui.label("Toggle wireframe");
                    ui.end_row();
                    ui.strong("N");
                    ui.label("Toggle nodes");
                    ui.end_row();
                    ui.strong("B");
                    ui.label("Toggle boundary conditions");
                    ui.end_row();
                    ui.strong("C");
                    ui.label("Toggle clipping plane");
                    ui.end_row();
                    ui.strong("G");
                    ui.label("Toggle grid");
                    ui.end_row();
                    ui.strong("I");
                    ui.label("Toggle stats overlay");
                    ui.end_row();
                    ui.strong("P");
                    ui.label("Toggle perspective/ortho");
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.heading("Camera");
            egui::Grid::new("shortcuts_camera")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("1");
                    ui.label("Front view");
                    ui.end_row();
                    ui.strong("2");
                    ui.label("Side view");
                    ui.end_row();
                    ui.strong("3");
                    ui.label("Top view");
                    ui.end_row();
                    ui.strong("0");
                    ui.label("Isometric view");
                    ui.end_row();
                    ui.strong("Home");
                    ui.label("Fit to mesh");
                    ui.end_row();
                    ui.strong("LMB drag");
                    ui.label("Rotate");
                    ui.end_row();
                    ui.strong("RMB drag");
                    ui.label("Pan");
                    ui.end_row();
                    ui.strong("Scroll");
                    ui.label("Zoom");
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.heading("Edit");
            egui::Grid::new("shortcuts_edit")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("Ctrl+Z");
                    ui.label("Undo");
                    ui.end_row();
                    ui.strong("Ctrl+Shift+Z");
                    ui.label("Redo");
                    ui.end_row();
                    ui.strong("Ctrl+S");
                    ui.label("Save settings");
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.heading("Clipping (when active)");
            egui::Grid::new("shortcuts_clip")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("X");
                    ui.label("Clip along X axis");
                    ui.end_row();
                    ui.strong("Y");
                    ui.label("Clip along Y axis");
                    ui.end_row();
                    ui.strong("Z");
                    ui.label("Clip along Z axis");
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.heading("Animation (with results)");
            egui::Grid::new("shortcuts_anim")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("Space");
                    ui.label("Play/pause animation");
                    ui.end_row();
                    ui.strong("</>");
                    ui.label("Previous/next frame");
                    ui.end_row();
                });

            ui.add_space(8.0);
            ui.heading("Panels");
            egui::Grid::new("shortcuts_panels")
                .num_columns(2)
                .spacing([20.0, 4.0])
                .show(ui, |ui| {
                    ui.strong("Ctrl+1");
                    ui.label("Mesh panel");
                    ui.end_row();
                    ui.strong("Ctrl+2");
                    ui.label("Setup panel");
                    ui.end_row();
                    ui.strong("Ctrl+3");
                    ui.label("Run panel");
                    ui.end_row();
                    ui.strong("Ctrl+4");
                    ui.label("Results panel");
                    ui.end_row();
                });

            ui.add_space(12.0);
            ui.horizontal(|ui| {
                if ui.button("Close").clicked() {
                    app.state.ui_state.show_shortcuts_help = false;
                }
            });
        });
}
