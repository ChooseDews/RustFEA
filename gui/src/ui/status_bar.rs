//! Bottom status bar

use eframe::egui;
use crate::app::FeaApp;

pub fn show(ctx: &egui::Context, app: &mut FeaApp) {
    egui::TopBottomPanel::bottom("status_bar")
        .resizable(false)
        .show(ctx, |ui| {
            ui.horizontal(|ui| {
                // Status indicator icon
                let status_icon = if app.state.is_running {
                    "◐"
                } else if app.state.results.is_some() {
                    "●"
                } else {
                    "○"
                };
                ui.label(status_icon);
                
                // Status message
                ui.label(&app.state.status_message);
                
                ui.separator();
                
                // Mesh info with tooltip
                if let Some(mesh) = app.state.current_mesh() {
                    let mesh_info = format!(
                        "Nodes: {} | Elements: {} | Groups: {}",
                        mesh.mesh.nodes.len(),
                        mesh.mesh.elements.len(),
                        mesh.mesh.node_groups.len() + mesh.mesh.element_groups.len()
                    );
                    ui.label(&mesh_info).on_hover_ui(|ui| {
                        ui.label("Mesh Statistics");
                        ui.separator();
                        egui::Grid::new("mesh_tooltip")
                            .num_columns(2)
                            .spacing([10.0, 2.0])
                            .show(ui, |ui| {
                                ui.label("Total Nodes:");
                                ui.label(format!("{}", mesh.mesh.nodes.len()));
                                ui.end_row();
                                ui.label("Total Elements:");
                                ui.label(format!("{}", mesh.mesh.elements.len()));
                                ui.end_row();
                                ui.label("Node Groups:");
                                ui.label(format!("{}", mesh.mesh.node_groups.len()));
                                ui.end_row();
                                ui.label("Element Groups:");
                                ui.label(format!("{}", mesh.mesh.element_groups.len()));
                                ui.end_row();
                                ui.label("Bodies:");
                                ui.label(format!("{}", mesh.mesh.bodies.len()));
                                ui.end_row();
                            });
                    });
                    
                    // Show clipping plane status
                    if app.state.ui_state.clipping_plane.enabled {
                        ui.separator();
                        ui.colored_label(egui::Color32::LIGHT_BLUE, "Clipped");
                    }
                }
                
                // Animation status
                if app.state.ui_state.playback_active {
                    ui.separator();
                        ui.colored_label(egui::Color32::LIGHT_GREEN, "Playing");
                }
                
                // Progress bar if running
                if app.state.is_running {
                    ui.separator();
                    ui.add(
                        egui::ProgressBar::new(app.state.progress)
                            .show_percentage()
                            .animate(true)
                    );
                }
                
                // Right-aligned items
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    // Camera info with tooltip
                    let cam = &app.state.ui_state.camera;
                    ui.label(format!(
                        "Cam: ({:.1}°, {:.1}°) d={:.1}",
                        cam.yaw.to_degrees(),
                        cam.pitch.to_degrees(),
                        cam.distance
                    )).on_hover_text("Camera: (yaw, pitch) distance\nLMB: Rotate | RMB: Pan | Scroll: Zoom");
                    
                    // Render stats (only show if mesh loaded)
                    if app.state.current_mesh().is_some() {
                        ui.separator();
                        let cache = &app.render_cache;
                        let efficiency = if cache.faces.len() > 0 {
                            (cache.last_rendered_faces as f32 / cache.faces.len() as f32 * 100.0) as u32
                        } else {
                            0
                        };
                        ui.label(format!(
                            "Faces: {}/{} ({}%) | Tris: {}",
                            cache.last_rendered_faces,
                            cache.faces.len(),
                            efficiency,
                            cache.last_rendered_triangles
                        )).on_hover_ui(|ui| {
                            ui.label("Render Statistics");
                            ui.separator();
                            ui.label(format!("Visible faces: {}", cache.last_rendered_faces));
                            ui.label(format!("Total faces: {}", cache.faces.len()));
                            ui.label(format!("Triangles drawn: {}", cache.last_rendered_triangles));
                            ui.label(format!("Culling efficiency: {}%", 100 - efficiency));
                            if cache.octree.is_some() {
                                ui.label("Octree enabled");
                            }
                        });
                    }
                });
            });
        });
}
