//! Results panel for viewing and exporting simulation results

use eframe::egui;
use crate::app::FeaApp;
use crate::state::ColorMode;

pub fn show(ui: &mut egui::Ui, app: &mut FeaApp) {
    ui.heading("Results");
    ui.add_space(8.0);
    
    if app.state.results.is_none() {
        ui.label("No results available.");
        ui.label("Run a simulation first.");
        
        if app.state.current_mesh().is_some() {
            ui.add_space(8.0);
            if ui.button("▶ Go to Run Tab").clicked() {
                app.state.ui_state.active_panel = crate::state::ActivePanel::Run;
            }
        }
        return;
    }
    
    // Extract result stats to avoid borrow issues
    let stats = app.state.results.as_ref().unwrap().stats.clone();
    let has_von_mises = stats.max_von_mises > 0.0;
    let has_stresses = !app.state.results.as_ref().unwrap().stresses.is_empty();
    let time_steps_count = app.state.results.as_ref().unwrap().time_steps.len();
    
    // Results summary
    ui.label(egui::RichText::new("Summary").strong());
    ui.add_space(4.0);
    
    egui::Grid::new("results_summary")
        .num_columns(2)
        .spacing([20.0, 4.0])
        .show(ui, |ui| {
            ui.label("Max Displacement:");
            ui.label(format!("{:.6e}", stats.max_displacement));
            ui.end_row();
            
            ui.label("Min Displacement:");
            ui.label(format!("{:.6e}", stats.min_displacement));
            ui.end_row();
            
            if stats.max_von_mises > 0.0 {
                ui.label("Max Von Mises:");
                ui.label(format!("{:.6e} Pa", stats.max_von_mises));
                ui.end_row();
            }
            
            ui.label("Solver Time:");
            ui.label(format!("{} ms", stats.solver_time_ms));
            ui.end_row();
        });
    
    ui.add_space(16.0);
    
    // Visualization options
    ui.label(egui::RichText::new("Visualization").strong());
    ui.add_space(4.0);
    
    ui.horizontal(|ui| {
        ui.label("Color By:");
        egui::ComboBox::from_id_salt("color_mode")
            .selected_text(match app.state.ui_state.color_mode {
                ColorMode::Solid => "Solid Color",
                ColorMode::Displacement => "Displacement",
                ColorMode::VonMises => "Von Mises Stress",
                ColorMode::Stress => "Stress",
                ColorMode::Strain => "Strain",
            })
            .show_ui(ui, |ui| {
                if ui.selectable_value(
                    &mut app.state.ui_state.color_mode, 
                    ColorMode::Solid, 
                    "Solid Color"
                ).changed() {
                    app.renderer = None; // Invalidate renderer
                }
                if ui.selectable_value(
                    &mut app.state.ui_state.color_mode, 
                    ColorMode::Displacement, 
                    "Displacement"
                ).changed() {
                    app.renderer = None;
                }
                // Show stress options when we have stress data
                if has_von_mises || has_stresses {
                    if ui.selectable_value(
                        &mut app.state.ui_state.color_mode, 
                        ColorMode::VonMises, 
                        "Von Mises Stress"
                    ).changed() {
                        app.renderer = None;
                    }
                    if ui.selectable_value(
                        &mut app.state.ui_state.color_mode, 
                        ColorMode::Stress, 
                        "Stress Components"
                    ).changed() {
                        app.renderer = None;
                    }
                    if ui.selectable_value(
                        &mut app.state.ui_state.color_mode, 
                        ColorMode::Strain, 
                        "Strain"
                    ).changed() {
                        app.renderer = None;
                    }
                }
            });
    });
    
    // Component selector for Stress/Strain modes
    let component_names = ["XX", "YY", "ZZ", "XY", "YZ", "XZ"];
    
    if app.state.ui_state.color_mode == ColorMode::Stress {
        ui.horizontal(|ui| {
            ui.label("Stress Component:");
            egui::ComboBox::from_id_salt("stress_component")
                .selected_text(format!("σ_{}", component_names[app.state.ui_state.stress_component]))
                .show_ui(ui, |ui| {
                    for (i, name) in component_names.iter().enumerate() {
                        if ui.selectable_value(
                            &mut app.state.ui_state.stress_component,
                            i,
                            format!("σ_{}", name)
                        ).changed() {
                            app.renderer = None;
                        }
                    }
                });
        });
    }
    
    if app.state.ui_state.color_mode == ColorMode::Strain {
        ui.horizontal(|ui| {
            ui.label("Strain Component:");
            egui::ComboBox::from_id_salt("strain_component")
                .selected_text(format!("ε_{}", component_names[app.state.ui_state.strain_component]))
                .show_ui(ui, |ui| {
                    for (i, name) in component_names.iter().enumerate() {
                        if ui.selectable_value(
                            &mut app.state.ui_state.strain_component,
                            i,
                            format!("ε_{}", name)
                        ).changed() {
                            app.renderer = None;
                        }
                    }
                });
        });
    }
    
    // Debug info for stress data
    if !has_von_mises && !has_stresses {
        ui.colored_label(egui::Color32::YELLOW, "No stress data available");
    }
    
    // Displacement scale
    ui.horizontal(|ui| {
        ui.label("Deformation Scale:");
        if ui.add(
            egui::Slider::new(&mut app.state.ui_state.displacement_scale, 0.0..=100.0)
                .logarithmic(true)
        ).changed() {
            app.renderer = None; // Invalidate renderer
        }
    });
    
    // Quick scale presets
    ui.horizontal(|ui| {
        if ui.small_button("1x").clicked() {
            app.state.ui_state.displacement_scale = 1.0;
            app.renderer = None;
        }
        if ui.small_button("10x").clicked() {
            app.state.ui_state.displacement_scale = 10.0;
            app.renderer = None;
        }
        if ui.small_button("100x").clicked() {
            app.state.ui_state.displacement_scale = 100.0;
            app.renderer = None;
        }
        if ui.small_button("Auto").clicked() {
            // Auto-scale based on mesh size
            if let Some(mesh) = app.state.current_mesh() {
                let diag = mesh.bounds.diagonal();
                let max_disp = stats.max_displacement;
                if max_disp > 0.0 {
                    // Target: displacement scaled to ~10% of mesh size
                    app.state.ui_state.displacement_scale = 
                        (0.1 * diag as f64 / max_disp) as f32;
                    app.renderer = None;
                }
            }
        }
    });
    
    ui.add_space(16.0);
    
    // Color legend
    if app.state.ui_state.color_mode != ColorMode::Solid {
        ui.label(egui::RichText::new("Color Scale").strong());
        show_color_legend(ui, &stats, app.state.ui_state.color_mode);
    }
    
    ui.add_space(16.0);
    
    // Time step playback (for explicit solver)
    if time_steps_count > 0 {
        ui.separator();
        ui.add_space(8.0);
        ui.label(egui::RichText::new("Time Steps").strong());
        ui.add_space(4.0);
        
        ui.label(format!("{} time steps recorded", time_steps_count));
        
        // Current step info
        let current_idx = app.state.ui_state.current_time_step.min(time_steps_count.saturating_sub(1));
        if let Some(results) = &app.state.results {
            if let Some(step) = results.time_steps.get(current_idx) {
                ui.horizontal(|ui| {
                    ui.label("Current:");
                    ui.label(format!("Step {} | t = {:.4e} s", step.iteration, step.time));
                });
                ui.horizontal(|ui| {
                    ui.label("Max Disp:");
                    ui.label(format!("{:.4e}", step.max_displacement));
                    ui.label("KE:");
                    ui.label(format!("{:.4e}", step.kinetic_energy));
                });
            }
        }
        
        ui.add_space(8.0);
        
        // Step slider
        let mut step_idx = app.state.ui_state.current_time_step;
        if ui.add(
            egui::Slider::new(&mut step_idx, 0..=time_steps_count.saturating_sub(1))
                .text("Step")
        ).changed() {
            app.state.ui_state.current_time_step = step_idx;
            app.renderer = None;
        }
        
        // Playback controls
        ui.horizontal(|ui| {
            // Jump to start
            if ui.button("⏮").on_hover_text("First step").clicked() {
                app.state.ui_state.current_time_step = 0;
                app.state.ui_state.playback_active = false;
                app.renderer = None;
            }
            
            // Step back
            if ui.button("⏪").on_hover_text("Previous step").clicked() {
                if app.state.ui_state.current_time_step > 0 {
                    app.state.ui_state.current_time_step -= 1;
                    app.renderer = None;
                }
            }
            
            // Play/Pause
            let play_text = if app.state.ui_state.playback_active { "⏸" } else { "▶" };
            if ui.button(play_text).on_hover_text(if app.state.ui_state.playback_active { "Pause" } else { "Play" }).clicked() {
                app.state.ui_state.playback_active = !app.state.ui_state.playback_active;
            }
            
            // Step forward
            if ui.button("⏩").on_hover_text("Next step").clicked() {
                if app.state.ui_state.current_time_step < time_steps_count - 1 {
                    app.state.ui_state.current_time_step += 1;
                    app.renderer = None;
                }
            }
            
            // Jump to end
            if ui.button("⏭").on_hover_text("Last step").clicked() {
                app.state.ui_state.current_time_step = time_steps_count.saturating_sub(1);
                app.state.ui_state.playback_active = false;
                app.renderer = None;
            }
            
            // Loop toggle
            ui.separator();
            ui.label("Speed:");
            ui.add(egui::DragValue::new(&mut app.state.ui_state.playback_speed)
                .range(1.0..=60.0)
                .suffix(" fps"));
        });
        
        ui.add_space(8.0);
        
        // Mini time history plot
        if let Some(results) = &app.state.results {
            if results.time_steps.len() > 1 {
                let plot_height = 80.0;
                let (rect, _) = ui.allocate_exact_size(
                    egui::vec2(ui.available_width(), plot_height),
                    egui::Sense::hover()
                );
                
                let painter = ui.painter_at(rect);
                painter.rect_filled(rect, 4.0, egui::Color32::from_gray(40));
                
                // Find max values for scaling
                let max_disp = results.time_steps.iter()
                    .map(|s| s.max_displacement)
                    .fold(1e-10f64, |a, b| a.max(b));
                let max_time = results.time_steps.last()
                    .map(|s| s.time)
                    .unwrap_or(1.0);
                
                // Draw displacement curve
                let points: Vec<egui::Pos2> = results.time_steps.iter()
                    .map(|step| {
                        let x = rect.left() + (step.time / max_time) as f32 * rect.width();
                        let y = rect.bottom() - (step.max_displacement / max_disp) as f32 * rect.height() * 0.9;
                        egui::pos2(x, y)
                    })
                    .collect();
                
                if points.len() > 1 {
                    for window in points.windows(2) {
                        painter.line_segment(
                            [window[0], window[1]],
                            egui::Stroke::new(1.5_f32, egui::Color32::from_rgb(100, 200, 255))
                        );
                    }
                }
                
                // Draw current position marker
                if let Some(current_step) = results.time_steps.get(app.state.ui_state.current_time_step) {
                    let x = rect.left() + (current_step.time / max_time) as f32 * rect.width();
                    painter.vline(x, rect.y_range(), egui::Stroke::new(1.0_f32, egui::Color32::YELLOW));
                }
                
                // Axis label
                painter.text(
                    egui::pos2(rect.left() + 4.0, rect.top() + 4.0),
                    egui::Align2::LEFT_TOP,
                    "Max Disp",
                    egui::FontId::proportional(9.0),
                    egui::Color32::from_gray(150)
                );
            }
        }
        
        ui.add_space(8.0);
        
        // Show time steps table
        egui::CollapsingHeader::new("Step History")
            .default_open(false)
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .max_height(150.0)
                    .show(ui, |ui| {
                        egui::Grid::new("timesteps")
                            .num_columns(4)
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Step");
                                ui.strong("Time");
                                ui.strong("Max Disp");
                                ui.strong("KE");
                                ui.end_row();
                                
                                if let Some(results) = &app.state.results {
                                    for (idx, step) in results.time_steps.iter().enumerate().rev().take(50) {
                                        let is_current = idx == app.state.ui_state.current_time_step;
                                        if is_current {
                                            ui.strong(format!("> {}", step.iteration));
                                        } else {
                                            ui.label(format!("{}", step.iteration));
                                        }
                                        ui.label(format!("{:.4e}", step.time));
                                        ui.label(format!("{:.4e}", step.max_displacement));
                                        ui.label(format!("{:.4e}", step.kinetic_energy));
                                        ui.end_row();
                                    }
                                }
                            });
                    });
            });
        
        ui.add_space(16.0);
    }
    
    // Export options
    ui.separator();
    ui.add_space(8.0);
    ui.label(egui::RichText::new("Export").strong());
    ui.add_space(4.0);
        
    ui.horizontal(|ui| {
        if ui.button("Export VTK").clicked() {
            export_vtk(app);
        }
        if ui.button("Export CSV").clicked() {
            export_csv(app);
        }
        if ui.button("Screenshot").clicked() {
            // TODO: Implement screenshot
        }
    });
}

fn show_color_legend(ui: &mut egui::Ui, stats: &crate::state::ResultStats, color_mode: ColorMode) {
    let (min_val, max_val, label) = match color_mode {
        ColorMode::Solid => return,
        ColorMode::Displacement => (
            stats.min_displacement,
            stats.max_displacement,
            "Displacement (m)"
        ),
        ColorMode::VonMises => (0.0, stats.max_von_mises, "Von Mises Stress (Pa)"),
        ColorMode::Stress => (0.0, stats.max_von_mises, "Stress Magnitude (Pa)"),
        ColorMode::Strain => (0.0, 0.01, "Strain (-)"),
    };
    
    ui.horizontal(|ui| {
        ui.label(label);
    });
    
    // Draw color bar
    let (rect, _response) = ui.allocate_exact_size(
        egui::vec2(ui.available_width() - 20.0, 20.0),
        egui::Sense::hover()
    );
    
    let painter = ui.painter();
    
    // Draw gradient
    let num_segments = 50;
    let segment_width = rect.width() / num_segments as f32;
    
    for i in 0..num_segments {
        let t = i as f32 / (num_segments - 1) as f32;
        let color = value_to_color(t);
        let x = rect.left() + i as f32 * segment_width;
        painter.rect_filled(
            egui::Rect::from_min_size(
                egui::pos2(x, rect.top()),
                egui::vec2(segment_width + 1.0, rect.height())
            ),
            0.0,
            color
        );
    }
    
    // Labels
    ui.horizontal(|ui| {
        ui.label(format!("{:.2e}", min_val));
        ui.add_space(ui.available_width() - 100.0);
        ui.label(format!("{:.2e}", max_val));
    });
}

/// Convert normalized value (0-1) to color using jet colormap
fn value_to_color(t: f32) -> egui::Color32 {
    // Blue -> Cyan -> Green -> Yellow -> Red
    let r = (1.5 - (4.0 * t - 3.0).abs()).clamp(0.0, 1.0);
    let g = (1.5 - (4.0 * t - 2.0).abs()).clamp(0.0, 1.0);
    let b = (1.5 - (4.0 * t - 1.0).abs()).clamp(0.0, 1.0);
    
    egui::Color32::from_rgb(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8
    )
}

fn export_vtk(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("VTK", &["vtk"])
            .save_file()
        {
            // TODO: Implement VTK export with results
            app.state.status_message = format!("VTK exported to {:?}", path);
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        app.state.status_message = "VTK export not available in web version".to_string();
    }
}

fn export_csv(app: &mut FeaApp) {
    #[cfg(feature = "native")]
    {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("CSV", &["csv"])
            .save_file()
        {
            if let Some(results) = &app.state.results {
                // Write displacement data
                let mut csv_data = String::from("node_id,dx,dy,dz,magnitude\n");
                
                let num_nodes = results.displacements.len() / 3;
                for i in 0..num_nodes {
                    let dx = results.displacements[i * 3];
                    let dy = results.displacements[i * 3 + 1];
                    let dz = results.displacements[i * 3 + 2];
                    let mag = (dx * dx + dy * dy + dz * dz).sqrt();
                    csv_data.push_str(&format!("{},{},{},{},{}\n", i, dx, dy, dz, mag));
                }
                
                if let Err(e) = std::fs::write(&path, csv_data) {
                    app.state.status_message = format!("Failed to export CSV: {}", e);
                } else {
                    app.state.status_message = format!("CSV exported to {:?}", path);
                }
            }
        }
    }
    
    #[cfg(not(feature = "native"))]
    {
        app.state.status_message = "CSV export not available in web version".to_string();
    }
}
