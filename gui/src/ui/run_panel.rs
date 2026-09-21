//! Run panel for executing simulations with detailed progress tracking

use eframe::egui;
use crate::app::FeaApp;
use crate::state::{SolverType, SolvePhaseCategory};

pub fn show(ui: &mut egui::Ui, app: &mut FeaApp) {
    ui.heading("Run Simulation");
    ui.add_space(8.0);
    
    // Pre-run checks
    let has_mesh = app.state.current_mesh().is_some();
    let has_bc = !app.state.simulation_config.boundary_conditions.is_empty();
    let is_running = app.state.is_running;
    
    // Status checks
    ui.label(egui::RichText::new("Pre-flight Check").strong());
    ui.add_space(4.0);
    
    show_check(ui, "Mesh loaded", has_mesh);
    show_check(ui, "Boundary conditions defined", has_bc);
    
    if !has_mesh {
        ui.label("Import a mesh in the Mesh tab");
    }
    if !has_bc {
        ui.label("Add boundary conditions in the Setup tab");
    }
    
    ui.add_space(16.0);
    
    // Simulation summary
    if let Some(mesh) = app.state.current_mesh() {
        ui.label(egui::RichText::new("Simulation Summary").strong());
        ui.add_space(4.0);
        
        egui::Grid::new("sim_summary")
            .num_columns(2)
            .spacing([20.0, 4.0])
            .show(ui, |ui| {
                ui.label("Mesh:");
                ui.label(&mesh.name);
                ui.end_row();
                
                ui.label("Nodes:");
                ui.label(format!("{}", mesh.mesh.nodes.len()));
                ui.end_row();
                
                ui.label("Elements:");
                ui.label(format!("{}", mesh.mesh.elements.len()));
                ui.end_row();
                
                ui.label("DOFs:");
                let dofs = app.state.simulation_config.dofs;
                let total_dofs = mesh.mesh.nodes.len() * dofs;
                ui.label(format!("{} per node ({} total)", dofs, total_dofs));
                ui.end_row();
                
                ui.label("Solver:");
                ui.label(match app.state.simulation_config.solver {
                    SolverType::Direct => "Direct (Linear)",
                    SolverType::Explicit => "Explicit (Time-stepping)",
                });
                ui.end_row();
                
                ui.label("BCs:");
                ui.label(format!("{}", app.state.simulation_config.boundary_conditions.len()));
                ui.end_row();
            });
    }
    
    ui.add_space(16.0);
    
    // Run controls
    ui.add_space(8.0);
    ui.label(egui::RichText::new("Controls").strong());
    ui.add_space(8.0);
    
    ui.horizontal(|ui| {
        let can_run = has_mesh && has_bc && !is_running;
        
        ui.add_enabled_ui(can_run, |ui| {
            if ui.button("▶ Run Simulation").clicked() {
                app.start_simulation();
            }
        });
        
        ui.add_enabled_ui(is_running, |ui| {
            if ui.button("■ Stop").clicked() {
                app.stop_simulation();
            }
        });
    });
    
    // Show detailed progress when running
    if is_running {
        ui.add_space(8.0);
        show_solve_progress(ui, app);
    }
    
    // Explicit solver options
    if app.state.simulation_config.solver == SolverType::Explicit {
        ui.add_space(16.0);
        ui.separator();
        ui.add_space(8.0);
        ui.label(egui::RichText::new("Explicit Solver Options").strong());
        ui.add_space(4.0);
        
        // Copy values to avoid borrow issues
        let mut time_steps = app.state.simulation_config.explicit_settings.time_steps as i32;
        let mut vtk_interval = app.state.simulation_config.explicit_settings.vtk_save_steps as i32;
        let mut state_interval = app.state.simulation_config.explicit_settings.state_save_steps as i32;
        
        let mut changed = false;
        
        egui::Grid::new("explicit_opts")
            .num_columns(2)
            .spacing([20.0, 4.0])
            .show(ui, |ui| {
                ui.label("Time Steps:");
                if ui.add(egui::DragValue::new(&mut time_steps).range(1..=10000000)).changed() {
                    changed = true;
                }
                ui.end_row();
                
                ui.label("VTK Save Interval:");
                if ui.add(egui::DragValue::new(&mut vtk_interval).range(1..=100000)).changed() {
                    changed = true;
                }
                ui.end_row();
                
                ui.label("State Save Interval:");
                if ui.add(egui::DragValue::new(&mut state_interval).range(1..=100000)).changed() {
                    changed = true;
                }
                ui.end_row();
            });
        
        if changed {
            app.state.simulation_config.explicit_settings.time_steps = time_steps as usize;
            app.state.simulation_config.explicit_settings.vtk_save_steps = vtk_interval as usize;
            app.state.simulation_config.explicit_settings.state_save_steps = state_interval as usize;
        }
        
        // Estimate
        ui.add_space(8.0);
        let estimated_outputs = time_steps / vtk_interval;
        ui.label(format!("Estimated VTK outputs: {}", estimated_outputs));
    }
    
    // Previous results info with phase timing
    if let Some(results) = &app.state.results {
        ui.add_space(16.0);
        ui.separator();
        ui.add_space(8.0);
        ui.label(egui::RichText::new("Last Results").strong());
        ui.add_space(4.0);
        
        egui::Grid::new("last_results")
            .num_columns(2)
            .spacing([20.0, 4.0])
            .show(ui, |ui| {
                ui.label("Max Displacement:");
                ui.label(format!("{:.6e}", results.stats.max_displacement));
                ui.end_row();
                
                ui.label("Total Solve Time:");
                ui.label(format_duration(results.stats.solver_time_ms));
                ui.end_row();
            });
        
        // Show phase timing breakdown if available
        if let Some(ref timing) = results.stats.phase_timing {
            ui.add_space(8.0);
            ui.collapsing("Phase Timing Breakdown", |ui| {
                show_phase_timing_breakdown(ui, timing);
            });
        }
        
        ui.add_space(8.0);
        if ui.button("View Results →").clicked() {
            app.state.ui_state.active_panel = crate::state::ActivePanel::Results;
        }
    }
}

/// Show detailed solve progress during simulation
fn show_solve_progress(ui: &mut egui::Ui, app: &FeaApp) {
    let progress = &app.state.solve_progress;
    
    // Overall progress bar
    ui.add(
        egui::ProgressBar::new(app.state.progress)
            .show_percentage()
            .animate(true)
    );
    
    // Elapsed time
    ui.horizontal(|ui| {
        ui.label("Elapsed:");
        ui.label(format_duration(progress.elapsed_ms));
        
        // Estimated remaining (if available)
        if let Some(remaining) = progress.estimated_remaining_ms {
            ui.separator();
            ui.label("ETA:");
            ui.label(format_duration(remaining));
        }
    });
    
    // Current phase with icon
    if let Some(ref phase) = progress.current_phase {
        ui.horizontal(|ui| {
            if let Some(cat) = progress.current_category {
                let color = category_color(cat);
                ui.label(egui::RichText::new(category_icon(cat)).color(color));
                ui.label(egui::RichText::new(phase).strong());
            } else {
                ui.label(egui::RichText::new("⏳").color(egui::Color32::YELLOW));
                ui.label(egui::RichText::new(phase).strong());
            }
        });
    }
    
    // Step progress for explicit solver
    if let Some((current, total)) = progress.step_progress {
        ui.add_space(4.0);
        let step_progress = current as f32 / total as f32;
        ui.add(
            egui::ProgressBar::new(step_progress)
                .text(format!("Step {}/{}", current, total))
        );
    }
    
    // Phase timeline (completed phases)
    if !progress.completed_phases.is_empty() {
        ui.add_space(8.0);
        ui.collapsing("Completed Phases", |ui| {
            egui::ScrollArea::vertical()
                .max_height(150.0)
                .show(ui, |ui| {
                    for entry in &progress.completed_phases {
                        show_phase_entry(ui, entry);
                    }
                });
        });
    }
}

/// Show a single phase entry in the log
fn show_phase_entry(ui: &mut egui::Ui, entry: &crate::state::SolvePhaseEntry) {
    ui.horizontal(|ui| {
        let color = category_color(entry.category);
        ui.label(egui::RichText::new(category_icon(entry.category)).color(color));
        ui.label(&entry.name);
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(egui::RichText::new(format_duration(entry.duration_ms)).monospace());
        });
    });
    if let Some(ref details) = entry.details {
        ui.indent("phase_detail", |ui| {
            ui.label(egui::RichText::new(details).small().weak());
        });
    }
}

/// Show phase timing breakdown from completed results
fn show_phase_timing_breakdown(ui: &mut egui::Ui, timing: &crate::state::SolvePhaseTimings) {
    // Waterfall/timeline visualization
    let total_ms = timing.total_ms.max(1) as f32;
    let available_width = ui.available_width() - 100.0; // Leave room for labels
    
    // Group by category for summary
    let mut category_totals: std::collections::HashMap<SolvePhaseCategory, u64> = std::collections::HashMap::new();
    for phase in &timing.phases {
        *category_totals.entry(phase.category).or_insert(0) += phase.duration_ms;
    }
    
    // Category summary bar
    ui.label("Time by Category:");
    let (rect, _response) = ui.allocate_exact_size(
        egui::vec2(available_width, 20.0),
        egui::Sense::hover()
    );
    
    let painter = ui.painter();
    let mut x_offset = rect.left();
    
    for cat in &[SolvePhaseCategory::Init, SolvePhaseCategory::Assembly, 
                 SolvePhaseCategory::Solve, SolvePhaseCategory::PostProcess] {
        if let Some(&cat_ms) = category_totals.get(cat) {
            let width = (cat_ms as f32 / total_ms) * available_width;
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
    
    // Legend
    ui.horizontal(|ui| {
        for cat in &[SolvePhaseCategory::Init, SolvePhaseCategory::Assembly, 
                     SolvePhaseCategory::Solve, SolvePhaseCategory::PostProcess] {
            if let Some(&cat_ms) = category_totals.get(cat) {
                let pct = (cat_ms as f32 / total_ms) * 100.0;
                let color = category_color(*cat);
                ui.colored_label(color, format!("{} {:.0}%", cat.name(), pct));
            }
        }
    });
    
    ui.add_space(8.0);
    
    // Detailed phase list
    ui.label("Phase Details:");
    egui::ScrollArea::vertical()
        .max_height(200.0)
        .show(ui, |ui| {
            egui::Grid::new("phase_details_grid")
                .num_columns(3)
                .spacing([10.0, 4.0])
                .striped(true)
                .show(ui, |ui| {
                    // Header
                    ui.label(egui::RichText::new("Phase").strong());
                    ui.label(egui::RichText::new("Duration").strong());
                    ui.label(egui::RichText::new("% of Total").strong());
                    ui.end_row();
                    
                    for phase in &timing.phases {
                        // Phase name with category icon
                        ui.horizontal(|ui| {
                            let color = category_color(phase.category);
                            ui.label(egui::RichText::new(category_icon(phase.category)).color(color));
                            ui.label(&phase.name);
                        });
                        
                        // Duration
                        ui.label(format_duration(phase.duration_ms));
                        
                        // Percentage
                        let pct = (phase.duration_ms as f32 / total_ms) * 100.0;
                        ui.label(format!("{:.1}%", pct));
                        ui.end_row();
                        
                        // Details (if any) on next row
                        if let Some(ref details) = phase.details {
                            ui.label("");
                            ui.add(egui::Label::new(
                                egui::RichText::new(details).small().weak()
                            ));
                            ui.label("");
                            ui.end_row();
                        }
                    }
                });
        });
}

fn show_check(ui: &mut egui::Ui, label: &str, passed: bool) {
    ui.horizontal(|ui| {
        if passed {
            ui.label(egui::RichText::new("✓").color(egui::Color32::GREEN));
            ui.label(label);
        } else {
            ui.label(egui::RichText::new("○").color(egui::Color32::GRAY));
            ui.label(egui::RichText::new(label).color(egui::Color32::GRAY));
        }
    });
}

/// Format duration in a human-readable way
fn format_duration(ms: u64) -> String {
    if ms < 1000 {
        format!("{}ms", ms)
    } else if ms < 60_000 {
        format!("{:.2}s", ms as f64 / 1000.0)
    } else {
        let secs = ms / 1000;
        let mins = secs / 60;
        let remaining_secs = secs % 60;
        format!("{}m {}s", mins, remaining_secs)
    }
}

/// Get icon for a phase category
fn category_icon(cat: SolvePhaseCategory) -> &'static str {
    match cat {
        SolvePhaseCategory::Init => "◆",
        SolvePhaseCategory::Assembly => "▤",
        SolvePhaseCategory::Solve => "⚡",
        SolvePhaseCategory::PostProcess => "◈",
        SolvePhaseCategory::Other => "○",
    }
}

/// Get color for a phase category
fn category_color(cat: SolvePhaseCategory) -> egui::Color32 {
    let rgba = cat.color();
    egui::Color32::from_rgba_unmultiplied(rgba[0], rgba[1], rgba[2], rgba[3])
}
