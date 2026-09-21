//! Left side panel with workflow tabs and content

use eframe::egui;
use crate::app::FeaApp;
use crate::state::ActivePanel;
use super::{mesh_panel, setup_panel, run_panel, results_panel};

pub fn show(ctx: &egui::Context, app: &mut FeaApp) {
    egui::SidePanel::left("side_panel")
        .resizable(true)
        .default_width(350.0)
        .min_width(280.0)
        .max_width(500.0)
        .show(ctx, |ui| {
            ui.add_space(8.0);
            
            // Workflow tabs with enhanced tooltips
            ui.horizontal(|ui| {
                // Mesh tab
                ui.selectable_value(
                    &mut app.state.ui_state.active_panel,
                    ActivePanel::Mesh,
                    "Mesh"
                ).on_hover_text("Import or create meshes (Ctrl+1)");
                
                // Setup tab
                ui.selectable_value(
                    &mut app.state.ui_state.active_panel,
                    ActivePanel::Setup,
                    "Setup"
                ).on_hover_text("Define materials and boundary conditions (Ctrl+2)");
                
                // Run tab
                let run_text = if app.state.is_running { "Run..." } else { "Run" };
                ui.selectable_value(
                    &mut app.state.ui_state.active_panel,
                    ActivePanel::Run,
                    run_text
                ).on_hover_text("Configure and run simulation (Ctrl+3)");
                
                // Results tab
                ui.selectable_value(
                    &mut app.state.ui_state.active_panel,
                    ActivePanel::Results,
                    "Results"
                ).on_hover_text("View and export results (Ctrl+4)");
            });
            
            ui.separator();
            
            // Panel content
            egui::ScrollArea::vertical()
                .auto_shrink([false; 2])
                .show(ui, |ui| {
                    match app.state.ui_state.active_panel {
                        ActivePanel::Mesh => mesh_panel::show(ui, app),
                        ActivePanel::Setup => setup_panel::show(ui, app),
                        ActivePanel::Run => run_panel::show(ui, app),
                        ActivePanel::Results => results_panel::show(ui, app),
                    }
                });
        });
}
