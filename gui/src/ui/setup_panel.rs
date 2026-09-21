//! Setup panel for simulation configuration (BCs, materials, solver settings)

use eframe::egui;
use crate::app::FeaApp;
use crate::state::{
    BoundaryConditionConfig, FixedBcConfig, LoadBcConfig, TorqueBcConfig, 
    ContactBcConfig, PressureBcConfig, TractionBcConfig, TractionTypeConfig,
    BodyForceBcConfig, BodyForceTypeConfig,
    MaterialConfig, SolverType, NewBcType
};

pub fn show(ui: &mut egui::Ui, app: &mut FeaApp) {
    ui.heading("Simulation Setup");
    ui.add_space(8.0);
    
    if app.state.current_mesh().is_none() {
        ui.label("⚠ No mesh loaded. Import a mesh first.");
        return;
    }
    
    // Solver settings
    ui.collapsing("Solver Settings", |ui| {
        show_solver_settings(ui, app);
    });
    
    ui.add_space(8.0);
    
    // Materials
    ui.collapsing("Materials", |ui| {
        show_materials(ui, app);
    });
    
    ui.add_space(8.0);
    
    // Boundary Conditions
    ui.collapsing("Boundary Conditions", |ui| {
        show_boundary_conditions(ui, app);
    });
    
    ui.add_space(8.0);
    
    // Output settings
    ui.collapsing("Output Settings", |ui| {
        show_output_settings(ui, app);
    });
}

fn show_solver_settings(ui: &mut egui::Ui, app: &mut FeaApp) {
    let config = &mut app.state.simulation_config;
    
    ui.horizontal(|ui| {
        ui.label("Solver Type:");
        egui::ComboBox::from_id_salt("solver_type")
            .selected_text(match config.solver {
                SolverType::Direct => "Direct",
                SolverType::Explicit => "Explicit",
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut config.solver, SolverType::Direct, "Direct");
                ui.selectable_value(&mut config.solver, SolverType::Explicit, "Explicit");
            });
    });
    
    ui.horizontal(|ui| {
        ui.label("DOFs per Node:");
        let mut dofs = config.dofs as i32;
        if ui.add(egui::DragValue::new(&mut dofs).range(1..=6)).changed() {
            config.dofs = dofs as usize;
        }
    });
    
    // Explicit solver settings
    if config.solver == SolverType::Explicit {
        ui.add_space(8.0);
        ui.label("Explicit Solver Settings:");
        
        ui.horizontal(|ui| {
            ui.label("Time Steps:");
            let mut steps = config.explicit_settings.time_steps as i32;
            if ui.add(egui::DragValue::new(&mut steps).range(1..=1000000)).changed() {
                config.explicit_settings.time_steps = steps as usize;
            }
        });
        
        ui.horizontal(|ui| {
            ui.label("VTK Save Interval:");
            let mut interval = config.explicit_settings.vtk_save_steps as i32;
            if ui.add(egui::DragValue::new(&mut interval).range(1..=100000)).changed() {
                config.explicit_settings.vtk_save_steps = interval as usize;
            }
        });
    }
}

fn show_materials(ui: &mut egui::Ui, app: &mut FeaApp) {
    let config = &mut app.state.simulation_config;
    
    // Add material button
    if ui.button("➕ Add Material").clicked() {
        let new_id = config.materials.len() + 1;
        config.materials.push(MaterialConfig {
            id: new_id,
            name: format!("Material {}", new_id),
            ..Default::default()
        });
    }
    
    ui.add_space(8.0);
    
    // Material list
    let mut remove_idx = None;
    let num_materials = config.materials.len();
    
    for idx in 0..num_materials {
        ui.push_id(idx, |ui| {
            let material = &mut config.materials[idx];
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.strong(&material.name);
                    if ui.small_button("🗑").clicked() && num_materials > 1 {
                        remove_idx = Some(idx);
                    }
                });
                
                egui::Grid::new("material_props")
                    .num_columns(2)
                    .spacing([10.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Name:");
                        ui.text_edit_singleline(&mut material.name);
                        ui.end_row();
                        
                        ui.label("Young's Modulus (Pa):");
                        ui.add(egui::DragValue::new(&mut material.youngs_modulus)
                            .speed(1e9)
                            .range(1e3..=1e15));
                        ui.end_row();
                        
                        ui.label("Poisson's Ratio:");
                        ui.add(egui::DragValue::new(&mut material.poissons_ratio)
                            .speed(0.01)
                            .range(0.0..=0.5));
                        ui.end_row();
                        
                        ui.label("Density (kg/m³):");
                        ui.add(egui::DragValue::new(&mut material.density)
                            .speed(100.0)
                            .range(1.0..=50000.0));
                        ui.end_row();
                    });
                
                // Preset buttons
                ui.horizontal(|ui| {
                    if ui.small_button("Steel").clicked() {
                        material.youngs_modulus = 200e9;
                        material.poissons_ratio = 0.3;
                        material.density = 7850.0;
                    }
                    if ui.small_button("Aluminum").clicked() {
                        material.youngs_modulus = 68.9e9;
                        material.poissons_ratio = 0.33;
                        material.density = 2700.0;
                    }
                    if ui.small_button("Titanium").clicked() {
                        material.youngs_modulus = 116e9;
                        material.poissons_ratio = 0.34;
                        material.density = 4500.0;
                    }
                });
            });
        });
        
        ui.add_space(4.0);
    }
    
    if let Some(idx) = remove_idx {
        config.materials.remove(idx);
    }
}

fn show_boundary_conditions(ui: &mut egui::Ui, app: &mut FeaApp) {
    // Get available node groups
    let node_groups: Vec<String> = app.state.current_mesh()
        .map(|m| m.mesh.node_groups.keys().cloned().collect())
        .unwrap_or_default();
    
    let element_groups: Vec<String> = app.state.current_mesh()
        .map(|m| m.mesh.element_groups.keys().cloned().collect())
        .unwrap_or_default();
    
    // Add BC controls
    ui.horizontal(|ui| {
        ui.label("Add BC:");
        egui::ComboBox::from_id_salt("new_bc_type")
            .selected_text(match app.state.ui_state.bc_editor.new_bc_type {
                NewBcType::Fixed => "Fixed",
                NewBcType::Load => "Load",
                NewBcType::Torque => "Torque",
                NewBcType::Contact => "Contact",
                NewBcType::Pressure => "Pressure",
                NewBcType::Traction => "Traction",
                NewBcType::BodyForce => "Body Force",
            })
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Fixed, "Fixed");
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Load, "Load");
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Torque, "Torque");
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Contact, "Contact");
                ui.separator();
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Pressure, "⬇ Pressure (surface)");
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::Traction, "↗ Traction (surface)");
                ui.selectable_value(&mut app.state.ui_state.bc_editor.new_bc_type, NewBcType::BodyForce, "🌍 Body Force (volume)");
            });
        
        if ui.button("➕").clicked() {
            let new_bc = match app.state.ui_state.bc_editor.new_bc_type {
                NewBcType::Fixed => BoundaryConditionConfig::Fixed(FixedBcConfig::default()),
                NewBcType::Load => BoundaryConditionConfig::Load(LoadBcConfig::default()),
                NewBcType::Torque => BoundaryConditionConfig::Torque(TorqueBcConfig::default()),
                NewBcType::Contact => BoundaryConditionConfig::Contact(ContactBcConfig::default()),
                NewBcType::Pressure => BoundaryConditionConfig::Pressure(PressureBcConfig::default()),
                NewBcType::Traction => BoundaryConditionConfig::Traction(TractionBcConfig::default()),
                NewBcType::BodyForce => BoundaryConditionConfig::BodyForce(BodyForceBcConfig::default()),
            };
            app.state.simulation_config.boundary_conditions.push(new_bc);
        }
    });
    
    ui.add_space(8.0);
    
    // BC list
    let mut remove_idx = None;
    let config = &mut app.state.simulation_config;
    
    for (idx, bc) in config.boundary_conditions.iter_mut().enumerate() {
        ui.push_id(idx, |ui| {
            ui.group(|ui| {
                match bc {
                    BoundaryConditionConfig::Fixed(cfg) => {
                        show_fixed_bc(ui, cfg, &node_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::Load(cfg) => {
                        show_load_bc(ui, cfg, &node_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::Torque(cfg) => {
                        show_torque_bc(ui, cfg, &node_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::Contact(cfg) => {
                        show_contact_bc(ui, cfg, &node_groups, &element_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::Pressure(cfg) => {
                        show_pressure_bc(ui, cfg, &element_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::Traction(cfg) => {
                        show_traction_bc(ui, cfg, &element_groups, &mut remove_idx, idx);
                    }
                    BoundaryConditionConfig::BodyForce(cfg) => {
                        show_body_force_bc(ui, cfg, &element_groups, &mut remove_idx, idx);
                    }
                }
            });
        });
        ui.add_space(4.0);
    }
    
    if let Some(idx) = remove_idx {
        config.boundary_conditions.remove(idx);
    }
}

fn show_fixed_bc(
    ui: &mut egui::Ui, 
    cfg: &mut FixedBcConfig, 
    node_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("🔒 Fixed BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    egui::Grid::new("fixed_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Node Group:");
            egui::ComboBox::from_id_salt("fixed_node_group")
                .selected_text(&cfg.node_group)
                .show_ui(ui, |ui| {
                    for group in node_groups {
                        ui.selectable_value(&mut cfg.node_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            // X constraint
            ui.label("X:");
            let mut x_fixed = cfg.constrain_x.is_some();
            ui.horizontal(|ui| {
                ui.checkbox(&mut x_fixed, "");
                if x_fixed {
                    let mut val = cfg.constrain_x.unwrap_or(0.0);
                    ui.add(egui::DragValue::new(&mut val).speed(0.001));
                    cfg.constrain_x = Some(val);
                } else {
                    cfg.constrain_x = None;
                    ui.label("Free");
                }
            });
            ui.end_row();
            
            // Y constraint
            ui.label("Y:");
            let mut y_fixed = cfg.constrain_y.is_some();
            ui.horizontal(|ui| {
                ui.checkbox(&mut y_fixed, "");
                if y_fixed {
                    let mut val = cfg.constrain_y.unwrap_or(0.0);
                    ui.add(egui::DragValue::new(&mut val).speed(0.001));
                    cfg.constrain_y = Some(val);
                } else {
                    cfg.constrain_y = None;
                    ui.label("Free");
                }
            });
            ui.end_row();
            
            // Z constraint
            ui.label("Z:");
            let mut z_fixed = cfg.constrain_z.is_some();
            ui.horizontal(|ui| {
                ui.checkbox(&mut z_fixed, "");
                if z_fixed {
                    let mut val = cfg.constrain_z.unwrap_or(0.0);
                    ui.add(egui::DragValue::new(&mut val).speed(0.001));
                    cfg.constrain_z = Some(val);
                } else {
                    cfg.constrain_z = None;
                    ui.label("Free");
                }
            });
            ui.end_row();
        });
}

fn show_load_bc(
    ui: &mut egui::Ui, 
    cfg: &mut LoadBcConfig, 
    node_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("⬇ Load BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    egui::Grid::new("load_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Node Group:");
            egui::ComboBox::from_id_salt("load_node_group")
                .selected_text(&cfg.node_group)
                .show_ui(ui, |ui| {
                    for group in node_groups {
                        ui.selectable_value(&mut cfg.node_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            ui.label("Force X (N):");
            ui.add(egui::DragValue::new(&mut cfg.force_x).speed(100.0));
            ui.end_row();
            
            ui.label("Force Y (N):");
            ui.add(egui::DragValue::new(&mut cfg.force_y).speed(100.0));
            ui.end_row();
            
            ui.label("Force Z (N):");
            ui.add(egui::DragValue::new(&mut cfg.force_z).speed(100.0));
            ui.end_row();
        });
}

fn show_torque_bc(
    ui: &mut egui::Ui, 
    cfg: &mut TorqueBcConfig, 
    node_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("🔄 Torque BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    egui::Grid::new("torque_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Node Group:");
            egui::ComboBox::from_id_salt("torque_node_group")
                .selected_text(&cfg.node_group)
                .show_ui(ui, |ui| {
                    for group in node_groups {
                        ui.selectable_value(&mut cfg.node_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            ui.label("Axis Point:");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut cfg.axis_point[0]).prefix("x: ").speed(0.1));
                ui.add(egui::DragValue::new(&mut cfg.axis_point[1]).prefix("y: ").speed(0.1));
                ui.add(egui::DragValue::new(&mut cfg.axis_point[2]).prefix("z: ").speed(0.1));
            });
            ui.end_row();
            
            ui.label("Axis Direction:");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut cfg.axis_direction[0]).prefix("x: ").speed(0.1));
                ui.add(egui::DragValue::new(&mut cfg.axis_direction[1]).prefix("y: ").speed(0.1));
                ui.add(egui::DragValue::new(&mut cfg.axis_direction[2]).prefix("z: ").speed(0.1));
            });
            ui.end_row();
            
            ui.label("Magnitude (N·m):");
            ui.add(egui::DragValue::new(&mut cfg.magnitude).speed(100.0));
            ui.end_row();
        });
}

fn show_contact_bc(
    ui: &mut egui::Ui, 
    cfg: &mut ContactBcConfig, 
    node_groups: &[String],
    element_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("👆 Contact BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    egui::Grid::new("contact_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Primary Surface:");
            egui::ComboBox::from_id_salt("primary_surface")
                .selected_text(&cfg.primary_surface)
                .show_ui(ui, |ui| {
                    for group in node_groups.iter().chain(element_groups.iter()) {
                        ui.selectable_value(&mut cfg.primary_surface, group.clone(), group);
                    }
                });
            ui.end_row();
            
            ui.label("Secondary Surface:");
            egui::ComboBox::from_id_salt("secondary_surface")
                .selected_text(&cfg.secondary_surface)
                .show_ui(ui, |ui| {
                    for group in node_groups.iter().chain(element_groups.iter()) {
                        ui.selectable_value(&mut cfg.secondary_surface, group.clone(), group);
                    }
                });
            ui.end_row();
        });
}

fn show_pressure_bc(
    ui: &mut egui::Ui, 
    cfg: &mut PressureBcConfig, 
    element_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("⬇ Pressure BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    ui.label("Applied normal to surface elements (F = ∫∫ p·n dA)");
    
    egui::Grid::new("pressure_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Element Group:");
            egui::ComboBox::from_id_salt("pressure_element_group")
                .selected_text(if cfg.element_group.is_empty() { "(select)" } else { &cfg.element_group })
                .show_ui(ui, |ui| {
                    for group in element_groups {
                        ui.selectable_value(&mut cfg.element_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            ui.label("Pressure (Pa):");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut cfg.pressure).speed(1e5));
                ui.label(if cfg.pressure > 0.0 { "(compression)" } else { "(tension)" });
            });
            ui.end_row();
        });
    
    // Preset buttons
    ui.horizontal(|ui| {
        if ui.small_button("1 kPa").clicked() { cfg.pressure = 1e3; }
        if ui.small_button("1 MPa").clicked() { cfg.pressure = 1e6; }
        if ui.small_button("100 MPa").clicked() { cfg.pressure = 100e6; }
    });
}

fn show_traction_bc(
    ui: &mut egui::Ui, 
    cfg: &mut TractionBcConfig, 
    element_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("↗ Traction BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    ui.label("Surface force per unit area (F = ∫∫ t dA)");
    
    egui::Grid::new("traction_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Element Group:");
            egui::ComboBox::from_id_salt("traction_element_group")
                .selected_text(if cfg.element_group.is_empty() { "(select)" } else { &cfg.element_group })
                .show_ui(ui, |ui| {
                    for group in element_groups {
                        ui.selectable_value(&mut cfg.element_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            // Traction type selector
            ui.label("Type:");
            let type_name = match &cfg.traction_type {
                TractionTypeConfig::Uniform { .. } => "Uniform",
                TractionTypeConfig::Normal { .. } => "Normal",
                TractionTypeConfig::Shear { .. } => "Shear",
            };
            egui::ComboBox::from_id_salt("traction_type")
                .selected_text(type_name)
                .show_ui(ui, |ui| {
                    if ui.selectable_label(matches!(&cfg.traction_type, TractionTypeConfig::Uniform { .. }), "Uniform").clicked() {
                        cfg.traction_type = TractionTypeConfig::Uniform { fx: 0.0, fy: 0.0, fz: -1e6 };
                    }
                    if ui.selectable_label(matches!(&cfg.traction_type, TractionTypeConfig::Normal { .. }), "Normal").clicked() {
                        cfg.traction_type = TractionTypeConfig::Normal { magnitude: 1e6 };
                    }
                    if ui.selectable_label(matches!(&cfg.traction_type, TractionTypeConfig::Shear { .. }), "Shear").clicked() {
                        cfg.traction_type = TractionTypeConfig::Shear { magnitude: 1e5 };
                    }
                });
            ui.end_row();
        });
    
    // Type-specific parameters
    match &mut cfg.traction_type {
        TractionTypeConfig::Uniform { fx, fy, fz } => {
            egui::Grid::new("traction_uniform_grid")
                .num_columns(2)
                .spacing([10.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Traction X (Pa):");
                    ui.add(egui::DragValue::new(fx).speed(1e4));
                    ui.end_row();
                    
                    ui.label("Traction Y (Pa):");
                    ui.add(egui::DragValue::new(fy).speed(1e4));
                    ui.end_row();
                    
                    ui.label("Traction Z (Pa):");
                    ui.add(egui::DragValue::new(fz).speed(1e4));
                    ui.end_row();
                });
        }
        TractionTypeConfig::Normal { magnitude } => {
            ui.horizontal(|ui| {
                ui.label("Magnitude (Pa):");
                ui.add(egui::DragValue::new(magnitude).speed(1e4));
                ui.label(if *magnitude > 0.0 { "(tension)" } else { "(compression)" });
            });
        }
        TractionTypeConfig::Shear { magnitude } => {
            ui.horizontal(|ui| {
                ui.label("Shear Magnitude (Pa):");
                ui.add(egui::DragValue::new(magnitude).speed(1e4));
            });
        }
    }
}

fn show_body_force_bc(
    ui: &mut egui::Ui, 
    cfg: &mut BodyForceBcConfig, 
    element_groups: &[String],
    remove_idx: &mut Option<usize>,
    idx: usize
) {
    ui.horizontal(|ui| {
        ui.strong("🌍 Body Force BC");
        if ui.small_button("🗑").clicked() {
            *remove_idx = Some(idx);
        }
    });
    
    ui.label("Distributed force per unit volume (F = ∫∫∫ b·N dV)");
    
    egui::Grid::new("bodyforce_bc_grid")
        .num_columns(2)
        .spacing([10.0, 4.0])
        .show(ui, |ui| {
            ui.label("Name:");
            ui.text_edit_singleline(&mut cfg.name);
            ui.end_row();
            
            ui.label("Element Group:");
            egui::ComboBox::from_id_salt("bodyforce_element_group")
                .selected_text(if cfg.element_group.is_empty() { "(all elements)" } else { &cfg.element_group })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut cfg.element_group, String::new(), "(all elements)");
                    for group in element_groups {
                        ui.selectable_value(&mut cfg.element_group, group.clone(), group);
                    }
                });
            ui.end_row();
            
            // Force type selector
            ui.label("Type:");
            let type_name = match &cfg.force_type {
                BodyForceTypeConfig::Gravity { .. } => "Gravity",
                BodyForceTypeConfig::Centrifugal { .. } => "Centrifugal",
                BodyForceTypeConfig::Uniform { .. } => "Uniform",
            };
            egui::ComboBox::from_id_salt("bodyforce_type")
                .selected_text(type_name)
                .show_ui(ui, |ui| {
                    if ui.selectable_label(matches!(&cfg.force_type, BodyForceTypeConfig::Gravity { .. }), "Gravity").clicked() {
                        cfg.force_type = BodyForceTypeConfig::Gravity { gx: 0.0, gy: -9.81, gz: 0.0 };
                    }
                    if ui.selectable_label(matches!(&cfg.force_type, BodyForceTypeConfig::Centrifugal { .. }), "Centrifugal").clicked() {
                        cfg.force_type = BodyForceTypeConfig::Centrifugal {
                            axis_point: [0.0, 0.0, 0.0],
                            axis_direction: [0.0, 1.0, 0.0],
                            angular_velocity: 100.0,
                        };
                    }
                    if ui.selectable_label(matches!(&cfg.force_type, BodyForceTypeConfig::Uniform { .. }), "Uniform").clicked() {
                        cfg.force_type = BodyForceTypeConfig::Uniform { fx: 0.0, fy: 0.0, fz: 0.0 };
                    }
                });
            ui.end_row();
        });
    
    // Type-specific parameters
    match &mut cfg.force_type {
        BodyForceTypeConfig::Gravity { gx, gy, gz } => {
            egui::Grid::new("gravity_grid")
                .num_columns(2)
                .spacing([10.0, 4.0])
                .show(ui, |ui| {
                    ui.label("g_x (m/s²):");
                    ui.add(egui::DragValue::new(gx).speed(0.1));
                    ui.end_row();
                    
                    ui.label("g_y (m/s²):");
                    ui.add(egui::DragValue::new(gy).speed(0.1));
                    ui.end_row();
                    
                    ui.label("g_z (m/s²):");
                    ui.add(egui::DragValue::new(gz).speed(0.1));
                    ui.end_row();
                });
            
            // Preset buttons
            ui.horizontal(|ui| {
                if ui.small_button("-Y (Earth)").clicked() { *gx = 0.0; *gy = -9.81; *gz = 0.0; }
                if ui.small_button("-Z (Earth)").clicked() { *gx = 0.0; *gy = 0.0; *gz = -9.81; }
                if ui.small_button("Moon").clicked() { *gx = 0.0; *gy = -1.62; *gz = 0.0; }
            });
        }
        BodyForceTypeConfig::Centrifugal { axis_point, axis_direction, angular_velocity } => {
            egui::Grid::new("centrifugal_grid")
                .num_columns(2)
                .spacing([10.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Axis Point:");
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(&mut axis_point[0]).prefix("x: ").speed(0.1));
                        ui.add(egui::DragValue::new(&mut axis_point[1]).prefix("y: ").speed(0.1));
                        ui.add(egui::DragValue::new(&mut axis_point[2]).prefix("z: ").speed(0.1));
                    });
                    ui.end_row();
                    
                    ui.label("Axis Direction:");
                    ui.horizontal(|ui| {
                        ui.add(egui::DragValue::new(&mut axis_direction[0]).prefix("x: ").speed(0.1));
                        ui.add(egui::DragValue::new(&mut axis_direction[1]).prefix("y: ").speed(0.1));
                        ui.add(egui::DragValue::new(&mut axis_direction[2]).prefix("z: ").speed(0.1));
                    });
                    ui.end_row();
                    
                    ui.label("ω (rad/s):");
                    ui.add(egui::DragValue::new(angular_velocity).speed(1.0));
                    ui.end_row();
                    
                    // Show RPM conversion
                    let rpm = *angular_velocity * 60.0 / (2.0 * std::f64::consts::PI);
                    ui.label("RPM:");
                    ui.label(format!("{:.1}", rpm));
                    ui.end_row();
                });
        }
        BodyForceTypeConfig::Uniform { fx, fy, fz } => {
            egui::Grid::new("uniform_body_grid")
                .num_columns(2)
                .spacing([10.0, 4.0])
                .show(ui, |ui| {
                    ui.label("Force X (N/m³):");
                    ui.add(egui::DragValue::new(fx).speed(1e3));
                    ui.end_row();
                    
                    ui.label("Force Y (N/m³):");
                    ui.add(egui::DragValue::new(fy).speed(1e3));
                    ui.end_row();
                    
                    ui.label("Force Z (N/m³):");
                    ui.add(egui::DragValue::new(fz).speed(1e3));
                    ui.end_row();
                });
        }
    }
}

fn show_output_settings(ui: &mut egui::Ui, app: &mut FeaApp) {
    let config = &mut app.state.simulation_config;
    
    ui.horizontal(|ui| {
        ui.label("VTK Output:");
        if let Some(path) = &config.output.vtk_path {
            ui.label(path.to_string_lossy().to_string());
        } else {
            ui.label("(not set)");
        }
        if ui.button("Browse...").clicked() {
            #[cfg(feature = "native")]
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("VTK", &["vtk"])
                .save_file()
            {
                config.output.vtk_path = Some(path);
            }
            
            #[cfg(not(feature = "native"))]
            {
                // VTK output not available on web
            }
        }
    });
    
    ui.checkbox(
        &mut config.output.save_stiffness_matrix, 
        "Save Stiffness Matrix"
    );
}
