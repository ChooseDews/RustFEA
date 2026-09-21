//! Project file serialization/deserialization
//! 
//! Handles saving and loading simulation configurations in a human-readable TOML format.

use serde::{Serialize, Deserialize};
use std::path::Path;
use crate::state::{
    SimulationConfig, SolverType, MaterialConfig, 
    BoundaryConditionConfig, FixedBcConfig, LoadBcConfig, TorqueBcConfig,
    ContactBcConfig, PressureBcConfig, TractionBcConfig, TractionTypeConfig,
    BodyForceBcConfig, BodyForceTypeConfig, ExplicitSettings,
};

/// Serializable project file format
#[derive(Serialize, Deserialize, Debug)]
pub struct ProjectFile {
    /// Project name
    pub name: String,
    
    /// Project description
    #[serde(default)]
    pub description: String,
    
    /// Mesh file path (relative to project file)
    pub mesh: Option<String>,
    
    /// Simulation configuration
    pub simulation: SimulationConfigSer,
}

/// Serializable simulation config
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SimulationConfigSer {
    /// Solver type: "direct" or "explicit"
    pub solver: String,
    
    /// DOFs per node
    #[serde(default = "default_dofs")]
    pub dofs: usize,
    
    /// Materials
    #[serde(default)]
    pub materials: Vec<MaterialConfigSer>,
    
    /// Boundary conditions
    #[serde(default)]
    pub boundary_conditions: Vec<BoundaryConditionSer>,
    
    /// Explicit solver settings (if applicable)
    #[serde(default)]
    pub explicit: Option<ExplicitSettingsSer>,
}

fn default_dofs() -> usize { 3 }

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MaterialConfigSer {
    pub id: usize,
    pub name: String,
    pub youngs_modulus: f64,
    pub poissons_ratio: f64,
    pub density: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ExplicitSettingsSer {
    pub time_steps: usize,
    #[serde(default)]
    pub time_step_override: Option<f64>,
    #[serde(default = "default_save_steps")]
    pub vtk_save_steps: usize,
    #[serde(default = "default_save_steps")]
    pub state_save_steps: usize,
}

fn default_save_steps() -> usize { 100 }

impl Default for ExplicitSettingsSer {
    fn default() -> Self {
        Self {
            time_steps: 1000,
            time_step_override: None,
            vtk_save_steps: 100,
            state_save_steps: 100,
        }
    }
}

/// Serializable boundary condition (tagged union)
#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "type")]
pub enum BoundaryConditionSer {
    #[serde(rename = "fixed")]
    Fixed {
        name: String,
        node_group: String,
        #[serde(default)]
        x: Option<f64>,
        #[serde(default)]
        y: Option<f64>,
        #[serde(default)]
        z: Option<f64>,
    },
    
    #[serde(rename = "load")]
    Load {
        name: String,
        node_group: String,
        #[serde(default)]
        force_x: f64,
        #[serde(default)]
        force_y: f64,
        #[serde(default)]
        force_z: f64,
    },
    
    #[serde(rename = "torque")]
    Torque {
        name: String,
        node_group: String,
        axis_point: [f64; 3],
        axis_direction: [f64; 3],
        magnitude: f64,
    },
    
    #[serde(rename = "contact")]
    Contact {
        name: String,
        primary_surface: String,
        secondary_surface: String,
    },
    
    #[serde(rename = "pressure")]
    Pressure {
        name: String,
        element_group: String,
        pressure: f64,
    },
    
    #[serde(rename = "traction")]
    Traction {
        name: String,
        element_group: String,
        traction_type: TractionTypeSer,
    },
    
    #[serde(rename = "body_force")]
    BodyForce {
        name: String,
        #[serde(default)]
        element_group: String,
        force_type: BodyForceTypeSer,
    },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "type")]
pub enum TractionTypeSer {
    #[serde(rename = "uniform")]
    Uniform { fx: f64, fy: f64, fz: f64 },
    #[serde(rename = "normal")]
    Normal { magnitude: f64 },
    #[serde(rename = "shear")]
    Shear { magnitude: f64 },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "type")]
pub enum BodyForceTypeSer {
    #[serde(rename = "gravity")]
    Gravity { gx: f64, gy: f64, gz: f64 },
    #[serde(rename = "centrifugal")]
    Centrifugal {
        axis_point: [f64; 3],
        axis_direction: [f64; 3],
        angular_velocity: f64,
    },
    #[serde(rename = "uniform")]
    Uniform { fx: f64, fy: f64, fz: f64 },
}

// ============================================================================
// Conversion functions
// ============================================================================

impl From<&SimulationConfig> for SimulationConfigSer {
    fn from(config: &SimulationConfig) -> Self {
        Self {
            solver: match config.solver {
                SolverType::Direct => "direct".to_string(),
                SolverType::Explicit => "explicit".to_string(),
            },
            dofs: config.dofs,
            materials: config.materials.iter().map(|m| MaterialConfigSer {
                id: m.id,
                name: m.name.clone(),
                youngs_modulus: m.youngs_modulus,
                poissons_ratio: m.poissons_ratio,
                density: m.density,
            }).collect(),
            boundary_conditions: config.boundary_conditions.iter().map(|bc| bc.into()).collect(),
            explicit: if config.solver == SolverType::Explicit {
                Some(ExplicitSettingsSer {
                    time_steps: config.explicit_settings.time_steps,
                    time_step_override: config.explicit_settings.time_step_override,
                    vtk_save_steps: config.explicit_settings.vtk_save_steps,
                    state_save_steps: config.explicit_settings.state_save_steps,
                })
            } else {
                None
            },
        }
    }
}

impl From<&BoundaryConditionConfig> for BoundaryConditionSer {
    fn from(bc: &BoundaryConditionConfig) -> Self {
        match bc {
            BoundaryConditionConfig::Fixed(cfg) => BoundaryConditionSer::Fixed {
                name: cfg.name.clone(),
                node_group: cfg.node_group.clone(),
                x: cfg.constrain_x,
                y: cfg.constrain_y,
                z: cfg.constrain_z,
            },
            BoundaryConditionConfig::Load(cfg) => BoundaryConditionSer::Load {
                name: cfg.name.clone(),
                node_group: cfg.node_group.clone(),
                force_x: cfg.force_x,
                force_y: cfg.force_y,
                force_z: cfg.force_z,
            },
            BoundaryConditionConfig::Torque(cfg) => BoundaryConditionSer::Torque {
                name: cfg.name.clone(),
                node_group: cfg.node_group.clone(),
                axis_point: cfg.axis_point,
                axis_direction: cfg.axis_direction,
                magnitude: cfg.magnitude,
            },
            BoundaryConditionConfig::Contact(cfg) => BoundaryConditionSer::Contact {
                name: cfg.name.clone(),
                primary_surface: cfg.primary_surface.clone(),
                secondary_surface: cfg.secondary_surface.clone(),
            },
            BoundaryConditionConfig::Pressure(cfg) => BoundaryConditionSer::Pressure {
                name: cfg.name.clone(),
                element_group: cfg.element_group.clone(),
                pressure: cfg.pressure,
            },
            BoundaryConditionConfig::Traction(cfg) => BoundaryConditionSer::Traction {
                name: cfg.name.clone(),
                element_group: cfg.element_group.clone(),
                traction_type: match &cfg.traction_type {
                    TractionTypeConfig::Uniform { fx, fy, fz } => 
                        TractionTypeSer::Uniform { fx: *fx, fy: *fy, fz: *fz },
                    TractionTypeConfig::Normal { magnitude } => 
                        TractionTypeSer::Normal { magnitude: *magnitude },
                    TractionTypeConfig::Shear { magnitude } => 
                        TractionTypeSer::Shear { magnitude: *magnitude },
                },
            },
            BoundaryConditionConfig::BodyForce(cfg) => BoundaryConditionSer::BodyForce {
                name: cfg.name.clone(),
                element_group: cfg.element_group.clone(),
                force_type: match &cfg.force_type {
                    BodyForceTypeConfig::Gravity { gx, gy, gz } => 
                        BodyForceTypeSer::Gravity { gx: *gx, gy: *gy, gz: *gz },
                    BodyForceTypeConfig::Centrifugal { axis_point, axis_direction, angular_velocity } =>
                        BodyForceTypeSer::Centrifugal {
                            axis_point: *axis_point,
                            axis_direction: *axis_direction,
                            angular_velocity: *angular_velocity,
                        },
                    BodyForceTypeConfig::Uniform { fx, fy, fz } =>
                        BodyForceTypeSer::Uniform { fx: *fx, fy: *fy, fz: *fz },
                },
            },
        }
    }
}

impl From<&SimulationConfigSer> for SimulationConfig {
    fn from(ser: &SimulationConfigSer) -> Self {
        Self {
            solver: match ser.solver.as_str() {
                "explicit" => SolverType::Explicit,
                _ => SolverType::Direct,
            },
            dofs: ser.dofs,
            materials: ser.materials.iter().map(|m| MaterialConfig {
                id: m.id,
                name: m.name.clone(),
                youngs_modulus: m.youngs_modulus,
                poissons_ratio: m.poissons_ratio,
                density: m.density,
            }).collect(),
            boundary_conditions: ser.boundary_conditions.iter().map(|bc| bc.into()).collect(),
            output: Default::default(),
            explicit_settings: ser.explicit.as_ref().map(|e| ExplicitSettings {
                time_steps: e.time_steps,
                time_step_override: e.time_step_override,
                vtk_save_steps: e.vtk_save_steps,
                state_save_steps: e.state_save_steps,
            }).unwrap_or_default(),
        }
    }
}

impl From<&BoundaryConditionSer> for BoundaryConditionConfig {
    fn from(ser: &BoundaryConditionSer) -> Self {
        match ser {
            BoundaryConditionSer::Fixed { name, node_group, x, y, z } => {
                BoundaryConditionConfig::Fixed(FixedBcConfig {
                    name: name.clone(),
                    node_group: node_group.clone(),
                    constrain_x: *x,
                    constrain_y: *y,
                    constrain_z: *z,
                })
            }
            BoundaryConditionSer::Load { name, node_group, force_x, force_y, force_z } => {
                BoundaryConditionConfig::Load(LoadBcConfig {
                    name: name.clone(),
                    node_group: node_group.clone(),
                    force_x: *force_x,
                    force_y: *force_y,
                    force_z: *force_z,
                })
            }
            BoundaryConditionSer::Torque { name, node_group, axis_point, axis_direction, magnitude } => {
                BoundaryConditionConfig::Torque(TorqueBcConfig {
                    name: name.clone(),
                    node_group: node_group.clone(),
                    axis_point: *axis_point,
                    axis_direction: *axis_direction,
                    magnitude: *magnitude,
                })
            }
            BoundaryConditionSer::Contact { name, primary_surface, secondary_surface } => {
                BoundaryConditionConfig::Contact(ContactBcConfig {
                    name: name.clone(),
                    primary_surface: primary_surface.clone(),
                    secondary_surface: secondary_surface.clone(),
                })
            }
            BoundaryConditionSer::Pressure { name, element_group, pressure } => {
                BoundaryConditionConfig::Pressure(PressureBcConfig {
                    name: name.clone(),
                    element_group: element_group.clone(),
                    pressure: *pressure,
                })
            }
            BoundaryConditionSer::Traction { name, element_group, traction_type } => {
                BoundaryConditionConfig::Traction(TractionBcConfig {
                    name: name.clone(),
                    element_group: element_group.clone(),
                    traction_type: match traction_type {
                        TractionTypeSer::Uniform { fx, fy, fz } =>
                            TractionTypeConfig::Uniform { fx: *fx, fy: *fy, fz: *fz },
                        TractionTypeSer::Normal { magnitude } =>
                            TractionTypeConfig::Normal { magnitude: *magnitude },
                        TractionTypeSer::Shear { magnitude } =>
                            TractionTypeConfig::Shear { magnitude: *magnitude },
                    },
                })
            }
            BoundaryConditionSer::BodyForce { name, element_group, force_type } => {
                BoundaryConditionConfig::BodyForce(BodyForceBcConfig {
                    name: name.clone(),
                    element_group: element_group.clone(),
                    force_type: match force_type {
                        BodyForceTypeSer::Gravity { gx, gy, gz } =>
                            BodyForceTypeConfig::Gravity { gx: *gx, gy: *gy, gz: *gz },
                        BodyForceTypeSer::Centrifugal { axis_point, axis_direction, angular_velocity } =>
                            BodyForceTypeConfig::Centrifugal {
                                axis_point: *axis_point,
                                axis_direction: *axis_direction,
                                angular_velocity: *angular_velocity,
                            },
                        BodyForceTypeSer::Uniform { fx, fy, fz } =>
                            BodyForceTypeConfig::Uniform { fx: *fx, fy: *fy, fz: *fz },
                    },
                })
            }
        }
    }
}

// ============================================================================
// File I/O
// ============================================================================

/// Save project to TOML file
pub fn save_project(
    path: &Path, 
    name: &str,
    mesh_path: Option<&Path>,
    config: &SimulationConfig
) -> Result<(), String> {
    let toml_str = serialize_project_toml(name, mesh_path, config)?;
    
    std::fs::write(path, toml_str)
        .map_err(|e| format!("Failed to write file: {}", e))?;
    
    Ok(())
}

/// Serialize project to TOML string (for WASM downloads)
pub fn serialize_project_toml(
    name: &str,
    mesh_path: Option<&Path>,
    config: &SimulationConfig
) -> Result<String, String> {
    let project = ProjectFile {
        name: name.to_string(),
        description: String::new(),
        mesh: mesh_path.map(|p| p.to_string_lossy().to_string()),
        simulation: config.into(),
    };
    
    toml::to_string_pretty(&project)
        .map_err(|e| format!("Failed to serialize project: {}", e))
}

/// Load project from TOML file
pub fn load_project(path: &Path) -> Result<(ProjectFile, SimulationConfig), String> {
    let contents = std::fs::read_to_string(path)
        .map_err(|e| format!("Failed to read file: {}", e))?;
    
    parse_project_toml(&contents)
}

/// Parse project from TOML string (for WASM file uploads)
pub fn parse_project_toml(contents: &str) -> Result<(ProjectFile, SimulationConfig), String> {
    let project: ProjectFile = toml::from_str(contents)
        .map_err(|e| format!("Failed to parse TOML: {}", e))?;
    
    let config: SimulationConfig = (&project.simulation).into();
    
    Ok((project, config))
}

/// Generate example project TOML string
pub fn generate_example_toml() -> String {
    r#"# RustFEA GUI Project File
name = "My Simulation"
description = "Example simulation project"

# Mesh file (relative path or absolute)
mesh = "meshes/my_mesh.inp"

[simulation]
solver = "direct"  # or "explicit"
dofs = 3

# Materials
[[simulation.materials]]
id = 1
name = "Steel"
youngs_modulus = 200e9
poissons_ratio = 0.3
density = 7850.0

# Boundary Conditions
[[simulation.boundary_conditions]]
type = "fixed"
name = "Fixed Support"
node_group = "fixed_nodes"
x = 0.0
y = 0.0
z = 0.0

[[simulation.boundary_conditions]]
type = "load"
name = "Applied Force"
node_group = "load_nodes"
force_x = 0.0
force_y = -10000.0
force_z = 0.0

# Explicit solver settings (only needed for explicit solver)
[simulation.explicit]
time_steps = 10000
vtk_save_steps = 100
state_save_steps = 100
"#.to_string()
}
