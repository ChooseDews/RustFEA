//! RustFEA GUI Library
//!
//! This library provides a graphical user interface for creating/importing meshes,
//! setting up simulations, running FEA analysis, and visualizing results.
//!
//! Built with egui for the UI and targeting both native (desktop) and web platforms.

pub mod app;
pub mod examples;
pub mod project_io;
pub mod render_cache;
pub mod renderer;
pub mod section_cut;
pub mod state;
pub mod ui;
pub mod web_file_io;

pub use app::FeaApp;
pub use state::AppState;
