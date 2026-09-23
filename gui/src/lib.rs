//! RustFEA GUI Library
//!
//! This library provides a graphical user interface for creating/importing meshes,
//! setting up simulations, running FEA analysis, and visualizing results.
//!
//! Built with egui for the UI and targeting both native (desktop) and web platforms.

pub mod app;
pub mod examples;
pub mod icons;
pub mod project_io;
pub mod render_cache;
pub mod renderer;
pub mod section_cut;
pub mod state;
pub mod ui;
pub mod web_file_io;

pub use app::FeaApp;
pub use state::AppState;

/// WASM entry point - called automatically when the module loads
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen::prelude::wasm_bindgen(start)]
pub fn wasm_main() {
    use wasm_bindgen::JsCast;
    
    // Redirect panics to console.error
    console_error_panic_hook::set_once();
    
    // Initialize logging
    console_log::init_with_level(log::Level::Debug).expect("Failed to initialize logger");
    
    let web_options = eframe::WebOptions::default();
    
    wasm_bindgen_futures::spawn_local(async {
        // Get the canvas element from the DOM
        let document = web_sys::window()
            .expect("No window")
            .document()
            .expect("No document");
        let canvas = document
            .get_element_by_id("rustfea_canvas")
            .expect("Failed to find canvas element")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("Element is not a canvas");
        
        eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(|cc| Ok(Box::new(FeaApp::new(cc)))),
            )
            .await
            .expect("Failed to start eframe");
    });
}
