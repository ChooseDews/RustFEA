//! RustFEA GUI - Main entry point
//!
//! A graphical user interface for creating/importing meshes, setting up simulations,
//! running FEA analysis, and visualizing results.

mod app;
mod examples;
mod project_io;
mod render_cache;
mod renderer;
mod section_cut;
mod ui;
mod state;
mod web_file_io;

use app::FeaApp;
#[cfg(not(target_arch = "wasm32"))]
use state::TestMode;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsCast;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result<()> {
    env_logger::init();
    
    // Parse command-line arguments for test mode
    let args: Vec<String> = std::env::args().collect();
    let test_mode = if args.iter().any(|a| a == "--auto-test") {
        let screenshot_path = args.iter()
            .position(|a| a == "--screenshot")
            .and_then(|i| args.get(i + 1))
            .map(|s| s.to_string())
            .unwrap_or_else(|| "test_screenshot.png".to_string());
        
        TestMode {
            enabled: true,
            screenshot_path: Some(screenshot_path),
        }
    } else {
        TestMode::default()
    };
    
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_min_inner_size([800.0, 600.0])
            .with_title("RustFEA - Finite Element Analysis"),
        ..Default::default()
    };
    
    eframe::run_native(
        "RustFEA",
        native_options,
        Box::new(move |cc| Ok(Box::new(FeaApp::new_with_test_mode(cc, test_mode.clone())))),
    )
}

#[cfg(target_arch = "wasm32")]
fn main() {
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
