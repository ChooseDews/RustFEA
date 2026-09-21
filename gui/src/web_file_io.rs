//! Web (WASM) file I/O operations using browser File API
//!
//! Provides upload (file picker) and download functionality for the browser environment.

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsCast;
#[cfg(target_arch = "wasm32")]
use web_sys::{Document, HtmlInputElement, HtmlAnchorElement, Blob, BlobPropertyBag, Url, FileReader};
#[cfg(target_arch = "wasm32")]
use std::cell::RefCell;
#[cfg(target_arch = "wasm32")]
use std::rc::Rc;

/// Pending file data from an upload operation
#[derive(Clone)]
pub struct PendingFile {
    pub name: String,
    pub data: Vec<u8>,
    pub file_type: FileType,
}

/// Type of file being loaded
#[derive(Clone, Copy, PartialEq)]
pub enum FileType {
    Mesh,       // .inp, .bin, .json mesh files
    Project,    // .toml project files
    Bundle,     // .rfea ZIP bundles (project + mesh)
}

/// Global state for pending file uploads (WASM only)
/// Uses thread-local storage since WASM is single-threaded
#[cfg(target_arch = "wasm32")]
thread_local! {
    static PENDING_FILE: RefCell<Option<PendingFile>> = RefCell::new(None);
}

/// Check if there's a pending file upload and take it
#[cfg(target_arch = "wasm32")]
pub fn take_pending_file() -> Option<PendingFile> {
    PENDING_FILE.with(|cell| cell.borrow_mut().take())
}

/// Check if there's a pending file (without taking it)
#[cfg(target_arch = "wasm32")]
pub fn has_pending_file() -> bool {
    PENDING_FILE.with(|cell| cell.borrow().is_some())
}

/// Open a file picker dialog for mesh files (.inp, .bin, .json)
#[cfg(target_arch = "wasm32")]
pub fn open_mesh_file_picker() {
    open_file_picker(".inp,.bin,.json,.xz", FileType::Mesh);
}

/// Open a file picker dialog for project files (.toml)
#[cfg(target_arch = "wasm32")]
pub fn open_project_file_picker() {
    open_file_picker(".toml", FileType::Project);
}

/// Internal: Create and trigger a file input element
#[cfg(target_arch = "wasm32")]
fn open_file_picker(accept: &str, file_type: FileType) {
    let window = match web_sys::window() {
        Some(w) => w,
        None => {
            log::error!("No window object");
            return;
        }
    };
    
    let document = match window.document() {
        Some(d) => d,
        None => {
            log::error!("No document object");
            return;
        }
    };
    
    // Create a hidden file input
    let input: HtmlInputElement = match document.create_element("input") {
        Ok(el) => match el.dyn_into::<HtmlInputElement>() {
            Ok(input) => input,
            Err(_) => {
                log::error!("Failed to cast to HtmlInputElement");
                return;
            }
        },
        Err(_) => {
            log::error!("Failed to create input element");
            return;
        }
    };
    
    input.set_type("file");
    input.set_accept(accept);
    input.style().set_property("display", "none").ok();
    
    // Append to body temporarily
    if let Some(body) = document.body() {
        body.append_child(&input).ok();
    }
    
    // Set up the change handler
    let input_clone = input.clone();
    let closure = Closure::wrap(Box::new(move |_event: web_sys::Event| {
        if let Some(files) = input_clone.files() {
            if files.length() > 0 {
                if let Some(file) = files.get(0) {
                    read_file(file, file_type);
                }
            }
        }
        // Clean up the input element
        if let Some(parent) = input_clone.parent_node() {
            parent.remove_child(&input_clone).ok();
        }
    }) as Box<dyn FnMut(_)>);
    
    input.set_onchange(Some(closure.as_ref().unchecked_ref()));
    closure.forget(); // Let JS garbage collector handle it
    
    // Trigger the file picker
    input.click();
}

/// Read a File object and store it in PENDING_FILE
#[cfg(target_arch = "wasm32")]
fn read_file(file: web_sys::File, file_type: FileType) {
    let file_name = file.name();
    
    let reader = match FileReader::new() {
        Ok(r) => r,
        Err(_) => {
            log::error!("Failed to create FileReader");
            return;
        }
    };
    
    let reader_clone = reader.clone();
    let file_name_clone = file_name.clone();
    
    let onload = Closure::wrap(Box::new(move |_event: web_sys::Event| {
        if let Ok(result) = reader_clone.result() {
            if let Some(array_buffer) = result.dyn_ref::<js_sys::ArrayBuffer>() {
                let uint8_array = js_sys::Uint8Array::new(array_buffer);
                let data = uint8_array.to_vec();
                
                let pending = PendingFile {
                    name: file_name_clone.clone(),
                    data,
                    file_type,
                };
                
                PENDING_FILE.with(|cell| {
                    *cell.borrow_mut() = Some(pending);
                });
                
                log::info!("File loaded: {} ({} bytes)", file_name_clone, uint8_array.length());
            }
        }
    }) as Box<dyn FnMut(_)>);
    
    reader.set_onload(Some(onload.as_ref().unchecked_ref()));
    onload.forget();
    
    // Start reading
    if let Err(e) = reader.read_as_array_buffer(&file) {
        log::error!("Failed to read file: {:?}", e);
    }
}

/// Download data as a file in the browser
#[cfg(target_arch = "wasm32")]
pub fn download_file(filename: &str, data: &[u8], mime_type: &str) {
    let window = match web_sys::window() {
        Some(w) => w,
        None => {
            log::error!("No window object");
            return;
        }
    };
    
    let document = match window.document() {
        Some(d) => d,
        None => {
            log::error!("No document object");
            return;
        }
    };
    
    // Create a Blob from the data
    let uint8_array = js_sys::Uint8Array::new_with_length(data.len() as u32);
    uint8_array.copy_from(data);
    
    let array = js_sys::Array::new();
    array.push(&uint8_array.buffer());
    
    let mut options = BlobPropertyBag::new();
    options.type_(mime_type);
    
    let blob = match Blob::new_with_u8_array_sequence_and_options(&array, &options) {
        Ok(b) => b,
        Err(_) => {
            log::error!("Failed to create Blob");
            return;
        }
    };
    
    // Create object URL
    let url = match Url::create_object_url_with_blob(&blob) {
        Ok(u) => u,
        Err(_) => {
            log::error!("Failed to create object URL");
            return;
        }
    };
    
    // Create anchor and trigger download
    let anchor: HtmlAnchorElement = match document.create_element("a") {
        Ok(el) => match el.dyn_into::<HtmlAnchorElement>() {
            Ok(a) => a,
            Err(_) => {
                log::error!("Failed to cast to HtmlAnchorElement");
                Url::revoke_object_url(&url).ok();
                return;
            }
        },
        Err(_) => {
            log::error!("Failed to create anchor element");
            Url::revoke_object_url(&url).ok();
            return;
        }
    };
    
    anchor.set_href(&url);
    anchor.set_download(filename);
    anchor.style().set_property("display", "none").ok();
    
    if let Some(body) = document.body() {
        body.append_child(&anchor).ok();
        anchor.click();
        body.remove_child(&anchor).ok();
    }
    
    // Clean up the URL (slight delay to ensure download starts)
    let url_clone = url.clone();
    let cleanup = Closure::wrap(Box::new(move || {
        Url::revoke_object_url(&url_clone).ok();
    }) as Box<dyn FnMut()>);
    
    window.set_timeout_with_callback_and_timeout_and_arguments_0(
        cleanup.as_ref().unchecked_ref(),
        100,
    ).ok();
    cleanup.forget();
}

/// Download a project as TOML
#[cfg(target_arch = "wasm32")]
pub fn download_project(filename: &str, toml_content: &str) {
    download_file(filename, toml_content.as_bytes(), "application/toml");
}

/// Download mesh as JSON
#[cfg(target_arch = "wasm32")]
pub fn download_mesh_json(filename: &str, json_content: &str) {
    download_file(filename, json_content.as_bytes(), "application/json");
}

/// Download results as VTK
#[cfg(target_arch = "wasm32")]
pub fn download_vtk(filename: &str, vtk_content: &str) {
    download_file(filename, vtk_content.as_bytes(), "text/plain");
}

/// Create and download a project bundle as a ZIP file
/// Contains: project.toml + mesh.json (serialized mesh)
#[cfg(target_arch = "wasm32")]
pub fn download_project_bundle(project_name: &str, toml_content: &str, mesh_json: Option<&str>) {
    use std::io::Write;
    use zip::write::SimpleFileOptions;
    use zip::ZipWriter;
    
    let mut buffer = Vec::new();
    
    {
        let mut zip = ZipWriter::new(std::io::Cursor::new(&mut buffer));
        let options = SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated);
        
        // Add project.toml
        if zip.start_file("project.toml", options).is_ok() {
            let _ = zip.write_all(toml_content.as_bytes());
        }
        
        // Add mesh.json if present
        if let Some(mesh) = mesh_json {
            if zip.start_file("mesh.json", options).is_ok() {
                let _ = zip.write_all(mesh.as_bytes());
            }
        }
        
        let _ = zip.finish();
    }
    
    let filename = format!("{}.rfea", project_name.replace(" ", "_"));
    download_file(&filename, &buffer, "application/zip");
}

/// Extracted project bundle contents
pub struct ProjectBundle {
    pub project_toml: Option<String>,
    pub mesh_json: Option<String>,
}

/// Extract a project bundle from ZIP data
pub fn extract_project_bundle(data: &[u8]) -> Result<ProjectBundle, String> {
    use std::io::Read;
    use zip::ZipArchive;
    
    let cursor = std::io::Cursor::new(data);
    let mut archive = ZipArchive::new(cursor)
        .map_err(|e| format!("Failed to read ZIP archive: {}", e))?;
    
    let mut bundle = ProjectBundle {
        project_toml: None,
        mesh_json: None,
    };
    
    for i in 0..archive.len() {
        let mut file = archive.by_index(i)
            .map_err(|e| format!("Failed to read ZIP entry: {}", e))?;
        
        let name = file.name().to_lowercase();
        
        if name == "project.toml" {
            let mut content = String::new();
            file.read_to_string(&mut content)
                .map_err(|e| format!("Failed to read project.toml: {}", e))?;
            bundle.project_toml = Some(content);
        } else if name == "mesh.json" {
            let mut content = String::new();
            file.read_to_string(&mut content)
                .map_err(|e| format!("Failed to read mesh.json: {}", e))?;
            bundle.mesh_json = Some(content);
        }
    }
    
    Ok(bundle)
}

/// Open a file picker for project bundles (.rfea ZIP files) and TOML projects
#[cfg(target_arch = "wasm32")]
pub fn open_bundle_file_picker() {
    let window = match web_sys::window() {
        Some(w) => w,
        None => return,
    };
    
    let document = match window.document() {
        Some(d) => d,
        None => return,
    };
    
    let input: HtmlInputElement = match document.create_element("input") {
        Ok(el) => match el.dyn_into::<HtmlInputElement>() {
            Ok(input) => input,
            Err(_) => return,
        },
        Err(_) => return,
    };
    
    input.set_type("file");
    input.set_accept(".rfea,.zip,.toml");
    input.style().set_property("display", "none").ok();
    
    let callback = {
        let input = input.clone();
        Closure::wrap(Box::new(move |_: web_sys::Event| {
            if let Some(files) = input.files() {
                if let Some(file) = files.get(0) {
                    let reader = FileReader::new().unwrap();
                    let reader_clone = reader.clone();
                    let file_name = file.name();
                    let is_toml = file_name.to_lowercase().ends_with(".toml");
                    
                    let onload = Closure::wrap(Box::new(move |_: web_sys::Event| {
                        if let Ok(result) = reader_clone.result() {
                            if let Some(array_buffer) = result.dyn_ref::<js_sys::ArrayBuffer>() {
                                let uint8_array = js_sys::Uint8Array::new(array_buffer);
                                let data: Vec<u8> = uint8_array.to_vec();
                                
                                // Determine file type based on extension
                                let file_type = if is_toml {
                                    FileType::Project
                                } else {
                                    FileType::Bundle
                                };
                                
                                PENDING_FILE.with(|cell| {
                                    *cell.borrow_mut() = Some(PendingFile {
                                        name: file_name.clone(),
                                        data,
                                        file_type,
                                    });
                                });
                            }
                        }
                    }) as Box<dyn FnMut(_)>);
                    
                    reader.set_onload(Some(onload.as_ref().unchecked_ref()));
                    onload.forget();
                    reader.read_as_array_buffer(&file).ok();
                }
            }
        }) as Box<dyn FnMut(_)>)
    };
    
    input.add_event_listener_with_callback("change", callback.as_ref().unchecked_ref()).ok();
    callback.forget();
    
    if let Some(body) = document.body() {
        body.append_child(&input).ok();
        input.click();
        body.remove_child(&input).ok();
    }
}

// Native stubs - these operations use rfd on native
#[cfg(not(target_arch = "wasm32"))]
pub fn take_pending_file() -> Option<PendingFile> { None }

#[cfg(not(target_arch = "wasm32"))]
pub fn has_pending_file() -> bool { false }

#[cfg(not(target_arch = "wasm32"))]
pub fn open_mesh_file_picker() {}

#[cfg(not(target_arch = "wasm32"))]
pub fn open_project_file_picker() {}

#[cfg(not(target_arch = "wasm32"))]
pub fn open_bundle_file_picker() {}

#[cfg(not(target_arch = "wasm32"))]
pub fn download_file(_filename: &str, _data: &[u8], _mime_type: &str) {}

#[cfg(not(target_arch = "wasm32"))]
pub fn download_project(_filename: &str, _toml_content: &str) {}

#[cfg(not(target_arch = "wasm32"))]
pub fn download_mesh_json(_filename: &str, _json_content: &str) {}

#[cfg(not(target_arch = "wasm32"))]
pub fn download_vtk(_filename: &str, _vtk_content: &str) {}

#[cfg(not(target_arch = "wasm32"))]
pub fn download_project_bundle(_project_name: &str, _toml_content: &str, _mesh_json: Option<&str>) {}
