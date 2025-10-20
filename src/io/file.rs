/// This module provides functions to read and write serialized data to and from files.
/// It supports JSON, Bincode, and optionally compressed files using the XZ format.

use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use serde::{Serialize, Deserialize};
use bincode;
use serde::de::DeserializeOwned;
use log::{debug, error};
use std::time::Instant;
/// Returns a writer for the specified file.
/// Supports `.json`, `.bin`, and `.xz` (compressed) file extensions.
///
/// # Arguments
///
/// * `filename` - A string slice that holds the name of the file.
///
/// # Panics
///
/// Panics if the file extension is unsupported or if the file cannot be created.
fn get_writer(filename: &str) -> Box<dyn Write> {
    let file_extension = Path::new(filename).extension().expect("Issue parsing file extension").to_str().unwrap();
    let file = File::create(filename).unwrap();
    match file_extension {
        "json" | "bin" => Box::new(BufWriter::new(file)),
        "xz" => {
            // Use level 1 for much faster compression (vs level 6)
            // Level 1: ~10x faster than level 6, still good compression
            Box::new(xz2::write::XzEncoder::new(file, 1))
        },
        "zst" => {
            // Zstandard compression - often faster than xz
            // Level 1: Very fast compression with good ratio
            // auto_finish() ensures internal buffer (~128kb) is flushed on drop
            Box::new(zstd::Encoder::new(file, 1).unwrap().auto_finish())
        },
        _ => panic!("Unsupported file extension: {}", file_extension),
    }
}

/// Returns a reader for the specified file.
/// Supports `.json`, `.bin`, and `.xz` (compressed) file extensions.
///
/// # Arguments
///
/// * `filename` - A string slice that holds the name of the file.
///
/// # Panics
///
/// Panics if the file extension is unsupported or if the file cannot be opened.
fn get_reader(filename: &str) -> Box<dyn Read> {
    let file_extension = Path::new(filename).extension().expect("Issue parsing file extension").to_str().unwrap();
    let file = File::open(filename).unwrap();
    match file_extension {
        "json" | "bin" => Box::new(BufReader::new(file)),
        "xz" => Box::new(xz2::read::XzDecoder::new(file)),
        "zst" => Box::new(zstd::Decoder::new(file).unwrap()),
        _ => panic!("Unsupported file extension: {}", file_extension),
    }
}

/// Reads and deserializes data from a file.
/// Supports `.json` and `.bin` file extensions.
///
/// # Arguments
///
/// * `filename` - A string slice that holds the name of the file.
///
/// # Returns
///
/// Returns the deserialized data of type `T`.
///
/// # Panics
///
/// Panics if the file extension is unsupported or if deserialization fails.
pub fn seralized_read<T: DeserializeOwned>(filename: &str) -> T {
    debug!("Reading serialized data from file: {}", filename);
    let mut reader: Box<dyn Read> = get_reader(filename);
    if filename.contains(".json") {
        serde_json::from_reader(&mut reader).expect("Failed to deserialize data from JSON")
    } else if filename.contains(".bin") {
        bincode::deserialize_from(&mut reader).expect("Failed to deserialize data from Bincode")
    } else {
        error!("Unsupported file extension for: {}", filename);
        panic!("Unsupported file extension for: {}", filename)
    }
}

/// Serializes and writes data to a file.
/// Supports `.json` and `.bin` file extensions.
///
/// # Arguments
///
/// * `filename` - A string slice that holds the name of the file.
/// * `data` - A reference to the data to be serialized.
///
/// # Panics
///
/// Panics if the file extension is unsupported or if serialization fails.
pub fn seralized_write<T: Serialize>(filename: &str, data: &T) {
    debug!("Writing serialized data to file: {}", filename);
    let start_time = Instant::now();
    if filename.contains(".json") {
        serde_json::to_writer(&mut get_writer(filename), &data).expect("Failed to serialize data to JSON")
    } else if filename.contains(".bin") {
        bincode::serialize_into(&mut get_writer(filename), &data).expect("Failed to serialize data to Bincode")
    } else {
        error!("Unsupported file extension for: {}", filename);
        panic!("Unsupported file extension for: {}", filename)
    }
    let end_time = Instant::now();
    let duration = end_time - start_time;
    debug!("Serialization time: {:.6} seconds", duration.as_secs_f64());
}