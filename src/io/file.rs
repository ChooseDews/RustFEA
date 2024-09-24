// read json config file
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use serde::Serialize;
use bincode;

fn get_writer(filename: &str) -> Box<dyn Write> {
    let file_extension = Path::new(filename).extension().expect("Issue parsing file extension").to_str().unwrap();
    let file = File::create(filename).unwrap();
    match file_extension {
        "json" | "bin" => Box::new(BufWriter::new(file)),
        "xz" => Box::new(xz2::write::XzEncoder::new(file, 6)),
        _ => panic!("Unsupported file extension: {}", file_extension),
    }
}

pub fn read_json_file(filename: &str) -> Value {
    let file_extension = Path::new(filename).extension().expect("Issue parsing file extension").to_str().unwrap();
    let file = File::open(filename).unwrap();
    let reader: Box<dyn Read> = if file_extension == "json" {
        Box::new(BufReader::new(file))
    } else if file_extension == "xz" {
        Box::new(xz2::read::XzDecoder::new(file))
    } else {
        panic!("Unsupported file extension: {}", file_extension);
    };
    serde_json::from_reader(reader).unwrap()
}

pub fn write_json_file(filename: &str, data: &Value) {
    let writer = get_writer(filename);
    serde_json::to_writer(writer, &data).expect("Failed to write JSON data");
    println!("Wrote JSON file to {:?}", filename);
}

pub fn seralized_write_json<T: Serialize>(filename: &str, data: &T) {
    let value = serde_json::to_value(data).expect("Failed to serialize data to JSON");
    write_json_file(filename, &value);
}

pub fn seralized_write_bincode<T: Serialize>(filename: &str, data: &T) { //also handle compression
    let encoded = bincode::serialize(data).expect("Failed to serialize data to Bincode");
    let mut writer = get_writer(filename);
    writer.write_all(&encoded).expect("Failed to write Bincode data");
    println!("Wrote Bincode file to {:?}", filename);
}





