use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use log::{info, debug, warn};

pub fn write_hashmap_sparse_matrix(filename: &str, sparse_matrix: &HashMap<(usize, usize), f64>) -> std::io::Result<()> {
    let nzn = sparse_matrix.len();
    info!("Writing Matrix Out - Number of Non Zero Values: {}", nzn);
    
    // Check size of matrix and reject if too large
    // This code is commented out, but you can bring it back if you need to.
    // if nzn > 100_000_000 {
    //     let error = format!("Matrix too large to write to file. Size: {}", sparse_matrix.len());
    //     warn!("Matrix too large to write to file. Size: {}", sparse_matrix.len());
    //     return Err(std::io::Error::new(std::io::ErrorKind::Other, error));
    // }

    debug!("Opening file for writing: {}", filename);
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    // Collect the filtered pairs into a Vec and sort
    let mut entries: Vec<_> = sparse_matrix.iter().filter(|&(&(row, col), _)| row <= col).collect();
    entries.sort_by_key(|&((row, col), _)| (row, col));

    for &((row, col), &value) in &entries {
        writeln!(writer, "{} {} {}", row, col, value)?;
    }

    debug!("Finished writing matrix to file");
    Ok(())
}

pub fn write_vector(filename: &str, vector: &Vec<f64>) -> std::io::Result<()> {
    debug!("Writing vector to file: {}", filename);
    // Open the file for writing and buffer the writes
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    for (i, value) in vector.iter().enumerate() {
        writeln!(writer, "{} {}", i, value)?;
    }

    // Ensure everything is written to the file
    writer.flush()?;

    debug!("Finished writing vector to file");
    Ok(())
}