use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;

pub fn write_hashmap_sparse_matrix(filename: &str, sparse_matrix: &HashMap<(usize, usize), f64>) -> std::io::Result<()> {
    let nzn = sparse_matrix.len();

    println!("Number of Non Zero Values: {}", nzn);
    // Check size of matrix and reject if too large
    // if nzn > 100_000_000 {
    //     let error = format!("Matrix too large to write to file. Size: {}", sparse_matrix.len());
    //     return Err(std::io::Error::new(std::io::ErrorKind::Other, error));
    // }

    // Open the file for writing directly
    let mut file = File::create(filename)?;

    // Sort hashmap keys and write to file directly
    let mut keys: Vec<(usize, usize)> = sparse_matrix.keys().cloned().collect();
    keys.sort();

    for key in keys {
        let value = sparse_matrix[&key];
        let (row, col) = key;

        // Only write upper right triangle
        if row <= col {
            writeln!(file, "{}, {}, {}", row, col, value)?;
        }
    }

    Ok(())
}