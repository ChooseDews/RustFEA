use std::fs::File;
use std::io::prelude::*;
use std::collections::HashMap;


pub fn write_hashmap_sparse_matrix(filename: &str, sparse_matrix: &HashMap<(usize, usize), f64>) -> std::io::Result<()> {

    //check size of matrix reject if too large
    let hashmap_length = sparse_matrix.len();
    if hashmap_length > 10_000_000 {
        let error = format!("Matrix too large to write to file. Size: {}", hashmap_length);
        return Err(std::io::Error::new(std::io::ErrorKind::Other, error));
    }
    
    //sort by row, col then write to file row, col, value
    let mut file = File::create(filename)?;
    //sort hashmap by key
    let mut keys: Vec<(usize, usize)> = sparse_matrix.keys().cloned().collect();
    keys.sort();
    //write to file
    for key in keys{
        let value = sparse_matrix[&key];
        let row = key.0;
        let col = key.1;
        //only write upper right triangle
        // if row > col {
        //     continue;
        // }
        let line = format!("{}, {}, {}\n", row, col, value);
        file.write_all(line.as_bytes())?;
    }
    //close file
    Ok(())
}