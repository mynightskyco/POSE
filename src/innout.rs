
extern crate serde_json;

use super::bodies;

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Main entry point into the init sequence
/// 
/// ## Argument
/// * 'file' - The name of the input file containing the bodies
pub fn parse_inpt(file: &str) -> Vec<bodies::SimobjT>{
    let mut sim_bodies: Vec<bodies::SimobjT> = Vec::new();

    let ser_objs = read_object_from_file(file).unwrap();

    //add objects to sim_bodies
    for elem in ser_objs.debris {
        //println!("id {}", elem.type_of());
        let p = Box::new(elem.clone());
        sim_bodies.push(p);
    }
    for elem in ser_objs.spacecraft {
        //println!("id {}", elem.type_of());
        let p = Box::new(elem.clone());
        sim_bodies.push(p);
    }


    return sim_bodies;
}

fn read_object_from_file<P: AsRef<Path>>(path: P) -> Result<bodies::Objects, Box<Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Objects`.
    let u = serde_json::from_reader(reader)?;

    // Return the `Objects`.
    Ok(u)
}