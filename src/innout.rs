
extern crate serde_json;

use super::bodies;

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/**
 * Main entry point into the init sequence
 *
 * ### Argument
 * * 'file' - The name of the input file containing the bodies
 *
 * ### Return
 *      A vector of bodies from the input file.
 */
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

    assign_id(&mut sim_bodies);

    return sim_bodies;
}

/**
 * Function responsable for handling opeing the file and conneting the
 * serde reader.
 *
 * ### Argument
 * * 'path' - The path to the input bodies json file.
 *
 * ### Return
 *      A result object loaded with an IO error on failure or the serde reader on
 *      success.
 */
fn read_object_from_file<P: AsRef<Path>>(path: P) -> Result<bodies::Objects, Box<Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Objects`.
    let u = serde_json::from_reader(reader)?;

    // Return the `Objects`.
    Ok(u)
}

/**
 * Adds an sequential id value to each of the simulation bodies.
 *
 * ### Argument
 * * 'sim_bodies' - A vector containing both debris and spacecraft objects.
 */
fn assign_id(sim_bodies: &mut Vec<bodies::SimobjT>) {
    let mut id_inc: u32 = 1;

    for body in sim_bodies {
        *body.id_mut() = id_inc;
        id_inc += 1;
    }
}