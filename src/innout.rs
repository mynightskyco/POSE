extern crate serde_json;

use super::bodies;
use crate::bodies::{SimobjT, Debris, Simobj};

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
/// Main entry point into the init sequence
/// 
/// # Argument
/// * 'file' - The name of the input file containing the bodies
pub fn parse_inpt(file: &str) -> Vec<bodies::SimobjT>{
    let mut sim_bodies: Vec<bodies::SimobjT> = Vec::new();

    // REMOVE THIS
    //println!("{file}", file=file);
//    let mut _test1 = bodies::Debris{id : 1,
//                              x_dis : 1f64,
//                              y_dis : 1f64,
//                              z_dis : 1f64,
//                              x_vel : 1f64,
//                              y_vel : 1f64,
//                              z_vel : 1f64};


    let mut _test: Vec<Debris> = read_object_from_file("test.json").unwrap();

    for elem in _test.iter() {
        //println!("{:?}", elem.id);
        let p = Box::new(elem.clone());
        sim_bodies.push(p);

    }

    for e in sim_bodies.iter() {
        println!("what is this {}", e.type_of());
    }

    //Testing for seeing what the output of the Debris looks like in JSON

//    let mut vec: Vec<Debris> = Vec::new();
//    vec.push(_test);
//    vec.push(_test2);
//    let str = serde_json::to_string_pretty(&vec);
//    println!("{}", str.unwrap());

    return sim_bodies;
}

fn read_object_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<Debris>, Box<Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Read the JSON contents of the file as an instance of `Debris`.
    let u = serde_json::from_reader(reader)?;

    // Return the `Debris`.
    Ok(u)
}
