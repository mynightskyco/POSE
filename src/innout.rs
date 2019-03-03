
use super::bodies;

/// Main entry point into the init sequence
/// 
/// # Argument
/// * 'file' - The name of the input file containing the bodies
pub fn parse_inpt(file: &str) -> Vec<bodies::SimobjT>{
    let sim_bodies: Vec<bodies::SimobjT> = Vec::new();

    // REMOVE THIS
    println!("{file}", file=file);
    let _test = bodies::Debris{id : 1, 
                              x_dis : 1f64, 
                              y_dis : 1f64,
                              z_dis : 1f64,
                              x_vel : 1f64,
                              y_vel : 1f64,
                              z_vel : 1f64};

    return sim_bodies;
}