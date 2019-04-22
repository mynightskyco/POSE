//!
//! POSE - Parallel Orbital Simulation Environment
//! TODO - Add more doc

#[macro_use]
mod macros;

extern crate clap;

mod innout;
mod bodies;

mod cli{


    ///Checks if value passed in to program argument is numeric. Returns a Result
    ///
    ///# Argument
    ///* 'strng' - The value passed by the user
    fn numeric_validator(strng: String) -> Result<(), String>{
        if strng.parse::<f64>().is_ok(){
            Ok(())
        } else {
            Err(String::from("Input is non-numeric"))
        }
    }

    /// Defines the argument structure for the pose simulation program
    /// Returns the result of user arguments passed over the cli
    pub fn check_cli() -> clap::ArgMatches<'static> {

        // Defines the input arguments from the cli
        let matches = clap::App::new("Parallel Orbital Simulation Environment (POSE)")
            .version("DEV0.1")
            .about("Simulation aimed to model the orbital environment around Earth for bodies at all magnitudes.")
            .args(&[
                clap::Arg::with_name("INPUT")
                    .help("json file containing information on bodies at initilzation.")
                    .required(true)
                    .index(1),
                clap::Arg::with_name("out")
                    .help("Main name of output files.")
                    .short("o")
                    .long("out")
                    .value_name("FILE_NAME")
                    .takes_value(true),
                clap::Arg::with_name("timei")
                    .help("Initial time for simulation start must, be in iso time format.")
                    .long("timei")
                    .value_name("INITAL_TIME")
                    .takes_value(true),
                clap::Arg::with_name("step")
                    .help("Simulation time step interval in seconds")
                    .short("s")
                    .long("step")
                    .value_name("STEP_INTERVAL")
                    .takes_value(true)
                    .default_value("0.1") // Simulation step time
                    .validator(numeric_validator)
            ])
            .get_matches();

        return matches;
    }
}

fn main() {
    let matches = cli::check_cli();
    let inpt_file = matches.value_of("INPUT").unwrap(); // Will always have INPUT
    let sim_bodies = innout::parse_inpt(inpt_file);

    // For testing --------------
    let day: f32 = 139f32; // Sat Jan 01 2000 01:01:01 GMT-0500 (EST)

    let mut planet_bodies = bodies::solar_system_objs(day);

    for e in sim_bodies.iter() {
        println!("{} with id {}", e.type_of(), e.get_id());
    }

    for i in 1..366 {
        bodies::update_solar_system_objs(&mut planet_bodies, i as f32);
        planet_bodies.get_mut(0).unwrap().mut_coords().to_meters();
        planet_bodies.get_mut(1).unwrap().mut_coords().to_meters();
        planet_bodies.get_mut(2).unwrap().mut_coords().to_meters();
        print!("Sun,{0},{1},{2},{3}\n",
               i,
               planet_bodies.get(0).unwrap().get_coords().xh,
               planet_bodies.get(0).unwrap().get_coords().yh,
               planet_bodies.get(0).unwrap().get_coords().zh,
        );
        print!("Earth,{0},{1},{2},{3}\n",
               i,
               planet_bodies.get(1).unwrap().get_coords().xh,
               planet_bodies.get(1).unwrap().get_coords().yh,
               planet_bodies.get(1).unwrap().get_coords().zh,
        );
        print!("Moon,{0},{1},{2},{3}\n",
               i,
               planet_bodies.get(2).unwrap().get_coords().xh,
               planet_bodies.get(2).unwrap().get_coords().yh,
               planet_bodies.get(2).unwrap().get_coords().zh,
        ); // REMOVE this
    }

    // ---------------------------

}
