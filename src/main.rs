
extern crate serde_json;
extern crate clap;

mod input;
mod bodies;

fn main() {

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
                .takes_value(true)
        ])
        .get_matches();

    let inpt_file = matches.value_of("INPUT").unwrap(); // Will always have INPUT
    input::load_inpt(inpt_file);
}
