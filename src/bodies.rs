use serde::{Deserialize, Serialize, Serializer};

pub type SimobjT = Box<dyn Simobj>;


#[derive(Serialize, Deserialize, Clone)]
pub struct Objects {
    pub Debris: Vec<Debris>,
    pub Spacecraft: Vec<Spacecraft>
}

pub trait Simobj {
    fn type_of(&self) -> String;
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Spacecraft{
    pub id: u32,
    pub x_dis: f64,
    pub y_dis: f64,
    pub z_dis: f64,
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64
    // Add more
}

impl Simobj for Spacecraft {
    fn type_of(&self) -> String {
        return String::from("Spacecraft");
    }
}
/// Struct for holding attributes relating to debris
#[derive(Serialize, Deserialize, Clone)]
pub struct Debris{
    pub id: u32,
    pub x_dis: f64,
    pub y_dis: f64,
    pub z_dis: f64,
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64
    // Add more
}

impl Simobj for Debris {
    fn type_of(&self) -> String {
        return String::from("Debris");
    }
}

pub enum Solarobj{
    Sun,
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune
}

pub struct LargeBody{
    pub solartype: Solarobj,
    pub mass: f64, // kg
    // meters
    pub radius: f64,
    pub x_dis: f64,
    pub y_dis: f64,
    pub z_dis: f64,
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64    
}