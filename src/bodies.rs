
pub trait Simobj {}

/// Struct for holding attributes relating to debris
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

impl Simobj for Debris {}

pub type SimobjT = Box<dyn Simobj>;

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