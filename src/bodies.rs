
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
}

impl Simobj for Debris {}

pub type SimobjT = Box<dyn Simobj>;