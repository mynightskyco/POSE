use serde::{Deserialize, Serialize};

pub trait Simobj {
    fn type_of(&self) -> String;
    fn get_id(&self) -> u32;
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
}

impl Simobj for Debris {

    fn type_of(&self) -> String {
        return String::from("Debris");
    }

    fn get_id(&self) -> u32 {
        return self.id;
    }
}

pub type SimobjT = Box<dyn Simobj>;