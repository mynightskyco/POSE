use serde::{Deserialize, Serialize};

pub type SimobjT = Box<dyn Simobj>;

#[derive(Serialize, Deserialize, Clone)]
pub struct Objects {
    pub debris: Vec<Debris>,
    pub spacecraft: Vec<Spacecraft>
}

pub trait Simobj {
    fn type_of(&self) -> String;
    fn get_id(&self) -> u32;
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
    fn get_id(&self) -> u32 { return self.id; }
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
    fn get_id(&self) -> u32 {
        return self.id;
    }
}

enum Solarobj{
    Sun{coords: CartesianCoords, attr: SolarAttr},
    Earth{coords: CartesianCoords, attr: SolarAttr},
    Moon{coords: CartesianCoords, attr: SolarAttr}
}

struct SolarAttr{
    radius: f64, // meters
    mass: f64 // kg
}

struct PlanetPL { // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    solartype: Solarobj, // Type enum of the solar obj
    n0: f32, nc: f32, // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
    i0: f32, ic: f32, // inclination to the ecliptic (deg)
    w0: f32, wc: f32, // argument of perihelion (deg)
    a0: f32, ac: f32, // semi-major axis, or mean distance from Sun (AU)
    e0: f32, ec: f32, // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
    m0: f32, mc: f32  // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change
}

pub struct CartesianCoords {
    xh: f64,
    yh: f64,
    zh: f64
}

mod kepler_utilities {
    use std::f32::{self, consts};
    use super::CartesianCoords;
    
    pub fn eccentric_anomaly(e: f32, m: f32) -> f32 {
        // TODO Create macro for sin and cos degrees
        let mut ecc: f32 = m + (e * sin_deg!(m) * (1f32 + (e * cos_deg!(m))));

        loop {
            let f: f32 = ecc - (ecc - (deg_from_rad!() * e * sin_deg!(ecc)) - m) / (1f32 - e * cos_deg!(ecc));
            let error = (f - ecc).abs();
            ecc = f;

            if error < 1.0e-8 { break; }
        };

        ecc
    }

    pub fn lunar_pertub(xh: f64, yh: f64, zh: f64, day: f64) -> CartesianCoords{

        // TODO add petrub code
        CartesianCoords{xh, yh, zh}
    }
}

trait KeplerModel{

    fn update_ecliptic_cartesian_coords(&self, day: f32);

    fn perturb(&self, xh: f64, yh: f64, zh: f64, day: f64) -> CartesianCoords{

        CartesianCoords{xh, yh, zh}
    }
}

impl KeplerModel for PlanetPL{

    fn update_ecliptic_cartesian_coords(&self, day: f32) {
        // Default impl
        let a = self.a0 + (day * self.ac);
        let e = self.e0 + (day * self.ec);
        let m_u = self.m0 + (day * self.mc);
        let n_u = self.n0 + (day * self.nc);
        let w = self.w0 + (day * self.wc);
        let i = self.i0 + (day * self.ic);
        let ecc = kepler_utilities::eccentric_anomaly(e, m_u);

        // TODO finsh this function

    }

    fn perturb(&self, xh: f64, yh: f64, zh: f64, day: f64) -> CartesianCoords {
        match &self.solartype {
            Solarobj::Moon{coords, attr} => CartesianCoords{xh, yh, zh},
            _ => CartesianCoords{xh, yh, zh}
        }
    }
}
