use serde::{Deserialize, Serialize};

pub type SimobjT = Box<dyn Simobj>;
pub type PlanetBody = Box<dyn KeplerModel>;

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

#[derive(Debug)]
pub enum Solarobj{
    Sun{attr: SolarAttr},
    Earth{attr: SolarAttr},
    Moon{attr: SolarAttr}
}

#[derive(Debug)]
pub struct SolarAttr{
    radius: f64, // meters
    mass: f64 // kg
}

#[derive(Debug)]
pub struct PlanetPL { // See  http://www.stjarnhimlen.se/comp/ppcomp.html#4
    solartype: Solarobj, // Type enum of the solar obj
    coords: CartesianCoords,
    n0: f32, nc: f32, // N0 = longitude of the ascending node (deg).  Nc = rate of change in deg/day
    i0: f32, ic: f32, // inclination to the ecliptic (deg)
    w0: f32, wc: f32, // argument of perihelion (deg)
    a0: f32, ac: f32, // semi-major axis, or mean distance from Sun (AU)
    e0: f32, ec: f32, // eccentricity (0=circle, 0..1=ellipse, 1=parabola)
    m0: f32, mc: f32  // M0 = mean anomaly  (deg) (0 at perihelion; increases uniformly with time).  Mc ("mean motion") = rate of change
}

#[derive(Debug)]
pub struct Earth {
    solartype: Solarobj,
    coords: CartesianCoords
}

#[derive(Debug)]
pub struct Sun {
    solartype: Solarobj,
    coords: CartesianCoords
}

#[derive(Debug)]
pub struct CartesianCoords {
    is_meters: bool,
    heliocentric: bool, // False if geocentric
    xh: f32, // X location in meters
    yh: f32, // Y location in meters
    zh: f32  // Z location in meters
}

impl CartesianCoords {
    fn to_meters(&mut self){
        if !self.is_meters {
            self.is_meters = true;
            self.xh *= 149600000000f32;
            self.zh *= 149600000000f32;
            self.yh *= 149600000000f32;
        }
    }

    fn to_au(&mut self){
        if self.is_meters {
            self.is_meters = false;
            self.xh /= 149600000000f32;
            self.zh /= 149600000000f32;
            self.yh /= 149600000000f32;
        }
    }
}

/// Provides utilities for calculating planetary bodies with a Kepler model
mod kepler_utilities {
    use std::f32::{self, consts};
    use super::CartesianCoords;
    
    /**
     *  Calculate the eccentric anomaly for a given body.
     * ### Arguments
     * * 'e' - TODO
     * * 'm' - TODO
     * 
     * ### Returns
     *      The eccentric anomaly for the provided input parameters.
     */
    pub fn eccentric_anomaly(e: f32, m: f32) -> f32 {

        let deg_from_rad = 180f32 / consts::PI;
        let mut ecc: f32 = m + (e * sin_deg!(m) * (1f32 + (e * cos_deg!(m))));

        loop {
            let f: f32 = ecc - (ecc - (deg_from_rad * e * sin_deg!(ecc)) - m) / (1f32 - e * cos_deg!(ecc));
            let error = (f - ecc).abs();
            ecc = f;

            if error < 1.0e-8 { break; }
        };

        ecc
    }

    pub fn lunar_pertub(xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords{

        // TODO add perturb code
        CartesianCoords{xh, yh, zh, is_meters: false, heliocentric: false}
    }
}

pub trait KeplerModel{

    fn update_ecliptic_cartesian_coords(&mut self, day: f32) -> ();

    fn perturb(&self, xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords{
        CartesianCoords{xh, yh, zh, is_meters: false, heliocentric: true}
    }
}

impl KeplerModel for PlanetPL{

    fn update_ecliptic_cartesian_coords(&mut self, day: f32) -> (){
        // Default impl
        let a = self.a0 + (day * self.ac);
        let e = self.e0 + (day * self.ec);
        let m_u = self.m0 + (day * self.mc);
        let n_u = self.n0 + (day * self.nc);
        let w = self.w0 + (day * self.wc);
        let i = self.i0 + (day * self.ic);
        let ecc = kepler_utilities::eccentric_anomaly(e, m_u);

        let xv = a * (cos_deg!(ecc) - e);
        let yv = a * ((1.0f32 - e*e).sqrt() * sin_deg!(ecc));

        let v = atan2_deg!(yv, xv); // True anomaly in degrees: the angle from perihelion of the body as seen by the Sun.
        let r = (xv*xv + yv*yv).sqrt(); // Distance from the Sun to the planet in AU

        let cos_n  = cos_deg!(n_u);
        let sin_n  = sin_deg!(n_u);
        let cosi  = cos_deg!(i);
        let sini  = sin_deg!(i);
        let cos_vw = cos_deg!(v + w);
        let sin_vw = sin_deg!(v + w);
    
        // Now we are ready to calculate (unperturbed) ecliptic cartesian heliocentric coordinates.
        let xh = r * (cos_n * cos_vw - sin_n * sin_vw * cosi);
        let yh = r * (sin_n * cos_vw + cos_n * sin_vw * cosi);
        let zh = r * sin_vw * sini;

        self.coords = self.perturb(xh, yh, zh, day);
        self.coords.to_meters();
    }

    /**
     * Calculates additional perturbations on top of main heliocentric position calculation.
     * Matches PlanetPL bodies using the type enum.
     *  
     * ### Arguments
     *  * 'xh' - X coord
     *  * 'yh' - Y coord
     *  * 'zh' - Z coord
     *  8 'day' - Day value
     * 
     * ### Returns
     *      Cartesian coords with the added perturbations.
     */
    fn perturb(&self, xh: f32, yh: f32, zh: f32, day: f32) -> CartesianCoords {
        match &self.solartype {
            Solarobj::Moon{attr} => kepler_utilities::lunar_pertub(xh, yh, zh, day),
            _ => CartesianCoords{xh, yh, zh, is_meters: false, heliocentric: true}
        }
    }
}

impl KeplerModel for Earth {

    /** 
     * Calculate the position of Earth relative to the Sun. 
     * See: http://cosinekitty.com/astronomy.js
     * 
     *  ### Arguments
     * * 'day' - Day as an f32
     * 
     * ### Return
     *     The coordinates of Earth at the provided time.
     */
    fn update_ecliptic_cartesian_coords(&mut self, day: f32) -> () {

        let d = day - 1.5;
        // Julian centuries since J2000.0
        let t = d / 36525.0;
        // Sun's mean longitude, in degrees
        let l_0 = 280.46645 + (36000.76983 * t) + (0.0003032 * t * t);
        // Sun's mean anomaly, in degrees
        let m_0 = 357.52910 + (35999.05030 * t) - (0.0001559 * t * t) - (0.00000048 * t * t * t);

        let c = // Sun's equation of center in degrees
            (1.914600 - 0.004817 * t - 0.000014 * t * t) * sin_deg!(m_0) +
            (0.01993 - 0.000101 * t) * sin_deg!(2f32 * m_0) +
            0.000290 * sin_deg!(3f32 * m_0);
        
        let ls = l_0 + c; // true elliptical longitude of Sun

         // The eccentricity of the Earth's orbit.
        let e = 0.016708617 - t * (0.000042037 + (0.0000001236 * t));
        // distance from Sun to Earth in astronomical units (AU)
        let distance_in_au = (1.000001018 * (1f32 - e * e)) / (1f32 + e * cos_deg!(m_0 + c));
        let x = -distance_in_au * cos_deg!(ls);
        let y = -distance_in_au * sin_deg!(ls);

        // the Earth's center is always on the plane of the ecliptic (z=0), by definition!
        self.coords = CartesianCoords{xh: x, yh: y, zh: 0f32, is_meters: false, heliocentric: true};
        self.coords.to_meters();
    }
}


impl KeplerModel for Sun {
    fn update_ecliptic_cartesian_coords(&mut self, day:f32) -> (){

    }
}

/**
 * Create the sun as a PlanetPL struct
 *
 * ### Return
 *      A newly crafted sun object
 */
pub fn make_sun() -> Sun {

    let solar_trait = Solarobj::Sun{attr: SolarAttr{radius: 6.95700e8, mass: 1.9891e30}};

    let sun_body = Sun{solartype: solar_trait, coords: CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, is_meters: true, heliocentric: true}};

    sun_body
}

pub fn make_earth(day: f32) -> Earth {

    // Completely not allowed, will cause wildly incorrect planetary calculations.
    if day < 0f32 {panic!("Provided day value is below 0.")}

    let solar_trait = Solarobj::Earth{attr: SolarAttr{radius: 6.3781e6, mass: 5.9722e24}};

    let mut earth_body = Earth{solartype: solar_trait,
                           coords: CartesianCoords{xh: 0f32, yh: 0f32, zh: 0f32, 
                                                   is_meters: false, heliocentric: true}};

    earth_body.update_ecliptic_cartesian_coords(day);

    earth_body
}
