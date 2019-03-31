
#[macro_export]
macro_rules! sin_deg {
    ($x: expr) => {
        $x.to_radians().sin().to_degrees()
    };
}

#[macro_export]
macro_rules! cos_deg {
    ($x: expr) => {
        $x.to_radians().cos().to_degrees()
    };
}

#[macro_export]
macro_rules! atan2_deg {
    ($x : expr, $y : expr) => {
        $x.to_radians().atan2($y.to_radians()).to_degrees()
    };
}