
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