use std::fmt;
use std::error::Error;
use std::thread::AccessError;
use std::fmt::{Debug,Display};

#[derive(Debug)]
/// Custom RStats Error
pub enum RError<T> where T:Sized+Debug {
    /// Error indicating that insufficient data has been supplied
    NoDataError(T),
    /// Error indicating that a wrong kind/size of data has been supplied
    DataError(T),
    /// Error indicating an invalid result, such as an attempt at division by zero
    ArithError(T),
    /// Other error converted to RError
    OtherError(T)
}

impl<T> Error for RError<T> where T:Sized+Debug+Display {}

impl<T> fmt::Display for RError<T> where T:Sized+Debug+Display {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RError::NoDataError(s) => write!(f,"Missing or insufficient data: {}",s),
            RError::DataError(s) => write!(f,"Wrong data: {}",s),
            RError::ArithError(s) => write!(f,"Arithmetic error: {}",s),
            RError::OtherError(s) => write!(f,"Converted from {}",s)
        }
    }
}

/// Example 'From' implementation for converting to RError
impl From<AccessError> for RError<& 'static str> {
    fn from(_: AccessError) -> Self {
        RError::OtherError("AccessError")
    }
}
