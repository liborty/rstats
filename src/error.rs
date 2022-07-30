use std::fmt;
use std::error::Error;
use std::thread::AccessError;
use std::fmt::{Debug,Display};

#[derive(Debug)]
/// Custom RStats Error
pub enum RError<T> where T:?Sized+Debug+'static {
    /// Error indicating that an insufficient data has been supplied
    NoDataError(&'static T),
    /// Error indicating that a wrong kind/size of data has been supplied
    DataError(&'static T),
    /// Error indicating an invalid result, such as division by zero
    ArithError(&'static T),
    /// Error obtained by converting from another error
    OtherError(&'static T)
}

impl<T> Error for RError<T> where T:Sized+Debug+Display {}

impl<T> fmt::Display for RError<T> where T:Sized+Debug+Display {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RError::NoDataError(s) => write!(f,"Missing or Insufficient Data: {}",s),
            RError::DataError(s) => write!(f,"Wrong Data: {}",s),
            RError::ArithError(s) => write!(f,"Arithmetic Error: {}",s),
            RError::OtherError(s) => write!(f,"Converted from {}",s)
        }
    }
}

impl From<AccessError> for RError<& 'static str> {
    fn from(_: AccessError) -> Self {
        RError::OtherError(&"AccessError")
    }
}
