use medians::error::MedError;
use ran::RanError;
use std::error::Error;
use std::fmt;
use std::fmt::{Debug, Display};
use std::thread::AccessError;

/// Shorthand type for returned errors with String (message) payload
pub type RE = RError<String>;

#[derive(Debug)]
/// Custom RStats Error
pub enum RError<T>
where
    T: Sized + Debug,
{
    /// Error indicating that insufficient data has been supplied
    NoDataError(T),
    /// Error indicating that a wrong kind/size of data has been supplied
    DataError(T),
    /// Error indicating an invalid result, such as an attempt at division by zero
    ArithError(T),
    /// Other error converted to RError
    OtherError(T),
}

impl<T> Error for RError<T> where T: Sized + Debug + Display {}

impl<T> fmt::Display for RError<T>
where
    T: Sized + Debug + Display,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RError::NoDataError(s) => write!(f, "Missing or insufficient data: {s}"),
            RError::DataError(s) => write!(f, "Wrong data: {s}"),
            RError::ArithError(s) => write!(f, "Arithmetic error: {s}"),
            RError::OtherError(s) => write!(f, "Converted from {s}"),
        }
    }
}
/// Automatically converting any RanError to RError::OtherError
impl From<RanError<String>> for RError<String> {
    fn from(e: RanError<String>) -> Self {
        RError::OtherError(format!("RanError: {e}"))
    }
}

/// Automatically converting any MedError to RError::OtherError
impl From<MedError<String>> for RError<String> {
    fn from(e: MedError<String>) -> Self {
        RError::OtherError(format!("MedError: {e}"))
    }
}

/// Example 'From' implementation for converting to RError
impl From<AccessError> for RError<String> {
    fn from(e: AccessError) -> Self {
        RError::OtherError(format!("AccessError: {e}"))
    }
}

// 'From' implementation for converting ioerror to RError
impl From<std::io::Error> for RError<String> {
    fn from(e: std::io::Error) -> Self {
        RError::OtherError(format!("IOError: {e}"))
    }
}

/// Convenience function for building RError<String>  
/// from short name and payload message, which can be either &str or String
pub fn re_error(kind: &str, msg: impl Into<String>) -> RE {
    match kind {
        "empty" => RError::NoDataError(msg.into()), 
        "size"  => RError::DataError(msg.into()), 
        "arith" => RError::ArithError(msg.into()),
        "other" => RError::OtherError(msg.into()),
        _ => RError::OtherError("Wrong error kind given to re_error".into())
    }
}
