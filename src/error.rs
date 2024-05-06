use medians::MedError;
use ran::RanError;
use std::error::Error;
use std::fmt;
use std::fmt::{Debug, Display};
use std::thread::AccessError;

/// Shorthand type for returned errors with String (message) payload
pub type RE = RError<String>;

#[derive(Debug)]
/// Custom RStats Error
/// Parameter <T> is future proofing, so that any type of argument may be returned.
/// Currently only messages of type <String> and <&str> are used
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
/// Convenience function for building RError::NoDataError(String)  
/// from payload message, which can be either literal `&str` or `String`.
/// `String` allows splicing in values of variables for debugging, using `format!`
pub fn nodata_error<T>(msg: impl Into<String>) -> Result<T,RError<String>> {
    Err(RError::NoDataError(msg.into()))
} 
/// Convenience function for building RError::DataError(String)  
/// from payload message, which can be either literal `&str` or `String`.
/// `String` allows splicing in values of variables for debugging, using `format!`
pub fn data_error<T>(msg: impl Into<String>) -> Result<T,RError<String>> {
    Err(RError::DataError(msg.into()))
}
/// Convenience function for building RError::ArithError(String)  
/// from payload message, which can be either literal `&str` or `String`.
/// `String` allows splicing in values of variables for debugging, using `format!`
pub fn arith_error<T>(msg: impl Into<String>) -> Result<T,RError<String>> {
    Err(RError::ArithError(msg.into()))
} 
/// Convenience function for building RError::ArithError(String)  
/// from payload message, which can be either literal `&str` or `String`.
/// `String` allows splicing in values of variables for debugging, using `format!`
pub fn other_error<T>(msg: impl Into<String>) -> Result<T,RError<String>> {
    Err(RError::OtherError(msg.into()))
} 
/*
/// Convenience function for building RError<String>  
/// from short name and payload message, which can be either literal `&str` or `String`.
/// `String` allows splicing in values of variables for debugging, using `format!`
pub fn re_error<T>(kind: &str, msg: impl Into<String>) -> Result<T,RError<String>> {
    match kind {
        "empty" => Err(RError::NoDataError(msg.into())), 
        "size"  => Err(RError::DataError(msg.into())), 
        "arith" => Err(RError::ArithError(msg.into())),
        "other" => Err(RError::OtherError(msg.into())),
        _ => Err(RError::OtherError("Wrong error kind given to re_error".into()))
    } 
}
*/
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
