use std::fmt;
use std::error::Error;
use std::thread::AccessError;

#[derive(Debug)]
/// Custom RStats Error
pub enum RError {
    /// Error indicating that an insufficient data has been supplied
    NoDataError,
    /// Error indicating that a wrong kind/size of data has been supplied
    DataError,
    /// Error indicating an invalid result, such as division by zero
    ArithError,
    /// Error obtained by converting from another error
    OtherError
}

impl Error for RError {}

impl fmt::Display for RError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RError::NoDataError => write!(f,"Missing or Insufficient Data"),
            RError::DataError => write!(f,"Wrong Data"),
            RError::ArithError => write!(f,"Arithmetic Error, such as abnormal f64 value"),
            RError::OtherError => write!(f,"Converted Error")
        }
    }
}

impl From<AccessError> for RError {
    fn from(_: AccessError) -> Self {
        RError::OtherError
    }
}
