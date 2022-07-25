use std::fmt;
use std::error::Error;
use std::thread::AccessError;

#[derive(Debug)]
/// Custom RStats Error
pub enum RError {
    /// Error usually indicating that an empty Vec has been supplied as an argument
    NoDataError,
    /// Error indicating that a wrong kind of data Vec has been supplied as an argument
    DataError,
    /// Error indicating that a division by zero would occur
    ArithError,
    /// Error obtained by converting from another error
    OtherError
}

impl Error for RError {}

impl fmt::Display for RError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RError::NoDataError => write!(f,"Missing Data"),
            RError::DataError => write!(f,"Wrong Data"),
            RError::ArithError => write!(f,"Arithmetic Error (abnormal f64 value)"),
            RError::OtherError => write!(f,"Other Error")
        }
    }
}

impl From<AccessError> for RError {
    fn from(_: AccessError) -> Self {
        RError::OtherError
    }
}
