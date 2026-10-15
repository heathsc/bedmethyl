use clap::{builder::PossibleValue, ValueEnum};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputType {
    NonconvConv,
    NonconvCov,
    MethCov,
    Meth,
}

impl ValueEnum for OutputType {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::NonconvConv,
            Self::NonconvCov,
            Self::MethCov,
            Self::Meth,
        ]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::NonconvConv => Some(PossibleValue::new("non_conv-conv")),
            Self::NonconvCov => Some(PossibleValue::new("non_conv-cov")),
            Self::MethCov => Some(PossibleValue::new("meth-cov")),
            Self::Meth => Some(PossibleValue::new("meth")),
        }
    }
}
