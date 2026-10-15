use clap::{ValueEnum, builder::PossibleValue};

#[derive(Debug, Clone, Copy)]
pub enum SimilarityType {
    Correlation,
    Covariance,
    Distance,
}

//                   .possible_values(["correlation", "covariance", "distance"])
impl ValueEnum for SimilarityType {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::Correlation,
            Self::Covariance,
            Self::Distance,
        ]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::Correlation => Some(PossibleValue::new("correlation")),
            Self::Covariance => Some(PossibleValue::new("space")),
            Self::Distance => Some(PossibleValue::new("comma")),
        }
    }
}
