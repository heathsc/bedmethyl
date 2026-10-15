use clap::{ValueEnum, builder::PossibleValue};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ValueDelim {
    Tab,
    Space,
    Comma,
    Semicolon,
    Slash,
    Split,
    None,
}

impl ValueEnum for ValueDelim {
    fn value_variants<'a>() -> &'a [Self] {
        &[
            Self::Tab,
            Self::Space,
            Self::Comma,
            Self::Semicolon,
            Self::Slash,
        ]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::Tab => Some(PossibleValue::new("tab")),
            Self::Space => Some(PossibleValue::new("space")),
            Self::Comma => Some(PossibleValue::new("comma")),
            Self::Semicolon => Some(PossibleValue::new("semicolon")),
            Self::Slash => Some(PossibleValue::new("slash")),
            _ => None,
        }
    }
}

impl ValueDelim {
	pub fn get_delim_char(&self) -> char {
		match self {
			ValueDelim::Space => ' ',
			ValueDelim::Comma => ',',
			ValueDelim::Semicolon => ';',
			_ => '\t',
		}
	}
}
