//! Functions used for the distributions.

// To verify that specific domains of parameters are respected.
// In fact this macro comes form: https://docs.rs/probability/latest/src/probability/lib.rs.html#1-35

#[macro_use]
pub mod macros {
    macro_rules! domain( // should you consider an error enum instead? What is faster? More understandable? One variant for param and one for quantile domain
        ($requirement:expr) => (debug_assert!($requirement));
        ($requirement:expr, $code:expr) => (debug_assert!($code, stringify!($requirement)));
    );
}

/// Distributional Quantity trait (i.e. each distribution will provide each of the following)
pub trait DistQuant {
    fn cdf(&self, x: f64) -> f64;      // cumulative density function (CDF)
    fn pdf(&self, x: f64) -> f64;      // probability density function (PDF)
    fn quantile(&self, x: f64) -> f64; // quantile function (i.e. inverse CDF)
    fn random(&self, seed: RandomSeed) -> f64;           // randomly generated value of the distribution
}

/// Seeding for the random generation of the distributions.
/// Can either be Empty (i.e. use random seed) or with a given u64 seed.
pub enum RandomSeed {
    Empty,
    Seed(u64),
}

/// Default seed is Empty (i.e. random seed)
impl Default for RandomSeed {
    fn default() -> Self {
        RandomSeed::Empty
    }
}

impl RandomSeed {
    /// To be able to retrieve the given seed, if ever needed.
    pub fn get_seed(&self) -> Option<u64> {
        match self {
            RandomSeed::Empty => None,
            RandomSeed::Seed(val) => Some(*val),
        }
    }
}
