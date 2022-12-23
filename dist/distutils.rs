//! Functions used for the distributions.

/// macro used to ensure that the given domain is valid.
#[macro_use]
pub mod macros {
    macro_rules! domain( // should you consider an error enum instead? What is faster? More understandable? One variant for param and one for quantile domain
        ($requirement:expr) => (debug_assert!($requirement));
        ($requirement:expr, $code:expr) => (debug_assert!($code, stringify!($requirement)));
    );
}

/// Distributional Quantity trait (i.e. each distribution will provide each of the following
/// quantities: the CDF, PDF, Quantile and random generation)
pub trait DistQuant {
    /// Cumulative Distribution Function (CDF)
    fn cdf(&self, x: f64) -> f64;
    /// Probability Density Function (PDF)
    fn pdf(&self, x: f64) -> f64;
    /// Quantile Function
    fn quantile(&self, x: f64) -> f64;
    /// Generate a random value from the distribution
    fn random(&self, seed: RandomSeed) -> f64;
}

/// Seeding for the random generation of the distributions.
/// Can either be Empty (i.e. use random seed) or with a given u64 seed.
pub enum RandomSeed {
    Empty,
    Seed(u64),
}

/// Default Trait
impl Default for RandomSeed {
    /// Default seed is Empty (i.e. random seed)
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
