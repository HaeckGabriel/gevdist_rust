//! The Gumbel Distribution.
use libm::{exp, log};

use crate::dist::distutils::*;

use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;
use rand::Rng;

//extern crate libm::exp;

/// Gumbel Dist. struct
#[derive(Clone, Copy)]
pub struct Gumbel {
    /// location parameter
    pub loc:   f64,
    /// scale parameter, must be positive
    pub scale: f64,
}

impl Gumbel {
    /// Create an instance of the Gumbel Distribution given location (loc) and scale parameter.
    /// The scale parameter must be larger than 0.
    #[inline]
    pub fn new(loc: f64, scale: f64) -> Self {
        domain!(scale > 0.0);
        Gumbel{loc, scale}
    }

    /// Obtain the location parameter
    #[inline(always)]
    pub fn loc(&self) -> f64 {
        self.loc
    }

    /// Obtain the scale parameter
    #[inline(always)]
    pub fn scale(&self) -> f64 {
        self.scale
    }

}

/// Distributional Quantities for the Gumbel Distribution.
impl DistQuant for Gumbel {

    /// CDF: $F(x) = \exp \left \{ - \exp \left \{- \frac{x - \loc}{\scale}  \right \} \right \} $
    /// for $x \in \mathbb{R}$
    fn cdf(&self, x: f64) -> f64 {
        let y: f64 = (x - self.loc) / self.scale; 
        exp(- exp(-y))
    }
    
    /// PDF of the Gumbel distribution.
    /// $f(x) = \frac{1}{\scale} \exp \left \{- \frac{x - \loc}{\scale} \right \} \exp \left \{- \exp \left \{ - \frac{x - \loc}{\scale} \right \} \right \}$
    fn pdf(&self, x: f64) -> f64 {
        let y: f64 = (x - self.loc) / self.scale;
        let constant: f64 = 1.0 / self.scale;
        constant * exp(- y) * exp(- exp(-y))
    }

    /// Quantile (inverse CDF) function.
    /// $F^{-1}(x) = \loc - \scale \log \left ( - \log \left ( x \right ) \right )$
    fn quantile(&self, x: f64) -> f64 {
        domain!(x >= 0.0 && x <= 1.0);
        self.loc - self.scale * log(-log(x))
    }

    /// Return a randomly generated value from the Gumbel distribution.
    fn random(&self, seed: RandomSeed) -> f64 {
        
        let mut rng = match seed {
            RandomSeed::Empty => ChaCha8Rng::from_entropy(),
            RandomSeed::Seed(val) => ChaCha8Rng::seed_from_u64(val), // ChaCha8Rng implements the SeedableRng trait
        };
        let rand_quant: f64 = rng.gen::<f64>(); // generate randomly from U(0,1)
        self.quantile(rand_quant) // then plug that random uniform into the quantile.
    }
}   

/// tests
#[cfg(test)]
mod tests {
    use super::*;

    // quick macro to create the instance of the Gumbel Distribution
    macro_rules! new_gumbel(
        ($loc:expr, $scale:expr) => (Gumbel::new($loc, $scale));
    );
    
    #[test]
    fn gumbel_cdf_test() {
        let gumb: Gumbel = new_gumbel!(0.5, 2.0);
        let ans: f64 = 0.6235249162568004;
        let cdf_gumb: f64 = gumb.cdf(2.0);
        assert_eq!(ans, cdf_gumb);
    }

    #[test]
    fn gumbel_pdf_test() {
        let gumb: Gumbel = new_gumbel!(0.5, 2.0);
        let ans: f64 = 0.14726615762017733;
        let pdf_gumb: f64 = gumb.pdf(2.0);
        assert_eq!(ans, pdf_gumb);
    }

    #[test]
    fn gumbel_quantile_test() {
        let gumb: Gumbel = new_gumbel!(0.5, 2.0);
        let ans: f64 = 2.5618608663174456;
        let gumb_quant: f64 = gumb.quantile(0.7);
        assert_eq!(ans, gumb_quant);
    }
}
