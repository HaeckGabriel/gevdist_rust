//! Fréchet Distribution
use libm::{exp, log, pow};

use crate::dist::distutils::*;

use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;
use rand::Rng;

/// Fréchet Dist. struct
#[derive(Clone, Copy)]
pub struct Frechet {
    pub loc:   f64, // location parameter, $\in \mathbb{R}$
    pub scale: f64, // scale parameter, $> 0$
    pub shape: f64, // shape parameter, $\in \mathbb{R}$
}

impl Frechet {
    /// Create Frechet Distribution given location (loc), scale and shape parameter.
    /// The scale and shape parameter must be larger than 0.
    #[inline]
    pub fn new(loc: f64, scale: f64, shape: f64) -> Self {
        domain!(scale > 0.0 && shape > 0.0);
        Frechet{loc, scale, shape}
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

    /// Obtain the shape parameter
    #[inline(always)]
    pub fn shape(&self) -> f64 {
        self.shape
    }

}

impl DistQuant for Frechet {
    /// CDF: $F(x) = \exp \left \{ - \left ( \frac{x - loc}{scale} \right)^{-shape} \right \} $
    /// for $x > loc$
    fn cdf(&self, x: f64) -> f64 {
        domain!(x > self.loc);
        let y: f64 = (x - self.loc) / self.scale;
        exp(- pow(y, - self.shape))
    }
    
    /// PDF of the Frechet distribution.
    /// $$f (x) = \frac{shape}{scale} \left(\frac{ x - loc }{scale}\right)^{-1 - shape} \exp \left \{ - \left( \frac{x - loc}{scale} \right)^{- shape}  \right \} $$
    fn pdf(&self, x: f64) -> f64 {
        domain!(x > self.loc);
        let y: f64 = (x - self.loc) / self.scale;
        let pow_const: f64 = self.shape / self.scale;
        pow_const * pow(y, -1.0 - self.shape) * exp(- pow(y, - self.shape))
    }

    /// Quantile (inverse CDF) function.
    /// $F^{-1}(x) = loc + scale \left(- \log x \right )^{- \frac{1}{shape}}$
    fn quantile(&self, x: f64) -> f64 {
        domain!(x >= 0.0 && x <= 1.0);
        self.loc + self.scale * pow(-log(x), - 1.0 / self.shape)
    }

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
    macro_rules! new_frechet(
        ($loc:expr, $scale:expr, $shape:expr) => (Frechet::new($loc, $scale, $shape));
    );
    
    #[test]
    fn frechet_cdf_test() {
        let frech: Frechet = new_frechet!(1.0, 0.1, 1.0);
        let ans: f64 = 0.951229424500714;
        let cdf_frechet: f64 = frech.cdf(3.0);
        assert_eq!(ans, cdf_frechet);
    }

    #[test]
    fn frechet_pdf_test() {
        let frech: Frechet = new_frechet!(1.0, 0.1, 1.0);
        let ans: f64 = 0.023780735612517853;
        let pdf_frechet: f64 = frech.pdf(3.0);
        assert_eq!(ans, pdf_frechet);
    }

    #[test]
    fn frechet_quantile_test() {
        let frech: Frechet = new_frechet!(1.0, 0.1, 1.0);
        let ans: f64 = 1.2803673252057128;
        let quant_frechet: f64 = frech.quantile(0.7);
        assert_eq!(ans, quant_frechet);
    }
}
