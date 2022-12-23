//! (Inverse) Weibull Distribution.
use libm::{exp, log, pow};

use crate::dist::distutils::*;

use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;
use rand::Rng;

/// Fréchet Dist. struct
#[derive(Clone, Copy)]
pub struct Weibull {
    pub loc:   f64, // location parameter, $\in \mathbb{R}$
    pub scale: f64, // scale parameter, $> 0$
    pub shape: f64, // shape parameter, $> 0$
}

impl Weibull {
    /// Create Weibull Distribution given location (loc), scale and shape parameter.
    /// The scale and shape parameter must be larger than 0.
    #[inline]
    pub fn new(loc: f64, scale: f64, shape: f64) -> Self {
        domain!(scale > 0.0 && shape > 0.0);
        Weibull{loc, scale, shape}
    }

    /// Obtain the location parameter
    #[inline(always)]
    pub fn loc(&self) -> f64 {
        self.loc
    }

    /// Obtain the scale parameter
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Obtain the shape parameter
    pub fn shape(&self) -> f64 {
        self.shape
    }

}

impl DistQuant for Weibull {
    /// CDF: $F(x) = \exp \left \{ - \left (  - \left ( \frac{x - loc}{ scale } \right) \right)^{shape}  \right \} $
    /// for $x < loc$, $loc \in \mathbb{R}$, $scale > 0$ and $shape > 0$.
    fn cdf(&self, x: f64) -> f64 {
        domain!(x < self.loc && self.scale > 0.0 && self.shape > 0.0);
        let y: f64 = (x - self.loc) / self.scale;
        exp(- pow(-y, self.shape))
    }
    
    /// PDF of the Weibull distribution.
    /// $$f(x) = \frac{shape}{scale} \left ( - \frac{x - loc}{scale} \right)^{shape -1} \cdot F(x) $$
    fn pdf(&self, x: f64) -> f64 {
        domain!(x < self.loc && self.scale > 0.0 && self.shape > 0.0);
        let y: f64 = (x - self.loc) / self.scale;
        let pow_const: f64 = self.shape / self.scale;
        pow_const * pow(-y, self.shape- 1.0 ) * exp(- pow(-y, self.shape))
    }

    /// Quantile (inverse CDF) function.
    /// $F^{-1}(x) = - scale \cdot \left(\log x  \right)^{\frac{1}{shape}} + loc$
    fn quantile(&self, x: f64) -> f64 {
        domain!(x >= 0.0 && x <= 1.0);
        self.loc - self.scale * pow(-log(x), 1.0 / self.shape)
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
    macro_rules! new_weibull(
        ($loc:expr, $scale:expr, $shape:expr) => (Weibull::new($loc, $scale, $shape));
    );
    
    #[test]
    fn weibull_cdf_test() {
        let weib: Weibull = new_weibull!(2.0, 2.0, 2.0);
        let ans: f64 = 0.7788007830714049;
        let cdf_weibull: f64 = weib.cdf(1.0);
        assert_eq!(ans, cdf_weibull);
    }

    #[test]
    fn weibull_pdf_test() {
        let weib: Weibull = new_weibull!(2.0, 2.0, 2.0);
        let ans: f64 = 0.38940039153570244;
        let pdf_weibull: f64 = weib.pdf(1.0);
        assert_eq!(ans, pdf_weibull);
    }

    #[test]
    fn weibull_quantile_test() {
        let weib: Weibull = new_weibull!(2.0, 2.0, 2.0);
        let ans: f64 = 0.8055546158342233;
        let quant_weibull: f64 = weib.quantile(0.7);
        assert_eq!(ans, quant_weibull);
    }
}
