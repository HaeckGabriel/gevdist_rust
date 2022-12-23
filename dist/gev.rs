//! GEV Distribution
use libm::{exp, log, pow};


use crate::dist::distutils::*;

use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;
use rand::Rng;

/// Fréchet Dist. struct
#[derive(Clone, Copy)]
pub struct GEV {
    loc:   f64, // location parameter, $\in \mathbb{R}$
    scale: f64, // scale parameter, $> 0$
    shape: f64, // shape parameter, $\in \mathbb{R}$
}

impl GEV {
    /// Create GEV Distribution given location (loc), scale and shape parameter.
    /// The scale parameter must be larger than 0.
    #[inline]
    pub fn new(loc: f64, scale: f64, shape: f64) -> Self {
        domain!(scale > 0.0);
        GEV{loc, scale, shape}
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

    /// t(x) function that depends on if the shape parameter is 0 or not.
    /// t(x) = \exp \left(x) = \left( 1 + \zeta \left( \frac{x - loc}{ scale} \right) \right)^{- \frac{1}{\zeta}}$$ if $\zeta \neq 0$,
    /// or $t(x) = \exp \left \{ - \frac{x - loc}{ scale}  \right \}$ if $\zeta = 0$t
    #[inline(always)]
    fn t_func(&self, x: f64) -> f64 {
        let y: f64 = (x - self.loc) / self.scale;
        if self.shape == 0.0 {
            exp(- y)
        } else {
            pow(1.0 + self.shape * y , - 1.0 / self.shape)
        }
    }

}

impl DistQuant for GEV {
     /// CDF: $F(x) = \exp \left \{ - t_func(x) \right \} $
    /// for $1 + shape \left( \frac{x - loc}{ scale} > 0$
    fn cdf(&self, x: f64) -> f64 {
        domain!(1.0 + self.shape * ( (x - self.loc ) / self.scale ) > 0.0 && self.scale > 0.0); // need  $1 + shape \left( \frac{x - loc}{ scale} > 0$ 
        let t_val: f64 = self.t_func(x);
        exp(- t_val)
    }
    
    /// PDF of the GEV distribution.
    /// $$ f(x) = \frac{1}{ scale } t_func(x)^{\zeta + 1} \cdot F(x)  $$
    fn pdf(&self, x: f64) -> f64 {
        domain!(1.0 + self.shape * ( (x - self.loc ) / self.scale ) > 0.0 && self.scale > 0.0); // need  $1 + shape \left( \frac{x - loc}{ scale} > 0$ 
        let mult_const: f64 = 1.0 / self.scale;
        let t_val: f64 = self.t_func(x);
        mult_const * pow(t_val, self.shape + 1.0) * exp(- t_val)
    }

    /// Quantile (inverse CDF) function.
    /// If $shape = 0$, $F^{-1}(x) = loc - scale * \log(- \log x)$
    /// o.w. we have $\frac{scale}{shape} * (- \log x)^{- shape} - \frac{scale}{shape} + loc$
    fn quantile(&self, x: f64) -> f64 {
        domain!(x >= 0.0 && x <= 1.0);
        if self.shape == 0.0 {
            - self.scale * log( - log(x)) + self.loc
        } else {
            let mult_const: f64 = self.scale / self.shape;
            mult_const * pow(- log(x) , - self.shape) - mult_const + self.loc
        }
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
    macro_rules! new_gev(
        ($loc:expr, $scale:expr, $shape:expr) => (GEV::new($loc, $scale, $shape));
    );
    
    #[test]
    fn gev_cdf_test_one() {
        let gev: GEV = new_gev!(2.0, 2.0, 2.0);
        let ans: f64 = 0.49306869139523984;
        let cdf_gev: f64 = gev.cdf(3.0);
        assert_eq!(ans, cdf_gev);
    }

    #[test]
    fn gev_pdf_test_one() {
        let gev: GEV = new_gev!(2.0, 2.0, 2.0);
        let ans: f64 = 0.08716305381908777;
        let pdf_gev: f64 = gev.pdf(3.0);
        assert_eq!(ans, pdf_gev);
    }

    #[test]
    fn gev_quantile_test_one() {
        let gev: GEV = new_gev!(2.0, 2.0, 2.0);
        let ans: f64 = 8.860583704300595;
        let quant_gev: f64 = gev.quantile(0.7);
        assert_eq!(ans, quant_gev);
    }

    /// now same series of test but with shape = 0.
    #[test]
    fn gev_cdf_test_two() {
        let gev: GEV = new_gev!(2.0, 2.0, 0.0);
        let ans: f64 = 0.545239211892605;
        let cdf_gev: f64 = gev.cdf(3.0);
        assert_eq!(ans, cdf_gev);
    }

    #[test]
    fn gev_pdf_test_two() {
        let gev: GEV = new_gev!(2.0, 2.0, 0.0);
        let ans: f64 = 0.16535214944520904;
        let pdf_gev: f64 = gev.pdf(3.0);
        assert_eq!(ans, pdf_gev);
    }

    #[test]
    fn gev_quantile_test_two() {
        let gev: GEV = new_gev!(2.0, 2.0, 0.0);
        let ans: f64 = 4.061860866317446;
        let quant_gev: f64 = gev.quantile(0.7);
        assert_eq!(ans, quant_gev);
    }

}
