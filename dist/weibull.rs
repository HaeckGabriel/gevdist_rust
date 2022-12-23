 
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

    /// Return a randomly generated value from the Weibull distribution.
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

    // quick macro to create the instance of the Weibull Distribution
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
