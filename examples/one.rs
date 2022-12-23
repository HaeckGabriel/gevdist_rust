use gevlib;

use gevlib::dist::gumbel::Gumbel;
use gevlib::dist::distutils::*;

fn main() {
    let gumbel_test: Gumbel = Gumbel {loc: 0.5, scale: 2.0};
    let cdf: f64 = gumbel_test.cdf(2.0);
    println!("CDF value at 2: {}", cdf);
}
