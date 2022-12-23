<h1 align="center"> Extreme Value Distributions (Rust) </h1>

<h4 align="center"> The family of Extreme Value Distributions in Rust. </h4>

<p align="center">
  <a href="https://crates.io/crates/gevlib">
    <img src="https://img.shields.io/badge/Version-0.1.4-blueviolet?style=for-the-badge&logo=Rust">
  </a>
  <a href="https://github.com/HaeckGabriel/gevdist_rust">
    <img src="https://img.shields.io/github/downloads/HaeckGabriel/gevdist_rust/total?label=Downloads&logo=Github&style=for-the-badge&color=blue">
  </a>
  <a href="https://crates.io/crates/gevlib">
    <img src="https://img.shields.io/crates/d/gevlib?label=Crates%20Downloads&logo=Rust&style=for-the-badge&color=9cf">
  </a>
</p>

Basic Distributional Quantities (CDF, PDF, Quantile and Random Generation) for the Gumbel, Fréchet, (inverse) Weibull and GEV Distributions.

<p align="center">
  <a href="#Installation">Installation</a> •
  <a href="#Details">Details</a>
</p>

## Installation
You can install the package with `cargo add gevlib` or adding the pacakge directly in the `Cargo.toml` file.

## Details

We quickly present the distributions in question.
You can read more about the family of GEV Distributions [here](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution). 

### Gumbel Distribution

The Gumbel distribution is defined, for location parameter $\mu$ and scale parameter $\theta > 0$, with the CDF
$$F(x) = \exp \left \{ - \exp  \{- \frac{x - \mu}{\theta}  \} \right \}, \quad x \in \mathbb{R}. $$

# To do
- [ ] add macros to create instances of each distribution.
- [ ] Clean code? 
- [ ] Add Rust readme file
- [ ] Write better documentation
