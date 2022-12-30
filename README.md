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
You can install the package with `cargo add gevlib` or adding the package directly in the `Cargo.toml` file.

## Details

We quickly present the distributions in question.
You can read more about the family of GEV Distributions [here](https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution). 

### Gumbel Distribution

The Gumbel distribution is defined, for location parameter $\mu \in \mathbb{R}$ and scale parameter $\sigma > 0$, with the CDF
$$F(x) = \exp \left ( - \exp  \left ( - \frac{x - \mu}{\sigma}   \right ) \right ), \quad x \in \mathbb{R}. $$

### Fréchet Distribution

The Fréchet distribution is defined, for location parameter $\mu \in \mathbb{R}$, scale and shape parameters $\sigma, \zeta >0$, with the CDF
$$F(x) = \exp \left ( - \left ( \frac{x - \mu}{\sigma} \right)^{-\zeta} \right ), \quad x \in \mathbb{R}. $$

### (Inverse) Weibull Distribution

The Weibull Distribution, which is in fact the Inverse Weibull distribution, is defined for location parameter $\mu \in \mathbb{R}$, scale and shape parameters $\sigma, \zeta >0$, with the CDF
$$F(x) = \exp \left ( - \left (  - \left ( \frac{x - \mu}{ \sigma } \right) \right)^{\zeta}  \right ), \quad x \in \mathbb{R}. $$

### GEV Distribution

The GEV Distribution generalizes all of the above distributions. It is defined for location parameter $\mu \in \mathbb{R}$, scale parameter $\sigma >0$ and shape parameters $\zeta \in \mathbb{R}$ with the CDF of $F(x) = \exp \left ( t(x) \right)$ for $1 + \zeta \left( \frac{x - \mu}{\sigma} \right ) > 0$, where
$$t(x) = \left( 1 + \zeta \left( \frac{x - \mu}{\sigma} \right) \right)^{- \frac{1}{\zeta}} \quad \text{if} \quad \zeta \neq 0, $$
and
$$t(x) = \exp \left ( - \frac{x - \mu}{\sigma}  \right ) \quad \text{if} \quad \zeta = 0.$$

In fact, for $\zeta = 0$ we recover the Gumbel distribution, for $\zeta > 0$ we recover the Fréchet distribution and for $\zeta < 0$ we have the Weibull distribution.

# To do
- [ ] add macros to create instances of each distribution.
- [ ] Clean code? 
- [ ] Add Rust readme file
- [ ] Write better documentation
