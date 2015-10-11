// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

//! Documentation is currently under construction.

// TODO: Make documentation not be "currently under construction." :P

#![feature(test)]
#![feature(step_by)]

#![allow(unused_imports)]
#![allow(non_snake_case)]
#![allow(unused_parens)]
#![deny(unused_must_use)]

extern crate num;
extern crate test;
extern crate rand;
extern crate bit_set;

pub use primes::PrimeTester;
pub use primes::PrimeSieve;
pub use primes::MillerRabinTester;
pub use factorizer::Factorizer;
pub use factorizer::TrialDivisionFactorizer;
pub use factorizer::StubbornFactorizer;
//pub use factorizer::FermatFactorizer;
pub use factorizer_dixon::DixonFactorizer;
pub use factorizer_pollard::PollardBrentFactorizer;
//pub use factorizer::GeneralFactorizer;
pub use factorizer::SafeFactorizer;
pub use factorizer::DefaultFactorizer;
pub use factorizer::ListFactorizer;
pub use factorization::Factorization;
pub use util::isqrt;
pub use util::gcd;

use util::literal;

use std::iter::FromIterator;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::Shr;

use num::{Zero,One,Integer};
use num::{FromPrimitive,ToPrimitive};

mod primes;
mod factorization;
mod factorizer;
mod factorizer_dixon;
mod factorizer_pollard;
mod util;


/// Factors a number using `DefaultFactorizer`.
pub fn factorize<T>(x: T) -> Factorization<T>
 where T: Eq + Clone + Zero + One + Integer + Shr<usize, Output=T> + Hash + ToPrimitive + FromPrimitive
{
	DefaultFactorizer.factorize(x)
}

/// Collects all primes up to a limit (inclusive)
///
/// # Example
///
/// ```
/// let v: Vec<u64> = primes_upto(17);
/// assert_eq!(v, vec![2, 3, 5, 7, 11, 13, 17]);
/// assert!(false);
/// ```
pub fn primes_upto<Out>(limit: usize) -> Out
 where Out: FromIterator<u64>
{
	let max = limit + 1;
	let sieve = primes::PrimeSieve::new(max);
	(2..max).step_by(2)
		.filter(|&i| sieve.is_prime(&i))
		.map(|i| i as u64)
		.collect()
}

