// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

//! Documentation is currently under construction.

// TODO: Make documentation not be "currently under construction." :P

#![allow(unused_imports)]
#![allow(non_snake_case)]
#![allow(unused_parens)]
#![deny(unused_must_use)]

extern crate num;
extern crate test;

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

use std::collections::hash_map::{HashMap,Hasher};
use std::num::{ToPrimitive,FromPrimitive};
use std::hash::Hash;
use std::ops::Shr;

use num::{Zero,One,Integer};

mod primes;
mod factorization;
mod factorizer;
mod factorizer_dixon;
mod factorizer_pollard;
mod util;


/// Factors a number using `DefaultFactorizer`.
pub fn factorize<T>(x: T) -> Factorization<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher>
{
	DefaultFactorizer.factorize(x)
}
