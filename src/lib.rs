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

extern crate num;
extern crate test;

use std::collections::hash_map::{HashMap,Hasher};
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;

pub use primes::PrimeTester;
pub use primes::PrimeSieve;
pub use factorizer::Factorizer;
pub use factorizer::TrialDivisionFactorizer;
//pub use factorizer::FermatFactorizer;
//pub use factorizer::DixonFactorizer;
//pub use factorizer::GeneralFactorizer;
pub use factorizer::SafeFactorizer;
pub use factorizer::DefaultFactorizer;
pub use factorizer::FactorStore;
pub use factorization::Factorization;

mod primes;
mod factorization;
mod factorizer;
mod util;

