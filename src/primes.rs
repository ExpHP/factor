// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;

use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::iter::range_step;

use num::{Zero,One,Integer};


pub trait PrimeTester<T>
 where T: Eq + Clone + Zero + One,
{
	fn is_prime(self: &Self, x: &T) -> bool;

	fn is_composite(self: &Self, x: &T) -> bool
	{
		if x == &Zero::zero() || x == &One::one() {
			false
		} else {
			!self.is_prime(x)
		}
	}
}

#[derive(Clone,Show)]
pub struct PrimeSieve
{
	sieve: Vec<bool>,
}

impl PrimeSieve
{
	/// Computes and constructs a Sieve of Eratosthenes capable of testing numbers
	///  below `limit` for primality.
	pub fn new(limit: usize) -> Self
	{
		PrimeSieve {
			sieve: compute_sieve_of_eratosthenes(limit),
		}
	}
}

impl<T> PrimeTester<T>
for PrimeSieve
 where T: Eq + Clone + ToPrimitive + FromPrimitive + Integer,
{
	//! All of the work involved in identifying prime numbers with a PrimeSieve
	//!  is done during the sieve's construction, so a PrimeSieve can test primality
	//!  and compositeness in nominal time.
	//!
	//! # Example
	//!
	//! ```
	//! use factor::{PrimeTester,PrimeSieve};
	//! 
	//! let tester: Box<PrimeTester<i32>> = Box::new(PrimeSieve::new(42));
	//! 
	//! assert!(!tester.is_prime(&0)  && !tester.is_composite(&0));
	//! assert!(!tester.is_prime(&1)  && !tester.is_composite(&1));
	//! assert!(tester.is_prime(&13)  && !tester.is_composite(&13));
	//! assert!(!tester.is_prime(&27) && tester.is_composite(&27));
	//! ```

	/// Tests a number for primality by indexing the underlying vector.
	///
	/// # Panics
	/// 
	/// Panics if ToPrimitive.to_uint() fails (for whatever reason), or if the index
	///  lies outside the array.
	// TODO: test panics
	#[inline]
	fn is_prime(self: &Self, x: &T) -> bool
	{
		self.sieve[x.to_uint().unwrap()]
	}
}

fn compute_sieve_of_eratosthenes(limit: usize) -> Vec<bool>
{
	let mut sieve: Vec<bool> = Vec::new();

	sieve.resize(limit, true); // not sure if idiomatic way to initialize?
	sieve[0] = false;
	sieve[1] = false;

	for p in (2..limit) {
		if sieve[p] {
			for multiple in range_step(p*p, limit, p) {
//			for multiple in (p*p..limit).step(p) {  // future form of above?
				sieve[multiple] = false;
			};
		};
	};
	sieve
}

