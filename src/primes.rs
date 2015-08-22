// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;
extern crate test;

use std::iter::range_step;
use std::ops::Shr;
use std::fmt::Debug;

use num::{Zero,One,Integer};
use num::BigUint;
use test::Bencher;

use util::literal;
use util::mod_pow;

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

#[derive(Clone,Debug)]
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
 where T: Eq + Clone + Integer,
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

pub struct MillerRabinTester<T>;

// TODO: implement non-deterministic miller-rabin
// TODO: implement different choices of witnesses for deterministic miller-rabin (and
//         have them validate their input against the max number for which they work)
impl<T> MillerRabinTester<T>
 where T: Eq + Clone + Integer,
{
	/// Produces some set of numbers in the half-open interval `[0,x)` to use
	///  as potential witnesses for the Miller Rabin Test.
	// TODO: a term like "witness" prolly isn't fit for public API without definition
	pub fn collect_witnesses(self: &Self, x: T) -> Vec<T> {
		// There's a couple of options here.
		// Several witnesses can be chosen at random for a strong heuristic test,
		//  or one can use a small set of witnesses that has been demonstrated
		//  to correctly identify all numbers up to some limit.

		// This set of witnesses is known to work for all n < 3,474,749,660,383

		let v = vec![literal(2), literal(3), literal(5), literal(7), literal(11), literal(13)];

		if x > literal(13) {
			v
		} else {
			v.into_iter().filter(|v| v < &x).collect()
		}
	}
}


impl<T> PrimeTester<T>
for MillerRabinTester<T>
 where T: Eq + Clone + Integer + Shr<usize, Output=T> + Debug,
{
	fn is_prime(self: &Self, x: &T) -> bool
	{
		if *x < literal(2) { return false; }
		if *x <= literal(3) { return *x >= literal(2); }

		let neg_one: T = x.clone() - literal(1);

		let (k,d) = decompose_pow2_odd(neg_one.clone());

		'witness: for witness in self.collect_witnesses(x.clone()).into_iter() {
			// TODO: This must be revised when witness production method can be selected

			let mut cur = mod_pow(witness.clone(), d.clone(), x.clone());

			// test witness ^ d
			if cur == One::one() || cur == neg_one {
				continue 'witness;  // not a witness
			}

			// test witness ^ 2d, witness ^ 4d, ..., witness ^ (2^(k-1) d)
			let mut counter: T = k.clone();
			while counter > literal(0) {  // turn this into a for loop, I dare you  (please?)
				counter = counter - literal(1);

				cur = cur.clone() * cur;
				cur = cur % x.clone();

				if cur == One::one() {
					return false;   // Composite
				}
				if cur == neg_one {
					continue 'witness;  // not a witness
				}
			}

			// We have a witness! (reached if no continues are encountered)
			return false;   // Composite
		}

		// We have no witnesses to the compositeness of x
		return true;   // Probable prime
	}
}



// Decomposes a number `x` into the form `2.pow(k) * d` where `d` is odd,
//  returning `(k,d)`.
fn decompose_pow2_odd<T>(x: T) -> (T, T)
 where T: Eq + Clone + Integer + Shr<usize, Output=T>,
{
	let mut remaining = x;
	let mut pow2: T = Zero::zero();
	while remaining.is_even() {
		pow2 = pow2 + literal(1);
		remaining = remaining >> literal(1);
	}
	return (pow2, remaining);
}

#[test]
fn test_miller_vs_sieve()
{
	// Test against prime sieve
	let limit = 20000usize;
	let expected = PrimeSieve::new(limit).sieve;

	for k in (0usize..limit) {
		if MillerRabinTester.is_prime(&k) != expected[k] {
			panic!("MillerRabinTester returned '{:?}' for {:?}", !expected[k], k);
		}
	}
}

#[test]
fn test_miller_mersenne()
{
	// Miller-Rabin should ALWAYS identify a prime as such, no matter how large.
	// (don't go too large, though, or the test will take a while :P)
	for p in vec![19usize, 31, 89, 107].into_iter() { // known p-values for some mersenne primes
		let one: BigUint = literal(1);
		let two: BigUint = literal(2);
		let mersenne: BigUint = num::pow(two,p) - one;
		assert!(MillerRabinTester.is_prime(&mersenne));
	}
}

#[bench]
fn bench_mersenne(b: &mut Bencher)
{
	let one: BigUint = literal(1);
	let two: BigUint = literal(2);
	let mersenne: BigUint = num::pow(two,107) - one;

	b.iter(|| {
		MillerRabinTester.is_prime(&mersenne)
	})
}

