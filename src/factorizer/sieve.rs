// Copyright 2015-2016 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use num::{Zero, One, Integer};
use ::TryFactor;
use ::Factored;
use util::MoreNumCast;

/// An object for efficiently factorizing many small numbers.
///
/// A `FactorSieve` is simply a precomputed array which stores the smallest
/// prime factor `p` for each integer `x`, up to some limit.
#[derive(Clone)]
pub struct FactorSieve {
	sieve: Vec<usize>,
}

impl FactorSieve
{
	/// Compute a factor sieve for numbers `< limit`.
	pub fn new<T>(limit: T) -> Self
	 where T: MoreNumCast + Integer,
	{
		assert!(limit >= T::zero());
		FactorSieve { sieve: factor_sieve_simple(limit.to_usize().unwrap()) }
	}
}
impl FactorSieve {
	pub fn into_vec(self) -> Vec<usize> { self.sieve }
	pub fn as_vec(&self) -> &Vec<usize> { &self.sieve }
	pub fn as_vec_mut(&mut self) -> &mut Vec<usize> { &mut self.sieve }
}

impl<T> TryFactor<T> for FactorSieve
 where T: Clone + Zero + One + Integer + MoreNumCast,
{
	fn try_factor_(&self, x: &T) -> T
	{ T::from_usize(self.sieve[x.to_usize().unwrap()]).unwrap() }

	fn factorize(&self, x: T) -> Factored<T>
	{ ::factorizer::helper::always_smallest_factorize(self, x) }
}

//------------------------------------------------------

fn factor_sieve_simple(limit: usize) -> Vec<usize>
{
	let mut sieve: Vec<_> = (0..limit).collect();

	// a very slight wheel optimization
	for x in (2..limit).step_by(2) { sieve[x] = 2; }

	for x in (3..).step_by(2) {
		// first multiple worth marking
		let first = x*x;
		if first > limit { break }

		// x is already marked
		if sieve[x] != x { continue }

		// x is prime; mark any multiples that aren't marked
		for m in (first..limit).step_by(2*x) {
			if sieve[m] == m { sieve[m] = x; }
		}
	}
	sieve
}

#[test]
fn test() {
	// a "definitely correct" sieve formed by trial division.
	// It is large enough to include values for which various certain
	// incorrectly-implemented optimizations could give the wrong prime
	// factor (e.g. a certain easy mistake results in giving 3 for 18,
	// 5 for 45, and 7 for 175; and the first incorrect value may be pushed
	// even further back if a large wheel optimization is used)
	let mut correct = vec![0,1];
	for i in 2..2000 {
		for k in 2..(i+1) {
			if i % k == 0 {
				correct.push(k);
				break;
			}
		}
	}
	for (i,x) in correct.clone().into_iter().enumerate() {
		if x == 3 {
			println!("{}", i);
		}
	}
	// a "definitely definitely correct" sieve to verify that one against,
	// up to a point
	assert_eq!(correct.len(), 2000);
	assert_eq!(&correct[..12], &vec![0,1,2,3,2,5,2,7,2,3,2,11][..]);

	// try various sizes, in particular:
	for len in vec![
		0, 1,   // Obvious edge cases (length 0, 1)
		2, 3,   // More subtle edge case (no multiples to cross off)
		9, 10, 11, // End in a prime square; a composite; and a prime
		correct.len(), // Catch erroneous values for larger x
	] {
		let actual = FactorSieve::new(len as usize).into_vec();
		assert_eq!(&actual[..], &correct[..len]);
	}
}
