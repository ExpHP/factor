// Copyright 2015-2016 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use num::{Zero, One, Integer, NumCast};
use ::std::hash::Hash;
use ::Factorizer;
use ::Factorization;

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
	 where T: NumCast + Integer,
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

impl<T> Factorizer<T> for FactorSieve
 where T: Eq + Clone + Zero + One + Integer + Hash + NumCast,
{
	fn get_factor(&self, x: &T) -> T
	{ T::from(self.sieve[x.to_usize().unwrap()]).unwrap() }

	// FIXME this is just a specialized implementation for factorizers which
	// always return the smallest factor. It should be pulled out and reused.
	fn factorize(&self, x: T) -> Factorization<T>
	{
		assert!(x >= T::zero());
		// special cases to make life easier
		if (x == T::zero()) { return Factorization::from_iter(vec![(x, 1)]); }
		if (x == T::one()) { return Factorization::from_iter(vec![]); }

		// the rest of this is basically trying to implement something like the
		// following in iterator-speak (but there's no grouping iterator in std,
		// and it is still painful to write an iterator in rust)
		//    iter.group_by(|p| p)
		//        .map(|(p, group)| (p, group.count()))
		// where iter is an imaginary iterator giving one prime at a time (with repeats)

		// begin first group
		let mut prev = self.get_factor(&x);
		let mut x = x/prev.clone();
		let mut count = 1;
		assert!(prev != T::one());

		let mut out = vec![];
		while x != T::one() {
			let p = self.get_factor(&x);
			if p != prev {
				// start new group
				out.push((prev, count));
				prev = p.clone();
				count = 0;
			}
			count += 1;
			x = x / p;
		}
		out.push((prev, count)); // final group
		Factorization::from_iter(out)
	}
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

