// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;
extern crate test;

use std::collections::hash_map::{HashMap,Hasher};
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;
use std::ops::{Shr,Rem};

use num::{Zero, One, Integer};

use util::isqrt;
use factorization::Factorization;


// Type-Synonyms, for semantic purposes.
// (alas, re-exports don't appear to take docstrings)

// /// A deterministic `Factorizer` that is guaranteed to work on any number, but may be fast.
pub use TrialDivisionFactorizer as SafeFactorizer;
// /// The `Factorizer` used by the `factor::factorize` method.
pub use TrialDivisionFactorizer as DefaultFactorizer;


//// XXX: Would reduce some of the pain in keeping impl type bounds consistent, but would also require the
////      trait to be explicitly implemented for all applicable types, which arguably sucks just as much.
//// Master trait of traits needed to implement Factorization<T>
//trait Factorable: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer { }


/// An interface for factorizing positive integers.
pub trait Factorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
{
	// TODO: Perhaps more sensible to return `1` for non-composite (so that division by the result
	//       is safe without needing to check for zero)
	/// Produces a single (not necessarily prime) factor of a number `x`.  Some implementations of
	///  `Factorizer` are deterministic and always produce the same factor for the same `x`, while
	///  others may have an element of randomness.
	///
	/// However, all Factorizers are expected to meet the following guarantees:
	///
	/// * `get_factor(x)` for any non-composite `x` (0, 1, or primes) returns `x`.
	/// * `get_factor(x)` for any composite `x` produces an arbitrary but non-trivial factor of `x`.
	/// 
	/// Keep in mind that, in addition to the value returned, another factor can be obtained by
	///  dividing `x` by the value.
	fn get_factor(self: &Self, x: &T) -> T;

	/// Builds a complete prime factorization of a number.  A default implementation is provided
	///  which calls `get_factor()` recursively on the factors produced.
	fn factorize(self: &Self, x: T) -> Factorization<T>
	{
		let a = self.get_factor(&x);

		// Non-composite point to themselves
		if a == x {
			return Factorization::from_factor(a);

		// Composite numbers
		} else {
			let b = x.clone() / a.clone();
			return self.factorize(a) * self.factorize(b);
		};
	}
}

/// Factors numbers using trial division.
///
/// Trial division is one of the most inefficient, yet most easily understood methods for
///  factorizing numbers.  `TrialDivisionFactorizer` tries dividing a number by successively
///  larger numbers until it encounters one that leaves a remainder of zero.
///
/// Despite its primitive nature, it can well outperform many of the more sophisticated methods
///  when factoring small numbers (TODO: of what magnitude?). However, it has trouble on numbers
///  with large prime factors.
pub struct TrialDivisionFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher>
;


impl<T> Factorizer<T>
for TrialDivisionFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher>,
{
	/// Produce a single factor of `x`.  TrialDivisionFactorizer is deterministic,
	///  and will always produce the smallest non-trivial factor of any composite number.
	///  Thus, the number it returns is also always prime.
	///
	/// The runtime scales linearly with the size of the smallest factor of `x`.
	fn get_factor(self: &Self, x: &T) -> T
	{
		// If you're reading this:
		//       I'm sorry.

		if x.is_zero() { return Zero::zero() };

		if x.is_even() { return FromPrimitive::from_uint(2).unwrap(); }

		let start: T = FromPrimitive::from_uint(3).unwrap();
		let stop:  T = isqrt(x.clone()) + FromPrimitive::from_uint(2).unwrap();
		let step:  T = FromPrimitive::from_uint(2).unwrap();
		
		//for odd in range_step(start, stop, step)  // needs trait Int
		let mut odd = start;
		while odd < stop {
			if x.is_multiple_of(&odd) {
				return odd;
			}
			odd = odd + step.clone();
		};

		// x is prime
		return x.clone();
	}

	// TODO: The default recursive algorithm for factorize() is poorly suited
	//       to this Factorizer, and should be overriden.
}

/// TODO
pub struct FermatFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
;

/// TODO
pub struct DixonFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
;

/// TODO
pub struct GeneralFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
;

/// Factors numbers by using results cached from another `Factorizer`.
///
/// A `Factorizer` which stores factors produced by another `Factorizer` for quick lookup.
/// Only a single non-trivial factor is stored for each composite number, from which
///  the full decomposition can be gathered recursively.
#[derive(Clone, Show)]
pub struct FactorStore<T> {
	// I wanted to call it FactorTree but it's really a DAG.  x_x

	factors: Vec<T>,
}

impl<T> FactorStore<T>
 where T: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + Integer
{
//	// construction method?
//	fn new(n: T, factorizer: Fn(T) -> T) { }
}

impl<T> Factorizer<T>
for FactorStore<T>
 where T: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + Zero + One + Integer + Shr<usize, Output=T>,
{
	/// Produces the factor stored for `x`.
	///
	/// # Panics
	///
	/// May panic. (TODO: When?)
	#[inline]
	fn get_factor(self: &Self, x: &T) -> T
	{
		self.factors[x.to_uint().unwrap()].clone()
	}
}


// Tests

#[cfg(test)]
mod tests {
	extern crate test;
	use super::*;
	use test::Bencher;
	use std::collections::hash_map::{HashMap,Hasher};
	use std::num::{ToPrimitive,FromPrimitive}; // and regret it
	use std::hash::Hash;
	use std::fmt::Show;

	use num::{BigUint, BigInt};
	use num::{Zero, One, Integer};

	//  A simple test to factorize 242 as an arbitrary data type using an arbitrary factorizer.
	fn test_242<T, U>(factorizer: U)
	 where T: Eq + Clone + Show + Hash<Hasher> + ToPrimitive + FromPrimitive + Integer,
	       U: Factorizer<T>,
	{
		// 242 = 2 * 11 * 11
		let x_t: T = FromPrimitive::from_uint(242).unwrap();

		let factors = factorizer.factorize(x_t.clone());

		// vector of expected powers
		let mut expected = vec![0us; 15];
		expected[2]  = 1;
		expected[11] = 2;

		// Check factorization against expected		
		for (k,v) in expected.iter().enumerate() {
			let k_t = FromPrimitive::from_uint(k).unwrap();  // cast k to T
			let v_t = FromPrimitive::from_uint(*v).unwrap(); // cast v to whatever the heck it is today
			assert_eq!(factors.get(&k_t), v_t);
		}
	}

	// Mix and match various factorizers and data types --- mostly to confirm that they compile
	//  (e.g. trait bounds are met, and whatnot)
	#[test]
	fn test_mix_and_match() {
		// differently sized integers
		test_242::<u32,_>(TrialDivisionFactorizer);
		test_242::<u64,_>(TrialDivisionFactorizer);
		test_242::<usize,_>(TrialDivisionFactorizer);

		// a signed type
		test_242::<isize,_>(TrialDivisionFactorizer);
	
		// my pain and suffering
		test_242::<BigUint,_>(TrialDivisionFactorizer);
		test_242::<BigInt,_>(TrialDivisionFactorizer);

		// TODO: Test other factorizers here when they exist
	}
}