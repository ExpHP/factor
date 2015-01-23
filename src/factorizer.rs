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

//// XXX: Would reduce some of the pain in keeping impl type bounds consistent, but would also require the
////      trait to be explicitly implemented for all applicable types, which arguably sucks just as much.
//// Master trait of traits needed to implement Factorization<T>
//trait Factorable: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer { }

/// An object which can produce a single factor or a complete Factorization object
///  for any given number.
///
pub trait Factorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
{
	/// `get_factor(x)` for any non-composite `x` (0, 1, or primes) should return `x`.
	/// `get_factor(x)` for any composite `x` will produce an arbitrary but non-trivial factor of `x`.
	fn get_factor(self: &Self, x: &T) -> T;

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

/// A deterministic factorizer that is guaranteed to produce a factor for any number,
///  but which may not be the fastest thing out there.
//type SafeFactorizer<T> = TrialDivisionFactorizer<T>;

pub struct TrialDivisionFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T>
;


impl<T> Factorizer<T>
for TrialDivisionFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher>,
{
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
}

pub struct FermatFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T>
;

/// A relatively compact way to store the factorizations of all numbers up to a limit.
/// Only a single non-trivial factor is stored for each composite number, from which
///  the full decomposition can be gathered recursively.
///
/// TODO: Despite the name, not strictly a tree
#[derive(Clone, Show)]
pub struct FactorTree<T> {
	factors: Vec<T>,
}

impl<T> FactorTree<T>
 where T: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + Integer
{
//	// construction method?
//	fn new(n: T, factorizer: Fn(T) -> T) { }

	/// Returns a single factor (not necessarily prime) for a number.  The same FactorTree
	/// will always produce the same factor for the same x, though different trees may
	/// produce different factors, depending on their construction.
	///
	/// # Panics
	///
	/// Panics if the index exceeds the maximum number 
	///
	/// FIXME TODO DERP DOCZ
	///
	/// get_factor(x) for any non-composite `x` (0, 1, or primes) will return `x`.
	/// get_factor(x) for any composite `x` will produce an arbitrary but non-trivial factor of `x`.
	///
	/// FIXME TODO DERP DOCZ
	#[inline]
	fn get_factor(self: &Self, x: &T) -> T
	{
		self.factors[x.to_uint().unwrap()].clone()
	}
}

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
