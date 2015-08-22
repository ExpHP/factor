// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;

use std::default::Default;
use std::ops::{Add,Sub,Mul,Div,Rem};
use std::cmp::{min,max};
use std::collections::HashMap;
use std::hash::Hash;

use num::{Zero,One,Integer,FromPrimitive};


/// TODO
#[derive(Eq,PartialEq,Debug,Clone)]
pub struct Factorization<K>
 where K: Eq + Clone + Integer + Hash,
{
	// powers are currently usize for consistency with the pow() functions
	//  in both std::num and the num crate.
	powers: HashMap<K, usize>,
}

// TODO: some of Factorization's methods need to be special cased for Zero

// TODO: While Factorization is intended to only store prime numbers (and zero)
//        as factors, it will currently accept anything.
//
//       Some alternatives under consideration:
//         * a Prime data type (possibly one form of a "known primality" enum
//             type) for keys
//         * allow non-prime keys.  The issue with this is that many of the
//           math methods (gcd, totient, etc.) are only worth implementing
//           for the all-prime case, and the user would have to be trusted to
//           use them responsibly.

impl<K> Factorization<K>
 where K: Eq + Clone + Integer + Hash,
{

	pub fn set(self: &mut Self, key: K, pwr: usize) {
		match pwr {
			0 => self.powers.remove(&key),
			_ => self.powers.insert(key, pwr),
		};
	}

	pub fn get(self: &Self, key: &K) -> usize {
		match self.powers.get(key) {
			None      => 0,
			Some(pwr) => *pwr,
		}
	}

	pub fn from_factor(factor: K) -> Self
	{
		let mut map = HashMap::new();
		map.insert(factor, 1usize);

		Factorization { powers: map }
	}

	// TODO: Make doctests less tedious-looking by using something like
	//        factorize(20), once such a function exists.

	/// Compute the integer represented by the `Factorization`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.product(), 20);
	/// ```
	pub fn product(self: &Self) -> K
	{
		let mut result: K = One::one();
		for (k,v) in self.powers.iter() {
			result = result * num::pow(k.clone(), v.clone());
		};
		result
	}

	/// Compute the total number of positive integers which evenly divide
	///  the number represented by the `Factorization` (including 1 and the
	///  number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.count_divisors(), 6);
	/// ```
	pub fn count_divisors(self: &Self) -> K
	 where K: FromPrimitive,
	{
		let mut result: K = One::one();
		for (_,v) in self.powers.iter() {
			result = result * FromPrimitive::from_usize(v.clone() + 1).unwrap();
		};
		result
	}

	/// Compute the sum of all positive integers which evenly divide the number
	///  represented by the `Factorization` (including 1 and the number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.sum_divisors(), 42);
	/// ```
	pub fn sum_divisors(self: &Self) -> K
	{
		// A geometric series for each factor
		let mut result: K = One::one();
		for (k,v) in self.powers.iter() {
			result = result * (num::pow(k.clone(), v.clone() + 1) - One::one());
			result = result / (k.clone() - One::one());
		};
		result
	}

	/// Compute the Euler totient `phi(x)`, the total number of
	///  positive integers less than `x` that are coprime to `x`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.totient(), 8);
	/// ```
	pub fn totient(self: &Self) -> K
	{
		let mut result: K = One::one();
		for (k,v) in self.powers.iter() {
			result = result * num::pow(k.clone(), v.clone() - 1);
			result = result * (k.clone() - One::one());
		};
		return result;
	}

	// TODO: pending unboxed abstract types 
	// fn iter_divisors(self: &Self) -> Iterator<K>

	// TODO: Docs
	// Returns None if any factors have odd powers
	/// Get the square root.
	///
	/// If the `Factorization` represents a perfect square, returns
	///  a `Some(Factorization)` representing the square root.
	///  Otherwise, returns `None`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let mut f = factorize(20u64);
	/// assert!(f.sqrt().is_none());    // 20 is not square!
	/// 
	/// f = f * factorize(5);
	/// assert_eq!(f.sqrt().unwrap().product(), 10);  // sqrt(100)
	/// ```
	pub fn sqrt(self: &Self) -> Option<Self>
	{
		let mut result = self.clone();
		for (_,v) in result.powers.iter_mut() {

			if *v % 2 == 1 {
				return None; // Odd power!  Abort!  Abort!
			};

			*v = *v / 2;
		};
		Some(result)
	}

	// TODO: Docs & test
	pub fn pow(self: &Self, pwr: usize) -> Self
	{
		let mut result = self.clone();
		for (_,v) in result.powers.iter_mut() {
			*v = *v * pwr;
		};
		result
	}

	// TODO: Docs & test
	pub fn gcd(self: &Self, other: &Self) -> Self
	{
		// Iterate through shorter list
		if other.powers.len() < self.powers.len() {
			return self.gcd(other);
		};

		let mut result: Self = One::one();
		for (k,v) in other.powers.clone().into_iter() {
			result.set(k.clone(), min(self.get(&k), v));
		};
		result
	}

	// TODO: Docs & test
	pub fn lcm(self: &Self, other: &Self) -> Self
	{
		let mut result = self.clone();
		for (k,v) in other.powers.clone().into_iter() {
			result.set(k.clone(), max(self.get(&k), v));
		};
		result
	}

	// TODO: Reasonable (?) methods to implement
	// pub fn is_multiple_of(self: &Self, other: &Self) -> bool
}


impl<K> One
for Factorization<K>
 where K: Eq + Clone + Hash + Integer,
{
	/// Creates a Factorization representing `1`.
	fn one() -> Self
	{ Factorization { powers: HashMap::new() } }
}


impl<K> Default
for Factorization<K>
 where K: Eq + Clone + Hash + Integer,
{
	/// Creates a Factorization representing `1`.
	fn default() -> Self
	{ One::one() }
}


impl<K> Mul<Factorization<K>>
for Factorization<K>
 where K: Eq + Clone + Hash + Integer,
{
	type Output = Self;
	fn mul(self: Self, other: Self) -> Self
	{
		let mut result = self;
		for (k,v) in other.powers.into_iter() {
			let current = result.get(&k); // owned proxy
			result.set(k.clone(), current + v);
		}
		result
	}
}
