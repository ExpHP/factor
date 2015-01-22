// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;

use std::collections::hash_map::{HashMap,Hasher};
use std::default::Default;
use std::ops::{Add,Sub,Mul,Div,Rem};
use std::cmp::{min,max};
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;

use num::{Zero,One,Integer};


/// TODO
#[derive(Show,Clone)]
pub struct Factorization<K>
 where K: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer,
{
	// powers are currently usize for consistency with the pow() functions
	//  in both std::num and the num crate.
	powers: HashMap<K, usize>,
}

// NOTE: ToPrimitive and FromPrimitive are unstable, and are currently included
//  in case I feel compelled to make a FactorTree using BigInt
impl<K> Factorization<K>
 where K: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer,
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

	// TODO: Make private
	pub fn from_factor(factor: K) -> Self
	{
		let mut map = HashMap::new();
		map.insert(factor, 1us);

		Factorization { powers: map }
	}

	// TODO: Make doctests less tedious-looking by using something like
	//        factorize(20), once such a function exists.

	/// The integer represented by this Factorization.
	///
	/// # Example:
	///
	/// ```
	/// use factor::Factorization;
	/// use std::default::Default;
	///
	/// let mut f: Factorization<u64> = Default::default();
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(5);
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

	/// The total number of positive integers which evenly divide
	///  the number represented by the Factorization (including 1
	///  and the number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::Factorization;
	/// use std::default::Default;
	///
	/// let mut f: Factorization<u64> = Default::default();
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(5);
	/// assert_eq!(f.count_divisors(), 6);
	/// ```
	pub fn count_divisors(self: &Self) -> K
	{
		let mut result: K = One::one();
		for (_,v) in self.powers.iter() {
			result = result * FromPrimitive::from_uint(v.clone() + 1).unwrap();
		};
		result
	}

	/// The sum of all positive integers which evenly divide
	///  the number represented by the Factorization (including 1
	///  and the number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::Factorization;
	/// use std::default::Default;
	///
	/// let mut f: Factorization<u64> = Default::default();
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(5);
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

	/// Euler totient function $\varphi(x)$, the total number of
	///  positive integers less than `x` that are coprime to `x`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::Factorization;
	/// use std::default::Default;
	///
	/// let mut f: Factorization<u64> = Default::default();
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(5);
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
	///
	/// # Example:
	///
	/// ```
	/// use factor::Factorization;
	/// use std::default::Default;
	///
	/// let mut f: Factorization<u64> = Default::default();
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(2);
	/// f = f * Factorization::from_factor(5);
	/// assert!(f.sqrt().is_none());  // 20: not square!
	/// 
	/// f = f * Factorization::from_factor(5);
	/// assert!(f.sqrt().is_some());  // 100: 10*10
	/// assert_eq!(f.sqrt().unwrap().product(), 10);
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

		let mut result: Self = Default::default();
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

}


impl<K> Default
for Factorization<K>
 where K: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer,
{
	/// Creates a Factorization representing 1.
	fn default() -> Self
	{ Factorization { powers: HashMap::new() } }
}


impl<K> Mul<Self> for Factorization<K>
 where K: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer,
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
