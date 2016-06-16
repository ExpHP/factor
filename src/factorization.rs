// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::default::Default;
use std::ops::{Add,Sub,Mul,Div,Rem};
use std::cmp::{min,max};
use std::collections::HashMap;
use std::collections::hash_map;
use std::hash::Hash;

use num;
use num::{One,FromPrimitive};
use num::Num;

// TODO After all this time I still don't like the name; it's too easy to
//  confuse with `Factorizer`... yet I can't think of a better alternative.
//  `PrimeFactorMap`?... bleh.
/// An integer stored as its prime factorization.
///
/// A `Factorization` represents a non-negative integer stored in the form
///
/// ```text
/// pow(p1, n1) * pow(p2, n2) * ... * pow(pN, nN)
/// ```
///
/// where `p1, p2, ...` are prime numbers (or zero) and `n1, n2, ...` are positive
/// integers. Some properties are easier to compute in this form, such as the euler
/// totient, or a sum of divisors. Also of note is that because the prime factorization
/// of any given number is unique, any two `Factorization`s constructed from the same
/// `x` are equal.
///
/// A special case is permitted to allow zero, defined as `{0: 1}`. Zero to any greater
/// power, or in combination with other factors is not permitted. Also worthy of note
/// is that the factorization of `1` is an empty map.
/// Negative numbers are not permitted as most of the special mathematical properties
/// are not defined for them.
///
/// One caveat:  It is not feasible for `Factorization` to enforce that the prime
///  factors it is built from are actually prime.  This is the job of the
///  `Factorizer` used to produce it!  The behavior of any mathematical method
///  on a `Factorization` with a composite factor is **undefined**. (in the sense
///  of nonsensical results; not in the sense that the safety guidelines of Rust
///  may be violated)
#[derive(Eq,PartialEq,Debug,Clone)]
pub struct Factorization<K>
 where K: Eq + Hash,
{
	// powers are currently usize for consistency with the pow() functions
	//  in both std::num and the num crate.
	powers: HashMap<K, usize>,
}

// TODO: some of Factorization's methods need to be special cased for Zero

impl<K> Factorization<K>
 where K: Eq + Hash,
{

	/// Set the power of a prime factor in the `Factorization`. `key` is
	///  assumed to be prime (or zero). Note that while composite numbers
	///  may be accepted, the behavior of the resulting `Factorization` is
	///  ill-defined.
	pub fn set(&mut self, key: K, pwr: usize)
	{
		match pwr {
			0 => self.powers.remove(&key),
			_ => self.powers.insert(key, pwr),
		};
	}

	/// Get the power of a prime factor in the `Factorization`. `key` is
	///  assumed to be prime (or zero). The value returned for non-prime
	///  keys is undefined.
	pub fn get(&self, key: &K) -> usize
	{
		match self.powers.get(key) {
			None      => 0,
			Some(pwr) => *pwr,
		}
	}

	/// Create a `Factorization` for the given number, which is assumed to be
	///  prime (or zero).  Note that while composite numbers may be accepted,
	///  the behavior of the resulting `Factorization` is ill-defined.
	pub fn from_prime(prime: K) -> Self
	{
		let mut map = HashMap::new();
		map.insert(prime, 1usize);

		Factorization { powers: map }
	}

	// TODO: Make doctests less tedious-looking by using something like
	//        factorize(20), once such a function exists.
}

#[inline]
fn clone2<A,B>((a,b): (&A,&B)) -> (A,B)
 where A: Clone, B: Clone { (a.clone(), b.clone()) }

impl<K> Factorization<K>
 where K: Eq + Hash,
{
	/// Recover the integer represented by the `Factorization`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.product(), 20);
	/// ```
	pub fn product(&self) -> K
	 where K: Clone + One,
	{ self.powers.iter().map(clone2).factor_product() }

	/// Compute the total number of positive integers which evenly divide
	///  a number from its prime decomposition.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.divisor_count(), 6);
	/// ```
	pub fn divisor_count(&self) -> K
	 where K: Clone + FromPrimitive + One,
	{ self.powers.iter().map(clone2).divisor_count() }

	/// Compute the sum of all positive integers which evenly divide the number
	///  represented by the `Factorization` (including 1 and the number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.divisor_sum(), 42);
	/// ```
	pub fn divisor_sum(&self) -> K
	 where K: Clone + One + Num,
	{ self.powers.iter().map(clone2).divisor_sum() }

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
	pub fn totient(&self) -> K
	 where K: Clone + One + Num,
	{ self.powers.iter().map(clone2).totient() }
}

// FIXME move
/// Prime-factorization-based functions implemented on arbitrary
///  iterables of `(prime, power)` pairs.
///
/// There are certain preconditions for using these methods:
///
///   * In each pair, `prime` must actually be a prime;
///     It is forbidden to be zero, one, negative, or composite.
///   * `power` must not be zero.
///   * The primes must be sorted, with each appearing only once.
///
/// A diagnostic for violations of these conditions (with the exception
///  of non-compositeness) may be provided in debug builds.
pub trait FactorExt<X>: Iterator<Item=(X,usize)>
{

	/// Recover an integer from its prime decomposition.
	///
	/// # Example:
	///
	/// ```
	/// use factor::{factorize, FactorExt};
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.into_iter().factor_product(), 20);
	/// ```
	fn factor_product(self) -> X
	 where X: Clone + One,
	       Self: Sized,
	{
		// FIXME why are we not using fold
		// product( pow(p[i], n[i]) )
		let mut result: X = One::one();
		for (k,v) in self {
			result = result * num::pow(k.clone(), v.clone());
		};
		result
	}

	/// Compute the total number of positive integers which evenly divide
	///  a number from its prime decomposition.
	///
	/// # Example:
	///
	/// ```
	/// use factor::{factorize, FactorExt};
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.into_iter().divisor_count(), 6);
	/// ```
	fn divisor_count(self) -> X
	 where X: Clone + FromPrimitive + One,
	       Self: Sized,
	{
		// FIXME why are we not using fold
		// product( n[i]+1 );  we're counting the possible ways to select a power
		//                      from 0 to n[i] (inclusive) for each prime p[i]
		let mut result: X = One::one();
		for (_,v) in self {
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
	/// use factor::{factorize, FactorExt};
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.into_iter().divisor_sum(), 42);
	/// ```
	fn divisor_sum(self) -> X
	 where X: Clone + One + Num,
	       Self: Sized,
	{
		// FIXME why are we not using fold
		// A geometric series for each factor
		let mut result: X = One::one();
		for (k,v) in self {
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
	/// use factor::{factorize, FactorExt};
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.into_iter().totient(), 8);
	/// ```
	fn totient(self) -> X
	 where X: Clone + One + Num,
	       Self: Sized,
	{
		// product( (p[i] - 1) * pow(p[i], n[i]-1) )
		let mut result: X = One::one();
		for (k,v) in self {
			result = result * num::pow(k.clone(), v.clone() - 1);
			result = result * (k.clone() - One::one());
		};
		return result;
	}

	// FIXME: sigma_k for other k?

	// FIXME for now I see no issue with just returning a ::std::vec::IntoIter.  Do that.
	// TODO: pending unboxed abstract types 
	// fn iter_divisors(&self) -> Iterator<K>
}

impl<I,X> FactorExt<X> for I
 where I: Iterator<Item=(X,usize)>,
{ }

impl<K> Factorization<K>
 where K: Eq + Hash,
{
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
	pub fn sqrt(&self) -> Option<Self>
	 where K: Clone,
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
	pub fn pow(&self, pwr: usize) -> Self
	 where K: Clone,
	{
		let mut result = self.clone();
		for (_,v) in result.powers.iter_mut() {
			*v = *v * pwr;
		};
		result
	}

	// TODO: Docs & test
	pub fn gcd(&self, other: &Self) -> Self
	 where K: Clone,
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
	pub fn lcm(&self, other: &Self) -> Self
	 where K: Clone,
	{
		let mut result = self.clone();
		for (k,v) in other.powers.clone().into_iter() {
			result.set(k.clone(), max(self.get(&k), v));
		};
		result
	}

	/// Iterate over unique prime factors.
	///
	/// Item type is `&K`.
	pub fn primes(&self) -> hash_map::Keys<K, usize> { self.powers.keys() }

	/// Iterate over `(prime, power)` pairs.
	///
	/// Item type is `(&K, &usize)`.
	pub fn iter(&self) -> hash_map::Iter<K, usize> { self.powers.iter() }

	/// Consume to iterate over the `(prime, power)` pairs.
	///
	/// Item type is `(K, usize)`.
	pub fn into_iter(self) -> hash_map::IntoIter<K, usize> { self.powers.into_iter() }

	/// Immutably borrow the underlying `HashMap`.
	pub fn as_hash_map(&self) -> &HashMap<K, usize> { &self.powers }
	/// Mutably borrow the underlying `HashMap`.
	pub fn as_mut_hash_map(&mut self) -> &mut HashMap<K, usize> { &mut self.powers }
	/// Consume to obtain the underlying `HashMap`.
	pub fn into_hash_map(self) -> HashMap<K, usize> { self.powers }

	/// Construct from an iterator of `(prime, power)` pairs.
	///
	/// The primes must be unique, else behavior is unspecified.
	///
	/// Each `prime` must be prime (or zero), or else the resulting `Factorization`
	///  will have ill-defined behavior.
	pub fn from_iter<I:IntoIterator<Item=(K,usize)>>(iter: I) -> Self {
		Factorization::from_hash_map(iter.into_iter().collect())
	}

	/// Construct from a `HashMap` of `prime => power`.
	///
	/// Each `prime` must be prime (or zero), or else the resulting `Factorization`
	///  will have ill-defined behavior.
	pub fn from_hash_map(powers: HashMap<K, usize>) -> Self {
		Factorization { powers: powers }
	}

	// TODO: Reasonable (?) methods to implement
	// pub fn is_multiple_of(&self, other: &Self) -> bool
}


impl<K> One
for Factorization<K>
 where K: Eq + Clone + Hash,
{
	/// Creates a Factorization representing `1`.
	fn one() -> Self
	{ Factorization { powers: HashMap::new() } }
}


impl<K> Default
for Factorization<K>
 where K: Eq + Clone + Hash,
{
	/// Creates a Factorization representing `1`.
	fn default() -> Self
	{ One::one() }
}


impl<K> Mul<Factorization<K>>
for Factorization<K>
 where K: Eq + Clone + Hash,
{
	type Output = Self;
	fn mul(self, other: Self) -> Self
	{
		let mut result = self;
		for (k,v) in other.powers.into_iter() {
			let current = result.get(&k); // owned proxy
			result.set(k.clone(), current + v);
		}
		result
	}
}
