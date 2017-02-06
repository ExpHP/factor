// Copyright (2015-2017) Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::default::Default;
use std::ops::Mul;
use std::cmp::{min,max};
use std::collections::BTreeMap;

use num;
use num::{One,FromPrimitive};
use num::Num;

use ::FactorExt;

/// An integer stored as its prime factorization.
///
/// A `Factored` represents a non-negative integer stored in the form
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
/// Negative numbers, zero, as well as negative powers (i.e. rationals) are not
/// permitted.
///
/// One caveat:  It is not feasible for `Factorization` to enforce that the prime
///  factors it is built from are actually prime.  This is the job of the
///  `Factorizer` used to produce it!  The behavior of any mathematical method
///  on a `Factorization` with a composite factor is **undefined**. (in the sense
///  of nonsensical results; not in the sense that the safety guidelines of Rust
///  may be violated)
#[derive(Eq,PartialEq,Debug,Clone)]
pub struct Factored<X>
{
	// powers are currently usize for consistency with the pow() function
	//  in the num crate.
	powers: BTreeMap<X, usize>,
}

// TODO: some of Factorization's methods need to be special cased for Zero

impl<X> Factored<X>
 where X: Ord,
{
	/// Set the power of a prime factor in the `Factorization`. `prime` is
	///  assumed to be prime. Note that while composite numbers
	///  may be accepted, the behavior of the resulting `Factorization` is
	///  ill-defined.
	// NOTE: No diagnostic is provided because it would require some sort
	//  of additional bound like Zero, One, or NumCast.
	pub fn set(&mut self, prime: X, exp: usize)
	{
		match exp {
			0 => self.powers.remove(&prime),
			_ => self.powers.insert(prime, exp),
		};
	}

	/// Get the power of a prime factor in the `Factorization`. `prime` is
	///  assumed to be prime. The value returned for non-prime
	///  keys is undefined.
	pub fn get(&self, prime: &X) -> usize
	{
		match self.powers.get(prime) {
			None       => 0,
			Some(&exp) => exp,
		}
	}


	/// Create a `Factored` from a single prime.
	///
	/// Do not use on a composite number (use `factor::factorize`)
	/// or one (use `One::one()`).
	pub fn from_prime(prime: X) -> Self
	{
		Factored { powers: ::std::iter::once((prime, 1)).collect() }
	}

	/// Recover the integer represented by the `Factored`.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.product(), 20);
	/// ```
	pub fn product(&self) -> X
	 where X: Clone + One,
	{ self.iter().map(clone2).factor_product() }

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
	pub fn divisor_count(&self) -> X
	 where X: Clone + FromPrimitive + One,
	{ self.iter().map(clone2).divisor_count() }

	/// Compute the sum of all positive integers which evenly divide the number
	///  represented by the `Factored` (including 1 and the number itself).
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let f = factorize(20u64);
	/// assert_eq!(f.divisor_sum(), 42);
	/// ```
	pub fn divisor_sum(&self) -> X
	 where X: Clone + One + Num,
	{ self.iter().map(clone2).divisor_sum() }

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
	pub fn totient(&self) -> X
	 where X: Clone + One + Num,
	{ self.iter().map(clone2).totient() }
}

impl<X> Factored<X>
 where X: Ord,
{
	/// Get the square root of a perfect square.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let mut f = factorize(20u64);
	/// assert!(f.clone().sqrt().is_none());    // 20 is not square!
	///
	/// f = f * factorize(5);
	/// assert_eq!(f.sqrt().unwrap().product(), 10);  // sqrt(100)
	/// ```
	pub fn sqrt(self) -> Option<Self>
	{ self.root(2) }

	/// Get the nth root of a perfect nth power.
	///
	/// # Example:
	///
	/// ```
	/// use factor::factorize;
	///
	/// let mut f = factorize(20u64);
	/// assert!(f.clone().root(3).is_none());    // 20 is not a cube!
	///
	/// f = f * factorize(50);
	/// assert_eq!(f.root(3).unwrap().product(), 10);  // cbrt(1000)
	/// ```
	pub fn root(mut self, n: usize) -> Option<Self>
	{
		for (_, v) in self.powers.iter_mut() {
			if *v % n == 1 {
				return None; // Not a root!  Abort!  Abort!
			};

			*v = *v / n;
		};
		Some(self)
	}

	// TODO: Docs & test
	pub fn pow(mut self, pwr: usize) -> Self
	 where X: Clone,
	{
		for (_, v) in self.powers.iter_mut() {
			*v = *v * pwr;
		};
		self
	}

/*
	// TODO: Docs & test
	pub fn gcd(mut self, other: mut Self) -> Self
	 where X: Clone,
	{
		// Iterate through shorter list
		if other.powers.len() < self.powers.len() {
			return self.gcd(other);
		};

		let mut result: Self = One::one();
		for (k,v) in other {
			self.set(k.clone(), min(self.get(&k), v));
		};
		self
	}

	// TODO: Docs & test
	pub fn lcm(&self, other: &Self) -> Self
	 where X: Clone,
	{
		let mut result = self.clone();
		for (k,v) in self.powers.clone().into_iter() {
			self.set(k.clone(), max(self.get(&k), v));
		};
		self
	}
*/

	/// Iterate over unique prime factors.
	///
	/// Item type is `&X`.
	pub fn primes(&self) -> Primes<X> { self.powers.keys() }

	/// Iterate over `(prime, power)` pairs.
	///
	/// Item type is `(&X, usize)`.
	pub fn iter(&self) -> Iter<X> { self.powers.iter() }

	/// Consume to iterate over the `(prime, power)` pairs.
	///
	/// Item type is `(X, usize)`.
	pub fn into_iter(self) -> IntoIter<X> { self.powers.into_iter() }

	/// Immutably borrow the underlying map.
	pub fn as_btree_map(&self) -> &BTreeMap<X, usize> { &self.powers }
	/// Mutably borrow the underlying map.
	///
	/// This is provided on the grounds of "eh, why the hell not."  After all,
	/// it's not exactly like `Factored` was doing a very good job at protecting
	/// any of its invariants to begin with.
	pub fn as_mut_btree_map(&mut self) -> &mut BTreeMap<X, usize> { &mut self.powers }
	/// Consume to obtain the underlying map.
	pub fn into_btree_map(self) -> BTreeMap<X, usize> { self.powers }

	/// Construct from a map of unique `(prime, power)` pairs.
	pub fn from_btree_map(map: BTreeMap<X, usize>) -> Self {
		Factored { powers: map }
	}

	/// Construct from an iterator of unique `(prime, power)` pairs.
	pub fn from_iter<I:IntoIterator<Item=(X,usize)>>(iter: I) -> Self {
		Factored::from_btree_map(iter.into_iter().collect())
	}

	// TODO: Reasonable (?) methods to implement
	// pub fn is_multiple_of(&self, other: &Self) -> bool
}

impl<X> One
for Factored<X>
 where X: Ord + Clone,
{
	/// Create a `Factored` representing `1`.
	fn one() -> Self { Self::from_btree_map(BTreeMap::new()) }
}

impl<X> Default
for Factored<X>
 where X: Ord + Clone,
{
	/// Create a `Factored` representing `1`.
	fn default() -> Self { One::one() }
}

impl<X> Mul<Factored<X>>
for Factored<X>
 where X: Ord + Clone,
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

//------------------------------------

pub type Iter<'a,T> = ::std::collections::btree_map::Iter<'a,T,usize>;
pub type IntoIter<T> = ::std::collections::btree_map::IntoIter<T,usize>;
pub type Primes<'a,T> = ::std::collections::btree_map::Keys<'a,T,usize>;

#[inline]
fn clone2<A,B>((a,b): (&A,&B)) -> (A,B)
 where A: Clone, B: Clone { (a.clone(), b.clone()) }

