// Copyright 2015-2016 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

//! Extension trait for factorization math on iterators.

use num;
use num::{One,FromPrimitive};
use num::Num;

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
	/// use factor::prelude::*;
	///
	/// let f = factor::factorize(20u64);
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
	/// use factor::prelude::*;
	///
	/// let f = factor::factorize(20u64);
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
	/// use factor::prelude::*;
	///
	/// let f = factor::factorize(20u64);
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
	/// use factor::prelude::*;
	///
	/// let f = factor::factorize(20u64);
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
