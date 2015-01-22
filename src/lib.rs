// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

//! Documentation is currently under construction.

// TODO: Make documentation not be "currently under construction." :P

#![allow(unused_imports)]

extern crate num;

use std::collections::hash_map::{HashMap,Hasher};
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;

pub use primes::{PrimeSieve,PrimeTester};
pub use factorization::Factorization;

mod primes;
mod factorization;

/// A relatively compact way to store the factorizations of all numbers up to a limit.
/// Only a single non-trivial factor is stored for each composite number, from which
///  the full decomposition can be gathered recursively.
///
/// TODO: Despite the name, not strictly a tree
#[derive(Clone, Show)]
pub struct FactorTree<T> {
	factors: Vec<T>,
}

// TODO: Test on bigint

//// XXX: Would reduce some of the pain in keeping impl type bounds consistent, but would also require the
////      trait to be explicitly implemented for all applicable types, which arguably sucks just as much.
//// Master trait of traits needed to implement Factorization<T>
//trait Factorable: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer { }

impl<T> FactorTree<T>
 where T: Eq + Clone + Hash<Hasher> + ToPrimitive + FromPrimitive + num::Integer
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

	/// Builds the complete prime factorization of a given number by traversing the structure as
	/// a tree.
	/// // TODO
	fn factorize(self: &Self, x: T) -> Factorization<T>
	{
		let a = self.get_factor(&x);

		// Non-composite numbers are leaves, pointing to themselves
		if a == x {
			// TODO: Construct new factorization
			return Factorization::from_factor(a);
		// Composite numbers
		} else {
			let b = x.clone() / a.clone();
			return self.factorize(a.clone()) * self.factorize(b.clone());
		};
	}

}

#[test]
fn it_works() {

}
