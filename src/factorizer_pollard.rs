// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

extern crate num;
extern crate test;

use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{Shr,Rem};
use std::mem::swap;

use std::cmp::min;

use num::{Zero, One, Integer};
use num::{FromPrimitive, ToPrimitive};
use rand::Rng;
use rand::weak_rng;
use rand::distributions::range::SampleRange;

use factorize;
use factorization::Factorization;
use factorizer::Factorizer;
use util::literal;
use util::gcd;
use util::mod_pow;

pub struct PollardBrentFactorizer;


impl<T> Factorizer<T>
for PollardBrentFactorizer
 where T: Eq + Clone + Zero + One + Integer + Shr<usize, Output=T> + Hash + SampleRange + FromPrimitive + ToPrimitive,
{
	/// Produce a single factor of `x`.  PollardBrentFactorizer is nondeterministic,
	/// and may sometimes fail to produce a non-trivial factor for composite `x`.
	fn get_factor(self: &Self, x: &T) -> T
	{
		// Adapted from https://comeoncodeon.wordpress.com/2010/09/18/pollard-rho-brent-integer-factorization/

		// In well-written functions, unused mutability is often indicative of a logic error.
		// So let it be known that this is NOT a well-written function:
		#![allow(unused_mut)]

		if x.is_even() { return literal(2); }
		if x.is_multiple_of(&literal(3)) { return literal(3); }
		if x < &literal(2) { return x.clone(); }

		let mut rng = weak_rng();
		let mut y: T = rng.gen_range(One::one(), x.clone()); // current value in the sequence:  y := y^2 + c (mod n)
		let mut c: T = rng.gen_range(One::one(), x.clone()); // parameter of y sequence
		let mut m: T = rng.gen_range(One::one(), x.clone()); // step size when multiplying crap together

		let mut g: T = One::one(); // contains the result
		let mut r: T = One::one(); // some kind of very coarse index
		let mut q: T = One::one(); // running product of `(z-y)` values

		let mut z: T = Zero::zero();      // Initial value of `y` for the current `r` value.
		let mut y_prev: T = Zero::zero(); // Initial value of `y` for the current `k` value.

		// Perform a coarse-grained search through the sequence of `y` values.
		while g == One::one() {
			z = y.clone();

			for _ in num::iter::range(Zero::zero(), r.clone()) {
				y = next_in_sequence(y, x.clone(), c.clone());
			}

			let mut k: T = Zero::zero();
			while k < r && g == One::one() {
				y_prev = y.clone();

				let niter = min(m.clone(), r.clone() - k.clone());

				// Multiply a bunch of (z-y) terms together (which may share factors with x)
				for _ in num::iter::range(Zero::zero(), niter) {
					y = next_in_sequence(y, x.clone(), c.clone());

					// Deviation from the source linked above, to support unsigned integers:
					//    abs(z-y) % x  --->  (x+z-y) % x
					// This is based on the notion that `gcd(+a % b, b) == gcd(-a % b, b)`,
					// so the absolute value isn't really necessary.
					q = q * (x.clone() + z.clone() - y.clone());
					q = q % x.clone();
				}

				g = gcd(x.clone(), q.clone());
				k = k + m.clone();
			}

			r = r * literal(2);
		} // end coarse-grained search

		// N.B. The following occurs when q == 0 (mod x).
		if &g == x {

			// Return to the beginning of this `k` step
			y = y_prev;

			loop {
				// Do a more fine grained search (computing the GCD every step)
				y = next_in_sequence(y, x.clone(), c.clone());
				g = gcd(x.clone(), x.clone() + z.clone() - y.clone()); // same deviation as noted above

				if g > One::one() { break; }
			}
		}

		// At this point, g is a nontrivial factor, or g == x.
		// In the latter case, g may still be composite (a "pseudoprime")
		return g;
	}
}

// computes (y**2 + c) % x
fn next_in_sequence<T>(y: T, x: T, c: T) -> T
	where T: Clone + Integer,
{
	let mut result = y.clone() * y;
	result = result % x.clone();
	result = result + c;
	result = result % x;
	return result;
}
