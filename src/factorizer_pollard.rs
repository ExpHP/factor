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
use std::collections::{Bitv,BitvSet};
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;
use std::ops::{Shr,Rem};
use std::rand::Rng;
use std::rand::weak_rng;
use std::rand::distributions::range::SampleRange;
use std::mem::swap;
use std::num::Int;

use std::cmp::min;

use num::{Zero, One, Integer};
use num::traits::Signed;

use factorize;
use factorization::Factorization;
use factorizer::Factorizer;
use util::literal;
use util::gcd;
use util::mod_pow;

pub struct PollardBrentFactorizer<T>;

impl<T> Factorizer<T>
for PollardBrentFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher> + SampleRange + Int + Signed,
{
	/// Produce a single factor of `x`.  PollardBrentFactorizer is nondeterministic,
	///  and will always produce the smallest non-trivial factor of any composite number.
	///  Thus, the number it returns is also always prime.
	///
	/// The runtime scales linearly with the size of the smallest factor of `x`.
	fn get_factor(self: &Self, x: &T) -> T
	{
		// Adapted from https://comeoncodeon.wordpress.com/2010/09/18/pollard-rho-brent-integer-factorization/
		if x.is_even() { return literal(2); }
		if x.is_multiple_of(&literal(3)) { return literal(3); }
		if x < &literal(2) { return x.clone(); }

		let mut rng = weak_rng();
		let mut y: T = rng.gen_range(One::one(), x.clone());
		let mut c: T = rng.gen_range(One::one(), x.clone());
		let mut m: T = rng.gen_range(One::one(), x.clone());

		let mut g: T = One::one();
		let mut r: T = One::one();
		let mut q: T = One::one();

		let mut y_prev: T = Zero::zero();

		let mut z: T = Zero::zero();

		while g == One::one() {
			let mut z = y.clone();

			for _ in num::iter::range(Zero::zero(), r) {
				y = next_in_sequence(y, x.clone(), c.clone());
			}

			// NOTE:  The purpose of this appears to be to avoid computing the
			//        gcd at every term in the sequence, and instead only computing
			//        it every m'th term.
			let mut k: T = Zero::zero();
			while k < r && g == One::one() {
				y_prev = y;

				let niter = min(m.clone(), r.clone() - k.clone());
				for _ in num::iter::range(Zero::zero(), niter) {
					y = next_in_sequence(y, x.clone(), c.clone());

					// FIXME: Not only will this fail for signed types, but it appears to be nonsense.
					//        Since when do absolute value and modular arithmetic EVER go together?
					q = q * (z.clone() - y.clone()).abs();
					q = q % x.clone();
				}

				g = gcd(x.clone(), q.clone());
				k = k + m.clone();
			}

			r = r * literal(2);
		}

		// FIXME:  The following block of code (which is in the original implementation from which
		//          this is adapted; see above) clearly never runs, as g = gcd(q, x), and the value
		//          of q is reduced mod x.
		//         Find out the original author's intent!
		if &g == x {
			y = y_prev;
			loop {
				y = next_in_sequence(y, x.clone(), c.clone());
				g = gcd(x.clone(), (z.clone() - y.clone()).abs());

				if g > One::one() { break; }
			}
		}

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
