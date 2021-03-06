// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::ops::{Shr,Rem};
use std::mem::swap;

use std::cmp::min;

use num;
use num::{Zero, One, Integer};
use num::{FromPrimitive, ToPrimitive};
use rand::{Rng, rngs::SmallRng};
use num::bigint::{BigInt, RandBigInt};
use rand::distributions::uniform::SampleUniform;
use rand::FromEntropy;

use crate::factorize;
use crate::Factored;
use crate::prelude::*;
use crate::util::literal;
use crate::util::gcd;
use crate::util::mod_pow;
use crate::util::MoreNumCast;

// PollardBrentBigInt is a hack because BigInt does not
//  satisfy SampleUniform (it has its own trait which takes by ref >_>)

pub struct PollardBrent;
pub struct PollardBrentBigInt;

impl<T> TryFactor<T> for PollardBrent
where
    T: Clone + Zero + One + Integer + Shr<usize, Output = T> + SampleUniform + MoreNumCast,
{
    /// Produce a single factor of `x`.  PollardBrent is nondeterministic,
    /// and may sometimes fail to produce a non-trivial factor for composite `x`.
    fn try_factor(&self, x: &T) -> Option<T> {
        let mut rng = SmallRng::from_entropy();
        do_pollard(x, |a, b| rng.gen_range(a.clone(), b.clone()))
    }
}

impl TryFactor<BigInt> for PollardBrentBigInt {
    fn try_factor(&self, x: &BigInt) -> Option<BigInt> {
        let mut rng = SmallRng::from_entropy();
        do_pollard(x, |a, b| rng.gen_bigint_range(a, b))
    }
}

/// Produce a single factor of `x`.  PollardBrent is nondeterministic,
/// and may sometimes fail to produce a non-trivial factor for composite `x`.
fn do_pollard<'c, T, F>(x: &'c T, mut rand_range: F) -> Option<T>
where
    T: Clone + Zero + One + Integer + Shr<usize, Output = T> + MoreNumCast,
    F: for<'a, 'b> FnMut(&'a T, &'b T) -> T,
{
    // Adapted from https://comeoncodeon.wordpress.com/2010/09/18/pollard-rho-brent-integer-factorization/

    // In well-written functions, unused mutability is often indicative of a logic error.
    // So let it be known that this is NOT a well-written function:
    #![allow(unused_mut)]

    if x < &T::zero() {
        panic!("x < 0");
    }
    if x.is_even() {
        return Some(literal(2));
    }
    if x.is_multiple_of(&literal(3)) {
        return Some(literal(3));
    }
    if x == &T::one() {
        return None;
    }

    let one = One::one();
    let mut y: T = rand_range(&one, &x); // current value in the sequence:  y := y^2 + c (mod n)
    let mut c: T = rand_range(&one, &x); // parameter of y sequence
    let mut m: T = rand_range(&one, &x); // step size when multiplying crap together

    let mut g: T = One::one(); // contains the result
    let mut r: T = One::one(); // some kind of very coarse index
    let mut q: T = One::one(); // running product of `(z-y)` values

    let mut z: T = Zero::zero(); // Initial value of `y` for the current `r` value.
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

            if g > One::one() {
                break;
            }
        }
    }

    // At this point, we may or may not have identified a nontrivial factor.
    debug_assert!(g != One::one());
    if &g == x {
        None
    } else {
        Some(g)
    }
}

// computes (y**2 + c) % x
fn next_in_sequence<T>(y: T, x: T, c: T) -> T
where
    T: Clone + Integer,
{
    let mut result = y.clone() * y;
    result = result % x.clone();
    result = result + c;
    result = result % x;
    return result;
}
