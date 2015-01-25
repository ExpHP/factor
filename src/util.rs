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
use std::default::Default;
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;
use std::ops::Shr;

use num::bigint::{ToBigUint,BigUint};
use num::{Zero,One,Integer};
use std::num::Float;

use test::Bencher;

/// Computes the greatest common divisor of two numbers using Euclid's method.
/// Behavior unspecified for negative numbers.
pub fn gcd<T>(a: T,  b: T) -> T
 where T: Clone + Zero + Integer,
{
	if (a == Zero::zero()) {
		b
	} else {
		gcd(b % a.clone(), a)
	}
}

/// Performs an integer square root, returning the largest integer whose square is not
///  greater than the argument.
pub fn isqrt<T>(n: T) -> T
 where T: Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T>,
{
	isqrt_fast(n.clone()).or_else(
		|| Some(isqrt_safe(n.clone()))
	).unwrap()
}

//-------------------------------
// isqrt helper methods

fn isqrt_fast<T>(x: T) -> Option<T>
 where T: FromPrimitive + ToPrimitive,
{
	x.to_f64().and_then(|f| {
		// Mantissa is 52 bits, and the square root takes half as many bits, so this
		//  may be a bit conservative.  The main concern is to avoid handling very
		//  large BigInts which may lose more than half of their precision.
		if f > 20f64.exp2() {
			None                 // Number too large, bail out!
		} else {
			FromPrimitive::from_f64(f.sqrt().floor())
		}
	})
}

fn isqrt_safe<T>(n: T) -> T
 where T: Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T>,
{
	// NOTE: while I'd like to remove the Shr bound, replacing '>> 1' with '/ 2' makes this
	//       algorithm take twice as long for BigInts :/
	if n.is_zero() { return Zero::zero(); }
	let mut x = n.clone();
	let mut y = (x.clone() + n.clone() / x.clone()) >> One::one();
	while y < x {
		x = y.clone();
		y = (x.clone() + n.clone() / x.clone()) >> One::one();
	}
	return x;
}

#[bench]
fn bench_fast(b: &mut Bencher) {
	b.iter(|| {
		range(0us, 1000us).map(|a| {isqrt_fast::<usize>(a).unwrap()}).collect::<Vec<usize>>()
	})
}

#[bench]
fn bench_safe(b: &mut Bencher) {
	b.iter(|| {
		range(0us, 1000us).map(isqrt_safe::<usize>).collect::<Vec<usize>>()
	})
}

#[bench]
fn bench_safe_bigint(b: &mut Bencher) {
	b.iter(|| {
		range(0us, 1000us).map(|a| isqrt_safe::<BigUint>(FromPrimitive::from_uint(a).unwrap())).collect::<Vec<BigUint>>()
	})
}

#[bench]
fn bench_safe_massive_bigint(b: &mut Bencher) {
	use std::rand::XorShiftRng;
	use num::bigint::RandBigInt;
	let mut r = XorShiftRng::new_unseeded();
	b.iter(|| {
		range(0us, 100us).map(|_| {
			isqrt_safe::<BigUint>(r.gen_biguint(100us))
		}).collect::<Vec<BigUint>>()
	})
}

#[test]
fn test_isqrt_consistency()
{
	(0us..1000).map(|x| {
		let bigX = x.to_biguint().unwrap();
		assert_eq!(isqrt_fast(x),            Some(isqrt_safe(x)));
		assert_eq!(isqrt_fast(bigX.clone()), Some(isqrt_safe(bigX.clone())));
	});
}

#[cfg(test)]
mod tests {
	extern crate test;
	use super::*;
	use std::num::FromPrimitive;
	use std::num::ToPrimitive;
	use num::bigint::{ToBigUint,BigUint};

	// need to test isqrt with BigInt more rigorously

}
