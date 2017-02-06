// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

use std::ops::Shr;
use std::fmt::Debug;

use num::bigint::{ToBigUint,BigUint};
use num::{Zero,One,Integer};
use num::{ToPrimitive,FromPrimitive};

/// Services more types than NumCast does, such as BigInt.
pub trait MoreNumCast: ToPrimitive + FromPrimitive { }
impl<T: ToPrimitive + FromPrimitive> MoreNumCast for T { }

use test::Bencher;

/// Computes the greatest common divisor of two numbers using Euclid's method.
/// Behavior unspecified for negative numbers.
pub fn gcd<T>(a: T,  b: T) -> T
 where T: Clone + Zero + Integer,
{
	let mut cur_a = a;
	let mut cur_b = b;
	while (!cur_a.is_zero()) {
		let old_b = cur_b;
		cur_b = cur_a.clone();
		cur_a = old_b % cur_a;
	}
	cur_b
}

#[bench]
fn bench_gcd(b: &mut Bencher) {
	use rand::weak_rng;
	use rand::Rng;
	let mut rng = weak_rng();
	b.iter(|| {
		let a = rng.gen_range(100000u32,1000000u32);
		let b = rng.gen_range(100000u32,1000000u32);
		gcd(a,b)
	})
}

/// Performs an integer square root, returning the largest integer whose square is not
///  greater than the argument.
pub fn isqrt<T>(n: T) -> T
 where T: Clone + Zero + Integer + Shr<usize, Output=T> + MoreNumCast
{
	isqrt_fast(n.clone()).or_else(
		|| Some(isqrt_safe(n.clone()))
	).unwrap()
}

/// Used to convert an integral literal into an arbitrary type.
/// For zero and one, `num::Zero::zero()` and `num::One::one()` is preferred when they
///  are used as the additive/multiplicative identity, and `literal` is used otherwise.
#[inline]
pub fn literal<T: MoreNumCast>(n: i32) -> T
 where T: MoreNumCast,
{
	T::from_i32(n).unwrap()
}

/// Computes `pow(x, power) % modulus` using exponentation by squaring.
pub fn mod_pow<T,P>(x: T, power: P, modulus: T) -> T
 where T: Eq + Clone + Integer,
       P: Eq + Clone + Integer + Shr<usize, Output=P>,
{
	let mut prod:T = One::one();
	let mut remaining = power;
	let mut cur = x;
	while remaining > Zero::zero() {
		if remaining.is_odd() {
			prod = prod * cur.clone();
			prod = prod % modulus.clone();
		}
		remaining = remaining >> 1;
		cur = cur.clone() * cur;
		cur = cur % modulus.clone();
	}
	prod
}

#[test]
fn test_mod_pow()
{
	assert_eq!(mod_pow(234u64, 0, 1259), 1);
	assert_eq!(mod_pow(234u64, 1, 1259), 234);
	assert_eq!(mod_pow(234u64, 2412, 1259), 1091);
}

//-------------------------------
// isqrt helper methods

fn isqrt_fast<T>(x: T) -> Option<T>
 where T: MoreNumCast,
{
	x.to_f64().and_then(|f| {
		// Mantissa is 52 bits, and the square root takes half as many bits, so this
		//  may be a bit conservative.  The main concern is to avoid handling very
		//  large BigInts which may lose more than half of their precision.
		if f > 20f64.exp2() {
			None                 // Number too large, bail out!
		} else {
			T::from_f64(f.sqrt().floor())
		}
	})
}

fn isqrt_safe<T>(n: T) -> T
 where T: Clone + Zero + Integer + Shr<usize, Output=T>,
{
	// NOTE: while I'd like to remove the Shr bound, replacing '>> 1' with '/ 2' makes this
	//       algorithm take twice as long for BigInts :/
	if n.is_zero() { return Zero::zero(); }
	let mut x = n.clone();
	let mut y = (x.clone() + n.clone() / x.clone()) >> literal(1);
	while y < x {
		x = y.clone();
		y = (x.clone() + n.clone() / x.clone()) >> literal(1);
	}
	return x;
}

#[bench]
fn bench_fast(b: &mut Bencher) {
	b.iter(|| {
		(0usize..1000).map(|a| {isqrt_fast::<usize>(a).unwrap()}).collect::<Vec<usize>>()
	})
}

#[bench]
fn bench_safe(b: &mut Bencher) {
	b.iter(|| {
		(0usize..1000).map(isqrt_safe::<usize>).collect::<Vec<usize>>()
	})
}

#[bench]
fn bench_safe_bigint(b: &mut Bencher) {
	b.iter(|| {
		(0usize..1000).map(|a| isqrt_safe::<BigUint>(T::from(a).unwrap())).collect::<Vec<BigUint>>()
	})
}

#[bench]
fn bench_safe_massive_bigint(b: &mut Bencher) {
	use rand::XorShiftRng;
	use num::bigint::RandBigInt;
	let mut r = XorShiftRng::new_unseeded();
	b.iter(|| {
		(0usize..100).map(|_| {
			isqrt_safe::<BigUint>(r.gen_biguint(100usize))
		}).collect::<Vec<BigUint>>()
	})
}

#[test]
fn test_isqrt_consistency()
{
	(0usize..1000).map(|x| {
		let bigX = x.to_biguint().unwrap();
		assert_eq!(isqrt_fast(x),            Some(isqrt_safe(x)));
		assert_eq!(isqrt_fast(bigX.clone()), Some(isqrt_safe(bigX.clone())));
	}).collect::<Vec<_>>();
}

#[cfg(test)]
mod tests {
	use super::*;
	use num::NumCast;
	use num::bigint::{ToBigUint,BigUint};

	// need to test isqrt with BigInt more rigorously

}
