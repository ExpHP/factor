
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
