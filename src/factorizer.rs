// Copyright 2015 Michael 'ExpHP' Lamparski
//
// Licensed under the terms of the MIT License, available at:
//  http://opensource.org/licenses/MIT
// and also included in the file COPYING at the root of this distribution.
// This file may not be copied, modified, or distributed except according
// to those terms.

// FIXME This file has a disabled test!

use std::ops::{Shr,Rem};
use std::fmt::Debug;
use std::marker::PhantomData;

use num;
use num::{Zero, One, Integer};
use num::{FromPrimitive, ToPrimitive};

use util::isqrt;
use util::literal;
use ::Factors;
use primes::PrimeTester;
use iter_ext::FactorExt;

// Type-Synonyms, for semantic purposes.
// (alas, re-exports don't appear to take docstrings)

// /// A deterministic `Factorizer` that is guaranteed to work on any number, but may be fast.
pub use TrialDivisionFactorizer as SafeFactorizer;
// /// The `Factorizer` used by the `factor::factorize` method.
pub use TrialDivisionFactorizer as DefaultFactorizer;

pub use factorizer_sieve::FactorSieve;
pub use factorizer_dixon::DixonFactorizer;
pub use factorizer_pollard::PollardBrentFactorizer;
pub use factorizer_pollard::PollardBrentFactorizerBigInt;

//// XXX: Would reduce some of the pain in keeping impl type bounds consistent, but would also require the
////      trait to be explicitly implemented for all applicable types, which arguably sucks just as much.
//// Master trait of traits needed to implement Factors<T>
//trait Factorable: Clone + num::Integer { }


// TODO: There should be a statically typed distinction for nondeterministic factorizers
//       which can fail to produce factors.  Alternatively, the library could simply not
//       publicly export any of these factorizers (exporting only reliable StubbornFactorizer
//       wrappers instead).

/// An interface for factorizing positive integers.
pub trait Factorizer<T>
 where T: Clone + Zero + One + Integer
{
	// TODO: Perhaps more sensible to return `1` for non-composite (so that division by the result
	//       is safe without needing to check for zero)
	/// Produces a single (not necessarily prime) factor of a number `x`.  Some implementations of
	///  `Factorizer` are deterministic and always produce the same factor for the same `x`, while
	///  others may have an element of randomness.
	///
	/// However, all Factorizers are expected to meet the following guarantees:
	///
	/// * `get_factor(x)` for any non-composite `x` (0, 1, or primes) returns `x`.
	/// * `get_factor(x)` for any composite `x` produces an arbitrary but non-trivial factor of `x`.
	///   (this factor is allowed to be composite)
	///
	/// Keep in mind that, in addition to the value returned, another factor can be obtained by
	///  dividing `x` by the value.
	fn get_factor(&self, x: &T) -> T;

	/// Builds a complete prime factorization of a number.  A default implementation is provided
	///  which calls `get_factor()` recursively on the factors produced.
	///
	/// Construction of a `Factorization` representing zero is not allowed, therefore
	/// `factorize(0)` is illegal.  Implementations are expected to check for this
	/// special case and emit an error.
	fn factorize(&self, x: T) -> Factors<T> where Self: Sized // FIXME why does this require Sized?
	{ self::helper::recursive_factorize(self, x) }
}

pub mod helper {
	use super::*;

	use num::{Zero, One, Integer};

	use ::Factors;

	/// The default implementation of `Factorizer::factorize`.
	///
	/// Factorizes a number by recursively factorizing both the
	///  produced factor and the remaining part.
	///
	/// Suitable for `Factorizers` which may return a composite number.
	pub fn recursive_factorize<F,T>(factorizer: &F, x: T) -> Factors<T>
	 where F: Factorizer<T>, T: Clone + Zero + One + Integer,
	{
		debug_assert!(x >= T::zero());
		if x == T::zero() { panic!("Zero has no factorization") }
		if x == T::one() { return One::one(); }

		let a = factorizer.get_factor(&x);

		// Non-composite point to themselves
		if a == x {
			return Factors::from_prime(a);

		// Composite numbers
		} else {
			let b = x.clone() / a.clone();
			let a_facs = recursive_factorize(factorizer, a);
			let b_facs = recursive_factorize(factorizer, b);
			return a_facs * b_facs;
		};
	}

	// NOTE: informal benchmarks indicate that the detectable speedup here is fairly
	// modest (a factor of ~2-3).
	/// A specialized implementation of `Factorizer::factorize` for certain types.
	///
	/// This is an optimized implementation for Factorizers which always produce
	///  the smallest nontrivial factor of any composite.
	pub fn always_smallest_factorize<F,T>(factorizer: &F, x: T) -> Factors<T>
	 where F: Factorizer<T>, T: Clone + Zero + One + Integer,
	{
		debug_assert!(x >= T::zero());
		if x == T::zero() { panic!("Zero has no factorization") }
		if x == T::one() { return One::one(); }

		// the rest of this is basically trying to implement something like the
		// following in iterator-speak (but there's no grouping iterator in std,
		// and it is still painful to write an iterator in rust)
		//    iter.group_by(|p| p)
		//        .map(|(p, group)| (p, group.count()))
		// where iter is an imaginary iterator giving one prime at a time (with repeats)

		// begin first group
		let mut prev = factorizer.get_factor(&x);
		let mut x = x/prev.clone();
		let mut count = 1;
		debug_assert!(prev != T::one());

		let mut out = vec![];
		while x != T::one() {
			let p = factorizer.get_factor(&x);
			if p != prev {
				debug_assert!(p > prev, "non-sorted primes in always_smallest_factorize");
				// start new group
				out.push((prev, count));
				prev = p.clone();
				count = 0;
			}
			count += 1;
			x = x / p;
		}
		out.push((prev, count)); // final group
		Factors::from_iter(out)
	}
}

/// Factors numbers using trial division.
///
/// Trial division is one of the most inefficient, yet most easily understood methods for
///  factorizing numbers.  `TrialDivisionFactorizer` tries dividing a number by successively
///  larger numbers until it encounters one that leaves a remainder of zero.
///
/// Despite its primitive nature, it can well outperform many of the more sophisticated methods
///  when factoring small numbers (TODO: of what magnitude?). However, it has trouble on numbers
///  with large prime factors.
pub struct TrialDivisionFactorizer;


impl<T> Factorizer<T>
for TrialDivisionFactorizer
 where T: Clone + Zero + One + Integer + Shr<usize, Output=T> + FromPrimitive + ToPrimitive,
{
	/// Produce a single factor of `x`.  TrialDivisionFactorizer is deterministic,
	///  and will always produce the smallest non-trivial factor of any composite number.
	///  Thus, the number it returns is also always prime.
	///
	/// The runtime scales linearly with the size of the smallest factor of `x`.
	fn get_factor(&self, x: &T) -> T
	{
		if x.is_zero() { return Zero::zero() };

		if x.is_even() { return literal(2); }

		let start: T = literal(3);
		let stop:  T = isqrt(x.clone()) + literal(2);
		let step:  T = literal(2);
		
		//for odd in range_step(start, stop, step)  // needs trait Int
		let mut odd = start;
		while odd < stop {
			if x.is_multiple_of(&odd) {
				return odd;
			}
			odd = odd + step.clone();
		};

		// x is prime
		return x.clone();
	}

	fn factorize(&self, x: T) -> Factors<T>
	{ self::helper::always_smallest_factorize(self, x) }
}

// TODO
//pub struct FermatFactorizer<T>
// where T: Clone + Zero + One + Integer
//;

// TODO
//pub struct GeneralFactorizer<T>
// where T: Clone + Zero + One + Integer
//;

// FIXME this shouldn't exist except maybe for testing purposes.
//       FactorSieve is the way to go.
/// Factors numbers by using results cached from another `Factorizer`.
///
/// Stores factors produced by another `Factorizer` for quick lookup. Only a single non-trivial
///  factor is stored for each composite number, from which the full decomposition can be
///  gathered recursively.
#[derive(Clone, Debug)]
pub struct ListFactorizer<T>
{
	factors: Vec<T>,
}

// TODO: Not sure how to do static dispatch here under the new orphan-checking rules.
//       Perhaps fix this up if/when unboxed abstract types make an appearance
impl<T> ListFactorizer<T>
 where T: Clone + Integer + FromPrimitive + ToPrimitive,
{
	/// Constructs a `ListFactorizer` containing factors for the numbers `0..n`, using the
	///  provided `factorizer` to generate them.  Be sure to wrap any nondeterministic
	///  factorizer in a `StubbornFactorizer` beforehand to ensure that only correct results
	///  get cached in the list.
	pub fn compute_new(n: T, factorizer: &Factorizer<T>) -> Self {
		ListFactorizer {
			factors: num::iter::range(literal(0), n).map(|x| factorizer.get_factor(&x)).collect(),
		}
	}
}

impl<T> Factorizer<T>
for ListFactorizer<T>
 where T: Clone + Integer + ToPrimitive,
{
	/// Produces the factor stored for `x`.
	///
	/// # Panics
	///
	/// May panic. (TODO: When?)
	#[inline]
	fn get_factor(&self, x: &T) -> T
	{
		self.factors[x.to_usize().unwrap()].clone()
	}
}

/// A `Factorizer` which doesn't take "no" for an answer.
///
/// It first tests the number for primality with its `PrimeTester` object.
/// If the number is not prime, the `StubbornFactorizer` will delegate to another `Factorizer`,
///  calling it repeatedly until a nontrivial factor is produced.
///
/// This makes it possible to factorize using nondeterministic `Factorizers` which can sometimes
///  fail to produce a nontrivial factor for composite numbers.
pub struct StubbornFactorizer<P,F,T>
 where P: PrimeTester<T>,
       F: Factorizer<T>,
       T: Clone + Zero + One + Integer,
{
	prime_tester: P,
	factorizer:   F,
	phantom:      PhantomData<T>,
}

impl<P,F,T> StubbornFactorizer<P,F,T>
 where P: PrimeTester<T>,
       F: Factorizer<T>,
       T: Clone + Zero + One + Integer,
{
	#[inline]
	pub fn new(prime_tester: P, factorizer: F) -> Self
	{
		StubbornFactorizer {
			prime_tester: prime_tester,
			factorizer:   factorizer,
			phantom:      PhantomData,
		}
	}
}

impl<P,F,T> Factorizer<T>
for StubbornFactorizer<P,F,T>
 where P: PrimeTester<T>,
       F: Factorizer<T>,
       T: Clone + Zero + One + Integer,
       T: Debug,
{
	fn get_factor(&self, x: &T) -> T
	{
		if self.prime_tester.is_composite(x) {
			// We are certain that x is composite, so keep trying to factor until we succeed
			loop {
				let factor = self.factorizer.get_factor(x);

				if &factor != x { return factor; };
			}

		} else {
			// x is 0, 1, or prime.  Or so says the prime tester, at least.
			x.clone()
		}
	}
}

// Tests

#[cfg(test)]
mod tests {
	use super::*;

	use std::fmt::Debug;

	use num;
	use num::{BigUint, BigInt};
	use num::{Zero, One, Integer};
	use num::{FromPrimitive, ToPrimitive};
	use test::Bencher;

	use util::literal;
	use primes::PrimeSieve;
	use primes::MillerRabinTester;

	//  A simple test to factorize 242 as an arbitrary data type using an arbitrary factorizer.
	fn test_242<T, U>(factorizer: U)
	 where T: Clone + Debug + Integer + FromPrimitive,
	       U: Factorizer<T>,
	{
		// 242 = 2 * 11 * 11
		let x_t: T = literal(242);

		let factors = factorizer.factorize(x_t.clone());

		// vector of expected powers
		let mut expected = vec![0usize; 15];
		expected[2]  = 1;
		expected[11] = 2;

		// Check factorization against expected		
		for (k,v) in expected.iter().enumerate() {
			let k_t = FromPrimitive::from_usize(k).unwrap();  // cast k to T
			let v_t = FromPrimitive::from_usize(*v).unwrap(); // cast v to whatever the heck it is today
			assert_eq!(factors.get(&k_t), v_t);
		}
	}

	// Mix and match various factorizers and data types --- mostly to confirm that they compile
	//  (e.g. trait bounds are met, and whatnot)
	#[test]
	fn test_mix_and_match() {
		// differently sized integers
		test_242::<u32,_>(TrialDivisionFactorizer);
		test_242::<u64,_>(TrialDivisionFactorizer);
		test_242::<usize,_>(TrialDivisionFactorizer);

		// a signed type
		test_242::<isize,_>(TrialDivisionFactorizer);
	
		// my pain and suffering
		test_242::<BigUint,_>(TrialDivisionFactorizer);
		test_242::<BigInt,_>(TrialDivisionFactorizer);

		// Non-deterministic factorizers: Test them using a StubbornFactorizer
		let primes = PrimeSieve::new(256);
		test_242::<u64,_>(StubbornFactorizer::new(primes.clone(), DixonFactorizer::new(vec![2,3,5])));
		test_242::<i64,_>(StubbornFactorizer::new(primes.clone(), DixonFactorizer::new(vec![2,3,5])));
		test_242::<u64,_>(StubbornFactorizer::new(primes.clone(), PollardBrentFactorizer));
		test_242::<i64,_>(StubbornFactorizer::new(primes.clone(), PollardBrentFactorizer));
		test_242::<BigInt,_>(StubbornFactorizer::new(primes.clone(), PollardBrentFactorizerBigInt));
	}

	fn make_list<F,T>(factorizer: F, limit: T) -> ListFactorizer<T>
	 where F: Factorizer<T>,
	       T: Clone + Debug + Integer + FromPrimitive + ToPrimitive,
	{
		return ListFactorizer::compute_new(limit, &factorizer);
	}

	fn make_list_stubborn<F,T>(factorizer: F, limit: T) -> ListFactorizer<T>
	 where F: Factorizer<T>,
	       T: Clone + Debug + Integer + FromPrimitive + ToPrimitive,
	{
		let primes = PrimeSieve::new(limit.to_usize().unwrap());
		let stubborn = StubbornFactorizer::new(primes, factorizer);
		return ListFactorizer::compute_new(limit, &stubborn);
	}

	// Builds a ListFactorizer up to some limit and verifies it against a list built
	//  using trial division. The factorizer provided to the macro will be wrapped
	//  in a StubbornFactorizer, so nondeterministic factorizers are OK.
	macro_rules! test_list_stubborn (
		($factorizer: expr, $limit: expr) => {
			{
				let expected = make_list(TrialDivisionFactorizer, $limit);
				let actual   = make_list_stubborn($factorizer, $limit);

				// We can't call factorize(0), but the value of get_factor(0) is still specified
				assert_eq!(expected.get_factor(&literal(0)), literal(0));

				// The elements of the two lists may differ, but the complete factorization
				//  of each number must agree:
				for i in num::iter::range(literal(1), $limit) {
					assert_eq!(expected.factorize(i.clone()), actual.factorize(i));
				};
			}
		}
	);

	/*
	// FIXME Disabled test!!!
	//  Reason:  DixonFactorizer is legit broken and I don't feel like fixing it.
	//           I'm going to make DixonFactorizer private instead.
	#[test]
	fn test_list_dixon() {
		test_list_stubborn!(DixonFactorizer::new(vec![2,3,5,7]), 100000u64);
	}
	*/

	#[test]
	fn test_list_pollard() {
		test_list_stubborn!(PollardBrentFactorizer, 100000u64);
	}

	/*
	// TODO: Make Dixon stop panicking so we can bench it.
	//       (also make it not suck)
	#[bench]
	fn bench_list_dixon(b: &mut Bencher) {
		b.iter(||
			make_list_stubborn(DixonFactorizer::new(vec![2,3,5]), 1000u64)
		);
	}
	*/

	#[bench]
	fn bench_list_trialdiv(b: &mut Bencher) {
		b.iter(||
			make_list(TrialDivisionFactorizer, 10000u64)
		);
	}

	#[bench]
	fn bench_factors_construction(b: &mut Bencher) {
		let sieve = ::FactorSieve::new(10000u64);
		b.iter(|| {
			(1..10000).map(|x| sieve.factorize(x).iter().count())
				.fold(0, |a,b| a+b)
		});
	}

	#[bench]
	fn bench_factors_product(b: &mut Bencher) {
		let sieve = ::FactorSieve::new(10000u64);
		b.iter(|| {
			(3000..4000).map(|x| sieve.factorize(x))
				.fold(::Factors::<u64>::one(), |a,b| a*b)
		});
	}

	#[bench]
	fn bench_list_pollard(b: &mut Bencher) {
		b.iter(||
			make_list_stubborn(PollardBrentFactorizer, 10000u64)
		);
	}

	// A "rough" number (no small factors) around 10^8
	const TEN_8_ROUGH: u64  = 99400891; // 9967 * 9973

	// A rough square around 10^8  (some algorithms have extra trouble with repeated factors)
	const TEN_8_SQUARE: u64 = 99341089; // 9967 * 9967

	#[bench]
	fn bench_ten_8_rough_trialdiv(b: &mut Bencher) {
		b.iter(|| {
			TrialDivisionFactorizer.get_factor(&TEN_8_ROUGH)
		});
	}

	#[bench]
	fn bench_ten_8_rough_pollard(b: &mut Bencher) {
		let factorizer = StubbornFactorizer::new(MillerRabinTester, PollardBrentFactorizer);
		b.iter(|| {
			factorizer.get_factor(&TEN_8_ROUGH)
		});
	}

	/*
	#[bench]
	fn bench_ten_8_rough_dixon(b: &mut Bencher) {
		let factorizer = StubbornFactorizer::new(MillerRabinTester, DixonFactorizer::new(vec![2,3,5]));
		b.iter(|| {
			let f = factorizer.get_factor(&TEN_8_ROUGH);
			println!("rough {:?}", f);
			f
		});
	}

	#[bench]
	fn bench_ten_8_square_dixon(b: &mut Bencher) {
		let factorizer = StubbornFactorizer::new(MillerRabinTester, DixonFactorizer::new(vec![2,3,5]));
		b.iter(|| {
			let f = factorizer.get_factor(&TEN_8_SQUARE);
			println!("square {:?}", f);
			f
		});
	}
	*/
}
