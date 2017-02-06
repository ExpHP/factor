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

use util::isqrt;
use util::literal;
use util::MoreNumCast;
use ::Factored;
use primes::PrimeTester;
use iter_ext::FactorExt;

mod sieve;
mod dixon;
mod pollard;
pub use self::sieve::FactorSieve;
pub use self::dixon::Dixon;
pub use self::pollard::PollardBrent;
pub use self::pollard::PollardBrentBigInt;

//// XXX: Would reduce some of the pain in keeping impl type bounds consistent, but would also require the
////      trait to be explicitly implemented for all applicable types, which arguably sucks just as much.
//// Master trait of traits needed to implement Factored<T>
//trait Factorable: Clone + num::Integer { }


// TODO: There should be a statically typed distinction for nondeterministic factorizers
//       which can fail to produce factors.  Alternatively, the library could simply not
//       publicly export any of these factorizers (exporting only reliable Stubborn
//       wrappers instead).

/// An interface for factorizing positive integers.
pub trait TryFactor<T>
 where T: Clone + Zero + One + Integer
{
	// TODO: Perhaps more sensible to return `1` for non-composite (so that division by the result
	//       is safe without needing to check for zero)
	/// Produces a single (not necessarily prime) factor of a number `x`.  Some implementations of
	///  `TryFactor` are deterministic and always produce the same factor for the same `x`, while
	///  others may have an element of randomness.
	///
	/// However, all Factorizers are expected to meet the following guarantees:
	///
	/// * `try_factor(x)` for any non-composite `x` (0, 1, or primes) returns `x`.
	/// * `try_factor(x)` for any composite `x` produces an arbitrary but non-trivial factor of `x`.
	///   (this factor is allowed to be composite)
	///
	/// Keep in mind that, in addition to the value returned, another factor can be obtained by
	///  dividing `x` by the value.
	fn try_factor(&self, x: &T) -> T;

	/// Builds a complete prime factorization of a number.  A default implementation is provided
	///  which calls `try_factor()` recursively on the factors produced.
	///
	/// Construction of a `Factorization` representing zero is not allowed, therefore
	/// `factorize(0)` is illegal.  Implementations are expected to check for this
	/// special case and emit an error.
	fn factorize(&self, x: T) -> Factored<T> where Self: Sized // FIXME why does this require Sized?
	{ self::helper::recursive_factorize(self, x) }
}

pub mod helper {
	use super::*;

	use num::{Zero, One, Integer};

	use ::Factored;

	/// The default implementation of `TryFactor::factorize`.
	///
	/// Factorizes a number by recursively factorizing both the
	///  produced factor and the remaining part.
	///
	/// Suitable for `Factorizers` which may return a composite number.
	pub fn recursive_factorize<F,T>(factorizer: &F, x: T) -> Factored<T>
	 where F: TryFactor<T>, T: Clone + Zero + One + Integer,
	{
		debug_assert!(x >= T::zero());
		if x == T::zero() { panic!("Zero has no factorization") }
		if x == T::one() { return One::one(); }

		let a = factorizer.try_factor(&x);

		// Non-composite point to themselves
		if a == x {
			return Factored::from_prime(a);

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
	/// A specialized implementation of `TryFactor::factorize` for certain types.
	///
	/// This is an optimized implementation for Factorizers which always produce
	///  the smallest nontrivial factor of any composite.
	pub fn always_smallest_factorize<F,T>(factorizer: &F, x: T) -> Factored<T>
	 where F: TryFactor<T>, T: Clone + Zero + One + Integer,
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
		let mut prev = factorizer.try_factor(&x);
		let mut x = x/prev.clone();
		let mut count = 1;
		debug_assert!(prev != T::one());

		let mut out = vec![];
		while x != T::one() {
			let p = factorizer.try_factor(&x);
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
		Factored::from_iter(out)
	}
}

/// Factors numbers using trial division.
///
/// Trial division is one of the most inefficient, yet most easily understood methods for
///  factorizing numbers.  `TrialDivision` tries dividing a number by successively
///  larger numbers until it encounters one that leaves a remainder of zero.
///
/// Despite its primitive nature, it can well outperform many of the more sophisticated methods
///  when factoring small numbers (TODO: of what magnitude?). However, it has trouble on numbers
///  with large prime factors.
pub struct TrialDivision;


impl<T> TryFactor<T>
for TrialDivision
 where T: Clone + Zero + One + Integer + Shr<usize, Output=T> + MoreNumCast,
{
	/// Produce a single factor of `x`.  TrialDivision is deterministic,
	///  and will always produce the smallest non-trivial factor of any composite number.
	///  Thus, the number it returns is also always prime.
	///
	/// The runtime scales linearly with the size of the smallest factor of `x`.
	fn try_factor(&self, x: &T) -> T
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

	fn factorize(&self, x: T) -> Factored<T>
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
/// Factors numbers by using results cached from another `TryFactor`.
///
/// Stores factors produced by another `TryFactor` for quick lookup. Only a single non-trivial
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
 where T: Clone + Integer + MoreNumCast,
{
	/// Constructs a `ListFactorizer` containing factors for the numbers `0..n`, using the
	///  provided `factorizer` to generate them.  Be sure to wrap any nondeterministic
	///  factorizer in a `Stubborn` beforehand to ensure that only correct results
	///  get cached in the list.
	pub fn compute_new(n: T, factorizer: &TryFactor<T>) -> Self {
		ListFactorizer {
			factors: num::iter::range(Zero::zero(), n).map(|x| factorizer.try_factor(&x)).collect(),
		}
	}
}

impl<T> TryFactor<T>
for ListFactorizer<T>
 where T: Clone + Integer + MoreNumCast,
{
	/// Produces the factor stored for `x`.
	///
	/// # Panics
	///
	/// May panic. (TODO: When?)
	#[inline]
	fn try_factor(&self, x: &T) -> T
	{
		self.factors[x.to_usize().unwrap()].clone()
	}
}

/// A `TryFactor` which doesn't take "no" for an answer.
///
/// It first tests the number for primality with its `PrimeTester` object.
/// If the number is not prime, the `Stubborn` will delegate to another `TryFactor`,
///  calling it repeatedly until a nontrivial factor is produced.
///
/// This makes it possible to factorize using nondeterministic `Factorizers` which can sometimes
///  fail to produce a nontrivial factor for composite numbers.
pub struct Stubborn<P,F,T>
 where P: PrimeTester<T>,
       F: TryFactor<T>,
       T: Clone + Zero + One + Integer,
{
	prime_tester: P,
	factorizer:   F,
	phantom:      PhantomData<T>,
}

impl<P,F,T> Stubborn<P,F,T>
 where P: PrimeTester<T>,
       F: TryFactor<T>,
       T: Clone + Zero + One + Integer,
{
	#[inline]
	pub fn new(prime_tester: P, factorizer: F) -> Self
	{
		Stubborn {
			prime_tester: prime_tester,
			factorizer:   factorizer,
			phantom:      PhantomData,
		}
	}
}

impl<P,F,T> TryFactor<T>
for Stubborn<P,F,T>
 where P: PrimeTester<T>,
       F: TryFactor<T>,
       T: Clone + Zero + One + Integer,
       T: Debug,
{
	fn try_factor(&self, x: &T) -> T
	{
		if self.prime_tester.is_composite(x) {
			// We are certain that x is composite, so keep trying to factor until we succeed
			loop {
				let factor = self.factorizer.try_factor(x);

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
	use num::{Zero, One, Integer, ToPrimitive};
	use test::Bencher;

	use util::literal;
	use util::MoreNumCast;
	use primes::PrimeSieve;
	use primes::MillerRabin;

	//  A simple test to factorize 242 as an arbitrary data type using an arbitrary factorizer.
	fn test_242<T, U>(factorizer: U)
	 where T: Clone + Debug + Integer + MoreNumCast,
	       U: TryFactor<T>,
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
			let k_t = T::from_usize(k).unwrap();  // cast k to T
			let v_t = T::from_usize(*v).unwrap(); // cast v to whatever the heck it is today
			assert_eq!(T::from_usize(factors.get(&k_t)).unwrap(), v_t);
		}
	}

	// Mix and match various factorizers and data types --- mostly to confirm that they compile
	//  (e.g. trait bounds are met, and whatnot)
	#[test]
	fn test_mix_and_match() {
		// differently sized integers
		test_242::<u32,_>(TrialDivision);
		test_242::<u64,_>(TrialDivision);
		test_242::<usize,_>(TrialDivision);

		// a signed type
		test_242::<isize,_>(TrialDivision);

		// my pain and suffering
		test_242::<BigUint,_>(TrialDivision);
		test_242::<BigInt,_>(TrialDivision);

		// Non-deterministic factorizers: Test them using a Stubborn
		let primes = PrimeSieve::new(256);
		test_242::<u64,_>(Stubborn::new(primes.clone(), Dixon::new(vec![2,3,5])));
		test_242::<i64,_>(Stubborn::new(primes.clone(), Dixon::new(vec![2,3,5])));
		test_242::<u64,_>(Stubborn::new(primes.clone(), PollardBrent));
		test_242::<i64,_>(Stubborn::new(primes.clone(), PollardBrent));
		test_242::<BigInt,_>(Stubborn::new(primes.clone(), PollardBrentBigInt));
	}

	fn make_list<F,T>(factorizer: F, limit: T) -> ListFactorizer<T>
	 where F: TryFactor<T>,
	       T: Clone + Debug + Integer + MoreNumCast,
	{
		return ListFactorizer::compute_new(limit, &factorizer);
	}

	fn make_list_stubborn<F,T>(factorizer: F, limit: T) -> ListFactorizer<T>
	 where F: TryFactor<T>,
	       T: Clone + Debug + Integer + MoreNumCast,
	{
		let primes = PrimeSieve::new(limit.to_usize().unwrap());
		let stubborn = Stubborn::new(primes, factorizer);
		return ListFactorizer::compute_new(limit, &stubborn);
	}

	// Builds a ListFactorizer up to some limit and verifies it against a list built
	//  using trial division. The factorizer provided to the macro will be wrapped
	//  in a Stubborn, so nondeterministic factorizers are OK.
	macro_rules! test_list_stubborn (
		($factorizer: expr, $limit: expr) => {
			{
				let expected = make_list(TrialDivision, $limit);
				let actual   = make_list_stubborn($factorizer, $limit);

				// We can't call factorize(0), but the value of try_factor(0) is still specified
				assert_eq!(expected.try_factor(&literal(0)), literal(0));

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
	//  Reason:  Dixon is legit broken and I don't feel like fixing it.
	//           I'm going to make Dixon private instead.
	#[test]
	fn test_list_dixon() {
		test_list_stubborn!(Dixon::new(vec![2,3,5,7]), 100000u64);
	}
	*/

	#[test]
	fn test_list_pollard() {
		test_list_stubborn!(PollardBrent, 100000u64);
	}

	/*
	// TODO: Make Dixon stop panicking so we can bench it.
	//       (also make it not suck)
	#[bench]
	fn bench_list_dixon(b: &mut Bencher) {
		b.iter(||
			make_list_stubborn(Dixon::new(vec![2,3,5]), 1000u64)
		);
	}
	*/

	#[bench]
	fn bench_list_trialdiv(b: &mut Bencher) {
		b.iter(||
			make_list(TrialDivision, 10000u64)
		);
	}

	#[bench]
	fn bench_factors_construction(b: &mut Bencher) {
		let sieve = ::factorizer::FactorSieve::new(10000u64);
		b.iter(|| {
			(1..10000).map(|x| sieve.factorize(x).iter().count())
				.fold(0, |a,b| a+b)
		});
	}

	#[bench]
	fn bench_factors_product(b: &mut Bencher) {
		let sieve = ::factorizer::FactorSieve::new(10000u64);
		b.iter(|| {
			(3000..4000).map(|x| sieve.factorize(x))
				.fold(::Factored::<u64>::one(), |a,b| a*b)
		});
	}

	#[bench]
	fn bench_list_pollard(b: &mut Bencher) {
		b.iter(||
			make_list_stubborn(PollardBrent, 10000u64)
		);
	}

	// A "rough" number (no small factors) around 10^8
	const TEN_8_ROUGH: u64  = 99400891; // 9967 * 9973

	// A rough square around 10^8  (some algorithms have extra trouble with repeated factors)
	const TEN_8_SQUARE: u64 = 99341089; // 9967 * 9967

	#[bench]
	fn bench_ten_8_rough_trialdiv(b: &mut Bencher) {
		b.iter(|| {
			TrialDivision.try_factor(&TEN_8_ROUGH)
		});
	}

	#[bench]
	fn bench_ten_8_rough_pollard(b: &mut Bencher) {
		let factorizer = Stubborn::new(MillerRabin, PollardBrent);
		b.iter(|| {
			factorizer.try_factor(&TEN_8_ROUGH)
		});
	}

	/*
	#[bench]
	fn bench_ten_8_rough_dixon(b: &mut Bencher) {
		let factorizer = Stubborn::new(MillerRabin, Dixon::new(vec![2,3,5]));
		b.iter(|| {
			let f = factorizer.try_factor(&TEN_8_ROUGH);
			println!("rough {:?}", f);
			f
		});
	}

	#[bench]
	fn bench_ten_8_square_dixon(b: &mut Bencher) {
		let factorizer = Stubborn::new(MillerRabin, Dixon::new(vec![2,3,5]));
		b.iter(|| {
			let f = factorizer.try_factor(&TEN_8_SQUARE);
			println!("square {:?}", f);
			f
		});
	}
	*/
}
