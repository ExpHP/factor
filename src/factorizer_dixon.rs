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
use std::default::Default;
use std::num::{ToPrimitive,FromPrimitive}; // and regret it
use std::hash::Hash;
use std::ops::{Shr,Rem};
use std::rand::Rng;
use std::rand::weak_rng;
use std::rand::distributions::range::SampleRange;
use std::mem::swap;

use num::{Zero, One, Integer};

use factorize;
use factorization::Factorization;
use factorizer::Factorizer;
use util::{isqrt,gcd};


/// TODO
pub struct DixonFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
{
	primes:       Vec<T>,
	extra_count:  usize,
	max_attempts: usize,
}

impl<T> DixonFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
{
	pub fn new(primes: Vec<T>) -> Self
	{
		DixonFactorizer {
			primes:       primes,
			extra_count:  3,
			max_attempts: 100,
		}
	}
}

impl<T> Factorizer<T>
for DixonFactorizer<T>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Shr<usize, Output=T> + Hash<Hasher> + SampleRange,
{
	/// Produce a single factor of `x`.  TrialDivisionFactorizer is deterministic,
	///  and will always produce the smallest non-trivial factor of any composite number.
	///  Thus, the number it returns is also always prime.
	///
	/// The runtime scales linearly with the size of the smallest factor of `x`.
	fn get_factor(self: &Self, x: &T) -> T
	{
		// Step 1: Collect congruences of the form a^2 = b (mod x), where b < x
		//          and b is smooth (composed only of small primes).
		let a_min = isqrt(x.clone()); // XXX: ceil? (ensure a^2 > x)

		let a_count = self.primes.len() + self.extra_count;
		let mut a_values: Vec<T> = Vec::new();
		let mut b_factorizations: Vec<Factorization<T>> = Vec::new();

		for _ in (0us..a_count) {
			for _ in (0us..self.max_attempts) {

				// a random number such that a < x and a*a > x.
				let a = weak_rng().gen_range(a_min.clone(), x.clone());
				let b = a.clone() * a.clone() % x.clone();

				// there are cases where x may divide directly into a*a, in which
				//  case b = 0.  In this case, gcd(a,x) must be a nontrivial
				//  factor (easy proof), and we can bail with an early result.
				if b.is_zero() {
					let candidate = gcd(a.clone(), x.clone());

					// Just to be sure...
					assert!(candidate != One::one());
					assert!(candidate != x.clone());
					return candidate;
				}

				// Try to factorize b using only small primes
				let b_factorization = factorize_limited(b, &self.primes);

				if b_factorization.is_some() {
					a_values.push(a);
					b_factorizations.push(b_factorization.unwrap());
				}
			}
		}

		// Step 2: Find products of b's which are square.

		// Use a bit array to represent each b's factorization mod 2
		let mut bitmatrix = bit_matrix_from_factorizations(&b_factorizations, &self.primes);

		// Put in row echelon form
		bit_matrix_to_ref(&mut bitmatrix);

		// Each row full of zeros in the matrix represents a set of b values that multiply
		//  together to form a square.
		for matrix_row in bitmatrix.into_rows().into_iter() {

			if matrix_row.is_all_zero() {

				let mut a_prod: T = One::one();
				let mut b_prod_factors: Factorization<T> = Default::default();

				for index in matrix_row.into_index_set().iter() {
					a_prod = a_prod * a_values[index].clone();
					b_prod_factors = b_prod_factors * b_factorizations[index].clone();
				}

				// we now have a congruence of squares (mod x) between a_prod^2 and b_prod
				let b_prodsqrt_factors = b_prod_factors.sqrt().unwrap();
				let b_prodsqrt = b_prodsqrt_factors.product();

				// a - sqrt(b) has a high chance of sharing a nontrivial factor in common with x
				let candidate = gcd(a_prod - b_prodsqrt,  x.clone());

				if candidate != One::one() && candidate != x.clone() {
					return candidate;
				}
			}
		}
		// temporary solution: panic on failure to find factor (FIXME)
		panic!("No factor found! :(");
	}
}



// Utility function that only returns a factorization if it can be constructed
//  *exclusively* from the given primes.
fn factorize_limited<T>(x: T, primes: &Vec<T>) -> Option<Factorization<T>>
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher> + SampleRange
{
	assert!(!x.is_zero());

	let mut f: Factorization<T> = Default::default();
	let mut c = x;
	for p in primes.iter() {

		let mut count = 0us;
		while c.is_multiple_of(p) {
			c = c / p.clone();
			count += 1;
		}

		f.set(p.clone(), count);
	}

	// Only report complete factorizations
	if c == One::one() {
		Some(f)
	} else {
		None
	}
}

#[test]
fn factorize_limited_test() {
	let primes = vec![2us, 5, 7];
	println!("{:?},{:?}",factorize_limited(1, &primes), Some(factorize(1us)));
	assert_eq!(factorize_limited(1, &primes), Some(factorize(1us)));
	assert_eq!(factorize_limited(2450, &primes), Some(factorize(2450us)));
	assert_eq!(factorize_limited(22, &primes), None);   // 11 not in prime list
	assert_eq!(factorize_limited(12, &primes), None);   // 3 not in prime list
}

//-------------------------------------------
// Private utility structs used in the "bit matrix to ref" algorithm, to separate
//  the underlying data representation from the general algorithm.
#[derive(Clone,Debug)]
struct DixonBitvec {
	elements: BitvSet, // represents the power of each prime modulo 2
	indices:  BitvSet, // indicates which rows have been xor'ed to make this row
}

#[derive(Clone,Debug)]
struct DixonBitmatrix {
	rows: Vec<DixonBitvec>,
}

impl DixonBitvec
{
	#[inline]
	fn is_all_zero(self: &Self) -> bool {
		self.elements.is_empty()
	}

	#[inline]
	fn into_index_set(self: Self) -> BitvSet {
		self.indices
	}

}

impl DixonBitmatrix
{
	// Matrix dimensions
	#[inline]
	fn nrows(self: &Self) -> usize { self.rows.len() }
	#[inline]
	fn ncols(self: &Self) -> usize { self.rows[0].elements.len() }

	// Index
	#[inline]
	fn get_elem(self: &Self, row: usize, col: usize) -> bool {
		self.rows[row].elements.contains(&col)
	}

	#[inline]
	fn swap_rows(self: &mut Self, i: usize, j: usize) {
		let temp = self.rows[i].clone();
		self.rows[i] = self.rows[j].clone();
		self.rows[j] = temp;
	}

	// Computes an XOR of rows src and dest, overwriting dest.
	#[inline]
	fn xor_update_row(self: &mut Self, src: usize, dest: usize) {
		let elems = self.rows[src].elements.clone();
		let inds  = self.rows[src].indices.clone();
		self.rows[dest].elements.symmetric_difference_with(&elems);
		self.rows[dest].indices.symmetric_difference_with(&inds);
	}

	#[inline]
	fn into_rows(self: Self) -> Vec<DixonBitvec> {
		self.rows
	}
}

// Produce matrix from initial input
fn bit_matrix_from_factorizations<T>(factorizations: &Vec<Factorization<T>>, primes: &Vec<T>) -> DixonBitmatrix
 where T: Eq + Clone + FromPrimitive + ToPrimitive + Zero + One + Integer + Hash<Hasher>
{
	let rows: Vec<DixonBitvec> = factorizations.iter().enumerate().map(|(row_index,fact)| {

		// set elements equal to powers in factorization, mod 2
		let elements: Bitv = primes.iter().map(|p| fact.get(p) % 2 == 1).collect();

		// indices initially contains just the index for this row
		let mut indices = BitvSet::new();
		indices.insert(row_index);

		// Construct the row
		DixonBitvec {
			elements: BitvSet::from_bitv(elements),
			indices:  indices,
		}
	}).collect();

	DixonBitmatrix { rows: rows }
}

// Manipulate bit matrix into row echelon form
fn bit_matrix_to_ref(matrix: &mut DixonBitmatrix) {

	let mut target_row = 0us;

	for col in (0us..matrix.ncols()) {

		// Look for a leading 1 in this column
		for source_row in (target_row..matrix.nrows()) {

			if matrix.get_elem(source_row, col) {

				// Move this row to its correct location (target_row)
				matrix.swap_rows(source_row, target_row);

				// Eliminate any remaining 1s below this one in the column by XORing
				for other_row in ((source_row+1)..matrix.nrows()) {
					if matrix.get_elem(other_row, col) {
						matrix.xor_update_row(target_row, other_row);
					}
				}

				target_row += 1; // done with this target_row
				break;           // also done with this column
			}
		}
		// (invariant: elements with row >= target_row,  col < leading_col are 0)
	}
}

