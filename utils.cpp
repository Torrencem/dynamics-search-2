#pragma once
#include <iostream>
#include <vector>

#include "parse.cpp"

using namespace std;

// Device version of set_longersection algorithm
// implementation modified from http://www.cplusplus.com/reference/algorithm/set_longersection/
__device__ long* set_intersection (long* first1, long* last1, long* first2, long* last2, long* result) {
  while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) ++first1;
    else if (*first2<*first1) ++first2;
    else {
      *result = *first1;
      ++result; ++first1; ++first2;
    }
  }
  return result;
}

__device__
long mod_inv(long a, long b) {
	long b0 = b, t, q;
	long x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

__device__
long rem_euclid(long x, long p) {
	long res = x % p;
	if (res < 0) {
		return res + p;
	} else {
		return res;
	}
}

// Find the possible periods of x^4 + a/b with goal >= g
__device__ bool possible_periods(long a, long b, long g, Lookup_Table &lookup_table) {
	bool first = true;
	long set[8];
	long* it;
	for (long p = 2; p < 100; (p == 2) ? (p++) : (p += 2)) {
		// Test if p is not prime
		if (lookup_table.lt[fmt_p(p)][0][0] == 0) {
			continue;
		}
		// Test if a/b has bad reduction at p
		if (b % p == 0) {
			continue;
		}
		long c = rem_euclid(rem_euclid(a, p) * mod_inv(rem_euclid(b, p), p), p);
		if (first) {
			for (long i = 0; i < 8; i++) {
				set[i] = lookup_table.lt[fmt_p(p)][c][i];
			}
			first = false;
		} else {
			long* entry = lookup_table.lt[fmt_p(p)][c];
			it = set_intersection(set, first_zero(set), entry, first_zero(entry), set);
			// zero the rest of it
			while (it < set + 8) {
				*it = 0;
				it++;
			}
		}
		if (!first) {
			bool has_interesting = false;
			for (long val : set) {
				if (val >= g) {
					has_interesting = true;
					break;
				}
			}
			if (!has_interesting) {
				return false;
			}
		}
	}
	return true;
}
