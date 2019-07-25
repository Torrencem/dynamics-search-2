#include<iostream>
#include <cstdlib>
#include <math.h>

#include "utils.cpp"
#include "cipolla.cpp"

using namespace std;

class gauss_int {
	private:
	public:
		long a, b;
		
		__device__ gauss_int(long a, long b) {
			this->a = a;
			this->b = b;
		}

		__device__ gauss_int() {
			this->a = 0;
			this->b = 0;
		}

		__device__ bool is_unit() const {
			if (abs(this->a) == 1 && this->b == 0) {
				return true;
			} else if (abs(this->b) == 1 && this->a == 0) {
				return true;
			}
			return false;
		}

		__device__ bool is_zero() const {
			return this->a == 0 && this->b == 0;
		}

		__device__ gauss_int operator+(gauss_int const &other) const {
			gauss_int ret(this->a + other.a, this->b + other.b);
			return ret;
		}

		__device__ gauss_int operator-(gauss_int const &other) const {
			gauss_int ret(this->a - other.a, this->b - other.b);
			return ret;
		}

		__device__ gauss_int operator*(gauss_int const &other) const {
			gauss_int ret(this->a * other.a - this->b * other.b, this->a * other.b + this->b * other.a);
			return ret;
		}

		__device__ gauss_int conjugate() const {
			gauss_int ret(this->a, -this->b);
			return ret;
		}

		__device__ long norm_sq() const {
			return this->a * this->a + this->b * this->b;
		}

		__device__ gauss_int operator/(gauss_int const &other) const {
			gauss_int numer = (*this) * other.conjugate();
			long normsq = other.norm_sq();
			double a = (double)(numer.a) / (double)(normsq);
			double b = (double)(numer.b) / (double)(normsq);
			gauss_int ret((long)round(a), (long)round(b));
			return ret;
		}

		__device__ gauss_int gcd(eis_int const &other) const {
			if (this->norm_sq() < other.norm_sq()) {
				return other.gcd(*this);
			}
			gauss_int a = *this;
			gauss_int b = other;
			gauss_int remainder;
			do {
				gauss_int q = a / b;
				remainder = a - (b * q);
				a = b;
				b = remainder;
			} while (!b.is_zero());
			return a;
		}

		__device__ long* reductions(long p) const {
			long d2   = mod_inv(2, p);
			cipolla_res s   = cipolla(p - 1, p);
			long w[2] = {(-1 + s.fst)*d2, (-1 + s.snd)*d2};
			long* ret;
			ret = (long*) malloc(2*sizeof(long));
			ret[0] = rem_euclid(this->a + w[0] * this->b, p);
			ret[1] = rem_euclid(this->a + w[1] * this->b, p);
			return ret;
		}

		__device__ float phase_angle(bool debug=false) {
			float cx = __ll2float_rd(this->a);
			float cy = __ll2float_rd(this->b); 
			float val = atan2f(cy, cx);
			if (val <= 0.0f) {
				val += 2*3.141593f;
			}
			return val;
		}
};

ostream &operator<<(ostream &os, eis_int const &m) {
    return os << m.a << " + w*" << m.b;
}

class qi_elem {
	private:
	public:
		gauss_int numer, denom;

		__device__ qw_elem(gauss_int numer, gauss_int denom) {
			this->numer = numer;
			this->denom = denom;
		}

		__device__ qw_elem() {
			this->numer = gauss_int();
			this->denom = gauss_int();
		}

		__device__ long* reductions(long p) {
			long* a = this->numer.reductions(p);
			long* b = this->denom.reductions(p);

			long* ret;
			ret = (long*) malloc(sizeof(long) * 2);

			if (b[0] == 0) {
				ret[0] = -1;
			} else {
				ret[0] = rem_euclid(a[0] * mod_inv(b[0], p), p);
			}

			if (b[1] == 0) {
				ret[1] = -1;
			} else {
				ret[1] = rem_euclid(a[1] * mod_inv(b[1], p), p);
			}
			free(a);
			free(b);
			return ret;
		}
};

__device__ bool possible_periods_qi(qi_elem &x, long g, Lookup_Table &lookup_table) {
        bool first = true;
        long set[14];
        long* it;
        for (long p = 2; p < 100; (p == 2) ? (p++) : (p += 2)) {
                // Test if p is not prime
                if (lookup_table.lt[fmt_p(p)][0][0] == 0) {
                        continue;
                }
                // Test if -3 is a quadratic residue mod p
                if (!has_qi_homomorphism(p)) {
			continue;
		}
                long* cs = x.reductions(p);
                if (first) {
                        if (cs[0] != -1) {
				for (int i = 0; i < 14; i++) {
					set[i] = lookup_table.lt[fmt_p(p)][cs[0]][i];
				}
				if (cs[1] != -1) {
					// There's another c reduction
					// intersect it
					long* entry = lookup_table.lt[fmt_p(p)][cs[1]];
					it = set_intersection(set, first_zero(set), entry, first_zero(entry), set);
					while (it < set + 14) {
						*it = 0;
						it++;
					}
				}
				first = false;
			} else if (cs[1] != -1) {
				for (int i = 0; i < 14; i++) {
					set[i] = lookup_table.lt[fmt_p(p)][cs[1]][i];
				}
				first = false;
			}
                } else {
			for (int cindex = 0; cindex < 2; cindex++) {
				if (cs[cindex] != -1) {
					// Intersect with cs[cindex] lookup
					long* entry = lookup_table.lt[fmt_p(p)][cs[cindex]];
					it = set_intersection(set, first_zero(set), entry, first_zero(entry), set);
					while (it < set + 14) {
						*it = 0;
						it++;
					}
				}
			}
                }
		free(cs);
                if (!first) {
                        bool has_interesting = false;
                        for (long val : set) {
                                if (val >= g) {
                                        has_interesting = true;
                                        break;
                                }
                        }
                        if (!has_interesting) {
				if (x.numer.a == 35 && x.numer.b == 36 && x.denom.a == 27 && x.denom.b == 0) {
					long* cs2 = x.reductions(p);
					printf("GOT HERE. Extra info p = %ld, reds = %ld, %ld\n", p, cs2[0], cs2[1]);
					free(cs2);
				}
                                return false;
                        }
                }
        }
        return true;
}
