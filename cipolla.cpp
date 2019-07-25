
#include "utils.cpp"

__device__ long uint64_log2(uint64_t n)
{
  //#define S(k) if (n >= (UINT64_C(1) << k)) { i += k; n >>= k; }

  //long i = -(n == 0); S(32); S(16); S(8); S(4); S(2); S(1); return i+1;

  //#undef S
  return 64 - __clzll(n) - 1;
}

__device__ long mod_power(long a, long b, long p) {
    long res = 1;
    for (int i = 0; i < b; i++) {
        res *= a;
        res %= p;
    }
    return rem_euclid(res, p);
}

struct cipolla_res {
    bool some;
    long fst, snd;
} CR;

__device__ long cipolla_multa(long a, long b, long c, long d, long w, long p) {
    return (a * c + b * d * w) % p;
}

__device__ long cipolla_multb(long a, long b, long c, long d, long w, long p) {
    return (a * d + b * c) % p;
}

__device__ cipolla_res cipolla(long n, long p, bool debug=false) {
    n = rem_euclid(n, p);
    cipolla_res res;
    if (n == 0 || n == 1) {
        res.some = true;
        res.fst = n;
        res.snd = rem_euclid(-n, p);
        return res;
    }
    long phi = p - 1;
    if (mod_power(n, phi/2, p) != 1) {
        res.some = false;
        return res;
    }
    if (p % 4 == 3) {
        long ans = mod_power(n, (p + 1)/4, p);
        res.some = true;
        res.fst = ans;
        res.snd = rem_euclid(-ans, p);
        return res;
    }
    long aa = 0;
    for (long i = 1; i < p; i++) {
        long temp = mod_power(rem_euclid(i*i - n, p), phi/2, p);
        if (temp == phi) {
            aa = i;
            break;
        }
    }
    long exponent =(p + 1)/2;
    long x1a = aa;
    long x1b = 1;
    long x2a = cipolla_multa(x1a, x1b, x1a, x1b, aa*aa - n, p);
    long x2b = cipolla_multb(x1a, x1b, x1a, x1b, aa*aa - n, p);
    long l = uint64_log2((unsigned long) exponent);
    for (long i = l-1; i >= 0; i--) {
        if (((exponent) & (1 << i)) == 0) {
            long x2a_ = cipolla_multa(x2a, x2b, x1a, x1b, aa*aa - n, p);
            long x2b_ = cipolla_multb(x2a, x2b, x1a, x1b, aa*aa - n, p);
            x2a = x2a_;
            x2b = x2b_;
            long x1a_ = cipolla_multa(x1a, x1b, x1a, x1b, aa*aa - n, p);
            long x1b_ = cipolla_multb(x1a, x1b, x1a, x1b, aa*aa - n, p);
            x1a = x1a_;
            x1b = x1b_;
        } else {
            long x1a_ = cipolla_multa(x1a, x1b, x2a, x2b, aa*aa - n, p);
            long x1b_ = cipolla_multb(x1a, x1b, x2a, x2b, aa*aa - n, p);
            x1a = x1a_;
            x1b = x1b_;
            long x2a_ = cipolla_multa(x2a, x2b, x2a, x2b, aa*aa - n, p);
            long x2b_ = cipolla_multb(x2a, x2b, x2a, x2b, aa*aa - n, p);
            x2a = x2a_;
            x2b = x2b_;
        }
    }
    res.some = true;
    res.fst = x1a;
    res.snd = rem_euclid(-x1a, p);
    return res;
}

__device__ bool has_qw_homomorphism(long p) {
    return p == 3 || mod_power(p-3, (p - 1)/2, p) == 1;
}
