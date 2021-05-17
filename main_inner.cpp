#include <iostream>
#include "main.h"
#include "src.h"

#ifndef CONFIG_SET
constexpr uint32_t BENCHES = 5;
constexpr uint32_t BENCHESMULT = 100000;

constexpr uint32_t n = 20;
constexpr uint32_t r = 30;
constexpr uint32_t m = 2;
constexpr uint32_t threads = 1;

// standard config.
// constexpr uint32_t s = 100;
// constexpr uint32_t t = 40;
// IMPORTANT. Make sure r is dividable by 3
constexpr uint32_t s = uint32_t(1) << r/3;
constexpr uint32_t t = uint32_t(1) << r/3;
#endif

static constexpr ConfigPrecompute configLegendre(n, r, s, t, m, threads);
static constexpr ConfigPrecomputeDLog configDLog(n, r, s, t, m, threads);
Precompute<configLegendre> leg;
PrecomputeDlog<configDLog> dlog;

void InnerFuntion() {
	double ltime = 0., dtime = 0.;
#ifndef BENCHMODE
	leg.info();
	dlog.info();
#endif

	leg.generate_instance ();
	dlog.generate_instance();

	// Some helper functions
	mpz_t y, tmp;
	mpz_init(y); mpz_init(tmp);
	gmp_randstate_t rstate;
	gmp_randinit_mt(rstate);
	mpz_urandomb(y, rstate, r);
	uint64_t as = leg.KeyedLegendreSequence(y, 0);

	for (uint32_t i = 0; i < BENCHESMULT; i++) {
		double t0 = ((double)clock()/CLOCKS_PER_SEC);

		mpz_set(tmp, y);
		leg.f(y, as);
		mpz_add(y, y, tmp);
		mpz_mod(y, y, leg.p);
		as = leg.KeyedLegendreSequence(y, 0);

		ltime += ((double)clock()/CLOCKS_PER_SEC) - t0;
	}

	for (uint32_t i = 0; i < BENCHESMULT; i++) {
		double t1 = ((double)clock()/CLOCKS_PER_SEC);
		mpz_powm(tmp, dlog.g, y, dlog.p);
		mpz_mul(y, tmp, dlog.h[0]);
		mpz_mod(y, y, dlog.p);
		dlog.f(tmp, y);

		dtime += ((double)clock()/CLOCKS_PER_SEC) - t1;
	}

#ifdef BENCHMODE
	std::cout << r << "," << BENCHESMULT << ","  << ltime << "," << dtime << "\n";
#else
	std::cout << "leg Time: " << ltime/BENCHESMULT << "\n";
	std::cout << "log Time: " << dtime/BENCHESMULT << "\n";
	std::cout << "coeff   : " << ltime/dtime << " leg slower than log\n";
#endif
	mpz_clear(y); mpz_clear(tmp);
	gmp_randclear(rstate);
}

int main() {
	InnerFuntion();
	return 0;
}
