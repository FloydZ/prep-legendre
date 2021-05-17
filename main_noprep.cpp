#include <iostream>

#define NOPREP
#include "main.h"

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

#include "src.h"

static constexpr ConfigPrecompute configLegendre(n, r, s, t, m, threads);
static constexpr ConfigPrecomputeDLog configDLog(n, r, s, t, m, threads);
Precompute<configLegendre> leg;
PrecomputeDlog<configDLog> dlog;

void LegendreAttack() {
	double time = 0;
	uint64_t loops = 0, loop, steps = 0, relations=0, too_longs=0;

	for (uint32_t i = 0; i < BENCHES; i++) {
		leg.generate_instance();

		double t0 = ((double)clock()/CLOCKS_PER_SEC);
		loop = leg.noprep_attack();
		if (loop == -1)
			continue;

		time += ((double)clock()/CLOCKS_PER_SEC) - t0;
		loops += loop;
		steps += leg.steps;
		relations += leg.rel_ctr;
		too_longs += leg.too_long;
	}

	std::cout << "l," << r << "," << t << "," << s << "," << m << ","  << time << "," << loops << "," << steps  << "," << relations << "," << too_longs << "," << BENCHES << "\n" << std::flush;
}
void DLogAttack() {
	double time = 0;
	uint64_t loops = 0, loop, steps = 0, relations=0, too_longs=0;

	for (uint32_t i = 0; i < BENCHES; i++) {
		dlog.generate_instance();

		double t0 = ((double)clock()/CLOCKS_PER_SEC);
		loop = dlog.noprep_attack();
		if (loop == -1)
			continue;

		time += ((double)clock()/CLOCKS_PER_SEC) - t0;
		loops += loop;
		steps += dlog.steps;
		relations += dlog.rel_ctr;
		too_longs += dlog.too_long;
	}

	std::cout << "d," << r << "," << t << "," << s << "," << m << ","  << time << "," << loops << "," << steps  << "," << relations << "," << too_longs << "," << BENCHES << "\n" << std::flush;
}

int main() {
	LegendreAttack();
	DLogAttack();
	return 0;
}
