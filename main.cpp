#include "main.h"
#include "src.h"
#include <iostream>

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
constexpr uint32_t s = uint32_t(1) << r / 3;
constexpr uint32_t t = uint32_t(1) << r / 3;
#endif

static constexpr ConfigPrecompute configLegendre(n, r, s, t, m, threads);
static constexpr ConfigPrecomputeDLog configDLog(n, r, s, t, m, threads);
Precompute<configLegendre> leg;
PrecomputeDlog<configDLog> dlog;

void LegendreAttack() {
	leg.generate_instance();
	leg.info();
	leg.attack();
}
void DLogAttack() {
	dlog.generate_instance();
	dlog.info();
	dlog.attack();
}

void LegendreTiming() {
	double time = 0;
	uint64_t loops = 0, loop, steps = 0;
#ifndef BENCHMODE
	leg.info();
#endif

	for (uint32_t i = 0; i < BENCHES; i++) {
		leg.generate_instance();

		double t0 = ((double) clock() / CLOCKS_PER_SEC);
		loop = leg.attack();
		if (loop == -1)
			continue;

		time += ((double) clock() / CLOCKS_PER_SEC) - t0;
		loops += loop;
		steps += leg.steps;

#ifndef BENCHMODE
		std::cout << i << " " << std::flush;
#endif
	}

#ifdef BENCHMODE
	std::cout << "l," << r << "," << t << "," << s << "," << m << "," << time << "," << time / BENCHES << "," << time / loops << "," << loops << "," << loops / BENCHES << "," << steps << "," << steps / loops << "\n"
	          << std::flush;
#else
	std::cout << "\nLegendre Time: " << time / BENCHES << ", time/loops: " << time / loops << ", loops: " << loops / BENCHES << ", steps: " << steps << "\n"
	          << std::flush;
#endif
}

void DlogTiming() {
	double time = 0;
	uint64_t loops = 0, loop, steps = 0;
#ifndef BENCHMODE
	dlog.info();
#endif

	for (uint32_t i = 0; i < BENCHES; i++) {
		dlog.generate_instance();

		double t0 = ((double) clock() / CLOCKS_PER_SEC);
		loop = dlog.attack();
		if (loop == -1)
			continue;

		time += ((double) clock() / CLOCKS_PER_SEC) - t0;
		loops += loop;
		steps += dlog.steps;
#ifndef BENCHMODE
		std::cout << i << " " << std::flush;
#endif
	}
#ifdef BENCHMODE
	std::cout << "d," << r << "," << t << "," << s << "," << m << "," << time << "," << time / BENCHES << "," << time / loops << "," << loops << "," << loops / BENCHES << "," << steps << "," << steps / loops << "\n"
	          << std::flush;
#else
	std::cout << "\nDLOG Time: " << time / BENCHES << ", time/loops: " << time / loops << ", loops: " << loops / BENCHES << ", steps: " << steps << "\n"
	          << std::flush;
#endif
}

void BenchMode() {
	LegendreTiming();
	DlogTiming();
}

int main() {
#ifdef BENCHMODE
	BenchMode();
#else
	//LegendreAttack();
	//DLogAttack();
	//LegendreTiming();
	//DlogTiming();
	InnerFuntion();
#endif
	return 0;
}
