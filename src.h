#ifndef CODE_SRC_H
#define CODE_SRC_H

#include <array>
#include <cassert>
#include <cmath>
#include <gmp.h>
#include <gmpxx.h>
#include <map>
#include <tuple>
#include <vector>

/// \param p
/// \return the msb set in p.
uint32_t msb(const mpz_t p) {
	uint32_t bit = 42;
	for (; bit > 0; bit--) {
		if (mpz_tstbit(p, bit))
			break;
	}

	return bit;
}


struct ConfigPrecompute {
public:
	const uint32_t n, r, s, t, m, threads;
	constexpr ConfigPrecompute(uint32_t n, uint32_t r, uint32_t s, uint32_t t, uint32_t m = 1, uint32_t threads = 1) : n(n), r(r), s(s), t(t), m(m), threads(threads) {}
};

template<const ConfigPrecompute &config>
class Precompute {
private:
	// Type definitions
	using LegInt = uint64_t;
	using ListElement = std::pair<LegInt, mpz_class>;

	// Problem instance variables
	constexpr static uint32_t n = config.n;
	constexpr static uint32_t r = config.r;// bits
	constexpr static uint32_t s = config.s;// size prep
	constexpr static uint32_t t = config.t;// time online
	constexpr static uint32_t m = config.m;// nr keys
	constexpr static uint32_t threads = config.threads;
	constexpr static uint32_t WalkLenOffline = m == 1 ? t / 2 : t / (2 * m);
	constexpr static uint32_t WalkLenOnline = m == 1 ? t : t / m;
	constexpr static double r_multiplier = 2;

public:
	mpz_t p;   // prime
	mpz_t k[m];// secret keys

	// PRNG
	gmp_randstate_t rstate;

	uint64_t steps = 0;   // total number of `f` invocations the programm did
	uint64_t rel_ctr = 0; // total number of realtions found in the `NOPREP` algorithm
	uint64_t too_long = 0;// total number of restarts in the `NOPRPEP` algorithm if one walk was too long.

#ifdef NOPREP
	std::map<LegInt, ListElement> L;
#else
	std::vector<ListElement> L;
#endif

	// constructor
	Precompute() {
#ifndef NOPREP
		L.resize(config.s);
#endif
		gmp_randinit_mt(rstate);

		// Init the rng
		srand((r + s + t + n) * time(NULL));
		for (int i = 0; i < rand() % 1000; ++i) {
			mpz_urandomb(p, rstate, r);
		}

		mpz_init(p);
		for (int i = 0; i < m; ++i) {
			mpz_init(k[i]);
		}

		// Implementation bound
		static_assert(r_multiplier * r <= 64);
		steps = 0;
	}

	~Precompute() {
		gmp_randclear(rstate);
		mpz_clear(p);
		for (int i = 0; i < m; ++i) {
			mpz_clear(k[i]);
		}
	}

	void f(mpz_t out, LegInt in) {
		mpz_set_ui(out, in);
		mpz_mod(out, out, p);
	}

	LegInt LegendreSequence(mpz_t in) {
		mpz_t tmp;
		mpz_init_set(tmp, in);
		LegInt seq = 0;
		for (uint32_t i = 0; i < r_multiplier * r; ++i) {
			unsigned char bit = (1 - mpz_legendre(tmp, p)) >> 1;
			seq |= (LegInt(bit) << i);
			mpz_add_ui(tmp, tmp, 1);
		}
		mpz_clear(tmp);
		return seq;
	}

	// Oracle access
	LegInt KeyedLegendreSequence(mpz_t in, const uint32_t l) {
		assert(l < m);

		mpz_t tmp;
		mpz_init_set(tmp, in);
		mpz_add(tmp, tmp, k[l]);
		mpz_mod(tmp, tmp, p);
		LegInt as = LegendreSequence(tmp);
		mpz_clear(tmp);
		return as;
	}

	void sort() {
		std::sort(L.begin(), L.end(),
		          [](const auto &e1, const auto &e2) {
			          return e1.first < e2.first;
		          });
	}

	uint64_t search(LegInt data) {
		static const mpz_class zero = 0;
		ListElement data_target(data, zero);

		auto res = std::lower_bound(L.begin(), L.end(), data_target,
		                            [](const auto &e1, const auto &e2) {
			                            return e1.first < e2.first;
		                            });

		if (res == L.end()) return -1;                // nothing found
		const uint64_t pos = distance(L.begin(), res);// calc the distance
		if (L[pos].first != data) return -1;
		assert(pos < s);

		return pos;
	}

	// compute a random walk W
	void offline_random_walk(const uint32_t i) {
		assert(i < s);

		LegInt as;
		mpz_t y_i, tmp;
		mpz_init(y_i);
		mpz_init(tmp);
		mpz_urandomb(y_i, rstate, r);
		mpz_mod(y_i, y_i, p);

		for (uint32_t j = 0; j < WalkLenOffline; ++j) {
			mpz_set(tmp, y_i);
			as = LegendreSequence(y_i);
			f(y_i, as);
			mpz_add(y_i, y_i, tmp);
			mpz_mod(y_i, y_i, p);
		}

		as = LegendreSequence(y_i);
		L[i] = ListElement(as, y_i);

		mpz_clear(y_i);
		mpz_clear(tmp);
	}

	// compute a random walk W
	bool online_random_walk(mpz_class &possible_k, const uint32_t i) {
		assert(i < m);

		mpz_t x_i, tmp;
		mpz_init(x_i);
		mpz_init(tmp);
		mpz_urandomb(x_i, rstate, r);
		mpz_mod(x_i, x_i, p);
		LegInt as = KeyedLegendreSequence(x_i, i);

		bool ret = false;
		for (uint32_t j = 0; j < WalkLenOnline; ++j) {
			uint64_t pos = search(as);
			if (pos != uint64_t(-1)) {
				possible_k = L[pos].second - mpz_class(x_i);
				ret = true;
				break;
			}

			mpz_set(tmp, x_i);
			f(x_i, as);
			mpz_add(x_i, x_i, tmp);
			mpz_mod(x_i, x_i, p);
			as = KeyedLegendreSequence(x_i, i);

			steps += 1;
		}

		mpz_clear(x_i);
		mpz_clear(tmp);
		return ret;
	}

	uint32_t check_key(const mpz_class &testk) {
		for (uint32_t i = 0; i < m; ++i) {
			if (testk == mpz_class(k[i])) {
				return i;
			}
		}

		return -1;
	}

	void generate_instance() {
		steps = 0;
		rel_ctr = 0;

		mpz_t tmp;
		mpz_init(tmp);

		do {
			mpz_rrandomb(tmp, rstate, r + 1);
			mpz_nextprime(p, tmp);
		} while (msb(p) != r);

		for (int i = 0; i < m; ++i) {
			mpz_urandomb(k[i], rstate, r);
			mpz_mod(k[i], k[i], p);
		}

		mpz_clear(tmp);
	}

	void listInfo() {
		for (int i = 0; i < s; ++i) {
			std::cout << "L(x): " << L[i].first << ", x: " << L[i].second << "\n";
		}
	}

	void info() {
		std::cout << "LEG Informations\n";
		std::cout << "r: " << this->r << ", s: " << s << ", t: " << t << ", m: " << m << ", threads: " << threads
		          << ", WalkLenOffline: " << WalkLenOffline << ", WalkLenOnline: " << WalkLenOnline << "\n";
		std::cout << "p: " << p << "\n";
		for (int i = 0; i < m; ++i) {
			std::cout << "k[" << i << "]: " << k[i] << "\n";
		}
	}

	/// TODO explain
	/// \param maxloops
	/// \return
	uint64_t _attack(const uint64_t maxloops) {
		uint64_t loops = 1;
		mpz_class _k;

		for (int i = 0; i < m; ++i) {
			while (true) {
				if (online_random_walk(_k, i)) {
					if (check_key(_k) != uint32_t(-1)) {
						break;
					}
				}

				loops += 1;

				if (loops > maxloops) {
					loops = -1;
					break;
				}
			}
		}

		return loops;
	}

	/// TODO explain
	/// \return
	uint64_t attack() {
		steps = 0;
		for (uint32_t j = 0; j < s; ++j) {
			offline_random_walk(j);
		}
		sort();


		return _attack(-1);
	}

	bool is_distinguished_point(const uint64_t a, const uint64_t mask_kD) {
		if ((a & mask_kD) == mask_kD)
			return true;
		return false;
	}

#ifdef NOPREP
	/// TODO explain
	/// \return
	uint64_t noprep_attack() {
		// constant algorithm factors
		mpz_t tmpt;
		mpz_init_set(tmpt, p);
		mpz_mul_ui(tmpt, tmpt, m);
		mpz_root(tmpt, tmpt, 2);
		const uint64_t tt = mpz_get_ui(tmpt);
		const double q = m / double(5 * tt);
		const uint64_t kD = ceil(log2(1. / q));
		const uint64_t mask_kD = (uint64_t(1) << kD) - 1;
		const uint64_t len = 8 * tt / m;
		const uint64_t reps = 4;

		//std::cout << "t: " << tt << ", q: " << q << ", k: " << kD << ", mkD: " << mask_kD << ", wlen: " << len << "\n";
		mpz_class _k;

		auto checkrelation = [this](const mpz_class &x1, const mpz_class &x2, const uint64_t i, const uint64_t j) -> bool {
			// checks if k[i] + x1 = k[j] + x2 %p
			// the only function with direct access to the secrets
			return (((mpz_class(k[i]) + x1) % mpz_class(p)) == ((mpz_class(k[j]) + x2) % mpz_class(p)));
		};

		// variables
		LegInt as;
		mpz_t x_l, tmp;
		mpz_init(x_l);
		mpz_init(tmp);

		steps = 0;
		rel_ctr = 0;
		too_long = 0;
		for (int l = 0; l < m; ++l) {
			for (int i = 0; i < reps; ++i) {
			noprep_restart:
				mpz_urandomb(x_l, rstate, r);
				mpz_mod(x_l, x_l, p);

				// random walk
				for (int j = 0; j < len; ++j) {
					mpz_set(tmp, x_l);
					as = KeyedLegendreSequence(x_l, l);
					steps += 1;
					if (is_distinguished_point(as, mask_kD)) {
						if (L.find(as) == L.end()) {
							//std::cout << "walk: " << (i*l)+i << " inserted disinguished point\n";
							L[as] = ListElement(l, x_l);
							goto relation_end;
						} else {
							if (checkrelation(mpz_class(x_l), L[as].second, l, L[as].first)) {
								//std::cout << "found relation: " << l << "," << L[as].first << "\n";
								rel_ctr += 1;
								goto relation_end;
							}
						}
						break;
					}

					f(x_l, as);
					mpz_add(x_l, x_l, tmp);
					mpz_mod(x_l, x_l, p);
				}
				too_long += 1;
				goto noprep_restart;
			relation_end:;
			}
		}

		//std ::cout << "found: " << rel_ctr << " relations\n";
		mpz_clear(x_l);
		mpz_clear(tmp);
		mpz_clear(tmpt);
		return 1;
	}
#endif
};


struct ConfigPrecomputeDLog {
public:
	const uint32_t n, r, s, t, m, threads;
	constexpr ConfigPrecomputeDLog(uint32_t n, uint32_t r, uint32_t s, uint32_t t, uint32_t m = 1, uint32_t threads = 1) : n(n), r(r), s(s), t(t), m(m), threads(threads) {}
};

template<const ConfigPrecomputeDLog &config>
class PrecomputeDlog {
private:
	// Problem instance variables
	constexpr static uint32_t n = config.n;
	constexpr static uint32_t r = config.r;
	constexpr static uint32_t s = config.s;
	constexpr static uint32_t t = config.t;
	constexpr static uint32_t m = config.m;
	constexpr static uint32_t threads = config.threads;
	constexpr static uint32_t WalkLenOffline = m == 1 ? t / 2 : t / (2 * m);
	constexpr static uint32_t WalkLenOnline = m == 1 ? t : t / m;

public:
	mpz_t q;   // safe prime bla blu bli
	mpz_t p;   // prime
	mpz_t p1;  // order
	mpz_t g;   // generator
	mpz_t k[m];// secret exponents
	mpz_t h[m];// g^k

	gmp_randstate_t rstate;// PRNG

	uint64_t steps = 0;
	uint64_t rel_ctr = 0;
	uint64_t too_long = 0;

#ifdef NOPREP
	//                              g_l         l           a_l
	using ListElement = std::tuple<mpz_class, uint64_t, mpz_class>;
	std::vector<ListElement> L;
#else
	using ListElement = std::pair<mpz_class, mpz_class>;
	std::vector<ListElement> L;
#endif

	PrecomputeDlog() {
#ifndef NOPREP
		L.resize(config.s);
#endif
		gmp_randinit_mt(rstate);
		mpz_init(p);

		// Init the rng
		srand((r + s + t + n) * time(NULL));
		for (int i = 0; i < rand() % 1000; ++i) {
			mpz_urandomb(p, rstate, r);
		}

		mpz_init(q);
		mpz_init(p1);
		mpz_init(g);
		for (int i = 0; i < m; ++i) {
			mpz_init(k[i]);
		}

		steps = 0;
	}

	~PrecomputeDlog() {
		gmp_randclear(rstate);
		mpz_clear(q);
		mpz_clear(p);
		mpz_clear(p1);
		mpz_clear(g);
		for (int i = 0; i < m; ++i) {
			mpz_clear(k[i]);
			mpz_clear(h[i]);
		}
	}

	void f(mpz_t out, const mpz_t in) {
		mpz_mod(out, in, p1);
	}


	void sort() {
		std::sort(L.begin(), L.end(),
		          [](const auto &e1, const auto &e2) {
			          return e1.first.get_ui() < e2.first.get_ui();
		          });
	}

	void generate_instance() {
		steps = 0;
		rel_ctr = 0;

		mpz_t tmp;
		mpz_init(tmp);


#ifdef SECUREPRIMEGEN
		mpz_rrandomb(q, rstate, r - 1);

		int is_prime;
		do {
			mpz_nextprime(q, q);

			mpz_mod_ui(tmp, q, 3);
			if (mpz_cmp_ui(tmp, 2) != 0)
				continue;

			mpz_mul_ui(p, q, 2);
			mpz_add_ui(p, p, 1);

			mpz_mod_ui(tmp, p, 3);
			if (mpz_cmp_ui(tmp, 2) != 0)
				continue;

			is_prime = mpz_probab_prime_p(p, 50);
			if (is_prime != 2)
				continue;
		} while (msb(p) < r);
#else
		mpz_rrandomb(q, rstate, r + 1);
		mpz_nextprime(p, q);
#endif
		mpz_set_ui(g, 2);

		mpz_set(p1, p);
		mpz_sub_ui(p1, p1, 1);

		for (uint32_t i = 0; i < m; ++i) {
			mpz_urandomb(k[i], rstate, r);
			mpz_mod(k[i], k[i], p1);
			mpz_powm(h[i], g, k[i], p);
		}

		mpz_clear(tmp);
	}

	void info() {
		std::cout << "DLOG Informations\n";
		std::cout << "r: " << r << ", s: " << s << ", t: " << t << ", m: " << m << ", threads: " << threads
		          << ", WalkLenOffline: " << WalkLenOffline << ", WalkLenOnline: " << WalkLenOnline << "\n";
		std::cout << "p: " << p << ", p-1: " << p1 << ", g: " << g << "\n";
		for (int i = 0; i < m; ++i) {
			std::cout << i << ": g**k_i=" << h[i] << ", k= " << k[i] << "\n";
		}
	}

#ifndef NOPREP

	uint64_t search(const mpz_t data) {
		static const mpz_class zero = 0;
		ListElement data_target(mpz_class(data), zero);

		auto res = std::lower_bound(L.begin(), L.end(), data_target,
		                            [](const auto &e1, const auto &e2) {
			                            return e1.first.get_ui() < e2.first.get_ui();
		                            });

		if (res == L.end()) return -1;                // nothing found
		const uint64_t pos = distance(L.begin(), res);// calc the distance
		if (L[pos].first != mpz_class(data)) return -1;
		assert(pos < s);

		return pos;
	}

	// compute a random walk W
	void offline_random_walk(const uint32_t i) {
		assert(i < s);

		mpz_t y_i, tmp;
		mpz_init(y_i);
		mpz_init(tmp);

		mpz_urandomb(y_i, rstate, r);
		mpz_mod(y_i, y_i, p1);

		for (uint64_t j = 0; j < WalkLenOffline; ++j) {
			mpz_powm(tmp, g, y_i, p);
			f(y_i, tmp);
		}

		mpz_powm(tmp, g, y_i, p);
		L[i] = ListElement(tmp, y_i);

		mpz_clear(y_i);
		mpz_clear(tmp);
	}

	// compute a random walk W
	bool online_random_walk(mpz_class &_k, const uint32_t i) {
		assert(i < m);

		mpz_t x_i, a_i, h_i, tmpg;
		mpz_init(x_i);
		mpz_init(a_i);
		mpz_init(h_i);
		mpz_init(tmpg);
		mpz_urandomb(a_i, rstate, r);
		mpz_mod(a_i, a_i, p1);
		mpz_set(h_i, h[i]);

		bool ret = false;
		for (uint32_t j = 0; j < WalkLenOnline; ++j) {
			mpz_powm(tmpg, g, a_i, p);
			mpz_mul(x_i, tmpg, h_i);

			mpz_mod(x_i, x_i, p);
			// search
			uint64_t pos = search(x_i);
			if (pos != uint64_t(-1)) {
				_k = L[pos].second - mpz_class(a_i);

				ret = true;
				break;
			}

			f(a_i, x_i);
			steps += 1;
		}

		mpz_clear(x_i);
		mpz_clear(a_i);
		mpz_clear(h_i);
		mpz_clear(tmpg);
		return ret;
	}

	// only function with access to the secret
	uint32_t check_key(const mpz_class &testk) {
		for (uint32_t i = 0; i < m; ++i) {
			if (testk == mpz_class(k[i])) {
				return i;
			}
		}

		return -1;
	}

	uint64_t attack() {
		steps = 0;
		mpz_class _k;
		for (int j = 0; j < s; ++j) {
			offline_random_walk(j);
		}
		sort();

		uint64_t loops = 1;
		for (uint32_t i = 0; i < m; ++i) {
			while (true) {
				if (online_random_walk(_k, i)) {
					if (check_key(_k) != uint32_t(-1)) {
						break;
					}
				}

				loops += 1;

				if (loops > uint64_t(-1)) {
					std::cout << "?";
					loops = -1;
					break;
				}
			}
		}

		return loops;
	}


#else
	bool is_distinguished_point(const mpz_t a, const uint64_t mask_kD) {
		if ((mpz_get_ui(a) & mask_kD) == mask_kD)
			return true;
		return false;
	}

	/// TODO explain
	/// \return
	uint64_t noprep_attack() {
		// constant algorithm factors
		mpz_t tmpt;
		mpz_init_set(tmpt, p);
		mpz_mul_ui(tmpt, tmpt, m);
		mpz_root(tmpt, tmpt, 2);
		const uint64_t tt = mpz_get_ui(tmpt);
		const double _q = m / double(5 * tt);
		const uint64_t kD = ceil(log2(1. / _q));
		const uint64_t mask_kD = (uint64_t(1) << kD) - 1;
		const uint64_t len = 8 * tt / m;
		const uint64_t reps = 4;

		//std::cout << "t: " << tt << ", q: " << _q << ", k: " << kD << ", mkD: " << mask_kD << ", wlen: " << len << "\n";

		mpz_class _k, _x_l;

		auto checkrelation = [this](const mpz_class &x1, const mpz_class &x2, const uint64_t i, const uint64_t j) -> bool {
			// checks if k[i] + x1 = k[j] + x2 %p
			// the only function with direct access to the secrets
			return (((mpz_class(k[i]) + x1) % mpz_class(p)) == ((mpz_class(k[j]) + x2) % mpz_class(p)));
		};

		// variables
		mpz_t x_l, a_l, tmpg, h_l;
		mpz_init(x_l);
		mpz_init(tmpg);
		mpz_init(a_l);
		mpz_init(h_l);
		steps = 0;
		rel_ctr = 0;
		too_long = 0;
		for (int l = 0; l < m; ++l) {
			mpz_set(h_l, h[l]);

			for (int i = 0; i < reps; ++i) {
			noprep_restart:
				mpz_urandomb(a_l, rstate, r);
				mpz_mod(a_l, a_l, p);

				// random walk
				for (int j = 0; j < len; ++j) {
					steps += 1;
					mpz_powm(tmpg, g, a_l, p);
					mpz_mul(x_l, tmpg, h_l);
					mpz_mod(x_l, x_l, p);

					if (is_distinguished_point(x_l, mask_kD)) {
						_x_l = mpz_class(x_l);

						// add to list of distinguished points
						auto it = std::find_if(L.begin(), L.end(),
						                       [_x_l](const auto &e1) {
							                       return std::get<0>(e1) == _x_l;
						                       });

						if (it == L.end()) {
							// not found
							L.template emplace_back(ListElement(_x_l, l, a_l));
							goto relation_end;
						} else {
							if (checkrelation(mpz_class(a_l), std::get<2>(*it), l, std::get<1>(*it))) {
								//std::cout << "found relation: " << l << "," << std::get<1>(*it) << "\n";
								rel_ctr += 1;
								goto relation_end;
							}
						}

						break;
					}


					f(a_l, x_l);
				}
				//std::cout << "too long\n";
				too_long += 1;
				goto noprep_restart;
			relation_end:;
			}
		}

		//std ::cout << "found: " << rel_ctr << " relations\n";
		mpz_clear(x_l);
		mpz_clear(a_l);
		mpz_clear(tmpg);
		mpz_clear(h_l);
		mpz_clear(tmpt);

		/*std::sort(L.begin(), L.end(),
		          [](const auto &e1, const auto &e2) {
			          return std::get<1>(e1) < std::get<1>(e2) ;
		          }
		);

		for (int i = 0; i < L.size(); ++i) {
			std::cout << std::get<0>(L[i]) << " " << std::get<1>(L[i]) << " " << std::get<2>(L[i]) << "\n";
		}*/

		return 1;
	}
#endif
};


#endif//CODE_SRC_H
