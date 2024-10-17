#pragma once

#include <random>
#include <iostream>
#include <iomanip>

class UURandom {
private:
    virtual double uniform_unit_random() = 0;

    unsigned parallel_seed(int myid, unsigned seed) {
        std::mt19937 rng(seed); // Standard Mersenne Twister RNG
        std::uniform_int_distribution<unsigned> dist(0, 0x7fffffff);

        for (int i = 0; i < myid; i++) {
            dist(rng); // Advance the RNG
        }
        return dist(rng); // Return a random seed
    }


public:
    virtual ~UURandom() {}
    virtual void init_random(int seed) = 0;
    void init_random_parallel(int myid, unsigned seed) {
        unsigned seed_parallel = parallel_seed(myid, seed);
        // std::cout << "random seed = " << std::setw(10) << seed_parallel << " with myid = " << myid << std::endl;
        init_random(seed_parallel);
    }

    double operator()() {
        return uniform_unit_random();
    }

    void operator()(VecReal& v) {
        for (int i = 0; i < v.size(); ++i) {
            v[i] = (*this)();
        }
    }

    void operator()(VecCmplx& v) {
        for (int i = 0; i < v.size(); ++i) {
            v[i] = Cmplx((*this)(), (*this)());
        }
    }
};

class UURandomMersenne : public UURandom {
private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;

    double uniform_unit_random() {
        return dist(rng);
    }

public:
    UURandomMersenne(int myid = 0, int seed = 1)
        : rng(seed), dist(0.0, 1.0) {
        init_random_parallel(myid, unsigned(seed));
    }

    void init_random(int seed) {
        rng.seed(seed);
    }
};

class UURandomSFMT : public UURandom {
private:
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> dist;

    double uniform_unit_random() {
        return dist(rng);
    }

public:
    UURandomSFMT(int myid = 0, int seed = 1)
        : rng(seed), dist(0.0, 1.0) {
        init_random_parallel(myid, unsigned(seed));
    }

    void init_random(int seed) {
        rng.seed(seed);
    }
};