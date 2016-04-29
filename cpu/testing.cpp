#include <iostream>
#include <random>
#include "vdb.hpp"
//#include <omp.h>

int main() {
    // Init RNG
    // From http://stackoverflow.com/questions/21102105/random-double-c11
    const double lower_bound = -255;
    const double upper_bound = 255;
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    std::random_device rand_dev;
    std::mt19937 rand_engine(rand_dev());

    VDB<3, 4, 5> vdb_square(1000 * 256 * 256, 2000, 10.0);
    // omp_set_num_threads(24);
    //#pragma omp parallel for collapse(3) shared(vdb_square)
    for (size_t i = 1; i < 20; ++i) {
        for (size_t j = 1; j < 24; ++j) {
            for (size_t k = 0; k < 23; ++k) {
                const Coord coord(i, j, k);
                double val = dist(rand_engine);
                vdb_square.random_insert(coord, val);
                double result = vdb_square.random_access(coord);
                if (val != result) {
                    std::cout << i << ' ' << j << ' ' << k << ' ' << result
                              << std::endl;
                }
            }
        }
    }
}
