#include "vdb.hpp"
#include <iostream>
#include <random>

int main() {
    // Init RNG
    // From http://stackoverflow.com/questions/21102105/random-double-c11
    const double lower_bound = -255;
    const double upper_bound = 255;
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    std::random_device rand_dev;
    std::mt19937 rand_engine(rand_dev());

    VDB<3, 4, 5> vdb_square(256*256*256, 2000, 10.0);
    // pragma omp parallal for collapse(3)
    for (size_t i = 0; i < 16; ++i) {
        for (size_t j = 0; j < 16; ++j) {
            for (size_t k = 0; k < 16; ++k) {
                const Coord coord(i, j, k);
                double val = dist(rand_engine);
                vdb_square.random_insert(coord, val);
                double result = vdb_square.random_access(coord);
                if (val != result) {
                    std::cout << i << ' ' << j << ' ' << k << ' ' << result << std::endl;
                }
            }
        }
    }
}
