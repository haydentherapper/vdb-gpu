#include <iostream>
#include <random>
#include "vdb.hpp"

#include <Kokkos_Core.hpp>

KOKKOS_INLINE_FUNCTION
void
convertSingleIndexToThree(const size_t pointIndex,
                          const size_t N_J,
                          const size_t N_K,
                          size_t* i,
                          size_t* j,
                          size_t* k) {
  *i = pointIndex / (N_J * N_K);
  *j = (pointIndex - *i * N_J * N_K) / N_K;
  *k = pointIndex % N_K;
}

int main(int argc, char** argv) {
    Kokkos::initialize(argc, argv);

    const size_t VDB_SIZE = 1000 * 256 * 256;
    const size_t VDB_HASH_MAP_SIZE = 2000;
    const double VDB_BACKGROUND_VALUE = 10.0;
    VDB<3, 4, 5> vdb_square(VDB_SIZE, VDB_HASH_MAP_SIZE, VDB_BACKGROUND_VALUE);

    const size_t N_I = 20;
    const size_t N_J = 21;
    const size_t N_K = 22;
    Kokkos::parallel_for(
        N_I * N_J * N_K,
        KOKKOS_LAMBDA (const size_t pointIndex) {
          size_t i;
          size_t j;
          size_t k;
          convertSingleIndexToThree(pointIndex, N_J, N_K,
                                    &i, &j, &k);
          const Coord coord(i, j, k);
          double val = i * 10000 + j * 100 + k;
          vdb_square.random_insert(coord, val);
          double result = vdb_square.random_access(coord);
          if (val != result) {
              printf("discrepancy at (%2zu, %2zu, %2zu), "
                     "expected %9.4lf but got %9.4lf\n",
                     i, j, k, val, result);
          }
        });

    printf("finished, finalizing, ignore the strange error message "
           "after this.\n");
    Kokkos::finalize();
}
