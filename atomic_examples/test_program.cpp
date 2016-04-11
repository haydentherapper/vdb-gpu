#include <cstdio>
#include <cstdint>

// Compile with: g++ -std=c++11 test_program.cpp -fopenmp

int main() {
    union InternalData {
        uint64_t index; //
        double tile_or_value; // tile or Value
    };

    InternalData test;
    test.index = 0;

    #pragma omp parallel for
    for (size_t i = 0; i < 1000000; ++i) {
        __sync_fetch_and_add(reinterpret_cast<uint64_t*>(&test), 1);
    }

    printf("%lu\n", test.index);
}
