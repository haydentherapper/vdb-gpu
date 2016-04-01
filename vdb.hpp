/**
 * \file vdb.hpp
 * \brief
 *
 * \author Hayden Blauzvern
 */


#ifndef VDB_HPP_INCLUDED
#define VDB_HPP_INCLUDED

#include <cmath>

#include "coord.hpp"

union InternalData {
    uint64_t index; //
    double tile_or_value; // tile or Value
};

template<uint64_t level0, uint64_t level1, uint64_t level2>
class VDB {
public:
    enum InternalNodeLevel { level2_node, level1_node };

    const size_t total_storage_size_;
    InternalData* vdb_storage_;

    /// Leaf node single-axis  log size (x, y, and z)
    const uint64_t LEVEL0_LOG2 = level0;
    const uint64_t LEVEL0_sLOG2 = level0;
    /// Leaf node voxel count
    const uint64_t LEVEL0_sSIZE = 1 << (LEVEL0_LOG2 * 3);

    /// Internal node (level 1) single-axis log size
    const uint64_t LEVEL1_LOG2 = level1;
    /// Internal node (level 1) log voxel-count
    const uint64_t LEVEL1_sLOG2 = LEVEL1_LOG2 + LEVEL0_sLOG2;
    /// Internal node (level 1) grid size
    const uint64_t LEVEL1_sSIZE = 1 << (LEVEL1_LOG2 * 3);

    /// Internal node (level 2) single-axis log size
    const uint64_t LEVEL2_LOG2 = level2;
    /// Internal node (level 2) log voxel-count
    const uint64_t LEVEL2_sLOG2 = LEVEL2_LOG2 + LEVEL1_sLOG2;
    /// Internal node (level 2) grid size
    const uint64_t LEVEL2_sSIZE = 1 << (LEVEL2_LOG2 * 3);

    /// Size of InternalData in bits
    const uint64_t INTERNAL_DATA_SIZE = sizeof(InternalData) * 8;
    /// Size of level 2 child and value mask
    const uint64_t LEVEL2_MASK_SIZE = LEVEL2_sSIZE / INTERNAL_DATA_SIZE;
    /// Size of level 1 child and value mask
    const uint64_t LEVEL1_MASK_SIZE = LEVEL1_sSIZE / INTERNAL_DATA_SIZE;
    /// Size of leaf node value mask
    const uint64_t LEVEL0_MASK_SIZE = LEVEL0_sSIZE / INTERNAL_DATA_SIZE;
    /// Total size of level 2 internal node
    const uint64_t LEVEL2_TOTAL_SIZE = LEVEL2_sSIZE + LEVEL2_MASK_SIZE * 2;
    /// Total size of level 1 internal node
    const uint64_t LEVEL1_TOTAL_SIZE = LEVEL1_sSIZE + LEVEL1_MASK_SIZE * 2;
    /// Total size of leaf node
    const uint64_t LEVEL0_TOTAL_SIZE = LEVEL0_sSIZE + LEVEL0_MASK_SIZE;

    /// Hashmap size
    const uint64_t HASHMAP_SIZE;
    const uint64_t HASHMAP_LOG_SIZE = uint64_t(log2(HASHMAP_SIZE));
    /// Index into vdb_storage where hashmap starts
    const uint64_t HASHMAP_START = 0;

    const double BACKGROUND_VALUE;

    uint64_t num_elements_;

    /// Disabled default constructor
    VDB() = delete;

    /// Constructor with memory size parameter
    VDB(uint64_t mem_size, uint64_t hashmap_size, double background);

    /// Destructor
    ~VDB();

    /// Disabled copy constructor
    VDB(VDB& other) = delete;

    /// Disabled assignment operator.
    VDB& operator=(VDB& rhs) = delete;

    /// Random access into VDB
    double random_access(const Coord xyz);

    /// Random insert into VDB
    bool random_insert(const Coord xyz, double value);

private:
    /// Initializes hashmap with uint64_t max
    void initialize_hashmap();

    /// Initializes vdb_storage_ with 0
    void initialize_vdb_storage();

    /// Gets index to internal node array from hashmap.
    /// Will return null if no node exists at coordinates
    InternalData get_node_from_hashmap(const Coord xyz);

    bool insert_node_into_hashmap(const Coord xyz, uint64_t node_index);

    void insert_internal_node(uint64_t index, InternalNodeLevel inl);

    void insert_leaf_node(uint64_t index);

    /// Returns index into internal node array
    uint64_t calculate_internal_offset(const Coord xyz, InternalNodeLevel inl);

    /// Returns index into leaf node array
    uint64_t calculate_leaf_offset(const Coord xyz);
};

#include "vdb-private.hpp"

#endif /* end of include guard: VDB_HPP_INCLUDED */