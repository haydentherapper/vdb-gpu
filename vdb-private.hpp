/**
 * \file vdb-private.hpp
 * \brief
 *
 * \author Hayden Blauzvern
 */

#ifndef VDB_PRIVATE_HPP_INCLUDED
#define VDB_PRIVATE_HPP_INCLUDED

#include <limits>
#include <array>
#include <iostream>
#include <cstdlib>
#include <cstring> // memset

#include "coord.hpp"
#include "utils.hpp"
#include "vdb.hpp"

template<uint64_t level0, uint64_t level1, uint64_t level2>
VDB<level0, level1, level2>::VDB(uint64_t mem_size,
                                 uint64_t hashmap_size,
                                 double background)
    : total_storage_size_(mem_size),
      vdb_storage_(new InternalData[mem_size]),
      HASHMAP_SIZE(hashmap_size),
      BACKGROUND_VALUE(background),
      num_elements_(HASHMAP_SIZE)
{
    initialize_vdb_storage();
    initialize_hashmap();
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
VDB<level0, level1, level2>::~VDB()
{
    delete [] vdb_storage_;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
void VDB<level0, level1, level2>::initialize_hashmap()
{
    for (uint64_t i = HASHMAP_START; i < HASHMAP_START + HASHMAP_SIZE; ++i) {
        vdb_storage_[i].index = std::numeric_limits<uint64_t>::max();
    }
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
void VDB<level0, level1, level2>::initialize_vdb_storage()
{
    memset(vdb_storage_, 0, total_storage_size_ * sizeof(InternalData));
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
double VDB<level0, level1, level2>::random_access(const Coord xyz)
{
    // Index to internal node or background
    InternalData node_from_hashmap = get_internal_node_from_hashmap(xyz);
    if (node_from_hashmap.tile_or_value == BACKGROUND_VALUE) {
        return node_from_hashmap.tile_or_value;
    }

    bool is_tile_level2(false);
    InternalData internal_node_level2 = get_internal_or_leaf_node(xyz, node_from_hashmap.index, &is_tile_level2, level2_node);
    if (internal_node_level2.tile_or_value == BACKGROUND_VALUE || is_tile_level2) {
        return internal_node_level2.tile_or_value;
    }

    bool is_tile_level1(false);
    InternalData internal_node_level1 = get_internal_or_leaf_node(xyz, internal_node_level2.index, &is_tile_level1, level1_node);
    if (internal_node_level1.tile_or_value == BACKGROUND_VALUE || is_tile_level1) {
        return internal_node_level1.tile_or_value;
    }

    double value = get_value_from_leaf_node(xyz, internal_node_level1.index);
    return value;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
bool VDB<level0, level1, level2>::random_insert(const Coord xyz, double value)
{
    // First, get index to internal node or background
    InternalData node_from_hashmap = get_internal_node_from_hashmap(xyz);
    if (node_from_hashmap.tile_or_value == BACKGROUND_VALUE) {
        // Node does not exist, insert into hashmap to update parent
        uint64_t result = insert_internal_node_into_hashmap(xyz);
        if (result == 0) {
            std::cout << "Not enough room in hashtable, exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        node_from_hashmap.index = result;
    }

    // Insert for level 2 internal node
    InternalData level2_index;
    level2_index.index = insert_internal_or_leaf_node(xyz, node_from_hashmap.index, level2_node);

    // Recurse for level 1 internal node
    InternalData level1_index;
    level1_index.index = insert_internal_or_leaf_node(xyz, level2_index.index, level1_node);

    // Set value in leaf node
    insert_value_into_leaf_node(xyz, level1_index.index, value);

    return true;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
InternalData VDB<level0, level1, level2>::get_internal_node_from_hashmap(const Coord xyz)
{
    std::array<uint64_t, 3> rootkey = xyz.get_rootkey(LEVEL2_sLOG2);
    uint64_t hash = Coord::hash(rootkey, HASHMAP_LOG_SIZE);
    uint64_t compressed_xyz = xyz.get_compressed_coord(LEVEL2_sLOG2);
    uint64_t size_of_axis = 20 - LEVEL2_sLOG2; // Likely 8
    uint64_t amount_to_pad = 64 - (size_of_axis * 3); // Likely 40
    for (uint64_t i = HASHMAP_START; i < HASHMAP_START + HASHMAP_SIZE; ++i) {
        uint64_t index = (hash + i * i) % HASHMAP_SIZE;
        if (vdb_storage_[index].index != std::numeric_limits<uint64_t>::max()) {
            // Compare to compressed_xyz to find correct hashmap value
            uint64_t compressed_xyz_index = vdb_storage_[index].index;
            if (((compressed_xyz_index >> amount_to_pad) << amount_to_pad) ==
                ((compressed_xyz >> amount_to_pad) << amount_to_pad)) {
                // TODO: Fix ugly hack
                uint64_t mask = (1LLU << amount_to_pad) - 1;
                while((vdb_storage_[index].index & mask) == mask);

                // Extract last 40 bits
                InternalData ret;
                ret.index = (vdb_storage_[index].index) & ((1LLU << amount_to_pad) - 1);
                return ret;
            }
        }
    }
    InternalData ret;
    ret.tile_or_value = BACKGROUND_VALUE;
    return ret;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
InternalData VDB<level0, level1, level2>::get_internal_or_leaf_node(const Coord xyz, uint64_t index, bool* is_tile, InternalNodeLevel inl)
{
    uint64_t size;
    uint64_t mask_size;
    switch (inl) {
        case level2_node:
            size = LEVEL2_sSIZE;
            mask_size = LEVEL2_MASK_SIZE;
            break;
        case level1_node:
            size = LEVEL1_sSIZE;
            mask_size = LEVEL1_MASK_SIZE;
            break;
    }

    // Pointer to first index into internal data array
    InternalData* internal_node_array_pointer = vdb_storage_ + index;
    uint64_t internal_offset = calculate_internal_offset(xyz, inl);

    // Spin until ready
    InternalData* lock_flag_array = internal_node_array_pointer +
                               // Length of node array
                               size +
                               // Length of value mask
                               mask_size +
                               // Length of child mask
                               mask_size;
    lock_flag_array += internal_offset / NUM_LOCK_FLAGS;
    uint64_t lock_offset = internal_offset % NUM_LOCK_FLAGS;
    while(*(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) != DONE);

    InternalData internal_node_index =
        *(internal_node_array_pointer + internal_offset);

    // Check value mask to see if node exists
    InternalData* value_mask = internal_node_array_pointer +
                               // Length of node array
                               size;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData value_mask_chunk = *value_mask;
    uint64_t value_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(value_mask_chunk.index, value_mask_offset)) {
        // If value_mask not set, return background
        InternalData ret;
        ret.tile_or_value = BACKGROUND_VALUE;
        return ret;
    }

    InternalData* child_mask = internal_node_array_pointer +
                               // Length of node array
                               size +
                               // Length of value mask
                               mask_size;
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData child_mask_chunk = *child_mask;

    uint64_t child_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(child_mask_chunk.index, child_mask_offset)) {
        // Returned value will be interpreted as tile
        *is_tile = true;
    }

    // Return index or tile, set with boolean flag
    return internal_node_index;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
double VDB<level0, level1, level2>::get_value_from_leaf_node(const Coord xyz, uint64_t index)
{
    InternalData* leaf_node_array_pointer = vdb_storage_ + index;
    uint64_t leaf_offset = calculate_leaf_offset(xyz);
    InternalData leaf_node_index = *(leaf_node_array_pointer + leaf_offset);
    // Value mask tracks active states
    InternalData* value_mask = leaf_node_array_pointer +
                               // Length of node array
                               LEVEL0_sSIZE;
    value_mask += leaf_offset / INTERNAL_DATA_SIZE;
    InternalData value_mask_chunk = *value_mask;
    uint64_t value_mask_offset = leaf_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(value_mask_chunk.index, value_mask_offset)) {
        return BACKGROUND_VALUE;
    }
    // Return voxel value
    return leaf_node_index.tile_or_value;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
uint64_t VDB<level0, level1, level2>::insert_internal_or_leaf_node(const Coord xyz, uint64_t index, InternalNodeLevel inl)
{
    uint64_t size;
    uint64_t mask_size;
    uint64_t total_size;
    switch (inl) {
        case level2_node:
            size = LEVEL2_sSIZE;
            mask_size = LEVEL2_MASK_SIZE;
            total_size = LEVEL1_TOTAL_SIZE;
            break;
        case level1_node:
            size = LEVEL1_sSIZE;
            mask_size = LEVEL1_MASK_SIZE;
            total_size = LEVEL0_TOTAL_SIZE;
            break;
    }

    // Pointer to first index into internal data array
    InternalData* internal_node_array_pointer = vdb_storage_ + index;
    uint64_t internal_offset = calculate_internal_offset(xyz, inl);
    InternalData& internal_node_index =
        *(internal_node_array_pointer + internal_offset);
    InternalData* value_mask = internal_node_array_pointer +
                               // Length of node array
                               size;
    InternalData* child_mask = internal_node_array_pointer +
                               // Length of node array
                               size +
                               // Length of value mask
                               mask_size;
   InternalData* lock_flag_array = internal_node_array_pointer +
                              // Length of node array
                              size +
                              // Length of value mask
                              mask_size +
                              // Length of child mask
                              mask_size;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    lock_flag_array += internal_offset / NUM_LOCK_FLAGS;

    InternalData child_mask_chunk = *child_mask;
    uint64_t mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    uint64_t lock_offset = internal_offset % NUM_LOCK_FLAGS;
    if (!extract_bit(child_mask_chunk.index, mask_offset)) {
        // No internal node child, add an internal node
        uint8_t old_lock_val = __atomic_exchange_n(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset, IP, __ATOMIC_SEQ_CST);
        if (old_lock_val == READY) {
            uint64_t old_num_elements = __atomic_fetch_add(&num_elements_, total_size, __ATOMIC_SEQ_CST);
            internal_node_index.index = old_num_elements;
            uint64_t constant = 1LLU << mask_offset;
            __atomic_fetch_or(reinterpret_cast<uint64_t*>(value_mask), constant, __ATOMIC_SEQ_CST);
            __atomic_fetch_or(reinterpret_cast<uint64_t*>(child_mask), constant, __ATOMIC_SEQ_CST);
            *(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) = DONE;
        } // Else, allocation already in progress
        if (old_lock_val == DONE) {
            __atomic_exchange_n(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset, DONE, __ATOMIC_SEQ_CST);
        }
    } // Else, spin until ready
    while(*(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) != DONE);

    // Return old_num_elements to be used in next iteration
    return internal_node_index.index;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
void VDB<level0, level1, level2>::insert_value_into_leaf_node(const Coord xyz, uint64_t index, double value)
{
    uint64_t leaf_pointer_index = index;
    InternalData* leaf_node_array_pointer = vdb_storage_ + leaf_pointer_index;
    uint64_t leaf_offset = calculate_leaf_offset(xyz);
    InternalData& leaf_node_index =
        *(leaf_node_array_pointer + leaf_offset);
    // Isn't currently used, tracks active states
    InternalData* value_mask = leaf_node_array_pointer +
                       // Length of node array
                       LEVEL0_sSIZE;
    value_mask += leaf_offset / INTERNAL_DATA_SIZE;
    uint64_t value_mask_offset = leaf_offset % INTERNAL_DATA_SIZE;
    // Multiple threads can update this, the last update wins
    leaf_node_index.tile_or_value = value;
    uint64_t constant = 1LLU << value_mask_offset;
    __atomic_fetch_or(reinterpret_cast<uint64_t*>(value_mask), constant, __ATOMIC_SEQ_CST);
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
uint64_t VDB<level0, level1, level2>::insert_internal_node_into_hashmap(const Coord xyz)
{
    std::array<uint64_t, 3> rootkey = xyz.get_rootkey(LEVEL2_sLOG2);
    uint64_t hash = Coord::hash(rootkey, HASHMAP_LOG_SIZE);
    uint64_t compressed_xyz = xyz.get_compressed_coord(LEVEL2_sLOG2);
    // Compress coord with num_elements
    uint64_t size_of_axis = 20 - LEVEL2_sLOG2; // Likely 8
    uint64_t amount_to_pad = 64 - (size_of_axis * 3); // Likely 40
    for (uint64_t i = HASHMAP_START; i < HASHMAP_START + HASHMAP_SIZE; ++i) {
        uint64_t index = (hash + i * i) % HASHMAP_SIZE;
        uint64_t expected(std::numeric_limits<uint64_t>::max());
        uint64_t old_val(std::numeric_limits<uint64_t>::max());
        // Use gcc builtin CAS. We must cast the union for the CAS to work.
        // Atomically update storage. We place a temporary value as a lock.
        __atomic_compare_exchange_n(reinterpret_cast<uint64_t*>(vdb_storage_+index), &old_val, compressed_xyz,
                                                false, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST);
        // Update index
        if (old_val == expected) {
            // Atomically update vdb_storage index
            uint64_t old_num_elements = __atomic_fetch_add(&num_elements_, LEVEL2_TOTAL_SIZE, __ATOMIC_SEQ_CST);
            // Atomically exchange storage value. compressed_xyz_index contains old_num_elements
            uint64_t compressed_xyz_index = ((compressed_xyz >> amount_to_pad) << amount_to_pad) + old_num_elements;
            __atomic_exchange_n(reinterpret_cast<uint64_t*>(vdb_storage_+index), compressed_xyz_index, __ATOMIC_SEQ_CST);
            return old_num_elements;
        }
        // Spin if we should be storing a value but another thread is.
        if (((old_val >> amount_to_pad) << amount_to_pad) ==
            ((compressed_xyz >> amount_to_pad) << amount_to_pad)) {
            uint64_t mask = (1LLU << amount_to_pad) - 1;
            while ((old_val & mask) == mask) {
                // TODO: Add thread fence
                old_val = vdb_storage_[index].index;
            }
            return (old_val & mask);
        } // Else, index is occupied and we should continue
    }
    // Error, return 0 and end program
    return 0;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
uint64_t VDB<level0, level1, level2>::calculate_internal_offset(const Coord xyz, InternalNodeLevel inl)
{
    uint64_t sLog2;
    uint64_t log2;
    uint64_t childSLog2;
    switch (inl) {
        case level2_node:
            sLog2 = LEVEL2_sLOG2;
            log2 = LEVEL2_LOG2;
            childSLog2 = LEVEL1_sLOG2;
            break;
        case level1_node:
            sLog2 = LEVEL1_sLOG2;
            log2 = LEVEL1_LOG2;
            childSLog2 = LEVEL0_sLOG2;
            break;
    }
    uint64_t internalOffset =
          (((xyz.x_ & ((1LLU << sLog2) - 1)) >> childSLog2) << (log2 + log2)) +
          (((xyz.y_ & ((1LLU << sLog2) - 1)) >> childSLog2) << log2) +
           ((xyz.z_ & ((1LLU << sLog2) - 1)) >> childSLog2);
    return internalOffset;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
uint64_t VDB<level0, level1, level2>::calculate_leaf_offset(const Coord xyz) {
    uint64_t leafOffset =
          ((xyz.x_ & ((1LLU << LEVEL0_sLOG2) - 1))
                << (LEVEL0_LOG2 + LEVEL0_LOG2)) +
          ((xyz.y_ & ((1LLU << LEVEL0_sLOG2) - 1)) << LEVEL0_LOG2) +
           (xyz.z_ & ((1LLU << LEVEL0_sLOG2) - 1));
    return leafOffset;
}

#endif /* end of include guard: VDB_PRIVATE_HPP_INCLUDED */
