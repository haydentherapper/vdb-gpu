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

// TODO: Make access safe with concurrent access and insertion by adding spin locks
template<uint64_t level0, uint64_t level1, uint64_t level2>
double VDB<level0, level1, level2>::random_access(const Coord xyz)
{
    // Index to internal node or background
    InternalData index = get_node_from_hashmap(xyz);
    if (index.tile_or_value == BACKGROUND_VALUE) {
        return index.tile_or_value;
    }

    // Pointer to first index into internal data array
    InternalData* internal_node_array_pointer = vdb_storage_ + index.index;
    uint64_t internal_offset = calculate_internal_offset(xyz, level2_node);
    InternalData internal_node_index =
        *(internal_node_array_pointer + internal_offset);

    // Check value mask to see if node exists
    InternalData* value_mask = internal_node_array_pointer +
                               // Length of node array
                               LEVEL2_sSIZE;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData value_mask_chunk = *value_mask;
    uint64_t value_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(value_mask_chunk.index, value_mask_offset)) {
        // If value_mask not set, return background
        return BACKGROUND_VALUE;
    }

    InternalData* child_mask = internal_node_array_pointer +
                               // Length of node array
                               LEVEL2_sSIZE +
                               // Length of value mask
                               LEVEL2_MASK_SIZE;
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData child_mask_chunk = *child_mask;

    uint64_t child_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(child_mask_chunk.index, child_mask_offset)) {
        // Return tile value
        return internal_node_index.tile_or_value;
    }

    // Recurse on level 1 internal node
    index.index = internal_node_index.index;
    internal_node_array_pointer = vdb_storage_ + index.index;
    internal_offset = calculate_internal_offset(xyz, level1_node);
    internal_node_index = *(internal_node_array_pointer + internal_offset);
    value_mask = internal_node_array_pointer +
                 // Length of node array
                 LEVEL1_sSIZE;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    value_mask_chunk = *value_mask;
    value_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(value_mask_chunk.index, value_mask_offset)) {
        // If value_mask not set, return background
        return BACKGROUND_VALUE;
    }
    child_mask = internal_node_array_pointer +
               // Length of node array
               LEVEL1_sSIZE +
               // Length of value mask
               LEVEL1_MASK_SIZE;
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask_chunk = *child_mask;
    child_mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(child_mask_chunk.index, child_mask_offset)) {
        // Return tile value
        return internal_node_index.tile_or_value;
    }

    // Extract index to leaf node
    uint64_t leaf_pointer_index = internal_node_index.index;
    InternalData* leaf_node_array_pointer = vdb_storage_ + leaf_pointer_index;
    uint64_t leaf_offset = calculate_leaf_offset(xyz);
    InternalData leaf_node_index =
        *(leaf_node_array_pointer + leaf_offset);
    // Isn't currently used, tracks active states
    value_mask = leaf_node_array_pointer +
                 // Length of node array
                 LEVEL0_sSIZE;
    value_mask += leaf_offset / INTERNAL_DATA_SIZE;
    value_mask_chunk = *value_mask;
    value_mask_offset = leaf_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(value_mask_chunk.index, value_mask_offset)) {
        return BACKGROUND_VALUE;
    }
    // Return voxel value
    return leaf_node_index.tile_or_value;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
bool VDB<level0, level1, level2>::random_insert(const Coord xyz, double value)
{
    // First, get index to internal node or background
    InternalData index = get_node_from_hashmap(xyz);
    if (index.tile_or_value == BACKGROUND_VALUE) {
        // Node does not exist, insert into hashmap to update parent
        uint64_t result = insert_node_into_hashmap(xyz);
        if (result == 0) {
            std::cout << "Not enough room in hashtable, exiting" << std::endl;
            exit(EXIT_FAILURE);
        }
        index.index = result;
    }

    // Pointer to first index into internal data array
    InternalData* internal_node_array_pointer = vdb_storage_ + index.index;
    uint64_t internal_offset = calculate_internal_offset(xyz, level2_node);
    InternalData& internal_node_index_2 =
        *(internal_node_array_pointer + internal_offset);
    InternalData* value_mask = internal_node_array_pointer +
                               // Length of node array
                               LEVEL2_sSIZE;
    InternalData* child_mask = internal_node_array_pointer +
                               // Length of node array
                               LEVEL2_sSIZE +
                               // Length of value mask
                               LEVEL2_MASK_SIZE;
   InternalData* lock_flag_array = internal_node_array_pointer +
                              // Length of node array
                              LEVEL2_sSIZE +
                              // Length of value mask
                              LEVEL2_MASK_SIZE +
                              // Length of child mask
                              LEVEL2_MASK_SIZE;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    lock_flag_array += internal_offset / NUM_LOCK_FLAGS;

    InternalData& child_mask_chunk_2 = *child_mask;
    uint64_t mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    uint64_t lock_offset = internal_offset % NUM_LOCK_FLAGS;
    if (!extract_bit(child_mask_chunk_2.index, mask_offset)) {
        // No internal node child, add an internal node
        uint8_t old_lock_val = __sync_lock_test_and_set(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset, IP);
        if (old_lock_val == READY) {
            uint64_t old_num_elements = __sync_fetch_and_add(&num_elements_, LEVEL1_TOTAL_SIZE);
            internal_node_index_2.index = old_num_elements;
            uint64_t constant = 1LLU << mask_offset;
            __sync_fetch_and_or(reinterpret_cast<uint64_t*>(value_mask), constant);
            __sync_fetch_and_or(reinterpret_cast<uint64_t*>(child_mask), constant);
            *(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) = READY;
        } // Else, allocation already in progress
    } // Else, spin until ready
    while(*(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) == IP);

    // Recurse on level 1 internal node
    index.index = internal_node_index_2.index;
    internal_node_array_pointer = vdb_storage_ + index.index;
    internal_offset = calculate_internal_offset(xyz, level1_node);
    InternalData& internal_node_index_1 = *(internal_node_array_pointer + internal_offset);
    value_mask = internal_node_array_pointer +
                               // Length of node array
                               LEVEL1_sSIZE;
    child_mask = internal_node_array_pointer +
               // Length of node array
               LEVEL1_sSIZE +
               // Length of value mask
               LEVEL1_MASK_SIZE;
   lock_flag_array = internal_node_array_pointer +
              // Length of node array
              LEVEL1_sSIZE +
              // Length of value mask
              LEVEL1_MASK_SIZE +
              // Length of child mask
              LEVEL1_MASK_SIZE;
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    lock_flag_array += internal_offset / NUM_LOCK_FLAGS;

    InternalData& child_mask_chunk_1 = *child_mask;
    mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    lock_offset = internal_offset % NUM_LOCK_FLAGS;
    if (!extract_bit(child_mask_chunk_1.index, mask_offset)) {
        // No leaf node child, add a leaf node
        uint8_t old_lock_val = __sync_lock_test_and_set(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset, IP);
        if (old_lock_val == READY) {
            uint64_t old_num_elements = __sync_fetch_and_add(&num_elements_, LEVEL0_TOTAL_SIZE);
            internal_node_index_1.index = old_num_elements;
            uint64_t constant = 1LLU << mask_offset;
            __sync_fetch_and_or(reinterpret_cast<uint64_t*>(value_mask), constant);
            __sync_fetch_and_or(reinterpret_cast<uint64_t*>(child_mask), constant);
            *(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) = READY;
        } // Else, allocation already in progress
    } // Else, spin until ready
    while(*(reinterpret_cast<uint8_t*>(lock_flag_array) + lock_offset) == IP);

    // Extract index to leaf node
    uint64_t leaf_pointer_index = internal_node_index_1.index;
    InternalData* leaf_node_array_pointer = vdb_storage_ + leaf_pointer_index;
    uint64_t leaf_offset = calculate_leaf_offset(xyz);
    InternalData& leaf_node_index =
        *(leaf_node_array_pointer + leaf_offset);
    // Isn't currently used, tracks active states
    value_mask = leaf_node_array_pointer +
                       // Length of node array
                       LEVEL0_sSIZE;
    value_mask += leaf_offset / INTERNAL_DATA_SIZE;
    uint64_t value_mask_offset = leaf_offset % INTERNAL_DATA_SIZE;
    // Multiple threads can update this, the last update wins
    leaf_node_index.tile_or_value = value;
    uint64_t constant = 1LLU << value_mask_offset;
    __sync_fetch_and_or(reinterpret_cast<uint64_t*>(value_mask), constant);
    return true;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
InternalData VDB<level0, level1, level2>::get_node_from_hashmap(const Coord xyz)
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
uint64_t VDB<level0, level1, level2>::insert_node_into_hashmap(const Coord xyz)
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
        // Use gcc builtin CAS. We must cast the union for the CAS to work.
        // Atomically update storage. We place a temporary value as a lock.
        uint64_t old_val = __sync_val_compare_and_swap(reinterpret_cast<uint64_t*>(vdb_storage_+index), expected, compressed_xyz);
        // Update index
        if (old_val == expected) {
            // Atomically update vdb_storage index
            uint64_t old_num_elements = __sync_fetch_and_add(&num_elements_, LEVEL2_TOTAL_SIZE);
            // TODO: Since other threads update num_elements_, passed num_elems is stale
            // Atomically exchange storage value. compressed_xyz_index contains old_num_elements
            uint64_t compressed_xyz_index = ((compressed_xyz >> amount_to_pad) << amount_to_pad) + old_num_elements;
            __sync_lock_test_and_set(reinterpret_cast<uint64_t*>(vdb_storage_+index), compressed_xyz_index);
            return old_num_elements;
        }
        // Spin if we should be storing a value but another thread is.
        if (((old_val >> amount_to_pad) << amount_to_pad) ==
            ((compressed_xyz >> amount_to_pad) << amount_to_pad)) {
            uint64_t mask = (1LLU << amount_to_pad) - 1;
            while ((old_val & mask) == mask) {
                // TODO: thread fence
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
        case level2_node: sLog2 = LEVEL2_sLOG2;
            log2 = LEVEL2_LOG2;
            childSLog2 = LEVEL1_sLOG2;
        case level1_node: sLog2 = LEVEL1_sLOG2;
            log2 = LEVEL1_LOG2;
            childSLog2 = LEVEL0_sLOG2;
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
