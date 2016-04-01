/**
 * \file vdb.cpp
 * \brief
 *
 * \author Hayden Blauzvern
 */

#ifndef VDB_PRIVATE_HPP_INCLUDED
#define VDB_PRIVATE_HPP_INCLUDED

#include <limits>
#include <array>
#include <iostream>

#include "coord.hpp"
#include "utils.hpp"
#include "vdb.hpp"

// TODO: Clang-format code

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
    memset(vdb_storage_, 0, total_storage_size_);
    // for (uint64_t i = HASHMAP_START + HASHMAP_SIZE; i < total_storage_size_; ++i) {
    //     vdb_storage_[i].index = 0;
    // }
}

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
    value_mask = internal_node_array_pointer +
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
        bool result = insert_node_into_hashmap(xyz, num_elements_);
        if (!result) {
            // TODO: Error inserting, do something
        }
        // Update child
        insert_internal_node(num_elements_, level2_node);
        index.index = num_elements_;
        // Update vdb_storage index
        num_elements_ += LEVEL2_TOTAL_SIZE;
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
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData& value_mask_chunk_2 = *value_mask;
    InternalData& child_mask_chunk_2 = *child_mask;
    // TODO: Change this to not be a reference

    uint64_t mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(child_mask_chunk_2.index, mask_offset)) {
        // No internal node child, add an internal node
        // Update child
        insert_internal_node(num_elements_, level1_node);
        // Update parent's node index, value mask, and child mask
        internal_node_index_2.index = num_elements_;
        update_bit(value_mask_chunk_2.index, mask_offset);
        update_bit(child_mask_chunk_2.index, mask_offset);
        // Update vdb_storage index
        num_elements_ += LEVEL1_TOTAL_SIZE;
    }

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
    value_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    child_mask += internal_offset / INTERNAL_DATA_SIZE; // vdb_storage index
    InternalData& value_mask_chunk_1 = *value_mask;
    InternalData& child_mask_chunk_1 = *child_mask;
    mask_offset = internal_offset % INTERNAL_DATA_SIZE;
    if (!extract_bit(child_mask_chunk_1.index, mask_offset)) {
        // No internal node child, add an internal node
        // Update child
        insert_leaf_node(num_elements_);
        // Update parent's node index, value mask, and child mask
        internal_node_index_1.index = num_elements_;
        update_bit(value_mask_chunk_1.index, mask_offset);
        update_bit(child_mask_chunk_1.index, mask_offset);
        // Update vdb_storage index
        num_elements_ += LEVEL0_TOTAL_SIZE;
    }

    // Extract index to leaf node
    uint64_t leaf_pointer_index = internal_node_index_1.index;
    InternalData* leaf_node_array_pointer = vdb_storage_ + leaf_pointer_index;
    uint64_t leaf_offset = calculate_leaf_offset(xyz);
    InternalData& leaf_node_index =
        *(leaf_node_array_pointer + leaf_offset);
    // Isn't currently used, tracks active states
    value_mask = internal_node_array_pointer +
                       // Length of node array
                       LEVEL0_sSIZE;
    value_mask += leaf_offset / INTERNAL_DATA_SIZE;
    InternalData& value_mask_chunk_0 = *value_mask;
    uint64_t value_mask_offset = leaf_offset % INTERNAL_DATA_SIZE;
    leaf_node_index.tile_or_value = value;
    update_bit(value_mask_chunk_0.index, value_mask_offset);
    return true;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
InternalData VDB<level0, level1, level2>::get_node_from_hashmap(const Coord xyz)
{
    std::array<int, 3> rootkey = xyz.get_rootkey(LEVEL2_sLOG2);
    for (uint64_t i = HASHMAP_START;
            i < HASHMAP_START + HASHMAP_SIZE; ++i) {
        uint64_t index = (Coord::hash(rootkey, HASHMAP_LOG_SIZE) + i * i)
                        % HASHMAP_SIZE;
        if (vdb_storage_[index].index != std::numeric_limits<uint64_t>::max()) {
            return vdb_storage_[index];
        }
    }
    InternalData ret;
    ret.tile_or_value = BACKGROUND_VALUE;
    return ret;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
bool VDB<level0, level1, level2>::insert_node_into_hashmap(const Coord xyz, uint64_t node_index)
{
    std::array<int, 3> rootkey = xyz.get_rootkey(LEVEL2_sLOG2);
    for (uint64_t i = HASHMAP_START;
            i < HASHMAP_START + HASHMAP_SIZE; ++i) {
        uint64_t index = (Coord::hash(rootkey, HASHMAP_LOG_SIZE) + i * i)
                        % HASHMAP_SIZE;
        if (vdb_storage_[index].index == std::numeric_limits<uint64_t>::max()) {
            vdb_storage_[index].index = node_index;
            return true;
        }
    }
    return false;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
void VDB<level0, level1, level2>::insert_internal_node(uint64_t index, InternalNodeLevel inl)
{
    uint64_t size;
    uint64_t mask_size;
    switch (inl) {
        case level2_node:
            size = LEVEL2_sSIZE;
            mask_size = LEVEL2_MASK_SIZE;
        case level1_node:
            size = LEVEL1_sSIZE;
            mask_size = LEVEL1_MASK_SIZE;
    }
    InternalData* internal_node_array_begin = vdb_storage_ + index;
    InternalData* value_mask = internal_node_array_begin +
                               // Length of node array
                               size;
    // memset(value_mask, 0, sizeof(InternalData));
    InternalData* child_mask = internal_node_array_begin +
                               // Length of node array
                               size +
                               // Length of value mask
                               mask_size;
    // memset(child_mask, 0, sizeof(InternalData));
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
void VDB<level0, level1, level2>::insert_leaf_node(uint64_t index)
{
    InternalData* internal_node_array_begin = vdb_storage_ + index;
    InternalData* value_mask = internal_node_array_begin +
                               // Length of node array
                               LEVEL0_sSIZE;
    // memset(value_mask, 0, sizeof(InternalData));
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
          (((xyz.x_ & (1 << sLog2) - 1) >> childSLog2) << (log2 + log2)) +
          (((xyz.y_ & (1 << sLog2) - 1) >> childSLog2) << log2) +
           ((xyz.z_ & (1 << sLog2) - 1) >> childSLog2);
    return internalOffset;
}

template<uint64_t level0, uint64_t level1, uint64_t level2>
uint64_t VDB<level0, level1, level2>::calculate_leaf_offset(const Coord xyz) {
    uint64_t leafOffset =
          ((xyz.x_ & (1 << LEVEL0_sLOG2) - 1)
                << (LEVEL0_LOG2 + LEVEL0_LOG2)) +
          ((xyz.y_ & (1 << LEVEL0_sLOG2) - 1) << LEVEL0_LOG2) +
           (xyz.z_ & (1 << LEVEL0_sLOG2) - 1);
    return leafOffset;
}

#endif /* end of include guard: VDB_PRIVATE_HPP_INCLUDED */

// TODO: Tile insert
// TODO: Refactor code with functions for repeated code in access/insert
