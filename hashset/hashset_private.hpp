/**
 * \file hashset-private.hpp
 * \brief Template declarations for HashSet.
 *
 * \author Hayden Blauzvern and Casey Chu
 */


#ifndef HASHSET_PRIVATE_HPP_INCLUDED
#define HASHSET_PRIVATE_HPP_INCLUDED

// Default constructor with templated size
template<typename T, size_t N>
HashSet<T,N>::HashSet()
    : size_(N), collisions_(0)
{
    table_.fill(nullptr);
}

// Returns the number of items stored in the set.
template<typename T, size_t N>
size_t HashSet<T,N>::size() const
{
    return size_;
}

// Adds x to the table.
template<typename T, size_t N>
void HashSet<T,N>::insert(const T& x)
{
    size_t size = size();
    for (size_t i = 0; i < size; ++i) {
        size_t index = (hash(x) + i * i) % size;
        if (table_[index] == nullptr) {
            table_[index] = x;
            return;
        }
        ++collisions_;
    }
}

// Returns true if x is present.
template<typename T, size_t N>
bool HashSet<T,N>::exists(const T& x) const
{
    size_t size = size();
    for (size_t i = 0; i < size; ++i) {
        size_t index = (hash(x) + i * i) % size;
        if (table_[index] != nullptr &&
            table_[index] == x) {
            return true;
        }
    }
    return false;
}

// Returns the number of times an insert has found a non-empty bucket.
template<typename T, size_t N>
size_t HashSet<T,N>::collisions() const
{
    return collisions_;
}

#endif
