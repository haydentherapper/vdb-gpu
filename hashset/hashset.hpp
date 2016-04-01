/**
 * \file hashset.hpp
 * \brief HashSet implements a set using separate chaining, thereby
 *        guaranteeing O(1) for insert, delete, and existence.
 * \author Hayden Blauzvern and Casey Chu
 */


#ifndef HASHSET_HPP_INCLUDED
#define HASHSET_HPP_INCLUDED

#include <array>

/// Implements a set using separate chaining, thereby
/// guaranteeing O(1) for insert, delete, and existence.
template<typename T, size_t N>
class HashSet {
public:
    /// Default constructor on templated size
    HashSet();

    /// Destructor.
    ~HashSet() = default;

    /// Returns the number of items stored in the set.
    size_t size() const;

    /// Adds x to the set.
    void insert(const T& x);

    /// Returns true if x is present in the set.
    bool exists(const T& x) const;

    /// Returns the number of times an insert has found a non-empty bucket.
    size_t collisions() const;

private:
    /// Disabled copy constructor.
    HashSet(HashSet<T>& other) = delete;

    /// Disabled assignment operator.
    HashSet<T>& operator=(HashSet<T>& rhs) = delete;

    /// The type of each bucket.
    typedef std::array<T, N> table_t;

    /// The table we'll put the values in.
    table_t table_;

    /// The number of elements in the hash table.
    size_t size_;

    /// The number of times an insert has found a non-empty bucket.
    size_t collisions_;
};

#include "hashset-private.hpp"
#endif
