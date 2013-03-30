/******************************************************************************
 * src/tools/lp_hash_table.h
 *
 * A size bounded hash table with linear probing.
 *
 * The complete structure uses only a limited amount of memory, thus inserted
 * key/value pairs may be discarded when the table gets too full. If it is
 * discarded, then for future gets the hash table will indicate that an pair
 * may have been discarded. This indication can have false positives.
 *
 ******************************************************************************
 * Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef LP_HASH_TABLE_H_
#define LP_HASH_TABLE_H_

#include <vector>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

namespace lp_hash_table {

template <typename KeyType>
inline size_t hash(const KeyType& key, unsigned int i);

template <typename KeyType>
inline size_t hash(const KeyType& key);

template <>
inline size_t hash(const uint32_t& key, unsigned int i)
{
    static const uint32_t factor[] = {
        2654435769u, // -phi
        3474701532u, // hat{phi}
        3373259426u, // pi
        2918732888u, // e
        3037000499u, // sqrt(2)
        3719550786u, // sqrt(3)
    };

    return (key * factor[i]);
}

template <>
inline size_t hash(const uint32_t& key)
{
    return hash(key,0);
}

template <>
inline size_t hash(const uint64_t& key, unsigned int i)
{
    static const uint64_t factor[] = {
        11400714819323198485llu, // -phi
        14923729446516375050llu, // hat{phi}
        14488038916154245684llu, // pi
        12535862302449814170llu, // e
        13043817825332782212llu, // sqrt(2)
        15975348984942515101llu, // sqrt(3)
    };

    return (key * factor[i]);
}

template <>
inline size_t hash(const uint64_t& key)
{
    return hash(key,0);
}

/**
 * Class for a simple Bloom fiter.
 */
template <typename KeyType>
class BloomFilter
{
public:
    /// the key type of the Bloom filter
    typedef KeyType     key_type;

public:

    // approximate false positive rate:
    // (1 - exp(- kn/m))^k

    // the optimal k for m/n is
    // k = ln(2) * m/n

    // with this k, bit array size m depending on probability p:
    // m  = - 2 n ln p

    /// bit array size m
    size_t      m_M;

    /// number of hash functions k
    int         m_K;

protected:

    /// the bit array, represented as size_ts
    std::vector<size_t> m_bits;

    /// set the b-th bit in the array
    inline void set_bit(unsigned int b)
    {
        unsigned int off = b / (8 * sizeof(size_t));
        unsigned int sht = b % (8 * sizeof(size_t));
        
        m_bits[off] |= (((size_t)1) << sht);
    }
    
    /// return the b-th bit in the array (for != 0 comparisons)
    inline size_t get_bit(unsigned int b) const
    {
        unsigned int off = b / (8 * sizeof(size_t));
        unsigned int sht = b % (8 * sizeof(size_t));
        
        return (m_bits[off] & (((size_t)1) << sht));
    }

public:

    /// Create Bloom filter with M bits and K hash functions, its memory usage
    /// is M/8 + O(1) bytes.
    explicit BloomFilter(size_t M, int K)
        : m_M(M),
          m_K(K),
          m_bits( ((M/8) + sizeof(size_t)-1) / sizeof(size_t), 0 )
    {
    }

    /// Clear Bloom filter
    inline void clear()
    {
        std::fill(m_bits.begin(), m_bits.end(), 0);
    }

    /// Mark the given key in the Bloom filter
    inline void set(const key_type& key)
    {
        for (int k = 0; k < m_K; ++k)
        {
            set_bit( hash(key,k) % m_M );
        }
    }

    /// Retreive whether the key was marked in the filter - with the probabilty
    /// of false positives!
    inline bool get(const key_type& key) const
    {
        for (int k = 0; k < m_K; ++k)
        {
            if (!get_bit( hash(key,k) % m_M )) return false;
        }
        return true;
    }
};

/// Enumeration indicating result of hash table lookup
enum result_type {
    FOUND = 1, NOT_FOUND = 0, UNKNOWN = -1
};

/**
 * A size bounded hash table with linear probing.
 *
 * The complete structure uses only a limited amount of memory, thus inserted
 * key/value pairs may be discarded when the table gets too full. If it is
 * discarded, then for future gets the hash table will indicate that an pair
 * may have been discarded. This indication can have false positives.
 *
 * By using timestamps, the clear() operations is very fast.
 */
template <typename KeyType, typename ValueType>
class LPHashTable
{
public:
    /// The key type
    typedef KeyType     key_type;

    /// The value type
    typedef ValueType   value_type;

protected:
    
    /// The type of all timestamps in the hash table
    typedef uint32_t    ts_type;

    /// Structure of the cells in the hash table.    
    struct Cell
    {
        key_type        key;
        value_type      value;
        ts_type         ts;
    };

    /// simple Bloom filter to indicate discarded keys (with false positives)
    BloomFilter<key_type>   m_bloom;

    /// The array of hash table cells
    std::vector<Cell>   m_cell;

    /// All cells with timestamp <= are considered invalid (set after last clear())
    ts_type             m_invalid_ts;

    /// Current timestamp of previous insertion
    ts_type             m_curr_ts;

    /// Length of the linear probing
    unsigned int        m_probe_length;

public:

    /// Count number of insertions
    size_t              m_count_insert;

    /// Count number of evictions
    size_t              m_count_eviction;

    /// Count number of successful lookups
    size_t              m_count_found;

    /// Count number of unsuccessful lookups
    size_t              m_count_notfound;

    /// Count number of unknown lookups
    size_t              m_count_unknown;

public:

    /// Construct a hash table with associated Bloom filter of maximum size M.
    /// Allocate 1/8 of M for the bloom filter, and the rest to the hash table cells.
    explicit LPHashTable(size_t M)
        : m_bloom( (M * 1/8) * 8, 1 ),                   // 1/8 of mem for bloom filter
          m_cell(  (M * 7/8) / sizeof(Cell) ),           // 7/8 of mem for cells
          m_invalid_ts(0), m_curr_ts(0),
          m_probe_length( sqrt(m_cell.size()) ),        // sqrt(n) probes
          m_count_insert(0), m_count_found(0),
          m_count_notfound(0), m_count_unknown(0)
    {
        // clear timestamps
        for (size_t i = 0; i < m_cell.size(); ++i)
            m_cell[i].ts = m_invalid_ts;
    }

    /// Clear the data structure by raising the minimum valid timestamp
    void clear()
    {
        m_invalid_ts = m_curr_ts;
        m_bloom.clear();
    }

    /// Insert a key->value pair in the hash table.
    void set(const key_type& key, const value_type& value)
    {
        ++m_count_insert;

        size_t idx = hash(key) % m_cell.size();

        if (m_curr_ts == std::numeric_limits<ts_type>::max())
        {
            // timestamp wrapped: set all to one and reset timestamp counters
            for (size_t i = 0; i < m_cell.size(); ++i)
                m_cell[i].ts = 1; // first valid timestamp
            m_invalid_ts = 0;
            m_curr_ts = 1;      // all timestamps 1 are "one insert"
        }

        ++m_curr_ts;    // increment timestamp for this insert

        ts_type chain_min_ts = m_curr_ts;
        size_t chain_min_idx = -1;

        unsigned int probes = 0;

        // look for identical key or a free place.
        while(
            // found a free cell in linear chain
            (m_cell[idx].ts > m_invalid_ts) &&
            // found key in linear chain
            (m_cell[idx].key != key)
            )
        {
            // sufficiently many chain cells probed -> displace min_idx
            if (++probes > m_probe_length)
            {
                assert(chain_min_idx != (size_t)-1);
                idx = chain_min_idx;

                assert(0);
                ++m_count_eviction;
                m_bloom.set( m_cell[idx].key );
                assert( m_bloom.get(m_cell[idx].key) );

                break;
            }

            // save minimum ts in chain
            if (m_cell[idx].ts < chain_min_ts) {
                chain_min_ts = m_cell[idx].ts;
                chain_min_idx = idx;
            }

            // continue to next entry in chain
            ++idx;
            idx %= m_cell.size();
        }

        // fill, update or displace idx
        m_cell[idx].key = key;
        m_cell[idx].value = value;
        m_cell[idx].ts = m_curr_ts;
    }

    /// Lookup a key, if found return 1 and copy the value. If never added,
    /// return -1. If unknown, return 0.
    result_type get(const key_type& key, value_type& value)
    {
        size_t idx = hash(key) % m_cell.size();
        unsigned int probes = 0;

        while(
            // found a free m_cell in linear chain
            (m_cell[idx].ts > m_invalid_ts) &&
            // found key in linear chain
            (m_cell[idx].key != key)
            )
        {
            // sufficiently many chain m_cells probed -> not found
            if (++probes > m_probe_length) 
                break;

            // continue to next entry in chain
            ++idx;
            idx %= m_cell.size();
        }

        if (m_cell[idx].key == key &&
            m_cell[idx].ts > m_invalid_ts)
        {
            value = m_cell[idx].value;
            ++m_count_found;
            return FOUND;
        }
        else
        {
            // bloom filter marks those keys, for which an entry was dropped. 
            if (m_bloom.get(key)) {
                ++m_count_unknown;
                return UNKNOWN;
            }
            else {
                ++m_count_notfound;
                return NOT_FOUND;
            }
        }
    }
};

static int testmain()
{
    LPHashTable<uint32_t,uint32_t> ht(2*1024);

    srand(2344234);

    for(int i = 0; i < 1000; ++i)
    {
        int k = i; // rand() % 2000;
        ht.set(k,i);
    }

    srand(2344234);

    for(int i = 0; i < 2000; ++i)
    {
        int k = i; // rand() % 2000;
        uint32_t out;
        result_type r = ht.get(k,out);
        std::cout << "i = " << i << " - " << r << "\n";
        if (i < 1000)
            assert( r == FOUND || r == UNKNOWN );
        else
            assert( r == NOT_FOUND || r == UNKNOWN );
    }

    return 0;
}

} // namespace lp_hash_table

#endif // LP_HASH_TABLE_H_
