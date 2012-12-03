/******************************************************************************
 * src/tools/radixsort.h
 *
 * Radix sort with a generic template interface.
 * 
 * Based on Kaerkkaeinen, J., & Rantala, T. (2009). Engineering radix sort for
 * strings. SPIRE'09.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
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

template <typename Iterator>
static inline void radixsort(Iterator begin, Iterator end, size_t depth = 0)
{
    if (end - begin < 128) {
        std::sort(begin, end);
        return;
    }

    const size_t K = begin->radixsort_limit();

    // bucketsize[0] is for finished strings, we sort by (ctype,char) in one round.
    size_t bucketsize[K];
    memset(bucketsize, 0, sizeof(bucketsize));
    // count character occurances
    for (Iterator i = begin; i != end; ++i) {
        ++bucketsize[ i->radixsort_index(depth) ];
    }
    // prefix sum
    ssize_t bucketindex[K];
    bucketindex[0] = bucketsize[0];
    size_t last_bucket_size = bucketsize[0];
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i];
        if (bucketsize[i]) last_bucket_size = bucketsize[i];
    }
    // in-place permutation
    for (Iterator i = begin, j; i < end - last_bucket_size; )
    {
        while ( (j = begin + --bucketindex[ i->radixsort_index(depth) ]) > i )
        {
            std::swap(*i, *j);
        }
        i += bucketsize[ i->radixsort_index(depth) ];
    }
    // recursion into bucket
    Iterator i = begin + bucketsize[0];
    for (size_t j = 1; j < K; i += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort(i, i + bucketsize[j], depth+1);
    }
}

template <typename Iterator>
static inline void radixsort_oracle(Iterator begin, Iterator end, size_t depth = 0)
{
    if (end - begin < 128) {
        std::sort(begin, end);
        return;
    }

    const size_t K = begin->radixsort_limit();

    // bucketsize[0] is for finished strings, we sort by (ctype,char) in one round.
    size_t bucketsize[K];
    memset(bucketsize, 0, sizeof(bucketsize));
    // fill oracle and count character occurances
    typedef typename Iterator::value_type::radixsort_oracle_type oracle_type;
    oracle_type* oracle = (oracle_type*)malloc((end - begin) * sizeof(oracle_type));
    size_t ic = 0, jc;
    for (Iterator i = begin; i != end; ++i, ++ic) {
        ++bucketsize[ oracle[ic] = i->radixsort_index(depth) ];
    }
    // prefix sum
    ssize_t bucketindex[K];
    bucketindex[0] = bucketsize[0];
    size_t last_bucket_size = bucketsize[0];
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i];
        if (bucketsize[i]) last_bucket_size = bucketsize[i];
    }
    // in-place permutation
    ic = 0;
    for (Iterator i = begin, j; i < end - last_bucket_size; )
    {
        while ( (jc = --bucketindex[ oracle[ic] ]) > ic )
        {
            j = begin + jc;
            assert( j > i );

            std::swap(*i, *j);
            std::swap(oracle[ic], oracle[jc]);
        }
        i  += bucketsize[ oracle[ic] ];
        ic += bucketsize[ oracle[ic] ];
    }
    free(oracle);
    // recursion into bucket
    Iterator i = begin + bucketsize[0];
    for (size_t j = 1; j < K; i += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort_oracle(i, i + bucketsize[j], depth+1);
    }
}

template <typename Transform, typename Iterator>
static inline void radixsort_transform(Iterator begin, Iterator end, size_t depth = 0)
{
    if (end - begin < 128 || depth >= Transform::maxdepth) {
        std::sort(begin, end);
        return;
    }

    const size_t K = Transform::limit(depth);

    // bucketsize[0] is for finished strings, we sort by (ctype,char) in one round.
    std::vector<size_t> bucketsize (K, 0);
    // count character occurances
    for (Iterator i = begin; i != end; ++i) {
        ++bucketsize[ Transform::index(*i, depth) ];
    }
    // prefix sum
    size_t* bucketindex = new size_t[K];
    bucketindex[0] = bucketsize[0];
    size_t last_bucket_size = bucketsize[0];
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i];
        if (bucketsize[i]) last_bucket_size = bucketsize[i];
    }
    // in-place permutation
    for (Iterator i = begin, j; i < end - last_bucket_size; )
    {
        while ( (j = begin + --bucketindex[ Transform::index(*i, depth) ]) > i )
        {
            std::swap(*i, *j);
        }
        i += bucketsize[ Transform::index(*i, depth) ];
    }
    delete [] bucketindex;
    // subsort equal bucket
    if (bucketsize[0] > 1)
        Transform::subsort(begin, begin + bucketsize[0]);
    // recursion into bucket
    Iterator i = begin + bucketsize[0];
    for (size_t j = 1; j < K; i += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort_transform<Transform,Iterator>(i, i + bucketsize[j], depth+1);
    }
}

template <typename Transform, typename Iterator>
size_t radixsort_transform_maxmemory(size_t n)
{
    // Calculate the maximum memory usage of radixsort_transform_oracle

    const size_t K = Transform::limit(0);
    typedef typename Transform::oracle_type oracle_type;

    return (Transform::maxdepth+1) * K * sizeof(size_t)  // bucketsize for each depth
        + n * sizeof(oracle_type);
}

template <typename Transform, typename Iterator>
static inline void radixsort_transform_oracle(Iterator begin, Iterator end, Transform transform, size_t depth = 0)
{
    if (end - begin < 32 || depth >= Transform::maxdepth) {
        std::sort(begin, end, transform);
        return;
    }

    const size_t K = Transform::limit(depth);

    // bucketsize[0] is for finished strings, we sort by (ctype,char) in one round.
    std::vector<size_t> bucketsize (K, 0);
    // fill oracle and count character occurances
    typedef typename Transform::oracle_type oracle_type;
    oracle_type* oracle = (oracle_type*)malloc((end - begin) * sizeof(oracle_type));
    size_t ic = 0, jc;
    for (Iterator i = begin; i != end; ++i, ++ic) {
        ++bucketsize[ oracle[ic] = transform.index(*i,depth) ];
    }
    // prefix sum
    size_t* bucketindex = (size_t*)malloc(K * sizeof(size_t));
    bucketindex[0] = bucketsize[0];
    size_t last_bucket_size = bucketsize[0];
    for (size_t i = 1; i < K; ++i) {
        bucketindex[i] = bucketindex[i-1] + bucketsize[i];
        if (bucketsize[i]) last_bucket_size = bucketsize[i];
    }
    // in-place permutation
    ic = 0;
    for (Iterator i = begin, j; i < end - last_bucket_size; )
    {
        while ( (jc = --bucketindex[ oracle[ic] ]) > ic )
        {
            j = begin + jc;
            assert( j > i );

            std::swap(*i, *j);
            std::swap(oracle[ic], oracle[jc]);
        }
        i  += bucketsize[ oracle[ic] ];
        ic += bucketsize[ oracle[ic] ];
    }
    free(bucketindex);
    free(oracle);
    // subsort equal bucket
    if (bucketsize[0] > 1)
        transform.subsort(begin, begin + bucketsize[0]);
    // recursion into bucket
    Iterator i = begin + bucketsize[0];
    for (size_t j = 1; j < K; i += bucketsize[j++]) {
        if (bucketsize[j] <= 1) continue;
        radixsort_transform_oracle<Transform,Iterator>(i, i + bucketsize[j], transform, depth+1);
    }
}
