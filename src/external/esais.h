// -*- mode: c++; mode: rebox; fill-column: 110; rebox-min-fill-column: 110 -*-

/******************************************************************************
 * src/external/esais.h
 *
 * eSAIS main algorithm
 *
 * Based on Bingmann, T., Fischer, J., Osipov, V. (2013).
 * Inducing suffix and LCP arrays in external memory. ALENEX'13.
 *
 ******************************************************************************
 * Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
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

namespace esais
{
/*************************************************************************************************************
 * Compilation Options                                                                                       *
 *************************************************************************************************************/

// run self-verifying parts of the code, very slow!
#define ESAIS_SELF_CHECK        0

// enable LCP calculation
#define ESAIS_LCP_CALC          LCP_CALC

// enable specific LCP calculation method
#if !defined(ESAIS_LCP_CALC_EXT) && !defined(ESAIS_LCP_CALC_INT)
#define ESAIS_LCP_CALC_EXT      0
#define ESAIS_LCP_CALC_INT      LCP_CALC
#endif

// enable discarding of sequences of unique characters
#ifndef ESAIS_DISCARD_UNIQUES
#define ESAIS_DISCARD_UNIQUES   0
#endif

// log PQ fill levels to file
#define ESAIS_LOG_PQFILL        0

// log estimate of wasted space in tuples
#ifndef ESAIS_COUNT_WASTED
#define ESAIS_COUNT_WASTED      0
#endif

// *******************************************************************************************************
// *** Debugging Switches

static const bool debug = true;

static const bool debug_substring_split = false;
static const bool debug_substring_decoder = false;
static const bool debug_substring_merge = false;
static const bool debug_lexnamepairs = false;
static const bool debug_rerank = false;
static const bool debug_recursive_input = false;
static const bool debug_recursive_output = false;
static const bool debug_sstarlcp = false;
static const bool debug_sstarlcp_merge = false;
static const bool debug_sstarlcp_result = false;
static const bool debug_sstarlcp_use = false;
static const bool debug_induce_input = false;
static const bool debug_induce_split = false;
static const bool debug_induce_arrays = false;
static const bool debug_induceL = false;
static const bool debug_induceL_lcp = false;
static const bool debug_induceS = false;
static const bool debug_induceS_lcp = false;
static const bool debug_output_write = false;
static const bool debug_output_merge = false;

/*************************************************************************************************************
 * Helper Routines                                                                                           *
 *************************************************************************************************************/

// {{{ Debug macros

#if ESAIS_LOG_PQFILL
#define LOG_SIZE(x)     do { (x); } while(0)
#else
#define LOG_SIZE(x)     do {} while(0)
#endif

#if ESAIS_COUNT_WASTED
#define LOG_WASTED(x)   do { g_wasted_iovolume += (x); } while(0)
#else
#define LOG_WASTED(x)   do {} while(0)
#endif

#if ESAIS_LCP_CALC
#define ESAIS_LCP_CALCX(x)      x
#else
#define ESAIS_LCP_CALCX(x)
#endif

#if ESAIS_LCP_CALC_INT
#define ESAIS_LCP_CALCX_INT(x)      x
#else
#define ESAIS_LCP_CALCX_INT(x)
#endif

#if ESAIS_COUNT_WASTED

/// count wasted io volumne, hacky global variable
static size_t   g_wasted_iovolume = 0;

#endif // ESAIS_COUNT_WASTED

/// count maximum size of internal memory RMQ structure, hacky global variable
static size_t   g_mainmemlcp = 0;

#if ESAIS_LCP_CALC_EXT && !defined(ESAIS_LCP_CALC_INT)
#define ESAIS_LCP_CALC_INT      0
#endif

#if ESAIS_LCP_CALC_INT && !defined(ESAIS_LCP_CALC_EXT)
#define ESAIS_LCP_CALC_EXT      0
#endif

// }}} Debug macros

// {{{ Variable stream writing with varints

//! \brief put some type as bytes at the end of a stream
template <typename Stream, typename Type>
static inline Stream& stream_put(Stream& s, const Type& t)
{
    assert( sizeof(t) % sizeof(typename Stream::block_type::type) == 0 );
    assert( sizeof(t) < Stream::block_type::size );

    for (size_t j = 0; j < sizeof(Type) / sizeof(typename Stream::block_type::type); ++j)
    {
        s << ( ((typename Stream::block_type::type*)&t)[j] );
    }
    return s;
}

//! \brief get some type as bytes from the end of a stream
template <typename Stream, typename Type>
static inline Stream& stream_get(Stream& s, Type& t)
{
    assert( sizeof(t) % sizeof(typename Stream::block_type::type) == 0 );
    assert( sizeof(t) < Stream::block_type::size );

    for (size_t j = 0; j < sizeof(Type) / sizeof(typename Stream::block_type::type); ++j)
    {
        s >> ( ((typename Stream::block_type::type*)&t)[j] );
    }
    return s;
}

//! \brief return number of bytes needed to encode a varint of given value
static unsigned int sizeof_varint(unsigned int x)
{
    if (x < 128) return 1;
    if (x < 128*128) return 2;
    if (x < 128*128*128) return 3;
    if (x < 128*128*128*128) return 4;
    abort();
    return 0;
}

//! \brief encode a varint and append it to a stream
template <typename Stream>
static inline Stream& stream_put_varint(Stream& s, size_t v)
{
    if (v < 128) {
        s << ((v >> 0) & 0x7F);
    }
    else if (v < 128*128) {
        s << (((v >> 0) & 0x7F) | 0x80);
        s << ((v >> 7) & 0x7F);
    }
    else if (v < 128*128*128) {
        s << (((v >> 0) & 0x7F) | 0x80);
        s << (((v >> 7) & 0x7F) | 0x80);
        s << ((v >> 14) & 0x7F);
    }
    else if (v < 128*128*128*128) {
        s << (((v >> 0) & 0x7F) | 0x80);
        s << (((v >> 7) & 0x7F) | 0x80);
        s << (((v >> 14) & 0x7F) | 0x80);
        s << ((v >> 21) & 0x7F);
    }
    else {
        abort(); // enlarge allowed range if ever needed.
    }
    return s;
}

//! \brief decode a varint from end of stream
template <typename Stream>
static inline Stream& stream_get_varint(Stream& s, unsigned int& v)
{
    unsigned char in;
    s >> in;

    if (in & 0x80) {
        v = (in & 0x7F);

        s >> in;
        if (in & 0x80) {
            v |= ((unsigned int)in & 0x7F) << 7;

            s >> in;
            if (in & 0x80) {
                v |= ((unsigned int)in & 0x7F) << 14;

                s >> in;
                if (in & 0x80) {
                    abort();
                }
                else {
                    v |= (unsigned int)in << 21;
                }
            }
            else {
                v |= (unsigned int)in << 14;
            }
        }
        else {
            v |= (unsigned int)in << 7;
        }
    }
    else {
        v = (in & 0x7F);
    }
    return s;
}

//! \brief decode a varint from end of stream
template <typename Stream>
static inline Stream& stream_get_varint(Stream& s, stxxl::uint40& v)
{
    // expand this proc to decode larger varints if ever needed.
    unsigned int x;
    stream_get_varint(s,x);
    v = x;
    return s;
}

// }}} Variable stream writing with varints

// {{{ STXXL Extensions

/**
 * Creates a stream object that reads from a stxxl::vector from back to front. This class is an adapter for
 * the reverse reader now in the STXXL; it matches the interface of recursive calls.
 */
template <class InputIterator>
class my_vector_bufreader_reverse : public stxxl::vector_bufreader_reverse<InputIterator>
{
public:

    typedef my_vector_bufreader_reverse self_type;

    typedef stxxl::vector_bufreader_reverse<InputIterator> super_type;

    my_vector_bufreader_reverse(InputIterator begin, InputIterator end, stxxl::unsigned_type nbuffers = 0)
        : super_type(begin, end, nbuffers)
    {
    }

    self_type& operator -- ()
    {
        super_type::operator++();
        return *this;
    }

    self_type& rewind(size_t /* mem */)
    {
        super_type::rewind();
        return *this;
    }

    //! deallocate stream reader. noop.
    self_type& finish()
    {
        return *this;
    }
};

/*************************************************************************************************************
 * Sorter for Variable Length Objects                                                                        *
 *************************************************************************************************************/

/**
 * The VarlengthSorter implements non-recursive external memory sorting for variable length objects. The
 * objects are passed by the caller in form of already-sorted runs of sizes theta(M), which are written to
 * disk. No aggregation of objects is done in this class. Once all runs are written, the sequences can be
 * streamed back, merged using a tournament tree and outputted in globally sorted order. Each objects is
 * serialized when written and deserialized when read using it's enclosed .serialize(stream) and
 * .deserialize(stream) methods.
 */
template <typename ValueType, int BlockSize>
class VarlengthSorter
{
public:

    // *** Template Arguments

    /// values serialized.
    typedef ValueType           value_type;

    /// block size of written runs
    static const int block_size = BlockSize;

    // *** Constructed Types

    typedef STXXL_DEFAULT_ALLOC_STRATEGY alloc_strategy_type;
    alloc_strategy_type     m_alloc_strategy;

    typedef stxxl::BID<block_size>   bid_type;

    typedef std::vector<bid_type>   bid_vector_type;

    typedef typename bid_vector_type::iterator   bid_iterator_type;

    /// create typed blocks containing characters
    typedef stxxl::typed_block<block_size, unsigned char> block_type;

    /// create byte-oriented buffered output streams
    typedef stxxl::buf_ostream<block_type, bid_iterator_type>       buf_ostream_type;

    /// create byte-oriented buffered input streams
    typedef stxxl::buf_istream<block_type, bid_iterator_type>       buf_istream_type;

    /// list of input streams, used when merging runs
    typedef std::vector<buf_istream_type*>       buf_istream_vector_type;

protected:

    /** Functor for tournament tree which feeds the input streams into the loser tree as required by its
     * interface. */
    struct StreamCompare
    {
        buf_istream_vector_type&        m_istream_vector;

        std::vector<value_type>         m_head;

        std::vector<bool>               m_exists;

        StreamCompare(buf_istream_vector_type& istream_vector)
            : m_istream_vector(istream_vector)
        {
        }

        void resize()
        {
            m_head.resize( m_istream_vector.size() );

            m_exists.clear();
            m_exists.resize( m_istream_vector.size(), true );
        }

        void fetch(unsigned int seq)
        {
            buf_istream_type& is = *m_istream_vector[seq];

            m_exists[seq] = m_head[seq].deserialize( is );
        }

        inline bool exists(unsigned int seq) const
        {
            return m_exists[seq];
        }

        inline int operator()(unsigned int seqa, unsigned int seqb) const
        {
            assert( exists(seqa) && exists(seqb) );
            return m_head[seqa].cmp( m_head[seqb] );
        }
    };

protected:

    /// block_manager of stxxl from which to allocate runs
    stxxl::block_manager*               m_bm;

    // *** Structures holding references to written runs

    std::vector<bid_vector_type>        m_block_seqs;

    // *** Structures when merging runs

    buf_istream_vector_type     m_istream_vector;

    StreamCompare               m_streamcmp;

    LoserTree3Way<StreamCompare>    m_losertree;

public:

    VarlengthSorter()
        : m_bm(stxxl::block_manager::get_instance()),
          m_streamcmp(m_istream_vector),
          m_losertree(m_streamcmp)
    {
    }

    ~VarlengthSorter()
    {
        clear();
    }

    // write out buffer of sorted objects into one run
    template <typename ConstIterator>
    void write(size_t totalsize, ConstIterator begin, ConstIterator end)
    {
        bid_vector_type newbids ( (totalsize + block_size - 1) / block_size );
        m_bm->new_blocks(m_alloc_strategy, newbids.begin(), newbids.end());

        buf_ostream_type os (newbids.begin(), 2);

        while( begin != end )
        {
            (*begin).serialize(os); ++begin;
        }

        os.fill(0);

        m_block_seqs.push_back(newbids);
    }

    // write out buffer of sorted objects into one run
    template <typename ConstIterator, typename Serializer>
    void write(size_t totalsize, ConstIterator begin, ConstIterator end, Serializer serializer)
    {
        bid_vector_type newbids ( (totalsize + block_size - 1) / block_size );
        m_bm->new_blocks(m_alloc_strategy, newbids.begin(), newbids.end());

        buf_ostream_type os (newbids.begin(), 2);

        while( begin != end )
        {
            serializer(os, *begin); ++begin;
        }

        os.fill(0);

        m_block_seqs.push_back(newbids);
    }

    /// finish input of runs and prepare tournament tree for output
    void sort()
    {
        DBG(debug, "VarlengthSorter merging " << m_block_seqs.size() << " sequences.");

        m_istream_vector.resize( m_block_seqs.size() );

        for (size_t i = 0; i < m_block_seqs.size(); ++i)
        {
            //DBG(debug, "  seq" << i << " : " << m_block_seqs[i].size() << " blocks of size " << blocksize);

            m_istream_vector[i] = new buf_istream_type( m_block_seqs[i].begin(), m_block_seqs[i].end(), 2 );
        }

        // create comparator for loser tree

        m_streamcmp.resize();

        for (size_t i = 0; i < m_istream_vector.size(); ++i)
        {
            m_streamcmp.fetch(i);
        }

        m_losertree.play_initial( m_block_seqs.size() );
    }

    /// Free all resources
    void clear()
    {
        for (size_t i = 0; i < m_block_seqs.size(); ++i)
        {
            delete m_istream_vector[i];

            m_bm->delete_blocks( m_block_seqs[i].begin(), m_block_seqs[i].end() );
        }

        m_istream_vector.clear();
        m_block_seqs.clear();
    }

    /// stream interface
    bool empty() const
    {
        return m_losertree.done();
    }

    /// stream interface
    value_type& operator*()
    {
        return m_streamcmp.m_head[ m_losertree.top() ];
    }

    /// stream interface
    VarlengthSorter& operator++()
    {
        // get next string in sequence and replay tree
        m_streamcmp.fetch( m_losertree.top() );
        m_losertree.replay();
        return *this;
    }

    /// True if the loser tree indicates that this value is equal to the preceding one. This is however only
    /// semi-definite; due to comparisons in the tree, an equal pair may go unnoticed.
    bool was_equal() const
    {
        return m_losertree.top_equal();
    }
};

// }}}

/*************************************************************************************************************
 * Top-Level Algorithm Class                                                                                 *
 *************************************************************************************************************/

template <typename AlphabetType, typename OffsetType, typename SizeType>
class eSAIS
{
public:

    // *******************************************************************************************************
    // *** Size parameters

    static const size_t block_size = ::block_size;

    static const size_t memsize = ram_use;

    // *******************************************************************************************************
    // *** Global typedefs and constants

    /// type used for character from the string
    typedef AlphabetType        alphabet_type;

    /// type used for offsets and indexes on storage
    typedef OffsetType          offset_type;

    /// type used for offsets and indexes in memory
    typedef SizeType            size_type;

    /// number of characters stored in a PQ tuple
#ifdef ESAIS_TUPLECHARLIMIT
    static const unsigned int   D = ESAIS_TUPLECHARLIMIT;
#else
    static const unsigned int   D = 3;
#endif

    /// type for classes of positions
    enum ctype_type { TYPE_L = 0, TYPE_S = 1 };

    /*********************************************************************************************************
     * General helpers depending on template parameters                                                      *
     *********************************************************************************************************/

    /// helper to print out readable characters
    static inline std::string strC(const alphabet_type& c)
    {
        std::ostringstream oss;
        if (c < alphabet_type(128) && isalnum((size_type)c)) oss << '\'' << (char)((size_type)c) << '\'';
        else oss << (size_type)c;
        return oss.str();
    }

    /// helper to format ctype_type values
    static inline const char* strT(const ctype_type& t)
    {
        if (t == TYPE_L) return "L";
        if (t == TYPE_S) return "S";
        return "?";
    }

    /// helper to format a sequence of characters
    template <typename Iterator>
    static inline std::string strS(Iterator begin, Iterator end)
    {
        std::ostringstream oss;
        if (begin == end) return oss.str();
        if (begin != end) oss << strC(*begin);
        for (Iterator s = begin+1; s != end; ++s)
        {
            oss << " " << strC(*s);
        }
        return oss.str();
    }

    /// helper to format a sequence of characters and calculate their type if the type of the last character
    /// is known.
    template <typename Iterator>
    static inline std::string strST(Iterator begin, Iterator end, ctype_type lasttype)
    {
        std::ostringstream oss;
        if (begin == end) return oss.str();

        std::vector<ctype_type> ctypearray (end - begin);
        {
            typename std::vector<ctype_type>::iterator cti = ctypearray.end();
            Iterator ai = end;

            --cti; --ai; *cti = lasttype;

            while( cti != ctypearray.begin() && ai != begin )
            {
                Iterator pi = ai-1;
                typename std::vector<ctype_type>::iterator pti = cti-1;

                *pti = (*pi < *ai || (*pi == *ai && *cti == TYPE_S)) ? TYPE_S : TYPE_L;

                ai = pi; cti = pti;
            }
        }

        typename std::vector<ctype_type>::iterator cti = ctypearray.begin();
        for (Iterator s = begin; s != end; ++s, ++cti)
        {
            if (s != begin) oss << " ";
            oss << strC(*s) << strT(*cti);
        }
        return oss.str();
    }

    /// helper to format a sequence of characters and calculate their type if all types of the characters are
    /// known.
    template <typename Iterator, typename IteratorCTypes>
    static inline std::string strST2(Iterator begin, Iterator end, IteratorCTypes begin_ctype)
    {
        std::ostringstream oss;
        if (begin == end) return oss.str();
        for (Iterator s = begin; s != end; ++s, ++begin_ctype)
        {
            if (s != begin) oss << " ";
            oss << strC(*s) << strT((ctype_type)*begin_ctype);
        }
        return oss.str();
    }

    /// comparator which compares characters and then the type of the last matching one if either string ends
    /// before the other.
    template <typename Iterator>
    static inline int lexcompare_type_3way(
        Iterator beginA, Iterator endA, ctype_type lasttypeA,
        Iterator beginB, Iterator endB, ctype_type lasttypeB
#if ESAIS_LCP_CALC
        , unsigned int endrepcountA, unsigned int endrepcountB
#endif // ESAIS_LCP_CALC
        )
    {
        assert( endA > beginA );        // no empty strings
        assert( endB > beginB );

        Iterator strA = beginA, strB = beginB;

        // compare characters until either string is empty
        for (; strA != endA && strB != endB; ++strA, ++strB)
        {
            if (*strA < *strB) return -1;
            if (*strB < *strA) return +1;
        }

        //std::cout << "strings [" << strST(beginA,endA,lasttypeA) << "] and [" << strST(beginB,endB,lasttypeB) << "] match chars\n" ;

        if (strA == endA && strB == endB)
        {
            // both strings end at same time.
#if ESAIS_LCP_CALC
            if (lasttypeA != lasttypeB)
                return (lasttypeA - lasttypeB);       // different end types
            else
            {
                // same end types: compare repcounts.
                if (lasttypeA == TYPE_L)
                    return (endrepcountA - endrepcountB);   // longer repcounts are lower
                else
                    return (endrepcountB - endrepcountA);   // longer repcounts are higher
            }
#else // !ESAIS_LCP_CALC
            return (lasttypeA - lasttypeB);       // different end types
#endif // !ESAIS_LCP_CALC
        }
        else if (strA == endA) // && (strB != endB)
        {
            // string A done, string B not yet ended, determine type of prev char of B to compare to lasttypeA

            assert( strB > beginB );
            --strB;

            // determine type of *strB
            while ( strB+1 != endB && *strB == *(strB+1) ) ++strB;

            ctype_type tb = ( strB+1 == endB ) ? lasttypeB : ( *strB < *(strB+1) ? TYPE_S : TYPE_L );

            //std::cout << "type of *strB = " << strT(tb) << "\n";

            if (lasttypeA != tb)
                return (lasttypeA - tb);

            return -1;  // types are equal, shorter string is smaller
        }
        else // (strA != endA) && (strB == endB)
        {
            // string B done, string A not yet ended, determine type of prev char of A to compare to lasttypeB

            assert( strA > beginA );
            --strA;

            // determine type of *strA
            while ( strA+1 != endA && *strA == *(strA+1) ) ++strA;

            ctype_type ta = ( strA+1 == endA ) ? lasttypeA : ( *strA < *(strA+1) ? TYPE_S : TYPE_L );

            //std::cout << "type of *strA = " << strT(ta) << "\n";

            if (ta != lasttypeB)
                return (ta - lasttypeB);

            return +1; // shorter string is smaller
        }
    }

#if ESAIS_LCP_CALC
    /// comparator which compares characters and then the type of the last matching one if either string ends
    /// before the other.
    template <typename Iterator>
    static inline int lexcompare_type_3way_lcp(Iterator beginA, Iterator endA, ctype_type lasttypeA, unsigned int endrepcountA,
                                               Iterator beginB, Iterator endB, ctype_type lasttypeB, unsigned int endrepcountB,
                                               offset_type& lcp)
    {
        assert( endA > beginA );        // no empty strings
        assert( endB > beginB );

        Iterator strA = beginA, strB = beginB;
        lcp = 0;

        // compare characters until either string is empty
        for (; strA != endA && strB != endB; ++strA, ++strB)
        {
            if (*strA < *strB) return -1;
            if (*strB < *strA) return +1;
            ++lcp;
        }

        //std::cout << "strings [" << strST(beginA,endA,lasttypeA) << "] and [" << strST(beginB,endB,lasttypeB) << "] match chars\n" ;

        if (strA == endA && strB == endB)
        {
            // both strings end at same time.

            lcp += std::min(endrepcountA,endrepcountB); // last char is equal, add repcount.

            if (lasttypeA != lasttypeB)
                return (lasttypeA - lasttypeB);       // different end types
            else
                return (endrepcountB - endrepcountA); // same end types: compare repcounts
        }
        else if (strA == endA) // && (strB != endB)
        {
            // string A done, string B not yet ended, must determine type of prev char of B to compare to lasttypeA

            assert( strB > beginB );
            --strB;

            // determine type of *strB and count chars equal to A's last char
            unsigned int repcount = 0;
            while ( strB+1 != endB && *strB == *(strB+1) ) ++strB, ++repcount;

            ctype_type tb = ( strB+1 == endB ) ? lasttypeB : ( *strB < *(strB+1) ? TYPE_S : TYPE_L );

            if (strB+1 == endB) repcount += endrepcountB; // add endrepcount if at end

            //std::cout << "type of *strB = " << strT(tb) << "\n";

            lcp += std::min(repcount, endrepcountA);

            if (lasttypeA != tb)
                return (lasttypeA - tb);

            return -1;  // shorter string is smaller
        }
        else // !(strA == endA) && (strB == endB)
        {
            assert( strA > beginA );
            --strA;

            // determine type of *strA
            unsigned int repcount = 0;
            while ( strA+1 != endA && *strA == *(strA+1) ) ++strA, ++repcount;

            ctype_type ta = ( strA+1 == endA ) ? lasttypeA : ( *strA < *(strA+1) ? TYPE_S : TYPE_L );

            if (strA+1 == endA) repcount += endrepcountA; // add endrepcount if at end

            //std::cout << "type of *strA = " << strT(ta) << "\n";

            lcp += std::min(repcount, endrepcountB);

            if (ta != lasttypeB)
                return (ta - lasttypeB);

            return +1; // shorter string is smaller
        }
    }

    static void test_assert_lexcompare(const std::string& strA, ctype_type lasttypeA, unsigned int endrepcountA,
                                       const std::string& strB, ctype_type lasttypeB, unsigned int endrepcountB,
                                       offset_type checklcp)
    {
        offset_type lcp;

        int c1 = lexcompare_type_3way_lcp(strA.begin(), strA.end(), lasttypeA, endrepcountA,
                                          strB.begin(), strB.end(), lasttypeB, endrepcountB, lcp);
        assert(c1 < 0); assert(lcp == checklcp);

        int c2 = lexcompare_type_3way_lcp(strB.begin(), strB.end(), lasttypeB, endrepcountB,
                                          strA.begin(), strA.end(), lasttypeA, endrepcountA, lcp);
        assert(c2 > 0); assert(lcp == checklcp);
    }

    static void test_lexcompare()
    {
        // (69013,['T'L 'T'L 'G'L 'G'L],L,2) (61227,['T'L 'T'L 'G'S],S,3)
        test_assert_lexcompare("TTGG", TYPE_L, 2,
                               "TTG",  TYPE_S, 3, 6);
    }
#endif // ESAIS_LCP_CALC

    /*********************************************************************************************************
     * Tuple Structures                                                                                      *
     *********************************************************************************************************/

    struct NameTuple
    {
        offset_type     index;          // start index of refered substring
        offset_type     name;           // lexname of substring
#if ESAIS_DISCARD_UNIQUES
        char            unique;         // uniqueness flag
#endif // ESAIS_DISCARD_UNIQUES

        NameTuple() {}

        NameTuple(const offset_type& i, const offset_type& n)
            : index(i), name(n) {}

        static std::ostream& description (std::ostream& os)
        {
            return os << "(index,name"
#if ESAIS_DISCARD_UNIQUES
                      << ",unique"
#endif // ESAIS_DISCARD_UNIQUES
                      << ")";
        }

        friend std::ostream& operator<< (std::ostream& os, const NameTuple& t)
        {
            return os << "(" << t.index << "," << t.name
#if ESAIS_DISCARD_UNIQUES
                      << "," << int(t.unique)
#endif // ESAIS_DISCARD_UNIQUES
                      << ")";
        }
    } __attribute__((packed));

    struct NameTupleOrder_IndexDesc
    {
        inline bool operator()(const NameTuple& a, const NameTuple& b) const
        {
            return a.index > b.index;
        }

        NameTuple min_value() const {
            return NameTuple(std::numeric_limits<offset_type>::max(), 0);
        }
        NameTuple max_value() const {
            return NameTuple(std::numeric_limits<offset_type>::min(), 0);
        }
    };

    struct RankTuple
    {
        offset_type     index;          // start index of refered substring
        offset_type     name;           // lexname or rank of substring

        RankTuple() {}

        RankTuple(const offset_type& i, const offset_type& r)
            : index(i), name(r) {}

        static std::ostream& description (std::ostream& os)
        {
            return os << "(index,name)";
        }

        friend std::ostream& operator<< (std::ostream& os, const RankTuple& t)
        {
            return os << "(" << t.index << "," << t.name << ")";
        }
    };

    struct RankTupleOrder_IndexDesc
    {
        inline bool operator()(const RankTuple& a, const RankTuple& b) const
        {
            return a.index > b.index;
        }

        RankTuple min_value() const {
            return RankTuple(std::numeric_limits<offset_type>::max(), 0);
        }
        RankTuple max_value() const {
            return RankTuple(std::numeric_limits<offset_type>::min(), 0);
        }
    };

#if ESAIS_SELF_CHECK
    struct RankTupleOrder_RankAsc
    {
        inline bool operator()(const RankTuple& a, const RankTuple& b) const
        {
            return a.name < b.name;
        }
    };
#endif // ESAIS_SELF_CHECK

#if ESAIS_DISCARD_UNIQUES
    struct RerankTuple
    {
        offset_type     index;          // start index of refered substring
        offset_type     name;           // lexname of substring
        offset_type     ISA;            // rank in recursive string

        RerankTuple() {}

        RerankTuple(const offset_type& i, const offset_type& n, const offset_type& r)
            : index(i), name(n), ISA(r) {}

        static std::ostream& description (std::ostream& os)
        {
            return os << "(index,name,ISA)";
        }

        friend std::ostream& operator<< (std::ostream& os, const RerankTuple& t)
        {
            return os << "(" << t.index << "," << t.name << "," << t.ISA << ")";
        }
    };

    struct RerankTuple_ISA
    {
        inline bool operator()(const RerankTuple& a, const RerankTuple& b) const
        {
            return a.ISA < b.ISA;
        }

        RerankTuple min_value() const {
            return RerankTuple(0,0,std::numeric_limits<offset_type>::min());
        }
        RerankTuple max_value() const {
            return RerankTuple(0,0,std::numeric_limits<offset_type>::max());
        }
    };
#endif

    /*********************************************************************************************************
     * Step 1 and 2: Sorting short S*-substrings                                                             *
     *********************************************************************************************************/

    // {{{ ShortStringSorter

    /**
     * ShortStringSorter splits the input into S*-substrings and creates Substring tuples, which are
     * serialized and sorted by a VarlengthSorter. The sorted output is then lexicographically named.
     */
    class ShortStringSorter
    {
    public:

        /// Prototype declaration
        struct Substring;

        /// Use the VarlengthSorter for Substrings, which are serialized.
        typedef VarlengthSorter<Substring, block_size>    ssorter_type;

        /// The VarlengthSorter
        ssorter_type            m_ssorter;

        /// Total input size
        size_type               m_inputsize;

        // ***************************************************************************************************
        // *** Substring classes with trailing type indicator

        struct Substring
        {
            offset_type                 index;
            std::vector<alphabet_type>  str;
            ctype_type                  lasttype;
            ESAIS_LCP_CALCX(offset_type endrepcount;)

            inline int cmp(const Substring& b) const
            {
                return lexcompare_type_3way(str.begin(), str.end(), lasttype,
                                            b.str.begin(), b.str.end(), b.lasttype
#if ESAIS_LCP_CALC
                                            ,endrepcount,b.endrepcount
#endif // ESAIS_LCP_CALC
                    );
            }

#if ESAIS_LCP_CALC
            inline int cmp(const Substring& b, offset_type& lcp) const
            {
                return lexcompare_type_3way_lcp(str.begin(), str.end(), lasttype, endrepcount,
                                                b.str.begin(), b.str.end(), b.lasttype, b.endrepcount, lcp);
            }
#endif // ESAIS_LCP_CALC

            inline bool operator< (const Substring& b) const
            {
                return cmp(b) < 0;
            }

            inline bool operator<= (const Substring& b) const
            {
                return cmp(b) <= 0;
            }

            friend std::ostream& operator<< (std::ostream& os, const Substring& s)
            {
                return os << "(" << s.index << ",[" << strST(s.str.begin(), s.str.end(), s.lasttype) << "]," << strT(s.lasttype)
                             ESAIS_LCP_CALCX(<< "," << s.endrepcount) << ")";
            }

            void swap(Substring& o)
            {
                std::swap(index, o.index);
                std::swap(str, o.str);
                std::swap(lasttype, o.lasttype);
                ESAIS_LCP_CALCX( std::swap(endrepcount, o.endrepcount); )
            }

            /// deserialize from byte input stream
            bool deserialize(typename ssorter_type::buf_istream_type& is)
            {
                // read offset_type as index

                stream_get( is, index );
                DBG(debug_substring_decoder, "index = " << index);

                // read varint length and lasttype ctype from sequence
                unsigned int len;
                stream_get_varint(is, len);

                lasttype = (ctype_type)(len & 1);
                len /= 2;

                // read varint S*-endrepcount
                ESAIS_LCP_CALCX( stream_get_varint(is,endrepcount); )

                DBG(debug_substring_decoder, "varint length = " << len << ", ctype = " << strT(lasttype)
                    ESAIS_LCP_CALCX(<< ", endrepcount = " << endrepcount));

                // read string

                str.resize(len);

                for (unsigned int i = 0; i < len; ++i)
                {
                    stream_get( is, str[i] );
                }

                if (debug_substring_decoder)
                {
                    std::cout << "string = ";
                    for (unsigned int i = 0; i < len; ++i)
                    {
                        std::cout << strC( str[i] ) << " ";
                    }
                    std::cout << "\n";
                }

                return ( str.size() != 0 );
            }
        };

        /**
         * When splitting the input into S*-substring, instead of copying characters we use this pointer
         * struct to mark offset and length. The Ptrs are sorted and serialized into the VarlengthSorter.
         */
        struct SubstringPtr
        {
            size_type           offset;
            uint16_t            len;    // no S*-substring is allowed > 64 KiB, otherwise it is split
            ESAIS_LCP_CALCX(offset_type endrepcount);

#if ESAIS_LCP_CALC
            SubstringPtr(size_type _offset, uint16_t _len, const offset_type& _endrepcount)
                : offset(_offset), len(_len), endrepcount(_endrepcount)
            {
                assert( len > 0 );
            }
#else
            SubstringPtr(size_type _offset, uint16_t _len, const offset_type&)
                : offset(_offset), len(_len)
            {
                assert( len > 0 );
            }
#endif
            friend std::ostream& operator<< (std::ostream& os, const SubstringPtr& s)
            {
                return os << "(" << s.offset << "," << s.len << ")";
            }

            Substring getSubstring(size_type baseoff, const alphabet_type* str, const unsigned char* ctype) const
            {
                Substring s;
                s.index = baseoff + offset;
                s.str.assign( str + offset, str + offset + len );
                s.lasttype = (ctype_type)ctype[ offset + len - 1 ];
                ESAIS_LCP_CALCX( s.endrepcount = endrepcount );
                return s;
            }
        };

        /**
         * Sorting Functional for SubstringPtr objects.
         */
        struct SubstringPtrSort
        {
            const alphabet_type*        str;
            const unsigned char*        ctype;

            SubstringPtrSort(const alphabet_type* _str, const unsigned char* _ctype)
                : str(_str), ctype(_ctype)
            {
            }

            // *** Interface for STL sort

            inline bool operator()(const SubstringPtr& a, const SubstringPtr& b) const
            {
                return lexcompare_type_3way(str + a.offset, str + a.offset + a.len, (ctype_type)ctype[ a.offset + a.len-1 ],
                                            str + b.offset, str + b.offset + b.len, (ctype_type)ctype[ b.offset + b.len-1 ]
#if ESAIS_LCP_CALC
                                            , a.endrepcount, b.endrepcount
#endif
                    ) < 0;
            }

            // *** Interface for generic radixsort

            typedef unsigned short oracle_type;

            static size_t limit(size_t /* depth */)
            {
                return 2 * (std::numeric_limits<unsigned char>::max() + 1);
            }

            static const size_t maxdepth = 16;

            oracle_type index(const SubstringPtr& sp, size_t depth) const
            {
                if (depth == sp.len) return 0;
                return 2 * oracle_type(str[sp.offset + depth]) + ctype[sp.offset + depth] + 1;
            }

            typedef typename std::vector<SubstringPtr>::iterator Iterator;

#if ESAIS_LCP_CALC
            static inline bool subsort_cmpL(const SubstringPtr& a, const SubstringPtr& b)
            {
                // for L-types longer repcounts are ranked lower
                return a.endrepcount < b.endrepcount;
            }

            static inline bool subsort_cmpS(const SubstringPtr& a, const SubstringPtr& b)
            {
                // for S-types longer repcounts are ranked higher
                return a.endrepcount > b.endrepcount;
            }
            // thereby the highest repcounts meet at the L/S-seam and correct LCPs are calculated

            // this is called by radixsort when all items in [begin,end) have same char and type, remains to
            // compare repcounts
            void subsort(Iterator begin, Iterator end) const
            {
#if ESAIS_SELF_CHECK
                for (Iterator x = begin+1; x != end; ++x) {
                    // all entries must have save type in the last character
                    assert( ctype[ begin->offset + begin->len-1 ] == ctype[ x->offset + x->len-1 ] );
                }
#endif // ESAIS_SELF_CHECK

                if ( ctype[ begin->offset + begin->len-1 ] == TYPE_L )
                    std::sort(begin, end, subsort_cmpL);
                else
                    std::sort(begin, end, subsort_cmpS);
            }
#else
            static void subsort(Iterator, Iterator) { }
#endif
        };

        /*
         * Serialization functional for serializing SubstringPtr into the VarlengthSorter. The serialization
         * is then read as a Substring object.
         */ 
        struct SubstringPtrWriter
        {
            size_type baseoffset;               // offset of all substrings
            alphabet_type* buffer;              // the current portion of the input
            unsigned char* buffer_types;        // the character types of the current portion

            inline void operator()(typename ssorter_type::buf_ostream_type& os, const SubstringPtr& sp) const
            {
                // write offset_type index
                offset_type index = baseoffset + sp.offset;
                stream_put( os, index );

                // prepare ctype bit of last character
                unsigned char ctype = buffer_types[ sp.offset + sp.len - 1 ];
                assert( ctype == 0 || ctype == 1 );

                // write varint string size: first character contains ctype as second-highest bit
                assert( sp.len != 0 );
                stream_put_varint( os, sp.len * 2 + ctype );

                // write varint S*-repcount
                ESAIS_LCP_CALCX( stream_put_varint( os, sp.endrepcount ); )

                // write string character-wise
                for (size_type i = sp.offset; i < sp.offset + sp.len; ++i)
                {
                    stream_put( os, buffer[i] );
                }
            }
        };

        /// Left-over from experiments with 8-bit radixsort
        struct RadixsortTransform1
        {
            static size_t limit(size_t depth)
            {
                return 2 * 0x100 + 1;
            }

            static inline unsigned int index(const SubstringPtr& p, size_t depth)
            {
                return p.radixsort_index(depth);
            }
        };

        /// Left-over from experiments with 16-bit radixsort
        struct RadixsortTransform2
        {
            static size_t limit(size_t depth)
            {
                return (depth % 2 == 1) ? (2 * 0x10000) + 1 : 0x10000 + 1;
            }

            static inline unsigned int index(const SubstringPtr& p, size_t depth)
            {
                if (depth / 2 == p.len) return 0;

                if (depth % 2 == 1)
                    return 2 * (p.str[depth / 2] & 0xFFFF) + p.ctypeptr[depth / 2] + 1;
                else
                    return ((p.str[depth / 2] >> 16) & 0xFFFF) + 1;
            }
        };

        // ***************************************************************************************************
        // *** Read input (in reverse) and split into S*- or split substrings.

#if ESAIS_LCP_CALC
        typedef stxxl::queue< uint16_t, block_size >         SStarSize_type;
        SStarSize_type          SStarSize;
#endif // ESAIS_LCP_CALC

        template <typename InputStreamReverse>
        void process_input(InputStreamReverse& inputrev, unsigned int depth)
        {
            // *** First scan the input back to front to find all S* positions

            stxxl::queue<offset_type,block_size>     sstar_positions;

            {
                inputrev.rewind(memsize);

                offset_type spos = (size_type)inputrev.size()-1;   // sentinel is not processed
                ctype_type prev_ctype = TYPE_L;         // last char is always L-type
                alphabet_type prev_char = *inputrev;

                DBG(debug_substring_split, "input[" << spos << "]: " << strC(prev_char) << " - " << strT(prev_ctype));

                --inputrev;

                while ( spos > offset_type(0) )
                {
                    assert( !inputrev.empty() );

                    ctype_type this_ctype = (*inputrev < prev_char || (*inputrev == prev_char && prev_ctype == TYPE_S)) ? TYPE_S : TYPE_L;

                    DBG(debug_substring_split, "input[" << spos-1 << "]: " << strC(*inputrev) << " - " << strT(this_ctype));

                    if (prev_ctype == TYPE_S && this_ctype == TYPE_L)
                    {
                        DBG(debug_substring_split, "S* at position " << spos);

                        sstar_positions.push(spos);
                    }

                    prev_ctype = this_ctype;
                    prev_char = *inputrev;

                    --spos;
                    --inputrev;
                }

                DBG(debug, "Total S*-indexes in input: " << sstar_positions.size());
                g_statscache >> "S*-indexes" >> depth << sstar_positions.size();
            }

            // *** Rescan input from back to front and build substring (pointer,size) pairs

            inputrev.rewind(memsize / 4);

            // value for for block splitting: at most 64kb lengths (so length can be 2 bytes)
            static const unsigned int split_size = std::min<unsigned int>(std::numeric_limits<uint16_t>::max(),block_size);

            // end pointer of current area
            size_type sarea_end = m_inputsize = inputrev.size();

            // for each position of buffer at most one char, 1 ctype and a half-StringPtr are needed.
            static const size_type per_substring_size = sizeof(alphabet_type) + sizeof(unsigned char) + sizeof(SubstringPtr) / 2;

            // for depth=0 we also need buffers for radixsort: K+1 * bucketsize + n * oracle_type
            static const size_type buffersize =
                depth == 0
                ? ( (memsize - (SubstringPtrSort::maxdepth+1) * SubstringPtrSort::limit(0) * sizeof(size_t))
                    / (per_substring_size + sizeof(typename SubstringPtrSort::oracle_type)) )
                : ( memsize * 3/4 / per_substring_size );

            DBGMEM("Prior to S*-substring buffer allocation");

            // temporary buffer for characters part of the current set of substrings
            stxxl::simple_vector<alphabet_type> buffer (buffersize);

            // temporary buffer for types of characters of the current set of substrings
            stxxl::simple_vector<unsigned char> buffer_types (buffersize);

            // pairs of (index,size) to sort split substrings in buffer
            std::vector< SubstringPtr > substrings;
            substrings.reserve( buffersize / 2 );

            DBGMEM("After buffer allocation");

            ctype_type prev_ctype = TYPE_L;     // start with sentinel, but here label it L to skip sentinel S*-substring
            alphabet_type prev_char = 0;
            size_type sarea_pos = sarea_end;

            size_type total_substrings = 0;     // total S*-substrings

            unsigned int overlap = 0;           // carry-over characters from last portion of input

#if ESAIS_LCP_CALC
            size_type repcount = 0;
#endif // ESAIS_LCP_CALC
            size_type prev_sstar_repcount = 0;

            while ( sarea_end > 0 && !sstar_positions.empty() )
            {
                // fetch new input portion of size M

                size_type sarea_begin = (sarea_end >= (size_type)buffersize ? sarea_end - buffersize : 0);
                size_type sarea_size = sarea_end - sarea_begin;

                // copy overlap to end of area
                if (overlap)
                {
                    DBG(debug_substring_split, "Copy overlap of " << overlap << " from beginning to end " << sarea_size);

                    std::copy( buffer.begin(), buffer.begin() + overlap, buffer.begin() + sarea_size - overlap );
                    std::copy( buffer_types.begin(), buffer_types.begin() + overlap, buffer_types.begin() + sarea_size - overlap );
                }

                DBG(debug, "working on M range = " << sarea_begin << " to " << sarea_end << ", overlap " << overlap);

                // length of current string
                size_type currstringlength = overlap;

                // total bytes that will be written in the end for this buffer
                size_type outputsize = 0;

                substrings.clear();

                for (size_type i = sarea_size - overlap; i > 0 && !sstar_positions.empty();)
                {
                    --i; --sarea_pos;

                    assert( sarea_pos == sarea_begin+i );

                    // read next character
                    buffer[i] = *inputrev;
                    --inputrev;

                    currstringlength++;

                    ctype_type this_ctype = (buffer[i] < prev_char || (buffer[i] == prev_char && prev_ctype == TYPE_S)) ? TYPE_S : TYPE_L;
                    buffer_types[i] = this_ctype;

#if ESAIS_LCP_CALC
                    repcount = (buffer[i] == prev_char) ? repcount+1 : 0;
#endif // ESAIS_LCP_CALC

                    DBG(debug_substring_split, "next input: " << strC( buffer[i] ) << " at " << sarea_begin + i << " of ctype " << strT(this_ctype)
                        ESAIS_LCP_CALCX(<< " with repcount " << repcount) << " and S*-distance: " << (sarea_pos - (size_type)sstar_positions.front()));

                    assert( sarea_pos >= (size_type)sstar_positions.front() );

                    if ( (sarea_pos - (size_type)sstar_positions.front()) % (split_size-1) == 0 ) // S* or block splitting boundary
                    {
                        // split at S*-position or block boundary

                        if (sarea_pos == (size_type)sstar_positions.front())
                        {
                            DBG1(debug_substring_split, "S*-substr: ");
                            sstar_positions.pop();
                        }
                        else
                        {
                            DBG1(debug_substring_split, "block-split: ");
                        }
#if ESAIS_LCP_CALC
                        assert(currstringlength > 0);
                        // save for LCP_S* calculation - deduct overlapping (usually S*) position at end,
                        // except for first tuple where the S* sentinel is only virtual.
                        SStarSize.push(currstringlength - (total_substrings == 0 ? 0 : 1));
#endif
                        assert( currstringlength <= std::numeric_limits<uint16_t>::max() );
                        substrings.push_back( SubstringPtr(i, currstringlength, prev_sstar_repcount) );

                        outputsize += sizeof(offset_type)
                            + sizeof_varint(currstringlength*2)
                            ESAIS_LCP_CALCX(+ sizeof_varint(prev_sstar_repcount))
                            + currstringlength * sizeof(alphabet_type);

                        ++total_substrings;

                        DBG2(debug_substring_split,
                             i << " = ["  << sarea_pos << "," << sarea_pos + currstringlength << ") = "
                             << substrings.back().getSubstring(sarea_pos,buffer.data(),buffer_types.data())
                             ESAIS_LCP_CALCX(<< " with sstar_repcount=" << prev_sstar_repcount) << "\n");

                        currstringlength = 1;
                        ESAIS_LCP_CALCX( prev_sstar_repcount = repcount; )
                    }

                    prev_char = buffer[i];
                    prev_ctype = this_ctype;
                }

                assert( substrings.size() <= buffersize/2 );

                DBG(debug_substring_split, "Remaining: ["  << sarea_begin << "," << sarea_begin + currstringlength << ") = "
                    << strS( &buffer[0], &buffer[currstringlength] ));

                // *** Sort SubstringPtr using radix or quicksort

                SubstringPtrSort sorter (buffer.data(), buffer_types.data());
                if ( sizeof(alphabet_type) == 1 )
                    radixsort_transform_oracle( substrings.begin(), substrings.end(), sorter );
                else
                    std::sort( substrings.begin(), substrings.end(), sorter );

#if !ESAIS_SELF_CHECK
                if (debug_substring_split)
#endif // !ESAIS_SELF_CHECK
                {
                    Substring prev = substrings[0].getSubstring(sarea_pos,buffer.data(),buffer_types.data());
                    DBG(debug_substring_split, "substr: " << prev);
                    for (size_t i = 1; i < substrings.size(); ++i)
                    {
                        Substring curr = substrings[i].getSubstring(sarea_pos,buffer.data(),buffer_types.data());
                        DBG(debug_substring_split, "substr: " << curr);
                        if( !(prev <= curr) ) DBG(debug, "WRONG substring order!");
                        prev = curr;
                    }
                }

                // at least one zero string at the end
                outputsize += sizeof(offset_type) + sizeof(unsigned char) + 4 * sizeof(alphabet_type);

                DBG(debug, "outputsize = " << outputsize);

                SubstringPtrWriter spw;
                spw.baseoffset = sarea_begin;
                spw.buffer = buffer.data();
                spw.buffer_types = buffer_types.data();
                m_ssorter.write( outputsize, substrings.begin(), substrings.end(), spw );

                DBG(debug_substring_split, "Reading new buffer of size M");

                overlap = currstringlength;
                sarea_end = sarea_begin + overlap;
            }

            //assert(inputrev.empty());
            assert(sstar_positions.empty());

            DBG(debug, "Total number of substrings " << total_substrings);

            inputrev.finish();
        }

    protected:

        /// flag whether every name was unique
        bool                    m_was_unique;

        /// total number of indexed substring positions
        size_type               m_totalsize;

        /// number of substrings that are duplicates
        size_type               m_duplicates;

    public:

#if ESAIS_LCP_CALC
        /// Save LCP_Names, the LCP of two consecutive lexnames, for answering RMQs later. One would expect
        /// this to be only uint16_t, but due to repcount being included in the LCP, we actually need an
        /// offset_type.
        typedef stxxl::sequence< offset_type, block_size >         lcp_seq_type;
        lcp_seq_type          lexname_lcp_seq;
#endif // ESAIS_LCP_CALC

        // ***************************************************************************************************
        // *** Read sorted sequences, merge using loser tree and create lexname pairs (index,name)

        template <typename OutputNameTupleStream>
        void merge(OutputNameTupleStream& outputnametuples)
        {
            m_ssorter.sort();

            size_type name = 0;
            m_was_unique = true;
            m_totalsize = 0;
            m_duplicates = 0;

            if (m_ssorter.empty()) return;

            // initialize prev variable with first item

            Substring prevstring = *m_ssorter;
            NameTuple prevtuple (prevstring.index, name);

#if ESAIS_DISCARD_UNIQUES
            prevtuple.unique = true;
#else // !ESAIS_DISCARD_UNIQUES
            bool prevunique = false;    // for counting duplicates and uniques
#endif

#if ESAIS_LCP_CALC
            offset_type prevlcp = 0, lcp = 0;
            lexname_lcp_seq.push_back(0);     // actually undefined
#endif // ESAIS_LCP_CALC

            ++m_totalsize;
            ++m_ssorter;

            while( !m_ssorter.empty() )
            {
                Substring& topstr = *m_ssorter;

                // check order of strings extracted from sorter
                assert( prevstring.cmp(topstr) <= 0 );

                // check if this string is different from last one.
                bool is_different = false;

                if ( m_ssorter.was_equal() ) {
                    // comparison not necessary: prev and this are equal by comparision within loser tree
                    is_different = false;
                    ESAIS_LCP_CALCX( lcp = topstr.str.size() + topstr.endrepcount );
                }
                else {
#if !ESAIS_LCP_CALC
                    is_different = ( prevstring.cmp(topstr) != 0 );
#else // ESAIS_LCP_CALC
                    is_different = ( prevstring.cmp(topstr, lcp) != 0 );
#endif // ESAIS_LCP_CALC
                }

                ESAIS_LCP_CALCX( lexname_lcp_seq.push_back(lcp) );       // save lcp in lexname order

                // select name or update unique flags
                if (is_different) {
#if ESAIS_DISCARD_UNIQUES
                    name = m_totalsize;            // select next name
#else
                    name++;
#endif
                }
                else {
                    m_was_unique = false;                       // global unique flag
#if ESAIS_DISCARD_UNIQUES
                    m_duplicates++;                             // this one is a duplicate
                    if (prevtuple.unique) ++m_duplicates;       // one unique less.
                    prevtuple.unique = false;                   // mark previous item as duplicate to this one
#else
                    m_duplicates++;                             // this one is a duplicate
                    if (prevunique) ++m_duplicates;             // one unique less.
#endif // ESAIS_DISCARD_UNIQUES
                }

                // push previous tuple
                DBG(debug_substring_merge, "tuple " << prevtuple << " - substring = " << prevstring
                    ESAIS_LCP_CALCX(<< " with lcp " << prevlcp));
                outputnametuples.push( prevtuple );

                // create next tuple from current string
                prevtuple = NameTuple(topstr.index, name);
#if ESAIS_DISCARD_UNIQUES
                prevtuple.unique = is_different;
#endif // ESAIS_DISCARD_UNIQUES
#if ESAIS_LCP_CALC
                prevlcp = lcp;
#endif // ESAIS_LCP_CALC

                // advance
                m_totalsize++;
                prevstring.swap(topstr);
#if !ESAIS_DISCARD_UNIQUES
                prevunique = is_different;
#endif // !ESAIS_DISCARD_UNIQUES

                ++m_ssorter;
            }

            if (m_totalsize != 0) {
                DBG(debug_substring_merge, "tuple " << prevtuple << " - substring = " << prevstring
                    ESAIS_LCP_CALCX(<< " with lcp " << prevlcp));
                outputnametuples.push( prevtuple );
            }

            m_ssorter.clear();  // free memory and allocated blocks
        }

        bool was_unique() const
        {
            return m_was_unique;
        }

        size_type totalsize() const
        {
            return m_totalsize;
        }

        size_type duplicates() const
        {
            return m_duplicates;
        }
    };

    // }}}

    // *******************************************************************************************************
    // *** Step 4: Inducing the SA using a PQ

    // {{{ Induce Result Structure: Form a Suffix Array using Queue and Stack

    /* eSAIS's induce step outputs all L-type chars followed by all S-type chars due to the different induce
      orders. The result object is used to merge these streams of different types back together into a
      consequtive suffix array. However, since there are no bucket counters as the alphabet can be very large,
      the two streams L-type and S-types are saved indepdendent from another into a queue and stack,
      respectively. Also recorded are the number of L-types of the same character, and which character that
      is. In the output phase, the queue and stack are merged by character yielding the consecutive suffix
      array. The lcp value is stored along each suffix array index if calculated.
     */
    struct Result
    {
    private:
        /// represents a "run" of suffix array indexes starting with the same character in the output.
        typedef std::pair<alphabet_type, offset_type>           result_run_type;

        typedef stxxl::queue< offset_type, block_size >         Lresult_queue_type;
        typedef stxxl::queue< result_run_type, block_size >     Lresult_run_queue_type;

        typedef typename stxxl::STACK_GENERATOR<offset_type, stxxl::external, stxxl::grow_shrink, 4, block_size>::result Sresult_stack_type;
        typedef typename stxxl::STACK_GENERATOR<result_run_type, stxxl::external, stxxl::grow_shrink, 4, block_size>::result Sresult_run_stack_type;

        Lresult_queue_type              Lqueue;
#if ESAIS_LCP_CALC
        Lresult_queue_type              Lqueue_lcp;
#endif // ESAIS_LCP_CALC
        Lresult_run_queue_type          Lrun_queue;

        Sresult_stack_type              Sstack;
#if ESAIS_LCP_CALC
        Sresult_stack_type              Sstack_lcp;
#endif // ESAIS_LCP_CALC
        Sresult_run_stack_type          Srun_stack;

        ctype_type                      LSstate;
        size_type                       LScounter;

        // ***************************************************************************************************
        // *** Input into Queue or Stack

    private:
        // Counter for output number
        size_type           output_count;
        alphabet_type       output_char;

    public:
        Result()
            : output_count(0), output_char( std::numeric_limits<alphabet_type>::min() )
        {
        }

        void output_Lentry(const alphabet_type& char0, const offset_type& index)
        {
            DBG(debug_output_write, "L-SAEntry: char " << strC(char0) << " index " << index);

            // save index into Lresult queue and if needed the size of the previous sequence
            if (output_char != char0)
            {
                if (output_count > 0)
                    Lrun_queue.push( result_run_type(output_char,output_count) );

                output_char = char0;
                output_count = 0;
            }

            Lqueue.push(index);
            output_count++;
        }

#if ESAIS_LCP_CALC
        void output_Lentry_lcp(const offset_type& lcp)
        {
            DBG(debug_output_write, "L-LCPEntry: lcp " << lcp);
            Lqueue_lcp.push(lcp);
        }
#endif // ESAIS_LCP_CALC

        void finish_Lsequence()
        {
            // no string can contain no L-types
            assert( output_count > 0 );
            Lrun_queue.push( result_run_type(output_char,output_count) );

            output_count = 0;
            output_char = std::numeric_limits<alphabet_type>::max();
        }

        void output_Sentry(const alphabet_type& char0, const offset_type& index)
        {
            DBG(debug_output_write, "S-SAEntry: char " << strC(char0) << " index " << index);

            // save index into Sresult stack and if needed the size of the previous sequence
            if (output_char != char0)
            {
                if (output_count > 0)
                    Srun_stack.push( result_run_type(output_char, output_count) );

                output_char = char0;
                output_count = 0;
            }

            Sstack.push(index);
            output_count++;
        }

#if ESAIS_LCP_CALC
        void output_Sentry_lcp(const offset_type& lcp)
        {
            DBG(debug_output_write, "S-LCPEntry: lcp " << lcp);
            Sstack_lcp.push(lcp);
        }
#endif // ESAIS_LCP_CALC

        void finish_Ssequence()
        {
            if ( output_count > 0 ) // A string can actually contain no S-chars: 'aaaa'
                Srun_stack.push( result_run_type(output_char,output_count) );

            // initialize output sequence
            LScounter = 0;
            operator++();
        }

    public:
        // ***************************************************************************************************
        // *** Output SA by alternating between Lqueue and Sstack

        Result& operator++()
        {
            if (LScounter > 0)   // remove current item
            {
                LScounter--;
#if ESAIS_LCP_CALC
                (LSstate == TYPE_L) ? (Lqueue.pop(), Lqueue_lcp.pop()) : (Sstack.pop(), Sstack_lcp.pop());
#else
                (LSstate == TYPE_L) ? Lqueue.pop() : Sstack.pop();
#endif // ESAIS_LCP_CALC
            }

            if (LScounter == 0)
            {
                if ( Lrun_queue.empty() && Srun_stack.empty() )
                {
                    assert( Lqueue.empty() );
                    assert( Sstack.empty() );
#if ESAIS_LCP_CALC
                    assert( Lqueue_lcp.empty() );
                    assert( Sstack_lcp.empty() );
#endif // ESAIS_LCP_CALC

                    DBG(debug_output_merge, "very empty");

                    // empty
                    return *this;
                }
                else if ( Srun_stack.empty() )
                {
                    DBG(debug_output_merge, "Finishing L-sequence of char " << strC(Lrun_queue.front().first)
                        << " of length " << Lrun_queue.front().second);

                    LScounter = (size_type)Lrun_queue.front().second;
                    Lrun_queue.pop();
                    LSstate = TYPE_L;
                }
                else if ( Lrun_queue.empty() )
                {
                    DBG(debug_output_merge, "Finishing S-sequence of char " << strC(Srun_stack.top().first)
                        << " of length " << Srun_stack.top().second);

                    LScounter = (size_type)Srun_stack.top().second;
                    Srun_stack.pop();
                    LSstate = TYPE_S;
                }
                else if ( Lrun_queue.front().first <= Srun_stack.top().first )
                {
                    DBG(debug_output_merge, "Switching to L-sequence of char " << strC(Lrun_queue.front().first)
                        << " of length " << Lrun_queue.front().second);

                    LScounter = (size_type)Lrun_queue.front().second;
                    Lrun_queue.pop();
                    LSstate = TYPE_L;
                }
                else
                {
                    DBG(debug_output_merge, "Switching to S-sequence of char " << strC(Srun_stack.top().first)
                        << " of length " << Srun_stack.top().second);

                    LScounter = (size_type)Srun_stack.top().second;
                    Srun_stack.pop();
                    LSstate = TYPE_S;
                }
            }

            return *this;
        }

        bool empty() const
        {
            return (Lqueue.empty() && Sstack.empty());
        }

        typedef offset_type     value_type;

        const value_type& operator*() const
        {
            assert(LScounter > 0);

            if (LSstate == TYPE_L)
            {
                DBG(debug_output_merge, "SA: " << Lqueue.front());
                return Lqueue.front();
            }
            else
            {
                DBG(debug_output_merge, "SA: " << Sstack.top());
                return Sstack.top();
            }
        }

#if ESAIS_LCP_CALC
        const value_type& lcp() const
        {
            assert(LScounter > 0);

            if (LSstate == TYPE_L)
            {
                DBG(debug_output_merge, "LCP: " << Lqueue_lcp.front());
                return Lqueue_lcp.front();
            }
            else
            {
                DBG(debug_output_merge, "LCP: " << Sstack_lcp.top());
                return Sstack_lcp.top();
            }
        }
#endif // ESAIS_LCP_CALC
    };

    // }}} Induce Result Structure: Array

    // {{{ Induce SA using PQ

    class Induce
    {
    public:
        // memsize divisors for induce() phase depends on extra LCP information:

#if ESAIS_LCP_CALC
        static const size_t mempartL = memsize / 8;
        static const size_t mempartS = memsize / 7;
#else
        static const size_t mempartL = memsize / 6;
        static const size_t mempartS = memsize / 5;
#endif // ESAIS_LCP_CALC

        static const offset_type lcp_unknown() { return std::numeric_limits<offset_type>::max(); }

        // ***************************************************************************************************
        // *** Tuples and Comparator for the L-PQ and S-PQ

        /// Continuation tuple -- both L-Array and S-Array contain CTuples
        struct CTuple
        {
            offset_type         index;          // i
            offset_type         repcount;       // repetition counter
            alphabet_type	chars[D];	// T[i,i-1,...,i-D+1]
            unsigned char       charfill : 7;
            unsigned char       continued : 1;

            static std::ostream& description (std::ostream& os)
            {
                return os << "(index,chars,continued,repcount)";
            }

            friend std::ostream& operator<< (std::ostream& os, const CTuple& t)
            {
                return os << "(" << t.index << ",[" << strS(t.chars, t.chars + t.charfill) << "],"
                          << (t.continued ? 'C' : '_') << "," << t.repcount << ")";
            }
        } __attribute__((packed));

        /// Order continuation tuples in L-Array by (chars[0] ASC, repcount ASC, index ASC)
        struct CTupleOrder_LArray
        {
            inline bool operator() (const CTuple& a, const CTuple& b) const
            {
                if (a.chars[0] == b.chars[0]) {
                    if (a.repcount == b.repcount) {
                        return (a.index < b.index);
                    }
                    return (a.repcount < b.repcount);
                }
                return (a.chars[0] < b.chars[0]);
            }

            static CTuple min_value()
            {
                CTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::min();
                t.repcount = std::numeric_limits<offset_type>::min();
                t.index = std::numeric_limits<offset_type>::min();
                return t;
            }
            static CTuple max_value()
            {
                CTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::max();
                t.repcount = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        /// Order continuation tuples in S-Array by (chars[0] DESC, repcount ASC, index ASC)
        struct CTupleOrder_SArray
        {
            inline bool operator() (const CTuple& a, const CTuple& b) const
            {
                if (a.chars[0] == b.chars[0]) {
                    if (a.repcount == b.repcount) {
                        return (a.index < b.index);
                    }
                    return (a.repcount < b.repcount);
                }
                return (a.chars[0] > b.chars[0]);
            }

            static CTuple min_value()
            {
                CTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::max();
                t.repcount = std::numeric_limits<offset_type>::min();
                t.index = std::numeric_limits<offset_type>::min();
                return t;
            }
            static CTuple max_value()
            {
                CTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::min();
                t.repcount = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        /// Tuples of the S*-Array
        struct STuple
        {
            offset_type         index;          // i
            offset_type         rank;		// ISA[i+1]
            alphabet_type       chars[D];       // T[i,i-1,...,i-D+1]
            unsigned char       charfill : 7;
            unsigned char       continued : 1;  // must look for continuation tuple when charfill=1
#if ESAIS_LCP_CALC
            offset_type         repcount;       // repetition count of char[0], the S*-index (zero for first char)
#endif
            static std::ostream& description (std::ostream& os)
            {
                return os << "(index,chars,rank,continued" ESAIS_LCP_CALCX(<< ",repcount") << ")";
            }

            friend std::ostream& operator<< (std::ostream& os, const STuple& t)
            {
                return os << "(" << t.index
                          << ",[" << strS(t.chars, t.chars + t.charfill) << "],"
                          << t.rank << "," << (t.continued ? 'C' : '_')
                          ESAIS_LCP_CALCX( << "," << t.repcount )
                          << ")";
            }

        } __attribute__((packed));

        /// Order seed tuples in S*-Array by (chars[0] ASC, rank ASC, index ASC)
        struct STupleOrder_SStarArray
        {
            inline bool operator() (const STuple& a, const STuple& b) const
            {
                if (a.chars[0] == b.chars[0])
                {
                    if (a.rank == b.rank)
                    {
                        return (a.index < b.index);
                    }
                    return (a.rank < b.rank);
                }
                return (a.chars[0] < b.chars[0]);
            }

            static STuple min_value()
            {
                STuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::min();
                t.rank = std::numeric_limits<offset_type>::min();
                t.index = std::numeric_limits<offset_type>::min();
                return t;
            }
            static STuple max_value()
            {
                STuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::max();
                t.rank = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        /// Tuples of the L-PQ and S-PQ, and of the L*-Array
        struct PQTuple
        {
            offset_type         index;          // i
            offset_type         rank;		// ISA[i+1]
            alphabet_type	chars[D];	// T[i,i-1,...,i-D+1]
            unsigned char       charfill : 7;
            unsigned char       continued : 1;  // must look for continuation tuple when charfill=1
#if ESAIS_LCP_CALC_INT
            offset_type         lcp;            // lcp precalculated or max() if unknown
#endif // ESAIS_LCP_CALC_INT

            static PQTuple fromCTuple(const CTuple& ct)
            {
                PQTuple pt;
                pt.index = ct.index;
#ifndef NDEBUG
                pt.rank = -1;
#endif
                memcpy(pt.chars, ct.chars, ct.charfill * sizeof(alphabet_type));
                pt.charfill = ct.charfill;
                pt.continued = ct.continued;
                return pt;
            }

            static PQTuple fromSTuple(const STuple& st)
            {
                PQTuple pt;
                pt.index = st.index;
                pt.rank = st.rank;
                memcpy(pt.chars, st.chars, st.charfill * sizeof(alphabet_type));
                pt.charfill = st.charfill;
                pt.continued = st.continued;
                ESAIS_LCP_CALCX_INT( pt.lcp = lcp_unknown(); )
                return pt;
            }

            static std::ostream& description (std::ostream& os)
            {
                return os << "(index,chars,rank,continued" ESAIS_LCP_CALCX_INT(<< ",lcp") << ")";
            }

            friend std::ostream& operator<< (std::ostream& os, const PQTuple& t)
            {
                os << "(" << t.index
                   << ",[" << strS(t.chars, t.chars + t.charfill) << "],"
                   << t.rank << ","
                   << (t.continued ? 'C' : '_');
#if ESAIS_LCP_CALC_INT
                if (t.lcp != lcp_unknown()) os << ",lcp=" << t.lcp; else os << ",nolcp";
#endif // ESAIS_LCP_CALC_INT
                return os << ")";
            }

            void decrease()
            {
                assert( charfill > 0 );

                --index;
                --charfill;
                for (unsigned int i = 0; i < charfill; ++i)
                {
                    chars[i] = chars[i+1];
                }
            }
        } __attribute__((packed));

        /// Order items in L-PQ by (chars[0] DESC, star DESC, rank DESC, index DESC) as PQ picks smallest element
        struct PQTupleOrder_LPQ
        {
            inline bool operator() (const PQTuple& a, const PQTuple& b) const
            {
                if (a.chars[0] == b.chars[0])
                {
                    if (a.rank == b.rank)
                    {
                        return (a.index > b.index);
                    }
                    return (a.rank > b.rank);
                }
                return (a.chars[0] > b.chars[0]);
            }

            PQTuple min_value() const
            {
                PQTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::max();
                t.rank = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        /// Order items in S-PQ by (chars[0] ASC, star DESC, rank DESC, index DESC) as PQ picks smallest element
        struct PQTupleOrder_SPQ
        {
            inline bool operator() (const PQTuple& a, const PQTuple& b) const
            {
                if (a.chars[0] == b.chars[0])
                {
                    if (a.rank == b.rank)
                    {
                        return (a.index > b.index);
                    }
                    return (a.rank > b.rank);
                }
                return (a.chars[0] < b.chars[0]);
            }

            PQTuple min_value() const
            {
                PQTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::min();
                t.rank = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        /// Order items in L*-Array by (chars[0] DESC, rank DESC, index ASC)
        struct PQTupleOrder_LStarArray
        {
            inline bool operator() (const PQTuple& a, const PQTuple& b) const
            {
                if (a.chars[0] == b.chars[0])
                {
                    if (a.rank == b.rank)
                    {
                        if (a.index == b.index)
                            return (a.charfill < b.charfill);

                        return (a.index < b.index);
                    }

                    return (a.rank > b.rank);
                }
                return (a.chars[0] > b.chars[0]);
            }

            PQTuple min_value() const
            {
                PQTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::max();
                t.rank = std::numeric_limits<offset_type>::min();
                t.index = std::numeric_limits<offset_type>::min();
                t.charfill = 0;
                return t;
            }

            PQTuple max_value() const
            {
                PQTuple t;
                t.chars[0] = std::numeric_limits<alphabet_type>::min();
                t.rank = std::numeric_limits<offset_type>::max();
                t.index = std::numeric_limits<offset_type>::max();
                t.charfill = 0x3F;
                return t;
            }
        };

        /// Tuples in continuation buffers, they are accumulated, sorted and merged with the L/S-Array
        struct CBufferTuple
        {
            offset_type         index;          // i
#ifndef NDEBUG
            alphabet_type       char0;          // chars[0] from discontinued tuple
#endif
            offset_type         rank;           // rank[i+1]

#if ESAIS_LCP_CALC_INT
            offset_type         lcp;            // have lcp precalculated - use max() as sentinel
#endif // ESAIS_LCP_CALC_INT

            static CBufferTuple fromPQTuple(const PQTuple& pt)
            {
                CBufferTuple ct;
                ct.index = pt.index;
#ifndef NDEBUG
                ct.char0 = pt.chars[0];
#endif
                ct.rank = pt.rank;
                ESAIS_LCP_CALCX_INT( ct.lcp = pt.lcp; )
                return ct;
            }

            static std::ostream& description (std::ostream& os)
            {
                return os << "(index,rank" << ESAIS_LCP_CALCX_INT(",lcp" <<) ")";
            }

            friend std::ostream& operator<< (std::ostream& os, const CBufferTuple& t)
            {
                return os << "(" << t.index << "," << t.rank << ESAIS_LCP_CALCX_INT("," << t.lcp <<) ")";
            }
        };

        /// Order continuation buffer tuples in L-MergeBuffer by (index ASC)
        struct CBufferTupleOrder_Index
        {
            inline bool operator() (const CBufferTuple& a, const CBufferTuple& b) const
            {
                return (a.index < b.index);
            }

            static CBufferTuple min_value()
            {
                CBufferTuple t;
                t.index = std::numeric_limits<offset_type>::min();
                return t;
            }
            static CBufferTuple max_value()
            {
                CBufferTuple t;
                t.index = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

#if ESAIS_LCP_CALC_EXT
        /// LCP-Tuples for creating LCP by solving a dynamic RMQ problem. The order of these types is
        /// important as processing order!
        enum LCPTupleType {
            LCP_LSTAR_MARKER,                                           // L* marker for target
            LCP_SETTER, LCP_SETTER_OUTPUT, LCP_SETTER_OUTPUT_LSTAR,     // Sets target to v1 and possibly output it
            LCP_LS_SEAM_QUERY,                                          // Sets target to L/S-seam repcount  of v2=char
            LCP_QUERY, LCP_QUERY_LSTAR                                  // RMQ(v1,target) -> v2
        };

        /// Tuples in LCP priority queue
        struct LCPTuple
        {
            unsigned char       type;           // tuple type
            offset_type         target;         // target, by which tuples are sorted: setter=position, query=right boundary
            offset_type         v1;             // setter=value, query=left boundary, L/S=repcount
            offset_type         v2;		// setter=nothing, query=target, L/S=char

            static std::ostream& description (std::ostream& os)
            {
                return os << "(type,target,v1,v2)";
            }

            friend std::ostream& operator<< (std::ostream& os, const LCPTuple& t)
            {
                if (t.type == LCP_LSTAR_MARKER)
                    return os << "(lcp[" << t.target << "] is L*)";
                else if (t.type == LCP_SETTER)
                    return os << "(lcp[" << t.target << "] := " << t.v1 << ")";
                else if (t.type == LCP_SETTER_OUTPUT)
                    return os << "(lcp[" << t.target << "] := " << t.v1 << " +output)";
                else if (t.type == LCP_SETTER_OUTPUT_LSTAR)
                    return os << "(lcp[" << t.target << "] := " << t.v1 << " +output +L*)";
                else if (t.type == LCP_QUERY)
                    return os << "(RMQ(" << t.v1 << "," << t.target << ") -> " << t.v2 << ")";
                else if (t.type == LCP_QUERY_LSTAR)
                    return os << "(RMQ(" << t.v1 << "," << t.target << ") -> " << t.v2 << " +L*)";
                else if (t.type == LCP_LS_SEAM_QUERY)
                    return os << "(L/S-seam repcount=" << t.v1 << ", match char=" << strC(alphabet_type(t.v2)) << " -> " << t.target << ")";
                else
                    return os << "(invalid LCPTuple)";
            }

        } __attribute__((packed));

        /// Order items in LCP PQ by target, which is setter-position or right-bound. Use (target DESC, type
        /// DESC) as PQ picks smallest element
        struct LCPTupleOrder
        {
            inline bool operator() (const LCPTuple& a, const LCPTuple& b) const
            {
                if (a.target == b.target)
                    return (a.type > b.type);
                else
                    return (a.target > b.target);
            }

            LCPTuple min_value() const
            {
                LCPTuple t;
                t.target = std::numeric_limits<offset_type>::max();
                t.type = std::numeric_limits<offset_type>::max();
                return t;
            }
        };

        typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<LCPTuple, LCPTupleOrder, mempartL, max_input_size/1024>::result LCP_Lpq_type;
        typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<LCPTuple, LCPTupleOrder, mempartS, max_input_size/1024>::result LCP_Spq_type;

        static const size_t lcprmq_slabsize = mempartL / (sizeof(offset_type)*2);       // size inside RMQ_Stack

#endif // ESAIS_LCP_CALC_EXT

        // ***************************************************************************************************
        // *** STXXL definitions to create L-PQ and S-PQ shared by the three steps

        typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<PQTuple, PQTupleOrder_LPQ,
                                                         mempartL, max_input_size/2/1024> Lpq_generator_type;

        typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<PQTuple, PQTupleOrder_SPQ,
                                                         mempartS, max_input_size/2/1024> Spq_generator_type;

        typedef typename Lpq_generator_type::result Lpq_type;
        typedef typename Spq_generator_type::result Spq_type;

        // input size
        offset_type     inputsize;

        // ***************************************************************************************************
        // *** STXXL typedefs to create L-Array, S-Array, S*-Array and L*-Array

        // Arrays containing continuation tuples are implemented as a runs_creator, which automatically
        // pre-sorts blocks before writing them to disk. The runs creator and mergers are encapsuled in the
        // stxxl::sorter.

        typedef stxxl::sorter< CTuple, CTupleOrder_LArray, block_size > LArray_sorter_type;
        typedef stxxl::sorter< CTuple, CTupleOrder_SArray, block_size > SArray_sorter_type;

        typedef stxxl::sorter< STuple, STupleOrder_SStarArray, block_size > SStarArray_sorter_type;
        typedef stxxl::sorter< PQTuple, PQTupleOrder_LStarArray, block_size > LStarArray_sorter_type;

        LArray_sorter_type        LArray;
        SArray_sorter_type        SArray;

        SStarArray_sorter_type    SStarArray;
        LStarArray_sorter_type*   LStarArray;   // created in induceL()

        typedef stxxl::sorter< CBufferTuple, CBufferTupleOrder_Index, block_size > CBufferTuple_sorter_type;

#if ESAIS_LCP_CALC
        typedef stxxl::sorter< offset_type, tuples::less<offset_type>, block_size > SStarLCPSkips_sorter_type;
        SStarLCPSkips_sorter_type       SStarLCPSkips;

        struct MaxRepcount
        {
            alphabet_type       charbkt;
            offset_type         maxrepcount;
            offset_type         lstar_repcount;
        };

        typedef typename stxxl::STACK_GENERATOR<MaxRepcount, stxxl::external, stxxl::grow_shrink, 4, block_size>::result MaxRepcountStack_type;

        MaxRepcountStack_type   MaxRepcountStack;

#endif // ESAIS_LCP_CALC

#if ESAIS_LCP_CALC_EXT

        typedef typename stxxl::STACK_GENERATOR<offset_type, stxxl::external, stxxl::grow_shrink, 4, block_size>::result LStarLCPStack_type;

        LStarLCPStack_type      LStarLCPStack;

#endif // ESAIS_LCP_CALC_EXT

#if ESAIS_LOG_PQFILL
        SizeLogger      Lpq_logger, Spq_logger;
        SizeLogger      LArray_logger, SArray_logger;
        SizeLogger      LStarArray_logger, SStarArray_logger;
#endif // ESAIS_LOG_PQFILL

        Result&                 m_result;

        size_type               m_mergecounter;

    public:

        Induce(Result* result)
            : LArray( CTupleOrder_LArray(), memsize / 6 ),
              SArray( CTupleOrder_SArray(), memsize / 6 ),
              SStarArray( STupleOrder_SStarArray(), memsize / 6 ),
#if ESAIS_LCP_CALC
              SStarLCPSkips( tuples::less<offset_type>(), memsize / 128 ),
#endif // ESAIS_LCP_CALC
#if ESAIS_LOG_PQFILL
              Lpq_logger("Lpq.log"), Spq_logger("Spq.log"),
              LArray_logger("LArray.log"), SArray_logger("SArray.log"),
              LStarArray_logger("LStarArray.log"), SStarArray_logger("SStarArray.log"),
#endif
              m_result(*result),
              m_mergecounter(0)
        {
        }

        // {{{ process_input(): create S* and L/S-Arrays

        // ***************************************************************************************************
        // *** Read input from back to front: iterating from S* to S* and storing continuation tuples if
        // *** needed.

        template <typename InputStreamReverse, typename SStarRankStream>
        void process_input(InputStreamReverse& inputrev, SStarRankStream& sstarrankstream)
        {
            DBGMEM("Induce process_input() 1");

            inputrev.rewind(memsize / 8);
            assert(!inputrev.empty());          // empty input

            if (inputrev.size() > std::numeric_limits<offset_type>::max())
            {
                std::cerr << "eSAIS Error: input size is larger than maximum index number.\n";
                abort();
            }

            inputsize = (size_type)inputrev.size();

            DBG(debug_induce_split, "inputsize = " << inputsize);

            // sentinel is always S* and last (real) character is always L-type, so start with a L-PQ-tuple

            offset_type spos = inputsize;

#if ESAIS_LCP_CALC_INT
            ctype_type pprev_ctype = TYPE_L; // need pre-previous items
            alphabet_type pprev_char = 0;
            size_type pprev_repcount = 0;
#endif
            ctype_type prev_ctype = TYPE_L;
            alphabet_type prev_char = 0;     // sentinel character '$', but is a special case and later ignored.
            size_type prev_repcount = 0;

            DBGMEM("Induce process_input() 2");

            while ( spos > offset_type(0) )
            {
                assert( !inputrev.empty() );
                assert( prev_ctype == TYPE_S || spos == inputsize );

#define CALC_THIS_CTYPE         ((*inputrev < prev_char || (*inputrev == prev_char && prev_ctype == TYPE_S)) ? TYPE_S : TYPE_L)

                ctype_type this_ctype = CALC_THIS_CTYPE;
                assert( this_ctype == TYPE_L || spos == inputsize );

                DBG(debug_induce_split, "input[" << spos-1 << "]: " << strC(*inputrev) << " - type " << strT(this_ctype) << " - prev_repcount " << prev_repcount);

                STuple     stuple;
                stuple.index = spos;
                stuple.chars[0] = prev_char;
                stuple.chars[1] = *inputrev;
                stuple.charfill = 2;
                stuple.continued = 0;
#if ESAIS_LCP_CALC
                stuple.repcount = prev_repcount;
#endif // ESAIS_LCP_CALC

                if (spos != inputsize)
                {
                    // Fetch rank by matching index in sorted SStarRankStream.

                    DBG(debug_induce_split, "Matching rank of S*-position " << spos << " in SStarRankStream");

                    assert( !sstarrankstream.empty() );

                    while( !sstarrankstream.empty() && sstarrankstream->index > spos ) {
                        DBG(debug, "Matching S*-position: advancing over forced split " << *sstarrankstream);
#if ESAIS_LCP_CALC
                        SStarLCPSkips.push( sstarrankstream->name );    // save forced split ranks
#endif // ESAIS_LCP_CALC
                        ++sstarrankstream;                              // advance over indexes from forced splits
                    }

                    assert( !sstarrankstream.empty() );
                    assert( sstarrankstream->index == spos );       // S*-indexes must match!

                    stuple.rank = sstarrankstream->name+1;          // names start with 0, but rank 0 is reserved for sentinel

                    DBG(debug_induce_split, "Found match to S*-position " << spos << " with rank " << stuple.rank);

                    if ( !sstarrankstream.empty() )
                        ++sstarrankstream;
                }
                else // spos == inputsize (this is the '$' S* position)
                {
                    stuple.rank = 0;
                }

                // advance to next character

#if ESAIS_LCP_CALC_INT

#define NEXTCHAR        pprev_repcount = prev_repcount; if (prev_char == *inputrev && spos != inputsize) ++prev_repcount; else prev_repcount = 0; \
                        pprev_ctype = prev_ctype, pprev_char = prev_char; \
                        prev_ctype = this_ctype, prev_char = *inputrev; \
                        --spos;

#else

#define NEXTCHAR        if (prev_char == *inputrev && spos != inputsize) ++prev_repcount; else prev_repcount = 0; \
                        prev_ctype = this_ctype, prev_char = *inputrev; \
                        --spos;

#endif
                NEXTCHAR;

                // prepare for continuation tuples

                CTuple          ctuple;
                ctype_type      ctuple_ctype = TYPE_L;  // arbitrary due to compiler warnings about uninitialized use
                ctuple.charfill = 0;

                while( spos > offset_type(0) &&
                       ( --inputrev, !((this_ctype = CALC_THIS_CTYPE) == TYPE_L && prev_ctype == TYPE_S) ) )
                {
                    assert( !inputrev.empty() );

                    DBG(debug_induce_split, "input[" << spos-1 << "]: " << strC(*inputrev) << " - type " << strT(this_ctype) << " - prev_repcount " << prev_repcount);

                    if (stuple.charfill < D)
                    {
                        // append this L-character to STuple
                        stuple.chars[stuple.charfill++] = *inputrev;

                        NEXTCHAR;
                    }
                    else
                    {
                        // STuple is full. Collect further L/S-characters into continuation tuple.

                        if (!stuple.continued)
                        {
                            // starting first L or S continuation tuple
                            stuple.continued = 1;
                        }
                        else
                        {
                            // starting next L or S continuation tuple
                            ctuple.continued = 1;

                            if (ctuple_ctype == TYPE_L) {
                                DBG(debug_induce_split, "Saving L-CTuple " << ctuple);
                                LArray.push(ctuple);
                                LOG_SIZE(LArray_logger << LArray.size());
                                LOG_WASTED(sizeof(alphabet_type) * (D - ctuple.charfill));
                            }
                            else {
                                DBG(debug_induce_split, "Saving S-CTuple " << ctuple);
                                SArray.push(ctuple);
                                LOG_SIZE(SArray_logger << SArray.size());
                                LOG_WASTED(sizeof(alphabet_type) * (D - ctuple.charfill));
                            }
                        }

#if ESAIS_LCP_CALC_INT
                        if (pprev_ctype == TYPE_L)   // continuation tuple for a L-position
#else
                        if (prev_ctype == TYPE_L)   // continuation tuple for a L-position
#endif
                        {
                            DBG(debug_induce_split, "STuple is full: " << stuple << " will be saved. Starting continuation L-CTuple.");

                            ctuple_ctype = TYPE_L;

#if ESAIS_LCP_CALC_INT
                            ctuple.index = spos+1;
                            ctuple.repcount = pprev_repcount;
                            ctuple.chars[0] = pprev_char;
                            ctuple.chars[1] = prev_char;
                            ctuple.chars[2] = *inputrev;
                            ctuple.charfill = 3;
#else
                            ctuple.index = spos;
                            ctuple.repcount = prev_repcount;
                            ctuple.chars[0] = prev_char;
                            ctuple.chars[1] = *inputrev;
                            ctuple.charfill = 2;
#endif
                            ctuple.continued = 0;


                            NEXTCHAR;

                            while( spos > offset_type(0) && ctuple.charfill < D &&
                                   ( --inputrev, !((this_ctype = CALC_THIS_CTYPE) == TYPE_L && prev_ctype == TYPE_S) ) )
                            {
                                assert( !inputrev.empty() );

                                // while this position is not an S*-position, fill up the continuation tuple

                                DBG(debug_induce_split, "input[" << spos-1 << "]: " << strC(*inputrev) << " - type " << strT(this_ctype) << " - prev_repcount " << prev_repcount);

                                // append this L/S-character to L-CTuple
                                ctuple.chars[ctuple.charfill++] = *inputrev;

                                NEXTCHAR;
                            }

                            DBG(debug_induce_split, "L-CTuple is full " << ctuple);
                        }
                        else
                        {
                            DBG(debug_induce_split, "STuple is full: " << stuple << " will be saved. Starting continuation S-CTuple.");

                            ctuple_ctype = TYPE_S;

#if ESAIS_LCP_CALC_INT
                            ctuple.index = spos+1;
                            ctuple.repcount = pprev_repcount;
                            ctuple.chars[0] = pprev_char;
                            ctuple.chars[1] = prev_char;
                            ctuple.chars[2] = *inputrev;
                            ctuple.charfill = 3;
#else
                            ctuple.index = spos;
                            ctuple.repcount = prev_repcount;
                            ctuple.chars[0] = prev_char;
                            ctuple.chars[1] = *inputrev;
                            ctuple.charfill = 2;
#endif
                            ctuple.continued = 0;

                            NEXTCHAR;

                            while( spos > offset_type(0) && ctuple.charfill < D &&
                                   ( --inputrev, !((this_ctype = CALC_THIS_CTYPE) == TYPE_L && prev_ctype == TYPE_S) ) )
                            {
                                assert( !inputrev.empty() );

                                // while this position is not an S*-position, fill up the continuation tuple

                                DBG(debug_induce_split, "input[" << spos-1 << "]: " << strC(*inputrev) << " - type " << strT(this_ctype) << " - prev_repcount " << prev_repcount);

                                // append this L/S-character to S-CTuple
                                ctuple.chars[ctuple.charfill++] = *inputrev;

                                NEXTCHAR;
                            }

                            DBG(debug_induce_split, "S-CTuple is full " << ctuple);
                        }

                        // immediately break (= end tuple) if this is S*-char before reading another char
                        if ((this_ctype = CALC_THIS_CTYPE) == TYPE_L && prev_ctype == TYPE_S)
                            break;
                    }
                }

                if (ctuple.charfill)
                {
                    if (ctuple_ctype == TYPE_L) {
                        DBG(debug_induce_split, "Saving L-CTuple " << ctuple);
                        LArray.push(ctuple);
                        LOG_SIZE(LArray_logger << LArray.size());
                        LOG_WASTED(sizeof(alphabet_type) * (D - ctuple.charfill));
                    }
                    else {
                        DBG(debug_induce_split, "Saving S-CTuple " << ctuple);
                        SArray.push(ctuple);
                        LOG_SIZE(SArray_logger << SArray.size());
                        LOG_WASTED(sizeof(alphabet_type) * (D - ctuple.charfill));
                    }
                }

                DBG(debug_induce_split, "Pushing " << stuple << " into S*-array, starting new STuple, spos = " << spos);

                SStarArray.push(stuple);
                LOG_SIZE(SStarArray_logger << SStarArray.size());
                LOG_WASTED(sizeof(alphabet_type) * (D - stuple.charfill));
            }
#undef CALC_THIS_CTYPE

            //assert( inputrev.empty() );
            assert( sstarrankstream.empty() );

            DBG_ST_ARRAY(debug_induce_arrays, "LArray", LArray);
            DBG_ST_ARRAY(debug_induce_arrays, "SArray", SArray);
            DBG_ST_ARRAY(debug_induce_arrays, "SStarArray", SStarArray);

            inputrev.finish();
            sstarrankstream.finish_clear();

            SArray.finish();
            LArray.finish();
            SStarArray.finish();
        }

        // }}}

#if ESAIS_LCP_CALC_INT
        // ***************************************************************************************************
        // *** Main Memory LCP calculation structure

        class MainMemLCP
        {
        private:

            RMQ_Stack<offset_type,offset_type> rmqstruct;

            //typedef std::map<alphabet_type,offset_type> bwtmap_type;
            //typedef typename bwtmap_type::const_iterator bwtiter_type;

            //typedef __gnu_pbds::cc_hash_table<alphabet_type,offset_type> bwtmap_type;
            //typedef typename bwtmap_type::const_point_iterator bwtiter_type;

            typedef std::tr1::unordered_map<alphabet_type,offset_type> bwtmap_type;
            typedef typename bwtmap_type::const_iterator bwtiter_type;

            bwtmap_type         prevbwt;

            alphabet_type       current_char;

        public:

            MainMemLCP()
                : current_char(0)
            {
            }

            ~MainMemLCP()
            {
                clear();
            }

            void clear()
            {
                g_mainmemlcp = std::max(g_mainmemlcp, rmqstruct.memsize() + prevbwt.size());

                rmqstruct.clear();
                prevbwt.clear();
            }

            void prepare_char(const alphabet_type& ch)
            {
                if (ch != current_char)
                {
                    clear();
                    current_char = ch;
                }
            }

            // save the BWT char of the extracted PQTuple or constructed S*/L* seed PQTuple
            void set_bwtchar(const PQTuple& t)
            {
                DBG(debug_induceL_lcp || debug_induceS_lcp, "Setting BWT at rank " << t.rank << " to bwtchar " << strC(t.chars[1]));

                prevbwt[ t.chars[1] ] = t.rank;
            }

            // save BWT and LCP without query
            void setL_noQuery(const PQTuple& t, const offset_type& thislcp)
            {
                DBG(debug_induceL_lcp, "LCP[" << t.rank << "] = " << thislcp);
                rmqstruct.set(t.rank, thislcp);
            }

            // queryL LCP RMQ: find prevbwt position and do RMQ (prevbwt_pos,this]
            void queryL(PQTuple& t, const offset_type& thislcp)
            {
                DBG(debug_induceL_lcp, "Setting L-LCP at rank " << t.rank << " to lcp " << thislcp << " and find prevbwt " << strC(t.chars[1]));

                DBG(debug_induceL_lcp, "LCP[" << t.rank << "] = " << thislcp);
                rmqstruct.set(t.rank, thislcp);

                bwtiter_type bwtiter = prevbwt.find(t.chars[1]);
                if (bwtiter == prevbwt.end())
                {
                    // no previous bwt occurance:
                    DBG(debug_induceL_lcp, "L-LCP bwtchar " << strC(t.chars[1]) << " has no previous occurance -> tenative lcp = 1");
                    t.lcp = 1;
                }
                else
                {
                    DBG1(debug_induceL_lcp, "L-LCP bwtchar " << strC(t.chars[1]) << " has previous occurance: RMQ(" << bwtiter->second+1 << " to " << t.rank << ")");
                    t.lcp = rmqstruct.query( bwtiter->second+1, t.rank ) + 1;
                    DBG3(debug_induceL_lcp, " -> result lcp = " << t.lcp);
                }

                set_bwtchar(t);
            }

            // queryLStar LCP: find RMQ to previous L*
            void queryL_lstar(PQTuple& t, const offset_type& thislcp, const offset_type& prevlstar)
            {
                DBG(debug_induceL_lcp, "L*-LCP: save lcp " << thislcp << " and find minimum to previous L*: RMQ(" << prevlstar+1 << "," << t.rank << ")");

                DBG(debug_induceL_lcp, "LCP[" << t.rank << "] = " << thislcp);
                rmqstruct.set(t.rank, thislcp);

                if (prevlstar != std::numeric_limits<offset_type>::max())
                {
                    t.lcp = rmqstruct.query(prevlstar+1, t.rank);
                    DBG(debug_induceL_lcp, "RMQ(" << prevlstar+1 << "," << t.rank << ") -> " << t.lcp);
                }
                else
                {
                    t.lcp = 0;  // no previous L* in this bucket
                    DBG(debug_induceL_lcp, "RMQ(" << prevlstar+1 << "," << t.rank << ") -> invalid = 0");
                }
            }

            // save BWT and LCP without query
            void setS_noQuery(const PQTuple& t, const offset_type& thislcp)
            {
                DBG(debug_induceS_lcp, "LCP[" << t.rank-1 << "] = " << thislcp);
                rmqstruct.set(t.rank-1, thislcp);
            }

            // queryS LCP RMQ for L*s: find prevbwt positions, do RMQ [prevbwt_pos,this) and then set L*-LCP
            void queryS_lstar(PQTuple& t, offset_type thislcp)  // no const& here, as thislcp must be copied.
            {
                DBG(debug_induceS_lcp, "S-LCP query at rank " << t.rank << " with bwtchar " << strC(t.chars[1]) << " for L*-LCP");

                bwtiter_type bwtiter = prevbwt.find(t.chars[1]);
                if (bwtiter == prevbwt.end())
                {
                    // no previous bwt occurance:
                    DBG(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has no previous occurance -> tenative lcp = 1");
                    t.lcp = 1;
                }
                else
                {
                    DBG1(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has previous occurance: RMQ(" << bwtiter->second << " to " << t.rank-1 << ")");
                    t.lcp = rmqstruct.query( bwtiter->second , t.rank - 1 ) + 1;
                    DBG3(debug_induceS_lcp, " -> result lcp = " << t.lcp);
                }

                DBG(debug_induceS_lcp, "Setting S-LCP at rank " << t.rank << " to lcp " << thislcp << " after prevbwt find.");

                DBG(debug_induceS_lcp, "LCP[" << t.rank << "] = " << thislcp);
                rmqstruct.set(t.rank, thislcp);

                set_bwtchar(t);
            }

            // queryS LCP RMQ without lcp labeling (for first positions in repbucket)
            void queryS_nolcp(PQTuple& t)
            {
                DBG(debug_induceS_lcp, "Query S-LCP at rank " << t.rank << " find prevbwt " << strC(t.chars[1]));

                bwtiter_type bwtiter = prevbwt.find(t.chars[1]);
                if (bwtiter == prevbwt.end())
                {
                    // no previous bwt occurance:
                    DBG(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has no previous occurance -> tenative lcp = 1");
                    t.lcp = 1;
                }
                else
                {
                    DBG1(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has previous occurance: RMQ(" << bwtiter->second << " to " << t.rank-1 << ")");
                    t.lcp = rmqstruct.query( bwtiter->second , t.rank - 1 ) + 1;
                    DBG3(debug_induceS_lcp, " -> result lcp = " << t.lcp);
                }

                set_bwtchar(t);
            }

            // queryS LCP RMQ: first set rank-1's LCP, find prevbwt positions, do RMQ [prevbwt_pos,this)
            void queryS(PQTuple& t, const offset_type& thislcp)
            {
                DBG(debug_induceS_lcp, "Setting S-LCP at rank " << t.rank << " to lcp " << thislcp << " and find prevbwt " << strC(t.chars[1]));

                DBG(debug_induceS_lcp, "LCP[" << t.rank-1 << "] = " << thislcp);
                rmqstruct.set(t.rank-1, thislcp);

                bwtiter_type bwtiter = prevbwt.find(t.chars[1]);
                if (bwtiter == prevbwt.end())
                {
                    // no previous bwt occurance:
                    DBG(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has no previous occurance -> tenative lcp = 1");
                    t.lcp = 1;
                }
                else
                {
                    DBG1(debug_induceS_lcp, "S-LCP bwtchar " << strC(t.chars[1]) << " has previous occurance: RMQ(" << bwtiter->second << " to " << t.rank-1 << ")");
                    t.lcp = rmqstruct.query( bwtiter->second , t.rank - 1 ) + 1;
                    DBG3(debug_induceS_lcp, " -> result lcp = " << t.lcp);
                }

                set_bwtchar(t);
            }

            // setRepbucketEnd LCP RMQ: set first position of an ending repbucket with repcount or L/S-seam result
            void setRepbucketEnd(const offset_type& rank, const offset_type& thislcp)
            {
                DBG(debug_induceS_lcp, "Setting S-LCP repbucket boundary at rank " << rank << " to lcp " << thislcp);

                DBG(debug_induceS_lcp, "LCP[" << rank << "] = " << thislcp);
                rmqstruct.set(rank, thislcp);
            }
        };

#endif // ESAIS_LCP_CALC_INT

#if ESAIS_LCP_CALC_EXT

        template <typename PQType>
        static inline void insertSplitRMQ(int type, const offset_type& left, const offset_type& right, const offset_type& target, PQType& pq)
        {
            if ((left / lcprmq_slabsize) == (right / lcprmq_slabsize)) {
                //std::cout << "Insert unsplit RMQ(" << left << "," << right << ")\n";
                LCPTuple l = { type, right, left, target };
                DBG(debug_induceL, "Query LCP-tuple: " << l);
                pq->push(l);
            }
            else {
                // calculate lower boundary of next slab after left
                size_type splitpos = (((size_type)(left / lcprmq_slabsize)+1) * lcprmq_slabsize);

                //std::cout << "Insert split RMQ(" << left << "," << right << ")"
                //          << " = RMQ(" << left << "," << splitpos-1 << ") & RMQ(" << splitpos << "," << right << ")\n";

                LCPTuple l1 = { type, splitpos-1, left, target };
                LCPTuple l2 = { type, right, splitpos, target };
                DBG(debug_induceL, "Query LCP-tuples: " << l1 << " and " << l2);
                pq->push(l1), pq->push(l2);
            }
        }
#endif // ESAIS_LCP_CALC_EXT

        // ***************************************************************************************************
        // *** Induce Ls

        void reinsertLTuple(PQTuple& t, Lpq_type& Lpq)
        {
            if ( t.chars[1] >= t.chars[0] )
            {
                // previous position is L-type, decrease and reinsert

                t.decrease();
                Lpq.push(t);

                DBG(debug_induceL, "reinserted " << t << " into L-PQ.");

                LOG_SIZE(Lpq_logger << Lpq.size());
                LOG_WASTED(sizeof(alphabet_type) * (D - t.charfill));
            }
            else // if (t.chars[1] < t.chars[0])
            {
                // previous position is S-type, so this is L*-type

                LStarArray->push(t);

                DBG(debug_induceL, "saved " << t << " into L*-Array.");

                LOG_SIZE(LStarArray_logger << LStarArray->size());
                LOG_WASTED(sizeof(alphabet_type) * (D - t.charfill));
            }
        }

        // induceL() input: S*-Array, L-Array, (S*-LCPs). work: L-PQ, PQ-pool, (LCP-PQ, PQ-pool), MergeBuffer. output: L*-Array.

        template <typename SStarLCPStream>
        void induceL(SStarLCPStream& sstarlcpstream)
        {
#if !ESAIS_LCP_CALC
            stxxl::STXXL_UNUSED(sstarlcpstream);
#endif // ESAIS_LCP_CALC

            DBG(debug_induceL, "=== start induceL() ============================================================");
            DBGMEM("induce initializing");

            offset_type relRank = 0;                        // rank for Ls (counts only Ls and S*s)

            // the L-PQ for inducing, allocated on heap due to it containing large buffer fields
            std::auto_ptr<Lpq_type> Lpq( new Lpq_type(mempartL / 2, mempartL / 2) );

            DBG(debug, "L-PQ parameters:"
                << " total_memory=" << Lpq->mem_cons()
                << " delete_buffer_size=" << Lpq->delete_buffer_size
                << " N=" << Lpq->N
                << " IntKMAX=" << Lpq->IntKMAX
                << " num_int_groups=" << Lpq->num_int_groups
                << " num_ext_groups=" << Lpq->num_ext_groups
                << " total_num_groups=" << Lpq->total_num_groups
                << " BlockSize=" << Lpq->BlockSize
                << " ExtKMAX=" << Lpq->ExtKMAX );

            // create L*-Array
            LStarArray = new LStarArray_sorter_type(PQTupleOrder_LStarArray(), mempartL);

#if ESAIS_LCP_CALC_EXT
            // the LCP L-PQ for LCP calculation, allocated on heap due to it containing large buffer fields
            std::auto_ptr<LCP_Lpq_type> LCPpq( new LCP_Lpq_type(mempartL / 4, mempartL / 4) );

            DBG(debug, "LCP-L-PQ parameters:"
                << " total_memory=" << LCPpq->mem_cons()
                << " delete_buffer_size=" << LCPpq->delete_buffer_size
                << " N=" << LCPpq->N
                << " IntKMAX=" << LCPpq->IntKMAX
                << " num_int_groups=" << LCPpq->num_int_groups
                << " num_ext_groups=" << LCPpq->num_ext_groups
                << " total_num_groups=" << LCPpq->total_num_groups
                << " BlockSize=" << LCPpq->BlockSize
                << " ExtKMAX=" << LCPpq->ExtKMAX );
#endif

#if ESAIS_LCP_CALC
            // initialize S*-LCP information
            sstarlcpstream.sort(mempartL / 2);

            SStarLCPSkips.sort(mempartL / 128);
            //DBG_ST_ARRAY(1, "SStarLCPSkips", SStarLCPSkips);

#endif // ESAIS_LCP_CALC

            // S*-Array input stream which contains the seed tuples
            SStarArray.sort(mempartL);

            // L-Array stream running through all continuation tuples in index order
            LArray.sort(mempartL);

            // buffer for tuples in need of continuation
            CBufferTuple_sorter_type    LMergeBuffer( CBufferTupleOrder_Index(), mempartL );

            size_type repcount = 0;
            alphabet_type prevCharLimit = std::numeric_limits<alphabet_type>::max();

#if ESAIS_LCP_CALC
            size_type lstar_repcount = 0;
#endif // ESAIS_LCP_CALC

#if ESAIS_LCP_CALC_INT
            MainMemLCP mmlcp;
            offset_type prev_lstar = std::numeric_limits<offset_type>::max();
#endif

            DBGMEM("induce starting");

            while( !Lpq->empty() ||                 // tuples left in PQ
                   !SStarArray.empty() )            // S*-positions left as seeds
            {
                // *** Seed from new S* positions into PQ if possible.

#if ESAIS_LCP_CALC
                bool firstsstar = true;
#endif // ESAIS_LCP_CALC

                while ( !SStarArray.empty() &&
                        ( Lpq->empty() || SStarArray->chars[0] < Lpq->top().chars[0] ) )
                {
                    PQTuple t = PQTuple::fromSTuple(*SStarArray);
#if ESAIS_LCP_CALC
                    offset_type sstarrepcount = SStarArray->repcount;
                    alphabet_type seedchar = t.chars[0];
#endif // ESAIS_LCP_CALC
                    ++SStarArray;

                    LOG_SIZE(SStarArray_logger << SStarArray.size());

                    DBG(debug_induceL, "Seeding from S*-tuple " << t);
                    assert( t.chars[1] > t.chars[0] || t.index == inputsize );

#if ESAIS_LCP_CALC_INT
                    mmlcp.prepare_char(t.chars[0]);
#endif

                    t.rank = relRank; ++relRank;

#if ESAIS_LCP_CALC
                    offset_type t_lcp;

                    if ( t.rank == offset_type(0) ) // first position has sentinel lcp 0
                    {
                        t_lcp = 0;
                        firstsstar = false;
                    }
                    else if (firstsstar) // first S* in a bucket must be compared to preceding L repcount
                    {
                        if (prevCharLimit == seedchar)
                        {
                            DBG(debug_induceL, "first S* in bucket: current L repcount = " << repcount << " S*-repcount " << sstarrepcount);
                            t_lcp = std::min<offset_type>(repcount, sstarrepcount) + 1;
                        }
                        else {
                            DBG(debug_induceL, "first S* in bucket: no preceding Ls -> 0");
                            t_lcp = 0;
                        }

                        while ( !sstarlcpstream.empty() && !SStarLCPSkips.empty() &&
                                sstarlcpstream.curr_index() == *SStarLCPSkips )
                        {
                            DBG(debug_sstarlcp_use, "Skipped S*-LCP[" << sstarlcpstream.curr_index() << "]: " << *sstarlcpstream);
                            ++sstarlcpstream, ++SStarLCPSkips;
                        }

                        assert( !sstarlcpstream.empty() );
                        DBG(debug_sstarlcp_use, "S*-LCP[" << sstarlcpstream.curr_index() << "]: " << *sstarlcpstream << " unused");
                        ++sstarlcpstream;       // skip unused S*-lcp

                        firstsstar = false;
                    }
                    else
                    {
                        assert( !sstarlcpstream.empty() );
                        t_lcp = *sstarlcpstream;

                        // if this position is not a real S*-index, but a forced-split, then calculate dumb
                        // RMQ to the next real S*.
                        while ( !sstarlcpstream.empty() && !SStarLCPSkips.empty() &&
                                sstarlcpstream.curr_index() == *SStarLCPSkips )
                        {
                            DBG(debug_sstarlcp_use, "Skipped S*-LCP[" << sstarlcpstream.curr_index() << "]: " << *sstarlcpstream );
                            ++sstarlcpstream, ++SStarLCPSkips;

                            assert( !sstarlcpstream.empty() );
                            t_lcp = std::min(t_lcp, *sstarlcpstream);
                        }

                        DBG(debug_sstarlcp_use, "S*-LCP[" << sstarlcpstream.curr_index() << "]: " << *sstarlcpstream);

                        ++sstarlcpstream;
                        DBG(debug_induceL, "seeding S*-LCP from LCP stream: " << t_lcp);
                    }

#if ESAIS_LCP_CALC_INT

                    mmlcp.queryL(t, t_lcp);

#elif ESAIS_LCP_CALC_EXT

                    LCPTuple l = { LCP_SETTER, t.rank, t_lcp, 0 };
                    DBG(debug_induceL, "Setter LCP-tuple: " << l);
                    LCPpq->push(l);
#endif

#endif // ESAIS_LCP_CALC

                    t.decrease();

                    Lpq->push(t);
                    LOG_SIZE(Lpq_logger << Lpq->size());
                    LOG_WASTED(sizeof(alphabet_type) * (D - t.charfill));

                    DBG(debug_induceL, "inserted L-tuple " << t << " into L-PQ");
                }
                // determine current limit char

                alphabet_type charLimit = Lpq->top().chars[0];
                offset_type rankLimit = relRank;

#if ESAIS_LCP_CALC
                if (charLimit != prevCharLimit)
                {
                    repcount = lstar_repcount = 0;
                    ESAIS_LCP_CALCX_INT( prev_lstar = std::numeric_limits<offset_type>::max(); )
                }
                else
                {
                    ++repcount;
                }

#if ESAIS_LCP_CALC_INT
                mmlcp.prepare_char(charLimit);
#endif
                offset_type prevSrcRank = std::numeric_limits<offset_type>::max();

#else // !ESAIS_LCP_CALC

                repcount = (charLimit == prevCharLimit) ? repcount+1 : 0;

#endif // ESAIS_LCP_CALC

                DBG(debug_induceL, "----------------------------------------------------------------------------------------------------");
                DBG(debug_induceL, "Extracting from L-PQ with repcount = " << repcount << ", charLimit = " << strC(charLimit) << " and rankLimit = " << rankLimit);

                while ( !Lpq->empty() &&
                        Lpq->top().chars[0] == charLimit && Lpq->top().rank < rankLimit )
                {
                    PQTuple t = Lpq->top(); Lpq->pop();

                    LOG_SIZE(Lpq_logger << Lpq->size());

                    DBG(debug_induceL, "Processing L-tuple " << t << " given index " << t.index << " relRank " << relRank);
                    DBG(debug_induceL, "--> SA " << t.index << " is next L-entry (relRank = " << relRank << ")");

                    m_result.output_Lentry(t.chars[0], t.index);

                    ESAIS_LCP_CALCX( offset_type srcRank = t.rank; )

                    t.rank = relRank; ++relRank;

#if ESAIS_LCP_CALC
                    bool isLStar = (t.charfill > 1 && t.chars[1] < t.chars[0]);

#if ESAIS_LCP_CALC_EXT
                    if (prevSrcRank == std::numeric_limits<offset_type>::max())
                    {
                        DBG(debug_induceL, "First entry in repbucket -> LCP = repcount = " << repcount);

                        LCPTuple l = { LCP_SETTER_OUTPUT + isLStar, t.rank, repcount, 0 };
                        DBG(debug_induceL, "Setter LCP-tuple: " << l);
                        LCPpq->push(l);
                    }
                    else
                    {
                        DBG(debug_induceL, "Calculate LCP between prevSrcRank = " << prevSrcRank << "+1 and " << srcRank << " relTarget = " << t.rank);
                        insertSplitRMQ(LCP_QUERY + isLStar, prevSrcRank+1, srcRank, t.rank, LCPpq);
                    }
#endif
#if ESAIS_LCP_CALC_INT
                    offset_type t_lcp = t.lcp;

                    if (prevSrcRank == std::numeric_limits<offset_type>::max())
                    {
                        DBG(debug_induceL_lcp, "L-LCP: first entry in repbucket -> override LCP = repcount = " << repcount);
                        t_lcp = repcount;
                    }

                    DBG(debug_induceL, "--> LCP " << t_lcp);
                    m_result.output_Lentry_lcp( t_lcp );

                    // if not L* -> precalculate LCP of induced L-type
                    if (!isLStar)
                    {
                        if (t.charfill > 1) {
                            mmlcp.queryL(t, t_lcp);
                        }
                        else {
                            mmlcp.setL_noQuery(t, t_lcp);
                        }
                    }
                    else { // for L*s (induced char is S-type) -> calculate LCP to the preceding L*
                        mmlcp.queryL_lstar(t, t_lcp, prev_lstar);
                        prev_lstar = t.rank;
                    }
#endif // ESAIS_LCP_CALC_INT

                    if (isLStar) lstar_repcount = repcount+1;
                    prevSrcRank = srcRank;

#endif // ESAIS_LCP_CALC

                    // decrease varlength tuple if possible and reinsert

#if !ESAIS_LCP_CALC_INT
                    if (t.charfill == 1)
                    {
                        if ( t.index == offset_type(0) )
                        {
                            DBG(debug_induceL, "First position processed and finished, no reinsert.");
                        }
                        else if (t.continued)
                        {
                            DBG(debug_induceL, "Process: charfill=1, put " << t << " into L-MergeBuffer.");

                            ++m_mergecounter;
                            LMergeBuffer.push( CBufferTuple::fromPQTuple(t) );
                        }
                    }
#else // ESAIS_LCP_CALC_INT
                    if (t.charfill == 1)
                    {
                        if ( t.index == offset_type(0) ) {
                            DBG(debug_induceL, "First position processed and finished, no reinsert.");
                        }
                        assert(!t.continued);
                    }
                    else if (t.charfill == 2 && t.continued)
                    {
                        DBG(debug_induceL, "Process: charfill=2, put " << t << " into L-MergeBuffer.");

                        ++m_mergecounter;
                        LMergeBuffer.push( CBufferTuple::fromPQTuple(t) );
                    }
#endif // ESAIS_LCP_CALC
                    else
                    {
                        reinsertLTuple(t, *Lpq);
                    }
                }

                if ( LMergeBuffer.size() )      // Some tuples need continuation
                {
                    // PQ has advanced as far as possible without merging in tuples that need more characters.

                    DBG(debug_induceL, "Must remerge tuples needing more characters.");

                    LMergeBuffer.sort_reuse();

                    DBG_ST_ARRAY(debug_induceL, "LMergeBuffer", LMergeBuffer);

#if ESAIS_SELF_CHECK && !defined(NDEBUG)
                    { // check that all starting chars of the merge buffer are the same

                        alphabet_type mergechar = LMergeBuffer->char0;
                        while ( !LMergeBuffer.empty() )
                        {
                            assert( LMergeBuffer->char0 == mergechar );
                            ++LMergeBuffer;
                        }
                        LMergeBuffer.rewind();
                    }
#endif
                    // merge array with tuples from L-Array

                    assert( !LArray.empty() );

                    while ( !LMergeBuffer.empty() )
                    {
                        assert( !LArray.empty() );
                        assert( LArray->index == LMergeBuffer->index );
                        assert( LArray->repcount == offset_type(repcount) );

                        // construct new PQ tuple from continuation tuple
                        PQTuple t = PQTuple::fromCTuple( *LArray );

                        t.rank = LMergeBuffer->rank;

#if ESAIS_LCP_CALC_INT
                        t.lcp = LMergeBuffer->lcp;
#endif
#if ESAIS_LCP_CALC_EXT
                        if (t.chars[1] < t.chars[0]) {
                            // this merged tuples is an L* tuple, so we missed the L*-star flag on the RMQ-tuple

                            LCPTuple l = { LCP_LSTAR_MARKER, LMergeBuffer->rank, 0, 0 };
                            DBG(debug_induceL, "Adding LCP-tuple marker tuple for missed L*: " << l);
                            LCPpq->push(l);

                            lstar_repcount = repcount+1;
                        }
#endif // ESAIS_LCP_CALC_EXT

                        DBG(debug_induceL, "Merged L-tuple " << t);

                        reinsertLTuple(t, *Lpq);

                        if ( !LArray.empty() ) {
                            ++LArray;
                            LOG_SIZE(LArray_logger << LArray.size());
                        }

                        ++LMergeBuffer;
                    }

                    LMergeBuffer.clear();
                }

#if ESAIS_LCP_CALC
                // *** Save repcount for L/S-seam calculation

                if ( Lpq->empty() || charLimit != Lpq->top().chars[0])
                {
                    DBG(debug_induceL, "RepcountStack: pushing repcount=" << repcount << " and L*-repcount=" << lstar_repcount << " for char " << strC(charLimit));
                    MaxRepcount mr = { charLimit, repcount+1, lstar_repcount };
                    MaxRepcountStack.push(mr);
                }
#endif

                prevCharLimit = charLimit;

                if (!Lpq->empty())
                    DBG(debug_induceL, "L-PQ.top after break: " << Lpq->top());
            }

            DBGMEM("induce-only finished");

            SStarArray.finish_clear();
            LArray.finish_clear();
            LMergeBuffer.finish_clear();

#if ESAIS_LCP_CALC_EXT
            DBGMEM("induce answering RMQs");

            DBG(debug_induceL, "====================================================================================================");

            RMQ_Stack_Blocked<offset_type,offset_type,lcprmq_slabsize> rmqstruct;

            bool lstarmarker = false;
            offset_type lstar_rangemin = 0;

            while ( !LCPpq->empty() )
            {
                LCPTuple t = LCPpq->top(); LCPpq->pop();

                if (t.target % 1000 == 0) {
                    //rmqstruct.print_size();
                }

                DBG(debug_induceL_lcp, "Extracted from LCP-PQ: " << t);

                if (t.type == LCP_LSTAR_MARKER)
                {
                    lstarmarker = true;
                }
                else if (t.type == LCP_SETTER || t.type == LCP_SETTER_OUTPUT || t.type == LCP_SETTER_OUTPUT_LSTAR)
                {
                    // collect more answers
                    while ( !LCPpq->empty() && LCPpq->top().target == t.target && LCPpq->top().type == t.type ) {
                        DBG(debug_induceL_lcp, "Min-merged with LCP-PQ: " << LCPpq->top());
                        t.v1 = std::min(t.v1, LCPpq->top().v1);
                        LCPpq->pop();
                    }

                    rmqstruct.set(t.target, t.v1);

                    lstar_rangemin = std::min(lstar_rangemin, t.v1);    // manual RMQ to last L*

                    if (t.type == LCP_SETTER_OUTPUT || t.type == LCP_SETTER_OUTPUT_LSTAR)
                    {
                        m_result.output_Lentry_lcp( t.v1 );
                    }
                    if (t.type == LCP_SETTER_OUTPUT_LSTAR || lstarmarker) {
                        // don't push the value itself, but the RMQ to the previous L*

                        DBG(debug_induceL_lcp, "Pushing into L*-LCP Stack: " << lstar_rangemin);
                        LStarLCPStack.push(lstar_rangemin);

                        lstar_rangemin = std::numeric_limits<offset_type>::max();
                        lstarmarker = false;
                    }
                }
                else if (t.type == LCP_QUERY || t.type == LCP_QUERY_LSTAR)
                {
                    const offset_type &left = t.v1, &right = t.target;
                    assert( left <= right );

                    const offset_type& lcp = rmqstruct.query(left, right) + 1;

                    LCPTuple l = { LCP_SETTER_OUTPUT + (t.type == LCP_QUERY_LSTAR ? 1 : 0), t.v2, lcp, 0 };
                    DBG(debug_induceL_lcp, "Reinserting answer LCP-tuple: " << l);
                    LCPpq->push(l);
                }
                else {
                    assert(!"Unknown RMQ-Tuple");
                }
            }

            rmqstruct.print_size();

#endif // ESAIS_LCP_CALC_EXT

            m_result.finish_Lsequence();

            DBGMEM("induce finished");
        }

        // ***************************************************************************************************
        // *** Induce Ss

        void reinsertSTuple(PQTuple& t, Spq_type& Spq)
        {
            if ( t.chars[1] <= t.chars[0] )
            {
                // previous position is S-type

                t.decrease();

                DBG(debug_induceS, "reinserted " << t << " into S-PQ.");

                Spq.push(t);
                LOG_SIZE(Spq_logger << Spq.size());
                LOG_WASTED(sizeof(alphabet_type) * (D - t.charfill));
            }
            else // if (t.chars[1] > t.chars[0])
            {
                // previous position is L-type, so this is S*-type
                assert(!"Impossible reinsertion of L-type as tuples should not contain L-types anymore!");

                DBG(debug_induceS, "not reinserted, finished.");
            }
        }

        // induceS() input: L*-Array, S-Array. work: S-PQ, PQ-pool, (LCP-PQ, PQ-pool), MergeBuffer. output: (only stack+queue).

        void induceS()
        {
            DBG(debug_induceS, "=== start induceS() ============================================================");
            DBGMEM("induce initializing");

            LStarArray_sorter_type& LStarArray = *this->LStarArray;

            DBG_ST_ARRAY(debug_induce_arrays, "LStarArray", LStarArray);

            offset_type relRank = 0;                // rank for Ss (in reverse order, counting only Ss)

            // the S-PQ for inducing, allocated on heap due to it containing large buffer fields
            std::auto_ptr<Spq_type> Spq( new Spq_type(mempartS / 2, mempartS / 2) );

            DBG(debug, "S-PQ parameters:"
                << " total_memory=" << Spq->mem_cons()
                << " delete_buffer_size=" << Spq->delete_buffer_size
                << " N=" << Spq->N
                << " IntKMAX=" << Spq->IntKMAX
                << " num_int_groups=" << Spq->num_int_groups
                << " num_ext_groups=" << Spq->num_ext_groups
                << " total_num_groups=" << Spq->total_num_groups
                << " BlockSize=" << Spq->BlockSize
                << " ExtKMAX=" << Spq->ExtKMAX );

            // L*-Array input stream which contains the seed tuples
            LStarArray.sort(mempartS);

#if ESAIS_LCP_CALC_EXT
            // the LCP S-PQ for LCP calculation, allocated on heap due to it containing large buffer fields
            std::auto_ptr<LCP_Spq_type> LCPpq( new LCP_Spq_type(mempartS / 2, mempartS / 2) );

            DBG(debug_induceS, "L*-LCPStack size " << LStarLCPStack.size() << " vs. L*-Array size " << LStarArray.size());
            assert((size_type)LStarLCPStack.size() == LStarArray.size());

            DBG(debug, "LCP-S-PQ parameters:"
                << " total_memory=" << LCPpq->mem_cons()
                << " delete_buffer_size=" << LCPpq->delete_buffer_size
                << " N=" << LCPpq->N
                << " IntKMAX=" << LCPpq->IntKMAX
                << " num_int_groups=" << LCPpq->num_int_groups
                << " num_ext_groups=" << LCPpq->num_ext_groups
                << " total_num_groups=" << LCPpq->total_num_groups
                << " BlockSize=" << LCPpq->BlockSize
                << " ExtKMAX=" << LCPpq->ExtKMAX );
#endif

#if ESAIS_LCP_CALC_INT
            MainMemLCP mmlcp;
#endif // ESAIS_LCP_CALC

            // buffer for tuples in need of continuation
            CBufferTuple_sorter_type    SMergeBuffer( CBufferTupleOrder_Index(), mempartS );

            // iterator of S-Array running through all continuation tuples in index order
            SArray.sort(mempartS);

            size_type repcount = 0;
            alphabet_type prevCharLimit = std::numeric_limits<alphabet_type>::min();

            DBGMEM("induce starting");

            while( !Spq->empty() ||                  // tuples left in PQ
                   !LStarArray.empty() )             // L*-positions left as seeds
            {
                // *** Seed from new L* positions into PQ if possible.

                while ( !LStarArray.empty() &&
                        ( Spq->empty() || LStarArray->chars[0] > Spq->top().chars[0] ) )
                {
                    PQTuple t = *LStarArray;
#if ESAIS_LCP_CALC_EXT
                    offset_type t_lcp = LStarLCPStack.top(); LStarLCPStack.pop();
#endif // ESAIS_LCP_CALC_EXT
                    ++LStarArray;
                    LOG_SIZE(LStarArray_logger << LStarArray.size());

                    DBG(debug_induceS, "Seeding from L*-tuple " << t);
                    assert( t.chars[1] < t.chars[0] );

                    t.rank = relRank; ++relRank;

#if ESAIS_LCP_CALC_INT

                    mmlcp.prepare_char(t.chars[0]);
                    mmlcp.queryS_lstar(t, t.lcp);

#endif
#if ESAIS_LCP_CALC_EXT

                    LCPTuple l = { LCP_SETTER, t.rank, t_lcp, 0 };
                    DBG(debug_induceS, "Setter LCP-tuple: " << l);
                    LCPpq->push(l);

#endif // ESAIS_LCP_CALC_EXT

                    t.decrease();

                    Spq->push(t);
                    LOG_SIZE(Spq_logger << Spq->size());
                    LOG_WASTED(sizeof(alphabet_type) * (D - t.charfill));

                    DBG(debug_induceS, "inserted S-tuple " << t << " into S-PQ");
                }

                // determine current limit char

                alphabet_type charLimit = Spq->top().chars[0];
                offset_type rankLimit = relRank;

#if ESAIS_LCP_CALC_INT
                mmlcp.prepare_char(charLimit);
#endif
#if ESAIS_LCP_CALC
                offset_type prevSrcRank = std::numeric_limits<offset_type>::max();
#endif

                repcount = (charLimit == prevCharLimit) ? repcount+1 : 0;

                DBG(debug_induceS, "----------------------------------------------------------------------------------------------------");
                DBG(debug_induceS, "Extracting from S-PQ with charLimit = " << strC(charLimit) << " and rankLimit = " << rankLimit);

                while ( !Spq->empty() &&
                        Spq->top().chars[0] == charLimit && Spq->top().rank < rankLimit)
                {
                    PQTuple t = Spq->top(); Spq->pop();

                    LOG_SIZE(Spq_logger << Spq->size());

                    DBG(debug_induceS, "Processing S-tuple " << t << " given index " << t.index << " relRank " << relRank);
                    DBG(debug_induceS, "--> SA " << t.index << " is next S-entry (relRank = " << relRank << ")");

                    m_result.output_Sentry(t.chars[0], t.index);

                    ESAIS_LCP_CALCX( offset_type srcRank = t.rank; )

                    t.rank = relRank; ++relRank;

#if ESAIS_LCP_CALC
                    if (prevSrcRank == std::numeric_limits<offset_type>::max())
                    {
                        // do nothing for very first entry from this repchar-sequence
#if ESAIS_LCP_CALC_INT
                        DBG(debug_induceS_lcp, "LCP: delay calculation for first entry of repbucket");
                        mmlcp.queryS_nolcp(t);
#endif
                    }
                    else
                    {
#if ESAIS_LCP_CALC_INT
                        DBG(debug_induceS, "--> LCP " << t.lcp);
                        m_result.output_Sentry_lcp( t.lcp );

                        if (t.charfill > 1) {
                            mmlcp.queryS(t, t.lcp);
                        }
                        else {
                            mmlcp.setS_noQuery(t, t.lcp);
                        }

#elif ESAIS_LCP_CALC_EXT
                        DBG(debug_induceS, "Calculate LCP between prevSrcRank = " << prevSrcRank << " and " << srcRank << "-1 relTarget = " << t.rank-1);
                        insertSplitRMQ(LCP_QUERY, prevSrcRank, srcRank-1, t.rank-1, LCPpq);
#endif
                    }

                    prevSrcRank = srcRank;
#endif // ESAIS_LCP_CALC

                    // decrease varlength tuple if possible and reinsert
#if !ESAIS_LCP_CALC_INT
                    if (t.charfill == 1)
                    {
                        if ( t.index == offset_type(0) )
                        {
                            DBG(debug_induceS, "first position finished. no reinsert.");
                        }
                        else if (t.continued)
                        {
                            DBG(debug_induceS, "charfill=1, put " << t << " into S-MergeBuffer");

                            ++m_mergecounter;
                            SMergeBuffer.push( CBufferTuple::fromPQTuple(t) );
                        }
                        else
                        {
                            DBG(debug_induceS, "no continuation, finished.");
                        }
                    }
#else // ESAIS_LCP_CALC_INT
                    if (t.charfill == 1)
                    {
                        if ( t.index == offset_type(0) ) {
                            DBG(debug_induceS, "first position finished. no reinsert.");
                        }
                        assert(!t.continued);
                    }
                    else if (t.charfill == 2 && t.continued)
                    {
                        DBG(debug_induceS, "charfill=2, put " << t << " into S-MergeBuffer");

                        ++m_mergecounter;
                        SMergeBuffer.push( CBufferTuple::fromPQTuple(t) );
                    }
#endif // ESAIS_LCP_CALC_EXT
                    else
                    {
                        reinsertSTuple(t, *Spq);
                    }
                }

                if ( SMergeBuffer.size() )      // Some tuples need continuation
                {
                    // PQ has advanced as far as possible without merging in
                    // tuples that need more characters.

                    DBG(debug_induceS, "Must remerge tuples needing more characters.");

                    SMergeBuffer.sort_reuse();

                    DBG_ST_ARRAY(debug_induceS, "SMergeBuffer", SMergeBuffer);

#if ESAIS_SELF_CHECK && !defined(NDEBUG)
                    { // check that all starting chars of the merge buffer are the same

                        alphabet_type mergechar = SMergeBuffer->char0;
                        while ( !SMergeBuffer.empty() )
                        {
                            assert( SMergeBuffer->char0 == mergechar );
                            ++SMergeBuffer;
                        }
                        SMergeBuffer.rewind();
                    }
#endif
                    // merge array with tuples from S-Array

                    assert( !SArray.empty() );

                    while ( !SMergeBuffer.empty() )
                    {
                        assert( !SArray.empty() );
                        assert( SArray->index == SMergeBuffer->index );
                        assert( SArray->repcount == offset_type(repcount) );

                        // construct new PQ tuple from continuation tuple
                        PQTuple t = PQTuple::fromCTuple( *SArray );

                        t.rank = SMergeBuffer->rank;
#if ESAIS_LCP_CALC_INT
                        t.lcp = SMergeBuffer->lcp;
#endif
                        DBG(debug_induceS, "Merged S-tuple " << t);

                        reinsertSTuple(t, *Spq);

                        if ( !SArray.empty() ) {
                            ++SArray;
                            LOG_SIZE(SArray_logger << SArray.size());
                        }

                        ++SMergeBuffer;
                    }

                    SMergeBuffer.clear();
                }

#if ESAIS_LCP_CALC_EXT
                // handle LCP for preceding position: either push repetition bucket number or L/S-seam tuples
                if (!Spq->empty() && charLimit == Spq->top().chars[0])
                {
                    // next position is still within the same bucket

                    assert(prevSrcRank != std::numeric_limits<offset_type>::max());

                    LCPTuple l = { LCP_SETTER_OUTPUT, relRank-1, repcount+1, 0 };
                    DBG(debug_induceS, "Pushed repbucket LCP-tuple: " << l << " with repcount " << repcount+1);
                    LCPpq->push(l);
                }
                else
                {
                    // create L/S-seam query tuple
                    LCPTuple l = { LCP_LS_SEAM_QUERY, relRank-1, repcount+1, offset_type(charLimit) };
                    DBG(debug_induceS, "Pushed seam LCP-tuple: " << l);
                    LCPpq->push(l);
                }
#endif // ESAIS_LCP_CALC_EXT

#if ESAIS_LCP_CALC_INT
                if (!Spq->empty() && charLimit == Spq->top().chars[0])
                {
                    // next position is still within the same bucket

                    assert(prevSrcRank != std::numeric_limits<offset_type>::max());

                    DBG(debug_induceS, "--> LCP " << repcount+1 << " (= repcount at repbucket end)");
                    m_result.output_Sentry_lcp( repcount+1 );

                    mmlcp.setRepbucketEnd( relRank-1, repcount+1 );
                }
                else
                {
                    offset_type m = repcount+1; // base value is the repcount.

                    while ( !MaxRepcountStack.empty() && MaxRepcountStack.top().charbkt > charLimit )   // pop all unneeded characters
                        MaxRepcountStack.pop();

                    if ( MaxRepcountStack.empty() )
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack is empty!");
                        m = 0;
                    }
                    else if ( MaxRepcountStack.top().charbkt != charLimit )
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack char " << strC(MaxRepcountStack.top().charbkt) << " != " << charLimit);
                        m = 0;
                    }
                    else
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack char " << strC(MaxRepcountStack.top().charbkt) << " matches, using max repcount value " << MaxRepcountStack.top().maxrepcount);
                        m = std::min(m, MaxRepcountStack.top().maxrepcount);
                    }

                    DBG(debug_induceS_lcp, "Stack max repcount for char " << strC(charLimit) << " -> result = " << m);
                    DBG(debug_induceS, "--> LCP " << m << " (= L/S-seam output)");

                    m_result.output_Sentry_lcp( m );

                    if ( !MaxRepcountStack.empty() && MaxRepcountStack.top().charbkt == charLimit )
                    {
                        DBG(debug_induceS_lcp, "Buckets match, applying minimum from L*-RepcountStack = " << MaxRepcountStack.top().lstar_repcount);
                        m = std::min(m, MaxRepcountStack.top().lstar_repcount);
                    }
                    else
                    {
                        m = 0;
                    }

                    DBG(debug_induceS, "++> LCP " << m << " (= L/S-seam L*-RMQ value)");
                    mmlcp.setRepbucketEnd( relRank-1, m );
                }
#endif // ESAIS_LCP_CALC_EXT

                prevCharLimit = charLimit;

                if (!Spq->empty())
                    DBG(debug_induceS, "S-PQ.top after break: " << Spq->top());
            }

            delete this->LStarArray;
            SArray.finish_clear();
            SMergeBuffer.finish_clear();

            DBGMEM("induce-only finished");

#if ESAIS_LCP_CALC_EXT
            DBGMEM("induce answering RMQs");

            DBG(debug_induceS, "====================================================================================================");

            RMQ_Stack_Blocked<offset_type,offset_type,lcprmq_slabsize> rmqstruct;

            while ( !LCPpq->empty() )
            {
                LCPTuple t = LCPpq->top(); LCPpq->pop();

                if (t.target % 1000 == 0) {
                    //rmqstruct.print_size();
                }

                DBG(debug_induceS_lcp, "Extracted from LCP-PQ: " << t);

                if (t.type == LCP_SETTER || t.type == LCP_SETTER_OUTPUT)
                {
                    // collect more answers
                    while ( !LCPpq->empty() && LCPpq->top().target == t.target && LCPpq->top().type == t.type ) {
                        DBG(debug_induceL_lcp, "Min-merged with LCP-PQ: " << LCPpq->top());
                        t.v1 = std::min(t.v1, LCPpq->top().v1);
                        LCPpq->pop();
                    }

                    rmqstruct.set(t.target, t.v1);

                    if (t.type == LCP_SETTER_OUTPUT)
                        m_result.output_Sentry_lcp( t.v1 );
                }
                else if (t.type == LCP_QUERY)
                {
                    const offset_type &left = t.v1, &right = t.target;
                    assert( left <= right );

                    const offset_type& lcp = rmqstruct.query(left,right) + 1;

                    LCPTuple l = { LCP_SETTER_OUTPUT + (t.type == LCP_QUERY_LSTAR ? 1 : 0), t.v2, lcp, 0 };
                    DBG(debug_induceS_lcp, "Reinserting answer LCP-tuple: " << l);
                    LCPpq->push(l);
                }
                else if (t.type == LCP_LS_SEAM_QUERY)
                {
                    alphabet_type a2 = t.v2;

                    offset_type m = t.v1;     // repcount from S-sequence

                    while ( !MaxRepcountStack.empty() && MaxRepcountStack.top().charbkt > a2 )
                        MaxRepcountStack.pop();

                    if ( MaxRepcountStack.empty() )
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack is empty!");
                        m = 0;
                    }
                    else if ( MaxRepcountStack.top().charbkt != a2 )
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack char " << strC(MaxRepcountStack.top().charbkt) << " != " << a2);
                        m = 0;
                    }
                    else
                    {
                        DBG(debug_induceS_lcp, "Top max repcount stack char " << strC(MaxRepcountStack.top().charbkt) << " matches, using max repcount value " << MaxRepcountStack.top().maxrepcount);
                        m = std::min(m, MaxRepcountStack.top().maxrepcount);
                    }

                    DBG(debug_induceS_lcp, "Stack max repcount for char " << strC(a2) << " -> result = " << m);

                    m_result.output_Sentry_lcp( m );

                    if ( !MaxRepcountStack.empty() && MaxRepcountStack.top().charbkt == a2 )
                    {
                        DBG(debug_induceS_lcp, "Buckets match, applying minimum from L*-RepcountStack = " << MaxRepcountStack.top().lstar_repcount);
                        m = std::min(m, MaxRepcountStack.top().lstar_repcount);
                    }
                    else
                    {
                        m = 0;
                    }

                    rmqstruct.set(t.target, m);
                }
                else {
                    assert(!"Unknown RMQ-Tuple");
                }
            }

            rmqstruct.print_size();

#endif // ESAIS_LCP_CALC_EXT

            m_result.finish_Ssequence();

            DBGMEM("induce finished");
        }

        // ***************************************************************************************************
        // *** Induce process in three steps

        /// Process the input stream in reverse together with the correctly ranked S*-positions. The
        /// SStarRankStream is instanciated with a runs_merger<NameTuple> for non-recursion cases and with
        /// RecursionInduceFilter for recursive calls, which returns NameTuples for each S*-position.
        template <typename InputStreamReverse, typename SStarRankStream, typename SStarLCPStream>
        void process(InputStreamReverse& inputrev, SStarRankStream& sstarrankstream, SStarLCPStream& sstarlcpstream, unsigned int depth)
        {
            DBGMEM("Induce process() started");

            process_input(inputrev, sstarrankstream);

            DBGMEM("Induce process() done");

            induceL(sstarlcpstream);
            induceS();

            DBG(debug, "Merge counter: " << m_mergecounter);
            g_statscache >> "tuplemerges" >> depth << m_mergecounter;
        }
    };

    // }}}

    /*********************************************************************************************************
     * Step 3: Recursion - input filter and output merger                                                    *
     *********************************************************************************************************/

    /// This stream class splits up the NameTuple of (index,name) prepared by run_merger into just a sequence
    /// of (name). This basically adapts the standard stxxl::sorter to the input sequence, which is read
    /// backwards.
    template <typename st_nametuple_type>
    class RecursionInputFilter
    {
    private:

        /// sorter NameTuple (index,name) stream
        st_nametuple_type&              m_st_nametuple;

        /// value type
        typedef offset_type             value_type;

    public:

        RecursionInputFilter(st_nametuple_type& st_nametuple)
            : m_st_nametuple(st_nametuple)
        {
        }

        RecursionInputFilter& rewind(size_t mem)
        {
            m_st_nametuple.sort(mem);
            return *this;
        }

        RecursionInputFilter& finish()
        {
            m_st_nametuple.finish();
            return *this;
        }

        value_type operator*() const
        {
            return m_st_nametuple->name;
        }

        value_type curr_index() const
        {
            return m_st_nametuple->index;
        }

        RecursionInputFilter& operator--()
        {
            ++m_st_nametuple;
            return *this;
        }

        bool empty() const
        {
            return m_st_nametuple.empty();
        }

        size_type size() const
        {
            return m_st_nametuple.size();
        }
    };

    template <typename InputPairs, typename InputISA>
    class RecursionInduceFilter
    {
    protected:

        InputPairs&     m_in1;

        InputISA&       m_in2;

        typedef RankTuple value_type;

        RankTuple       m_current;

    public:

        RecursionInduceFilter(InputPairs& in1, InputISA& in2)
            : m_in1(in1), m_in2(in2)
        {
            m_current = RankTuple(m_in1.curr_index(), m_in2->second);
        }

        bool empty() const
        {
            return (m_in1.empty() || m_in2.empty());
        }

        void finish_clear()
        {
            m_in2.finish_clear();      // deallocate ISA array
        }

        RecursionInduceFilter& operator++()
        {
            --m_in1; ++m_in2;
            if (!empty())
                m_current = RankTuple(m_in1.curr_index(), m_in2->second);
            return *this;
        }

        const value_type& operator*() const
        {
            return m_current;
        }

        const value_type* operator->() const
        {
            return &m_current;
        }
    };

    class SeqLCPStream : public ShortStringSorter::lcp_seq_type::stream
    {
    private:
        offset_type     index;

    public:
        SeqLCPStream(const typename ShortStringSorter::lcp_seq_type& sequence)
            : ShortStringSorter::lcp_seq_type::stream(sequence),
              index(0)
        {
        }

        void sort(size_t)
        {
        }

        SeqLCPStream& operator++() {
            ++index;
            ShortStringSorter::lcp_seq_type::stream::operator++();
            return *this;
        }

        const offset_type& curr_index() const { return index; }
    };

    /// Dummy LCP stream returning only zeros as needed for the base-case where all lcps are zero.
    class ZeroLCPStream
    {
    public:
        typedef offset_type     value_type;

    protected:
        value_type      m_current;

    public:
        ZeroLCPStream() : m_current(0) {}

        bool empty() const { return false; }

        ZeroLCPStream& operator++() { return *this; }

        const value_type& operator*() const { return m_current; }

        void sort(size_t) {}

        const offset_type& curr_index() const { return 0; }
    };

    template <typename SorterLCPrangesumAnswer, typename SorterRMQresult>
    class SStarLCPStream
    {
    public:
        typedef offset_type     value_type;

    protected:
        bool            m_empty;
        value_type      m_current;
        offset_type     m_index;

        SorterLCPrangesumAnswer&  LCPrangesumAnswer;
        SorterRMQresult&        RMQresult;

    public:
        SStarLCPStream(SorterLCPrangesumAnswer& _LCPrangesumAnswer, SorterRMQresult& _RMQresult)
            : m_empty(false),
              m_current(0),
              m_index(0),               // first lcp is sentinel zero
              LCPrangesumAnswer(_LCPrangesumAnswer),
              RMQresult(_RMQresult)
        {
        }

        void sort(const size_t& memsize)
        {
            LCPrangesumAnswer.sort(memsize / 2);
            RMQresult.sort(memsize / 2);
        }

        bool empty() const
        {
            return m_empty;
        }

        SStarLCPStream& operator++()
        {
            if (RMQresult.empty()) {
                m_empty = true;
                return *this;
            }

            m_index = RMQresult->first;
            m_current = RMQresult->second;
            ++RMQresult;

            DBG(debug_sstarlcp_merge, "Next S*-lcp at position " << m_index << " current minimum " << m_current);

            while ( !RMQresult.empty() && RMQresult->first == m_index )
            {
                m_current = std::min( m_current, RMQresult->second );
                ++RMQresult;
                DBG(debug_sstarlcp_merge, "After another RMQanswer: " << m_current);
            }

            if ( !LCPrangesumAnswer.empty() && LCPrangesumAnswer->first == m_index )
            {
                offset_type rangesum1 = LCPrangesumAnswer->second;

                ++LCPrangesumAnswer;
                assert( !LCPrangesumAnswer.empty() );   // range sums always come in pairs
                assert( LCPrangesumAnswer->first == m_index );

                offset_type rangesum2 = LCPrangesumAnswer->second;
                ++LCPrangesumAnswer;

                offset_type rangesum = (rangesum1 < rangesum2) ? (rangesum2 - rangesum1) : (rangesum1 - rangesum2);

                m_current += rangesum;

                DBG(debug_sstarlcp_merge, "Added range sum = " << rangesum);
            }

            DBG(debug_sstarlcp_merge, "Final S*-lcp = " << m_current);

            return *this;
        }

        const value_type& operator*() const
        {
            return m_current;
        }

        const offset_type& curr_index() const
        {
            return m_index;
        }
    };

#if ESAIS_DISCARD_UNIQUES
    /// This stream class splits up the sequences of (index,name) prepared by run_merger into just a sequence
    /// of (name). Therewhile discarding sequences of unique substrings, leaving only the first unique name.
    template <typename st_nametuple_type>
    class RecursionInputFilterDiscard
    {
    private:

        /// sorter NameTuple (index,name) stream
        st_nametuple_type&              m_st_nametuple;

        /// value type
        typedef offset_type             value_type;

        /// the current value of this stream, while the runs_merger is always one step ahead.
        typename st_nametuple_type::value_type        m_current;

        /// end of stream, extra flag needed because it does not coincide with runs_merger.empty()
        bool                            m_empty;

        /// total size
        size_type                       m_totalsize;

    public:

        RecursionInputFilterDiscard(st_nametuple_type& st_nametuple)
            : m_st_nametuple(st_nametuple)
        {
            m_empty = false;
            operator--();

            m_totalsize = 0;    // count number of characters after discarding uniques
            while( !empty() ) {
                ++m_totalsize;
                operator--();
            }

            DBG(debug, "Total size after discarding: " << m_totalsize
                << " - original: " << m_st_nametuple.size() << " - "
                << (m_totalsize * 100.0 / m_st_nametuple.size()) << " percent");
        }

        RecursionInputFilterDiscard& rewind(size_t mem)
        {
            m_st_nametuple.sort(mem);
            m_empty = false;
            operator--();

            return *this;
        }

        RecursionInputFilterDiscard& finish()
        {
            m_st_nametuple.finish();
            return *this;
        }

        value_type operator*() const
        {
            return m_current.name;
        }

        value_type curr_index() const
        {
            return m_current.index;
        }

        RecursionInputFilterDiscard& operator--()
        {
            assert(!m_empty);

            if (m_st_nametuple.empty()) {
                m_empty = true;
                return *this;
            }

            m_current = *m_st_nametuple;
            ++m_st_nametuple;

            while ( !m_st_nametuple.empty() && m_current.unique && m_st_nametuple->unique ) {
                m_current = *m_st_nametuple;
                ++m_st_nametuple;
            }
            return *this;
        }

        bool empty() const
        {
            return m_empty;
        }

        size_type size() const
        {
            return m_totalsize;
        }
    };
#endif

     /// Sorts left RMQ tuples boundaries by (i DIV S, tgt) with S the RMQ slab size
    template <int S>
    struct indexpair_RMQ
    {
        typedef tuples::pair<offset_type,offset_type>      value_type;

        bool operator()(const value_type& a, const value_type& b) const {
            if ( (a.first / S) == (b.first / S) )
                return a.second < b.second;
            else
                return (a.first / S) < (b.first / S);
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

     /// Sorts right RMQ tuples boundaries by (i DIV S, tgt) with S the RMQ slab size
    template <typename TripleType, int S>
    struct indextriple_RMQ
    {
        typedef TripleType      value_type;

        bool operator()(const value_type& a, const value_type& b) const {
            if ( (a.first / S) == (b.first / S) )
                return a.second < b.second;
            else
                return (a.first / S) < (b.first / S);
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /*********************************************************************************************************
     * eSAIS Main                                                                                            *
     *********************************************************************************************************/

    /// instance of inducing process's result which delivers the output stream
    Result*             m_result;

public:

    eSAIS()
        : m_result(NULL)
    {
        assert( sizeof(alphabet_type) <= sizeof(offset_type) );
    }

    ~eSAIS()
    {
        if (m_result) delete m_result;
    }

    template <typename InputStreamReverse>
    unsigned int run(InputStreamReverse& inputrev, unsigned int depth = 0)
    {
        TimerSeries timer;
        timer.record("start");

        unsigned int maxdepth;

        DBGMEM("run()");

        DBG(debug, "inputsize=" << inputrev.size() << ", " <<
            "sizeof(CTuple)=" << sizeof(typename Induce::CTuple) << ", "
            "sizeof(PQTuple)=" << sizeof(typename Induce::PQTuple));

        // ***************************************************************************************************
        // *** Typedefs of sorter for index pairs of (index,rank)

        // pair of (index,name) created by lexnaming in substring sorter
        typedef tuples::pair<offset_type,offset_type> indexpair_type;

        // comparator type that will sort by first component: index
        typedef tuples::pair_less1st<indexpair_type> indexpair_less1st_type;
        typedef tuples::pair_greater1st<indexpair_type> indexpair_greater1st_type;

        indexpair_greater1st_type indexpair_greater1st;

        // complete sorter for indexpair_type
        typedef stxxl::sorter<indexpair_type, indexpair_greater1st_type, block_size> st_indexpair_byindex_type;

        // ***************************************************************************************************
        // *** Read input from back to front, split at S*-positions or block boundaries and feed substrings to
        // *** sorter

        ShortStringSorter sss;

        DBGMEM("ShortStringSorter start");

        sss.process_input( inputrev, depth );

        DBGMEM("ShortStringSorter process_input() done");
        timer.record("sssort").printlast();

        // ***************************************************************************************************
        // *** Get sorted NameTuples (name,index,uniqueflag,length) via merge of sorted sequences and sort
        // *** output by index

        // sorter for NameTuples (index,name,uniqueflag,length)
        typedef stxxl::sorter<NameTuple, NameTupleOrder_IndexDesc, block_size> st_nametuple_byindex_type;

        st_nametuple_byindex_type    st_nametuple(NameTupleOrder_IndexDesc(), memsize, memsize / 8);

        // run substring merger: merge sorted sequences do lexnaming and push (index,name,...) into runs creator
        sss.merge( st_nametuple );

        DBG_ST_ARRAY(debug_lexnamepairs, "NameTuple list", st_nametuple);

        DBGMEM("ShortStringSorter merge() done");
        timer.record("ssname").printlast();

        DBG(debug, "Substring tuples: " << sss.totalsize()
            << " total with " << sss.duplicates() << " duplicates and " << sss.totalsize() - sss.duplicates() << " uniques,"
            << " ratio = " << sss.totalsize() * 1.0 / sss.duplicates());

        g_statscache >> "S*ss-duplicates" >> depth << sss.duplicates();
        g_statscache >> "S*ss-uniques" >> depth << (sss.totalsize() - sss.duplicates());

        // ***************************************************************************************************
        // *** Depending of whether names were unique: go into recursion or use names as ranks.

        if ( sss.was_unique() )
        {
            // base case with unique names: use (index,rank) directly for PQ-creation

            DBG(debug, " No recursion needed, depth " << depth << " *************************************************************");

            if (debug_recursive_input || debug_induce_input)
            {
                // sorter yields (index,rank) order by index DESC.
                DBG_ST_ARRAY(debug_recursive_input || debug_induce_input, "Base-case Input for Induce (in reverse)", st_nametuple);

#if ESAIS_LCP_CALC
                typename ShortStringSorter::lcp_seq_type::stream lcpstream ( sss.lexname_lcp_seq );

                DBG(debug_recursive_input || debug_induce_input, "Base-case LCP Stream for Induce (in SA-order): ");

                while ( !lcpstream.empty() )
                {
                    DBG(debug_recursive_input || debug_induce_input, *lcpstream);
                    ++lcpstream;
                }
#endif // ESAIS_LCP_CALC
            }

            // sorter yields (index,rank) order by index DESC.
            st_nametuple.sort(memsize / 4);

#if ESAIS_LCP_CALC
            SeqLCPStream lcpstream ( sss.lexname_lcp_seq );
#else // !ESAIS_LCP_CALC
            ZeroLCPStream lcpstream;
#endif // ESAIS_LCP_CALC

            m_result = new Result;

            std::auto_ptr<Induce> induce (new Induce(m_result));
            induce->process( inputrev, st_nametuple, lcpstream, depth );

            DBGMEM("Induce process() done");
            timer.record("done").printall();

            // result of induce can be read by caller via stream interface.
            return depth;
        }
#if ESAIS_DISCARD_UNIQUES
        else if (sss.totalsize() / sss.duplicates() >= 3)     // unique discarding
        {
            // ***********************************************************************************************
            // *** recursion case with duplicate names: create recursion string by removing second component
            // *** and discarding sequences of unique names

            DBG(debug, "*** Entering Recursion with unique discarding, depth " << depth << " ************************************************************");

#if ESAIS_LCP_CALC
#error "LCP Calculation not implemented with discarding."
#endif // ESAIS_LCP_CALC

            // recursion filter: extracts sequences of (name) from pairs and indexes in a queue

            st_nametuple.sort(memsize / 8);
            typedef RecursionInputFilterDiscard<st_nametuple_byindex_type> rec_input_type;
            rec_input_type rec_input( st_nametuple );

            ESAIS_LCP_CALCX( assert(sss.SStarSize.size() == rec_input.size()) );

            DBGMEM("Built recursion filter");

            if (debug_recursive_input)
            {
                DBG(debug_recursive_input, "Recursive Input (in reverse):");
                while ( ! rec_input.empty() )
                {
                    DBG(debug_recursive_input, *rec_input << " (idx " << rec_input.curr_index() << ")");
                    --rec_input;
                }
                DBG(debug_recursive_input, "End");

                rec_input.rewind(memsize / 8);
            }

            // --- run recursion

            DBGMEM("Enter recursion");
            timer.record("recstart").printlast();

            eSAIS<offset_type, offset_type, size_type> esais_recurse;
            maxdepth = esais_recurse.run( rec_input, depth+1 );

            DBG(debug, "*** Exiting Recursion, depth " << depth << " *************************************************************");

            timer.record("recdone").printlast();
            DBGMEM("Exit recursion");

            // --- process recursion result

            // sorted_runs to contain result SA and build inverse suffix array
            st_indexpair_byindex_type   st_ISA(indexpair_greater1st, memsize / 2);

            DBGMEM("Created result sorters");

#if ESAIS_SELF_CHECK
            std::vector<offset_type>    recoutput;
#endif

            for ( size_type i = 0; !esais_recurse.empty(); ++i, ++esais_recurse )
            {
                DBG(debug_recursive_output, "recursive SA[" << i << "] = " << *esais_recurse ESAIS_LCP_CALCX(<< " with lcp[" << i << "] = " << esais_recurse.lcp()));

                st_ISA.push( indexpair_type(*esais_recurse, i) );

#if ESAIS_SELF_CHECK
                recoutput.push_back( *esais_recurse );
#endif
            }

            st_ISA.finish();

#if ESAIS_SELF_CHECK
            {
                // read recursive input into memory
                rec_input.rewind(memsize / 4);
                std::vector<offset_type> recinput(rec_input.size());
                for(size_type i = recinput.size(); !rec_input.empty(); --rec_input)
                    recinput[--i] = *rec_input;

                if (debug_recursive_output)
                {
                    std::cout << "\nResulting suffix array: \n";

                    for (size_t i = 0; i < recoutput.size(); ++i)
                    {
                        std::cout << i << " : " << recoutput[i] << " : ";

                        for (size_t j = 0; recoutput[i]+j < recinput.size() && j < 20; ++j)
                            std::cout << strC(recinput[recoutput[i]+j]) << " ";

                        std::cout << "\n";
                    }
                }

                (std::cout << " checking Recursion ").flush();

                if (! sacheck::check_sa_vectors(recinput, recoutput) )
                    std::cout << "failed!" << std::endl;
                else
                    std::cout << "ok" << std::endl;
            }
#endif // ESAIS_SELF_CHECK

            DBGMEM("Iterating over ISA build RerankTuples");

            // final sorter yielding (index,rank) by index of all S*-positions, used by Induce
            stxxl::sorter<RankTuple, RankTupleOrder_IndexDesc, block_size> st_sstarranks (RankTupleOrder_IndexDesc(), memsize / 4);

            // *** Merge (i,ISA) from recursion with original recstr (i,name[i]). Those that were discarded
            // *** are put into the final sorter, others are merged with recISA into a RerankTuple

            // runs creator for RerankTuple (index,rank,ISA)
            stxxl::sorter<RerankTuple, RerankTuple_ISA, block_size> st_reranktuple (RerankTuple_ISA(), memsize / 4);

            {
                // runs merger yields (i,ISA[i]) order by index DESC.
                st_ISA.sort(memsize / 4);

                // runs merger yields (index,rank) order by index DESC.
                st_nametuple.sort(memsize / 4);

                NameTuple current = *st_nametuple;
                ++st_nametuple;

                while (!st_nametuple.empty())
                {
                    if (current.unique && st_nametuple->unique)
                    {
                        // current was discarded
                        st_sstarranks.push( RankTuple(current.index, current.name) );
                        DBG(debug_rerank, "Discarded from recursion: (" << current.index << "," << current.name << ")");
                    }
                    else
                    {
                        assert( !st_ISA.empty() );

                        RerankTuple rt(current.index, current.name, st_ISA->second);
                        st_reranktuple.push(rt);
                        ++st_ISA;

                        DBG(debug_rerank, "Reranking tuple " << rt);
                    }

                    current = *st_nametuple;
                    ++st_nametuple;
                }

                // process last tuple
                assert( !st_ISA.empty() );

                RerankTuple rt(current.index, current.name, st_ISA->second);
                st_reranktuple.push(rt);
                ++st_ISA;

                assert( st_ISA.empty() );
                DBG(debug_rerank, "Reranking tuple " << rt);
            }

            st_nametuple.finish_clear();
            st_ISA.finish_clear();

            DBGMEM("Reranking Tuples");

            st_reranktuple.sort(memsize / 2);

            {
                // *** Iterate over RerankTuples (index,name,ISA) and rerank names into ranks ascending with ISA.

                offset_type lastrank = std::numeric_limits<offset_type>::max();
                size_type rankcounter = 0;

                while (!st_reranktuple.empty())
                {
                    DBG(debug_rerank, "Rerank tuple: " << *st_reranktuple);

                    if (st_reranktuple->name != lastrank) {
                        lastrank = st_reranktuple->name;
                        rankcounter = 0;
                    }
                    else {
                        rankcounter++;
                    }

                    DBG(debug_rerank, "--- reranking to " << lastrank+rankcounter);
                    st_sstarranks.push( RankTuple(st_reranktuple->index,lastrank+rankcounter) );

                    ++st_reranktuple;
                }
            }

            st_reranktuple.finish_clear();
            st_sstarranks.sort(memsize / 8);

            DBG_ST_ARRAY(debug_rerank, "SStarRanks", st_sstarranks);

            DBGMEM("Induce created");
            timer.record("rerank").printlast();

            m_result = new Result;
            ZeroLCPStream lcpstream;

            std::auto_ptr<Induce> induce (new Induce(m_result));
            induce->process( inputrev, st_sstarranks, lcpstream, depth );

            DBGMEM("Induce process() done");
            timer.record("done").printall();

            // result of induce can be read by caller via stream interface.
            return maxdepth;
        }
#endif // ESAIS_DISCARD_UNIQUES
        else
        {
            // ***********************************************************************************************
            // *** recursion case with duplicate names: create recursion string by removing second component.

            DBG(debug, "*** Entering Recursion, depth " << depth << " ************************************************************");

            // recursion filter: extracts sequences of (name) from pairs and indexes in a queue

            st_nametuple.sort(memsize / 8);
            typedef RecursionInputFilter<st_nametuple_byindex_type> rec_input_type;
            rec_input_type rec_input( st_nametuple );

            DBGMEM("Built recursion filter");

            if (debug_recursive_input)
            {
                DBG(debug_recursive_input, "Recursive Input (in reverse):");
                while ( ! rec_input.empty() )
                {
                    DBG(debug_recursive_input, *rec_input << " (idx " << rec_input.curr_index() << ")");
                    --rec_input;
                }
                DBG(debug_recursive_input, "End");

                rec_input.rewind(memsize / 4);
            }

            // --- run recursion

            DBGMEM("Enter recursion");
            timer.record("recstart").printlast();

            eSAIS<offset_type, offset_type, size_type> esais_recurse;
            maxdepth = esais_recurse.run( rec_input, depth+1 );

            DBG(debug, "*** Exiting Recursion, depth " << depth << " *************************************************************");

            timer.record("recdone").printlast();
            DBGMEM("Exit recursion");

            std::cout << "g_mainmemlcp: " << g_mainmemlcp << "\n";

            st_nametuple.finish();

            // --- process recursion result

#if !ESAIS_LCP_CALC
            // sorter to read result SA and build inverse suffix array
            st_indexpair_byindex_type   st_ISA(indexpair_greater1st, memsize / 2);
#else
            // sorter to read result SA and build inverse suffix array
            st_indexpair_byindex_type   st_ISA(indexpair_greater1st, memsize / 8);

            // sorter to contain the range queries needed for sum term in LCP
            st_indexpair_byindex_type   st_LCPrangesum(indexpair_greater1st, memsize / 8);

            // sorter to query the ISA array for specific positions, which are later transformed to RMQ
            st_indexpair_byindex_type   st_ISAqueryLeft(indexpair_greater1st, memsize / 8);
            st_indexpair_byindex_type   st_ISAqueryRight(indexpair_greater1st, memsize / 8);

            // direct queries into lcp_names for recursive lcp = 0

            typedef stxxl::sequence< offset_type, block_size >         LCPdirect_seq_type;
            LCPdirect_seq_type        dq_LCPDirect;

#endif // ESAIS_LCP_CALC

            DBGMEM("Created result sorters");

#if ESAIS_SELF_CHECK
            std::vector<offset_type>    recoutput;
#if ESAIS_LCP_CALC
            std::vector<offset_type>    recoutput_lcp;
#endif // ESAIS_LCP_CALC
#endif // ESAIS_SELF_CHECK

#if ESAIS_LCP_CALC
            offset_type prevSA = 0;
#endif // ESAIS_LCP_CALC

            for ( size_type i = 0; !esais_recurse.empty(); ++i, ++esais_recurse )
            {
                DBG(debug_recursive_output, "recursive SA[" << i << "] = " << *esais_recurse ESAIS_LCP_CALCX(<< " with lcp[" << i << "] = " << esais_recurse.lcp()));

                st_ISA.push( indexpair_type(*esais_recurse, i) );

#if ESAIS_LCP_CALC
                if (i != 0)
                {
                    if (esais_recurse.lcp() != offset_type(0))
                    {
                        // add LCP range sum query tuples
                        DBG(debug_sstarlcp, "range sum query pos " << *esais_recurse << " target " << i);
                        DBG(debug_sstarlcp, "range sum query pos " << *esais_recurse + esais_recurse.lcp() << " target " << i);

                        st_LCPrangesum.push(indexpair_type( *esais_recurse, i )); // left position
                        st_LCPrangesum.push(indexpair_type( *esais_recurse + esais_recurse.lcp(), i )); // right position

                        // add ISA query tuples
                        DBG(debug_sstarlcp, "ISA query left " << prevSA + esais_recurse.lcp() << " target " << i);
                        DBG(debug_sstarlcp, "ISA query right " << *esais_recurse + esais_recurse.lcp() << " target " << i);

                        st_ISAqueryLeft.push(indexpair_type( prevSA + esais_recurse.lcp(), i ));
                        st_ISAqueryRight.push(indexpair_type( *esais_recurse + esais_recurse.lcp(), i ));
                    }
                    else
                    {
                        // don't need to query ISA: direct lookup of adjacent lexname lcps
                        DBG(debug_sstarlcp, "LCP=0 direct query into LCP_Names for target " << i);
                        dq_LCPDirect.push_back(i);
                    }
                }

                prevSA = *esais_recurse;
#endif // ESAIS_LCP_CALC

#if ESAIS_SELF_CHECK
                recoutput.push_back( *esais_recurse );
#if ESAIS_LCP_CALC
                recoutput_lcp.push_back( esais_recurse.lcp() );
#endif // ESAIS_LCP_CALC
#endif // ESAIS_SELF_CHECK
            }

            st_ISA.finish();

#if ESAIS_LCP_CALC
            st_LCPrangesum.finish();
            st_ISAqueryLeft.finish();
            st_ISAqueryRight.finish();
#endif // ESAIS_LCP_CALC

            DBGMEM("Done reading recursive input");

#if ESAIS_SELF_CHECK
            {
                // read recursive input into memory
                rec_input.rewind(memsize / 4);
                std::vector<offset_type> recinput(rec_input.size());
                for(size_type i = recinput.size(); !rec_input.empty(); --rec_input)
                    recinput[--i] = *rec_input;

                if (debug_recursive_output)
                {
                    std::cout << "\nResulting suffix array: \n";

                    for (size_t i = 0; i < recoutput.size(); ++i)
                    {
                        std::cout << i << " : " << recoutput[i] << " : ";

#if ESAIS_LCP_CALC
                        std::cout << i << " lcp " << recoutput_lcp[i] << " : ";
#endif // ESAIS_LCP_CALC

                        for (size_t j = 0; recoutput[i]+j < recinput.size() && j < 20; ++j)
                            std::cout << strC(recinput[recoutput[i]+j]) << " ";

                        std::cout << "\n";
                    }
                }

                (std::cout << " checking recursive SA ").flush();

                if (! sacheck::check_sa_vectors(recinput, recoutput) )
                    std::cout << "failed!" << std::endl;
                else
                    std::cout << "ok" << std::endl;

#if ESAIS_LCP_CALC
                (std::cout << " checking recursive LCP ").flush();

                bool good = true;
                for (size_t i = 1; i < recoutput_lcp.size(); ++i)
                {
                    size_type lcp = 0;
                    for (unsigned int j = 0; std::max(recoutput[i-1],recoutput[i]) + j < recinput.size(); ++j)
                    {
                        if (recinput[ recoutput[i-1] + j ] != recinput[ recoutput[i] + j ])
                            break;
                        ++lcp;
                    }

                    if (lcp != recoutput_lcp[i]) {
                        DBG(debug, "recursive LCP incorrect: lcp[" << i << "] = " << recoutput_lcp[i] << " which should be " << lcp);
                        good = false;
                    }
                }
                std::cout << (good ? "good" : "failed!") << "\n";
#endif // ESAIS_LCP_CALC
            }
#endif // ESAIS_SELF_CHECK

#if ESAIS_LCP_CALC

            // sorter to contain the range query answers needed for sum term in LCP
            typedef stxxl::sorter<indexpair_type, indexpair_less1st_type, block_size> st_LCPrangesumAnswer_type;
            indexpair_less1st_type indexpair_less1st;
            st_LCPrangesumAnswer_type st_LCPrangesumAnswer (indexpair_less1st, memsize / 2);

            {
                st_LCPrangesum.sort(memsize / 2);

                size_type suffixsum = 0;        // reading sizes backwards, so this is a suffix sum (instead of prefix sum)
                offset_type index = sss.SStarSize.size();

                while ( !sss.SStarSize.empty() && !st_LCPrangesum.empty() )
                {
                    --index;
                    suffixsum += sss.SStarSize.front();

                    DBG(debug_sstarlcp, index << " - suffixsum = " << suffixsum);

                    while ( !st_LCPrangesum.empty() && st_LCPrangesum->first == index )
                    {
                        DBG(debug_sstarlcp, "answering query pos " << st_LCPrangesum->first << ", target " << st_LCPrangesum->second << " with " << suffixsum);

                        st_LCPrangesumAnswer.push( indexpair_type(st_LCPrangesum->second, suffixsum) );

                        ++st_LCPrangesum;
                    }

                    sss.SStarSize.pop();
                }

                st_LCPrangesum.finish_clear();
                st_LCPrangesumAnswer.finish();

                DBG_ST_ARRAY(debug_sstarlcp, "st_LCPrangesumAnswer", st_LCPrangesumAnswer);
            }

            DBGMEM("Done LCP range sum query");

            st_ISAqueryLeft.sort(memsize / 4);
            st_ISAqueryRight.sort(memsize / 4);
            st_ISA.sort(memsize / 4);

            // *** step 1: answer queries on ISA and regroup RMQs as (target,boundary)

            stxxl::sorter<indexpair_type, indexpair_less1st_type, block_size> st_ISAresult (indexpair_less1st, memsize / 4);

            while ( !st_ISA.empty() && (!st_ISAqueryLeft.empty() || !st_ISAqueryRight.empty()) )
            {
                while ( !st_ISAqueryLeft.empty() && st_ISAqueryLeft->first == st_ISA->first )
                {
                    DBG(debug_sstarlcp, "matching ISA-query pos " << st_ISAqueryLeft->first << " left boundary to " << st_ISA->second << " tgt " << st_ISAqueryLeft->second);
                    st_ISAresult.push( indexpair_type(st_ISAqueryLeft->second, st_ISA->second + 1) );
                    ++st_ISAqueryLeft;
                }

                while ( !st_ISAqueryRight.empty() && st_ISAqueryRight->first == st_ISA->first )
                {
                    DBG(debug_sstarlcp, "matching ISA-query pos " << st_ISAqueryRight->first << " right boundary to " << st_ISA->second << " tgt " << st_ISAqueryRight->second);
                    st_ISAresult.push( indexpair_type(st_ISAqueryRight->second, st_ISA->second) );
                    ++st_ISAqueryRight;
                }

                ++st_ISA;
            }

            st_ISA.finish();
            st_ISAqueryLeft.finish_clear();
            st_ISAqueryRight.finish_clear();

            DBGMEM("Done ISA queries");

            // *** step 2: generate RMQ tuples from the ISA values, depending on the size of the RMQ interval

            static const unsigned int slabsize = memsize / 4 / sizeof(offset_type);

            typedef tuples::triple<offset_type,offset_type,offset_type> indextriple_type;

            stxxl::sorter<indexpair_type, indexpair_RMQ<slabsize>, block_size> st_RMQleft (indexpair_RMQ<slabsize>(), memsize / 3);
            stxxl::sorter<indextriple_type, indextriple_RMQ<indextriple_type, slabsize>, block_size> st_RMQright (indextriple_RMQ<indextriple_type, slabsize>(), memsize / 3);

            st_ISAresult.sort(memsize / 3);

            while ( !st_ISAresult.empty() )
            {
                offset_type target = st_ISAresult->first;
                offset_type left = st_ISAresult->second;

                ++st_ISAresult;
                assert( st_ISAresult->first == target );

                offset_type right = st_ISAresult->second;
                ++st_ISAresult;

                assert( left <= right );

                if ( (left / slabsize) == (right / slabsize) )
                {
                    // case 1: query is within the same slab -> create one right tuple

                    st_RMQright.push( indextriple_type(right, left, target) );

                    DBG(debug_sstarlcp, "Create one right tuple for enclosed RMQ (" << left << "," << right << ") -> target " << target);
                }
                else
                {
                    // case 2: query spans blocks and must be answered using two queries

                    st_RMQleft.push( indexpair_type(left, target) );
                    st_RMQright.push( indextriple_type(right, left, target) );

                    DBG(debug_sstarlcp, "Create two tuples for spanning RMQ (" << left << "," << right << ") -> target " << target);
                }
            }

            st_ISAresult.finish_clear();

            DBGMEM("Done RMQ queries create");

            st_RMQleft.sort(memsize / 8);
            st_RMQright.sort(memsize / 8);

            DBG_ST_ARRAY(debug_sstarlcp, "st_RMQleft", st_RMQleft);
            DBG_ST_ARRAY(debug_sstarlcp, "st_RMQright", st_RMQright);

            // *** step 3: create blocked RMQ structures and answer queries within those blocks

            // sorter to contain the RMQ results for each position (real RMQ is minimum of matches)
            typedef stxxl::sorter<indexpair_type, indexpair_less1st_type, block_size> st_RMQresult_type;
            st_RMQresult_type st_RMQresult (indexpair_less1st, memsize / 4);

            {
                typename ShortStringSorter::lcp_seq_type::stream lcpstream ( sss.lexname_lcp_seq );

                typename LCPdirect_seq_type::stream lcpdirect ( dq_LCPDirect );

                size_type start = 0;
                stxxl::simple_vector<offset_type> buffer (slabsize);

                RMQ_Stack<uint32_t, offset_type> slabrmq;

                for (size_type thisslab = 0; !lcpstream.empty(); ++thisslab)
                {
                    size_type thissize = 0;
                    assert( start / slabsize == thisslab );

                    while ( thissize < slabsize && !lcpstream.empty() )
                    {
                        DBG(debug_sstarlcp, "buffer[" << thissize << "] = " << *lcpstream);
                        buffer[thissize++] = *lcpstream;
                        ++lcpstream;
                    }

                    DBGMEM("Prior to building succinct RMQ data structure");

                    DBG(debug_sstarlcp, "Building RMQ on range [" << start << "," << start+thissize << ")");
                    RMQ_succinct<offset_type, uint32_t> RMQ(buffer.data(), thissize);
                    DBG(debug_sstarlcp, "Done building RMQ");

                    DBGMEM("After building succinct RMQ data structure");

                    while ( !st_RMQleft.empty() && (st_RMQleft->first / slabsize) == thisslab )
                    {
                        // case 1: left RMQ boundary in this slab, right is beyond the border
                        DBG(debug_sstarlcp, "RMQ left in slab: (" << st_RMQleft->first << ",inf) target " << st_RMQleft->second);

                        uint32_t resultpos = RMQ.query( st_RMQleft->first - start, thissize-1 );
                        assert( resultpos < thissize );
                        const offset_type& result = buffer[resultpos];

                        DBG(debug_sstarlcp, "RMQ result " << result);
                        st_RMQresult.push( indexpair_type(st_RMQleft->second, result) );

                        ++st_RMQleft;
                    }

                    while ( !st_RMQright.empty() && (st_RMQright->first / slabsize) == thisslab )
                    {
                        if ( (st_RMQright->second / slabsize) == thisslab )
                        {
                            // case 2: both RMQ boundaries are for this slab and the same target
                            DBG(debug_sstarlcp, "RMQ pair in slab: (" << st_RMQright->second << "," << st_RMQright->first << ") target " << st_RMQright->third);

                            uint32_t resultpos = RMQ.query( st_RMQright->second - start, st_RMQright->first - start );
                            assert( resultpos < thissize );
                            const offset_type& result = buffer[resultpos];

                            DBG(debug_sstarlcp, "RMQ result: " << result);
                            st_RMQresult.push( indexpair_type(st_RMQright->third, result) );

                            ++st_RMQright;
                        }
                        else
                        {
                            // case 3: right RMQ boundary in this slab, left is beyond the border
                            DBG(debug_sstarlcp, "RMQ right in slab: (-inf," << st_RMQright->first << ") target " << st_RMQright->third);

                            uint32_t resultpos = RMQ.query( 0, st_RMQright->first - start );
                            assert( resultpos < thissize );
                            offset_type result = buffer[resultpos];

                            assert( (st_RMQright->second / slabsize) < thisslab );
                            if ( (st_RMQright->second / slabsize) < thisslab-1 ) {
                                const offset_type& min2 = slabrmq.query( (st_RMQright->second / slabsize)+1,  thisslab-1 );
                                result = std::min(result, min2);
                            }

                            DBG(debug_sstarlcp, "RMQ result: " << result);
                            st_RMQresult.push( indexpair_type(st_RMQright->third, result) );

                            ++st_RMQright;
                        }
                    }

                    // Copy LCP_Names directly into RMQresult for all positions with LCP_Recursive = 0
                    while ( !lcpdirect.empty() && *lcpdirect < offset_type(start + thissize) )
                    {
                        DBG(debug_sstarlcp, "Direct LCP_Names copy value " << buffer[*lcpdirect - start] << " for target " << *lcpdirect);

                        st_RMQresult.push( indexpair_type(*lcpdirect, buffer[*lcpdirect - start]) );

                        ++lcpdirect;
                    }

                    // Add global minimum to in-memory stack-structure
                    {
                        uint32_t gminpos = RMQ.query(0,thissize-1);
                        offset_type gmin = buffer[gminpos];
                        DBG(debug_sstarlcp, "Saving slab global minimum " << gmin);
                        slabrmq.set( thisslab, gmin );
                    }

                    start += thissize;
                }
            }

            DBG_ST_ARRAY(debug_sstarlcp, "st_RMQresult", st_RMQresult);

            if (debug_sstarlcp_result)
            {
                SStarLCPStream<st_LCPrangesumAnswer_type,st_RMQresult_type> sstarlcpstream (st_LCPrangesumAnswer, st_RMQresult);
                sstarlcpstream.sort(memsize / 8);

                while ( !sstarlcpstream.empty() )
                {
                    DBG(debug_sstarlcp_result, "S*-LCP[" << sstarlcpstream.curr_index() << "] = " << *sstarlcpstream);
                    ++sstarlcpstream;
                }

                st_LCPrangesumAnswer.finish();
            }

#if ESAIS_SELF_CHECK
            { // check S*-LCPs by reading S*-index list and calculating LCP naively

                (std::cout << " checking S*-LCPs ").flush();

                // read this reclevel's input into memory
                inputrev.rewind(memsize / 4);
                std::vector<alphabet_type> input(inputrev.size());
                for (size_type i = input.size(); !inputrev.empty(); --inputrev)
                    input[--i] = *inputrev;

                // read induce input to find ordered S*-indexes
                std::vector<RankTuple> sstars;

                st_ISA.sort(); rec_input.rewind(memsize / 8);
                RecursionInduceFilter<rec_input_type, st_indexpair_byindex_type> induceinput (rec_input, st_ISA);

                for ( ; !induceinput.empty(); ++induceinput )
                    sstars.push_back( *induceinput );

                std::sort( sstars.begin(), sstars.end(), RankTupleOrder_RankAsc()); // sort by rank

                // naively calculate lcp between S*-indexes in rank order

                SStarLCPStream<st_LCPrangesumAnswer_type,st_RMQresult_type> sstarlcpstream (st_LCPrangesumAnswer, st_RMQresult);
                sstarlcpstream.sort(memsize / 8);

                assert( *sstarlcpstream == offset_type(0) ); // first is zero
                ++sstarlcpstream;

                bool good = true;
                for (unsigned int i = 1; i < sstars.size(); ++i)
                {
                    size_type lcp = 0;
                    for (unsigned int j = 0; std::max(sstars[i].index,sstars[i-1].index) + j < input.size(); ++j)
                    {
                        if (input[ sstars[i-1].index + j ] !=
                            input[ sstars[i].index + j ])
                            break;

                        ++lcp;
                    }

                    if (lcp != *sstarlcpstream) {
                        good = false;
                        DBG(debug, "LCP of S*s " << sstars[i-1] << " and " << sstars[i] << " incorrectly calculated " << *sstarlcpstream << " correct is " << lcp);
                    }

                    ++sstarlcpstream;
                }
                if (good) DBG(debug, "good");
                else assert(!"bad");
            }
#endif // ESAIS_SELF_CHECK

            st_RMQresult.finish();

            DBGMEM("Done with S*-LCP calculation");
#endif // ESAIS_LCP_CALC

            if (debug_induce_input)
            {
                st_ISA.sort();

                rec_input.rewind(memsize / 4);
                RecursionInduceFilter<rec_input_type, st_indexpair_byindex_type> induceinput (rec_input, st_ISA);

                DBG(debug_induce_input, "Input for Induce:");
                for ( ; !induceinput.empty(); ++induceinput )
                    DBG(debug_induce_input, *induceinput);
            }

            // Induce::process_input() memory: inputrev (1/8) + induceinput (2*1/8) + rec_input(1/8) + 3 output arrays (1/6 each)
            st_nametuple.sort(memsize / 8);
            st_ISA.sort(memsize / 8);

            rec_input.rewind(memsize / 8);
            RecursionInduceFilter<rec_input_type, st_indexpair_byindex_type> induceinput (rec_input, st_ISA);

            DBGMEM("Induce created");
            timer.record("isadone").printlast();

            m_result = new Result;

#if ESAIS_LCP_CALC
            SStarLCPStream<st_LCPrangesumAnswer_type,st_RMQresult_type> lcpstream (st_LCPrangesumAnswer, st_RMQresult);
#else
            ZeroLCPStream lcpstream;
#endif // ESAIS_LCP_CALC

            std::auto_ptr<Induce> induce (new Induce(m_result));
            induce->process( inputrev, induceinput, lcpstream, depth );

            DBGMEM("Induce process() done");
            timer.record("done").printall();

            // result of induce can be read by caller via stream interface.
            return maxdepth;
        }
    }

    template <typename StringContainer, typename SuffixArrayContainer>
    void run_vector(const StringContainer& string, SuffixArrayContainer& suffixarray)
    {
        my_vector_bufreader_reverse<typename StringContainer::const_iterator> vectstream(string.begin(), string.end());

        unsigned int maxdepth = run(vectstream);
        g_statscache >> "maxdepth" << maxdepth;

        suffixarray.resize( string.size() );
        stxxl::stream::materialize(*this, suffixarray.begin(), suffixarray.end());
    }

#if ESAIS_LCP_CALC
    template <typename StringContainer, typename SuffixArrayContainer, typename LCPArrayContainer>
    void run_vector_lcp(const StringContainer& string, SuffixArrayContainer& suffixarray, LCPArrayContainer& lcparray)
    {
        my_vector_bufreader_reverse<typename StringContainer::const_iterator> vectstream(string.begin(), string.end());

        unsigned int maxdepth = run(vectstream);
        g_statscache >> "maxdepth" << maxdepth;

        suffixarray.resize( string.size() );
        lcparray.resize( string.size() );

        stxxl::vector_bufwriter<typename SuffixArrayContainer::iterator> sa_writer (suffixarray.begin());
        stxxl::vector_bufwriter<typename LCPArrayContainer::iterator> lcp_writer (lcparray.begin());

        for (size_type i = 0; !empty(); operator++(), ++i)
        {
            //std::cout << "SA[" << i << "] = " << operator*() << " - LCP[" << i << "] = " << lcp() << "\n";
            sa_writer << operator*();
            lcp_writer << lcp();
        }
    }
#endif // ESAIS_LCP_CALC

    /*********************************************************************************************************
     * Stream Interface of eSAIS class                                                                       *
     *********************************************************************************************************/

public:

    typedef offset_type         value_type;

    const offset_type& operator*() const
    {
        assert(m_result);
        return *(*m_result);
    }

#if ESAIS_LCP_CALC
    const offset_type& lcp() const
    {
        assert(m_result);
        return m_result->lcp();
    }
#endif // ESAIS_LCP_CALC

    bool empty() const
    {
        assert(m_result);
        return m_result->empty();
    }

    eSAIS& operator++()
    {
        assert(m_result);
        ++(*m_result);
        return *this;
    }
};


template <typename StringContainer, typename SuffixArrayContainer>
struct SACA
{
    std::string name()
    {
#if ESAIS_LCP_CALC_EXT
        return "eSAISlcpext";
#elif ESAIS_LCP_CALC_INT
        return "eSAISlcpint";
#elif ESAIS_DISCARD_UNIQUES
        return "eSAISdisc";
#else
        return "eSAIS";
#endif // ESAIS_DISCARD_UNIQUES
    }

    void prepare(const StringContainer& /* string */, SuffixArrayContainer& /* suffixarray */, unsigned int /* K */)
    {
    }

    void run(const StringContainer& string, SuffixArrayContainer& suffixarray, unsigned int /* K */)
    {
        eSAIS<typename StringContainer::value_type, typename SuffixArrayContainer::value_type, stxxl::uint64> algo;

        g_statscache >> "tuplecharlimit" << int(algo.D);

        algo.run_vector(string, suffixarray);

#if ESAIS_COUNT_WASTED
        g_statscache >> "wasted_iovolume" << g_wasted_iovolume;
#endif
    }

#if ESAIS_LCP_CALC
    template <typename LCPArrayContainer>
    void run_lcp(const StringContainer& string, SuffixArrayContainer& suffixarray, LCPArrayContainer& lcparray, unsigned int /* K */)
    {
        eSAIS<typename StringContainer::value_type, typename SuffixArrayContainer::value_type, stxxl::uint64> algo;

        g_statscache >> "tuplecharlimit" << int(algo.D);

#if ESAIS_SELF_CHECK
        algo.test_lexcompare();
#endif // ESAIS_SELF_CHECK

        algo.run_vector_lcp(string, suffixarray, lcparray);

#if ESAIS_COUNT_WASTED
        g_statscache >> "wasted_iovolume" << g_wasted_iovolume;
#endif
        g_statscache >> "mainmemlcp" << g_mainmemlcp;
        std::cout << "g_mainmemlcp: " << g_mainmemlcp << "\n";
    }
#endif // ESAIS_LCP_CALC
};

} // namespace esais
