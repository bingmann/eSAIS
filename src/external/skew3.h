// -*- mode: c++; mode: rebox; fill-column: 100; rebox-min-fill-column: 100 -*-

/******************************************************************************
 * src/external/skew3.h
 *
 * Algorithm implementing skew3 aka DC3 in external memory.
 *
 * Originally implemented in STXXL by Jens Mehnert (2004)
 * LCP construction with RMQs and cleanup by Timo Bingmann (2012)
 *
 * Dementiev, R., Kaerkkaeinen, J., Mehnert, J., & Sanders, P. (2008).
 * Better external memory suffix array construction. JEA, Volume 12.
 *
 ******************************************************************************
 * Copyright (C) 2004 Jens Mehnert
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

namespace skew {

using namespace tuples;
using namespace stxxl;

template <typename offset_type>
class skew
{
public:
    /***********************************************************************************************
     * Skew debug configurations                                                                   *
     ***********************************************************************************************/

    static const bool debug = true;

    static const bool debug_input = false;
    static const bool debug_lexname = false;
    static const bool debug_isa = false;
    static const bool debug_merge = false;
    static const bool debug_mcreate = false;
    static const bool debug_S012 = false;
    static const bool debug_merge_result = false;
    static const bool debug_merge_lcp = false;
    static const bool debug_lcp_expand = false;
    static const bool debug_lcp_build = false;
    static const bool debug_rmq = false;

    /***********************************************************************************************
     * Skew numbering overview                                                                     *
     ***********************************************************************************************/

    /*
     * Numbering and alphabet translation must be defined correctly calculate LCPs. The current
     * schema translates in the following steps:
     * a) input -> quads - add +1 to characters when placed into quads to free 0 for termination
     * b) lexnaming of quads returns lexnames starting with 0.
     * c) lexnaming create LCP_N, also based at 0, thus LCP_N[1] = lcp_T(name_0,name_1)
     * d) SA is generated zero-based, as well as ISA.
     * e) during build_sa +1 is added to the ranks from ISA, to free 0 for termination ranks.
     *
     */

    /***********************************************************************************************
     * Skew algorithm configurations                                                               *
     ***********************************************************************************************/

    /** Types */

    typedef tuples::pair<offset_type, offset_type> skew_pair_type;
    typedef tuples::triple<offset_type, offset_type, offset_type> skew_triple_type;
    typedef tuples::quad<offset_type, offset_type, offset_type, offset_type> skew_quad_type;
    typedef tuples::quint<offset_type, offset_type, offset_type, offset_type, offset_type> skew_quint_type;

    /* Perhaps use a templateparameter for the following */

    static const size_t vector_block_size = ::block_size;

    typedef typename stxxl::VECTOR_GENERATOR<offset_type, 8, 2, vector_block_size>::result offset_array_type;
    typedef stream::vector_iterator2stream <typename offset_array_type::iterator> offset_array_it_rg;

#if LCP_CALC
    typedef typename stxxl::VECTOR_GENERATOR<skew_pair_type, 8, 2, vector_block_size>::result skew_pair_array_type;
    typedef stream::vector_iterator2stream <typename skew_pair_array_type::iterator> skew_pair_array_it_rg;

    typedef typename stxxl::VECTOR_GENERATOR<uint8_t, 8, 2, vector_block_size>::result uint8_array_type;
    typedef stream::vector_iterator2stream <typename uint8_array_type::iterator> uint8_array_it_rg;
#endif // LCP_CALC

    /***********************************************************************************************
     * Helper Functions needed by the skew algorithm.                                              *
     ***********************************************************************************************/

    /** comparism function for the mod0 tuples. */

    struct less_mod0
    {
        typedef skew_quint_type value_type;

        bool operator() (const value_type& a, const value_type& b) const
        {
            if (a.second == b.second) {
                return a.fourth < b.fourth;
            }
            else {
                return a.second < b.second;
            }
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /** compares with respect to the first and second component. */

    typedef tuples::pair_less1st<skew_pair_type> less_pair_1st;
    typedef tuples::pair_less2nd<skew_pair_type> less_pair_2nd;

    typedef tuples::triple_less1st<skew_triple_type> less_triple_1st;

    /** compares quadruples with respect to the first component. */

    struct less_quad_2nd
    {
        typedef skew_quad_type value_type;

        bool operator() (const skew_quad_type& a, const skew_quad_type& b) const
        {
            return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /** compares five-tuples with respect to the second component. */

    struct less_quint_2nd
    {
        typedef skew_quint_type value_type;

        bool operator() (const value_type& a, const value_type& b) const
        {
            return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /**
     * Put the (0 mod 2) [which are the 1,2 mod 3 tuples] tuples at the begin.
     */
    struct less_skew
    {
        typedef skew_pair_type value_type;

        bool operator() (const value_type& a, const value_type& b) const
        {
            if ((a.first & 1) == (b.first & 1)) {
                return a.first < b.first;
            }
            else
            {
                return (a.first & 1) < (b.first & 1);
            }
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /** Sort skew_quad datatype. */
    template <typename alphabet_type>
    struct less_quad
    {
        typedef typename tuples::quad<offset_type, alphabet_type, alphabet_type, alphabet_type> value_type;

        bool operator() (const value_type& a, const value_type& b) const
        {
            if (a.second == b.second) {
                if (a.third == b.third) {
                    return a.fourth < b.fourth;
                }
                else {
                    return a.third < b.third;
                }
            }
            else {
                return a.second < b.second;
            }
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /***********************************************************************************************
     * Skew naming.                                                                                *
     ***********************************************************************************************/

    /**
     * Check, if last two components of two quads are equal.
     */

    template <class quad_type>
    static inline bool quad_eq(const quad_type& a, const quad_type& b)
    {
        return (a.second == b.second) && (a.third == b.third) && (a.fourth == b.fourth);
    }

    template <class quad_type>
    static inline uint8_t quad_lcp(const quad_type& a, const quad_type& b)
    {
        if (a.second != b.second) return 0;
        if (a.third != b.third) return 1;
        if (a.fourth != b.fourth) return 2;
        return 3;
    }

    /**
     * Naming pipe for the conventional skew algorithm without discarding.
     */
    template <class Input>
    class naming
    {
    public:
        typedef typename Input::value_type quad_type;

        typedef skew_pair_type  value_type;

    private:
        Input& A;

#if LCP_CALC
        uint8_array_type& lcp_names;
#endif // LCP_CALC

        bool & unique;
        offset_type lexname;
        quad_type prev;
        skew_pair_type result;

    public:

        naming(Input& A_,
#if LCP_CALC
               uint8_array_type& lcp_names_,
#endif // LCP_CALC
               bool & unique_)
            : A(A_),
#if LCP_CALC
              lcp_names(lcp_names_),
#endif // LCP_CALC
              unique(unique_),
              lexname(0)
        {
            assert(!A.empty());
            unique = true;

            prev = *A;
            result.first = prev.first;
            result.second = lexname;

            DBG(debug_lexname, "lexname[" << lexname << "]: " << prev << " and sentinel lcp 0");

#if LCP_CALC
            lcp_names.push_back(0);    // lcp to $ is zero
#endif // LCP_CALC
        }

        const value_type& operator*() const
        {
            return result;
        }

        naming& operator++()
        {
            assert(!A.empty());

            ++A;
            if (A.empty()) return *this;

            quad_type curr = *A;
#if !LCP_CALC
            if (!quad_eq(prev, curr)) {
                ++lexname;
            }
#else // LCP_CALC
            uint8_t lcp;
            if ( (lcp = quad_lcp(prev, curr)) != 3 ) {
                ++lexname;
            }
#endif // LCP_CALC
            else {
                if (!A.empty() && curr.second != offset_type(0)) {
                    unique = false;
                }
            }

#if !LCP_CALC
            DBG(debug_lexname, "lexname[" << lexname << "]: " << curr << " and prev " << prev);
#else // LCP_CALC
            DBG(debug_lexname, "lexname[" << lexname << "]: " << curr << " and prev " << prev << " have lcp " << int(lcp));
            lcp_names.push_back(lcp);
#endif // LCP_CALC

            result.first = curr.first;
            result.second = lexname;

            prev = curr;
            return *this;
        }

        bool empty() const
        {
            return A.empty();
        }
    };

    /***********************************************************************************************
     * A skew algorithm specific namespace.                                                        *
     ***********************************************************************************************/

    /**
     * Create tuples until one of the input streams are empty.
     */
    template<class InputA, class InputB, const size_t add_alphabet = 0>
    class make_pairs
    {
    public:
        // typedef tuples::pair<typename InputA::value_type, offset_type> value_type;
        typedef skew_pair_type value_type;

    private:
        InputA& A;
        InputB& B;

        value_type result;

    public:

        make_pairs(InputA& a, InputB& b)
            : A(a), B(b)
        {
            assert(!A.empty()); assert(!B.empty());
            if (!empty()) {
                result = value_type(*A, *B + add_alphabet);
            }
        }

        const value_type & operator*() const
        {
            return result;
        }

        make_pairs & operator++()
        {
            assert(!A.empty()); assert(!B.empty());

            ++A; ++B;

            if (!A.empty() && !B.empty()) {
                result = value_type(*A, *B + add_alphabet);
            }

            return *this;
        }

        bool empty() const
        {
            return (A.empty() || B.empty());
        }
    };

    /**
     * Collect three characters beginning at the index. Since we need at least one unique
     * endcaracter, we free the first characters.
     * Maps (t_i) -> (i,t_i,t_{i+1},t_{i+2})
     */
    template <class Input, typename alphabet_type, const size_t add_alphabet = 0>
    class make_quads
    {
    public:
        typedef tuples::quad<offset_type, alphabet_type, alphabet_type, alphabet_type> value_type;

    private:
        Input & A;
        value_type current;
        offset_type counter;
        unsigned int z3z;	// = counter mod 3, ("+",Z/3Z) is cheaper than %
        bool finished;

        offset_array_type & backup;

    public:
        make_quads(Input & data_in_, offset_array_type & backup_)
            : A(data_in_),
              current(0,0,0,0),
              counter(0),
              z3z(0),
              finished(false),
              backup(backup_)
        {
            assert(!A.empty());

            current.first = counter;
            current.second = (*A).second + add_alphabet;
            ++A;

            if (!A.empty()) {
                current.third = (*A).second + add_alphabet;
                ++A;
            }
            else {
                current.third = 0;
                current.fourth = 0;
            }

            if (!A.empty()) {
                current.fourth = (*A).second + add_alphabet;
                //++A;
            }
            else {
                current.fourth = 0;
            }
        }

        const value_type & operator*() const
        {
            return current;
        }

        make_quads & operator++()
        {
            assert(!A.empty() || !finished);

            if (current.second != offset_type(0)) {
                DBG(debug_input, "Input[" << counter << "] = " << current.second - add_alphabet);
                backup.push_back(current.second);
            }

            // Calculate module
            if (++z3z == 3) z3z = 0;

            current.first = ++counter;
            current.second = current.third;
            current.third = current.fourth;

            if (!A.empty())
                ++A;

            if (!A.empty()) {
                current.fourth = (*A).second + add_alphabet;
            }
            else {
                current.fourth = 0;
            }

            // Inserts a dummy tuple for input sizes of n%3==1

            if ((current.second == offset_type(0)) && (z3z != 1)) {
                finished = true;
            }

            return *this;
        }

        bool empty() const
        {
            return (A.empty() && finished);
        }
    };

    /**
     * Drop 1/3 of the input. More exactly the offsets at positions (0 mod 3). Index begins with 0.
     */
    template <class Input>
    class extract_mod12
    {
    public:
        typedef typename Input::value_type value_type;

    private:
        Input & A;
        offset_type counter;
        offset_type output_counter;
        value_type result;

    public:

        extract_mod12(Input & A_)
            : A(A_),
              counter(0),
              output_counter(0)
        {
            assert(!A.empty());
            ++A, ++counter; // skip 0 = mod0 offset
            if (!A.empty()) {
                result = *A;
                result.set_1(output_counter);
            }
        }

        const value_type & operator*() const
        {
            return result;
        }

        extract_mod12 & operator++()
        {
            assert(!A.empty());

            ++A, ++counter, ++output_counter;

            if (!A.empty() && (counter % 3) == 0) {
                // skip mod0 offsets
                ++A, ++counter;
            }
            if (!A.empty()) {
                result = *A;
                result.set_1(output_counter);
            }

            return *this;
        }

        bool empty() const
        {
            return A.empty();
        }
    };

#if LCP_CALC
    /**
     * RMQ answering structure
     */

    class RMQ
    {
    private:
        // memsize_compute is constant due to slabsize calculation
        static const size_t memsize_compute = ram_use / 5;

        // slab size of RMQs
        static const size_t slabsize = memsize_compute / 4 / sizeof(offset_type);

        // left query: (left boundary,target) sorted by (left boundary / slabsize)
        typedef tuples::pair<offset_type,offset_type> RMQLeft;

        // right query: (right boundary,left boundary,target) sorted by (right boundary / slabsize)
        typedef tuples::triple<offset_type,offset_type,offset_type> RMQRight;

        // answer tuple: (target,rm answer) sorted by target
        typedef tuples::pair<offset_type,offset_type> RMQAnswer;

        // comparator for RMQ tuples
        template <typename RMQTuple>
        struct RMQTupleOrder
        {
            typedef RMQTuple value_type;

            inline bool operator()(const value_type& a, const value_type& b) const
            {
                return (a.first / slabsize) < (b.first / slabsize);
            }

            static value_type min_value() { return value_type::min_value(); }
            static value_type max_value() { return value_type::max_value(); }
        };

        typedef stxxl::sorter<RMQLeft, RMQTupleOrder<RMQLeft>, block_size> st_RMQLeft_type;
        typedef stxxl::sorter<RMQRight, RMQTupleOrder<RMQRight>, block_size> st_RMQRight_type;

        typedef stxxl::sorter<RMQAnswer, tuples::pair_less1st<RMQAnswer>, block_size> st_RMQAnswer_type;

    private:
        st_RMQLeft_type    st_RMQLeft;
        st_RMQRight_type   st_RMQRight;

        std::auto_ptr<st_RMQAnswer_type>  st_RMQAnswer;

        RMQAnswer          current;
        bool               computed, done;

    public:

        RMQ(size_t memsize_input)
            : st_RMQLeft(RMQTupleOrder<RMQLeft>(), memsize_input / 2),
              st_RMQRight(RMQTupleOrder<RMQRight>(), memsize_input / 2),
              computed(false)
        {
        }

        void query(const offset_type& left, const offset_type& right, const offset_type& target)
        {
            if ( (left / slabsize) == (right / slabsize) )
            {
                // one RMQ tuple
                st_RMQRight.push( RMQRight(right,left,target) );
            }
            else
            {
                // two RMQ tuples: one for left end and one for spanning and right end
                st_RMQLeft.push( RMQLeft(left,target) );
                st_RMQRight.push( RMQRight(right,left,target) );
            }
        }

        template <typename DataStream>
        void compute(DataStream& stream)
        {
            assert(!computed); computed = true;

            st_RMQLeft.sort( memsize_compute / 4 );
            st_RMQRight.sort( memsize_compute / 4 );
            st_RMQAnswer = std::auto_ptr<st_RMQAnswer_type>
                ( new st_RMQAnswer_type(tuples::pair_less1st<RMQAnswer>(), memsize_compute / 4) );

            // buffer takes fourth part of memory

            size_t start = 0;
            offset_type* buffer = new offset_type[slabsize];

            RMQ_Stack<offset_type,offset_type> spanrmq;

            for (size_t thisslab = 0; !stream.empty(); ++thisslab)
            {
                size_t thissize = 0;
                assert( start / slabsize == thisslab );

                // read next slab into buffer
                while ( thissize < slabsize && !stream.empty() )
                {
                    buffer[thissize++] = *stream;
                    ++stream;
                }

                // build rmq structure for slab
                RMQ_succinct<offset_type, uint32_t> slabrmq(buffer, thissize);

                // answer all left queries for this slab
                while ( !st_RMQLeft.empty() && (st_RMQLeft->first / slabsize) == thisslab )
                {
                    // case 1: left RMQ boundary in this slab, right is beyond the border
                    DBG(debug_rmq,"RMQ left in slab: (" << st_RMQLeft->first << ",inf) target " << st_RMQLeft->second);

                    offset_type result = slabrmq.query( st_RMQLeft->first - start, thissize-1 );
                    assert( result < offset_type(thissize) );
                    result = buffer[result];

                    DBG(debug_rmq,"RMQ result " << result);
                    st_RMQAnswer->push( RMQAnswer(st_RMQLeft->second, result) );

                    ++st_RMQLeft;
                }

                // answer all right queries for this slab
                while ( !st_RMQRight.empty() && (st_RMQRight->first / slabsize) == thisslab )
                {
                    if ( (st_RMQRight->second / slabsize) == thisslab )
                    {
                        // case 2: both RMQ boundaries are for this slab and the same target
                        DBG(debug_rmq,"RMQ pair in slab: (" << st_RMQRight->second << "," << st_RMQRight->first << ") target " << st_RMQRight->third);

                        offset_type result;
                        if ( st_RMQRight->first < offset_type(start + thissize) ) {
                            result = slabrmq.query( st_RMQRight->second - start, st_RMQRight->first - start );
                            assert( result < offset_type(thissize) );
                            result = buffer[result];
                        }
                        else {
                            result = 0;
                        }

                        DBG(debug_rmq,"RMQ result: " << result);
                        st_RMQAnswer->push( RMQAnswer(st_RMQRight->third, result) );

                        ++st_RMQRight;
                    }
                    else
                    {
                        // case 3: right RMQ boundary in this slab, left is beyond the border
                        DBG(debug_rmq,"RMQ right in slab: (-inf," << st_RMQRight->first << ") target " << st_RMQRight->third);

                        offset_type result;
                        if ( st_RMQRight->first < offset_type(start + thissize) ) {
                            result = slabrmq.query( 0, st_RMQRight->first - start );
                            assert( result < offset_type(thissize) );
                            result = buffer[result];
                        }
                        else {
                            result = 0;
                        }

                        assert( (st_RMQRight->second / slabsize) < thisslab );
                        if ( (st_RMQRight->second / slabsize) < thisslab-1 ) {
                            const offset_type& min2 = spanrmq.query( (st_RMQRight->second / slabsize)+1,  thisslab-1 );
                            result = std::min(result, min2);
                        }

                        DBG(debug_rmq,"RMQ result: " << result);
                        st_RMQAnswer->push( RMQAnswer(st_RMQRight->third, result) );

                        ++st_RMQRight;
                    }
                }

                // calculate spanning slab rm and save to RMQ stack
                {
                    offset_type gmin = slabrmq.query(0,thissize-1);
                    gmin = buffer[gmin];
                    DBG(debug_rmq,"Saving slab global minimum " << gmin);
                    spanrmq.set( thisslab, gmin );
                }

                start += thissize;
            }

            delete [] buffer;

            st_RMQLeft.finish_clear();
            st_RMQRight.finish_clear();
            st_RMQAnswer->finish();
        }

        void start_output(size_t memsize_output)
        {
            st_RMQAnswer->sort(memsize_output);

            done = false;
            operator++();
        }

        // *** Stream interface to RMQ answers

        typedef RMQAnswer       value_type;

        const value_type & operator*() const
        {
            return current;
        }

        const value_type* operator->() const
        {
            return &current;
        }

        RMQ& operator++()
        {
            if (st_RMQAnswer->empty()) {
                done = true;
                return *this;
            }

            // first answer tuple
            current = *(*st_RMQAnswer);
            ++(*st_RMQAnswer);

            // optional second answer tuple
            if (!st_RMQAnswer->empty() && (*st_RMQAnswer)->first == current.first) {
                current.second = std::min(current.second, (*st_RMQAnswer)->second);
                ++(*st_RMQAnswer);
            }

            return *this;
        }

        bool empty() const
        {
            return done;
        }
    };
#endif // LCP_CALC

    /**
     * Create the suffix array from the current sub problem
     */
    template <class Mod0, class Mod1, class Mod2>
    class merge_sa
    {
    public:
        typedef offset_type value_type;

    private:
        Mod0 & A;
        Mod1 & B;
        Mod2 & C;

        skew_quint_type s0;
        skew_quad_type s1;
        skew_quint_type s2;

        int selected;

        bool exists[3];

        offset_type index;
        offset_type merge_result;

#if LCP_CALC
        // previous winning tuple
        int prev_selected;

        skew_quint_type p0;
        skew_quad_type  p1;
        skew_quint_type p2;

        uint8_array_type& lcp_base;
        RMQ& lcp_rmq;
        skew_pair_array_type& lcp_rmq_ranks;
#endif // LCP_CALC

        /**
         * Constellation of chars and ranks in streams:
         * Mod0 quint: <index, character, character+1, ~SA12+1, ~SA12+2>
         * Mod1 quad : <index, ~SA12 name, character, ~SA12+1 name>
         * Mod2 quint: <index ,~SA12 name, character, character+1, ~SA12+2 name>
         */

        bool cmp_mod12()
        {
            return s1.get_2() < s2.get_2();
        }

        bool cmp_mod02()
        {
            if (s0.get_2() == s2.get_3()) {
                if (s0.get_3() == s2.get_4()) {
                    return s0.get_5() < s2.get_5();
                }
                else {
                    return s0.get_3() < s2.get_4();
                }
            }
            else {
                return s0.get_2() < s2.get_3();
            }
        }

        bool cmp_mod01()
        {
            if (s0.get_2() == s1.get_3()) {
                return s0.get_4() < s1.get_4();
            }
            else {
                return s0.get_2() < s1.get_3();
            }
        }

        void get012()
        {
            DBG1(debug_merge, "0" << s0 << " 1" << s1 << " 2" << s2 << " -> ");

            if (cmp_mod01()) {
                if (cmp_mod02()) {
                    selected = 0;
                    merge_result = s0.get_1();
                }
                else {
                    selected = 2;
                    merge_result = s2.get_1();
                }
            }
            else {
                if (cmp_mod12()) {
                    selected = 1;
                    merge_result = s1.get_1();
                }
                else {
                    selected = 2;
                    merge_result = s2.get_1();
                }
            }

            DBG3(debug_merge, "selected " << selected << " - index " << merge_result);
        }

        void get01()
        {
            DBG1(debug_merge, "0" << s0 << " 1" << s1 << " -> ");

            if (cmp_mod01()) {
                selected = 0;
                merge_result = s0.get_1();
            }
            else {
                selected = 1;
                merge_result = s1.get_1();
            }

            DBG3(debug_merge, "selected " << selected << " - index " << merge_result);
        }

        void get12()
        {
            DBG1(debug_merge, "1" << s1 << " 2" << s2 << " -> ");

            if (cmp_mod12()) {
                selected = 1;
                merge_result = s1.get_1();
            }
            else {
                selected = 2;
                merge_result = s2.get_1();
            }

            DBG3(debug_merge, "selected " << selected << " - index " << merge_result);
        }

        void get02()
        {
            DBG1(debug_merge, "0" << s0 << " 2" << s2 << " -> ");

            if (cmp_mod02()) {
                selected = 0;
                merge_result = s0.get_1();
            }
            else {
                selected = 2;
                merge_result = s2.get_1();
            }

            DBG3(debug_merge, "selected " << selected << " - index " << merge_result);
        }

        void solve()
        {
#if LCP_CALC
            prev_selected = selected;
#endif // LCP_CALC

            if (exists[0]) {
                if (exists[1]) {
                    if (exists[2]) {
                        get012();
                    }
                    else {
                        get01();
                    }
                }
                else {
                    if (exists[2]) {
                        get02();
                    }
                    else {
                        selected = 0;
                        merge_result = s0.get_1();
                        DBG(debug_merge, "get0(): 0" << s0 << " -> selected " << selected << " - index " << merge_result);
                    }
                }
            }
            else {
                if (exists[1]) {
                    if (exists[2]) {
                        get12();
                    }
                    else {
                        selected = 1;
                        merge_result = s1.get_1();
                        DBG(debug_merge, "get1(): 1" << s1 << " -> selected " << selected << " - index " << merge_result);
                    }
                }
                else {
                    if (exists[2]) {
                        selected = 2;
                        merge_result = s2.get_1();
                        DBG(debug_merge, "get2(): 2" << s2 << " -> selected " << selected << " - index " << merge_result);
                    }
                    else {
                        assert(false);
                    }

                }
            }

            DBG(debug_merge_result, merge_result);

#if LCP_CALC
            // calculate lcp of previous tuple to next winner

            if (prev_selected < 0)
            {
                lcp_base.push_back(0); // first lcp is zero
            }
            else if (prev_selected == 0 && selected == 0)
            {
                calc_lcp_char2(p0.get_2(), s0.get_2(),
                               p0.get_3(), s0.get_3(),
                               p0.get_5(), s0.get_5());
            }
            else if (prev_selected == 0 && selected == 1)
            {
                calc_lcp_char1(p0.get_2(), s1.get_3(),
                               p0.get_4(), s1.get_4());
            }
            else if (prev_selected == 0 && selected == 2)
            {
                calc_lcp_char2(p0.get_2(), s2.get_3(),
                               p0.get_3(), s2.get_4(),
                               p0.get_5(), s2.get_5());
            }
            else if (prev_selected == 1 && selected == 0)
            {
                calc_lcp_char1(p1.get_3(), s0.get_2(),
                               p1.get_4(), s0.get_4());
            }
            else if (prev_selected == 1 && selected == 1)
            {
                calc_lcp_char0(p1.get_2(), s1.get_2());
            }
            else if (prev_selected == 1 && selected == 2)
            {
                calc_lcp_char0(p1.get_2(), s2.get_2());
            }
            else if (prev_selected == 2 && selected == 0)
            {
                calc_lcp_char2(p2.get_3(), s0.get_2(),
                               p2.get_4(), s0.get_3(),
                               p2.get_5(), s0.get_5());
            }
            else if (prev_selected == 2 && selected == 1)
            {
                calc_lcp_char0(p2.get_2(), s1.get_2());
            }
            else if (prev_selected == 2 && selected == 2)
            {
                calc_lcp_char0(p2.get_2(), s2.get_2());
            }
            else {
                assert(!"Impossible case");
            }

            if (selected == 0)
                p0 = s0;
            else if (selected == 1)
                p1 = s1;
            else if (selected == 2)
                p2 = s2;
            else
                assert(!"Impossible case");
#endif // LCP_CALC
        }

#if LCP_CALC
        void save_rmq(const offset_type& rankA, const offset_type& rankB)
        {
            if (rankA == offset_type(0)) return;
            assert(rankA > offset_type(0)); assert(rankB > offset_type(0));

            lcp_rmq.query(rankA-1+1, rankB-1, index);
            lcp_rmq_ranks.push_back( skew_pair_type(rankA-1, rankB-1) );
        }

        void calc_lcp_char0(const offset_type& rankA, const offset_type& rankB)
        {
            DBG(debug_merge_lcp, "lcp[" << index << "] = 0 + rmq_lcp(" << rankA-1 << "+1," << rankB-1 << ")");
            assert(rankB > offset_type(0));
            lcp_base.push_back(0);
            save_rmq(rankA, rankB);
        }

        void calc_lcp_char1(const offset_type& chA0, const offset_type& chB0,
                            const offset_type& rankA, const offset_type& rankB)
        {
            if (chA0 != chB0) {
                DBG(debug_merge_lcp, "lcp[" << index << "] = 0");
                lcp_base.push_back(0);
            }
            else {
                DBG(debug_merge_lcp, "lcp[" << index << "] = 1 + rmq_lcp(" << rankA-1 << "+1," << rankB-1 << ")");
                assert(rankB > offset_type(0));
                lcp_base.push_back(1);
                save_rmq(rankA, rankB);
            }
        }

        void calc_lcp_char2(const offset_type& chA0, const offset_type& chB0,
                            const offset_type& chA1, const offset_type& chB1,
                            const offset_type& rankA, const offset_type& rankB)
        {
            if (chA0 != chB0) {
                DBG(debug_merge_lcp, "lcp[" << index << "] = 0");
                lcp_base.push_back(0);
            }
            else if (chA1 != chB1) {
                DBG(debug_merge_lcp, "lcp[" << index << "] = 1");
                lcp_base.push_back(1);
            }
            else {
                DBG(debug_merge_lcp, "lcp[" << index << "] = 2 + rmq_lcp(" << rankA-1 << "+1," << rankB-1 << ")");
                assert(rankB > offset_type(0));
                lcp_base.push_back(2);
                save_rmq(rankA, rankB);
            }
        }
#endif // LCP_CALC

    public:

        bool empty() const
        {
            return (A.empty() && B.empty() && C.empty());
        }

        merge_sa(Mod0 & x1, Mod1 & x2, Mod2 & x3
#if LCP_CALC
                 , uint8_array_type& lcp_base_, RMQ& lcp_rmq_, skew_pair_array_type& lcp_rmq_ranks_
#endif // LCP_CALC
            )
            : A(x1), B(x2), C(x3), selected(-1), index(0)
#if LCP_CALC
              , lcp_base(lcp_base_), lcp_rmq(lcp_rmq_), lcp_rmq_ranks(lcp_rmq_ranks_)
#endif // LCP_CALC
        {
            assert(!A.empty());
            assert(!B.empty());
            assert(!C.empty());
            exists[0] = true;
            exists[1] = true;
            exists[2] = true;
            s0 = *A;
            s1 = *B;
            s2 = *C;

            solve();
        }

        const value_type & operator* () const
        {
            return merge_result;
        }

        merge_sa & operator++ ()
        {
            if (selected == 0) {
                assert(!A.empty());
                ++A;
                if (!A.empty())
                    s0 = *A;
                else
                    exists[0] = false;
            }
            else if (selected == 1) {
                assert(!B.empty());
                ++B;
                if (!B.empty())
                    s1 = *B;
                else
                    exists[1] = false;
            }
            else {
                assert(!C.empty());
                assert(selected == 2);
                ++C;
                if (!C.empty())
                    s2 = *C;
                else
                    exists[2] = false;
            }

            ++index;
            if (!empty()) solve();

            return *this;
        }
    };

    /**
     * Helper function for computing the size of the 2/3 subproblem
     */
    static inline size_t subp_size(size_t n)
    {
        return (n / 3) * 2 + ((n % 3) == 2);
    }

    /***********************************************************************************************
     * A skew algorithm specific namespace.                                                        *
     ***********************************************************************************************/

    /**
     * The merging step of the skew algorithm.
     * @param S input string pipe type.
     * @param Mod1 mod1 tuples input pipe type.
     * @param Mod2 mod2 tuples input pipe type.
     * @param block_size block size of the used sorter.
     * @param memsize memory requirement of one sorter.
     */
    template <class S, class Mod1, class Mod2,
              const int block_size, const size_t memsize>
    class build_sa
    {
    public:
        typedef offset_type value_type;

        static const size_t add_rank = 1;  // free first rank to mark ranks beyond end of input

    private:

        // Mod1 types

        typedef typename stream::use_push < skew_quad_type > mod1_push_type;
        typedef typename stream::runs_creator < mod1_push_type, less_quad_2nd, block_size, RC > mod1_runs_type;
        typedef typename mod1_runs_type::sorted_runs_type sorted_mod1_runs_type;
        typedef typename stream::runs_merger < sorted_mod1_runs_type, less_quad_2nd > mod1_rm_type;

        // Mod2 types

        typedef typename stream::use_push < skew_quint_type > mod2_push_type;
        typedef typename stream::runs_creator < mod2_push_type, less_quint_2nd, block_size, RC > mod2_runs_type;
        typedef typename mod2_runs_type::sorted_runs_type sorted_mod2_runs_type;
        typedef typename stream::runs_merger < sorted_mod2_runs_type, less_quint_2nd > mod2_rm_type;

        // Mod0 types

        typedef typename stream::use_push < skew_quint_type > mod0_push_type;
        typedef typename stream::runs_creator < mod0_push_type, less_mod0, block_size, RC > mod0_runs_type;
        typedef typename mod0_runs_type::sorted_runs_type sorted_mod0_runs_type;
        typedef typename stream::runs_merger < sorted_mod0_runs_type, less_mod0 > mod0_rm_type;

        // Merge type

        typedef merge_sa < mod0_rm_type, mod1_rm_type, mod2_rm_type > merge_sa_type;

        // Functions

        less_mod0 c0;
        less_quad_2nd c1;
        less_quint_2nd c2;

        // Runs merger

        mod1_rm_type *mod1_result;
        mod2_rm_type *mod2_result;
        mod0_rm_type *mod0_result;

        // Merger

        merge_sa_type *vmerge_sa;

        // Input

        S & source;
        Mod1 & mod_1;
        Mod2 & mod_2;

        // Tmp variables

        offset_type t[3];
        offset_type old_t2;
        offset_type old_mod2;
        bool exists[3];
        offset_type mod_one;
        offset_type mod_two;

        size_t  index;

        // Empty_flag

        bool ready;

        // Result

        value_type result;

    public:

        build_sa(S & source_, Mod1 & mod_1_, Mod2 & mod_2_, size_t a_size
#if LCP_CALC
                 , uint8_array_type& lcp_base, RMQ& lcp_rmq1, skew_pair_array_type& lcp_rmq1_ranks
#endif // LCP_CALC
            )
            : source(source_), mod_1(mod_1_), mod_2(mod_2_), index(0), ready(false)
        {
            assert(!source_.empty());

            // Runs storage

            DBGMEM("creating mod0,1,2");

            // input: ISA_1,2 from previous level
            mod0_runs_type mod0_runs(c0, memsize / 3);
            mod1_runs_type mod1_runs(c1, memsize / 3);
            mod2_runs_type mod2_runs(c2, memsize / 3);

            while (!source.empty())
            {
                exists[0] = false;
                exists[1] = false;
                exists[2] = false;

                if (!source.empty()) {
                    t[0] = *source;
                    ++source;
                    exists[0] = true;
                }

                if (!source.empty()) {
                    assert(!mod_1.empty());
                    t[1] = *source;
                    ++source;
                    mod_one = *mod_1 + add_rank;
                    ++mod_1;
                    exists[1] = true;
                }

                if (!source.empty()) {
                    assert(!mod_2.empty());
                    t[2] = *source;
                    ++source;
                    mod_two = *mod_2 + add_rank;
                    ++mod_2;
                    exists[2] = true;
                }

                // Check special cases in the middle of "source"
                // Cases are cx|xc cxx|cxx and cxxc|xxc

                assert(t[0] != offset_type(0));
                assert(t[1] != offset_type(0));
                assert(t[2] != offset_type(0));

                // Mod 0 : (index0,char0,char1,mod1,mod2)
                // Mod 1 : (index1,mod1,char1,mod2)
                // Mod 2 : (index2,mod2)

                if (exists[2]) {	// Nothing is missed
                    DBG(debug_mcreate, "nomiss");
                    DBG(debug_mcreate, "Mod0 : (" << index << "|" << t[0] << "," << t[1] << " - " << mod_one << "," << mod_two << ")");
                    DBG(debug_mcreate, "Mod1 : (" << (index + 1) << "|" << t[1] << " - " << mod_one << "," << mod_two << ")");

                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, mod_two));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], mod_two));

                    if (index != offset_type(0)) {
                        DBG(debug_mcreate, "Mod2 : (" << (index - 1) << "|" << old_t2 << "," << t[0] << " - " << old_mod2 << "," << mod_one << ")");
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                    }
                }
                else if (exists[1]) {	// Last element missed
                    DBG(debug_mcreate, "nolast");
                    DBG(debug_mcreate, "Mod0 : (" << index << "|" << t[0] << "," << t[1] << " - " << mod_one << "," << 0 << ")");
                    DBG(debug_mcreate, "Mod1 : (" << (index + 1) << "|" << t[1] << " - " << mod_one << "," << 0 << ")");

                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, 0));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], 0));

                    if (index != offset_type(0)) {
                        DBG(debug_mcreate, "Mod2 : (" << (index - 1) << "|" << old_t2 << "," << t[0] << " - " << old_mod2 << "," << mod_one << ")");
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                    }
                }
                else {	// Only one element left
                    assert(exists[0]);
                    DBG(debug_mcreate, "twomiss");
                    DBG(debug_mcreate, "Mod0 : (" << index << "|" << t[0] << "," << 0 << " - " << 0 << "," << 0 << ")");

                    mod0_runs.push(skew_quint_type(index, t[0], 0, 0, 0));

                    if (index != offset_type(0)) {
                        DBG(debug_mcreate, "Mod2 : (" << (index - 1) << "|" << old_t2 << "," << t[0] << " - " << old_mod2 << "," << 0 << ")");
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], 0));
                    }
                }

                old_mod2 = mod_two;
                old_t2 = t[2];
                index += 3;
            }

            if ((a_size % 3) == 0) {	// changed
                if (index != offset_type(0)) {
                    DBG(debug_mcreate, "Mod2 : (" << (index - 1) << "|" << old_t2 << "," << 0 << " - " << old_mod2 << "," << 0 << ")");
                    mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, 0, 0));
                }
            }

            if (debug_S012)
            {
                offset_type i = 0;
                for(mod0_rm_type mod0_result(mod0_runs.result(), less_mod0(), memsize / 3);
                    !mod0_result.empty(); ++mod0_result, ++i)
                    DBG(debug_S012, "S0[" << i << "]: " << *mod0_result);

                i = 0;
                for(mod1_rm_type mod1_result(mod1_runs.result(), less_quad_2nd(), memsize / 3);
                    !mod1_result.empty(); ++mod1_result, ++i)
                    DBG(debug_S012, "S1[" << i << "]: " << *mod1_result);

                i = 0;
                for(mod2_rm_type mod2_result(mod2_runs.result(), less_quint_2nd(), memsize / 3);
                    !mod2_result.empty(); ++mod2_result, ++i)
                    DBG(debug_S012, "S2[" << i << "]: " << *mod2_result);
            }

            mod0_runs.deallocate();
            mod1_runs.deallocate();
            mod2_runs.deallocate();

            // Prepare for merging

            mod0_result = new mod0_rm_type(mod0_runs.result(), less_mod0(), memsize / 3);
            mod1_result = new mod1_rm_type(mod1_runs.result(), less_quad_2nd(), memsize / 3);
            mod2_result = new mod2_rm_type(mod2_runs.result(), less_quint_2nd(), memsize / 3);
            // output: ISA_1,2 for next level

#if !LCP_CALC
            vmerge_sa = new merge_sa_type(*mod0_result, *mod1_result, *mod2_result);
#else // LCP_CALC
            vmerge_sa = new merge_sa_type(*mod0_result, *mod1_result, *mod2_result, lcp_base, lcp_rmq1, lcp_rmq1_ranks);
#endif // LCP_CALC

            result = *(*vmerge_sa);
        }

        const value_type & operator*() const
        {
            return result;
        }

        build_sa & operator++()
        {
            assert(vmerge_sa != 0 && !vmerge_sa->empty());

            ++(*vmerge_sa);
            if (!vmerge_sa->empty()) {
                result = *(*vmerge_sa);
            }
            else {		// cleaning up
                assert(vmerge_sa->empty());
                ready = true;

                assert(vmerge_sa != NULL);
                delete vmerge_sa; vmerge_sa = NULL;

                assert(mod0_result != NULL && mod1_result != NULL && mod2_result != NULL);
                delete mod0_result; mod0_result = NULL;
                delete mod1_result; mod1_result = NULL;
                delete mod2_result; mod2_result = NULL;
            }

            return *this;
        }

        bool empty() const
        {
            return ready;
        }
    };

#if LCP_CALC
    class build_lcp
    {
    public:

        typedef offset_type  value_type;

    private:

        std::auto_ptr<offset_array_type> lcp_base;

        offset_array_it_rg       lcp_baseiter;

        std::auto_ptr<RMQ>      lcp_rmq;

        offset_type             index;

        value_type              result;

        void calc_next_lcp()
        {
            if (empty()) return;

            result = *lcp_baseiter;
            DBG(debug_lcp_build, "lcp[" << index << "] = " << result);

            // add lcp of associated RMQ, if found.
            if (!lcp_rmq->empty())
            {
                if ((*lcp_rmq)->first == index)
                {
                    result += (*lcp_rmq)->second;
                    DBG(debug_lcp_build, "lcp[" << index << "] += " << (*lcp_rmq)->second << " -> " << result);
                    ++(*lcp_rmq);
                }
                assert((*lcp_rmq).empty() || (*lcp_rmq)->first > index);
            }
        }

    public:

        build_lcp(offset_array_type* _lcp_base, RMQ* _lcp_rmq, uint8_array_type& lcp_names)
            : lcp_base(_lcp_base),
              lcp_baseiter(lcp_base->begin(), lcp_base->end()),
              lcp_rmq(_lcp_rmq),
              index(0)
        {
            DBGMEM("computing RMQs on LCP_N");

            // some RMQ problem on LCP_names and prepare answers for iteration
            uint8_array_it_rg lcp_names_stream(lcp_names.begin(), lcp_names.end());
            lcp_rmq->compute(lcp_names_stream);

            DBGMEM("done computing RMQs on LCP_N");
        }

        void start_output(size_t memsize)
        {
            lcp_rmq->start_output(memsize);
            calc_next_lcp();
        }

        const value_type & operator*() const
        {
            return result;
        }

        build_lcp& operator++ ()
        {
            ++lcp_baseiter, ++index;
            calc_next_lcp();
            return *this;
        }

        bool empty() const
        {
            return lcp_baseiter.empty();
        }
    };
#endif // LCP_CALC

    /***********************************************************************************************
     * The skew algorithm                                                                          *
     ***********************************************************************************************/

    /**
     * The skew algorithm (in its pure form)
     * @param Input type of the input pipe.
     * @param block_size block size of an external sorter we use.
     * @param quiet if true we omit informations about the computation process.
     */

    template <class Input, const int block_size = ::block_size>
    class algorithm
    {
#if LCP_CALC
        static const int memdivA        = 5;        // memory divisor of loop
#else
        static const int memdivA        = 4;
#endif

    public:
        typedef offset_type value_type;

        typedef typename Input::value_type       alphabet_type;

    private:
        // All needed types

        // *** Types for alphabet_type input - make input stream compatible to loop functors.

        typedef tuples::counter <offset_type> counter_stream_type;
        // (t_i) -> (i,t_i)
        typedef make_pairs <counter_stream_type, Input> make_pairs_input_type;

        // All needed types for recursive name construction loop.

        // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
        typedef make_quads <make_pairs_input_type, offset_type, 1> make_quads_input_type;
        // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
        typedef extract_mod12 <make_quads_input_type> mod12_quads_input_type;

        // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
        typedef typename stream::sort<mod12_quads_input_type, less_quad<offset_type>, block_size> sort_mod12_input_type;

        // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
        typedef naming <sort_mod12_input_type> naming_input_type;

        // (i,n_i) sort by (i mod 3, i div 3) to create recursive text
        typedef typename stream::runs_creator <naming_input_type, less_skew, block_size> rc_order_type;
        typedef typename rc_order_type::sorted_runs_type sr_order_type;
        typedef typename stream::runs_merger <sr_order_type, less_skew> rm_order_type;

        // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with recursive text
        typedef make_quads <rm_order_type, offset_type, 1> make_quads_loop_type;

        // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
        typedef extract_mod12 <make_quads_loop_type> mod12_quads_loop_type;

        // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
        typedef typename stream::sort <mod12_quads_loop_type, less_quad<offset_type>, block_size> sort_loop_type;

        // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
        typedef naming <sort_loop_type> naming_loop_type;

        typedef typename stream::runs_creator <naming_loop_type, less_skew, block_size> rc_loop_type;

        // Types for sub problem loop

        // build suffix array in base case from saved ISA (unique names recursion string)
        typedef build_sa <offset_array_it_rg, offset_array_it_rg, offset_array_it_rg, block_size, ram_use * 3/memdivA> build_sa_type;

        // Inverse sa types

        typedef make_pairs < counter_stream_type, build_sa_type > prepare_isa_type;

        // Seperate mod1 from mod2 tuples

        typedef typename stream::use_push < skew_pair_type > isa1_push_type;
        typedef typename stream::runs_creator < isa1_push_type, less_pair_2nd, block_size, RC > isa1_rc_type;
        typedef typename isa1_rc_type::sorted_runs_type isa1_sr_type;
        typedef typename stream::runs_merger < isa1_sr_type, less_pair_2nd > isa1_rm_type;

        typedef typename stream::use_push < skew_pair_type > isa2_push_type;
        typedef typename stream::runs_creator < isa2_push_type, less_pair_2nd, block_size, RC > isa2_rc_type;
        typedef typename isa2_rc_type::sorted_runs_type isa2_sr_type;
        typedef typename stream::runs_merger < isa2_sr_type, less_pair_2nd > isa2_rm_type;

        // Loop types

        typedef stream::choose <isa1_rm_type,1> isa1_first_type;
        typedef stream::choose <isa2_rm_type,1> isa2_first_type;
        typedef build_sa < offset_array_it_rg, isa1_first_type, isa2_first_type, block_size, ram_use * 3/memdivA> build_sa_loop_type;
        typedef make_pairs < counter_stream_type, build_sa_loop_type > prepare_isa_loop_type;

        typedef typename stream::use_push < skew_pair_type > isa_loop_push_type;
        typedef typename stream::runs_creator < isa_loop_push_type, less_pair_2nd, block_size, RC > isa_loop_rc_type;

        // Variables

#if !LCP_CALC
        build_sa_loop_type *out_sa;
#else
        offset_array_it_rg* out_sa;
        offset_array_type* out_sa_array;

        build_lcp* out_lcp;
#endif

        bool finished;

        template <class Stream>
        static inline void print(Stream& s)
        {
            while (!s.empty())
            {
                std::cout << *s << "\n";
                ++s;
            }
        }

    protected:

#if LCP_CALC
        static build_lcp* expand_build_lcp(build_lcp* prev_lcp, offset_array_type* prev_sa, isa1_sr_type &isa1_sr, isa2_sr_type& isa2_sr,
                                           uint8_array_type& lcp_names,
                                           uint8_array_type* lcp_base8, RMQ* lcp_rmq1_ptr, skew_pair_array_type* lcp_rmq1_ranks  )
        {
            DBGMEM("computing RMQs on LCP_R");

            // compute RMQs on recursive LCP array
            RMQ& lcp_rmq1 = *lcp_rmq1_ptr;
            prev_lcp->start_output(ram_use / 2);
            lcp_rmq1.compute(*prev_lcp);
            delete prev_lcp;

            // read RMQ results: create expanded lcp_next and RMQs on the LCP_lexnames array

            skew_pair_array_it_rg rmq1_ranks(lcp_rmq1_ranks->begin(), lcp_rmq1_ranks->end());
            offset_array_type* lcp_next = new offset_array_type;

            lcp_rmq1.start_output(ram_use / 3);
            RMQ* lcp_rmq2 = new RMQ(ram_use / 3);

            stxxl::sorter<skew_triple_type,less_triple_1st,block_size> SAqueries (less_triple_1st(), ram_use / 3);

            // memsize: input stream, rmq1, output: rmq2 and SAqueries
            {
                uint8_array_it_rg lcp_base8_stream(lcp_base8->begin(), lcp_base8->end());

                for (offset_type index = 0; !lcp_base8_stream.empty(); ++index, ++lcp_base8_stream)
                {
                    offset_type lcp = *lcp_base8_stream;

                    if ( !lcp_rmq1.empty() && lcp_rmq1->first == index ) // RMQ available?
                    {
                        assert( !rmq1_ranks.empty() );

                        if (lcp_rmq1->second == offset_type(0)) // RMQ result was 0 -> no lookup on SA/ISA needed
                        {
                            DBG(debug_lcp_expand, "next_lcp[" << index << "] = " << lcp << " + direct RMQ on ranks " << *rmq1_ranks);
                            assert(rmq1_ranks->first <= rmq1_ranks->second);

                            lcp_rmq2->query( rmq1_ranks->first+1, rmq1_ranks->second, index );
                        }
                        else
                        {
                            lcp += 3 * lcp_rmq1->second;        // add 3 for each matching triple

                            DBG(debug_lcp_expand, "next_lcp[" << index << "] = " << lcp << " + far RMQ on ISA[SA[" << *rmq1_ranks << "] + " << lcp_rmq1->second << "]");

                            SAqueries.push( skew_triple_type(rmq1_ranks->first, lcp_rmq1->second, index) );
                            SAqueries.push( skew_triple_type(rmq1_ranks->second, lcp_rmq1->second, index) );
                        }

                        ++lcp_rmq1, ++rmq1_ranks;
                    }
                    else
                    {
                        DBG(debug_lcp_expand, "next_lcp[" << index << "] = " << lcp);
                    }

                    lcp_next->push_back(lcp);
                }

                delete lcp_base8;
                delete lcp_rmq1_ptr;
                delete lcp_rmq1_ranks;
            }

            // *** Answer queries on SA

            // following processes run with lcp_rmq2 keeping it's memory!
            SAqueries.sort(ram_use / 3);

            if (SAqueries.size())
            {
                // *** Scan saved SA_R and answer queries for SA positions

                stxxl::sorter<skew_pair_type,less_pair_1st,block_size> ISAqueries (less_pair_1st(), ram_use / 3);

                offset_array_it_rg sa(prev_sa->begin(), prev_sa->end());

                for (offset_type index = 0; !sa.empty() && !SAqueries.empty(); ++index, ++sa)
                {
                    DBG(debug_lcp_expand, "SA[" << index << "] = " << *sa);
                    while ( !SAqueries.empty() && SAqueries->first == index )
                    {
                        DBG(debug_lcp_expand, "Matching SA[" << index << "] query with " << *sa << " + " << SAqueries->second << " for target " << SAqueries->third);
                        ISAqueries.push( skew_pair_type(*sa + SAqueries->second, SAqueries->third) );
                        ++SAqueries;
                    }
                }

                assert( SAqueries.empty() );
                SAqueries.finish_clear();
                delete prev_sa;

                // *** answer queries on ISA by reading SA again

                stxxl::sorter<skew_pair_type,less_pair_1st,block_size> ISAanswers (less_pair_1st(), ram_use * 2/9);
                ISAqueries.sort(ram_use * 2/9);

                for (isa1_rm_type isa1_rm(isa1_sr, less_pair_2nd(), ram_use * 2/9);
                     !isa1_rm.empty() && !ISAqueries.empty(); ++isa1_rm)
                {
                    DBG(debug_lcp_expand, "ISA1: " << *isa1_rm);

                    while ( !ISAqueries.empty() && ISAqueries->first == isa1_rm->second )
                    {
                        DBG(debug_lcp_expand, "Matching ISA[" << ISAqueries->first << "] query for target " << ISAqueries->second << " with " << isa1_rm->first);
                        ISAanswers.push( skew_pair_type(ISAqueries->second, isa1_rm->first) );
                        ++ISAqueries;
                    }
                }

                for (isa2_rm_type isa2_rm(isa2_sr, less_pair_2nd(), ram_use * 2/9);
                     !isa2_rm.empty() && !ISAqueries.empty(); ++isa2_rm)
                {
                    DBG(debug_lcp_expand, "ISA2: " << *isa2_rm);

                    while ( !ISAqueries.empty() && ISAqueries->first == isa2_rm->second )
                    {
                        DBG(debug_lcp_expand, "Matching ISA[" << ISAqueries->first << "] query for target " << ISAqueries->second << " with " << isa2_rm->first);
                        ISAanswers.push( skew_pair_type(ISAqueries->second, isa2_rm->first) );
                        ++ISAqueries;
                    }
                }

                while ( !ISAqueries.empty() )
                {
                    DBG(debug_lcp_expand, "Out-of-bounds ISA[" << ISAqueries->first << "] query for target " << ISAqueries->second << " with ");
                    ++ISAqueries;
                }

                assert( ISAqueries.empty() );
                ISAqueries.finish_clear();

                // *** collect answers from ISA and prepare RMQ on LCP_N

                ISAanswers.sort(ram_use / 3);

                while ( !ISAanswers.empty() )
                {
                    skew_pair_type p1 = *ISAanswers; ++ISAanswers;

                    if ( ISAanswers.empty() || ISAanswers->first != p1.first )
                    {
                        DBG(debug_lcp_expand, "Discarding out-of-bounds RMQ(" << p1.second << "+1,?) -> target " << p1.first);
                        continue;
                    }

                    skew_pair_type p2 = *ISAanswers; ++ISAanswers;
                    assert( p1.first == p2.first );

                    if (p1.second > p2.second) std::swap(p1.second, p2.second);

                    DBG(debug_lcp_expand, "Matching for RMQ(" << p1.second << "+1," << p2.second << ") -> target " << p1.first);

                    lcp_rmq2->query( p1.second+1, p2.second, p1.first );
                }
            }

            DBGMEM("done");

            // Prepare LCP calculation from RMQs and LCP_lexnames

            // build_lcp takes ownership of lcp_next and lcp_rmq2
            return new build_lcp(lcp_next, lcp_rmq2, lcp_names);
        }
#endif // LCP_CALC

    public:
        algorithm(Input & data_in)
            : finished(false)
        {
            // store intermediate STXXL-vectors
            typename std::vector<offset_array_type *> vectors;

#if LCP_CALC
            typename std::vector<uint8_array_type *> lcp_name_vectors;
#endif // LCP_CALC

            /**
             * Steps for building 2/3 arrays for the subproblem:
             * 1) Build quads <index,*,*,*>. The last three components are the characters, beginning at index.
             *    Store the incomming elements into a vector. Do not forget special cases.
             * 2) Drop the elements at positions (index modulo 3 = 0).
             * 3) Sort the elements in respect to the last three components (lex. order).
             * 4) Scann throught the elements and give to each a discarding compatible names and check for uniqueness.
             * 5) Sort the elements in respect to the index.
             * 6) Do the 1-5 within a loop, until the names are unique.
             * -mehnert
             */

            vectors.push_back(new offset_array_type);           // to store current S
#if LCP_CALC
            lcp_name_vectors.push_back(new uint8_array_type);       // to store lcp array of lexname
#endif // LCP_CALC

            DBG(debug, "Lexnaming input in depth 0");

            // (t_i) -> (i,t_i)
            counter_stream_type dummy;
            make_pairs_input_type pairs_input(dummy, data_in);

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
            make_quads_input_type quads_input(pairs_input, *vectors.back());

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
            mod12_quads_input_type mod12_quads_input(quads_input);

            // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
            sort_mod12_input_type sort_mod12_input(mod12_quads_input, less_quad<offset_type>(), ram_use, ram_use / 2);

            // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
            bool unique = false;        // ist the current quad array unique ?
#if !LCP_CALC
            naming_input_type names_input(sort_mod12_input, unique);
#else // LCP_CALC
            naming_input_type names_input(sort_mod12_input, *lcp_name_vectors.back(), unique);
#endif // LCP_CALC

            rc_order_type *resume_order_rc = new rc_order_type(names_input, less_skew(), ram_use / 2);
            rm_order_type *resume_order_rm = new rm_order_type(resume_order_rc->result(), less_skew(), ram_use / 2);
            rc_loop_type *loop_order_rc = NULL;

            DBGMEM("First naming done");

            /** loop for building the skew data */

            size_t depth = 0;

            while (!unique)
            {
                DBG(debug, "Lexnaming " << resume_order_rm->size() << " triples in depth " << ++depth);

                vectors.push_back(new offset_array_type);           // to store current S
#if LCP_CALC
                lcp_name_vectors.push_back(new uint8_array_type);   // to store lcp array of lexname
#endif // LCP_CALC

                // If we cross the border from mod1 tuples to mod2 tuples give the last two tuples
                // the right name !!!!!!!!!!!! -mehnert

                // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with recursive text
                make_quads_loop_type quads_loop(*resume_order_rm, *vectors.back());

                // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
                mod12_quads_loop_type mod12_quads_loop(quads_loop);

                // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
                sort_loop_type sort_quads_loop(mod12_quads_loop, less_quad<offset_type>(), ram_use / 2);

                if (resume_order_rc != NULL) {
                    delete resume_order_rc; resume_order_rc = NULL;
                }
                if (loop_order_rc != NULL) {
                    delete loop_order_rc; loop_order_rc = NULL;
                }

                // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
#if !LCP_CALC
                naming_loop_type name_quads_loop(sort_quads_loop, unique);
#else // LCP_CALC
                naming_loop_type name_quads_loop(sort_quads_loop, *lcp_name_vectors.back(), unique);
#endif // LCP_CALC

                assert(resume_order_rm != NULL);
                delete resume_order_rm; resume_order_rm = NULL;

                loop_order_rc = new rc_loop_type(name_quads_loop, less_skew(), ram_use / 2);
                resume_order_rm = new rm_order_type(loop_order_rc->result(), less_skew(), ram_use / 2);

                DBGMEM("Next naming level done");

                if (unique) // Push the subsuffixarray on the list
                {
                    vectors.push_back( new offset_array_type(resume_order_rm->size()) );	// to store last S

                    stream::choose<rm_order_type,2> merger(*resume_order_rm);
                    stream::materialize(merger, vectors.back()->begin(), vectors.back()->end());
                }

            } // while (!unique)

            DBG(debug, "Found unique lexnames.");
            g_statscache >> "maxdepth" << depth;

            // clean garbage

            if (loop_order_rc != NULL) {
                delete loop_order_rc; loop_order_rc = NULL;
            }
            if (resume_order_rc != NULL) {
                delete resume_order_rc; resume_order_rc = NULL;
            }
            if (resume_order_rm != NULL) {
                delete resume_order_rm; resume_order_rm = NULL;
            }

            DBGMEM("Unique names found.");

#if SKEW_DEBUG
            for (size_t i = 0; i < vectors.size(); ++i)
            {
                offset_array_type & test = *vectors[i];
                for (size_t e = 0; e < test.size(); ++e)
                    std::cout << test[e] << " ";
                if ((test.size() & 1) == 1)
                    std::cout << " Break: " << (test.size() >> 1) + 1;
                else
                    std::cout << " Break: " << (test.size() >> 1);
                std::cout << " Nextbreak: " << subp_size(test.size());
                std::cout << " Size: " << test.size();
                std::cout << std::endl;
            }
#endif

            // *************************************************************************************
            // *** Base case suffix array creation

            /**
             * Simple suffix array creation in base case with unique names
             * Mod0 quint: <index, character, character+1, ~SA12+1, ~SA12+2>
             * Mod1 quad : <index, ~SA12 name, character, ~SA12+1 name>
             * Mod2 quint: <index, ~SA12 name, character, character+1, ~SA12+2 name>
             * -mehnert
             */

            assert(vectors.size() > 2);	// We do not want all names finished at this time ...

            isa1_sr_type isa1_sr;
            isa2_sr_type isa2_sr;
#if LCP_CALC
            offset_array_type* prev_sa;
            build_lcp* prev_lcp;
#endif // LCP_CALC

            { // base case suffix array

                DBG(debug, "Merging from base case in depth " << depth);

                offset_array_type& base_source = *vectors.back(); // contains ISA of base case -tb
                offset_array_type& prev_array = *vectors[vectors.size() - 2];

                // calculate start of mod2 area

                typename offset_array_type::iterator start_mod2 = base_source.begin();
                start_mod2 = start_mod2 + ((prev_array.size() / 3) + ((prev_array.size() % 3) == 2));

                // Special case, if the injected tuple had affected the next subproblem

                size_t special = (subp_size(prev_array.size()) != base_source.size()) ? 1 : 0;

                // Generate streams from vector

                offset_array_it_rg source(prev_array.begin(), prev_array.end());
                offset_array_it_rg mod1(base_source.begin(), start_mod2);
                offset_array_it_rg mod2(start_mod2 + special, base_source.end());

                // Build suffix array
#if !LCP_CALC
                build_sa_type get_sa(source, mod1, mod2, prev_array.size());
#else // LCP_CALC
                uint8_array_type lcp_base8;
                RMQ* lcp_rmq = new RMQ(ram_use / memdivA);
                skew_pair_array_type lcp_rmq1_ranks;

                build_sa_type get_sa(source, mod1, mod2, prev_array.size(), lcp_base8, *lcp_rmq, lcp_rmq1_ranks);
#endif // LCP_CALC

                // Prepare inverse suffix array

                counter_stream_type isa_index;
                prepare_isa_type isa_tuples(isa_index, get_sa);

                // Seperate isa

                isa1_rc_type isa1_rc(less_pair_2nd(), ram_use / memdivA/2);
                isa2_rc_type isa2_rc(less_pair_2nd(), ram_use / memdivA/2);

                offset_array_type *last_array = vectors[vectors.size() - 3];
                special = (prev_array.size() != subp_size(last_array->size())); // Special case (see above)
                offset_type mod2_pos = (subp_size(last_array->size()) >> 1) + (subp_size(last_array->size()) & 1) + special;
                DBG(0, "mod2_pos: " << mod2_pos << " special:" << special << " minus2: " << last_array->size() << " prev2: " << prev_array.size());

#if LCP_CALC
                prev_sa = new offset_array_type;
#endif // LCP_CALC

                while (!isa_tuples.empty()) { // (first,second) = (idx,sa[idx])
                    const skew_pair_type& tmp = *isa_tuples; // (first,second) = (i,sa[i])

#if LCP_CALC
                    prev_sa->push_back(tmp.second);
#endif // LCP_CALC

                    if (tmp.second < mod2_pos) {
#if !LCP_CALC
                        if (tmp.second + special < mod2_pos) // else: special sentinel tuple is dropped -tb
#endif
                            isa1_rc.push(tmp);
                    }
                    else {
                        isa2_rc.push(tmp);
                    }
                    ++isa_tuples;
                }

                isa1_rc.deallocate(); // deallocate run_creator's buffers
                isa2_rc.deallocate();

#if LCP_CALC
                // Prepare LCP calculation from RMQs and LCP_lexnames

                assert(base_source.size() == lcp_name_vectors.back()->size());

                // transform uint8_t lcps to offset_type
                offset_array_type* lcp_base = new offset_array_type(lcp_base8.size());
                uint8_array_it_rg lcp_base8_iter(lcp_base8.begin(), lcp_base8.end());
                stream::materialize(lcp_base8_iter, lcp_base->begin(), lcp_base->end());

                // build_lcp takes ownership of lcp_base and lcp_rmq
                prev_lcp = new build_lcp(lcp_base, lcp_rmq, *lcp_name_vectors.back());

                delete lcp_name_vectors.back(); // in build_lcp constructor the lcp_name RMQs are
                                                // answered, so we can delete LCP_names here
                lcp_name_vectors.pop_back();
#endif // LCP_CALC

                // Save ISA runs for next subproblem
                isa1_sr = isa1_rc.result();
                isa2_sr = isa2_rc.result();

                // debug ISA
                if (debug_isa)
                {
                    for (isa1_rm_type isa1_rm(isa1_sr, less_pair_2nd(), ram_use / 8);
                         !isa1_rm.empty(); ++isa1_rm)
                        DBG(debug_isa, "ISA1: " << *isa1_rm);

                    for (isa2_rm_type isa2_rm(isa2_sr, less_pair_2nd(), ram_use / 8);
                         !isa2_rm.empty(); ++isa2_rm)
                        DBG(debug_isa, "ISA2: " << *isa2_rm);
                }

            } // end base case suffix array

            // *************************************************************************************
            // *** Loop head

            // delete last subproblem
            delete vectors.back();
            vectors.pop_back();

#if SKEW_DEBUG
            for(size_t i = 0; i < vectors.size(); ++i) {
                offset_array_type& test = *vectors[i];
                for(size_t e = 0; e < test.size(); ++e)
                    std::cout << test[e] << " ";
                std::cout.flush();
                if ((test.size() & 1) == 1)
                    std::cout << " Break: " << (test.size() >> 1) + 1;
                else
                    std::cout << " Break: " << (test.size() >> 1);
                std::cout << " Nextbreak: " << subp_size(test.size());
                std::cout << " Size: " << test.size();
                std::cout << std::endl;
            }
#endif

            offset_array_type *prev_array = vectors[vectors.size() - 2];

            while (vectors.size() > 2)
            {
                DBG(debug, "Exiting recursion depth " << --depth);

                DBGMEM("Loop start");

                isa1_rm_type isa1_rm(isa1_sr, less_pair_2nd(), ram_use / memdivA/2);
                isa2_rm_type isa2_rm(isa2_sr, less_pair_2nd(), ram_use / memdivA/2);

                isa1_first_type isa1_first(isa1_rm);
                isa2_first_type isa2_first(isa2_rm);

                prev_array = vectors[vectors.size() - 2];

                offset_array_it_rg source(prev_array->begin(), prev_array->end());

                // The runs merger will be emptied after suffix array construction in the following line
#if !LCP_CALC

                build_sa_loop_type get_sa(source, isa1_first, isa2_first, prev_array->size());

#else // LCP_CALC

                uint8_array_type* lcp_base8 = new uint8_array_type;
                RMQ* lcp_rmq1 = new RMQ(ram_use / memdivA);
                skew_pair_array_type* lcp_rmq1_ranks = new skew_pair_array_type;

                build_sa_loop_type get_sa(source, isa1_first, isa2_first, prev_array->size(),
                                          *lcp_base8, *lcp_rmq1, *lcp_rmq1_ranks);

                offset_array_type* suffixarray = new offset_array_type();

#endif // LCP_CALC

                isa1_rm.deallocate(); isa2_rm.deallocate();
#if !LCP_CALC
                isa1_sr->clear(); isa2_sr->clear();       // for LCP calculation, prev_isa is still needed.
#endif

                // *** Read subsuffixarray and generate ISA from it. This also generates the LCP requests

                counter_stream_type isa_loop_index;
                prepare_isa_loop_type isa_loop_tuples(isa_loop_index, get_sa);

                // build runs creator

                isa_loop_rc_type isa1_loop_rc(less_pair_2nd(), ram_use / memdivA/2);
                isa_loop_rc_type isa2_loop_rc(less_pair_2nd(), ram_use / memdivA/2);

                offset_array_type* last_array = vectors[vectors.size() - 3];

                size_t special = (prev_array->size() != subp_size(last_array->size())) ? 1 : 0;
                offset_type mod2_loop_pos = (subp_size(last_array->size()) >> 1) + (subp_size(last_array->size()) & 1) + special;

                DBG(0, "Mod2_Pos: " << mod2_loop_pos << " Special:" << special << " Minus2: " << last_array->size() << " Prev2: " << prev_array->size() << " Vsize:"  << vectors.size() ); // -mehnert

                while (!isa_loop_tuples.empty()) {
                    const skew_pair_type& tmp = *isa_loop_tuples; // (first,second) = (i,sa[i])
#if LCP_CALC
                    suffixarray->push_back(tmp.second);
#endif // LCP_CALC
                    if (tmp.second < mod2_loop_pos) {
#if !LCP_CALC
                        if (tmp.second + special < mod2_loop_pos) // else: special sentinel tuple is dropped -tb
#endif
                            isa1_loop_rc.push(tmp);
                    }
                    else {
                        isa2_loop_rc.push(tmp);
                    }
                    ++isa_loop_tuples;
                }

                isa1_rm.deallocate();
                isa2_rm.deallocate();
                isa1_loop_rc.deallocate(); // deallocate run_creator's buffers
                isa2_loop_rc.deallocate();

#if LCP_CALC
                // Prepare LCP calculation from RMQs and LCP_lexnames

                assert( lcp_base8->size() == suffixarray->size() );
                assert( vectors.back()->size() == lcp_name_vectors.back()->size() );

                prev_lcp = expand_build_lcp(prev_lcp, prev_sa, isa1_sr, isa2_sr,
                                            *lcp_name_vectors.back(),
                                            lcp_base8, lcp_rmq1, lcp_rmq1_ranks);

                delete lcp_name_vectors.back(); // in build_lcp constructor the lcp_name RMQs are
                                                // answered, so we can delete LCP_names here
                lcp_name_vectors.pop_back();
                prev_sa = suffixarray;

#endif // LCP_CALC

                // save ISA runs for next level -- frees previous result (smart pointer)
                isa1_sr = isa1_loop_rc.result();
                isa2_sr = isa2_loop_rc.result();

                // debug ISA
                if (debug_isa)
                {
                    for (isa1_rm_type isa1_rm(isa1_sr, less_pair_2nd(), ram_use / 8);
                         !isa1_rm.empty(); ++isa1_rm)
                        DBG(debug_isa, "ISA1: " << *isa1_rm);

                    for (isa2_rm_type isa2_rm(isa2_sr, less_pair_2nd(), ram_use / 8);
                         !isa2_rm.empty(); ++isa2_rm)
                        DBG(debug_isa, "ISA2: " << *isa2_rm);
                }

                // delete last subproblem
                delete vectors.back();
                vectors.pop_back();
            }

            // *************************************************************************************
            // *** Final stage returning top suffix array

            DBG(debug, "Merging final suffix array.");

            {
                isa1_rm_type isa1_rm(isa1_sr, less_pair_2nd(), ram_use / memdivA/2);
                isa2_rm_type isa2_rm(isa2_sr, less_pair_2nd(), ram_use / memdivA/2);

                isa1_first_type isa1_first(isa1_rm);
                isa2_first_type isa2_first(isa2_rm);

                prev_array = vectors[vectors.size() - 2]; // this is the plain input string

                offset_array_it_rg source(prev_array->begin(), prev_array->end());

                DBG(debug, "Final build_sa: source " << prev_array->size() << " - mod1 " << isa1_rm.size() << " - mod2 " << isa2_rm.size());

#if !LCP_CALC

                out_sa = new build_sa_loop_type(source, isa1_first, isa2_first, prev_array->size());

#else // LCP_CALC

                uint8_array_type* lcp_base8 = new uint8_array_type;
                RMQ* lcp_rmq1 = new RMQ(ram_use / memdivA);
                skew_pair_array_type* lcp_rmq1_ranks = new skew_pair_array_type;

                out_sa_array = new offset_array_type(prev_array->size());

                {
                    build_sa_loop_type get_sa(source, isa1_first, isa2_first, prev_array->size(),
                                              *lcp_base8, *lcp_rmq1, *lcp_rmq1_ranks);

                    // read suffixarray into vector, this also creates the lcp_base all LCP-RMQs on
                    // the recursive LCP array.
                    stream::materialize(get_sa, out_sa_array->begin(), out_sa_array->end());
                }

                isa1_rm.deallocate();
                isa2_rm.deallocate();

                assert( lcp_base8->size() == out_sa_array->size() );
                assert( vectors.back()->size() == lcp_name_vectors.back()->size() );

                prev_lcp = expand_build_lcp(prev_lcp, prev_sa, isa1_sr, isa2_sr,
                                            *lcp_name_vectors.back(),
                                            lcp_base8, lcp_rmq1, lcp_rmq1_ranks);
                prev_lcp->start_output(ram_use / 4);

                delete lcp_name_vectors.back(); // in build_lcp constructor the lcp_name RMQs are
                                                // answered, so we can delete LCP_names here
                lcp_name_vectors.pop_back();

                out_sa = new offset_array_it_rg(out_sa_array->begin(), out_sa_array->end());
                out_lcp = prev_lcp;
#endif
            }

            DBG(debug, "Done with final merge.");

            DBGMEM("done");
        }

        const value_type & operator*() const
        {
            return *(*out_sa);
        }

#if LCP_CALC
        const value_type& lcp() const
        {
            return *(*out_lcp);
        }
#endif // LCP_CALC

        algorithm & operator++()
        {
            assert(!out_sa->empty());

            ++(*out_sa);
#if LCP_CALC
            ++(*out_lcp);
#endif // LCP_CALC

            if ((out_sa != NULL) && (out_sa->empty())) {
                finished = true;
                delete out_sa; out_sa = NULL;
#if LCP_CALC
                delete out_sa_array; out_sa_array = NULL;
                delete out_lcp; out_lcp = NULL;
#endif // LCP_CALC
            }
            return *this;
        }

        bool empty() const
        {
            return finished;
        }
    };
};

template <typename StringContainer, typename SuffixArrayContainer>
struct SACA
{
    std::string name()
    {
#if !LCP_CALC
        return "skew3";
#else // LCP_CALC
        return "skew3-lcp";
#endif // LCP_CALC
    }

    void prepare(const StringContainer& string, SuffixArrayContainer& suffixarray, unsigned int K)
    {
        stxxl::STXXL_UNUSED(string);
        stxxl::STXXL_UNUSED(suffixarray);
        stxxl::STXXL_UNUSED(K);
    }

    void run(const StringContainer& string, SuffixArrayContainer& suffixarray, unsigned int K)
    {
        stxxl::STXXL_UNUSED(K);

        typedef typename StringContainer::value_type alphabet_type;
        typedef typename SuffixArrayContainer::value_type offset_type;

        typedef stxxl::stream::vector_iterator2stream< typename StringContainer::const_iterator > input_type;
        typedef typename skew<offset_type>::template algorithm<input_type> saca_type;

        input_type input = input_type(string.begin(), string.end());
        saca_type saca(input);
        stxxl::stream::materialize(saca, suffixarray.begin(), suffixarray.end());
    }

#if LCP_CALC
    template <typename LCPArrayContainer>
    void run_lcp(const StringContainer& string, SuffixArrayContainer& suffixarray, LCPArrayContainer& lcparray, unsigned int /* K */)
    {
        typedef typename StringContainer::value_type alphabet_type;
        typedef typename SuffixArrayContainer::value_type offset_type;

        typedef stxxl::stream::vector_iterator2stream< typename StringContainer::const_iterator > input_type;
        typedef typename skew<offset_type>::template algorithm<input_type> saca_type;

        input_type input = input_type(string.begin(), string.end());
        saca_type saca(input);

        stxxl::vector_bufwriter<SuffixArrayContainer> sa_writer (suffixarray.begin(), suffixarray.end());
        stxxl::vector_bufwriter<LCPArrayContainer> lcp_writer (lcparray.begin(), lcparray.end());

        for (size_t i = 0; !saca.empty(); ++saca, ++i)
        {
            //std::cout << "SA[" << i << "] = " << operator*() << " - LCP[" << i << "] = " << lcp() << "\n";
            sa_writer << *saca;
            lcp_writer << saca.lcp();
        }
    }
#endif // LCP_CALC
};

} // namespace skew
