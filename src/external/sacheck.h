/******************************************************************************
 * src/external/sacheck.h
 *
 * Algorithm to check a suffix array. Loosely based on the ideas of
 * Kaerkkaeinen und Burghardt.
 *
 * Reimplemented using triples by Timo Bingmann (2012)
 * Originally implemented in STXXL by Jens Mehnert (2004)
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

/* Definitions:
 * T   = The original text, from which the suffix array was build
 * SA  = Perhaps suffixarray from T
 * ISA = The inverse of SA
 */

namespace sacheck
{
/**
 * Build tuples with index: (SA[i]) -> (SA[i],i)
 */
template < class InputSA >
class tuplize_sa_index
{
    typedef typename InputSA::value_type offset_type;

    typedef tuples::pair<offset_type, offset_type>	tuple_type;

    InputSA&		m_input;

    /// current position in input
    offset_type         m_counter;

    bool		m_finished;

    tuple_type		m_result;

public:

    typedef tuple_type	value_type;

    tuplize_sa_index(InputSA& input)
        : m_input(input), m_counter(0), m_finished(false)
    {
        assert(!m_input.empty());
        m_result = value_type(m_counter, *m_input);
    }

    const value_type & operator*() const
    {
        return m_result;
    }

    tuplize_sa_index & operator++()
    {
        assert(!m_input.empty());

        ++m_input;
        ++m_counter;

        if (!m_input.empty())
            m_result = value_type(m_counter, *m_input);
        else
            m_finished = true;

        return *this;
    }

    bool empty()
    {
        return m_finished;
    }
};

template <typename InputT, typename InputSA, typename SizeType>
bool check_sa(InputT& inputT, InputSA& inputSA)
{
    using namespace stxxl;
    using namespace tuples;

    typedef typename InputT::value_type         alphabet_type;
    typedef typename InputSA::value_type        offset_type;
    typedef SizeType                            size_type;

    typedef tuples::pair<offset_type, offset_type>	          pair_type;
    typedef tuples::triple<offset_type,alphabet_type,offset_type> triple_type;

    // *** Pipeline Declaration ***

    // build (SA[i], i) and sort to (i, ISA[i])

    typedef pair_less2nd<pair_type>    pair_less_type;

    typedef tuplize_sa_index<InputSA>   tuplize_sa_index_type;
    typedef typename stream::sort<tuplize_sa_index_type, pair_less_type, block_size> build_isa_type;

    tuplize_sa_index_type	tuplize_sa_index(inputSA);
    build_isa_type		build_isa(tuplize_sa_index, pair_less_type(), ram_use / 3);

    // build (ISA[i], T[i], ISA[i+1]) and sort to (i, T[SA[i]], ISA[SA[i]+1])

    typedef triple_less1st<triple_type>	triple_less_type;	// comparison relation

    typedef typename stream::use_push<triple_type> triple_push_type;	// indicator use push()

    typedef typename stream::runs_creator<triple_push_type, triple_less_type, block_size, RC>	triple_rc_type;
    typedef typename stream::runs_merger<typename triple_rc_type::sorted_runs_type, triple_less_type> triple_rm_type;

    triple_rc_type	triple_rc(triple_less_type(), ram_use / 3);

    // *** Process ***

    // ********************************************************************************
    // loop 1: read ISA and check for a permutation. Simultaneously create runs
    // of triples by iterating ISA and T.

    size_type totalSize;

    {
        offset_type prev_isa = (*build_isa).first;
        offset_type counter = 0;

        while (!build_isa.empty())
        {
            if ((*build_isa).second != counter) {
                std::cout << "Error: suffix array is not a permutation of 0..n-1.\n";
                return false;
            }

            ++counter;
            ++build_isa;		// ISA is one in front of T

            if (!build_isa.empty())
            {
                //std::cout << "(" << prev_isa << "," << *inputT << "," << (*build_isa).first << ")" << std::endl;

                triple_rc.push( triple_type(prev_isa, *inputT, (*build_isa).first) );

                prev_isa = (*build_isa).first;
            }

            ++inputT;
        }

        totalSize = counter;
    }

    if (totalSize == 1) return true;

    // ********************************************************************************
    // loop 2: read triples (i,T[SA[i]],ISA[SA[i]+1]) and check for correct ordering.

    triple_rm_type	triple_rm(triple_rc.result(), triple_less_type(), ram_use / 3);

    {
        triple_type	prev_triple = *triple_rm;
        size_type	counter = 0;

        ++triple_rm;

        while (!triple_rm.empty())
        {
            const triple_type& this_triple = *triple_rm;

            //(*triple_rm).to_cout();

            if (prev_triple.second > this_triple.second)
            {
                // simple check of first character of suffix
                std::cout << "Error: suffix array position " << counter  << " ordered incorrectly.\n";
                return false;
            }
            else if (prev_triple.second == this_triple.second)
            {
                if ( this_triple.third == (offset_type)totalSize ) {
                    // last suffix of string must be first among those with same
                    // first character
                    std::cout << "Error: suffix array position " << counter << " ordered incorrectly.\n";
                    return false;
                }
                if ( prev_triple.third != (offset_type)totalSize &&
                     prev_triple.third > this_triple.third )
                {
                    // positions SA[i] and SA[i-1] has same first character but
                    // their suffixes are ordered incorrectly: the suffix
                    // position of SA[i] is given by ISA[SA[i]]
                    std::cout << "Error: suffix array position " << counter << " ordered incorrectly.\n";
                    return false;
                }
            }

            prev_triple = this_triple;

            ++triple_rm;
            ++counter;
        }
    }

    return true;
}

template <typename InputT, typename InputSA>
bool check_sa_vectors(InputT& inputT, InputSA& inputSA)
{
    using namespace tuples;

    typedef it_rg<typename InputT::iterator> InputT_it_rg;
    InputT_it_rg string_input = create_ls(inputT.begin(), inputT.end());

    typedef it_rg<typename InputSA::iterator> InputSA_it_rg;
    InputSA_it_rg sa_input = create_ls(inputSA.begin(), inputSA.end());

    return check_sa<InputT_it_rg, InputSA_it_rg, typename InputSA::size_type>(string_input, sa_input);
}

} // namespace sacheck
