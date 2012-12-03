/******************************************************************************
 * src/external/lcp.h
 *
 * Algorithm to calculate the LCP array from a string and suffix array in
 * linear time. Based on the ideas of Kasai et al. (2001)
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

namespace lcp {

template <typename InputStream>
struct sa_index_stream_type
{
    typedef typename InputStream::value_type offset_type;

    typedef tuples::triple<offset_type, offset_type, offset_type> value_type;

private:
    size_t      m_counter;

    InputStream&  m_input;

    value_type  m_curr;

public:
    sa_index_stream_type(InputStream& input)
        : m_counter(0), m_input(input)
    {
        if (!m_input.empty())
            m_curr = value_type(*m_input, m_counter++, 0);
    }

    const value_type& operator* () const
    {
        return m_curr;
    }

    sa_index_stream_type& operator++ ()
    {
        ++m_input;
        if (!m_input.empty())
            m_curr = value_type(*m_input, m_counter++, m_curr.first);
        return *this;
    }

    bool empty() const
    {
        return m_input.empty();
    }
};

template <typename InputStream>
struct pair_2nd_stream_type
{
    typedef typename InputStream::value_type pair_type;

    typedef typename pair_type::first_type value_type;

private:
    InputStream&  m_input;

public:
    pair_2nd_stream_type(InputStream& input)
        : m_input(input)
    {
    }

    const value_type& operator* () const
    {
        return m_input->second;
    }

    pair_2nd_stream_type& operator++ ()
    {
        ++m_input;
        return *this;
    }

    bool empty() const
    {
        return m_input.empty();
    }
};

template <typename StringContainer, typename SAContainer, typename LCPContainer>
void lcparray_stxxl_kasai(const StringContainer& string, const SAContainer& SA, LCPContainer& lcp)
{
    assert( string.size() == SA.size() );

    // *** Generate ISA

    typedef typename StringContainer::value_type     alphabet_type;
    typedef typename SAContainer::value_type         offset_type;
    typedef typename SAContainer::size_type          size_type;

    // pair and triples of offsets
    typedef tuples::pair<offset_type,offset_type>               offset_pair_type;
    typedef tuples::triple<offset_type,offset_type,offset_type> offset_triple_type;

    // define stream iterating over a STXXL vector
    typedef stxxl::stream::iterator2stream< typename SAContainer::const_iterator > sa_stream_type;
    sa_stream_type sa_stream (SA.begin(), SA.end());

    sa_index_stream_type<sa_stream_type> sa_index_stream(sa_stream);

    // comparator for ISA sort: 1st component ascending
    typedef tuples::triple_less1st<offset_triple_type> offset_triple_less1st_type;
    offset_triple_less1st_type offset_triple_less1st;

    // runs creator for ISA stream
    typedef stxxl::stream::sort<sa_index_stream_type<sa_stream_type>, offset_triple_less1st_type, block_size> isa_sort_type;
    isa_sort_type isa_sort(sa_index_stream, offset_triple_less1st, ram_use / 8, ram_use / 8);

    // *** Output sorter: 

    typedef tuples::pair_less1st<offset_pair_type> offset_pair_less1st_type;
    offset_pair_less1st_type offset_pair_less1st;

    typedef stxxl::sorter<offset_pair_type, offset_pair_less1st_type, block_size> lcp_sorter_type;
    lcp_sorter_type lcp_sorter(offset_pair_less1st, ram_use / 8);

    // *** Kasai algorithm: iterate over ISA

    {
        size_type h = 0;  // current height
        size_type i = 0;  // ISA index counter;

        size_type N10 = string.size() / 10;

        lcp_sorter.push( offset_pair_type(0,0) );

        std::vector<alphabet_type> stringRAM (string.begin(), string.end());

        while ( !isa_sort.empty() )
        {
            offset_type k = isa_sort->second;     // k = ISA[i]

            if ( (i % N10) == 0 ) (std::cout << (i / N10 * 10) << "%...").flush();

            if (k > offset_type(0))
            {
                size_type j = isa_sort->third;   // j = SA[k-1];

                while(i + h < stringRAM.size() && j + h < stringRAM.size() &&
                      stringRAM[i+h] == stringRAM[j+h])
                    h++;

                lcp_sorter.push( offset_pair_type(k,h) );
            }
            if (h > 0) h--;

            ++isa_sort, ++i;
        }
        std::cout << "\n";
    }

    // *** Collect output

    if (0) // simple way
    {
        lcp.resize(string.size());
        lcp_sorter.sort();

        offset_type i = 0;
        while( !lcp_sorter.empty() )
        {
            //std::cout << "(" << lcp_sorter->first << "," << lcp_sorter->second << ")\n";

            assert( lcp_sorter->first == i );
            lcp[i] = lcp_sorter->second;
            ++lcp_sorter, ++i;
        }
    }
    else // efficient way
    {
        pair_2nd_stream_type<lcp_sorter_type> pair_2nd_stream(lcp_sorter);

        lcp.resize(string.size());
        lcp_sorter.sort();

        stxxl::stream::materialize(pair_2nd_stream, lcp.begin(), lcp.end());
    }
}

} // namespace lcp
