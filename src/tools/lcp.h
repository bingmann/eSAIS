/******************************************************************************
 * src/tools/lcp.h
 *
 * Algorithm to calculate the LCP array from a string and suffix array in
 * linear time.
 *
 * Based on Kasai, T., Lee, G., Arimura, H., Arikawa, S., & Park, K. (2001).
 * Linear-time longest-common-prefix computation in suffix arrays and its
 * applications. 2001. CPM'01 (pp. 181-192). Springer
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

template <typename alphabet_type, typename offset_type>
void lcparray(const alphabet_type* string, const offset_type* SA, std::vector<offset_type>& lcp, unsigned int N)
{
    // calculate ISA
    std::vector<offset_type> ISA (N);

    for(unsigned int i = 0; i < N; i++)
	ISA[SA[i]] = i;

    // allocate lcp vector and iterate over ISA
    lcp.resize(N);

    unsigned int h = 0; // current height

    lcp[0] = 0;
    for(unsigned int i = 0; i < N; i++)
    {
	unsigned int k = ISA[i];
 	if (k > 0) {
	    unsigned int j = SA[k-1];
	    while(i + h < N && j + h < N &&
		  string[i+h] == string[j+h])
		h++;
	    lcp[k] = h;
	}
	if (h > 0) h--;
    }
}

template <typename alphabet_type, typename offset_type>
void lcparray(const std::vector<alphabet_type>& string, const std::vector<offset_type>& SA, std::vector<offset_type>& lcp)
{
    assert(string.size() == SA.size());

    lcparray(string.data(), SA.data(), lcp, string.size());
}

template <typename StringContainer, typename SAContainer, typename LCPContainer>
void lcparray_naive(const StringContainer& string, const SAContainer& SA, LCPContainer& lcp)
{
    assert( string.size() == SA.size() );
    lcp.resize(string.size());

    lcp[0] = 0;
    for (size_t i = 1; i < string.size(); ++i)
    {
        typename StringContainer::const_iterator
            pa = string.begin() + SA[i-1],
            pb = string.begin() + SA[i];

        size_t h = 0;

        while (pa != string.end() && pb != string.end() && *pa == *pb)
            ++pa, ++pb, ++h;

        lcp[i] = h;
    }
}

} // namespace lcp
