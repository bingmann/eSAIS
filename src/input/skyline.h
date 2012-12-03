/******************************************************************************
 * src/input/skyline.h
 *
 * Algorithm to generate worst-case input for SAIS and divsufsort as a stream.
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

namespace skyline {

template <typename Container>
void generate(Container& out)
{
    if (out.size() == 0) return;

    size_t len = out.size();

    size_t log2len = ceil( log(len+1) / log(2) );

    for (size_t i = 0; i < len; ++i)
    {
        // count number of least-significant binary zeros in i+1
        size_t j = i+1;

        unsigned int zc = 0;
        while ( (j & 1) == 0 ) {
            zc++; j >>= 1;
            assert( j != 0 );
        }

        zc = std::min<size_t>(zc, log2len - 1); // repeat a if size is not a power of 2

        out[i] = 'a' + (log2len - zc - 1);
        //std::cout << i << " : " << (out[i] - 'a') << " - " << out[i] << "\n";
    }
}

} // namespace skyline
