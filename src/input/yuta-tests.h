/******************************************************************************
 * src/input/yuta-tests.h
 *
 * Algorithm to generate test input according to Yuta Mori generator2
 *
 * Adapted from http://homepage3.nifty.com/wpage/junk/generator2.c.bz2 and
 * enhanced to generate larger strings.
 *
 ******************************************************************************
 * Copyright (C) Yuta Mori
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

namespace yuta_tests {

/* 0x00, 0x00, 0x00, 0x01, ..., 0x00, 0xfe, 0x00, 0xff */
template <typename Container>
void generate1(Container& out)
{
    for (size_t i = 0; i+1 < out.size();)
    {
	for(size_t j = 0x00; j <= 0xff && i+1 < out.size(); ++j) {
	    out[i++] = 0x00;
	    out[i++] = j;
	}
    }
}

/* 0x00, 0xff, 0x00, 0xfe, ..., 0x00, 0x01, 0x00, 0x00 */
template <typename Container>
void generate2(Container& out)
{
    for (size_t i = 0; i+1 < out.size();)
    {
	for(ssize_t j = 0xff; 0x00 <= j && i+1 < out.size(); --j) {
	    out[i++] = 0x00;
	    out[i++] = j;
	}
    }
}

/* 0x00, 0x00, 0x01, 0x00, 0x02, ..., 0x00, 0xfe, 0x00, 0xff,
   0x01, 0x01, 0x02, 0x01, 0x03, ..., 0x01, 0xfe, 0x01, 0xff,
   ...,
   0xfd, 0xfd, 0xfe, 0xfd, 0xff,
   0xfe, 0xfe, 0xff,
   0xff */
template <typename Container>
void generate3(Container& out)
{
    for (size_t i = 0; i+1 < out.size();)
    {
	for(size_t j = 0x00; j <= 0xff && i < out.size(); ++j) {
	    out[i++] = j;
	    for(size_t k = j + 1; k <= 0xff && i+1 < out.size(); ++k) {
		out[i++] = j;
		out[i++] = k;
	    }
	}
    }
}

/* 0xff, 0xff, 0xfe, 0xff, 0xfd, ..., 0xff, 0x01, 0xff, 0x00,
   0xfe, 0xfe, 0xfd, 0xfe, 0xfc, ..., 0xfe, 0x01, 0xfe, 0x00,
   ...,
   0x02, 0x02, 0x01, 0x02, 0x00,
   0x01, 0x01, 0x00,
   0x00 */
template <typename Container>
void generate4(Container& out)
{
    for (size_t i = 0; i+1 < out.size();)
    {
	for(ssize_t j = 0xff; 0x00 <= j && i < out.size(); --j) {
	    out[i++] = j;
	    for(ssize_t k = j - 1; 0x00 <= k && i+1 < out.size(); --k) {
		out[i++] = j;
		out[i++] = k;
	    }
	}
    }
}

template <typename Container>
void generate5(Container& out)
{
    int fibsize = 1 << 8;
    std::vector<typename Container::value_type> fibstr(fibsize);
    fibstr[0] = 0;
    fibstr[1] = 1;
    for(int i = 0, j = 0, k = 2; k < fibsize; ++i, ++k) {
	fibstr[k] = fibstr[i];
	if (i == j) {
	    i = -1;
	    j = k - 1;
	}
    }

    for(size_t i = 0, j = 0; i < out.size(); ++i) {
	out[i] = (fibstr[j] + (i / fibsize)) % 256;
	j = (j+1) % fibsize;
    }
}

} // namespace yutatests
