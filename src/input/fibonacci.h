/******************************************************************************
 * src/input/fibonacci.h
 *
 * Algorithm to generate a Fibonacci word
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

namespace fibonacci {

template <typename Container>
void generate(Container& out)
{
    if (out.size() == 0) return;

    out[0] = 'a';

    size_t j = 1;
    for (size_t i = 0; j < out.size(); ++i)
    {
	if (out[i] == 'a')
	{
	    out[j++] = 'b';
	    if (j >= out.size()) break;

	    out[j++] = 'a';
	}
	else
	{
	    out[j++] = 'a';
	}
    }
}

} // namespace fibonacci
