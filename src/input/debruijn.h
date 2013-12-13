/******************************************************************************
 * src/input/debruijn.h
 *
 * Algorithm to generate De Bruijn input sequences.
 * extracted from FXT algorithm library http://www.jjj.de/fxt/fxtpage.html
 *
 ******************************************************************************
 * Original copyrighted by FXT algorithm library authors.
 * Modifications copyright (C) 2012 Timo Bingmann <tb@panthema.net>
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

namespace debruijn {

// Generate necklaces, iterative algorithm.
class necklace
{
public:
    ulong *a_;  // the string, NOTE: one-based
    ulong *dv_; // delta sequence of divisors of n
    ulong n_;   // length of strings
    ulong m1_;  // m-ary strings, m1=m-1
    ulong j_;   // period of the word (if necklaces)

public:
    necklace(ulong m, ulong n)
    {
        n_ = ( n ? n : 1 );  // at least 1
        m1_ = ( m>1 ? m-1 : 1); // at least 2
        a_ = new ulong[n_+1];
        dv_ = new ulong[n_+1];
        for (ulong j=1; j<=n; ++j)  dv_[j] = ( 0==(n_%j ) );  // divisors
        first();
    }

    ~necklace()
    {
        delete [] a_;
        delete [] dv_;
    }

    void first()
    {
        for (ulong j=0; j<=n_; ++j)  a_[j] = 0;
        j_ = 1;
    }

    const ulong * data()  const  { return  a_ + 1; }

    ulong next_pre()  // next pre-necklace
    // return j (zero when finished)
    {
        // Find rightmost digit that can be incremented:
        ulong j = n_;
        while ( a_[j] == m1_ )  { --j; }

        // Increment:
        // if ( 0==j_ )   return 0;  // last
        ++a_[j];

        // Copy periodically:
        for (ulong k=j+1; k<=n_; ++k)  a_[k] = a_[k-j];

        j_ = j;
        return  j;
    }

    bool is_necklace()  const
    {
        return ( 0!=dv_[j_] );  // whether j divides n
    }

    bool is_lyn()  const
    {
        return ( j_==n_ );  // whether j equals n
    }

    ulong next()  // next necklace
    {
        do
        {
            next_pre();
            if ( 0==j_ )  return 0;
        }
        while ( 0==dv_[j_] );  // until j divides n
        return j_;
    }

    ulong next_lyn()  // next Lyndon word
    {
        do
        {
            next_pre();
            if ( 0==j_ )  return 0;
        }
        while ( j_==n_ );  // until j equals n
        return j_;  // == n
    }
};

// Lexicographic minimal De Bruijn sequence.
class debruijn : public necklace
{
public:
    ulong i_;   // position of current digit in current string

public:
    debruijn(ulong m, ulong n)
        : necklace(m, n)
    { first_string(); }

    ~debruijn()  { ; }

    ulong first_string()
    {
        necklace::first();
        i_ = 1;
        return j_;
    }

    ulong next_string()  // make new string, return its length
    {
        necklace::next();
        i_ = (j_ != 0);
        return j_;
    }

    ulong next_digit()
    // Return current digit and move to next digit.
    // Return m if previous was last.
    {
        if ( i_ == 0 )  return necklace::m1_ + 1;
        ulong d = a_[ i_ ];
        if ( i_ == j_ )  next_string();
        else  ++i_;
        return d;
    }

    ulong first_digit()
    {
        first_string();
        return next_digit();
    }
};

template <typename Container>
void generate(int k, int n, Container& out)
{
    debruijn deb(k,n);

    if (pow(k,n) != out.size()) {
	std::cout << "Length of de Bruijn sequence k=" << k << ", n = " << n << " is " << (int)pow(k,n) << "\n";
	abort();
    }

    int i, o = 0;

    while( (i = deb.next_digit()) != k )
    {
	out[o++] = (i - 1 + '0');
    }
}

} // namespace debruijn
