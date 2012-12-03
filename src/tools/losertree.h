/******************************************************************************
 * src/tools/losertree.h
 *
 * Loser tree implementation with improved caching of equality.
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

template <typename Comparator>
class LoserTree3Way
{
private:

    /// marker inside loser tree to save equality comparisons
    static const uint32_t       EqualMark = 0x80000000;

    /// marker inside loser tree for finished players
    static const uint32_t       DoneMark = 0x7FFFFFFF;

    /// the comparator object of this tree
    const Comparator&	m_cmp;

    /// number of players
    unsigned int	m_size;

    /// the tree of size n-1
    std::vector<uint32_t> m_tree;

    /// saves the equality mark of top after extraction
    bool		m_top_equal;

    /// all players are done
    bool                m_done;

public:

    /// Initialize loser tree and play initial games.
    LoserTree3Way(const Comparator& cmp)
        : m_cmp(cmp)
    {
    }

    /// Play initial round of comparison games.
    bool play_initial(unsigned int size)
    {
        m_size = size;
        m_top_equal = false;
        m_done = false;

        // initialize loser tree of appropriate size (round up to power of two)
        unsigned int treesize = (1 << (int)(log2(size - 1) + 2)) - 1;
        if (treesize == 0) treesize = 1;
        m_tree.resize(treesize, DoneMark);

        // fill in lowest level: all ascending player numbers
        int levelput = m_tree.size() / 2;
        for (unsigned int i = 0; i < m_size; ++i)
        {
            m_tree[levelput + i] = m_cmp.exists(i) ? i : DoneMark;
        }

        unsigned int levelsize = levelput + 1;

        // construct higher levels iteratively from bottom up
        while ( levelsize > 1 )
        {
            levelsize = levelsize / 2;
            int levelget = levelput;
            levelput /= 2;

            for (unsigned int i = 0; i < levelsize; ++i)
            {
                if ( m_tree[levelget + 2*i + 1] == DoneMark )
                    m_tree[levelput + i] = m_tree[levelget + 2*i] & ~EqualMark;
                else if ( m_tree[levelget + 2*i] == DoneMark )
                    m_tree[levelput + i] = m_tree[levelget + 2*i + 1] & ~EqualMark;
                else
                {
                    int cmp = m_cmp( m_tree[levelget + 2*i] & ~EqualMark, m_tree[levelget + 2*i + 1] & ~EqualMark );
                    if (cmp < 0)
                        m_tree[levelput + i] = m_tree[levelget + 2*i] & ~EqualMark;
                    else if (cmp == 0)
                        m_tree[levelput + i] = m_tree[levelget + 2*i] | EqualMark;
                    else
                        m_tree[levelput + i] = m_tree[levelget + 2*i + 1] & ~EqualMark;
                }
            }
        }

        m_top_equal = false;
        return (m_done = (m_tree[0] == DoneMark));
    }

    /// Get the top element of the tree: index of the losing player.
    unsigned int top() const
    {
        return m_tree[0] & ~EqualMark;
    }

    /// Get flag whether the current losing player's key is equal to the previous one's. If this returns true,
    /// the previous key is _always_ equal to this one. However, this flag _need_not_ be set if the previous
    /// keys are equal. Due to optimization in the comparison tree, it might not have been necessary to
    /// compare the two player's keys directly.
    bool top_equal() const
    {
        return m_top_equal;
    }

    /// Extract the currently losing player and replay the tree with the new updated values for the losing
    /// player. (The updates are performed outside of the tree implementation.)
    bool replay()
    {
        m_top_equal = (m_tree[0] & EqualMark);	// save equality marker for query

        int top = m_tree[0] & ~EqualMark;

        int p = (m_tree.size() / 2) + top;

        if (!m_cmp.exists(top)) top = DoneMark;	// mark player as done

        while( p > 0 )
        {
            m_tree[p] = top;

            p -= (p+1) % 2;	// round down to left node position

            if (m_tree[p] == DoneMark)
                top = m_tree[p+1] & ~EqualMark;
            else if (m_tree[p+1] == DoneMark)
                top = m_tree[p] & ~EqualMark;
            else
                if (m_tree[p/2] & EqualMark)
                {
                    // do nothing, the previous comparison of the parent key was equal -> this key can only be
                    // greater-equal than that one.
                    top = m_tree[p+1] & ~EqualMark;
                }
                else
                {
                    int cmp = m_cmp( m_tree[p] & ~EqualMark, m_tree[p+1] & ~EqualMark );
                    if ( cmp < 0 )
                        top = m_tree[p] & ~EqualMark;
                    else if (cmp == 0)
                        top = m_tree[p] | EqualMark;
                    else
                        top = m_tree[p+1] & ~EqualMark;
                }

            p /= 2;
        }
        assert(p == 0);

        m_tree[0] = top;
        return (m_done = (m_tree[0] == DoneMark));
    }

    /// Return true if all players are done
    bool done() const
    {
        return m_done;
    }

    void print() const
    {
        unsigned int levelsize = 1;
        unsigned int j = 0;

        for (unsigned int i = 0; i < m_tree.size(); ++i)
        {
            if (i >= j + levelsize) {
                std::cout << "\n";
                j = i; levelsize *= 2;
            }
            std::cout << (m_tree[i] & ~EqualMark) << ((m_tree[i] & EqualMark) ? 'e' : ' ') << " ";
        }
        std::cout << "\n";
    }
};
