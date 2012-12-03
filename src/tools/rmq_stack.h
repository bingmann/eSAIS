/******************************************************************************
 * src/tools/rmq_stack.h
 *
 * Simplistic, dynamic and succinct RMQ data structure that only allows querys
 * to the left of the current end of the array.
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

template <typename offset_type, typename value_type>
class RMQ_Stack
{
private:

    /// Pairs (index,minimum) saved on stack with most index and minimum ascending.
    struct Item {
        offset_type     index;
        value_type      minimum;

        Item(const offset_type& i, const value_type& m)
            : index(i), minimum(m) {}
    };

    /// Use a std::vector as stack-like structure
    typedef typename std::vector<Item> itemlist_type;
    typedef typename itemlist_type::iterator itemiter_type;
    typedef typename itemlist_type::const_iterator itemciter_type;

    /// Only for debugging purposes:
    offset_type         m_right;

    /// The current stack of minimum pairs
    itemlist_type       m_stack;

    static bool ItemCmpIndex(const Item& i, const offset_type& index) {
        return i.index < index;
    }

    static bool ItemCmpMinimum(const Item& i, const value_type& m) {
        return i.minimum < m;
    }

public:

    /// Constructor, does nothing.
    RMQ_Stack()
        : m_right(0)
    {
    }

    /// Append new value at the end and move right boundary
    void set(const offset_type& index, const value_type& val)
    {
        //assert(index == m_right);
        assert( m_right <= index );

        itemiter_type p = std::lower_bound(m_stack.begin(), m_stack.end(), val, ItemCmpMinimum);

        m_stack.erase(p, m_stack.end());
        m_stack.push_back( Item(index,val) );

        //++m_right;
        m_right = index+1;

        //print();
    }

    /// Perform left-ward range minimum query
    const value_type& query(const offset_type& left, const offset_type& right) const
    {
        assert(right+1 == m_right);
        assert(left <= m_right);
        stxxl::STXXL_UNUSED(right);

        itemciter_type p = std::lower_bound(m_stack.begin(), m_stack.end(), left, ItemCmpIndex);

        //std::cout << "rmq(" << left << "," << right << ") = " << p->minimum << "\n";

        return p->minimum;
    }

    /// Return the minimum over all current item
    const value_type& global_minimum() const
    {
        assert(m_stack.size());
        return m_stack.front().minimum;
    }

    /// Clear data structure
    void clear()
    {
        m_right = 0;
        m_stack.clear();
    }

    /// Return memory size of data structure
    size_t memsize() const
    {
        return m_stack.capacity() * sizeof(m_stack[0]) + sizeof(*this);
    }

    /// Print stack as tuples
    void print()
    {
        std::cout << "Stack: ";
        for (itemiter_type i = m_stack.begin(); i != m_stack.end(); ++i)
        {
            std::cout << "(" << i->index << "," << i->minimum << ") ";
        }
        std::cout << "\n";
    }
};

/**
 * Another simplistic dynamic RMQ data structure that only allows querys to the
 * left of the current end of the array. Furthermore, this structue also splits
 * the area into slabs: the left boundary can either be any query <= previous
 * slab boundary, or it must be exactly a slab-boundary. The class is used to
 * implement blocked RMQs
 */
template <typename offset_type, typename value_type, size_t slabsize>
class RMQ_Stack_Blocked
{
private:

    // RMQ stack for near querys <= previous slab boundary
    RMQ_Stack<offset_type,value_type>   rmq_near;

    // RMQ stack for slab spanning queries
    RMQ_Stack<offset_type,value_type>   rmq_slabs;

    // previous slab boundary: always a multiple of slabsize
    offset_type         thisslab;

public:

    RMQ_Stack_Blocked()
        : thisslab(0)
    {
    }

    /// Append new value at the end and move right boundary
    void set(const offset_type& index, const value_type& val)
    {
        assert(index >= thisslab);
        if (index >= offset_type(thisslab + slabsize))
        {
            // save global minimum of passed slab
            thisslab = (offset_type)(index / slabsize) * slabsize;      // new boundary
            //std::cout << "save global min: " << thisslab-1 << " - " << rmq_near.global_minimum() << "\n";
            rmq_slabs.set(thisslab-1, rmq_near.global_minimum());
            rmq_near.clear();
        }

        rmq_near.set(index, val);
    }

    /// Perform left-ward range minimum query
    const value_type& query(const offset_type& left, const offset_type& right) const
    {
        if (left >= thisslab)
        {
            return rmq_near.query(left,right);
        }
        else
        {
            //std::cout << "split query: " << "(" << left << "," << thisslab << ") and (" << thisslab << "," << right << ")\n";
            assert( left % slabsize == 0 );
            const value_type& min1 = rmq_slabs.query(left, thisslab-1);
            const value_type& min2 = rmq_near.query(thisslab, right);
            return (min1 < min2) ? min1 : min2;
        }
    }

    /// Print size of data structure
    void print_size()
    {
        std::cout << "RMQ_near: " << rmq_near.memsize() << " + RMQ_slabs: " << rmq_slabs.memsize() << "\n";
    }
};
