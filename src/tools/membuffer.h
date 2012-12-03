/******************************************************************************
 * src/tools/losertree.h
 *
 * Dumb memory buffer using malloc which does not initialize the values.
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

template <typename Type>
class membuffer
{
private:
    /// pointer to allocated memory area
    Type*       m_ptr;
    
    /// size of allocated memory
    size_t      m_size;

    /// protect copy-constructor
    membuffer(const membuffer& b);

    /// protect assignment operator
    membuffer& operator=(const membuffer& b);

public:

    /// Allocate memory buffer
    inline membuffer(size_t size)
        : m_ptr( (Type*)malloc(size * sizeof(Type)) ),
          m_size(size)
    {
    }

    /// Deallocate memory buffer
    inline ~membuffer()
    {
        free(m_ptr);
    }

    /// Accessor to elements
    inline Type& operator[](size_t n) const
    {
        return m_ptr[n];
    }

    inline Type* data() { return m_ptr; }
    inline const Type* data() const { return m_ptr; }

    inline Type* begin() { return m_ptr; }
    inline Type* end() { return m_ptr + m_size; }

    inline const Type* begin() const { return m_ptr; }
    inline const Type* end() const { return m_ptr + m_size; }
};
