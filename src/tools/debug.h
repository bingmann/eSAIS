/******************************************************************************
 * src/tools/debug.h
 *
 * Debugging macros and memory debugging utility
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

#define DBGX(dbg,X)   do { if (dbg) { std::cout << X; } } while(0)

#define DBG(dbg,X)    DBGX(dbg, __FUNCTION__ << "() " << X << std::endl)

#define DBG1(dbg,X)   DBGX(dbg, __FUNCTION__ << "() " << X)
#define DBG2(dbg,X)   DBGX(dbg, X)
#define DBG3(dbg,X)   DBGX(dbg, X << std::endl)

#define DBG_ST_ARRAY(dbg,text,X)  do {                                          \
    if (dbg)									\
    {										\
        X.sort();                                                               \
        DBG1(dbg, text << " line " << __LINE__ << " : i - " #X "[i]");          \
        if (!X.empty()) DBG2(dbg," " << (*X).description); DBG2(dbg,std::endl); \
        for (size_t i = 0; !X.empty(); ++i, ++X)                                \
            std::cout << i << " : " << *X << "\n";                              \
        X.rewind();                                                             \
    }										\
} while(0)

#ifdef MALLOCCOUNT

// prototypes to call the functions in malloccount.c
extern "C" size_t my_memory_current();
extern "C" size_t my_memory_maximum();

static size_t dmlastsize = 0;
static size_t dklastsize = 0;
stxxl::stats_data g_stats;

template <typename T>
static inline std::string numFormat(const T& v) {
    std::ostringstream ss;
    ss.imbue(std::locale(""));
    ss << v;
    return ss.str();
}

#define DBGMEM(text) do {                                                       \
    stxxl::stats_data nst = *stxxl::stats::get_instance();              \
    stxxl::stats_data st = nst - g_stats;                               \
    std::cout << text << " - MEM: total= " << numFormat(my_memory_current()) << " delta= " << numFormat((ssize_t)my_memory_current() - (ssize_t)dmlastsize) << " maxalloc= " << numFormat(my_memory_maximum()) << "\n"; \
    std::cout << text << " - DISK: total= " << numFormat(stxxl::block_manager::get_instance()->current_allocated()) << " delta= " << numFormat((long)stxxl::block_manager::get_instance()->current_allocated() - (long)dklastsize) << " bytes maxalloc= " << numFormat(stxxl::block_manager::get_instance()->max_allocated()) << " io volume = " << numFormat(st.get_read_volume() + st.get_written_volume()) << "\n"; \
    dmlastsize = my_memory_current();                                           \
    dklastsize = stxxl::block_manager::get_instance()->current_allocated();     \
} while (0)

#else

#define DBGMEM(text) do { } while (0)

#endif
