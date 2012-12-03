/******************************************************************************
 * src/malloccount.cc
 *
 * malloc() allocation counter based on http://ozlabs.org/~jk/code/
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

#if 0
base=malloccount
cmd="gcc -Wall -Wstrict-prototypes -fPIC -D_GNU_SOURCE -ldl \
    -shared -o lib$base.so $base.c"
echo $cmd
$cmd
if [ $? != 0 ]; then echo compile failed; return; fi
echo -e "\nrun with:\nLD_PRELOAD=$PWD/lib$base.so <program>"
return
#endif

static const int log_operations = 0;

/*
  gcc -Wall -Wstrict-prototypes -fPIC -D_GNU_SOURCE -ldl -shared -o libmalloccount.so malloccount.c
*/

#include <stdlib.h>
#include <dlfcn.h>
#include <string.h>
#include <stdio.h>
#include <locale.h>

static void *(*real_malloc)(size_t) = NULL;
static void (*real_free)(void *) = NULL;
static void *(*real_realloc)(void *, size_t) = NULL;

#define INIT_HEAP_SIZE 1024*1024
static char init_heap[INIT_HEAP_SIZE];
static size_t init_heap_use = 0;

#define min(a,b) (a) < (b) ? (a) : (b)

long long peak = 0, curr = 0, total = 0;

static inline void inc_count(size_t inc)
{
    if ((curr += inc) > peak) peak = curr;
    total += inc;
}

static inline void dec_count(size_t dec)
{
    curr -= dec;
}

extern size_t my_memory_current(void)
{
    return curr;
}

extern size_t my_memory_maximum(void)
{
    return peak;
}

void *malloc(size_t size)
{
    void *ret;

    if (!size) return NULL;

    if (real_malloc)
    {
        ret = (*real_malloc)(size + sizeof(size_t));

        inc_count(size);
        if (log_operations && size > 1124*1024) {
            printf("malloc(%lld) = %p   (curr %'lld)\n", (long long)size, ret, curr);
        }

        *((size_t *)ret) = size;

        return ret + sizeof(size_t);
    }
    else
    {
        if (init_heap_use + sizeof(size_t) + size > INIT_HEAP_SIZE) {
            fprintf(stderr, "init heap full!\n");
            exit(EXIT_FAILURE);
        }

        ret = init_heap + init_heap_use;
        init_heap_use += sizeof(size_t) + size;
        *((size_t *)ret) = size;

        //printf("malloc(%lld) = %p   on init heap\n", (long long)size, ret + sizeof(size_t));

        return ret + sizeof(size_t);
    }
}

void free(void *ptr)
{
    long long size;

    if (!ptr) return;

    if ((char*)ptr >= init_heap &&
        (char*)ptr <= init_heap + init_heap_use)
    {
        //printf("free(%p) = on init heap\n", ptr);
        return;
    }

    if (!real_free) {
        printf("free(%p) = without real_free !!!\n", ptr);
        return;
    }

    ptr -= sizeof(size_t);

    size = *(size_t *)ptr;
    dec_count(size);

    if (log_operations && size > 1124*1024) {
        printf("  free(%lld) <- %p   (curr %'lld)\n", size, ptr, curr);
    }

    (*real_free)(ptr);
}

void *calloc(size_t nmemb, size_t size)
{
    void *ret;
    size *= nmemb;
    if (!size)
        return NULL;
    ret = malloc(size);
    memset(ret, 0, size);
    return ret;
}

void *realloc(void *ptr, size_t size)
{
    void *newptr;
    size_t oldsize;

    if (ptr >= (void *)init_heap &&
        ptr <= (void *)init_heap + init_heap_use)
    {
        //printf("realloc(%p) = on init heap\n", ptr);

        oldsize = *(size_t *)(ptr - sizeof(size_t));
        if (oldsize <= size) {
            *(size_t *)(ptr - sizeof(size_t)) = size;
            return ptr;
        }
        else {
            newptr = malloc(size);
            memcpy(newptr, ptr, min(oldsize, size));
            free(ptr);
            return newptr;
        }
    }

    if (!size) {
        free(ptr);
        return NULL;
    }

    if (!ptr) return malloc(size);

    ptr -= sizeof(size_t);
    oldsize = *(size_t *)ptr;

    dec_count(oldsize);
    inc_count(size);

    newptr = (*real_realloc)(ptr, size + sizeof(size_t));

    if (log_operations && size > 1124*1024) {
        if (newptr == ptr)
            printf("realloc(%lld -> %lld) = %p   (curr %'lld)\n", (long long)oldsize, (long long)size, newptr, curr);
        else
            printf("realloc(%lld -> %lld) = %p -> %p   (curr %'lld)\n", (long long)oldsize, (long long)size, ptr, newptr, curr);
    }

    *((size_t *)newptr) = size;

    return newptr + sizeof(size_t);
}

__attribute__((constructor)) void init(void)
{
    char *error;

    setlocale(LC_NUMERIC, "");

    dlerror();

    real_malloc = dlsym(RTLD_NEXT, "malloc");
    if ((error = dlerror()) != NULL) {
        fprintf(stderr, "malloccount: %s\n", error);
        exit(EXIT_FAILURE);
    }

    real_realloc = dlsym(RTLD_NEXT, "realloc");
    if ((error = dlerror()) != NULL) {
        fprintf(stderr, "malloccount: %s\n", error);
        exit(EXIT_FAILURE);
    }

    real_free = dlsym(RTLD_NEXT, "free");
    if ((error = dlerror()) != NULL) {
        fprintf(stderr, "malloccount: %s\n", error);
        exit(EXIT_FAILURE);
    }
}

__attribute__((destructor)) void fini(void)
{
    printf("MALLOCCOUNT ### bytes allocated: total: %'lld, peak: %'lld, curr: %'lld\n",
           total, peak, curr);
}
