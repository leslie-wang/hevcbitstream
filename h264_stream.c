/* 
 * h264bitstream - a library for reading and writing H.264 video
 * Copyright (C) 2005-2007 Auroras Entertainment, LLC
 * Copyright (C) 2008-2011 Avail-TVN
 * Copyright (C) 2012 Alex Izvorski
 *
 * Written by Alex Izvorski <aizvorski@gmail.com> and Alex Giladi <alex.giladi@gmail.com>
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "bs.h"
#include "h264_stream.h"
#include "h264_sei.h"

FILE* h264_dbgfile = NULL;

#define printf(...) fprintf((h264_dbgfile == NULL ? stdout : h264_dbgfile), __VA_ARGS__)

/** 
 Calculate the log base 2 of the argument, rounded up. 
 Zero or negative arguments return zero 
 Idea from http://www.southwindsgames.com/blog/2009/01/19/fast-integer-log2-function-in-cc/
 */
int intlog2(int x)
{
    int log = 0;
    if (x < 0) { x = 0; }
    while ((x >> log) > 0)
    {
        log++;
    }
    if (log > 0 && x == 1<<(log-1)) { log--; }
    return log;
}

int is_slice_type(int slice_type, int cmp_type)
{
    if (slice_type >= 5) { slice_type -= 5; }
    if (cmp_type >= 5) { cmp_type -= 5; }
    if (slice_type == cmp_type) { return 1; }
    else { return 0; }
}

int more_rbsp_data(bs_t* bs)
{
    // TODO this version handles reading only. writing version?

    // no more data
    if (bs_eof(bs)) { return 0; }

    // no rbsp_stop_bit yet
    if (bs_peek_u1(bs) == 0) { return 1; }

    // next bit is 1, is it the rsbp_stop_bit? only if the rest of bits are 0
    bs_t bs_tmp;
    bs_clone(&bs_tmp, bs);
    bs_skip_u1(&bs_tmp);
    while(!bs_eof(&bs_tmp))
    {
        // A later bit was 1, it wasn't the rsbp_stop_bit
        if (bs_read_u1(&bs_tmp) == 1) { return 1; }
    }

    // All following bits were 0, it was the rsbp_stop_bit
    return 0;
}

int more_rbsp_trailing_data(bs_t* b) { return !bs_eof(b); }

int _read_ff_coded_number(bs_t* b)
{
    int n1 = 0;
    int n2;
    do 
    {
        n2 = bs_read_u8(b);
        n1 += n2;
    } while (n2 == 0xff);
    return n1;
}

void _write_ff_coded_number(bs_t* b, int n)
{
    while (1)
    {
        if (n > 0xff)
        {
            bs_write_u8(b, 0xff);
            n -= 0xff;
        }
        else
        {
            bs_write_u8(b, n);
            break;
        }
    }
}

void debug_bytes(uint8_t* buf, int len)
{
    int i;
    for (i = 0; i < len; i++)
    {
        printf("%02X ", buf[i]);
        if ((i+1) % 16 == 0) { printf ("\n"); }
    }
    printf("\n");
}

//7.3.2.11 RBSP trailing bits syntax
void read_rbsp_trailing_bits(bs_t* b)
{
    /* rbsp_stop_one_bit */ bs_skip_u(b, 1);
    
    while( !bs_byte_aligned(b) )
    {
        /* rbsp_alignment_zero_bit */ bs_skip_u(b, 1);
    }
}


