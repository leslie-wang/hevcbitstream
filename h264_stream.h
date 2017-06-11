/* 
 * h264bitstream - a library for reading and writing H.264 video
 * Copyright (C) 2005-2007 Auroras Entertainment, LLC
 * Copyright (C) 2008-2011 Avail-TVN
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

#ifndef _H264_STREAM_H
#define _H264_STREAM_H        1

#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "bs.h"

#ifdef __cplusplus
extern "C" {
#endif

//Appendix E. Table E-1  Meaning of sample aspect ratio indicator
#define SAR_Unspecified  0           // Unspecified
#define SAR_1_1        1             //  1:1
#define SAR_12_11      2             // 12:11
#define SAR_10_11      3             // 10:11
#define SAR_16_11      4             // 16:11
#define SAR_40_33      5             // 40:33
#define SAR_24_11      6             // 24:11
#define SAR_20_11      7             // 20:11
#define SAR_32_11      8             // 32:11
#define SAR_80_33      9             // 80:33
#define SAR_18_11     10             // 18:11
#define SAR_15_11     11             // 15:11
#define SAR_64_33     12             // 64:33
#define SAR_160_99    13             // 160:99
                                     // 14..254           Reserved
#define SAR_Extended      255        // Extended_SAR

extern int find_nal_unit(uint8_t* buf, int size, int* nal_start, int* nal_end);

extern int rbsp_to_nal(const uint8_t* rbsp_buf, const int* rbsp_size, uint8_t* nal_buf, int* nal_size);
extern int nal_to_rbsp(const uint8_t* nal_buf, int* nal_size, uint8_t* rbsp_buf, int* rbsp_size);

extern int more_rbsp_data(bs_t* bs);
extern int more_rbsp_trailing_data(bs_t* b);

extern int _read_ff_coded_number(bs_t* b);
extern void _write_ff_coded_number(bs_t* b, int n);

extern void debug_bytes(uint8_t* buf, int len);
    
void read_rbsp_trailing_bits(bs_t* b);


// file handle for debug output
extern FILE* h264_dbgfile;

#ifdef __cplusplus
}
#endif

#endif
