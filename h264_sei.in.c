/*
 * h264bitstream - a library for reading and writing H.264 video
 * Copyright (C) 2005-2007 Auroras Entertainment, LLC
 * Copyright (C) 2008-2011 Avail-TVN
 * Copyright (C) 2012 Alex Izvorski
 *
 * This file is written by Leslie Wang <wqyuwss@gmail.com>
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

#include <stdio.h>
#include <stdlib.h> // malloc
#include <string.h> // memset

sei_t* sei_new()
{
    sei_t* s = (sei_t*)calloc(1, sizeof(sei_t));
    memset(s, 0, sizeof(sei_t));
    s->data = NULL;
    return s;
}

void sei_free(sei_t* s)
{
    switch( s->payloadType ) {
        default:
            if ( s->data != NULL ) free(s->data);
    }
    free(s);
}

void read_sei_end_bits(bs_t* b )
{
    // if the message doesn't end at a byte border
    if ( !bs_byte_aligned( b ) )
    {
        if ( !bs_read_u1( b ) ) fprintf(stderr, "WARNING: bit_equal_to_one is 0!!!!\n");
        while ( ! bs_byte_aligned( b ) )
        {
            if ( bs_read_u1( b ) ) fprintf(stderr, "WARNING: bit_equal_to_zero is 1!!!!\n");
        }
    }
    
    read_rbsp_trailing_bits(b);
}

#end_preamble

#function_declarations

// D.1 SEI payload syntax
void structure(sei_payload)( sei_t* s, bs_t* b )
{
    int i;
    switch( s->payloadType )
    {
        default:
            if( is_reading )
            {
                s->data = (uint8_t*)calloc(1, s->payloadSize);
            }
            
            for ( i = 0; i < s->payloadSize; i++ )
                value( s->data[i], u8 );
    }
    
    //if( is_reading )
    //    read_sei_end_bits(h, b);
}
