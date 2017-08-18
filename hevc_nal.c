/* 
 * hevcbitstream - a library for reading and writing H.264 video
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

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "bs.h"
#include "hevc_stream.h"

/**
 Create a new HEVC stream object.  Allocates all structures contained within it.
 @return    the stream object
 */
hevc_stream_t* hevc_new()
{
    hevc_stream_t* h = (hevc_stream_t*)calloc(1, sizeof(hevc_stream_t));

    h->nal = (hevc_nal_t*)calloc(1, sizeof(hevc_nal_t));
    
    // initialize tables
    for ( int i = 0; i < 32; i++ ) { h->sps_table[i] = (hevc_sps_t*)calloc(1, sizeof(hevc_sps_t)); }
    for ( int i = 0; i < 256; i++ ) { h->pps_table[i] = (hevc_pps_t*)calloc(1, sizeof(hevc_pps_t)); }

    h->vps = (hevc_vps_t*)calloc(1, sizeof(hevc_vps_t));
    h->sps = (hevc_sps_t*)calloc(1, sizeof(hevc_sps_t));
    h->pps = (hevc_pps_t*)calloc(1, sizeof(hevc_pps_t));
    h->aud = (hevc_aud_t*)calloc(1, sizeof(hevc_aud_t));
    /* TODO: support SEI
    h->num_seis = 0;
    h->seis = NULL;
    h->sei = NULL;  //This is a TEMP pointer at whats in h->seis...
     */
    h->sh = (hevc_slice_header_t*)calloc(1, sizeof(hevc_slice_header_t));
    h->slice_data = (hevc_slice_data_rbsp_t*)calloc(1, sizeof(hevc_slice_data_rbsp_t));

    return h;
}


/**
 Free an existing HEVC stream object.  Frees all contained structures.
 @param[in,out] h   the stream object
 */
void hevc_free(hevc_stream_t* h)
{
    free(h->nal);

    for ( int i = 0; i < 32; i++ ) { free( h->sps_table[i] ); }
    for ( int i = 0; i < 256; i++ ) { free( h->pps_table[i] ); }

    /* TODO: support SEI
    if(h->seis != NULL)
    {
        for( int i = 0; i < h->num_seis; i++ )
        {
            sei_t* sei = h->seis[i];
            sei_free(sei);
        }
        free(h->seis);
    }
     */
    free(h->slice_data);
    free(h->sh);
    
    free(h->aud);
    free(h->pps);
    free(h->sps);
    free(h->vps);
    
    free(h);
}

/**
 Read only the NAL headers (enough to determine unit type) from a byte buffer.
 @return unit type if read successfully, or -1 if this doesn't look like a nal
*/
int peek_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    bs_t* b = bs_new(buf, size);

    /* forbidden_zero_bit */ bs_skip_u(b, 1);
    nal->nal_unit_type = bs_read_u(b, 6);
    nal->nal_layer_id = bs_read_u(b, 6);
    nal->nal_temporal_id_plus1 = bs_read_u(b, 3);

    bs_free(b);

    // basic verification, per 7.4.1
    if ( nal->nal_unit_type <= 0 || nal->nal_unit_type > MAX_HEVC_VAL_UNIT_TYPE ) { return -1; }

    return nal->nal_unit_type;
}
