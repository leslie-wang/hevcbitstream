/*
  hevc_stream.c
  Created by leslie_qiwa@gmail.com on 6/8/17
 */


#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "bs.h"
#include "h264_stream.h"
#include "hevc_stream.h"
#include "h264_sei.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static void init_slice_hevc(hevc_stream_t* h) 
{
    memset(h->sh, 0, sizeof(hevc_slice_header_t));
    
    h->sh->collocated_from_l0_flag = 1;
}

static int NumDeltaPocs[MAX_NUM_SHORT_TERM_REF_PICS];
static int NumNegativePics[MAX_NUM_NEGATIVE_PICS];
static int NumPositivePics[MAX_NUM_POSITIVE_PICS];
static int DeltaPocS0[MAX_NUM_REF_PICS_L0][MAX_NUM_NEGATIVE_PICS];
static int UsedByCurrPicS0[MAX_NUM_REF_PICS_L0][MAX_NUM_NEGATIVE_PICS];
static int DeltaPocS1[MAX_NUM_REF_PICS_L1][MAX_NUM_POSITIVE_PICS];
static int UsedByCurrPicS1[MAX_NUM_REF_PICS_L1][MAX_NUM_POSITIVE_PICS];

//TODO: return right one
static int getNumPicTotalCurr(hevc_sps_t* sps, hevc_slice_header_t* sh)
{
    int i;
    int NumPicTotalCurr = 0;
    int CurrRpsIdx = sps->num_short_term_ref_pic_sets;
    
    if( sh->short_term_ref_pic_set_sps_flag )
        CurrRpsIdx = sh->short_term_ref_pic_set_idx;
    for( i = 0; i < NumNegativePics[ CurrRpsIdx ]; i++ )
        if( UsedByCurrPicS0[ CurrRpsIdx ][ i ] ) 
            NumPicTotalCurr++;
    for( i = 0; i < NumPositivePics[ CurrRpsIdx ]; i++) 
        if( UsedByCurrPicS1[ CurrRpsIdx ][ i ] )
            NumPicTotalCurr++;
    for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
        int UsedByCurrPicLt = 0;
        if( i < sh->num_long_term_sps )
            UsedByCurrPicLt = sps->used_by_curr_pic_lt_sps_flag[ sh->lt_idx_sps[ i ] ];
        else
            UsedByCurrPicLt = sh->used_by_curr_pic_lt_flag[ i ];
        if( UsedByCurrPicLt ) 
            NumPicTotalCurr++;
    }
    return NumPicTotalCurr;
}

static void updateNumDeltaPocs(hevc_st_ref_pic_set_t *st_ref_pic_set, int stRpsIdx) {
    int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
    
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        int i, j, dPoc;
        int deltaRps = ( 1 - 2 * st_ref_pic_set->delta_rps_sign ) * ( st_ref_pic_set->abs_delta_rps_minus1 + 1 );
        i = 0;
        for( j = NumPositivePics[ RefRpsIdx ] - 1; j >= 0; j-- ) {
            dPoc = DeltaPocS1[ RefRpsIdx ][ j ] + deltaRps;
            if( dPoc < 0 && st_ref_pic_set->use_delta_flag[ NumNegativePics[ RefRpsIdx ] + j ] ) {
                DeltaPocS0[ stRpsIdx ][ i ] = dPoc;
                UsedByCurrPicS0[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ NumNegativePics[ RefRpsIdx ] + j ];
            }
        }
        if( deltaRps < 0 && st_ref_pic_set->use_delta_flag[ NumDeltaPocs[ RefRpsIdx ] ] ) {
            DeltaPocS0[ stRpsIdx ][ i ] = deltaRps;
            UsedByCurrPicS0[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ NumDeltaPocs[ RefRpsIdx ] ];
        }
        for( j = 0; j < NumNegativePics[ RefRpsIdx ]; j++ ) {
            dPoc = DeltaPocS0[ RefRpsIdx ][ j ] + deltaRps;
            if( dPoc < 0 && st_ref_pic_set->use_delta_flag[ j ] ) {
                DeltaPocS0[ stRpsIdx ][ i ] = dPoc;
                UsedByCurrPicS0[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ j ];
            }
        }
        NumNegativePics[ stRpsIdx ] = i;
        i = 0;
        for( j = NumNegativePics[ RefRpsIdx ] - 1; j >= 0; j-- ) {
            dPoc = DeltaPocS0[ RefRpsIdx ][ j ] + deltaRps;
            if( dPoc > 0 && st_ref_pic_set->use_delta_flag[ j ] ) {
                DeltaPocS1[ stRpsIdx ][ i ] = dPoc;
                UsedByCurrPicS1[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ j ];
            }
        }
        if( deltaRps > 0 && st_ref_pic_set->use_delta_flag[ NumDeltaPocs[ RefRpsIdx ] ] ) {
            DeltaPocS1[ stRpsIdx ][ i ] = deltaRps;
            UsedByCurrPicS1[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ NumDeltaPocs[ RefRpsIdx ] ];
        }
        for( j = 0; j < NumPositivePics[ RefRpsIdx ]; j++) {
            dPoc = DeltaPocS1[ RefRpsIdx ][ j ] + deltaRps;
            if( dPoc > 0 && st_ref_pic_set->use_delta_flag[ NumNegativePics[ RefRpsIdx ] + j ] ) {
                DeltaPocS1[ stRpsIdx ][ i ] = dPoc;
                UsedByCurrPicS1[ stRpsIdx ][ i++ ] = st_ref_pic_set->used_by_curr_pic_flag[ NumNegativePics[ RefRpsIdx ] + j ];
            }
        }
        NumPositivePics[ stRpsIdx ] = i;
    } else {
        NumNegativePics[ stRpsIdx ] = st_ref_pic_set->num_negative_pics;
        NumPositivePics[ stRpsIdx ] = st_ref_pic_set->num_positive_pics;
    }
    
    NumDeltaPocs[ stRpsIdx ] = NumNegativePics[ stRpsIdx ] + NumPositivePics[ stRpsIdx ];
}

static int getSliceSegmentAddressBitLength(hevc_sps_t *sps) {
    int MinCbLog2SizeY = sps->log2_min_luma_coding_block_size_minus3 + 3;
    int CtbLog2SizeY = MinCbLog2SizeY + sps->log2_diff_max_min_luma_coding_block_size;
    int CtbSizeY = 1 << CtbLog2SizeY;
    int PicWidthInCtbsY = ceil( sps->pic_width_in_luma_samples * 1.0f / CtbSizeY );
    int PicHeightInCtbsY = ceil( sps->pic_height_in_luma_samples * 1.0f / CtbSizeY );
    int PicSizeInCtbsY = PicWidthInCtbsY * PicHeightInCtbsY;
    return ceil( log2( PicSizeInCtbsY ) ) ;
}

#end_preamble

#function_declarations


//7.3.1 NAL unit syntax
int structure(hevc_nal_unit)(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    int nal_size = size;
    int rbsp_size = size;
    uint8_t* rbsp_buf = (uint8_t*)calloc(1, rbsp_size);

    if( is_reading )
    {
        int rc = nal_to_rbsp(buf, &nal_size, rbsp_buf, &rbsp_size);

        if (rc < 0) { free(rbsp_buf); return -1; } // handle conversion error
    }

    if( is_writing )
    {
        rbsp_size = size*3/4; // NOTE this may have to be slightly smaller (3/4 smaller, worst case) in order to be guaranteed to fit
    }

    bs_t* b = bs_new(rbsp_buf, rbsp_size);
    value( forbidden_zero_bit, f(1, 0) );
    value( nal->nal_unit_type, u(6) );
    value( nal->nal_layer_id,  u(6) );
    value( nal->nal_temporal_id_plus1, u(3) );

    switch ( nal->nal_unit_type )
    {
        case HEVC_NAL_UNIT_TYPE_TRAIL_N:
        case HEVC_NAL_UNIT_TYPE_TRAIL_R:  
        case HEVC_NAL_UNIT_TYPE_TSA_N:
        case HEVC_NAL_UNIT_TYPE_TSA_R:
        case HEVC_NAL_UNIT_TYPE_STSA_N:
        case HEVC_NAL_UNIT_TYPE_STSA_R:
        case HEVC_NAL_UNIT_TYPE_RADL_N:
        case HEVC_NAL_UNIT_TYPE_RADL_R:
        case HEVC_NAL_UNIT_TYPE_RASL_N:
        case HEVC_NAL_UNIT_TYPE_RASL_R:
        case HEVC_NAL_UNIT_TYPE_BLA_W_LP:
        case HEVC_NAL_UNIT_TYPE_BLA_W_RADL:
        case HEVC_NAL_UNIT_TYPE_BLA_N_LP:
        case HEVC_NAL_UNIT_TYPE_IDR_W_RADL:
        case HEVC_NAL_UNIT_TYPE_IDR_N_LP:
        case HEVC_NAL_UNIT_TYPE_CRA_NUT:
            
            structure(hevc_slice_layer_rbsp)(h, b);
            break;

#ifdef HAVE_SEI
        case NAL_UNIT_TYPE_SEI:
            structure(hevc_sei_rbsp)(h, b);
            break;
#endif

        case HEVC_NAL_UNIT_TYPE_VPS_NUT: 
            structure(hevc_video_parameter_set_rbsp)(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_SPS_NUT: 
            structure(hevc_seq_parameter_set_rbsp)(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_PPS_NUT:   
            structure(hevc_pic_parameter_set_rbsp)(h, b);
            break;

        default:
            return -1;
    }

    if (bs_overrun(b)) { bs_free(b); free(rbsp_buf); return -1; }

    if( is_writing )
    {
        // now get the actual size used
        rbsp_size = bs_pos(b);

        int rc = rbsp_to_nal(rbsp_buf, &rbsp_size, buf, &nal_size);
        if (rc < 0) { bs_free(b); free(rbsp_buf); return -1; }
    }

    bs_free(b);
    free(rbsp_buf);

    return nal_size;
}

//7.3.2.1 Sequence parameter set RBSP syntax
void structure(hevc_video_parameter_set_rbsp)(hevc_stream_t* h, bs_t* b)
{
    int i, j;

    hevc_vps_t* vps = h->vps;
    if( is_reading )
    {
        memset(vps, 0, sizeof(hevc_vps_t));
    }
 
    value( vps->vps_video_parameter_set_id,               u(4) );
    value( vps->vps_base_layer_internal_flag,             u1 );
    value( vps->vps_base_layer_available_flag,            u1 );
    value( vps->vps_max_layers_minus1,                    u(6) );
    value( vps->vps_max_sub_layers_minus1,                u(3) );
    value( vps->vps_temporal_id_nesting_flag,             u1 );
    value( vps_reserved_0xffff_16bits,                    f(16, 0xffff) );
    
    structure(hevc_profile_tier_level)(&vps->ptl, b, 1, vps->vps_max_sub_layers_minus1); 
    
    value( vps->vps_sub_layer_ordering_info_present_flag, u1 );
    for( i = ( vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1 ); 
            i <= vps->vps_max_sub_layers_minus1; i++ ) {
        value( vps->vps_max_dec_pic_buffering_minus1[ i ],ue );        
        value( vps->vps_max_num_reorder_pics[ i ],        ue );        
        value( vps->vps_max_latency_increase_plus1[ i ],  ue );        
    }
    value( vps->vps_max_layer_id,                         u(6) );
    value( vps->vps_num_layer_sets_minus1,                ue );
    for( i = 1; i <= vps->vps_num_layer_sets_minus1; i++ )
        for( j = 0; j <= vps->vps_max_layer_id; j++ ) {
            value( vps->layer_id_included_flag[ i ][ j ], u1 );
        }
    value( vps->vps_timing_info_present_flag,             u1 );
    if( vps->vps_timing_info_present_flag ) {
        value( vps->vps_num_units_in_tick,                u(32) );
        value( vps->vps_time_scale,                       u(32) );
        value( vps->vps_poc_proportional_to_timing_flag,  u1 );
        if( vps->vps_poc_proportional_to_timing_flag ) {
            value( vps->vps_num_ticks_poc_diff_one_minus1,ue );
        }
        value( vps->vps_num_hrd_parameters,               ue );
        for( i = 0; i < vps->vps_num_hrd_parameters; i++ ) {
            value( vps->hrd_layer_set_idx[ i ],           ue );
            if (i > 0) {
                value( vps->cprms_present_flag[ i ],      u1 );
            }
            structure(hevc_hrd_parameters)(&vps->hrd[i], b,
                                           vps->cprms_present_flag[ i ],
                                           vps->vps_max_sub_layers_minus1);
        }
    }
    value( vps->vps_extension_flag,                       u1 );
    //TODO: support extension data
    //if (vps->vps_extension_flag)    

    structure(hevc_rbsp_trailing_bits)(b);
}

//7.3.2.2 Sequence parameter set RBSP syntax
void structure(hevc_seq_parameter_set_rbsp)(hevc_stream_t* h, bs_t* b)
{
    int i;

    hevc_sps_t* sps = h->sps;
    if( is_reading )
    {
        memset(sps, 0, sizeof(hevc_sps_t));
    }
 
    value( sps->sps_video_parameter_set_id,               u(4) );
    value( sps->sps_max_sub_layers_minus1,                u(3) );
    value( sps->sps_temporal_id_nesting_flag,             u1 );
    structure(hevc_profile_tier_level)(&sps->ptl, b, 1, sps->sps_max_sub_layers_minus1); 
    value( sps->sps_seq_parameter_set_id,                 ue );
    value( sps->chroma_format_idc,                        ue );
    if( sps->chroma_format_idc == 3 ) {
        value( sps->separate_colour_plane_flag,           u1 );
    }
    value( sps->pic_width_in_luma_samples,                ue );
    value( sps->pic_height_in_luma_samples,               ue );
    value( sps->conformance_window_flag,                  u1 );
    if( sps->conformance_window_flag ) {
        value( sps->conf_win_left_offset,                 ue );
        value( sps->conf_win_right_offset,                ue );
        value( sps->conf_win_top_offset,                  ue );
        value( sps->conf_win_bottom_offset,               ue );
    }
    value( sps->bit_depth_luma_minus8,                    ue );
    value( sps->bit_depth_chroma_minus8,                  ue );
    value( sps->log2_max_pic_order_cnt_lsb_minus4,        ue );
    value( sps->sps_sub_layer_ordering_info_present_flag, u1 );
    for( i = ( sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1 ); 
            i <= sps->sps_max_sub_layers_minus1; i++ ) {
        value( sps->sps_max_dec_pic_buffering_minus1 [ i ], ue );
        value( sps->sps_max_num_reorder_pics [ i ],         ue );
        value( sps->sps_max_latency_increase_plus1 [ i ],   ue );
    }
    value( sps->log2_min_luma_coding_block_size_minus3,      ue );
    value( sps->log2_diff_max_min_luma_coding_block_size,    ue );
    value( sps->log2_min_luma_transform_block_size_minus2,   ue );
    value( sps->log2_diff_max_min_luma_transform_block_size, ue );
    value( sps->max_transform_hierarchy_depth_inter,         ue );
    value( sps->max_transform_hierarchy_depth_intra,         ue );
    value( sps->scaling_list_enabled_flag,                   u1 );
    
    if( sps->scaling_list_enabled_flag ) {
        value( sps->sps_scaling_list_data_present_flag,      u1 );
        if( sps->sps_scaling_list_data_present_flag ) {
            structure(hevc_scaling_list_data)(&sps->scaling_list_data, b); 
        }
    }
    
    value( sps->amp_enabled_flag,                                 u1 );
    value( sps->sample_adaptive_offset_enabled_flag,              u1 );
    value( sps->pcm_enabled_flag,                                 u1 );
    if( sps->pcm_enabled_flag ) {
        value( sps->pcm_sample_bit_depth_luma_minus1,             u(4) );
        value( sps->pcm_sample_bit_depth_chroma_minus1,           u(4) );
        value( sps->log2_min_pcm_luma_coding_block_size_minus3,   ue );
        value( sps->log2_diff_max_min_pcm_luma_coding_block_size, ue );
        value( sps->pcm_loop_filter_disabled_flag,                u1 );
    }
    value( sps->num_short_term_ref_pic_sets, ue );
    for( i = 0; i < sps->num_short_term_ref_pic_sets; i++) {
        structure(hevc_st_ref_pic_set)(&sps->st_ref_pic_set[i], b, i, sps->num_short_term_ref_pic_sets);
    }
    
    value( sps->long_term_ref_pics_present_flag, u1 );
    if( sps->long_term_ref_pics_present_flag ) {
        value( sps->num_long_term_ref_pics_sps, ue );
        for( i = 0; i < sps->num_long_term_ref_pics_sps; i++ ) {
            value( sps->lt_ref_pic_poc_lsb_sps[ i ], u( sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ) );
            value( sps->used_by_curr_pic_lt_sps_flag[ i ], u1 );
        }
    }
    value( sps->sps_temporal_mvp_enabled_flag, u1 );
    value( sps->strong_intra_smoothing_enabled_flag, u1 );
    value( sps->vui_parameters_present_flag, u1 );
    if( sps->vui_parameters_present_flag ) {
        structure(hevc_vui_parameters)(sps, b);
    }
    value( sps->sps_extension_present_flag, u1 );
    
    if( sps->sps_extension_present_flag ) {
        value( sps->sps_range_extension_flag, u1 );
        value( sps->sps_multilayer_extension_flag, u1 );
        value( sps->sps_3d_extension_flag, u1 );
        value( sps->sps_extension_5bits, u(5) );
    }
    if( sps->sps_range_extension_flag ) {
        structure(hevc_sps_range_extension)( &sps->sps_range_ext, b);
    }
    
    if( is_reading )
    {
        memcpy(h->sps_table[sps->sps_seq_parameter_set_id], h->sps, sizeof(hevc_sps_t));
    }
}

//7.3.2.2.2 Sequence parameter set range extension syntax
void structure(hevc_sps_range_extension)(hevc_sps_range_ext_t* sps_range_ext, bs_t* b)
{
    value( sps_range_ext->transform_skip_rotation_enabled_flag, u1);
    value( sps_range_ext->transform_skip_context_enabled_flag, u1);
    value( sps_range_ext->implicit_rdpcm_enabled_flag, u1);
    value( sps_range_ext->explicit_rdpcm_enabled_flag, u1);
    value( sps_range_ext->extended_precision_processing_flag, u1);
    value( sps_range_ext->intra_smoothing_disabled_flag, u1);
    value( sps_range_ext->high_precision_offsets_enabled_flag, u1);
    value( sps_range_ext->persistent_rice_adaptation_enabled_flag, u1);
    value( sps_range_ext->cabac_bypass_alignment_enabled_flag, u1);
}


//7.3.2.3 Picture parameter set RBSP syntax
void structure(hevc_pic_parameter_set_rbsp)(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_pps_t* pps = h->pps;
    if( is_reading )
    {
        memset(pps, 0, sizeof(hevc_pps_t));
    }

    value( pps->pic_parameter_set_id, ue);
    value( pps->seq_parameter_set_id, ue );
    value( pps->dependent_slice_segments_enabled_flag, u1 );
    value( pps->output_flag_present_flag, u1 );
    value( pps->num_extra_slice_header_bits, u( 3 ) );
    value( pps->sign_data_hiding_enabled_flag, u1 );
    value( pps->cabac_init_present_flag, u1 );
    value( pps->num_ref_idx_l0_default_active_minus1, ue );
    value( pps->num_ref_idx_l1_default_active_minus1, ue );
    value( pps->init_qp_minus26, se );
    value( pps->constrained_intra_pred_flag, u1 );
    value( pps->transform_skip_enabled_flag, u1 );
    value( pps->cu_qp_delta_enabled_flag, u1 );
    if( pps->cu_qp_delta_enabled_flag ) {
        value( pps->diff_cu_qp_delta_depth, ue);
    }
    value( pps->pps_cb_qp_offset, se );
    value( pps->pps_cr_qp_offset, se );
    value( pps->pps_slice_chroma_qp_offsets_present_flag, u1 );
    value( pps->weighted_pred_flag, u1 );
    value( pps->weighted_bipred_flag, u1 );
    value( pps->transquant_bypass_enabled_flag, u1 );
    value( pps->tiles_enabled_flag, u1 );
    value( pps->entropy_coding_sync_enabled_flag, u1 );
    if( pps->tiles_enabled_flag ) {
        value( pps->num_tile_columns_minus1, ue );
        value( pps->num_tile_rows_minus1, ue );
        value( pps->uniform_spacing_flag, u1 );
        if( !pps->uniform_spacing_flag ) {
            for( i = 0; i < pps->num_tile_columns_minus1; i++ ) {
                value( pps->column_width_minus1[ i ], ue );
            }
            for( i = 0; i < pps->num_tile_rows_minus1; i++ ) {
                value( pps->row_height_minus1[ i ], ue );
            }
        }
        value( pps->loop_filter_across_tiles_enabled_flag, u1 );
    }
    value( pps->pps_loop_filter_across_slices_enabled_flag, u1 );
    value( pps->deblocking_filter_control_present_flag, u1 );
    if( pps->deblocking_filter_control_present_flag ) {
        value( pps->deblocking_filter_override_enabled_flag, u1 );
        value( pps->pps_deblocking_filter_disabled_flag, u1 );
        if( pps->pps_deblocking_filter_disabled_flag ) {
            value( pps->pps_beta_offset_div2, se );
            value( pps->pps_tc_offset_div2, se );
        }
    }
    value( pps->pps_scaling_list_data_present_flag, u1 );
    if( pps->pps_scaling_list_data_present_flag ) {
        structure(hevc_scaling_list_data)(&pps->scaling_list_data, b);
    }
    value( pps->lists_modification_present_flag, u1 );
    value( pps->log2_parallel_merge_level_minus2, ue );
    value( pps->slice_segment_header_extension_present_flag, u1 );
    value( pps->pps_extension_present_flag, u1 );
    if( pps->pps_extension_present_flag ) {
        value( pps->pps_range_extension_flag, u1 );
        value( pps->pps_multilayer_extension_flag, u1 );
        value( pps->pps_3d_extension_flag, u1 );
        value( pps->pps_extension_5bits, u1 );
    }
    if( pps->pps_range_extension_flag ) {
        structure(hevc_pps_range_extension)( pps, b);
    }

    structure(hevc_rbsp_trailing_bits)(b);

    if( is_reading )
    {
        memcpy(h->pps_table[pps->pic_parameter_set_id], h->pps, sizeof(hevc_pps_t));
    }
}

//7.3.2.3.2 Picture parameter set range extension syntax
void structure(hevc_pps_range_extension)(hevc_pps_t* pps, bs_t* b)
{
    hevc_pps_range_ext_t *pps_range_ext = &pps->pps_range_ext;;
    if( pps->transform_skip_enabled_flag ) {
        value( pps_range_ext->log2_max_transform_skip_block_size_minus2, ue);
    }
    value( pps_range_ext->cross_component_prediction_enabled_flag, u1);
    value( pps_range_ext->chroma_qp_offset_list_enabled_flag, u1);
    if( pps_range_ext->chroma_qp_offset_list_enabled_flag ) {
        value( pps_range_ext->diff_cu_chroma_qp_offset_depth, ue);
        value( pps_range_ext->chroma_qp_offset_list_len_minus1, ue);
        for( int i = 0; i <= pps_range_ext->chroma_qp_offset_list_len_minus1; i++ ) {
            value( pps_range_ext->cb_qp_offset_list[ i ], se);
            value( pps_range_ext->cr_qp_offset_list[ i ], se);
        }
    }
    value( pps_range_ext->log2_sao_offset_scale_luma, ue);
    value( pps_range_ext->log2_sao_offset_scale_chroma, ue);
}

#ifdef HAVE_SEI
//7.3.2.4 Supplemental enhancement information RBSP syntax
void structure(sei_rbsp)(hevc_stream_t* h, bs_t* b)
{
    if( is_reading )
    {
        for( int i = 0; i < h->num_seis; i++ ) {
            sei_free(h->seis[i]);
        }
    
        h->num_seis = 0;
        do {
            h->num_seis++;
            h->seis = (sei_t**)realloc(h->seis, h->num_seis * sizeof(sei_t*));
            h->seis[h->num_seis - 1] = sei_new();
            h->sei = h->seis[h->num_seis - 1];
            structure(sei_message)(h, b);
        } while( more_rbsp_data(h, b) );
    }

    if( is_writing )
    {
        for (int i = 0; i < h->num_seis; i++) {
            h->sei = h->seis[i];
            structure(sei_message)(h, b);
        }
        h->sei = NULL;
    }

    structure(hevc_rbsp_trailing_bits)(b);
}

//7.3.5 Supplemental enhancement information message syntax
void structure(sei_message)(hevc_stream_t* h, bs_t* b)
{
    if( is_writing )
    {
        _write_ff_coded_number(b, h->sei->payloadType);
        _write_ff_coded_number(b, h->sei->payloadSize);
    }
    if( is_reading )
    {
        h->sei->payloadType = _read_ff_coded_number(b);
        h->sei->payloadSize = _read_ff_coded_number(b);
    }
    structure(sei_payload)( h, b, h->sei->payloadType, h->sei->payloadSize );
}
#endif

//7.3.2.5 Access unit delimiter RBSP syntax
void structure(hevc_access_unit_delimiter_rbsp)(hevc_stream_t* h, bs_t* b)
{
    value( h->aud->primary_pic_type, u(3) );
    structure(hevc_rbsp_trailing_bits)(b);
}

//7.3.2.6 End of sequence RBSP syntax
void structure(hevc_end_of_seq_rbsp)()
{
}

//7.3.2.7 End of bitstream RBSP syntax
void structure(end_of_bitstream_rbsp)()
{
}

//7.3.2.8 Filler data RBSP syntax
void structure(filler_data_rbsp)(bs_t* b)
{
    while( bs_next_bits(b, 8) == 0xFF )
    {
        value( ff_byte, f(8, 0xFF) );
    }
    structure(hevc_rbsp_trailing_bits)(b);
}

//7.3.2.9 Slice segment layer RBSP syntax
void structure(hevc_slice_layer_rbsp)(hevc_stream_t* h,  bs_t* b)
{
    structure(hevc_slice_header)(h, b);
    hevc_slice_data_rbsp_t* slice_data = h->slice_data;

    if ( slice_data != NULL )
    {
        if ( slice_data->rbsp_buf != NULL ) free( slice_data->rbsp_buf ); 
        uint8_t *sptr = b->p + (!!b->bits_left); // CABAC-specific: skip alignment bits, if there are any
        slice_data->rbsp_size = b->end - sptr;
        
        slice_data->rbsp_buf = (uint8_t*)malloc(slice_data->rbsp_size);
        memcpy( slice_data->rbsp_buf, sptr, slice_data->rbsp_size );
    }

    //structure(hevc_slice_data)(h, b); /* all categories of slice_data( ) syntax */
    structure(hevc_rbsp_slice_trailing_bits)( b );
}

//7.3.2.10 RBSP slice trailing bits syntax
void structure(hevc_rbsp_slice_trailing_bits)(bs_t* b)
{
    structure(hevc_rbsp_trailing_bits)(b);
    //while( more_rbsp_trailing_data(b) )
    //{
    //    value( cabac_zero_word, f(16, 0x0000) );
    //}
}

//7.3.2.11 RBSP trailing bits syntax
void structure(hevc_rbsp_trailing_bits)(bs_t* b)
{
    value( rbsp_stop_one_bit, f(1, 1) );

    while( !bs_byte_aligned(b) )
    {
        value( rbsp_alignment_zero_bit, f(1, 0) );
    }
}

//7.3.2.12 Byte alignment syntax
void structure(hevc_byte_alignment)(bs_t* b)
{
    value( alignment_bit_equal_to_one, f(1, 1) );

    while( !bs_byte_aligned(b) )
    {
        value( alignment_bit_equal_to_zero, f(1, 0) );
    }
}

//7.3.3 Profile, tier and level syntax
void structure(hevc_profile_tier_level)(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1)
{
    int i, j;
    if( profilePresentFlag ) {
        value( ptl->general_profile_space, u(2) );
        value( ptl->general_tier_flag,     u1 );
        value( ptl->general_profile_idc,   u(5) );
        for( i = 0; i < 32; i++ ) {
            value( ptl->general_profile_compatibility_flag[ i ],  u1 );
        }
        value( ptl->general_progressive_source_flag,            u1 );
        value( ptl->general_interlaced_source_flag,             u1 );
        value( ptl->general_non_packed_constraint_flag,         u1 );
        value( ptl->general_frame_only_constraint_flag,         u1 );
        if( ptl->general_profile_idc == 4 || ptl->general_profile_compatibility_flag[ 4 ] || 
            ptl->general_profile_idc == 5 || ptl->general_profile_compatibility_flag[ 5 ] || 
            ptl->general_profile_idc == 6 || ptl->general_profile_compatibility_flag[ 6 ] || 
            ptl->general_profile_idc == 7 || ptl->general_profile_compatibility_flag[ 7 ] ) {
                
            value( ptl->general_max_12bit_constraint_flag,        u1 );
            value( ptl->general_max_10bit_constraint_flag,        u1 );
            value( ptl->general_max_8bit_constraint_flag,         u1 );
            value( ptl->general_max_422chroma_constraint_flag,    u1 );
            value( ptl->general_max_420chroma_constraint_flag,    u1 );
            value( ptl->general_max_monochrome_constraint_flag,   u1 );
            value( ptl->general_intra_constraint_flag,            u1 );
            value( ptl->general_one_picture_only_constraint_flag, u1 );
            value( ptl->general_lower_bit_rate_constraint_flag,   u1 );
            value( general_reserved_zero_34bits,                  f(34, 0) );
        } else {
            value( general_reserved_zero_43bits,                  f(43, 0) );
        }
        if( ( ptl->general_profile_idc >= 1 && ptl->general_profile_idc <= 5 ) ||
              ptl->general_profile_compatibility_flag[ 1 ] ||
              ptl->general_profile_compatibility_flag[ 2 ] ||
              ptl->general_profile_compatibility_flag[ 3 ] ||
              ptl->general_profile_compatibility_flag[ 4 ] ||
              ptl->general_profile_compatibility_flag[ 5 ] ) {

            value( ptl->general_inbld_flag,                       u1 );
        } else {
            value( general_reserved_zero_bit,                     f(1, 0) );
        }
        value( ptl->general_level_idc, u8 );
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) {
            value( ptl->sub_layer_profile_present_flag[ i ],      u1 );
            value( ptl->sub_layer_level_present_flag[ i ],        u1 );
        }
        if( maxNumSubLayersMinus1 > 0 ) {
            for( i = maxNumSubLayersMinus1; i < 8; i++ ) {
                value( reserved_zero_xxbits,                      f(2, 0) );
            }
        }
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) { 
            if( ptl->sub_layer_profile_present_flag[ i ] ) {
                value( ptl->sub_layer_profile_space[ i ],         u(2) );
                value( ptl->sub_layer_tier_flag[ i ],             u1 );
                value( ptl->sub_layer_profile_idc[ i ],           u(5) );
                for( j = 0; j < 32; j++ ) {
                    value( ptl->sub_layer_profile_compatibility_flag[ i ][ j ], u(1));
                }
                value( ptl->sub_layer_progressive_source_flag[ i ],             u1 );
                value( ptl->sub_layer_interlaced_source_flag[ i ],              u1 );
                value( ptl->sub_layer_non_packed_constraint_flag[ i ],          u1 );
                value( ptl->sub_layer_frame_only_constraint_flag[ i ],          u1 );
                if( ptl->sub_layer_profile_idc[ i ] == 4 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 4 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 5 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 5 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 6 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 6 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 7 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 7 ] ) {
                    value( ptl->sub_layer_max_12bit_constraint_flag[ i ],           u1 );
                    value( ptl->sub_layer_max_10bit_constraint_flag[ i ],           u1 );
                    value( ptl->sub_layer_max_8bit_constraint_flag[ i ],            u1 );
                    value( ptl->sub_layer_max_422chroma_constraint_flag[ i ],       u1 );
                    value( ptl->sub_layer_max_420chroma_constraint_flag[ i ],       u1 );
                    value( ptl->sub_layer_max_monochrome_constraint_flag[ i ],      u1 );
                    value( ptl->sub_layer_intra_constraint_flag[ i ],               u1 );
                    value( ptl->sub_layer_one_picture_only_constraint_flag[ i ],    u1 );
                    value( ptl->sub_layer_lower_bit_rate_constraint_flag[ i ],      u1 );
                    value( sub_layer_reserved_zero_34bits,                          f(34, 0) );
                } else {
                    value( sub_layer_reserved_zero_43bits,                          f(43, 0) );
                }
            
                if( ( ptl->sub_layer_profile_idc[ i ] >= 1 && ptl->sub_layer_profile_idc[ i ] <= 5 ) ||
                   ptl->sub_layer_profile_compatibility_flag[ 1 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 2 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 3 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 4 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 5 ] ) {
                    value( ptl->sub_layer_inbld_flag[ i ],                          u1 );
                } else {
                    value( sub_layer_reserved_zero_bit,                             f(1, 0) );
                }
            }
            if( ptl->sub_layer_level_present_flag[ i ] ) {
                value( ptl->sub_layer_level_idc[ i ],                               u8 );
            }
        }
    }
}

//7.3.4 Scaling list data syntax
void structure(hevc_scaling_list_data)( hevc_scaling_list_data_t *sld, bs_t* b )
{
    int nextCoef, coefNum;
    for( int sizeId = 0; sizeId < 4; sizeId++ )
        for( int matrixId = 0; matrixId < 6; matrixId += ( sizeId == 3 ) ? 3 : 1 ) {
            value( sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ], u1 );
            if( !sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] ) {
                value( sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ], ue );
            } else {
                nextCoef = 8;
                coefNum=MIN(64, (1 << (4+(sizeId << 1))));
                if( sizeId > 1 ) {
                    value( sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ], se );
                }
 
                for( int i = 0; i < coefNum; i++) {
                    value( sld->scaling_list_delta_coef[ sizeId ][ matrixId ], se );
                }
            }
 
        }
}

//7.3.6 Slice header syntax
void structure(hevc_slice_header)(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    if( is_reading )
    {
        init_slice_hevc(h);
    }

    hevc_nal_t* nal = h->nal;

    value( sh->first_slice_segment_in_pic_flag, u1 );
    if( nal->nal_unit_type >= HEVC_NAL_UNIT_TYPE_BLA_W_LP && 
       nal->nal_unit_type <= HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23) {
            value( sh->no_output_of_prior_pics_flag, u1 );
    }
    value( sh->pic_parameter_set_id, ue );
    
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];

    //set default value
    sh->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    sh->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;

    if( !sh->first_slice_segment_in_pic_flag ) {
        if( pps->dependent_slice_segments_enabled_flag ) {
            value( sh->dependent_slice_segment_flag, u1 );
        }
        value( sh->slice_segment_address, u( getSliceSegmentAddressBitLength( sps ) ) );
    }
    
    if( !sh->dependent_slice_segment_flag ) {
        for( i = 0; i < pps->num_extra_slice_header_bits; i++ ) {
            value( slice_reserved_flag,                              f(1, 1) );
        }
        value( sh->slice_type, ue );
        if( pps->output_flag_present_flag ) {
            value( sh->pic_output_flag, u1 );
        }
        if( sps->separate_colour_plane_flag == 1 ) {
            value( sh->colour_plane_id, u(2) );
        }
        if( nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_W_RADL &&
            nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_N_LP) {
            value( sh->slice_pic_order_cnt_lsb, u( sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ) );
            value( sh->short_term_ref_pic_set_sps_flag, u1 );
            if( !sh->short_term_ref_pic_set_sps_flag ) {
                structure(hevc_st_ref_pic_set)( &sh->st_ref_pic_set, b, sps->num_short_term_ref_pic_sets, sps->num_short_term_ref_pic_sets );
            } else if( sps->num_short_term_ref_pic_sets > 1 ) {
                value( sh->short_term_ref_pic_set_idx, u( ceil( log2( sps->num_short_term_ref_pic_sets ) ) ) );
            }
            if( sps->long_term_ref_pics_present_flag ) {
                if( sps->num_long_term_ref_pics_sps > 0 ) {
                    value( sh->num_long_term_sps, ue );
                }
                value( sh->num_long_term_pics, ue );
                for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
                    if( i < sh->num_long_term_sps ) {
                        if( sps->num_long_term_ref_pics_sps > 1 ) {
                            value( sh->lt_idx_sps[ i ], u( ceil( log2( sps->num_long_term_ref_pics_sps ) ) ) );
                        }
                    } else {
                        value( sh->poc_lsb_lt[ i ], u( sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ) );
                        value( sh->used_by_curr_pic_lt_flag[ i ], u1 );
                    }
                    value( sh->delta_poc_msb_present_flag[ i ], u1 );
                    if( sh->delta_poc_msb_present_flag[ i ]) {
                        value( sh->delta_poc_msb_cycle_lt[ i ], ue );
                    }
                }
            }
            if( sps->sps_temporal_mvp_enabled_flag ) {
                value( sh->slice_temporal_mvp_enabled_flag, u1 );
            }
        }
        if( sps->sample_adaptive_offset_enabled_flag ) {
            value( sh->slice_sao_luma_flag, u1 );
            int ChromaArrayType = 0;
            if( sps->separate_colour_plane_flag == 0) {
                ChromaArrayType = sps->chroma_format_idc;
            }
            if( ChromaArrayType != 0 ) {
                value( sh->slice_sao_chroma_flag, u1 );
            }
        }
        if( sh->slice_type == HEVC_SLICE_TYPE_P || sh->slice_type == HEVC_SLICE_TYPE_B ){
            value( sh->num_ref_idx_active_override_flag, u1 );
            if( sh->num_ref_idx_active_override_flag ) {
                value( sh->num_ref_idx_l0_active_minus1, ue );
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    value( sh->num_ref_idx_l1_active_minus1, ue );
                }
            }
            if( pps->lists_modification_present_flag && getNumPicTotalCurr( sps, sh ) > 1 ) {
                structure(hevc_ref_pic_lists_modification)( h, b );
            }
            
            if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                value( sh->mvd_l1_zero_flag, u1 );
            }
            if( pps->cabac_init_present_flag ) {
                value( sh->cabac_init_flag, u1 );
            }
            if( sh->slice_temporal_mvp_enabled_flag ) {
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    value( sh->collocated_from_l0_flag, u1 );
                }
                if( ( sh->collocated_from_l0_flag && sh->num_ref_idx_l0_active_minus1 > 0 ) ||
                   ( !sh->collocated_from_l0_flag && sh->num_ref_idx_l1_active_minus1 > 0 ) ) {
                    value( sh->collocated_ref_idx, ue );
                }
            }
            if( ( pps->weighted_pred_flag && sh->slice_type == HEVC_SLICE_TYPE_P ) || 
                ( pps->weighted_bipred_flag && sh->slice_type == HEVC_SLICE_TYPE_B ) ) {
                structure(hevc_pred_weight_table)( h, b );
            }
            value( sh->five_minus_max_num_merge_cand, ue );
        }
        value( sh->slice_qp_delta, se );
        if( pps->pps_slice_chroma_qp_offsets_present_flag ) {
            value( sh->slice_cb_qp_offset, se );
            value( sh->slice_cr_qp_offset, se );
        }
        if( pps->pps_range_ext.chroma_qp_offset_list_enabled_flag ) {
            value( sh->cu_chroma_qp_offset_enabled_flag, u1 );
        }
        if( pps->deblocking_filter_override_enabled_flag ) {
            value( sh->deblocking_filter_override_flag, u1 );
        }
        if( sh->deblocking_filter_override_flag ) {
            value( sh->slice_deblocking_filter_disabled_flag, u1 );
            if( !sh->slice_deblocking_filter_disabled_flag ) {
                value( sh->slice_beta_offset_div2, se );
                value( sh->slice_tc_offset_div2, se );
            }
        }
        if( pps->pps_loop_filter_across_slices_enabled_flag &&
           ( sh->slice_sao_luma_flag || sh->slice_sao_chroma_flag || !sh->slice_deblocking_filter_disabled_flag ) ) {
            value( sh->slice_loop_filter_across_slices_enabled_flag, u1 );
        }
    }
    if( pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag ) {
        value( sh->num_entry_point_offsets, ue );
        if( sh->num_entry_point_offsets > 0 ) {
            value( sh->offset_len_minus1, ue );
            for( i = 0; i < sh->num_entry_point_offsets; i++ ) {
                value( sh->entry_point_offset_minus1[ i ], u( sh->offset_len_minus1 + 1 ) );
            }
        }
    }
    if( pps->slice_segment_header_extension_present_flag ) {
        value( sh->slice_segment_header_extension_length, ue );
        //TODO: support header extension,
        for( i = 0; i < sh->slice_segment_header_extension_length; i++) {
            value( slice_segment_header_extension_data_byte, f(8, 0) );
        }
    }
    structure(hevc_byte_alignment)( b );
}

//7.3.6.2 Reference picture list reordering syntax
void structure(hevc_ref_pic_lists_modification)(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    
    value( sh->rpld.ref_pic_list_modification_flag_l0, u1 );
    if( sh->rpld.ref_pic_list_modification_flag_l0 ) {
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            value( sh->rpld.list_entry_l0[ i ], u( ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) ) );
        }
    }
    
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        value( sh->rpld.ref_pic_list_modification_flag_l1, 1 );
        if( sh->rpld.ref_pic_list_modification_flag_l1 ) {
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                value( sh->rpld.list_entry_l1[ i ], u( ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) ) );
            }
        }
    }
}

//7.3.6.3 Prediction weight table syntax
void structure(hevc_pred_weight_table)(hevc_stream_t* h, bs_t* b)
{
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    hevc_pred_weight_table_t *pwt = &sh->pwt;

    int i, j;

    value( pwt->luma_log2_weight_denom, ue );
    
    int ChromaArrayType = 0;
    if( sps->separate_colour_plane_flag == 0) {
        ChromaArrayType = sps->chroma_format_idc;
    }

    if( ChromaArrayType != 0 ) {
        value( pwt->delta_chroma_log2_weight_denom, se );
    }
    
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        value( pwt->luma_weight_l0_flag[i], u1 );
    }
    if( ChromaArrayType != 0 )
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            value( pwt->chroma_weight_l0_flag[i], u1 );
        }
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        if( pwt->luma_weight_l0_flag[i] ) {
            value( pwt->delta_luma_weight_l0[i], se );
            value( pwt->luma_offset_l0[i], se );
        }
        if( pwt->chroma_weight_l0_flag[i] ) {
            for( j =0; j < 2; j++ ) {
                value( pwt->delta_chroma_weight_l0[i][j], se );
                value( pwt->delta_chroma_offset_l0[i][j], se );
            }
        }
    }
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            value( pwt->luma_weight_l1_flag[i], u1 );
        }
        if( ChromaArrayType != 0 )
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                value( pwt->chroma_weight_l1_flag[i], u1 );
            }
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            if( pwt->luma_weight_l1_flag[i] ) {
                value( pwt->delta_luma_weight_l1[i], se );
                value( pwt->luma_offset_l1[i], se );
            }
            if( pwt->chroma_weight_l1_flag[i] ) {
                for( j =0; j < 2; j++ ) {
                    value( pwt->delta_chroma_weight_l1[i][j], se );
                    value( pwt->delta_chroma_offset_l1[i][j], se );
                }
            }
        }        
    }
}

//7.3.7 Short-term reference picture set syntax
void structure(hevc_st_ref_pic_set)( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets )
{
    int i, j;
    
    if( stRpsIdx != 0 ) {
        value( st_ref_pic_set->inter_ref_pic_set_prediction_flag,     u1);
    }
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        if( stRpsIdx == num_short_term_ref_pic_sets ) {
            value( st_ref_pic_set->delta_idx_minus1,                  ue );
        }
        value( st_ref_pic_set->delta_rps_sign,                        u1 );
        value( st_ref_pic_set->abs_delta_rps_minus1,                  ue );
        
        int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
        
        for( j = 0; j <= NumDeltaPocs[ RefRpsIdx ]; j++ ) {
            value( st_ref_pic_set->used_by_curr_pic_flag[ j ],        u1 );
            if( !st_ref_pic_set->used_by_curr_pic_flag[ j ] ) {
                value( st_ref_pic_set->use_delta_flag[ j ],           u1 );
            }
        }
    } else {
        value( st_ref_pic_set->num_negative_pics,                     ue );
        value( st_ref_pic_set->num_positive_pics,                     ue );
        for( i = 0; i < st_ref_pic_set->num_negative_pics; i++ ) {
            value( st_ref_pic_set->delta_poc_s0_minus1[ i ],          ue );
            value( st_ref_pic_set->used_by_curr_pic_s0_flag[ i ],     u1 );
            
            //update derived field
            UsedByCurrPicS0[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s0_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS0[ stRpsIdx ][ i ] = -1 * ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            } else {
                DeltaPocS0[ stRpsIdx ][ i ] = DeltaPocS0[ stRpsIdx ][ i - 1 ] - ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            }
        }
        for( i = 0; i < st_ref_pic_set->num_positive_pics; i++ ) {
            value( st_ref_pic_set->delta_poc_s1_minus1[ i ],          ue );
            value( st_ref_pic_set->used_by_curr_pic_s1_flag[ i ],     u1 );
        
            //update derived field
            UsedByCurrPicS1[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s1_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS1[ stRpsIdx ][ i ] = st_ref_pic_set->delta_poc_s1_minus1[ i ] + 1;
            } else {
                DeltaPocS1[ stRpsIdx ][ i ] = DeltaPocS1[ stRpsIdx ][ i - 1 ] + ( st_ref_pic_set->delta_poc_s1_minus1[ i ] + 1 );
            }
        }
    }
    updateNumDeltaPocs( st_ref_pic_set, stRpsIdx);
}

//Appendix E.2.1 VUI parameters syntax
void structure(hevc_vui_parameters)(hevc_sps_t* sps, bs_t* b)
{
    hevc_vui_t* vui = &sps->vui;
    value( vui->aspect_ratio_info_present_flag, u1 );
    if( vui->aspect_ratio_info_present_flag )
    {
        value( vui->aspect_ratio_idc, u8 );
        if( vui->aspect_ratio_idc == SAR_Extended )
        {
            value( vui->sar_width, u(16) );
            value( vui->sar_height, u(16) );
        }
    }
    value( vui->overscan_info_present_flag, u1 );
    if( vui->overscan_info_present_flag ) {
        value( vui->overscan_appropriate_flag, u1 );
    }
    value( vui->video_signal_type_present_flag, u1 );
    if( vui->video_signal_type_present_flag ) {
        value( vui->video_format, u(3) );
        value( vui->video_full_range_flag, u1 );
        value( vui->colour_description_present_flag, u1 );
        if( vui->colour_description_present_flag ) {
            value( vui->colour_primaries, u8 );
            value( vui->transfer_characteristics, u8 );
            value( vui->matrix_coefficients, u8 );
        }
    }
    value( vui->chroma_loc_info_present_flag, u1 );
    if( vui->chroma_loc_info_present_flag ) {
        value( vui->chroma_sample_loc_type_top_field, ue );
        value( vui->chroma_sample_loc_type_bottom_field, ue );
    }
    
    value( vui->neutral_chroma_indication_flag, u1 );
    value( vui->field_seq_flag, u1 );
    value( vui->frame_field_info_present_flag, u1 );
    value( vui->default_display_window_flag, u1 );
    if( vui->default_display_window_flag ) {
        value( vui->def_disp_win_left_offset, ue );
        value( vui->def_disp_win_right_offset, ue );
        value( vui->def_disp_win_top_offset, ue );
        value( vui->def_disp_win_bottom_offset, ue );
    }
    value( vui->vui_timing_info_present_flag, u1 );
    if( vui->vui_timing_info_present_flag ) {
        value( vui->vui_num_units_in_tick, u(32) );
        value( vui->vui_time_scale, u(32) );
        value( vui->vui_poc_proportional_to_timing_flag, u1 );
        if( vui->vui_poc_proportional_to_timing_flag ) {
            value( vui->vui_num_ticks_poc_diff_one_minus1, ue );
        }
        value( vui->vui_hrd_parameters_present_flag, u1 );
        if( vui->vui_hrd_parameters_present_flag ) {
            structure(hevc_hrd_parameters)( &vui->hrd, b, 1, sps->sps_max_sub_layers_minus1 );
        }
    }
    value( vui->bitstream_restriction_flag, u1 );
    if( vui->bitstream_restriction_flag )
    {
        value( vui->tiles_fixed_structure_flag, u1 );
        value( vui->motion_vectors_over_pic_boundaries_flag, u1 );
        value( vui->restricted_ref_pic_lists_flag, u1 );
        value( vui->min_spatial_segmentation_idc, ue );
        value( vui->max_bytes_per_pic_denom, ue );
        value( vui->max_bits_per_min_cu_denom, ue );
        value( vui->log2_max_mv_length_horizontal, ue );
        value( vui->log2_max_mv_length_vertical, ue );
    }
}

//Appendix E.2.2 HRD parameters syntax
void structure(hevc_hrd_parameters)(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1)
{
    if( commonInfPresentFlag ) {
        value( hrd->nal_hrd_parameters_present_flag, u1 );
        value( hrd->vcl_hrd_parameters_present_flag, u1 );
        if( hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag ){
            value( hrd->sub_pic_hrd_params_present_flag, u1 );
            if( hrd->sub_pic_hrd_params_present_flag ) {
                value( hrd->tick_divisor_minus2, u8 );
                value( hrd->du_cpb_removal_delay_increment_length_minus1, u(5) );
                value( hrd->sub_pic_cpb_params_in_pic_timing_sei_flag, u1 );
                value( hrd->dpb_output_delay_du_length_minus1, u(5) );
            }
            value( hrd->bit_rate_scale, u(4) );
            value( hrd->cpb_size_scale, u(4) );
            if( hrd->sub_pic_hrd_params_present_flag ) {
                value( hrd->cpb_size_du_scale, u(4) );
            }
            value( hrd->initial_cpb_removal_delay_length_minus1, u(5) );
            value( hrd->au_cpb_removal_delay_length_minus1, u(5) );
            value( hrd->dpb_output_delay_length_minus1, u(5) );
        }
    }
    
    for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
        value( hrd->fixed_pic_rate_general_flag[ i ], u1 );
        if( !hrd->fixed_pic_rate_general_flag[ i ] ) {
            value( hrd->fixed_pic_rate_within_cvs_flag[ i ], u1 );
        }
        if( hrd->fixed_pic_rate_within_cvs_flag[ i ] ) {
            value( hrd->elemental_duration_in_tc_minus1[ i ], ue );
        } else {
            value( hrd->low_delay_hrd_flag[ i ], u1 );
        }
        if( hrd->low_delay_hrd_flag[ i ] ) {
            value( hrd->cpb_cnt_minus1[ i ], ue );
        }
        if( hrd->nal_hrd_parameters_present_flag ) {
            structure(hevc_sub_layer_hrd_parameters)(&hrd->sub_layer_hrd_nal[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
        if( hrd->vcl_hrd_parameters_present_flag ) {
            structure(hevc_sub_layer_hrd_parameters)(&hrd->sub_layer_hrd_vcl[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
    }
}

//Appendix E.2.3 Sub-layer HRD parameters syntax
void structure(hevc_sub_layer_hrd_parameters)(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag)
{
    for( int i = 0; i <= CpbCnt; i++ ) {
        value( sub_layer_hrd->bit_rate_value_minus1[i], ue );
        value( sub_layer_hrd->cpb_size_value_minus1[i], ue );
        if( sub_pic_hrd_params_present_flag ) {
            value( sub_layer_hrd->cpb_size_du_value_minus1[i], ue );
            value( sub_layer_hrd->bit_rate_du_value_minus1[i], ue );
        }
        value( sub_layer_hrd->cbr_flag[i], u1 );
    }
}
