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



void read_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b);
void read_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b);
void read_sei_rbsp(hevc_stream_t* h, bs_t* b);
void read_sei_message(hevc_stream_t* h, bs_t* b);
void read_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b);
void read_hevc_end_of_seq_rbsp();
void read_end_of_bitstream_rbsp();
void read_filler_data_rbsp(bs_t* b);
void read_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b);
void read_hevc_rbsp_slice_trailing_bits(bs_t* b);
void read_hevc_rbsp_trailing_bits(bs_t* b);
void read_hevc_byte_alignment(bs_t* b);
void read_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1);
void read_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b );
void read_hevc_slice_header(hevc_stream_t* h, bs_t* b);
void read_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b);
void read_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b);
void read_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets );
void read_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b);
void read_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1);
void read_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag);



//7.3.1 NAL unit syntax
int read_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    int nal_size = size;
    int rbsp_size = size;
    uint8_t* rbsp_buf = (uint8_t*)calloc(1, rbsp_size);

    if( 1 )
    {
        int rc = nal_to_rbsp(buf, &nal_size, rbsp_buf, &rbsp_size);

        if (rc < 0) { free(rbsp_buf); return -1; } // handle conversion error
    }

    if( 0 )
    {
        rbsp_size = size*3/4; // NOTE this may have to be slightly smaller (3/4 smaller, worst case) in order to be guaranteed to fit
    }

    bs_t* b = bs_new(rbsp_buf, rbsp_size);
    /* forbidden_zero_bit */ bs_skip_u(b, 1);
    nal->nal_unit_type = bs_read_u(b, 6);
    nal->nal_layer_id = bs_read_u(b, 6);
    nal->nal_temporal_id_plus1 = bs_read_u(b, 3);

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
            
            read_hevc_slice_layer_rbsp(h, b);
            break;

#ifdef HAVE_SEI
        case NAL_UNIT_TYPE_SEI:
            read_hevc_sei_rbsp(h, b);
            break;
#endif

        case HEVC_NAL_UNIT_TYPE_VPS_NUT: 
            read_hevc_video_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_SPS_NUT: 
            read_hevc_seq_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_PPS_NUT:   
            read_hevc_pic_parameter_set_rbsp(h, b);
            break;

        default:
            return -1;
    }

    if (bs_overrun(b)) { bs_free(b); free(rbsp_buf); return -1; }

    if( 0 )
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
void read_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i, j;

    hevc_vps_t* vps = h->vps;
    if( 1 )
    {
        memset(vps, 0, sizeof(hevc_vps_t));
    }
 
    vps->vps_video_parameter_set_id = bs_read_u(b, 4);
    vps->vps_base_layer_internal_flag = bs_read_u1(b);
    vps->vps_base_layer_available_flag = bs_read_u1(b);
    vps->vps_max_layers_minus1 = bs_read_u(b, 6);
    vps->vps_max_sub_layers_minus1 = bs_read_u(b, 3);
    vps->vps_temporal_id_nesting_flag = bs_read_u1(b);
    /* vps_reserved_0xffff_16bits */ bs_skip_u(b, 16);
    
    read_hevc_profile_tier_level(&vps->ptl, b, 1, vps->vps_max_sub_layers_minus1); 
    
    vps->vps_sub_layer_ordering_info_present_flag = bs_read_u1(b);
    for( i = ( vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1 ); 
            i <= vps->vps_max_sub_layers_minus1; i++ ) {
        vps->vps_max_dec_pic_buffering_minus1[ i ] = bs_read_ue(b);        
        vps->vps_max_num_reorder_pics[ i ] = bs_read_ue(b);        
        vps->vps_max_latency_increase_plus1[ i ] = bs_read_ue(b);        
    }
    vps->vps_max_layer_id = bs_read_u(b, 6);
    vps->vps_num_layer_sets_minus1 = bs_read_ue(b);
    for( i = 1; i <= vps->vps_num_layer_sets_minus1; i++ )
        for( j = 0; j <= vps->vps_max_layer_id; j++ ) {
            vps->layer_id_included_flag[ i ][ j ] = bs_read_u1(b);
        }
    vps->vps_timing_info_present_flag = bs_read_u1(b);
    if( vps->vps_timing_info_present_flag ) {
        vps->vps_num_units_in_tick = bs_read_u(b, 32);
        vps->vps_time_scale = bs_read_u(b, 32);
        vps->vps_poc_proportional_to_timing_flag = bs_read_u1(b);
        if( vps->vps_poc_proportional_to_timing_flag ) {
            vps->vps_num_ticks_poc_diff_one_minus1 = bs_read_ue(b);
        }
        vps->vps_num_hrd_parameters = bs_read_ue(b);
        for( i = 0; i < vps->vps_num_hrd_parameters; i++ ) {
            vps->hrd_layer_set_idx[ i ] = bs_read_ue(b);
            if (i > 0) {
                vps->cprms_present_flag[ i ] = bs_read_u1(b);
            }
            read_hevc_hrd_parameters(&vps->hrd[i], b,
                                           vps->cprms_present_flag[ i ],
                                           vps->vps_max_sub_layers_minus1);
        }
    }
    vps->vps_extension_flag = bs_read_u1(b);
    //TODO: support extension data
    //if (vps->vps_extension_flag)    

    read_hevc_rbsp_trailing_bits(b);
}

//7.3.2.2 Sequence parameter set RBSP syntax
void read_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;

    hevc_sps_t* sps = h->sps;
    if( 1 )
    {
        memset(sps, 0, sizeof(hevc_sps_t));
    }
 
    sps->sps_video_parameter_set_id = bs_read_u(b, 4);
    sps->sps_max_sub_layers_minus1 = bs_read_u(b, 3);
    sps->sps_temporal_id_nesting_flag = bs_read_u1(b);
    read_hevc_profile_tier_level(&sps->ptl, b, 1, sps->sps_max_sub_layers_minus1); 
    sps->sps_seq_parameter_set_id = bs_read_ue(b);
    sps->chroma_format_idc = bs_read_ue(b);
    if( sps->chroma_format_idc == 3 ) {
        sps->separate_colour_plane_flag = bs_read_u1(b);
    }
    sps->pic_width_in_luma_samples = bs_read_ue(b);
    sps->pic_height_in_luma_samples = bs_read_ue(b);
    sps->conformance_window_flag = bs_read_u1(b);
    if( sps->conformance_window_flag ) {
        sps->conf_win_left_offset = bs_read_ue(b);
        sps->conf_win_right_offset = bs_read_ue(b);
        sps->conf_win_top_offset = bs_read_ue(b);
        sps->conf_win_bottom_offset = bs_read_ue(b);
    }
    sps->bit_depth_luma_minus8 = bs_read_ue(b);
    sps->bit_depth_chroma_minus8 = bs_read_ue(b);
    sps->log2_max_pic_order_cnt_lsb_minus4 = bs_read_ue(b);
    sps->sps_sub_layer_ordering_info_present_flag = bs_read_u1(b);
    for( i = ( sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1 ); 
            i <= sps->sps_max_sub_layers_minus1; i++ ) {
        sps->sps_max_dec_pic_buffering_minus1 [ i ] = bs_read_ue(b);
        sps->sps_max_num_reorder_pics [ i ] = bs_read_ue(b);
        sps->sps_max_latency_increase_plus1 [ i ] = bs_read_ue(b);
    }
    sps->log2_min_luma_coding_block_size_minus3 = bs_read_ue(b);
    sps->log2_diff_max_min_luma_coding_block_size = bs_read_ue(b);
    sps->log2_min_luma_transform_block_size_minus2 = bs_read_ue(b);
    sps->log2_diff_max_min_luma_transform_block_size = bs_read_ue(b);
    sps->max_transform_hierarchy_depth_inter = bs_read_ue(b);
    sps->max_transform_hierarchy_depth_intra = bs_read_ue(b);
    sps->scaling_list_enabled_flag = bs_read_u1(b);
    
    if( sps->scaling_list_enabled_flag ) {
        sps->sps_scaling_list_data_present_flag = bs_read_u1(b);
        if( sps->sps_scaling_list_data_present_flag ) {
            read_hevc_scaling_list_data(&sps->scaling_list_data, b); 
        }
    }
    
    sps->amp_enabled_flag = bs_read_u1(b);
    sps->sample_adaptive_offset_enabled_flag = bs_read_u1(b);
    sps->pcm_enabled_flag = bs_read_u1(b);
    if( sps->pcm_enabled_flag ) {
        sps->pcm_sample_bit_depth_luma_minus1 = bs_read_u(b, 4);
        sps->pcm_sample_bit_depth_chroma_minus1 = bs_read_u(b, 4);
        sps->log2_min_pcm_luma_coding_block_size_minus3 = bs_read_ue(b);
        sps->log2_diff_max_min_pcm_luma_coding_block_size = bs_read_ue(b);
        sps->pcm_loop_filter_disabled_flag = bs_read_u1(b);
    }
    sps->num_short_term_ref_pic_sets = bs_read_ue(b);
    for( i = 0; i < sps->num_short_term_ref_pic_sets; i++) {
        read_hevc_st_ref_pic_set(&sps->st_ref_pic_set[i], b, i, sps->num_short_term_ref_pic_sets);
    }
    
    sps->long_term_ref_pics_present_flag = bs_read_u1(b);
    if( sps->long_term_ref_pics_present_flag ) {
        sps->num_long_term_ref_pics_sps = bs_read_ue(b);
        for( i = 0; i < sps->num_long_term_ref_pics_sps; i++ ) {
            sps->lt_ref_pic_poc_lsb_sps[ i ] = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 );
            sps->used_by_curr_pic_lt_sps_flag[ i ] = bs_read_u1(b);
        }
    }
    sps->sps_temporal_mvp_enabled_flag = bs_read_u1(b);
    sps->strong_intra_smoothing_enabled_flag = bs_read_u1(b);
    sps->vui_parameters_present_flag = bs_read_u1(b);
    if( sps->vui_parameters_present_flag ) {
        read_hevc_vui_parameters(sps, b);
    }
    sps->sps_extension_present_flag = bs_read_u1(b);
    
    if( sps->sps_extension_present_flag ) {
        sps->sps_range_extension_flag = bs_read_u1(b);
        sps->sps_multilayer_extension_flag = bs_read_u1(b);
        sps->sps_3d_extension_flag = bs_read_u1(b);
        sps->sps_extension_5bits = bs_read_u(b, 5);
    }
    if( sps->sps_range_extension_flag ) {
        read_hevc_sps_range_extension( &sps->sps_range_ext, b);
    }
    
    if( 1 )
    {
        memcpy(h->sps_table[sps->sps_seq_parameter_set_id], h->sps, sizeof(hevc_sps_t));
    }
}

//7.3.2.2.2 Sequence parameter set range extension syntax
void read_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b)
{
    sps_range_ext->transform_skip_rotation_enabled_flag = bs_read_u1(b);
    sps_range_ext->transform_skip_context_enabled_flag = bs_read_u1(b);
    sps_range_ext->implicit_rdpcm_enabled_flag = bs_read_u1(b);
    sps_range_ext->explicit_rdpcm_enabled_flag = bs_read_u1(b);
    sps_range_ext->extended_precision_processing_flag = bs_read_u1(b);
    sps_range_ext->intra_smoothing_disabled_flag = bs_read_u1(b);
    sps_range_ext->high_precision_offsets_enabled_flag = bs_read_u1(b);
    sps_range_ext->persistent_rice_adaptation_enabled_flag = bs_read_u1(b);
    sps_range_ext->cabac_bypass_alignment_enabled_flag = bs_read_u1(b);
}


//7.3.2.3 Picture parameter set RBSP syntax
void read_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_pps_t* pps = h->pps;
    if( 1 )
    {
        memset(pps, 0, sizeof(hevc_pps_t));
    }

    pps->pic_parameter_set_id = bs_read_ue(b);
    pps->seq_parameter_set_id = bs_read_ue(b);
    pps->dependent_slice_segments_enabled_flag = bs_read_u1(b);
    pps->output_flag_present_flag = bs_read_u1(b);
    pps->num_extra_slice_header_bits = bs_read_u(b,  3 );
    pps->sign_data_hiding_enabled_flag = bs_read_u1(b);
    pps->cabac_init_present_flag = bs_read_u1(b);
    pps->num_ref_idx_l0_default_active_minus1 = bs_read_ue(b);
    pps->num_ref_idx_l1_default_active_minus1 = bs_read_ue(b);
    pps->init_qp_minus26 = bs_read_se(b);
    pps->constrained_intra_pred_flag = bs_read_u1(b);
    pps->transform_skip_enabled_flag = bs_read_u1(b);
    pps->cu_qp_delta_enabled_flag = bs_read_u1(b);
    if( pps->cu_qp_delta_enabled_flag ) {
        pps->diff_cu_qp_delta_depth = bs_read_ue(b);
    }
    pps->pps_cb_qp_offset = bs_read_se(b);
    pps->pps_cr_qp_offset = bs_read_se(b);
    pps->pps_slice_chroma_qp_offsets_present_flag = bs_read_u1(b);
    pps->weighted_pred_flag = bs_read_u1(b);
    pps->weighted_bipred_flag = bs_read_u1(b);
    pps->transquant_bypass_enabled_flag = bs_read_u1(b);
    pps->tiles_enabled_flag = bs_read_u1(b);
    pps->entropy_coding_sync_enabled_flag = bs_read_u1(b);
    if( pps->tiles_enabled_flag ) {
        pps->num_tile_columns_minus1 = bs_read_ue(b);
        pps->num_tile_rows_minus1 = bs_read_ue(b);
        pps->uniform_spacing_flag = bs_read_u1(b);
        if( !pps->uniform_spacing_flag ) {
            for( i = 0; i < pps->num_tile_columns_minus1; i++ ) {
                pps->column_width_minus1[ i ] = bs_read_ue(b);
            }
            for( i = 0; i < pps->num_tile_rows_minus1; i++ ) {
                pps->row_height_minus1[ i ] = bs_read_ue(b);
            }
        }
        pps->loop_filter_across_tiles_enabled_flag = bs_read_u1(b);
    }
    pps->pps_loop_filter_across_slices_enabled_flag = bs_read_u1(b);
    pps->deblocking_filter_control_present_flag = bs_read_u1(b);
    if( pps->deblocking_filter_control_present_flag ) {
        pps->deblocking_filter_override_enabled_flag = bs_read_u1(b);
        pps->pps_deblocking_filter_disabled_flag = bs_read_u1(b);
        if( pps->pps_deblocking_filter_disabled_flag ) {
            pps->pps_beta_offset_div2 = bs_read_se(b);
            pps->pps_tc_offset_div2 = bs_read_se(b);
        }
    }
    pps->pps_scaling_list_data_present_flag = bs_read_u1(b);
    if( pps->pps_scaling_list_data_present_flag ) {
        read_hevc_scaling_list_data(&pps->scaling_list_data, b);
    }
    pps->lists_modification_present_flag = bs_read_u1(b);
    pps->log2_parallel_merge_level_minus2 = bs_read_ue(b);
    pps->slice_segment_header_extension_present_flag = bs_read_u1(b);
    pps->pps_extension_present_flag = bs_read_u1(b);
    if( pps->pps_extension_present_flag ) {
        pps->pps_range_extension_flag = bs_read_u1(b);
        pps->pps_multilayer_extension_flag = bs_read_u1(b);
        pps->pps_3d_extension_flag = bs_read_u1(b);
        pps->pps_extension_5bits = bs_read_u1(b);
    }
    if( pps->pps_range_extension_flag ) {
        read_hevc_pps_range_extension( pps, b);
    }

    read_hevc_rbsp_trailing_bits(b);

    if( 1 )
    {
        memcpy(h->pps_table[pps->pic_parameter_set_id], h->pps, sizeof(hevc_pps_t));
    }
}

//7.3.2.3.2 Picture parameter set range extension syntax
void read_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b)
{
    hevc_pps_range_ext_t *pps_range_ext = &pps->pps_range_ext;;
    if( pps->transform_skip_enabled_flag ) {
        pps_range_ext->log2_max_transform_skip_block_size_minus2 = bs_read_ue(b);
    }
    pps_range_ext->cross_component_prediction_enabled_flag = bs_read_u1(b);
    pps_range_ext->chroma_qp_offset_list_enabled_flag = bs_read_u1(b);
    if( pps_range_ext->chroma_qp_offset_list_enabled_flag ) {
        pps_range_ext->diff_cu_chroma_qp_offset_depth = bs_read_ue(b);
        pps_range_ext->chroma_qp_offset_list_len_minus1 = bs_read_ue(b);
        for( int i = 0; i <= pps_range_ext->chroma_qp_offset_list_len_minus1; i++ ) {
            pps_range_ext->cb_qp_offset_list[ i ] = bs_read_se(b);
            pps_range_ext->cr_qp_offset_list[ i ] = bs_read_se(b);
        }
    }
    pps_range_ext->log2_sao_offset_scale_luma = bs_read_ue(b);
    pps_range_ext->log2_sao_offset_scale_chroma = bs_read_ue(b);
}

#ifdef HAVE_SEI
//7.3.2.4 Supplemental enhancement information RBSP syntax
void read_sei_rbsp(hevc_stream_t* h, bs_t* b)
{
    if( 1 )
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
            read_sei_message(h, b);
        } while( more_rbsp_data(h, b) );
    }

    if( 0 )
    {
        for (int i = 0; i < h->num_seis; i++) {
            h->sei = h->seis[i];
            read_sei_message(h, b);
        }
        h->sei = NULL;
    }

    read_hevc_rbsp_trailing_bits(b);
}

//7.3.5 Supplemental enhancement information message syntax
void read_sei_message(hevc_stream_t* h, bs_t* b)
{
    if( 0 )
    {
        _write_ff_coded_number(b, h->sei->payloadType);
        _write_ff_coded_number(b, h->sei->payloadSize);
    }
    if( 1 )
    {
        h->sei->payloadType = _read_ff_coded_number(b);
        h->sei->payloadSize = _read_ff_coded_number(b);
    }
    read_sei_payload( h, b, h->sei->payloadType, h->sei->payloadSize );
}
#endif

//7.3.2.5 Access unit delimiter RBSP syntax
void read_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b)
{
    h->aud->primary_pic_type = bs_read_u(b, 3);
    read_hevc_rbsp_trailing_bits(b);
}

//7.3.2.6 End of sequence RBSP syntax
void read_hevc_end_of_seq_rbsp()
{
}

//7.3.2.7 End of bitstream RBSP syntax
void read_end_of_bitstream_rbsp()
{
}

//7.3.2.8 Filler data RBSP syntax
void read_filler_data_rbsp(bs_t* b)
{
    while( bs_next_bits(b, 8) == 0xFF )
    {
        /* ff_byte */ bs_skip_u(b, 8);
    }
    read_hevc_rbsp_trailing_bits(b);
}

//7.3.2.9 Slice segment layer RBSP syntax
void read_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b)
{
    read_hevc_slice_header(h, b);
    hevc_slice_data_rbsp_t* slice_data = h->slice_data;

    if ( slice_data != NULL )
    {
        if ( slice_data->rbsp_buf != NULL ) free( slice_data->rbsp_buf ); 
        uint8_t *sptr = b->p + (!!b->bits_left); // CABAC-specific: skip alignment bits, if there are any
        slice_data->rbsp_size = b->end - sptr;
        
        slice_data->rbsp_buf = (uint8_t*)malloc(slice_data->rbsp_size);
        memcpy( slice_data->rbsp_buf, sptr, slice_data->rbsp_size );
    }

    //read_hevc_slice_data(h, b); /* all categories of slice_data( ) syntax */
    read_hevc_rbsp_slice_trailing_bits( b );
}

//7.3.2.10 RBSP slice trailing bits syntax
void read_hevc_rbsp_slice_trailing_bits(bs_t* b)
{
    read_hevc_rbsp_trailing_bits(b);
    //while( more_rbsp_trailing_data(b) )
    //{
    //    value( cabac_zero_word, f(16, 0x0000) );
    //}
}

//7.3.2.11 RBSP trailing bits syntax
void read_hevc_rbsp_trailing_bits(bs_t* b)
{
    /* rbsp_stop_one_bit */ bs_skip_u(b, 1);

    while( !bs_byte_aligned(b) )
    {
        /* rbsp_alignment_zero_bit */ bs_skip_u(b, 1);
    }
}

//7.3.2.12 Byte alignment syntax
void read_hevc_byte_alignment(bs_t* b)
{
    /* alignment_bit_equal_to_one */ bs_skip_u(b, 1);

    while( !bs_byte_aligned(b) )
    {
        /* alignment_bit_equal_to_zero */ bs_skip_u(b, 1);
    }
}

//7.3.3 Profile, tier and level syntax
void read_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1)
{
    int i, j;
    if( profilePresentFlag ) {
        ptl->general_profile_space = bs_read_u(b, 2);
        ptl->general_tier_flag = bs_read_u1(b);
        ptl->general_profile_idc = bs_read_u(b, 5);
        for( i = 0; i < 32; i++ ) {
            ptl->general_profile_compatibility_flag[ i ] = bs_read_u1(b);
        }
        ptl->general_progressive_source_flag = bs_read_u1(b);
        ptl->general_interlaced_source_flag = bs_read_u1(b);
        ptl->general_non_packed_constraint_flag = bs_read_u1(b);
        ptl->general_frame_only_constraint_flag = bs_read_u1(b);
        if( ptl->general_profile_idc == 4 || ptl->general_profile_compatibility_flag[ 4 ] || 
            ptl->general_profile_idc == 5 || ptl->general_profile_compatibility_flag[ 5 ] || 
            ptl->general_profile_idc == 6 || ptl->general_profile_compatibility_flag[ 6 ] || 
            ptl->general_profile_idc == 7 || ptl->general_profile_compatibility_flag[ 7 ] ) {
                
            ptl->general_max_12bit_constraint_flag = bs_read_u1(b);
            ptl->general_max_10bit_constraint_flag = bs_read_u1(b);
            ptl->general_max_8bit_constraint_flag = bs_read_u1(b);
            ptl->general_max_422chroma_constraint_flag = bs_read_u1(b);
            ptl->general_max_420chroma_constraint_flag = bs_read_u1(b);
            ptl->general_max_monochrome_constraint_flag = bs_read_u1(b);
            ptl->general_intra_constraint_flag = bs_read_u1(b);
            ptl->general_one_picture_only_constraint_flag = bs_read_u1(b);
            ptl->general_lower_bit_rate_constraint_flag = bs_read_u1(b);
            /* general_reserved_zero_34bits */ bs_skip_u(b, 34);
        } else {
            /* general_reserved_zero_43bits */ bs_skip_u(b, 43);
        }
        if( ( ptl->general_profile_idc >= 1 && ptl->general_profile_idc <= 5 ) ||
              ptl->general_profile_compatibility_flag[ 1 ] ||
              ptl->general_profile_compatibility_flag[ 2 ] ||
              ptl->general_profile_compatibility_flag[ 3 ] ||
              ptl->general_profile_compatibility_flag[ 4 ] ||
              ptl->general_profile_compatibility_flag[ 5 ] ) {

            ptl->general_inbld_flag = bs_read_u1(b);
        } else {
            /* general_reserved_zero_bit */ bs_skip_u(b, 1);
        }
        ptl->general_level_idc = bs_read_u8(b);
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) {
            ptl->sub_layer_profile_present_flag[ i ] = bs_read_u1(b);
            ptl->sub_layer_level_present_flag[ i ] = bs_read_u1(b);
        }
        if( maxNumSubLayersMinus1 > 0 ) {
            for( i = maxNumSubLayersMinus1; i < 8; i++ ) {
                /* reserved_zero_xxbits */ bs_skip_u(b, 2);
            }
        }
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) { 
            if( ptl->sub_layer_profile_present_flag[ i ] ) {
                ptl->sub_layer_profile_space[ i ] = bs_read_u(b, 2);
                ptl->sub_layer_tier_flag[ i ] = bs_read_u1(b);
                ptl->sub_layer_profile_idc[ i ] = bs_read_u(b, 5);
                for( j = 0; j < 32; j++ ) {
                    ptl->sub_layer_profile_compatibility_flag[ i ][ j ] = bs_read_u(b, 1);
                }
                ptl->sub_layer_progressive_source_flag[ i ] = bs_read_u1(b);
                ptl->sub_layer_interlaced_source_flag[ i ] = bs_read_u1(b);
                ptl->sub_layer_non_packed_constraint_flag[ i ] = bs_read_u1(b);
                ptl->sub_layer_frame_only_constraint_flag[ i ] = bs_read_u1(b);
                if( ptl->sub_layer_profile_idc[ i ] == 4 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 4 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 5 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 5 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 6 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 6 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 7 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 7 ] ) {
                    ptl->sub_layer_max_12bit_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_max_10bit_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_max_8bit_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_max_422chroma_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_max_420chroma_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_max_monochrome_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_intra_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_one_picture_only_constraint_flag[ i ] = bs_read_u1(b);
                    ptl->sub_layer_lower_bit_rate_constraint_flag[ i ] = bs_read_u1(b);
                    /* sub_layer_reserved_zero_34bits */ bs_skip_u(b, 34);
                } else {
                    /* sub_layer_reserved_zero_43bits */ bs_skip_u(b, 43);
                }
            
                if( ( ptl->sub_layer_profile_idc[ i ] >= 1 && ptl->sub_layer_profile_idc[ i ] <= 5 ) ||
                   ptl->sub_layer_profile_compatibility_flag[ 1 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 2 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 3 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 4 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 5 ] ) {
                    ptl->sub_layer_inbld_flag[ i ] = bs_read_u1(b);
                } else {
                    /* sub_layer_reserved_zero_bit */ bs_skip_u(b, 1);
                }
            }
            if( ptl->sub_layer_level_present_flag[ i ] ) {
                ptl->sub_layer_level_idc[ i ] = bs_read_u1(b);
            }
        }
    }
}

//7.3.4 Scaling list data syntax
void read_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b )
{
    int nextCoef, coefNum;
    for( int sizeId = 0; sizeId < 4; sizeId++ )
        for( int matrixId = 0; matrixId < 6; matrixId += ( sizeId == 3 ) ? 3 : 1 ) {
            sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] = bs_read_u1(b);
            if( !sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] ) {
                sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ] = bs_read_ue(b);
            } else {
                nextCoef = 8;
                coefNum=MIN(64, (1 << (4+(sizeId << 1))));
                if( sizeId > 1 ) {
                    sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ] = bs_read_se(b);
                }
 
                for( int i = 0; i < coefNum; i++) {
                    sld->scaling_list_delta_coef[ sizeId ][ matrixId ] = bs_read_se(b);
                }
            }
 
        }
}

//7.3.6 Slice header syntax
void read_hevc_slice_header(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    if( 1 )
    {
        init_slice_hevc(h);
    }

    hevc_nal_t* nal = h->nal;

    sh->first_slice_segment_in_pic_flag = bs_read_u1(b);
    if( nal->nal_unit_type >= HEVC_NAL_UNIT_TYPE_BLA_W_LP && 
       nal->nal_unit_type <= HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23) {
            sh->no_output_of_prior_pics_flag = bs_read_u1(b);
    }
    sh->pic_parameter_set_id = bs_read_ue(b);
    
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];

    //set default value
    sh->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    sh->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;

    if( !sh->first_slice_segment_in_pic_flag ) {
        if( pps->dependent_slice_segments_enabled_flag ) {
            sh->dependent_slice_segment_flag = bs_read_u1(b);
        }
        sh->slice_segment_address = bs_read_u(b,  getSliceSegmentAddressBitLength( sps ) );
    }
    
    if( !sh->dependent_slice_segment_flag ) {
        for( i = 0; i < pps->num_extra_slice_header_bits; i++ ) {
            /* slice_reserved_flag */ bs_skip_u(b, 1);
        }
        sh->slice_type = bs_read_ue(b);
        if( pps->output_flag_present_flag ) {
            sh->pic_output_flag = bs_read_u1(b);
        }
        if( sps->separate_colour_plane_flag == 1 ) {
            sh->colour_plane_id = bs_read_u(b, 2);
        }
        if( nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_W_RADL &&
            nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_N_LP) {
            sh->slice_pic_order_cnt_lsb = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 );
            sh->short_term_ref_pic_set_sps_flag = bs_read_u1(b);
            if( !sh->short_term_ref_pic_set_sps_flag ) {
                read_hevc_st_ref_pic_set( &sh->st_ref_pic_set, b, sps->num_short_term_ref_pic_sets, sps->num_short_term_ref_pic_sets );
            } else if( sps->num_short_term_ref_pic_sets > 1 ) {
                sh->short_term_ref_pic_set_idx = bs_read_u(b,  ceil( log2( sps->num_short_term_ref_pic_sets ) ) );
            }
            if( sps->long_term_ref_pics_present_flag ) {
                if( sps->num_long_term_ref_pics_sps > 0 ) {
                    sh->num_long_term_sps = bs_read_ue(b);
                }
                sh->num_long_term_pics = bs_read_ue(b);
                for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
                    if( i < sh->num_long_term_sps ) {
                        if( sps->num_long_term_ref_pics_sps > 1 ) {
                            sh->lt_idx_sps[ i ] = bs_read_u(b,  ceil( log2( sps->num_long_term_ref_pics_sps ) ) );
                        }
                    } else {
                        sh->poc_lsb_lt[ i ] = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 );
                        sh->used_by_curr_pic_lt_flag[ i ] = bs_read_u1(b);
                    }
                    sh->delta_poc_msb_present_flag[ i ] = bs_read_u1(b);
                    if( sh->delta_poc_msb_present_flag[ i ]) {
                        sh->delta_poc_msb_cycle_lt[ i ] = bs_read_ue(b);
                    }
                }
            }
            if( sps->sps_temporal_mvp_enabled_flag ) {
                sh->slice_temporal_mvp_enabled_flag = bs_read_u1(b);
            }
        }
        if( sps->sample_adaptive_offset_enabled_flag ) {
            sh->slice_sao_luma_flag = bs_read_u1(b);
            int ChromaArrayType = 0;
            if( sps->separate_colour_plane_flag == 0) {
                ChromaArrayType = sps->chroma_format_idc;
            }
            if( ChromaArrayType != 0 ) {
                sh->slice_sao_chroma_flag = bs_read_u1(b);
            }
        }
        if( sh->slice_type == HEVC_SLICE_TYPE_P || sh->slice_type == HEVC_SLICE_TYPE_B ){
            sh->num_ref_idx_active_override_flag = bs_read_u1(b);
            if( sh->num_ref_idx_active_override_flag ) {
                sh->num_ref_idx_l0_active_minus1 = bs_read_ue(b);
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    sh->num_ref_idx_l1_active_minus1 = bs_read_ue(b);
                }
            }
            if( pps->lists_modification_present_flag && getNumPicTotalCurr( sps, sh ) > 1 ) {
                read_hevc_ref_pic_lists_modification( h, b );
            }
            
            if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                sh->mvd_l1_zero_flag = bs_read_u1(b);
            }
            if( pps->cabac_init_present_flag ) {
                sh->cabac_init_flag = bs_read_u1(b);
            }
            if( sh->slice_temporal_mvp_enabled_flag ) {
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    sh->collocated_from_l0_flag = bs_read_u1(b);
                }
                if( ( sh->collocated_from_l0_flag && sh->num_ref_idx_l0_active_minus1 > 0 ) ||
                   ( !sh->collocated_from_l0_flag && sh->num_ref_idx_l1_active_minus1 > 0 ) ) {
                    sh->collocated_ref_idx = bs_read_ue(b);
                }
            }
            if( ( pps->weighted_pred_flag && sh->slice_type == HEVC_SLICE_TYPE_P ) || 
                ( pps->weighted_bipred_flag && sh->slice_type == HEVC_SLICE_TYPE_B ) ) {
                read_hevc_pred_weight_table( h, b );
            }
            sh->five_minus_max_num_merge_cand = bs_read_ue(b);
        }
        sh->slice_qp_delta = bs_read_se(b);
        if( pps->pps_slice_chroma_qp_offsets_present_flag ) {
            sh->slice_cb_qp_offset = bs_read_se(b);
            sh->slice_cr_qp_offset = bs_read_se(b);
        }
        if( pps->pps_range_ext.chroma_qp_offset_list_enabled_flag ) {
            sh->cu_chroma_qp_offset_enabled_flag = bs_read_u1(b);
        }
        if( pps->deblocking_filter_override_enabled_flag ) {
            sh->deblocking_filter_override_flag = bs_read_u1(b);
        }
        if( sh->deblocking_filter_override_flag ) {
            sh->slice_deblocking_filter_disabled_flag = bs_read_u1(b);
            if( !sh->slice_deblocking_filter_disabled_flag ) {
                sh->slice_beta_offset_div2 = bs_read_se(b);
                sh->slice_tc_offset_div2 = bs_read_se(b);
            }
        }
        if( pps->pps_loop_filter_across_slices_enabled_flag &&
           ( sh->slice_sao_luma_flag || sh->slice_sao_chroma_flag || !sh->slice_deblocking_filter_disabled_flag ) ) {
            sh->slice_loop_filter_across_slices_enabled_flag = bs_read_u1(b);
        }
    }
    if( pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag ) {
        sh->num_entry_point_offsets = bs_read_ue(b);
        if( sh->num_entry_point_offsets > 0 ) {
            sh->offset_len_minus1 = bs_read_ue(b);
            for( i = 0; i < sh->num_entry_point_offsets; i++ ) {
                sh->entry_point_offset_minus1[ i ] = bs_read_u(b,  sh->offset_len_minus1 + 1 );
            }
        }
    }
    if( pps->slice_segment_header_extension_present_flag ) {
        sh->slice_segment_header_extension_length = bs_read_ue(b);
        //TODO: support header extension,
        for( i = 0; i < sh->slice_segment_header_extension_length; i++) {
            /* slice_segment_header_extension_data_byte */ bs_skip_u(b, 8);
        }
    }
    read_hevc_byte_alignment( b );
}

//7.3.6.2 Reference picture list reordering syntax
void read_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    
    sh->rpld.ref_pic_list_modification_flag_l0 = bs_read_u1(b);
    if( sh->rpld.ref_pic_list_modification_flag_l0 ) {
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            sh->rpld.list_entry_l0[ i ] = bs_read_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) );
        }
    }
    
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        // ERROR: value( sh->rpld.ref_pic_list_modification_flag_l1, 1 );
        if( sh->rpld.ref_pic_list_modification_flag_l1 ) {
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                sh->rpld.list_entry_l1[ i ] = bs_read_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) );
            }
        }
    }
}

//7.3.6.3 Prediction weight table syntax
void read_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b)
{
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    hevc_pred_weight_table_t *pwt = &sh->pwt;

    int i, j;

    pwt->luma_log2_weight_denom = bs_read_ue(b);
    
    int ChromaArrayType = 0;
    if( sps->separate_colour_plane_flag == 0) {
        ChromaArrayType = sps->chroma_format_idc;
    }

    if( ChromaArrayType != 0 ) {
        pwt->delta_chroma_log2_weight_denom = bs_read_se(b);
    }
    
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        pwt->luma_weight_l0_flag[i] = bs_read_u1(b);
    }
    if( ChromaArrayType != 0 )
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            pwt->chroma_weight_l0_flag[i] = bs_read_u1(b);
        }
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        if( pwt->luma_weight_l0_flag[i] ) {
            pwt->delta_luma_weight_l0[i] = bs_read_se(b);
            pwt->luma_offset_l0[i] = bs_read_se(b);
        }
        if( pwt->chroma_weight_l0_flag[i] ) {
            for( j =0; j < 2; j++ ) {
                pwt->delta_chroma_weight_l0[i][j] = bs_read_se(b);
                pwt->delta_chroma_offset_l0[i][j] = bs_read_se(b);
            }
        }
    }
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            pwt->luma_weight_l1_flag[i] = bs_read_u1(b);
        }
        if( ChromaArrayType != 0 )
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                pwt->chroma_weight_l1_flag[i] = bs_read_u1(b);
            }
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            if( pwt->luma_weight_l1_flag[i] ) {
                pwt->delta_luma_weight_l1[i] = bs_read_se(b);
                pwt->luma_offset_l1[i] = bs_read_se(b);
            }
            if( pwt->chroma_weight_l1_flag[i] ) {
                for( j =0; j < 2; j++ ) {
                    pwt->delta_chroma_weight_l1[i][j] = bs_read_se(b);
                    pwt->delta_chroma_offset_l1[i][j] = bs_read_se(b);
                }
            }
        }        
    }
}

//7.3.7 Short-term reference picture set syntax
void read_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets )
{
    int i, j;
    
    if( stRpsIdx != 0 ) {
        st_ref_pic_set->inter_ref_pic_set_prediction_flag = bs_read_u1(b);
    }
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        if( stRpsIdx == num_short_term_ref_pic_sets ) {
            st_ref_pic_set->delta_idx_minus1 = bs_read_ue(b);
        }
        st_ref_pic_set->delta_rps_sign = bs_read_u1(b);
        st_ref_pic_set->abs_delta_rps_minus1 = bs_read_ue(b);
        
        int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
        
        for( j = 0; j <= NumDeltaPocs[ RefRpsIdx ]; j++ ) {
            st_ref_pic_set->used_by_curr_pic_flag[ j ] = bs_read_u1(b);
            if( !st_ref_pic_set->used_by_curr_pic_flag[ j ] ) {
                st_ref_pic_set->use_delta_flag[ j ] = bs_read_u1(b);
            }
        }
    } else {
        st_ref_pic_set->num_negative_pics = bs_read_ue(b);
        st_ref_pic_set->num_positive_pics = bs_read_ue(b);
        for( i = 0; i < st_ref_pic_set->num_negative_pics; i++ ) {
            st_ref_pic_set->delta_poc_s0_minus1[ i ] = bs_read_ue(b);
            st_ref_pic_set->used_by_curr_pic_s0_flag[ i ] = bs_read_u1(b);
            
            //update derived field
            UsedByCurrPicS0[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s0_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS0[ stRpsIdx ][ i ] = -1 * ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            } else {
                DeltaPocS0[ stRpsIdx ][ i ] = DeltaPocS0[ stRpsIdx ][ i - 1 ] - ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            }
        }
        for( i = 0; i < st_ref_pic_set->num_positive_pics; i++ ) {
            st_ref_pic_set->delta_poc_s1_minus1[ i ] = bs_read_ue(b);
            st_ref_pic_set->used_by_curr_pic_s1_flag[ i ] = bs_read_u1(b);
        
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
void read_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b)
{
    hevc_vui_t* vui = &sps->vui;
    vui->aspect_ratio_info_present_flag = bs_read_u1(b);
    if( vui->aspect_ratio_info_present_flag )
    {
        vui->aspect_ratio_idc = bs_read_u8(b);
        if( vui->aspect_ratio_idc == SAR_Extended )
        {
            vui->sar_width = bs_read_u(b, 16);
            vui->sar_height = bs_read_u(b, 16);
        }
    }
    vui->overscan_info_present_flag = bs_read_u1(b);
    if( vui->overscan_info_present_flag ) {
        vui->overscan_appropriate_flag = bs_read_u1(b);
    }
    vui->video_signal_type_present_flag = bs_read_u1(b);
    if( vui->video_signal_type_present_flag ) {
        vui->video_format = bs_read_u(b, 3);
        vui->video_full_range_flag = bs_read_u1(b);
        vui->colour_description_present_flag = bs_read_u1(b);
        if( vui->colour_description_present_flag ) {
            vui->colour_primaries = bs_read_u8(b);
            vui->transfer_characteristics = bs_read_u8(b);
            vui->matrix_coefficients = bs_read_u8(b);
        }
    }
    vui->chroma_loc_info_present_flag = bs_read_u1(b);
    if( vui->chroma_loc_info_present_flag ) {
        vui->chroma_sample_loc_type_top_field = bs_read_ue(b);
        vui->chroma_sample_loc_type_bottom_field = bs_read_ue(b);
    }
    
    vui->neutral_chroma_indication_flag = bs_read_u1(b);
    vui->field_seq_flag = bs_read_u1(b);
    vui->frame_field_info_present_flag = bs_read_u1(b);
    vui->default_display_window_flag = bs_read_u1(b);
    if( vui->default_display_window_flag ) {
        vui->def_disp_win_left_offset = bs_read_ue(b);
        vui->def_disp_win_right_offset = bs_read_ue(b);
        vui->def_disp_win_top_offset = bs_read_ue(b);
        vui->def_disp_win_bottom_offset = bs_read_ue(b);
    }
    vui->vui_timing_info_present_flag = bs_read_u1(b);
    if( vui->vui_timing_info_present_flag ) {
        vui->vui_num_units_in_tick = bs_read_u(b, 32);
        vui->vui_time_scale = bs_read_u(b, 32);
        vui->vui_poc_proportional_to_timing_flag = bs_read_u1(b);
        if( vui->vui_poc_proportional_to_timing_flag ) {
            vui->vui_num_ticks_poc_diff_one_minus1 = bs_read_ue(b);
        }
        vui->vui_hrd_parameters_present_flag = bs_read_u1(b);
        if( vui->vui_hrd_parameters_present_flag ) {
            read_hevc_hrd_parameters( &vui->hrd, b, 1, sps->sps_max_sub_layers_minus1 );
        }
    }
    vui->bitstream_restriction_flag = bs_read_u1(b);
    if( vui->bitstream_restriction_flag )
    {
        vui->tiles_fixed_structure_flag = bs_read_u1(b);
        vui->motion_vectors_over_pic_boundaries_flag = bs_read_u1(b);
        vui->restricted_ref_pic_lists_flag = bs_read_u1(b);
        vui->min_spatial_segmentation_idc = bs_read_ue(b);
        vui->max_bytes_per_pic_denom = bs_read_ue(b);
        vui->max_bits_per_min_cu_denom = bs_read_ue(b);
        vui->log2_max_mv_length_horizontal = bs_read_ue(b);
        vui->log2_max_mv_length_vertical = bs_read_ue(b);
    }
}

//Appendix E.2.2 HRD parameters syntax
void read_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1)
{
    if( commonInfPresentFlag ) {
        hrd->nal_hrd_parameters_present_flag = bs_read_u1(b);
        hrd->vcl_hrd_parameters_present_flag = bs_read_u1(b);
        if( hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag ){
            hrd->sub_pic_hrd_params_present_flag = bs_read_u1(b);
            if( hrd->sub_pic_hrd_params_present_flag ) {
                hrd->tick_divisor_minus2 = bs_read_u8(b);
                hrd->du_cpb_removal_delay_increment_length_minus1 = bs_read_u(b, 5);
                hrd->sub_pic_cpb_params_in_pic_timing_sei_flag = bs_read_u1(b);
                hrd->dpb_output_delay_du_length_minus1 = bs_read_u(b, 5);
            }
            hrd->bit_rate_scale = bs_read_u(b, 4);
            hrd->cpb_size_scale = bs_read_u(b, 4);
            if( hrd->sub_pic_hrd_params_present_flag ) {
                hrd->cpb_size_du_scale = bs_read_u(b, 4);
            }
            hrd->initial_cpb_removal_delay_length_minus1 = bs_read_u(b, 5);
            hrd->au_cpb_removal_delay_length_minus1 = bs_read_u(b, 5);
            hrd->dpb_output_delay_length_minus1 = bs_read_u(b, 5);
        }
    }
    
    for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
        hrd->fixed_pic_rate_general_flag[ i ] = bs_read_u1(b);
        if( !hrd->fixed_pic_rate_general_flag[ i ] ) {
            hrd->fixed_pic_rate_within_cvs_flag[ i ] = bs_read_u1(b);
        }
        if( hrd->fixed_pic_rate_within_cvs_flag[ i ] ) {
            hrd->elemental_duration_in_tc_minus1[ i ] = bs_read_ue(b);
        } else {
            hrd->low_delay_hrd_flag[ i ] = bs_read_u1(b);
        }
        if( hrd->low_delay_hrd_flag[ i ] ) {
            hrd->cpb_cnt_minus1[ i ] = bs_read_ue(b);
        }
        if( hrd->nal_hrd_parameters_present_flag ) {
            read_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_nal[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
        if( hrd->vcl_hrd_parameters_present_flag ) {
            read_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_vcl[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
    }
}

//Appendix E.2.3 Sub-layer HRD parameters syntax
void read_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag)
{
    for( int i = 0; i <= CpbCnt; i++ ) {
        sub_layer_hrd->bit_rate_value_minus1[i] = bs_read_ue(b);
        sub_layer_hrd->cpb_size_value_minus1[i] = bs_read_ue(b);
        if( sub_pic_hrd_params_present_flag ) {
            sub_layer_hrd->cpb_size_du_value_minus1[i] = bs_read_ue(b);
            sub_layer_hrd->bit_rate_du_value_minus1[i] = bs_read_ue(b);
        }
        sub_layer_hrd->cbr_flag[i] = bs_read_u1(b);
    }
}


void write_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b);
void write_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b);
void write_sei_rbsp(hevc_stream_t* h, bs_t* b);
void write_sei_message(hevc_stream_t* h, bs_t* b);
void write_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b);
void write_hevc_end_of_seq_rbsp();
void write_end_of_bitstream_rbsp();
void write_filler_data_rbsp(bs_t* b);
void write_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b);
void write_hevc_rbsp_slice_trailing_bits(bs_t* b);
void write_hevc_rbsp_trailing_bits(bs_t* b);
void write_hevc_byte_alignment(bs_t* b);
void write_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1);
void write_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b );
void write_hevc_slice_header(hevc_stream_t* h, bs_t* b);
void write_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b);
void write_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b);
void write_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets );
void write_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b);
void write_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1);
void write_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag);



//7.3.1 NAL unit syntax
int write_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    int nal_size = size;
    int rbsp_size = size;
    uint8_t* rbsp_buf = (uint8_t*)calloc(1, rbsp_size);

    if( 0 )
    {
        int rc = nal_to_rbsp(buf, &nal_size, rbsp_buf, &rbsp_size);

        if (rc < 0) { free(rbsp_buf); return -1; } // handle conversion error
    }

    if( 1 )
    {
        rbsp_size = size*3/4; // NOTE this may have to be slightly smaller (3/4 smaller, worst case) in order to be guaranteed to fit
    }

    bs_t* b = bs_new(rbsp_buf, rbsp_size);
    /* forbidden_zero_bit */ bs_write_u(b, 1, 0);
    bs_write_u(b, 6, nal->nal_unit_type);
    bs_write_u(b, 6, nal->nal_layer_id);
    bs_write_u(b, 3, nal->nal_temporal_id_plus1);

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
            
            write_hevc_slice_layer_rbsp(h, b);
            break;

#ifdef HAVE_SEI
        case NAL_UNIT_TYPE_SEI:
            write_hevc_sei_rbsp(h, b);
            break;
#endif

        case HEVC_NAL_UNIT_TYPE_VPS_NUT: 
            write_hevc_video_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_SPS_NUT: 
            write_hevc_seq_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_PPS_NUT:   
            write_hevc_pic_parameter_set_rbsp(h, b);
            break;

        default:
            return -1;
    }

    if (bs_overrun(b)) { bs_free(b); free(rbsp_buf); return -1; }

    if( 1 )
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
void write_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i, j;

    hevc_vps_t* vps = h->vps;
    if( 0 )
    {
        memset(vps, 0, sizeof(hevc_vps_t));
    }
 
    bs_write_u(b, 4, vps->vps_video_parameter_set_id);
    bs_write_u1(b, vps->vps_base_layer_internal_flag);
    bs_write_u1(b, vps->vps_base_layer_available_flag);
    bs_write_u(b, 6, vps->vps_max_layers_minus1);
    bs_write_u(b, 3, vps->vps_max_sub_layers_minus1);
    bs_write_u1(b, vps->vps_temporal_id_nesting_flag);
    /* vps_reserved_0xffff_16bits */ bs_write_u(b, 16, 0xffff);
    
    write_hevc_profile_tier_level(&vps->ptl, b, 1, vps->vps_max_sub_layers_minus1); 
    
    bs_write_u1(b, vps->vps_sub_layer_ordering_info_present_flag);
    for( i = ( vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1 ); 
            i <= vps->vps_max_sub_layers_minus1; i++ ) {
        bs_write_ue(b, vps->vps_max_dec_pic_buffering_minus1[ i ]);        
        bs_write_ue(b, vps->vps_max_num_reorder_pics[ i ]);        
        bs_write_ue(b, vps->vps_max_latency_increase_plus1[ i ]);        
    }
    bs_write_u(b, 6, vps->vps_max_layer_id);
    bs_write_ue(b, vps->vps_num_layer_sets_minus1);
    for( i = 1; i <= vps->vps_num_layer_sets_minus1; i++ )
        for( j = 0; j <= vps->vps_max_layer_id; j++ ) {
            bs_write_u1(b, vps->layer_id_included_flag[ i ][ j ]);
        }
    bs_write_u1(b, vps->vps_timing_info_present_flag);
    if( vps->vps_timing_info_present_flag ) {
        bs_write_u(b, 32, vps->vps_num_units_in_tick);
        bs_write_u(b, 32, vps->vps_time_scale);
        bs_write_u1(b, vps->vps_poc_proportional_to_timing_flag);
        if( vps->vps_poc_proportional_to_timing_flag ) {
            bs_write_ue(b, vps->vps_num_ticks_poc_diff_one_minus1);
        }
        bs_write_ue(b, vps->vps_num_hrd_parameters);
        for( i = 0; i < vps->vps_num_hrd_parameters; i++ ) {
            bs_write_ue(b, vps->hrd_layer_set_idx[ i ]);
            if (i > 0) {
                bs_write_u1(b, vps->cprms_present_flag[ i ]);
            }
            write_hevc_hrd_parameters(&vps->hrd[i], b,
                                           vps->cprms_present_flag[ i ],
                                           vps->vps_max_sub_layers_minus1);
        }
    }
    bs_write_u1(b, vps->vps_extension_flag);
    //TODO: support extension data
    //if (vps->vps_extension_flag)    

    write_hevc_rbsp_trailing_bits(b);
}

//7.3.2.2 Sequence parameter set RBSP syntax
void write_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;

    hevc_sps_t* sps = h->sps;
    if( 0 )
    {
        memset(sps, 0, sizeof(hevc_sps_t));
    }
 
    bs_write_u(b, 4, sps->sps_video_parameter_set_id);
    bs_write_u(b, 3, sps->sps_max_sub_layers_minus1);
    bs_write_u1(b, sps->sps_temporal_id_nesting_flag);
    write_hevc_profile_tier_level(&sps->ptl, b, 1, sps->sps_max_sub_layers_minus1); 
    bs_write_ue(b, sps->sps_seq_parameter_set_id);
    bs_write_ue(b, sps->chroma_format_idc);
    if( sps->chroma_format_idc == 3 ) {
        bs_write_u1(b, sps->separate_colour_plane_flag);
    }
    bs_write_ue(b, sps->pic_width_in_luma_samples);
    bs_write_ue(b, sps->pic_height_in_luma_samples);
    bs_write_u1(b, sps->conformance_window_flag);
    if( sps->conformance_window_flag ) {
        bs_write_ue(b, sps->conf_win_left_offset);
        bs_write_ue(b, sps->conf_win_right_offset);
        bs_write_ue(b, sps->conf_win_top_offset);
        bs_write_ue(b, sps->conf_win_bottom_offset);
    }
    bs_write_ue(b, sps->bit_depth_luma_minus8);
    bs_write_ue(b, sps->bit_depth_chroma_minus8);
    bs_write_ue(b, sps->log2_max_pic_order_cnt_lsb_minus4);
    bs_write_u1(b, sps->sps_sub_layer_ordering_info_present_flag);
    for( i = ( sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1 ); 
            i <= sps->sps_max_sub_layers_minus1; i++ ) {
        bs_write_ue(b, sps->sps_max_dec_pic_buffering_minus1 [ i ]);
        bs_write_ue(b, sps->sps_max_num_reorder_pics [ i ]);
        bs_write_ue(b, sps->sps_max_latency_increase_plus1 [ i ]);
    }
    bs_write_ue(b, sps->log2_min_luma_coding_block_size_minus3);
    bs_write_ue(b, sps->log2_diff_max_min_luma_coding_block_size);
    bs_write_ue(b, sps->log2_min_luma_transform_block_size_minus2);
    bs_write_ue(b, sps->log2_diff_max_min_luma_transform_block_size);
    bs_write_ue(b, sps->max_transform_hierarchy_depth_inter);
    bs_write_ue(b, sps->max_transform_hierarchy_depth_intra);
    bs_write_u1(b, sps->scaling_list_enabled_flag);
    
    if( sps->scaling_list_enabled_flag ) {
        bs_write_u1(b, sps->sps_scaling_list_data_present_flag);
        if( sps->sps_scaling_list_data_present_flag ) {
            write_hevc_scaling_list_data(&sps->scaling_list_data, b); 
        }
    }
    
    bs_write_u1(b, sps->amp_enabled_flag);
    bs_write_u1(b, sps->sample_adaptive_offset_enabled_flag);
    bs_write_u1(b, sps->pcm_enabled_flag);
    if( sps->pcm_enabled_flag ) {
        bs_write_u(b, 4, sps->pcm_sample_bit_depth_luma_minus1);
        bs_write_u(b, 4, sps->pcm_sample_bit_depth_chroma_minus1);
        bs_write_ue(b, sps->log2_min_pcm_luma_coding_block_size_minus3);
        bs_write_ue(b, sps->log2_diff_max_min_pcm_luma_coding_block_size);
        bs_write_u1(b, sps->pcm_loop_filter_disabled_flag);
    }
    bs_write_ue(b, sps->num_short_term_ref_pic_sets);
    for( i = 0; i < sps->num_short_term_ref_pic_sets; i++) {
        write_hevc_st_ref_pic_set(&sps->st_ref_pic_set[i], b, i, sps->num_short_term_ref_pic_sets);
    }
    
    bs_write_u1(b, sps->long_term_ref_pics_present_flag);
    if( sps->long_term_ref_pics_present_flag ) {
        bs_write_ue(b, sps->num_long_term_ref_pics_sps);
        for( i = 0; i < sps->num_long_term_ref_pics_sps; i++ ) {
            bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sps->lt_ref_pic_poc_lsb_sps[ i ]);
            bs_write_u1(b, sps->used_by_curr_pic_lt_sps_flag[ i ]);
        }
    }
    bs_write_u1(b, sps->sps_temporal_mvp_enabled_flag);
    bs_write_u1(b, sps->strong_intra_smoothing_enabled_flag);
    bs_write_u1(b, sps->vui_parameters_present_flag);
    if( sps->vui_parameters_present_flag ) {
        write_hevc_vui_parameters(sps, b);
    }
    bs_write_u1(b, sps->sps_extension_present_flag);
    
    if( sps->sps_extension_present_flag ) {
        bs_write_u1(b, sps->sps_range_extension_flag);
        bs_write_u1(b, sps->sps_multilayer_extension_flag);
        bs_write_u1(b, sps->sps_3d_extension_flag);
        bs_write_u(b, 5, sps->sps_extension_5bits);
    }
    if( sps->sps_range_extension_flag ) {
        write_hevc_sps_range_extension( &sps->sps_range_ext, b);
    }
    
    if( 0 )
    {
        memcpy(h->sps_table[sps->sps_seq_parameter_set_id], h->sps, sizeof(hevc_sps_t));
    }
}

//7.3.2.2.2 Sequence parameter set range extension syntax
void write_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b)
{
    bs_write_u1(b, sps_range_ext->transform_skip_rotation_enabled_flag);
    bs_write_u1(b, sps_range_ext->transform_skip_context_enabled_flag);
    bs_write_u1(b, sps_range_ext->implicit_rdpcm_enabled_flag);
    bs_write_u1(b, sps_range_ext->explicit_rdpcm_enabled_flag);
    bs_write_u1(b, sps_range_ext->extended_precision_processing_flag);
    bs_write_u1(b, sps_range_ext->intra_smoothing_disabled_flag);
    bs_write_u1(b, sps_range_ext->high_precision_offsets_enabled_flag);
    bs_write_u1(b, sps_range_ext->persistent_rice_adaptation_enabled_flag);
    bs_write_u1(b, sps_range_ext->cabac_bypass_alignment_enabled_flag);
}


//7.3.2.3 Picture parameter set RBSP syntax
void write_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_pps_t* pps = h->pps;
    if( 0 )
    {
        memset(pps, 0, sizeof(hevc_pps_t));
    }

    bs_write_ue(b, pps->pic_parameter_set_id);
    bs_write_ue(b, pps->seq_parameter_set_id);
    bs_write_u1(b, pps->dependent_slice_segments_enabled_flag);
    bs_write_u1(b, pps->output_flag_present_flag);
    bs_write_u(b,  3 , pps->num_extra_slice_header_bits);
    bs_write_u1(b, pps->sign_data_hiding_enabled_flag);
    bs_write_u1(b, pps->cabac_init_present_flag);
    bs_write_ue(b, pps->num_ref_idx_l0_default_active_minus1);
    bs_write_ue(b, pps->num_ref_idx_l1_default_active_minus1);
    bs_write_se(b, pps->init_qp_minus26);
    bs_write_u1(b, pps->constrained_intra_pred_flag);
    bs_write_u1(b, pps->transform_skip_enabled_flag);
    bs_write_u1(b, pps->cu_qp_delta_enabled_flag);
    if( pps->cu_qp_delta_enabled_flag ) {
        bs_write_ue(b, pps->diff_cu_qp_delta_depth);
    }
    bs_write_se(b, pps->pps_cb_qp_offset);
    bs_write_se(b, pps->pps_cr_qp_offset);
    bs_write_u1(b, pps->pps_slice_chroma_qp_offsets_present_flag);
    bs_write_u1(b, pps->weighted_pred_flag);
    bs_write_u1(b, pps->weighted_bipred_flag);
    bs_write_u1(b, pps->transquant_bypass_enabled_flag);
    bs_write_u1(b, pps->tiles_enabled_flag);
    bs_write_u1(b, pps->entropy_coding_sync_enabled_flag);
    if( pps->tiles_enabled_flag ) {
        bs_write_ue(b, pps->num_tile_columns_minus1);
        bs_write_ue(b, pps->num_tile_rows_minus1);
        bs_write_u1(b, pps->uniform_spacing_flag);
        if( !pps->uniform_spacing_flag ) {
            for( i = 0; i < pps->num_tile_columns_minus1; i++ ) {
                bs_write_ue(b, pps->column_width_minus1[ i ]);
            }
            for( i = 0; i < pps->num_tile_rows_minus1; i++ ) {
                bs_write_ue(b, pps->row_height_minus1[ i ]);
            }
        }
        bs_write_u1(b, pps->loop_filter_across_tiles_enabled_flag);
    }
    bs_write_u1(b, pps->pps_loop_filter_across_slices_enabled_flag);
    bs_write_u1(b, pps->deblocking_filter_control_present_flag);
    if( pps->deblocking_filter_control_present_flag ) {
        bs_write_u1(b, pps->deblocking_filter_override_enabled_flag);
        bs_write_u1(b, pps->pps_deblocking_filter_disabled_flag);
        if( pps->pps_deblocking_filter_disabled_flag ) {
            bs_write_se(b, pps->pps_beta_offset_div2);
            bs_write_se(b, pps->pps_tc_offset_div2);
        }
    }
    bs_write_u1(b, pps->pps_scaling_list_data_present_flag);
    if( pps->pps_scaling_list_data_present_flag ) {
        write_hevc_scaling_list_data(&pps->scaling_list_data, b);
    }
    bs_write_u1(b, pps->lists_modification_present_flag);
    bs_write_ue(b, pps->log2_parallel_merge_level_minus2);
    bs_write_u1(b, pps->slice_segment_header_extension_present_flag);
    bs_write_u1(b, pps->pps_extension_present_flag);
    if( pps->pps_extension_present_flag ) {
        bs_write_u1(b, pps->pps_range_extension_flag);
        bs_write_u1(b, pps->pps_multilayer_extension_flag);
        bs_write_u1(b, pps->pps_3d_extension_flag);
        bs_write_u1(b, pps->pps_extension_5bits);
    }
    if( pps->pps_range_extension_flag ) {
        write_hevc_pps_range_extension( pps, b);
    }

    write_hevc_rbsp_trailing_bits(b);

    if( 0 )
    {
        memcpy(h->pps_table[pps->pic_parameter_set_id], h->pps, sizeof(hevc_pps_t));
    }
}

//7.3.2.3.2 Picture parameter set range extension syntax
void write_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b)
{
    hevc_pps_range_ext_t *pps_range_ext = &pps->pps_range_ext;;
    if( pps->transform_skip_enabled_flag ) {
        bs_write_ue(b, pps_range_ext->log2_max_transform_skip_block_size_minus2);
    }
    bs_write_u1(b, pps_range_ext->cross_component_prediction_enabled_flag);
    bs_write_u1(b, pps_range_ext->chroma_qp_offset_list_enabled_flag);
    if( pps_range_ext->chroma_qp_offset_list_enabled_flag ) {
        bs_write_ue(b, pps_range_ext->diff_cu_chroma_qp_offset_depth);
        bs_write_ue(b, pps_range_ext->chroma_qp_offset_list_len_minus1);
        for( int i = 0; i <= pps_range_ext->chroma_qp_offset_list_len_minus1; i++ ) {
            bs_write_se(b, pps_range_ext->cb_qp_offset_list[ i ]);
            bs_write_se(b, pps_range_ext->cr_qp_offset_list[ i ]);
        }
    }
    bs_write_ue(b, pps_range_ext->log2_sao_offset_scale_luma);
    bs_write_ue(b, pps_range_ext->log2_sao_offset_scale_chroma);
}

#ifdef HAVE_SEI
//7.3.2.4 Supplemental enhancement information RBSP syntax
void write_sei_rbsp(hevc_stream_t* h, bs_t* b)
{
    if( 0 )
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
            write_sei_message(h, b);
        } while( more_rbsp_data(h, b) );
    }

    if( 1 )
    {
        for (int i = 0; i < h->num_seis; i++) {
            h->sei = h->seis[i];
            write_sei_message(h, b);
        }
        h->sei = NULL;
    }

    write_hevc_rbsp_trailing_bits(b);
}

//7.3.5 Supplemental enhancement information message syntax
void write_sei_message(hevc_stream_t* h, bs_t* b)
{
    if( 1 )
    {
        _write_ff_coded_number(b, h->sei->payloadType);
        _write_ff_coded_number(b, h->sei->payloadSize);
    }
    if( 0 )
    {
        h->sei->payloadType = _read_ff_coded_number(b);
        h->sei->payloadSize = _read_ff_coded_number(b);
    }
    write_sei_payload( h, b, h->sei->payloadType, h->sei->payloadSize );
}
#endif

//7.3.2.5 Access unit delimiter RBSP syntax
void write_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b)
{
    bs_write_u(b, 3, h->aud->primary_pic_type);
    write_hevc_rbsp_trailing_bits(b);
}

//7.3.2.6 End of sequence RBSP syntax
void write_hevc_end_of_seq_rbsp()
{
}

//7.3.2.7 End of bitstream RBSP syntax
void write_end_of_bitstream_rbsp()
{
}

//7.3.2.8 Filler data RBSP syntax
void write_filler_data_rbsp(bs_t* b)
{
    while( bs_next_bits(b, 8) == 0xFF )
    {
        /* ff_byte */ bs_write_u(b, 8, 0xFF);
    }
    write_hevc_rbsp_trailing_bits(b);
}

//7.3.2.9 Slice segment layer RBSP syntax
void write_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b)
{
    write_hevc_slice_header(h, b);
    hevc_slice_data_rbsp_t* slice_data = h->slice_data;

    if ( slice_data != NULL )
    {
        if ( slice_data->rbsp_buf != NULL ) free( slice_data->rbsp_buf ); 
        uint8_t *sptr = b->p + (!!b->bits_left); // CABAC-specific: skip alignment bits, if there are any
        slice_data->rbsp_size = b->end - sptr;
        
        slice_data->rbsp_buf = (uint8_t*)malloc(slice_data->rbsp_size);
        memcpy( slice_data->rbsp_buf, sptr, slice_data->rbsp_size );
    }

    //write_hevc_slice_data(h, b); /* all categories of slice_data( ) syntax */
    write_hevc_rbsp_slice_trailing_bits( b );
}

//7.3.2.10 RBSP slice trailing bits syntax
void write_hevc_rbsp_slice_trailing_bits(bs_t* b)
{
    write_hevc_rbsp_trailing_bits(b);
    //while( more_rbsp_trailing_data(b) )
    //{
    //    value( cabac_zero_word, f(16, 0x0000) );
    //}
}

//7.3.2.11 RBSP trailing bits syntax
void write_hevc_rbsp_trailing_bits(bs_t* b)
{
    /* rbsp_stop_one_bit */ bs_write_u(b, 1, 1);

    while( !bs_byte_aligned(b) )
    {
        /* rbsp_alignment_zero_bit */ bs_write_u(b, 1, 0);
    }
}

//7.3.2.12 Byte alignment syntax
void write_hevc_byte_alignment(bs_t* b)
{
    /* alignment_bit_equal_to_one */ bs_write_u(b, 1, 1);

    while( !bs_byte_aligned(b) )
    {
        /* alignment_bit_equal_to_zero */ bs_write_u(b, 1, 0);
    }
}

//7.3.3 Profile, tier and level syntax
void write_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1)
{
    int i, j;
    if( profilePresentFlag ) {
        bs_write_u(b, 2, ptl->general_profile_space);
        bs_write_u1(b, ptl->general_tier_flag);
        bs_write_u(b, 5, ptl->general_profile_idc);
        for( i = 0; i < 32; i++ ) {
            bs_write_u1(b, ptl->general_profile_compatibility_flag[ i ]);
        }
        bs_write_u1(b, ptl->general_progressive_source_flag);
        bs_write_u1(b, ptl->general_interlaced_source_flag);
        bs_write_u1(b, ptl->general_non_packed_constraint_flag);
        bs_write_u1(b, ptl->general_frame_only_constraint_flag);
        if( ptl->general_profile_idc == 4 || ptl->general_profile_compatibility_flag[ 4 ] || 
            ptl->general_profile_idc == 5 || ptl->general_profile_compatibility_flag[ 5 ] || 
            ptl->general_profile_idc == 6 || ptl->general_profile_compatibility_flag[ 6 ] || 
            ptl->general_profile_idc == 7 || ptl->general_profile_compatibility_flag[ 7 ] ) {
                
            bs_write_u1(b, ptl->general_max_12bit_constraint_flag);
            bs_write_u1(b, ptl->general_max_10bit_constraint_flag);
            bs_write_u1(b, ptl->general_max_8bit_constraint_flag);
            bs_write_u1(b, ptl->general_max_422chroma_constraint_flag);
            bs_write_u1(b, ptl->general_max_420chroma_constraint_flag);
            bs_write_u1(b, ptl->general_max_monochrome_constraint_flag);
            bs_write_u1(b, ptl->general_intra_constraint_flag);
            bs_write_u1(b, ptl->general_one_picture_only_constraint_flag);
            bs_write_u1(b, ptl->general_lower_bit_rate_constraint_flag);
            /* general_reserved_zero_34bits */ bs_write_u(b, 34, 0);
        } else {
            /* general_reserved_zero_43bits */ bs_write_u(b, 43, 0);
        }
        if( ( ptl->general_profile_idc >= 1 && ptl->general_profile_idc <= 5 ) ||
              ptl->general_profile_compatibility_flag[ 1 ] ||
              ptl->general_profile_compatibility_flag[ 2 ] ||
              ptl->general_profile_compatibility_flag[ 3 ] ||
              ptl->general_profile_compatibility_flag[ 4 ] ||
              ptl->general_profile_compatibility_flag[ 5 ] ) {

            bs_write_u1(b, ptl->general_inbld_flag);
        } else {
            /* general_reserved_zero_bit */ bs_write_u(b, 1, 0);
        }
        bs_write_u8(b, ptl->general_level_idc);
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) {
            bs_write_u1(b, ptl->sub_layer_profile_present_flag[ i ]);
            bs_write_u1(b, ptl->sub_layer_level_present_flag[ i ]);
        }
        if( maxNumSubLayersMinus1 > 0 ) {
            for( i = maxNumSubLayersMinus1; i < 8; i++ ) {
                /* reserved_zero_xxbits */ bs_write_u(b, 2, 0);
            }
        }
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) { 
            if( ptl->sub_layer_profile_present_flag[ i ] ) {
                bs_write_u(b, 2, ptl->sub_layer_profile_space[ i ]);
                bs_write_u1(b, ptl->sub_layer_tier_flag[ i ]);
                bs_write_u(b, 5, ptl->sub_layer_profile_idc[ i ]);
                for( j = 0; j < 32; j++ ) {
                    bs_write_u(b, 1, ptl->sub_layer_profile_compatibility_flag[ i ][ j ]);
                }
                bs_write_u1(b, ptl->sub_layer_progressive_source_flag[ i ]);
                bs_write_u1(b, ptl->sub_layer_interlaced_source_flag[ i ]);
                bs_write_u1(b, ptl->sub_layer_non_packed_constraint_flag[ i ]);
                bs_write_u1(b, ptl->sub_layer_frame_only_constraint_flag[ i ]);
                if( ptl->sub_layer_profile_idc[ i ] == 4 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 4 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 5 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 5 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 6 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 6 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 7 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 7 ] ) {
                    bs_write_u1(b, ptl->sub_layer_max_12bit_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_max_10bit_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_max_8bit_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_max_422chroma_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_max_420chroma_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_max_monochrome_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_intra_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_one_picture_only_constraint_flag[ i ]);
                    bs_write_u1(b, ptl->sub_layer_lower_bit_rate_constraint_flag[ i ]);
                    /* sub_layer_reserved_zero_34bits */ bs_write_u(b, 34, 0);
                } else {
                    /* sub_layer_reserved_zero_43bits */ bs_write_u(b, 43, 0);
                }
            
                if( ( ptl->sub_layer_profile_idc[ i ] >= 1 && ptl->sub_layer_profile_idc[ i ] <= 5 ) ||
                   ptl->sub_layer_profile_compatibility_flag[ 1 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 2 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 3 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 4 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 5 ] ) {
                    bs_write_u1(b, ptl->sub_layer_inbld_flag[ i ]);
                } else {
                    /* sub_layer_reserved_zero_bit */ bs_write_u(b, 1, 0);
                }
            }
            if( ptl->sub_layer_level_present_flag[ i ] ) {
                bs_write_u1(b, ptl->sub_layer_level_idc[ i ]);
            }
        }
    }
}

//7.3.4 Scaling list data syntax
void write_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b )
{
    int nextCoef, coefNum;
    for( int sizeId = 0; sizeId < 4; sizeId++ )
        for( int matrixId = 0; matrixId < 6; matrixId += ( sizeId == 3 ) ? 3 : 1 ) {
            bs_write_u1(b, sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ]);
            if( !sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] ) {
                bs_write_ue(b, sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ]);
            } else {
                nextCoef = 8;
                coefNum=MIN(64, (1 << (4+(sizeId << 1))));
                if( sizeId > 1 ) {
                    bs_write_se(b, sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ]);
                }
 
                for( int i = 0; i < coefNum; i++) {
                    bs_write_se(b, sld->scaling_list_delta_coef[ sizeId ][ matrixId ]);
                }
            }
 
        }
}

//7.3.6 Slice header syntax
void write_hevc_slice_header(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    if( 0 )
    {
        init_slice_hevc(h);
    }

    hevc_nal_t* nal = h->nal;

    bs_write_u1(b, sh->first_slice_segment_in_pic_flag);
    if( nal->nal_unit_type >= HEVC_NAL_UNIT_TYPE_BLA_W_LP && 
       nal->nal_unit_type <= HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23) {
            bs_write_u1(b, sh->no_output_of_prior_pics_flag);
    }
    bs_write_ue(b, sh->pic_parameter_set_id);
    
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];

    //set default value
    sh->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    sh->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;

    if( !sh->first_slice_segment_in_pic_flag ) {
        if( pps->dependent_slice_segments_enabled_flag ) {
            bs_write_u1(b, sh->dependent_slice_segment_flag);
        }
        bs_write_u(b,  getSliceSegmentAddressBitLength( sps ) , sh->slice_segment_address);
    }
    
    if( !sh->dependent_slice_segment_flag ) {
        for( i = 0; i < pps->num_extra_slice_header_bits; i++ ) {
            /* slice_reserved_flag */ bs_write_u(b, 1, 1);
        }
        bs_write_ue(b, sh->slice_type);
        if( pps->output_flag_present_flag ) {
            bs_write_u1(b, sh->pic_output_flag);
        }
        if( sps->separate_colour_plane_flag == 1 ) {
            bs_write_u(b, 2, sh->colour_plane_id);
        }
        if( nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_W_RADL &&
            nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_N_LP) {
            bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sh->slice_pic_order_cnt_lsb);
            bs_write_u1(b, sh->short_term_ref_pic_set_sps_flag);
            if( !sh->short_term_ref_pic_set_sps_flag ) {
                write_hevc_st_ref_pic_set( &sh->st_ref_pic_set, b, sps->num_short_term_ref_pic_sets, sps->num_short_term_ref_pic_sets );
            } else if( sps->num_short_term_ref_pic_sets > 1 ) {
                bs_write_u(b,  ceil( log2( sps->num_short_term_ref_pic_sets ) ) , sh->short_term_ref_pic_set_idx);
            }
            if( sps->long_term_ref_pics_present_flag ) {
                if( sps->num_long_term_ref_pics_sps > 0 ) {
                    bs_write_ue(b, sh->num_long_term_sps);
                }
                bs_write_ue(b, sh->num_long_term_pics);
                for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
                    if( i < sh->num_long_term_sps ) {
                        if( sps->num_long_term_ref_pics_sps > 1 ) {
                            bs_write_u(b,  ceil( log2( sps->num_long_term_ref_pics_sps ) ) , sh->lt_idx_sps[ i ]);
                        }
                    } else {
                        bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sh->poc_lsb_lt[ i ]);
                        bs_write_u1(b, sh->used_by_curr_pic_lt_flag[ i ]);
                    }
                    bs_write_u1(b, sh->delta_poc_msb_present_flag[ i ]);
                    if( sh->delta_poc_msb_present_flag[ i ]) {
                        bs_write_ue(b, sh->delta_poc_msb_cycle_lt[ i ]);
                    }
                }
            }
            if( sps->sps_temporal_mvp_enabled_flag ) {
                bs_write_u1(b, sh->slice_temporal_mvp_enabled_flag);
            }
        }
        if( sps->sample_adaptive_offset_enabled_flag ) {
            bs_write_u1(b, sh->slice_sao_luma_flag);
            int ChromaArrayType = 0;
            if( sps->separate_colour_plane_flag == 0) {
                ChromaArrayType = sps->chroma_format_idc;
            }
            if( ChromaArrayType != 0 ) {
                bs_write_u1(b, sh->slice_sao_chroma_flag);
            }
        }
        if( sh->slice_type == HEVC_SLICE_TYPE_P || sh->slice_type == HEVC_SLICE_TYPE_B ){
            bs_write_u1(b, sh->num_ref_idx_active_override_flag);
            if( sh->num_ref_idx_active_override_flag ) {
                bs_write_ue(b, sh->num_ref_idx_l0_active_minus1);
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    bs_write_ue(b, sh->num_ref_idx_l1_active_minus1);
                }
            }
            if( pps->lists_modification_present_flag && getNumPicTotalCurr( sps, sh ) > 1 ) {
                write_hevc_ref_pic_lists_modification( h, b );
            }
            
            if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                bs_write_u1(b, sh->mvd_l1_zero_flag);
            }
            if( pps->cabac_init_present_flag ) {
                bs_write_u1(b, sh->cabac_init_flag);
            }
            if( sh->slice_temporal_mvp_enabled_flag ) {
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    bs_write_u1(b, sh->collocated_from_l0_flag);
                }
                if( ( sh->collocated_from_l0_flag && sh->num_ref_idx_l0_active_minus1 > 0 ) ||
                   ( !sh->collocated_from_l0_flag && sh->num_ref_idx_l1_active_minus1 > 0 ) ) {
                    bs_write_ue(b, sh->collocated_ref_idx);
                }
            }
            if( ( pps->weighted_pred_flag && sh->slice_type == HEVC_SLICE_TYPE_P ) || 
                ( pps->weighted_bipred_flag && sh->slice_type == HEVC_SLICE_TYPE_B ) ) {
                write_hevc_pred_weight_table( h, b );
            }
            bs_write_ue(b, sh->five_minus_max_num_merge_cand);
        }
        bs_write_se(b, sh->slice_qp_delta);
        if( pps->pps_slice_chroma_qp_offsets_present_flag ) {
            bs_write_se(b, sh->slice_cb_qp_offset);
            bs_write_se(b, sh->slice_cr_qp_offset);
        }
        if( pps->pps_range_ext.chroma_qp_offset_list_enabled_flag ) {
            bs_write_u1(b, sh->cu_chroma_qp_offset_enabled_flag);
        }
        if( pps->deblocking_filter_override_enabled_flag ) {
            bs_write_u1(b, sh->deblocking_filter_override_flag);
        }
        if( sh->deblocking_filter_override_flag ) {
            bs_write_u1(b, sh->slice_deblocking_filter_disabled_flag);
            if( !sh->slice_deblocking_filter_disabled_flag ) {
                bs_write_se(b, sh->slice_beta_offset_div2);
                bs_write_se(b, sh->slice_tc_offset_div2);
            }
        }
        if( pps->pps_loop_filter_across_slices_enabled_flag &&
           ( sh->slice_sao_luma_flag || sh->slice_sao_chroma_flag || !sh->slice_deblocking_filter_disabled_flag ) ) {
            bs_write_u1(b, sh->slice_loop_filter_across_slices_enabled_flag);
        }
    }
    if( pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag ) {
        bs_write_ue(b, sh->num_entry_point_offsets);
        if( sh->num_entry_point_offsets > 0 ) {
            bs_write_ue(b, sh->offset_len_minus1);
            for( i = 0; i < sh->num_entry_point_offsets; i++ ) {
                bs_write_u(b,  sh->offset_len_minus1 + 1 , sh->entry_point_offset_minus1[ i ]);
            }
        }
    }
    if( pps->slice_segment_header_extension_present_flag ) {
        bs_write_ue(b, sh->slice_segment_header_extension_length);
        //TODO: support header extension,
        for( i = 0; i < sh->slice_segment_header_extension_length; i++) {
            /* slice_segment_header_extension_data_byte */ bs_write_u(b, 8, 0);
        }
    }
    write_hevc_byte_alignment( b );
}

//7.3.6.2 Reference picture list reordering syntax
void write_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    
    bs_write_u1(b, sh->rpld.ref_pic_list_modification_flag_l0);
    if( sh->rpld.ref_pic_list_modification_flag_l0 ) {
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            bs_write_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) , sh->rpld.list_entry_l0[ i ]);
        }
    }
    
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        // ERROR: value( sh->rpld.ref_pic_list_modification_flag_l1, 1 );
        if( sh->rpld.ref_pic_list_modification_flag_l1 ) {
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                bs_write_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) , sh->rpld.list_entry_l1[ i ]);
            }
        }
    }
}

//7.3.6.3 Prediction weight table syntax
void write_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b)
{
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    hevc_pred_weight_table_t *pwt = &sh->pwt;

    int i, j;

    bs_write_ue(b, pwt->luma_log2_weight_denom);
    
    int ChromaArrayType = 0;
    if( sps->separate_colour_plane_flag == 0) {
        ChromaArrayType = sps->chroma_format_idc;
    }

    if( ChromaArrayType != 0 ) {
        bs_write_se(b, pwt->delta_chroma_log2_weight_denom);
    }
    
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        bs_write_u1(b, pwt->luma_weight_l0_flag[i]);
    }
    if( ChromaArrayType != 0 )
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            bs_write_u1(b, pwt->chroma_weight_l0_flag[i]);
        }
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        if( pwt->luma_weight_l0_flag[i] ) {
            bs_write_se(b, pwt->delta_luma_weight_l0[i]);
            bs_write_se(b, pwt->luma_offset_l0[i]);
        }
        if( pwt->chroma_weight_l0_flag[i] ) {
            for( j =0; j < 2; j++ ) {
                bs_write_se(b, pwt->delta_chroma_weight_l0[i][j]);
                bs_write_se(b, pwt->delta_chroma_offset_l0[i][j]);
            }
        }
    }
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            bs_write_u1(b, pwt->luma_weight_l1_flag[i]);
        }
        if( ChromaArrayType != 0 )
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                bs_write_u1(b, pwt->chroma_weight_l1_flag[i]);
            }
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            if( pwt->luma_weight_l1_flag[i] ) {
                bs_write_se(b, pwt->delta_luma_weight_l1[i]);
                bs_write_se(b, pwt->luma_offset_l1[i]);
            }
            if( pwt->chroma_weight_l1_flag[i] ) {
                for( j =0; j < 2; j++ ) {
                    bs_write_se(b, pwt->delta_chroma_weight_l1[i][j]);
                    bs_write_se(b, pwt->delta_chroma_offset_l1[i][j]);
                }
            }
        }        
    }
}

//7.3.7 Short-term reference picture set syntax
void write_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets )
{
    int i, j;
    
    if( stRpsIdx != 0 ) {
        bs_write_u1(b, st_ref_pic_set->inter_ref_pic_set_prediction_flag);
    }
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        if( stRpsIdx == num_short_term_ref_pic_sets ) {
            bs_write_ue(b, st_ref_pic_set->delta_idx_minus1);
        }
        bs_write_u1(b, st_ref_pic_set->delta_rps_sign);
        bs_write_ue(b, st_ref_pic_set->abs_delta_rps_minus1);
        
        int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
        
        for( j = 0; j <= NumDeltaPocs[ RefRpsIdx ]; j++ ) {
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_flag[ j ]);
            if( !st_ref_pic_set->used_by_curr_pic_flag[ j ] ) {
                bs_write_u1(b, st_ref_pic_set->use_delta_flag[ j ]);
            }
        }
    } else {
        bs_write_ue(b, st_ref_pic_set->num_negative_pics);
        bs_write_ue(b, st_ref_pic_set->num_positive_pics);
        for( i = 0; i < st_ref_pic_set->num_negative_pics; i++ ) {
            bs_write_ue(b, st_ref_pic_set->delta_poc_s0_minus1[ i ]);
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_s0_flag[ i ]);
            
            //update derived field
            UsedByCurrPicS0[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s0_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS0[ stRpsIdx ][ i ] = -1 * ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            } else {
                DeltaPocS0[ stRpsIdx ][ i ] = DeltaPocS0[ stRpsIdx ][ i - 1 ] - ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            }
        }
        for( i = 0; i < st_ref_pic_set->num_positive_pics; i++ ) {
            bs_write_ue(b, st_ref_pic_set->delta_poc_s1_minus1[ i ]);
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_s1_flag[ i ]);
        
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
void write_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b)
{
    hevc_vui_t* vui = &sps->vui;
    bs_write_u1(b, vui->aspect_ratio_info_present_flag);
    if( vui->aspect_ratio_info_present_flag )
    {
        bs_write_u8(b, vui->aspect_ratio_idc);
        if( vui->aspect_ratio_idc == SAR_Extended )
        {
            bs_write_u(b, 16, vui->sar_width);
            bs_write_u(b, 16, vui->sar_height);
        }
    }
    bs_write_u1(b, vui->overscan_info_present_flag);
    if( vui->overscan_info_present_flag ) {
        bs_write_u1(b, vui->overscan_appropriate_flag);
    }
    bs_write_u1(b, vui->video_signal_type_present_flag);
    if( vui->video_signal_type_present_flag ) {
        bs_write_u(b, 3, vui->video_format);
        bs_write_u1(b, vui->video_full_range_flag);
        bs_write_u1(b, vui->colour_description_present_flag);
        if( vui->colour_description_present_flag ) {
            bs_write_u8(b, vui->colour_primaries);
            bs_write_u8(b, vui->transfer_characteristics);
            bs_write_u8(b, vui->matrix_coefficients);
        }
    }
    bs_write_u1(b, vui->chroma_loc_info_present_flag);
    if( vui->chroma_loc_info_present_flag ) {
        bs_write_ue(b, vui->chroma_sample_loc_type_top_field);
        bs_write_ue(b, vui->chroma_sample_loc_type_bottom_field);
    }
    
    bs_write_u1(b, vui->neutral_chroma_indication_flag);
    bs_write_u1(b, vui->field_seq_flag);
    bs_write_u1(b, vui->frame_field_info_present_flag);
    bs_write_u1(b, vui->default_display_window_flag);
    if( vui->default_display_window_flag ) {
        bs_write_ue(b, vui->def_disp_win_left_offset);
        bs_write_ue(b, vui->def_disp_win_right_offset);
        bs_write_ue(b, vui->def_disp_win_top_offset);
        bs_write_ue(b, vui->def_disp_win_bottom_offset);
    }
    bs_write_u1(b, vui->vui_timing_info_present_flag);
    if( vui->vui_timing_info_present_flag ) {
        bs_write_u(b, 32, vui->vui_num_units_in_tick);
        bs_write_u(b, 32, vui->vui_time_scale);
        bs_write_u1(b, vui->vui_poc_proportional_to_timing_flag);
        if( vui->vui_poc_proportional_to_timing_flag ) {
            bs_write_ue(b, vui->vui_num_ticks_poc_diff_one_minus1);
        }
        bs_write_u1(b, vui->vui_hrd_parameters_present_flag);
        if( vui->vui_hrd_parameters_present_flag ) {
            write_hevc_hrd_parameters( &vui->hrd, b, 1, sps->sps_max_sub_layers_minus1 );
        }
    }
    bs_write_u1(b, vui->bitstream_restriction_flag);
    if( vui->bitstream_restriction_flag )
    {
        bs_write_u1(b, vui->tiles_fixed_structure_flag);
        bs_write_u1(b, vui->motion_vectors_over_pic_boundaries_flag);
        bs_write_u1(b, vui->restricted_ref_pic_lists_flag);
        bs_write_ue(b, vui->min_spatial_segmentation_idc);
        bs_write_ue(b, vui->max_bytes_per_pic_denom);
        bs_write_ue(b, vui->max_bits_per_min_cu_denom);
        bs_write_ue(b, vui->log2_max_mv_length_horizontal);
        bs_write_ue(b, vui->log2_max_mv_length_vertical);
    }
}

//Appendix E.2.2 HRD parameters syntax
void write_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1)
{
    if( commonInfPresentFlag ) {
        bs_write_u1(b, hrd->nal_hrd_parameters_present_flag);
        bs_write_u1(b, hrd->vcl_hrd_parameters_present_flag);
        if( hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag ){
            bs_write_u1(b, hrd->sub_pic_hrd_params_present_flag);
            if( hrd->sub_pic_hrd_params_present_flag ) {
                bs_write_u8(b, hrd->tick_divisor_minus2);
                bs_write_u(b, 5, hrd->du_cpb_removal_delay_increment_length_minus1);
                bs_write_u1(b, hrd->sub_pic_cpb_params_in_pic_timing_sei_flag);
                bs_write_u(b, 5, hrd->dpb_output_delay_du_length_minus1);
            }
            bs_write_u(b, 4, hrd->bit_rate_scale);
            bs_write_u(b, 4, hrd->cpb_size_scale);
            if( hrd->sub_pic_hrd_params_present_flag ) {
                bs_write_u(b, 4, hrd->cpb_size_du_scale);
            }
            bs_write_u(b, 5, hrd->initial_cpb_removal_delay_length_minus1);
            bs_write_u(b, 5, hrd->au_cpb_removal_delay_length_minus1);
            bs_write_u(b, 5, hrd->dpb_output_delay_length_minus1);
        }
    }
    
    for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
        bs_write_u1(b, hrd->fixed_pic_rate_general_flag[ i ]);
        if( !hrd->fixed_pic_rate_general_flag[ i ] ) {
            bs_write_u1(b, hrd->fixed_pic_rate_within_cvs_flag[ i ]);
        }
        if( hrd->fixed_pic_rate_within_cvs_flag[ i ] ) {
            bs_write_ue(b, hrd->elemental_duration_in_tc_minus1[ i ]);
        } else {
            bs_write_u1(b, hrd->low_delay_hrd_flag[ i ]);
        }
        if( hrd->low_delay_hrd_flag[ i ] ) {
            bs_write_ue(b, hrd->cpb_cnt_minus1[ i ]);
        }
        if( hrd->nal_hrd_parameters_present_flag ) {
            write_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_nal[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
        if( hrd->vcl_hrd_parameters_present_flag ) {
            write_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_vcl[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
    }
}

//Appendix E.2.3 Sub-layer HRD parameters syntax
void write_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag)
{
    for( int i = 0; i <= CpbCnt; i++ ) {
        bs_write_ue(b, sub_layer_hrd->bit_rate_value_minus1[i]);
        bs_write_ue(b, sub_layer_hrd->cpb_size_value_minus1[i]);
        if( sub_pic_hrd_params_present_flag ) {
            bs_write_ue(b, sub_layer_hrd->cpb_size_du_value_minus1[i]);
            bs_write_ue(b, sub_layer_hrd->bit_rate_du_value_minus1[i]);
        }
        bs_write_u1(b, sub_layer_hrd->cbr_flag[i]);
    }
}


void read_debug_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b);
void read_debug_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b);
void read_debug_sei_rbsp(hevc_stream_t* h, bs_t* b);
void read_debug_sei_message(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_end_of_seq_rbsp();
void read_debug_end_of_bitstream_rbsp();
void read_debug_filler_data_rbsp(bs_t* b);
void read_debug_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b);
void read_debug_hevc_rbsp_slice_trailing_bits(bs_t* b);
void read_debug_hevc_rbsp_trailing_bits(bs_t* b);
void read_debug_hevc_byte_alignment(bs_t* b);
void read_debug_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1);
void read_debug_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b );
void read_debug_hevc_slice_header(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b);
void read_debug_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets );
void read_debug_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b);
void read_debug_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1);
void read_debug_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag);



//7.3.1 NAL unit syntax
int read_debug_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    int nal_size = size;
    int rbsp_size = size;
    uint8_t* rbsp_buf = (uint8_t*)calloc(1, rbsp_size);

    if( 1 )
    {
        int rc = nal_to_rbsp(buf, &nal_size, rbsp_buf, &rbsp_size);

        if (rc < 0) { free(rbsp_buf); return -1; } // handle conversion error
    }

    if( 0 )
    {
        rbsp_size = size*3/4; // NOTE this may have to be slightly smaller (3/4 smaller, worst case) in order to be guaranteed to fit
    }

    bs_t* b = bs_new(rbsp_buf, rbsp_size);
    printf("%ld.%d: ", b->p - b->start, b->bits_left); int forbidden_zero_bit = bs_read_u(b, 1); printf("forbidden_zero_bit: %d \n", forbidden_zero_bit); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); nal->nal_unit_type = bs_read_u(b, 6); printf("nal->nal_unit_type: %d \n", nal->nal_unit_type); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); nal->nal_layer_id = bs_read_u(b, 6); printf("nal->nal_layer_id: %d \n", nal->nal_layer_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); nal->nal_temporal_id_plus1 = bs_read_u(b, 3); printf("nal->nal_temporal_id_plus1: %d \n", nal->nal_temporal_id_plus1); 

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
            
            read_debug_hevc_slice_layer_rbsp(h, b);
            break;

#ifdef HAVE_SEI
        case NAL_UNIT_TYPE_SEI:
            read_debug_hevc_sei_rbsp(h, b);
            break;
#endif

        case HEVC_NAL_UNIT_TYPE_VPS_NUT: 
            read_debug_hevc_video_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_SPS_NUT: 
            read_debug_hevc_seq_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_PPS_NUT:   
            read_debug_hevc_pic_parameter_set_rbsp(h, b);
            break;

        default:
            return -1;
    }

    if (bs_overrun(b)) { bs_free(b); free(rbsp_buf); return -1; }

    if( 0 )
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
void read_debug_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i, j;

    hevc_vps_t* vps = h->vps;
    if( 1 )
    {
        memset(vps, 0, sizeof(hevc_vps_t));
    }
 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_video_parameter_set_id = bs_read_u(b, 4); printf("vps->vps_video_parameter_set_id: %d \n", vps->vps_video_parameter_set_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_base_layer_internal_flag = bs_read_u1(b); printf("vps->vps_base_layer_internal_flag: %d \n", vps->vps_base_layer_internal_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_base_layer_available_flag = bs_read_u1(b); printf("vps->vps_base_layer_available_flag: %d \n", vps->vps_base_layer_available_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_layers_minus1 = bs_read_u(b, 6); printf("vps->vps_max_layers_minus1: %d \n", vps->vps_max_layers_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_sub_layers_minus1 = bs_read_u(b, 3); printf("vps->vps_max_sub_layers_minus1: %d \n", vps->vps_max_sub_layers_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_temporal_id_nesting_flag = bs_read_u1(b); printf("vps->vps_temporal_id_nesting_flag: %d \n", vps->vps_temporal_id_nesting_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); int vps_reserved_0xffff_16bits = bs_read_u(b, 16); printf("vps_reserved_0xffff_16bits: %d \n", vps_reserved_0xffff_16bits); 
    
    read_debug_hevc_profile_tier_level(&vps->ptl, b, 1, vps->vps_max_sub_layers_minus1); 
    
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_sub_layer_ordering_info_present_flag = bs_read_u1(b); printf("vps->vps_sub_layer_ordering_info_present_flag: %d \n", vps->vps_sub_layer_ordering_info_present_flag); 
    for( i = ( vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1 ); 
            i <= vps->vps_max_sub_layers_minus1; i++ ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_dec_pic_buffering_minus1[ i ] = bs_read_ue(b); printf("vps->vps_max_dec_pic_buffering_minus1[ i ]: %d \n", vps->vps_max_dec_pic_buffering_minus1[ i ]);         
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_num_reorder_pics[ i ] = bs_read_ue(b); printf("vps->vps_max_num_reorder_pics[ i ]: %d \n", vps->vps_max_num_reorder_pics[ i ]);         
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_latency_increase_plus1[ i ] = bs_read_ue(b); printf("vps->vps_max_latency_increase_plus1[ i ]: %d \n", vps->vps_max_latency_increase_plus1[ i ]);         
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_max_layer_id = bs_read_u(b, 6); printf("vps->vps_max_layer_id: %d \n", vps->vps_max_layer_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_num_layer_sets_minus1 = bs_read_ue(b); printf("vps->vps_num_layer_sets_minus1: %d \n", vps->vps_num_layer_sets_minus1); 
    for( i = 1; i <= vps->vps_num_layer_sets_minus1; i++ )
        for( j = 0; j <= vps->vps_max_layer_id; j++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->layer_id_included_flag[ i ][ j ] = bs_read_u1(b); printf("vps->layer_id_included_flag[ i ][ j ]: %d \n", vps->layer_id_included_flag[ i ][ j ]); 
        }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_timing_info_present_flag = bs_read_u1(b); printf("vps->vps_timing_info_present_flag: %d \n", vps->vps_timing_info_present_flag); 
    if( vps->vps_timing_info_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_num_units_in_tick = bs_read_u(b, 32); printf("vps->vps_num_units_in_tick: %d \n", vps->vps_num_units_in_tick); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_time_scale = bs_read_u(b, 32); printf("vps->vps_time_scale: %d \n", vps->vps_time_scale); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_poc_proportional_to_timing_flag = bs_read_u1(b); printf("vps->vps_poc_proportional_to_timing_flag: %d \n", vps->vps_poc_proportional_to_timing_flag); 
        if( vps->vps_poc_proportional_to_timing_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_num_ticks_poc_diff_one_minus1 = bs_read_ue(b); printf("vps->vps_num_ticks_poc_diff_one_minus1: %d \n", vps->vps_num_ticks_poc_diff_one_minus1); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_num_hrd_parameters = bs_read_ue(b); printf("vps->vps_num_hrd_parameters: %d \n", vps->vps_num_hrd_parameters); 
        for( i = 0; i < vps->vps_num_hrd_parameters; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->hrd_layer_set_idx[ i ] = bs_read_ue(b); printf("vps->hrd_layer_set_idx[ i ]: %d \n", vps->hrd_layer_set_idx[ i ]); 
            if (i > 0) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->cprms_present_flag[ i ] = bs_read_u1(b); printf("vps->cprms_present_flag[ i ]: %d \n", vps->cprms_present_flag[ i ]); 
            }
            read_debug_hevc_hrd_parameters(&vps->hrd[i], b,
                                           vps->cprms_present_flag[ i ],
                                           vps->vps_max_sub_layers_minus1);
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vps->vps_extension_flag = bs_read_u1(b); printf("vps->vps_extension_flag: %d \n", vps->vps_extension_flag); 
    //TODO: support extension data
    //if (vps->vps_extension_flag)    

    read_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.2 Sequence parameter set RBSP syntax
void read_debug_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;

    hevc_sps_t* sps = h->sps;
    if( 1 )
    {
        memset(sps, 0, sizeof(hevc_sps_t));
    }
 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_video_parameter_set_id = bs_read_u(b, 4); printf("sps->sps_video_parameter_set_id: %d \n", sps->sps_video_parameter_set_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_max_sub_layers_minus1 = bs_read_u(b, 3); printf("sps->sps_max_sub_layers_minus1: %d \n", sps->sps_max_sub_layers_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_temporal_id_nesting_flag = bs_read_u1(b); printf("sps->sps_temporal_id_nesting_flag: %d \n", sps->sps_temporal_id_nesting_flag); 
    read_debug_hevc_profile_tier_level(&sps->ptl, b, 1, sps->sps_max_sub_layers_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_seq_parameter_set_id = bs_read_ue(b); printf("sps->sps_seq_parameter_set_id: %d \n", sps->sps_seq_parameter_set_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->chroma_format_idc = bs_read_ue(b); printf("sps->chroma_format_idc: %d \n", sps->chroma_format_idc); 
    if( sps->chroma_format_idc == 3 ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->separate_colour_plane_flag = bs_read_u1(b); printf("sps->separate_colour_plane_flag: %d \n", sps->separate_colour_plane_flag); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pic_width_in_luma_samples = bs_read_ue(b); printf("sps->pic_width_in_luma_samples: %d \n", sps->pic_width_in_luma_samples); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pic_height_in_luma_samples = bs_read_ue(b); printf("sps->pic_height_in_luma_samples: %d \n", sps->pic_height_in_luma_samples); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->conformance_window_flag = bs_read_u1(b); printf("sps->conformance_window_flag: %d \n", sps->conformance_window_flag); 
    if( sps->conformance_window_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->conf_win_left_offset = bs_read_ue(b); printf("sps->conf_win_left_offset: %d \n", sps->conf_win_left_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->conf_win_right_offset = bs_read_ue(b); printf("sps->conf_win_right_offset: %d \n", sps->conf_win_right_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->conf_win_top_offset = bs_read_ue(b); printf("sps->conf_win_top_offset: %d \n", sps->conf_win_top_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->conf_win_bottom_offset = bs_read_ue(b); printf("sps->conf_win_bottom_offset: %d \n", sps->conf_win_bottom_offset); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->bit_depth_luma_minus8 = bs_read_ue(b); printf("sps->bit_depth_luma_minus8: %d \n", sps->bit_depth_luma_minus8); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->bit_depth_chroma_minus8 = bs_read_ue(b); printf("sps->bit_depth_chroma_minus8: %d \n", sps->bit_depth_chroma_minus8); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_max_pic_order_cnt_lsb_minus4 = bs_read_ue(b); printf("sps->log2_max_pic_order_cnt_lsb_minus4: %d \n", sps->log2_max_pic_order_cnt_lsb_minus4); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_sub_layer_ordering_info_present_flag = bs_read_u1(b); printf("sps->sps_sub_layer_ordering_info_present_flag: %d \n", sps->sps_sub_layer_ordering_info_present_flag); 
    for( i = ( sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1 ); 
            i <= sps->sps_max_sub_layers_minus1; i++ ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_max_dec_pic_buffering_minus1 [ i ] = bs_read_ue(b); printf("sps->sps_max_dec_pic_buffering_minus1 [ i ]: %d \n", sps->sps_max_dec_pic_buffering_minus1 [ i ]); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_max_num_reorder_pics [ i ] = bs_read_ue(b); printf("sps->sps_max_num_reorder_pics [ i ]: %d \n", sps->sps_max_num_reorder_pics [ i ]); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_max_latency_increase_plus1 [ i ] = bs_read_ue(b); printf("sps->sps_max_latency_increase_plus1 [ i ]: %d \n", sps->sps_max_latency_increase_plus1 [ i ]); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_min_luma_coding_block_size_minus3 = bs_read_ue(b); printf("sps->log2_min_luma_coding_block_size_minus3: %d \n", sps->log2_min_luma_coding_block_size_minus3); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_diff_max_min_luma_coding_block_size = bs_read_ue(b); printf("sps->log2_diff_max_min_luma_coding_block_size: %d \n", sps->log2_diff_max_min_luma_coding_block_size); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_min_luma_transform_block_size_minus2 = bs_read_ue(b); printf("sps->log2_min_luma_transform_block_size_minus2: %d \n", sps->log2_min_luma_transform_block_size_minus2); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_diff_max_min_luma_transform_block_size = bs_read_ue(b); printf("sps->log2_diff_max_min_luma_transform_block_size: %d \n", sps->log2_diff_max_min_luma_transform_block_size); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->max_transform_hierarchy_depth_inter = bs_read_ue(b); printf("sps->max_transform_hierarchy_depth_inter: %d \n", sps->max_transform_hierarchy_depth_inter); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->max_transform_hierarchy_depth_intra = bs_read_ue(b); printf("sps->max_transform_hierarchy_depth_intra: %d \n", sps->max_transform_hierarchy_depth_intra); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->scaling_list_enabled_flag = bs_read_u1(b); printf("sps->scaling_list_enabled_flag: %d \n", sps->scaling_list_enabled_flag); 
    
    if( sps->scaling_list_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_scaling_list_data_present_flag = bs_read_u1(b); printf("sps->sps_scaling_list_data_present_flag: %d \n", sps->sps_scaling_list_data_present_flag); 
        if( sps->sps_scaling_list_data_present_flag ) {
            read_debug_hevc_scaling_list_data(&sps->scaling_list_data, b); 
        }
    }
    
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->amp_enabled_flag = bs_read_u1(b); printf("sps->amp_enabled_flag: %d \n", sps->amp_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sample_adaptive_offset_enabled_flag = bs_read_u1(b); printf("sps->sample_adaptive_offset_enabled_flag: %d \n", sps->sample_adaptive_offset_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pcm_enabled_flag = bs_read_u1(b); printf("sps->pcm_enabled_flag: %d \n", sps->pcm_enabled_flag); 
    if( sps->pcm_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pcm_sample_bit_depth_luma_minus1 = bs_read_u(b, 4); printf("sps->pcm_sample_bit_depth_luma_minus1: %d \n", sps->pcm_sample_bit_depth_luma_minus1); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pcm_sample_bit_depth_chroma_minus1 = bs_read_u(b, 4); printf("sps->pcm_sample_bit_depth_chroma_minus1: %d \n", sps->pcm_sample_bit_depth_chroma_minus1); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_min_pcm_luma_coding_block_size_minus3 = bs_read_ue(b); printf("sps->log2_min_pcm_luma_coding_block_size_minus3: %d \n", sps->log2_min_pcm_luma_coding_block_size_minus3); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->log2_diff_max_min_pcm_luma_coding_block_size = bs_read_ue(b); printf("sps->log2_diff_max_min_pcm_luma_coding_block_size: %d \n", sps->log2_diff_max_min_pcm_luma_coding_block_size); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->pcm_loop_filter_disabled_flag = bs_read_u1(b); printf("sps->pcm_loop_filter_disabled_flag: %d \n", sps->pcm_loop_filter_disabled_flag); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->num_short_term_ref_pic_sets = bs_read_ue(b); printf("sps->num_short_term_ref_pic_sets: %d \n", sps->num_short_term_ref_pic_sets); 
    for( i = 0; i < sps->num_short_term_ref_pic_sets; i++) {
        read_debug_hevc_st_ref_pic_set(&sps->st_ref_pic_set[i], b, i, sps->num_short_term_ref_pic_sets);
    }
    
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->long_term_ref_pics_present_flag = bs_read_u1(b); printf("sps->long_term_ref_pics_present_flag: %d \n", sps->long_term_ref_pics_present_flag); 
    if( sps->long_term_ref_pics_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->num_long_term_ref_pics_sps = bs_read_ue(b); printf("sps->num_long_term_ref_pics_sps: %d \n", sps->num_long_term_ref_pics_sps); 
        for( i = 0; i < sps->num_long_term_ref_pics_sps; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->lt_ref_pic_poc_lsb_sps[ i ] = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ); printf("sps->lt_ref_pic_poc_lsb_sps[ i ]: %d \n", sps->lt_ref_pic_poc_lsb_sps[ i ]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->used_by_curr_pic_lt_sps_flag[ i ] = bs_read_u1(b); printf("sps->used_by_curr_pic_lt_sps_flag[ i ]: %d \n", sps->used_by_curr_pic_lt_sps_flag[ i ]); 
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_temporal_mvp_enabled_flag = bs_read_u1(b); printf("sps->sps_temporal_mvp_enabled_flag: %d \n", sps->sps_temporal_mvp_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->strong_intra_smoothing_enabled_flag = bs_read_u1(b); printf("sps->strong_intra_smoothing_enabled_flag: %d \n", sps->strong_intra_smoothing_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->vui_parameters_present_flag = bs_read_u1(b); printf("sps->vui_parameters_present_flag: %d \n", sps->vui_parameters_present_flag); 
    if( sps->vui_parameters_present_flag ) {
        read_debug_hevc_vui_parameters(sps, b);
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_extension_present_flag = bs_read_u1(b); printf("sps->sps_extension_present_flag: %d \n", sps->sps_extension_present_flag); 
    
    if( sps->sps_extension_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_range_extension_flag = bs_read_u1(b); printf("sps->sps_range_extension_flag: %d \n", sps->sps_range_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_multilayer_extension_flag = bs_read_u1(b); printf("sps->sps_multilayer_extension_flag: %d \n", sps->sps_multilayer_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_3d_extension_flag = bs_read_u1(b); printf("sps->sps_3d_extension_flag: %d \n", sps->sps_3d_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sps->sps_extension_5bits = bs_read_u(b, 5); printf("sps->sps_extension_5bits: %d \n", sps->sps_extension_5bits); 
    }
    if( sps->sps_range_extension_flag ) {
        read_debug_hevc_sps_range_extension( &sps->sps_range_ext, b);
    }
    
    if( 1 )
    {
        memcpy(h->sps_table[sps->sps_seq_parameter_set_id], h->sps, sizeof(hevc_sps_t));
    }
}

//7.3.2.2.2 Sequence parameter set range extension syntax
void read_debug_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b)
{
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->transform_skip_rotation_enabled_flag = bs_read_u1(b); printf("sps_range_ext->transform_skip_rotation_enabled_flag: %d \n", sps_range_ext->transform_skip_rotation_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->transform_skip_context_enabled_flag = bs_read_u1(b); printf("sps_range_ext->transform_skip_context_enabled_flag: %d \n", sps_range_ext->transform_skip_context_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->implicit_rdpcm_enabled_flag = bs_read_u1(b); printf("sps_range_ext->implicit_rdpcm_enabled_flag: %d \n", sps_range_ext->implicit_rdpcm_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->explicit_rdpcm_enabled_flag = bs_read_u1(b); printf("sps_range_ext->explicit_rdpcm_enabled_flag: %d \n", sps_range_ext->explicit_rdpcm_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->extended_precision_processing_flag = bs_read_u1(b); printf("sps_range_ext->extended_precision_processing_flag: %d \n", sps_range_ext->extended_precision_processing_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->intra_smoothing_disabled_flag = bs_read_u1(b); printf("sps_range_ext->intra_smoothing_disabled_flag: %d \n", sps_range_ext->intra_smoothing_disabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->high_precision_offsets_enabled_flag = bs_read_u1(b); printf("sps_range_ext->high_precision_offsets_enabled_flag: %d \n", sps_range_ext->high_precision_offsets_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->persistent_rice_adaptation_enabled_flag = bs_read_u1(b); printf("sps_range_ext->persistent_rice_adaptation_enabled_flag: %d \n", sps_range_ext->persistent_rice_adaptation_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sps_range_ext->cabac_bypass_alignment_enabled_flag = bs_read_u1(b); printf("sps_range_ext->cabac_bypass_alignment_enabled_flag: %d \n", sps_range_ext->cabac_bypass_alignment_enabled_flag); 
}


//7.3.2.3 Picture parameter set RBSP syntax
void read_debug_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_pps_t* pps = h->pps;
    if( 1 )
    {
        memset(pps, 0, sizeof(hevc_pps_t));
    }

    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pic_parameter_set_id = bs_read_ue(b); printf("pps->pic_parameter_set_id: %d \n", pps->pic_parameter_set_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->seq_parameter_set_id = bs_read_ue(b); printf("pps->seq_parameter_set_id: %d \n", pps->seq_parameter_set_id); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->dependent_slice_segments_enabled_flag = bs_read_u1(b); printf("pps->dependent_slice_segments_enabled_flag: %d \n", pps->dependent_slice_segments_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->output_flag_present_flag = bs_read_u1(b); printf("pps->output_flag_present_flag: %d \n", pps->output_flag_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->num_extra_slice_header_bits = bs_read_u(b,  3 ); printf("pps->num_extra_slice_header_bits: %d \n", pps->num_extra_slice_header_bits); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->sign_data_hiding_enabled_flag = bs_read_u1(b); printf("pps->sign_data_hiding_enabled_flag: %d \n", pps->sign_data_hiding_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->cabac_init_present_flag = bs_read_u1(b); printf("pps->cabac_init_present_flag: %d \n", pps->cabac_init_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->num_ref_idx_l0_default_active_minus1 = bs_read_ue(b); printf("pps->num_ref_idx_l0_default_active_minus1: %d \n", pps->num_ref_idx_l0_default_active_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->num_ref_idx_l1_default_active_minus1 = bs_read_ue(b); printf("pps->num_ref_idx_l1_default_active_minus1: %d \n", pps->num_ref_idx_l1_default_active_minus1); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->init_qp_minus26 = bs_read_se(b); printf("pps->init_qp_minus26: %d \n", pps->init_qp_minus26); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->constrained_intra_pred_flag = bs_read_u1(b); printf("pps->constrained_intra_pred_flag: %d \n", pps->constrained_intra_pred_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->transform_skip_enabled_flag = bs_read_u1(b); printf("pps->transform_skip_enabled_flag: %d \n", pps->transform_skip_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->cu_qp_delta_enabled_flag = bs_read_u1(b); printf("pps->cu_qp_delta_enabled_flag: %d \n", pps->cu_qp_delta_enabled_flag); 
    if( pps->cu_qp_delta_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->diff_cu_qp_delta_depth = bs_read_ue(b); printf("pps->diff_cu_qp_delta_depth: %d \n", pps->diff_cu_qp_delta_depth); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_cb_qp_offset = bs_read_se(b); printf("pps->pps_cb_qp_offset: %d \n", pps->pps_cb_qp_offset); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_cr_qp_offset = bs_read_se(b); printf("pps->pps_cr_qp_offset: %d \n", pps->pps_cr_qp_offset); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_slice_chroma_qp_offsets_present_flag = bs_read_u1(b); printf("pps->pps_slice_chroma_qp_offsets_present_flag: %d \n", pps->pps_slice_chroma_qp_offsets_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->weighted_pred_flag = bs_read_u1(b); printf("pps->weighted_pred_flag: %d \n", pps->weighted_pred_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->weighted_bipred_flag = bs_read_u1(b); printf("pps->weighted_bipred_flag: %d \n", pps->weighted_bipred_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->transquant_bypass_enabled_flag = bs_read_u1(b); printf("pps->transquant_bypass_enabled_flag: %d \n", pps->transquant_bypass_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->tiles_enabled_flag = bs_read_u1(b); printf("pps->tiles_enabled_flag: %d \n", pps->tiles_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->entropy_coding_sync_enabled_flag = bs_read_u1(b); printf("pps->entropy_coding_sync_enabled_flag: %d \n", pps->entropy_coding_sync_enabled_flag); 
    if( pps->tiles_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->num_tile_columns_minus1 = bs_read_ue(b); printf("pps->num_tile_columns_minus1: %d \n", pps->num_tile_columns_minus1); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->num_tile_rows_minus1 = bs_read_ue(b); printf("pps->num_tile_rows_minus1: %d \n", pps->num_tile_rows_minus1); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->uniform_spacing_flag = bs_read_u1(b); printf("pps->uniform_spacing_flag: %d \n", pps->uniform_spacing_flag); 
        if( !pps->uniform_spacing_flag ) {
            for( i = 0; i < pps->num_tile_columns_minus1; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->column_width_minus1[ i ] = bs_read_ue(b); printf("pps->column_width_minus1[ i ]: %d \n", pps->column_width_minus1[ i ]); 
            }
            for( i = 0; i < pps->num_tile_rows_minus1; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->row_height_minus1[ i ] = bs_read_ue(b); printf("pps->row_height_minus1[ i ]: %d \n", pps->row_height_minus1[ i ]); 
            }
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->loop_filter_across_tiles_enabled_flag = bs_read_u1(b); printf("pps->loop_filter_across_tiles_enabled_flag: %d \n", pps->loop_filter_across_tiles_enabled_flag); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_loop_filter_across_slices_enabled_flag = bs_read_u1(b); printf("pps->pps_loop_filter_across_slices_enabled_flag: %d \n", pps->pps_loop_filter_across_slices_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->deblocking_filter_control_present_flag = bs_read_u1(b); printf("pps->deblocking_filter_control_present_flag: %d \n", pps->deblocking_filter_control_present_flag); 
    if( pps->deblocking_filter_control_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->deblocking_filter_override_enabled_flag = bs_read_u1(b); printf("pps->deblocking_filter_override_enabled_flag: %d \n", pps->deblocking_filter_override_enabled_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_deblocking_filter_disabled_flag = bs_read_u1(b); printf("pps->pps_deblocking_filter_disabled_flag: %d \n", pps->pps_deblocking_filter_disabled_flag); 
        if( pps->pps_deblocking_filter_disabled_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_beta_offset_div2 = bs_read_se(b); printf("pps->pps_beta_offset_div2: %d \n", pps->pps_beta_offset_div2); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_tc_offset_div2 = bs_read_se(b); printf("pps->pps_tc_offset_div2: %d \n", pps->pps_tc_offset_div2); 
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_scaling_list_data_present_flag = bs_read_u1(b); printf("pps->pps_scaling_list_data_present_flag: %d \n", pps->pps_scaling_list_data_present_flag); 
    if( pps->pps_scaling_list_data_present_flag ) {
        read_debug_hevc_scaling_list_data(&pps->scaling_list_data, b);
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->lists_modification_present_flag = bs_read_u1(b); printf("pps->lists_modification_present_flag: %d \n", pps->lists_modification_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->log2_parallel_merge_level_minus2 = bs_read_ue(b); printf("pps->log2_parallel_merge_level_minus2: %d \n", pps->log2_parallel_merge_level_minus2); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->slice_segment_header_extension_present_flag = bs_read_u1(b); printf("pps->slice_segment_header_extension_present_flag: %d \n", pps->slice_segment_header_extension_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_extension_present_flag = bs_read_u1(b); printf("pps->pps_extension_present_flag: %d \n", pps->pps_extension_present_flag); 
    if( pps->pps_extension_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_range_extension_flag = bs_read_u1(b); printf("pps->pps_range_extension_flag: %d \n", pps->pps_range_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_multilayer_extension_flag = bs_read_u1(b); printf("pps->pps_multilayer_extension_flag: %d \n", pps->pps_multilayer_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_3d_extension_flag = bs_read_u1(b); printf("pps->pps_3d_extension_flag: %d \n", pps->pps_3d_extension_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps->pps_extension_5bits = bs_read_u1(b); printf("pps->pps_extension_5bits: %d \n", pps->pps_extension_5bits); 
    }
    if( pps->pps_range_extension_flag ) {
        read_debug_hevc_pps_range_extension( pps, b);
    }

    read_debug_hevc_rbsp_trailing_bits(b);

    if( 1 )
    {
        memcpy(h->pps_table[pps->pic_parameter_set_id], h->pps, sizeof(hevc_pps_t));
    }
}

//7.3.2.3.2 Picture parameter set range extension syntax
void read_debug_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b)
{
    hevc_pps_range_ext_t *pps_range_ext = &pps->pps_range_ext;;
    if( pps->transform_skip_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->log2_max_transform_skip_block_size_minus2 = bs_read_ue(b); printf("pps_range_ext->log2_max_transform_skip_block_size_minus2: %d \n", pps_range_ext->log2_max_transform_skip_block_size_minus2); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->cross_component_prediction_enabled_flag = bs_read_u1(b); printf("pps_range_ext->cross_component_prediction_enabled_flag: %d \n", pps_range_ext->cross_component_prediction_enabled_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->chroma_qp_offset_list_enabled_flag = bs_read_u1(b); printf("pps_range_ext->chroma_qp_offset_list_enabled_flag: %d \n", pps_range_ext->chroma_qp_offset_list_enabled_flag); 
    if( pps_range_ext->chroma_qp_offset_list_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->diff_cu_chroma_qp_offset_depth = bs_read_ue(b); printf("pps_range_ext->diff_cu_chroma_qp_offset_depth: %d \n", pps_range_ext->diff_cu_chroma_qp_offset_depth); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->chroma_qp_offset_list_len_minus1 = bs_read_ue(b); printf("pps_range_ext->chroma_qp_offset_list_len_minus1: %d \n", pps_range_ext->chroma_qp_offset_list_len_minus1); 
        for( int i = 0; i <= pps_range_ext->chroma_qp_offset_list_len_minus1; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->cb_qp_offset_list[ i ] = bs_read_se(b); printf("pps_range_ext->cb_qp_offset_list[ i ]: %d \n", pps_range_ext->cb_qp_offset_list[ i ]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->cr_qp_offset_list[ i ] = bs_read_se(b); printf("pps_range_ext->cr_qp_offset_list[ i ]: %d \n", pps_range_ext->cr_qp_offset_list[ i ]); 
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->log2_sao_offset_scale_luma = bs_read_ue(b); printf("pps_range_ext->log2_sao_offset_scale_luma: %d \n", pps_range_ext->log2_sao_offset_scale_luma); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); pps_range_ext->log2_sao_offset_scale_chroma = bs_read_ue(b); printf("pps_range_ext->log2_sao_offset_scale_chroma: %d \n", pps_range_ext->log2_sao_offset_scale_chroma); 
}

#ifdef HAVE_SEI
//7.3.2.4 Supplemental enhancement information RBSP syntax
void read_debug_sei_rbsp(hevc_stream_t* h, bs_t* b)
{
    if( 1 )
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
            read_debug_sei_message(h, b);
        } while( more_rbsp_data(h, b) );
    }

    if( 0 )
    {
        for (int i = 0; i < h->num_seis; i++) {
            h->sei = h->seis[i];
            read_debug_sei_message(h, b);
        }
        h->sei = NULL;
    }

    read_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.5 Supplemental enhancement information message syntax
void read_debug_sei_message(hevc_stream_t* h, bs_t* b)
{
    if( 0 )
    {
        _write_ff_coded_number(b, h->sei->payloadType);
        _write_ff_coded_number(b, h->sei->payloadSize);
    }
    if( 1 )
    {
        h->sei->payloadType = _read_ff_coded_number(b);
        h->sei->payloadSize = _read_ff_coded_number(b);
    }
    read_debug_sei_payload( h, b, h->sei->payloadType, h->sei->payloadSize );
}
#endif

//7.3.2.5 Access unit delimiter RBSP syntax
void read_debug_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b)
{
    printf("%ld.%d: ", b->p - b->start, b->bits_left); h->aud->primary_pic_type = bs_read_u(b, 3); printf("h->aud->primary_pic_type: %d \n", h->aud->primary_pic_type); 
    read_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.6 End of sequence RBSP syntax
void read_debug_hevc_end_of_seq_rbsp()
{
}

//7.3.2.7 End of bitstream RBSP syntax
void read_debug_end_of_bitstream_rbsp()
{
}

//7.3.2.8 Filler data RBSP syntax
void read_debug_filler_data_rbsp(bs_t* b)
{
    while( bs_next_bits(b, 8) == 0xFF )
    {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); int ff_byte = bs_read_u(b, 8); printf("ff_byte: %d \n", ff_byte); 
    }
    read_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.9 Slice segment layer RBSP syntax
void read_debug_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b)
{
    read_debug_hevc_slice_header(h, b);
    hevc_slice_data_rbsp_t* slice_data = h->slice_data;

    if ( slice_data != NULL )
    {
        if ( slice_data->rbsp_buf != NULL ) free( slice_data->rbsp_buf ); 
        uint8_t *sptr = b->p + (!!b->bits_left); // CABAC-specific: skip alignment bits, if there are any
        slice_data->rbsp_size = b->end - sptr;
        
        slice_data->rbsp_buf = (uint8_t*)malloc(slice_data->rbsp_size);
        memcpy( slice_data->rbsp_buf, sptr, slice_data->rbsp_size );
    }

    //read_debug_hevc_slice_data(h, b); /* all categories of slice_data( ) syntax */
    read_debug_hevc_rbsp_slice_trailing_bits( b );
}

//7.3.2.10 RBSP slice trailing bits syntax
void read_debug_hevc_rbsp_slice_trailing_bits(bs_t* b)
{
    read_debug_hevc_rbsp_trailing_bits(b);
    //while( more_rbsp_trailing_data(b) )
    //{
    //    value( cabac_zero_word, f(16, 0x0000) );
    //}
}

//7.3.2.11 RBSP trailing bits syntax
void read_debug_hevc_rbsp_trailing_bits(bs_t* b)
{
    printf("%ld.%d: ", b->p - b->start, b->bits_left); int rbsp_stop_one_bit = bs_read_u(b, 1); printf("rbsp_stop_one_bit: %d \n", rbsp_stop_one_bit); 

    while( !bs_byte_aligned(b) )
    {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); int rbsp_alignment_zero_bit = bs_read_u(b, 1); printf("rbsp_alignment_zero_bit: %d \n", rbsp_alignment_zero_bit); 
    }
}

//7.3.2.12 Byte alignment syntax
void read_debug_hevc_byte_alignment(bs_t* b)
{
    printf("%ld.%d: ", b->p - b->start, b->bits_left); int alignment_bit_equal_to_one = bs_read_u(b, 1); printf("alignment_bit_equal_to_one: %d \n", alignment_bit_equal_to_one); 

    while( !bs_byte_aligned(b) )
    {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); int alignment_bit_equal_to_zero = bs_read_u(b, 1); printf("alignment_bit_equal_to_zero: %d \n", alignment_bit_equal_to_zero); 
    }
}

//7.3.3 Profile, tier and level syntax
void read_debug_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1)
{
    int i, j;
    if( profilePresentFlag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_profile_space = bs_read_u(b, 2); printf("ptl->general_profile_space: %d \n", ptl->general_profile_space); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_tier_flag = bs_read_u1(b); printf("ptl->general_tier_flag: %d \n", ptl->general_tier_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_profile_idc = bs_read_u(b, 5); printf("ptl->general_profile_idc: %d \n", ptl->general_profile_idc); 
        for( i = 0; i < 32; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_profile_compatibility_flag[ i ] = bs_read_u1(b); printf("ptl->general_profile_compatibility_flag[ i ]: %d \n", ptl->general_profile_compatibility_flag[ i ]); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_progressive_source_flag = bs_read_u1(b); printf("ptl->general_progressive_source_flag: %d \n", ptl->general_progressive_source_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_interlaced_source_flag = bs_read_u1(b); printf("ptl->general_interlaced_source_flag: %d \n", ptl->general_interlaced_source_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_non_packed_constraint_flag = bs_read_u1(b); printf("ptl->general_non_packed_constraint_flag: %d \n", ptl->general_non_packed_constraint_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_frame_only_constraint_flag = bs_read_u1(b); printf("ptl->general_frame_only_constraint_flag: %d \n", ptl->general_frame_only_constraint_flag); 
        if( ptl->general_profile_idc == 4 || ptl->general_profile_compatibility_flag[ 4 ] || 
            ptl->general_profile_idc == 5 || ptl->general_profile_compatibility_flag[ 5 ] || 
            ptl->general_profile_idc == 6 || ptl->general_profile_compatibility_flag[ 6 ] || 
            ptl->general_profile_idc == 7 || ptl->general_profile_compatibility_flag[ 7 ] ) {
                
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_12bit_constraint_flag = bs_read_u1(b); printf("ptl->general_max_12bit_constraint_flag: %d \n", ptl->general_max_12bit_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_10bit_constraint_flag = bs_read_u1(b); printf("ptl->general_max_10bit_constraint_flag: %d \n", ptl->general_max_10bit_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_8bit_constraint_flag = bs_read_u1(b); printf("ptl->general_max_8bit_constraint_flag: %d \n", ptl->general_max_8bit_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_422chroma_constraint_flag = bs_read_u1(b); printf("ptl->general_max_422chroma_constraint_flag: %d \n", ptl->general_max_422chroma_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_420chroma_constraint_flag = bs_read_u1(b); printf("ptl->general_max_420chroma_constraint_flag: %d \n", ptl->general_max_420chroma_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_max_monochrome_constraint_flag = bs_read_u1(b); printf("ptl->general_max_monochrome_constraint_flag: %d \n", ptl->general_max_monochrome_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_intra_constraint_flag = bs_read_u1(b); printf("ptl->general_intra_constraint_flag: %d \n", ptl->general_intra_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_one_picture_only_constraint_flag = bs_read_u1(b); printf("ptl->general_one_picture_only_constraint_flag: %d \n", ptl->general_one_picture_only_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_lower_bit_rate_constraint_flag = bs_read_u1(b); printf("ptl->general_lower_bit_rate_constraint_flag: %d \n", ptl->general_lower_bit_rate_constraint_flag); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); int general_reserved_zero_34bits = bs_read_u(b, 34); printf("general_reserved_zero_34bits: %d \n", general_reserved_zero_34bits); 
        } else {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); int general_reserved_zero_43bits = bs_read_u(b, 43); printf("general_reserved_zero_43bits: %d \n", general_reserved_zero_43bits); 
        }
        if( ( ptl->general_profile_idc >= 1 && ptl->general_profile_idc <= 5 ) ||
              ptl->general_profile_compatibility_flag[ 1 ] ||
              ptl->general_profile_compatibility_flag[ 2 ] ||
              ptl->general_profile_compatibility_flag[ 3 ] ||
              ptl->general_profile_compatibility_flag[ 4 ] ||
              ptl->general_profile_compatibility_flag[ 5 ] ) {

            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_inbld_flag = bs_read_u1(b); printf("ptl->general_inbld_flag: %d \n", ptl->general_inbld_flag); 
        } else {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); int general_reserved_zero_bit = bs_read_u(b, 1); printf("general_reserved_zero_bit: %d \n", general_reserved_zero_bit); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->general_level_idc = bs_read_u8(b); printf("ptl->general_level_idc: %d \n", ptl->general_level_idc); 
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_profile_present_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_profile_present_flag[ i ]: %d \n", ptl->sub_layer_profile_present_flag[ i ]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_level_present_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_level_present_flag[ i ]: %d \n", ptl->sub_layer_level_present_flag[ i ]); 
        }
        if( maxNumSubLayersMinus1 > 0 ) {
            for( i = maxNumSubLayersMinus1; i < 8; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); int reserved_zero_xxbits = bs_read_u(b, 2); printf("reserved_zero_xxbits: %d \n", reserved_zero_xxbits); 
            }
        }
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) { 
            if( ptl->sub_layer_profile_present_flag[ i ] ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_profile_space[ i ] = bs_read_u(b, 2); printf("ptl->sub_layer_profile_space[ i ]: %d \n", ptl->sub_layer_profile_space[ i ]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_tier_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_tier_flag[ i ]: %d \n", ptl->sub_layer_tier_flag[ i ]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_profile_idc[ i ] = bs_read_u(b, 5); printf("ptl->sub_layer_profile_idc[ i ]: %d \n", ptl->sub_layer_profile_idc[ i ]); 
                for( j = 0; j < 32; j++ ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_profile_compatibility_flag[ i ][ j ] = bs_read_u(b, 1); printf("ptl->sub_layer_profile_compatibility_flag[ i ][ j ]: %d \n", ptl->sub_layer_profile_compatibility_flag[ i ][ j ]); 
                }
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_progressive_source_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_progressive_source_flag[ i ]: %d \n", ptl->sub_layer_progressive_source_flag[ i ]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_interlaced_source_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_interlaced_source_flag[ i ]: %d \n", ptl->sub_layer_interlaced_source_flag[ i ]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_non_packed_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_non_packed_constraint_flag[ i ]: %d \n", ptl->sub_layer_non_packed_constraint_flag[ i ]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_frame_only_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_frame_only_constraint_flag[ i ]: %d \n", ptl->sub_layer_frame_only_constraint_flag[ i ]); 
                if( ptl->sub_layer_profile_idc[ i ] == 4 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 4 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 5 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 5 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 6 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 6 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 7 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 7 ] ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_12bit_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_12bit_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_12bit_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_10bit_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_10bit_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_10bit_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_8bit_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_8bit_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_8bit_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_422chroma_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_422chroma_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_422chroma_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_420chroma_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_420chroma_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_420chroma_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_max_monochrome_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_max_monochrome_constraint_flag[ i ]: %d \n", ptl->sub_layer_max_monochrome_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_intra_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_intra_constraint_flag[ i ]: %d \n", ptl->sub_layer_intra_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_one_picture_only_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_one_picture_only_constraint_flag[ i ]: %d \n", ptl->sub_layer_one_picture_only_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_lower_bit_rate_constraint_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_lower_bit_rate_constraint_flag[ i ]: %d \n", ptl->sub_layer_lower_bit_rate_constraint_flag[ i ]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); int sub_layer_reserved_zero_34bits = bs_read_u(b, 34); printf("sub_layer_reserved_zero_34bits: %d \n", sub_layer_reserved_zero_34bits); 
                } else {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); int sub_layer_reserved_zero_43bits = bs_read_u(b, 43); printf("sub_layer_reserved_zero_43bits: %d \n", sub_layer_reserved_zero_43bits); 
                }
            
                if( ( ptl->sub_layer_profile_idc[ i ] >= 1 && ptl->sub_layer_profile_idc[ i ] <= 5 ) ||
                   ptl->sub_layer_profile_compatibility_flag[ 1 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 2 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 3 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 4 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 5 ] ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_inbld_flag[ i ] = bs_read_u1(b); printf("ptl->sub_layer_inbld_flag[ i ]: %d \n", ptl->sub_layer_inbld_flag[ i ]); 
                } else {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); int sub_layer_reserved_zero_bit = bs_read_u(b, 1); printf("sub_layer_reserved_zero_bit: %d \n", sub_layer_reserved_zero_bit); 
                }
            }
            if( ptl->sub_layer_level_present_flag[ i ] ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); ptl->sub_layer_level_idc[ i ] = bs_read_u1(b); printf("ptl->sub_layer_level_idc[ i ]: %d \n", ptl->sub_layer_level_idc[ i ]); 
            }
        }
    }
}

//7.3.4 Scaling list data syntax
void read_debug_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b )
{
    int nextCoef, coefNum;
    for( int sizeId = 0; sizeId < 4; sizeId++ )
        for( int matrixId = 0; matrixId < 6; matrixId += ( sizeId == 3 ) ? 3 : 1 ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] = bs_read_u1(b); printf("sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ]: %d \n", sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ]); 
            if( !sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ] = bs_read_ue(b); printf("sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ]: %d \n", sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ]); 
            } else {
                nextCoef = 8;
                coefNum=MIN(64, (1 << (4+(sizeId << 1))));
                if( sizeId > 1 ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ] = bs_read_se(b); printf("sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ]: %d \n", sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ]); 
                }
 
                for( int i = 0; i < coefNum; i++) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sld->scaling_list_delta_coef[ sizeId ][ matrixId ] = bs_read_se(b); printf("sld->scaling_list_delta_coef[ sizeId ][ matrixId ]: %d \n", sld->scaling_list_delta_coef[ sizeId ][ matrixId ]); 
                }
            }
 
        }
}

//7.3.6 Slice header syntax
void read_debug_hevc_slice_header(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    if( 1 )
    {
        init_slice_hevc(h);
    }

    hevc_nal_t* nal = h->nal;

    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->first_slice_segment_in_pic_flag = bs_read_u1(b); printf("sh->first_slice_segment_in_pic_flag: %d \n", sh->first_slice_segment_in_pic_flag); 
    if( nal->nal_unit_type >= HEVC_NAL_UNIT_TYPE_BLA_W_LP && 
       nal->nal_unit_type <= HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->no_output_of_prior_pics_flag = bs_read_u1(b); printf("sh->no_output_of_prior_pics_flag: %d \n", sh->no_output_of_prior_pics_flag); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->pic_parameter_set_id = bs_read_ue(b); printf("sh->pic_parameter_set_id: %d \n", sh->pic_parameter_set_id); 
    
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];

    //set default value
    sh->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    sh->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;

    if( !sh->first_slice_segment_in_pic_flag ) {
        if( pps->dependent_slice_segments_enabled_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->dependent_slice_segment_flag = bs_read_u1(b); printf("sh->dependent_slice_segment_flag: %d \n", sh->dependent_slice_segment_flag); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_segment_address = bs_read_u(b,  getSliceSegmentAddressBitLength( sps ) ); printf("sh->slice_segment_address: %d \n", sh->slice_segment_address); 
    }
    
    if( !sh->dependent_slice_segment_flag ) {
        for( i = 0; i < pps->num_extra_slice_header_bits; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); int slice_reserved_flag = bs_read_u(b, 1); printf("slice_reserved_flag: %d \n", slice_reserved_flag); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_type = bs_read_ue(b); printf("sh->slice_type: %d \n", sh->slice_type); 
        if( pps->output_flag_present_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->pic_output_flag = bs_read_u1(b); printf("sh->pic_output_flag: %d \n", sh->pic_output_flag); 
        }
        if( sps->separate_colour_plane_flag == 1 ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->colour_plane_id = bs_read_u(b, 2); printf("sh->colour_plane_id: %d \n", sh->colour_plane_id); 
        }
        if( nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_W_RADL &&
            nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_N_LP) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_pic_order_cnt_lsb = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ); printf("sh->slice_pic_order_cnt_lsb: %d \n", sh->slice_pic_order_cnt_lsb); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->short_term_ref_pic_set_sps_flag = bs_read_u1(b); printf("sh->short_term_ref_pic_set_sps_flag: %d \n", sh->short_term_ref_pic_set_sps_flag); 
            if( !sh->short_term_ref_pic_set_sps_flag ) {
                read_debug_hevc_st_ref_pic_set( &sh->st_ref_pic_set, b, sps->num_short_term_ref_pic_sets, sps->num_short_term_ref_pic_sets );
            } else if( sps->num_short_term_ref_pic_sets > 1 ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->short_term_ref_pic_set_idx = bs_read_u(b,  ceil( log2( sps->num_short_term_ref_pic_sets ) ) ); printf("sh->short_term_ref_pic_set_idx: %d \n", sh->short_term_ref_pic_set_idx); 
            }
            if( sps->long_term_ref_pics_present_flag ) {
                if( sps->num_long_term_ref_pics_sps > 0 ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_long_term_sps = bs_read_ue(b); printf("sh->num_long_term_sps: %d \n", sh->num_long_term_sps); 
                }
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_long_term_pics = bs_read_ue(b); printf("sh->num_long_term_pics: %d \n", sh->num_long_term_pics); 
                for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
                    if( i < sh->num_long_term_sps ) {
                        if( sps->num_long_term_ref_pics_sps > 1 ) {
                            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->lt_idx_sps[ i ] = bs_read_u(b,  ceil( log2( sps->num_long_term_ref_pics_sps ) ) ); printf("sh->lt_idx_sps[ i ]: %d \n", sh->lt_idx_sps[ i ]); 
                        }
                    } else {
                        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->poc_lsb_lt[ i ] = bs_read_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 ); printf("sh->poc_lsb_lt[ i ]: %d \n", sh->poc_lsb_lt[ i ]); 
                        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->used_by_curr_pic_lt_flag[ i ] = bs_read_u1(b); printf("sh->used_by_curr_pic_lt_flag[ i ]: %d \n", sh->used_by_curr_pic_lt_flag[ i ]); 
                    }
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->delta_poc_msb_present_flag[ i ] = bs_read_u1(b); printf("sh->delta_poc_msb_present_flag[ i ]: %d \n", sh->delta_poc_msb_present_flag[ i ]); 
                    if( sh->delta_poc_msb_present_flag[ i ]) {
                        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->delta_poc_msb_cycle_lt[ i ] = bs_read_ue(b); printf("sh->delta_poc_msb_cycle_lt[ i ]: %d \n", sh->delta_poc_msb_cycle_lt[ i ]); 
                    }
                }
            }
            if( sps->sps_temporal_mvp_enabled_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_temporal_mvp_enabled_flag = bs_read_u1(b); printf("sh->slice_temporal_mvp_enabled_flag: %d \n", sh->slice_temporal_mvp_enabled_flag); 
            }
        }
        if( sps->sample_adaptive_offset_enabled_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_sao_luma_flag = bs_read_u1(b); printf("sh->slice_sao_luma_flag: %d \n", sh->slice_sao_luma_flag); 
            int ChromaArrayType = 0;
            if( sps->separate_colour_plane_flag == 0) {
                ChromaArrayType = sps->chroma_format_idc;
            }
            if( ChromaArrayType != 0 ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_sao_chroma_flag = bs_read_u1(b); printf("sh->slice_sao_chroma_flag: %d \n", sh->slice_sao_chroma_flag); 
            }
        }
        if( sh->slice_type == HEVC_SLICE_TYPE_P || sh->slice_type == HEVC_SLICE_TYPE_B ){
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_ref_idx_active_override_flag = bs_read_u1(b); printf("sh->num_ref_idx_active_override_flag: %d \n", sh->num_ref_idx_active_override_flag); 
            if( sh->num_ref_idx_active_override_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_ref_idx_l0_active_minus1 = bs_read_ue(b); printf("sh->num_ref_idx_l0_active_minus1: %d \n", sh->num_ref_idx_l0_active_minus1); 
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_ref_idx_l1_active_minus1 = bs_read_ue(b); printf("sh->num_ref_idx_l1_active_minus1: %d \n", sh->num_ref_idx_l1_active_minus1); 
                }
            }
            if( pps->lists_modification_present_flag && getNumPicTotalCurr( sps, sh ) > 1 ) {
                read_debug_hevc_ref_pic_lists_modification( h, b );
            }
            
            if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->mvd_l1_zero_flag = bs_read_u1(b); printf("sh->mvd_l1_zero_flag: %d \n", sh->mvd_l1_zero_flag); 
            }
            if( pps->cabac_init_present_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->cabac_init_flag = bs_read_u1(b); printf("sh->cabac_init_flag: %d \n", sh->cabac_init_flag); 
            }
            if( sh->slice_temporal_mvp_enabled_flag ) {
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->collocated_from_l0_flag = bs_read_u1(b); printf("sh->collocated_from_l0_flag: %d \n", sh->collocated_from_l0_flag); 
                }
                if( ( sh->collocated_from_l0_flag && sh->num_ref_idx_l0_active_minus1 > 0 ) ||
                   ( !sh->collocated_from_l0_flag && sh->num_ref_idx_l1_active_minus1 > 0 ) ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->collocated_ref_idx = bs_read_ue(b); printf("sh->collocated_ref_idx: %d \n", sh->collocated_ref_idx); 
                }
            }
            if( ( pps->weighted_pred_flag && sh->slice_type == HEVC_SLICE_TYPE_P ) || 
                ( pps->weighted_bipred_flag && sh->slice_type == HEVC_SLICE_TYPE_B ) ) {
                read_debug_hevc_pred_weight_table( h, b );
            }
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->five_minus_max_num_merge_cand = bs_read_ue(b); printf("sh->five_minus_max_num_merge_cand: %d \n", sh->five_minus_max_num_merge_cand); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_qp_delta = bs_read_se(b); printf("sh->slice_qp_delta: %d \n", sh->slice_qp_delta); 
        if( pps->pps_slice_chroma_qp_offsets_present_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_cb_qp_offset = bs_read_se(b); printf("sh->slice_cb_qp_offset: %d \n", sh->slice_cb_qp_offset); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_cr_qp_offset = bs_read_se(b); printf("sh->slice_cr_qp_offset: %d \n", sh->slice_cr_qp_offset); 
        }
        if( pps->pps_range_ext.chroma_qp_offset_list_enabled_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->cu_chroma_qp_offset_enabled_flag = bs_read_u1(b); printf("sh->cu_chroma_qp_offset_enabled_flag: %d \n", sh->cu_chroma_qp_offset_enabled_flag); 
        }
        if( pps->deblocking_filter_override_enabled_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->deblocking_filter_override_flag = bs_read_u1(b); printf("sh->deblocking_filter_override_flag: %d \n", sh->deblocking_filter_override_flag); 
        }
        if( sh->deblocking_filter_override_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_deblocking_filter_disabled_flag = bs_read_u1(b); printf("sh->slice_deblocking_filter_disabled_flag: %d \n", sh->slice_deblocking_filter_disabled_flag); 
            if( !sh->slice_deblocking_filter_disabled_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_beta_offset_div2 = bs_read_se(b); printf("sh->slice_beta_offset_div2: %d \n", sh->slice_beta_offset_div2); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_tc_offset_div2 = bs_read_se(b); printf("sh->slice_tc_offset_div2: %d \n", sh->slice_tc_offset_div2); 
            }
        }
        if( pps->pps_loop_filter_across_slices_enabled_flag &&
           ( sh->slice_sao_luma_flag || sh->slice_sao_chroma_flag || !sh->slice_deblocking_filter_disabled_flag ) ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_loop_filter_across_slices_enabled_flag = bs_read_u1(b); printf("sh->slice_loop_filter_across_slices_enabled_flag: %d \n", sh->slice_loop_filter_across_slices_enabled_flag); 
        }
    }
    if( pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->num_entry_point_offsets = bs_read_ue(b); printf("sh->num_entry_point_offsets: %d \n", sh->num_entry_point_offsets); 
        if( sh->num_entry_point_offsets > 0 ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->offset_len_minus1 = bs_read_ue(b); printf("sh->offset_len_minus1: %d \n", sh->offset_len_minus1); 
            for( i = 0; i < sh->num_entry_point_offsets; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->entry_point_offset_minus1[ i ] = bs_read_u(b,  sh->offset_len_minus1 + 1 ); printf("sh->entry_point_offset_minus1[ i ]: %d \n", sh->entry_point_offset_minus1[ i ]); 
            }
        }
    }
    if( pps->slice_segment_header_extension_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->slice_segment_header_extension_length = bs_read_ue(b); printf("sh->slice_segment_header_extension_length: %d \n", sh->slice_segment_header_extension_length); 
        //TODO: support header extension,
        for( i = 0; i < sh->slice_segment_header_extension_length; i++) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); int slice_segment_header_extension_data_byte = bs_read_u(b, 8); printf("slice_segment_header_extension_data_byte: %d \n", slice_segment_header_extension_data_byte); 
        }
    }
    read_debug_hevc_byte_alignment( b );
}

//7.3.6.2 Reference picture list reordering syntax
void read_debug_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    
    printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->rpld.ref_pic_list_modification_flag_l0 = bs_read_u1(b); printf("sh->rpld.ref_pic_list_modification_flag_l0: %d \n", sh->rpld.ref_pic_list_modification_flag_l0); 
    if( sh->rpld.ref_pic_list_modification_flag_l0 ) {
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->rpld.list_entry_l0[ i ] = bs_read_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) ); printf("sh->rpld.list_entry_l0[ i ]: %d \n", sh->rpld.list_entry_l0[ i ]); 
        }
    }
    
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); // ERROR: value( sh->rpld.ref_pic_list_modification_flag_l1, 1 ); printf("sh->rpld.ref_pic_list_modification_flag_l1: %d \n", sh->rpld.ref_pic_list_modification_flag_l1); 
        if( sh->rpld.ref_pic_list_modification_flag_l1 ) {
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); sh->rpld.list_entry_l1[ i ] = bs_read_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) ); printf("sh->rpld.list_entry_l1[ i ]: %d \n", sh->rpld.list_entry_l1[ i ]); 
            }
        }
    }
}

//7.3.6.3 Prediction weight table syntax
void read_debug_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b)
{
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    hevc_pred_weight_table_t *pwt = &sh->pwt;

    int i, j;

    printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->luma_log2_weight_denom = bs_read_ue(b); printf("pwt->luma_log2_weight_denom: %d \n", pwt->luma_log2_weight_denom); 
    
    int ChromaArrayType = 0;
    if( sps->separate_colour_plane_flag == 0) {
        ChromaArrayType = sps->chroma_format_idc;
    }

    if( ChromaArrayType != 0 ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_chroma_log2_weight_denom = bs_read_se(b); printf("pwt->delta_chroma_log2_weight_denom: %d \n", pwt->delta_chroma_log2_weight_denom); 
    }
    
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->luma_weight_l0_flag[i] = bs_read_u1(b); printf("pwt->luma_weight_l0_flag[i]: %d \n", pwt->luma_weight_l0_flag[i]); 
    }
    if( ChromaArrayType != 0 )
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->chroma_weight_l0_flag[i] = bs_read_u1(b); printf("pwt->chroma_weight_l0_flag[i]: %d \n", pwt->chroma_weight_l0_flag[i]); 
        }
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        if( pwt->luma_weight_l0_flag[i] ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_luma_weight_l0[i] = bs_read_se(b); printf("pwt->delta_luma_weight_l0[i]: %d \n", pwt->delta_luma_weight_l0[i]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->luma_offset_l0[i] = bs_read_se(b); printf("pwt->luma_offset_l0[i]: %d \n", pwt->luma_offset_l0[i]); 
        }
        if( pwt->chroma_weight_l0_flag[i] ) {
            for( j =0; j < 2; j++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_chroma_weight_l0[i][j] = bs_read_se(b); printf("pwt->delta_chroma_weight_l0[i][j]: %d \n", pwt->delta_chroma_weight_l0[i][j]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_chroma_offset_l0[i][j] = bs_read_se(b); printf("pwt->delta_chroma_offset_l0[i][j]: %d \n", pwt->delta_chroma_offset_l0[i][j]); 
            }
        }
    }
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->luma_weight_l1_flag[i] = bs_read_u1(b); printf("pwt->luma_weight_l1_flag[i]: %d \n", pwt->luma_weight_l1_flag[i]); 
        }
        if( ChromaArrayType != 0 )
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->chroma_weight_l1_flag[i] = bs_read_u1(b); printf("pwt->chroma_weight_l1_flag[i]: %d \n", pwt->chroma_weight_l1_flag[i]); 
            }
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            if( pwt->luma_weight_l1_flag[i] ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_luma_weight_l1[i] = bs_read_se(b); printf("pwt->delta_luma_weight_l1[i]: %d \n", pwt->delta_luma_weight_l1[i]); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->luma_offset_l1[i] = bs_read_se(b); printf("pwt->luma_offset_l1[i]: %d \n", pwt->luma_offset_l1[i]); 
            }
            if( pwt->chroma_weight_l1_flag[i] ) {
                for( j =0; j < 2; j++ ) {
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_chroma_weight_l1[i][j] = bs_read_se(b); printf("pwt->delta_chroma_weight_l1[i][j]: %d \n", pwt->delta_chroma_weight_l1[i][j]); 
                    printf("%ld.%d: ", b->p - b->start, b->bits_left); pwt->delta_chroma_offset_l1[i][j] = bs_read_se(b); printf("pwt->delta_chroma_offset_l1[i][j]: %d \n", pwt->delta_chroma_offset_l1[i][j]); 
                }
            }
        }        
    }
}

//7.3.7 Short-term reference picture set syntax
void read_debug_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets )
{
    int i, j;
    
    if( stRpsIdx != 0 ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->inter_ref_pic_set_prediction_flag = bs_read_u1(b); printf("st_ref_pic_set->inter_ref_pic_set_prediction_flag: %d \n", st_ref_pic_set->inter_ref_pic_set_prediction_flag); 
    }
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        if( stRpsIdx == num_short_term_ref_pic_sets ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->delta_idx_minus1 = bs_read_ue(b); printf("st_ref_pic_set->delta_idx_minus1: %d \n", st_ref_pic_set->delta_idx_minus1); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->delta_rps_sign = bs_read_u1(b); printf("st_ref_pic_set->delta_rps_sign: %d \n", st_ref_pic_set->delta_rps_sign); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->abs_delta_rps_minus1 = bs_read_ue(b); printf("st_ref_pic_set->abs_delta_rps_minus1: %d \n", st_ref_pic_set->abs_delta_rps_minus1); 
        
        int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
        
        for( j = 0; j <= NumDeltaPocs[ RefRpsIdx ]; j++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->used_by_curr_pic_flag[ j ] = bs_read_u1(b); printf("st_ref_pic_set->used_by_curr_pic_flag[ j ]: %d \n", st_ref_pic_set->used_by_curr_pic_flag[ j ]); 
            if( !st_ref_pic_set->used_by_curr_pic_flag[ j ] ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->use_delta_flag[ j ] = bs_read_u1(b); printf("st_ref_pic_set->use_delta_flag[ j ]: %d \n", st_ref_pic_set->use_delta_flag[ j ]); 
            }
        }
    } else {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->num_negative_pics = bs_read_ue(b); printf("st_ref_pic_set->num_negative_pics: %d \n", st_ref_pic_set->num_negative_pics); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->num_positive_pics = bs_read_ue(b); printf("st_ref_pic_set->num_positive_pics: %d \n", st_ref_pic_set->num_positive_pics); 
        for( i = 0; i < st_ref_pic_set->num_negative_pics; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->delta_poc_s0_minus1[ i ] = bs_read_ue(b); printf("st_ref_pic_set->delta_poc_s0_minus1[ i ]: %d \n", st_ref_pic_set->delta_poc_s0_minus1[ i ]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->used_by_curr_pic_s0_flag[ i ] = bs_read_u1(b); printf("st_ref_pic_set->used_by_curr_pic_s0_flag[ i ]: %d \n", st_ref_pic_set->used_by_curr_pic_s0_flag[ i ]); 
            
            //update derived field
            UsedByCurrPicS0[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s0_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS0[ stRpsIdx ][ i ] = -1 * ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            } else {
                DeltaPocS0[ stRpsIdx ][ i ] = DeltaPocS0[ stRpsIdx ][ i - 1 ] - ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            }
        }
        for( i = 0; i < st_ref_pic_set->num_positive_pics; i++ ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->delta_poc_s1_minus1[ i ] = bs_read_ue(b); printf("st_ref_pic_set->delta_poc_s1_minus1[ i ]: %d \n", st_ref_pic_set->delta_poc_s1_minus1[ i ]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); st_ref_pic_set->used_by_curr_pic_s1_flag[ i ] = bs_read_u1(b); printf("st_ref_pic_set->used_by_curr_pic_s1_flag[ i ]: %d \n", st_ref_pic_set->used_by_curr_pic_s1_flag[ i ]); 
        
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
void read_debug_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b)
{
    hevc_vui_t* vui = &sps->vui;
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->aspect_ratio_info_present_flag = bs_read_u1(b); printf("vui->aspect_ratio_info_present_flag: %d \n", vui->aspect_ratio_info_present_flag); 
    if( vui->aspect_ratio_info_present_flag )
    {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->aspect_ratio_idc = bs_read_u8(b); printf("vui->aspect_ratio_idc: %d \n", vui->aspect_ratio_idc); 
        if( vui->aspect_ratio_idc == SAR_Extended )
        {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->sar_width = bs_read_u(b, 16); printf("vui->sar_width: %d \n", vui->sar_width); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->sar_height = bs_read_u(b, 16); printf("vui->sar_height: %d \n", vui->sar_height); 
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->overscan_info_present_flag = bs_read_u1(b); printf("vui->overscan_info_present_flag: %d \n", vui->overscan_info_present_flag); 
    if( vui->overscan_info_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->overscan_appropriate_flag = bs_read_u1(b); printf("vui->overscan_appropriate_flag: %d \n", vui->overscan_appropriate_flag); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->video_signal_type_present_flag = bs_read_u1(b); printf("vui->video_signal_type_present_flag: %d \n", vui->video_signal_type_present_flag); 
    if( vui->video_signal_type_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->video_format = bs_read_u(b, 3); printf("vui->video_format: %d \n", vui->video_format); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->video_full_range_flag = bs_read_u1(b); printf("vui->video_full_range_flag: %d \n", vui->video_full_range_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->colour_description_present_flag = bs_read_u1(b); printf("vui->colour_description_present_flag: %d \n", vui->colour_description_present_flag); 
        if( vui->colour_description_present_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->colour_primaries = bs_read_u8(b); printf("vui->colour_primaries: %d \n", vui->colour_primaries); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->transfer_characteristics = bs_read_u8(b); printf("vui->transfer_characteristics: %d \n", vui->transfer_characteristics); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->matrix_coefficients = bs_read_u8(b); printf("vui->matrix_coefficients: %d \n", vui->matrix_coefficients); 
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->chroma_loc_info_present_flag = bs_read_u1(b); printf("vui->chroma_loc_info_present_flag: %d \n", vui->chroma_loc_info_present_flag); 
    if( vui->chroma_loc_info_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->chroma_sample_loc_type_top_field = bs_read_ue(b); printf("vui->chroma_sample_loc_type_top_field: %d \n", vui->chroma_sample_loc_type_top_field); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->chroma_sample_loc_type_bottom_field = bs_read_ue(b); printf("vui->chroma_sample_loc_type_bottom_field: %d \n", vui->chroma_sample_loc_type_bottom_field); 
    }
    
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->neutral_chroma_indication_flag = bs_read_u1(b); printf("vui->neutral_chroma_indication_flag: %d \n", vui->neutral_chroma_indication_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->field_seq_flag = bs_read_u1(b); printf("vui->field_seq_flag: %d \n", vui->field_seq_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->frame_field_info_present_flag = bs_read_u1(b); printf("vui->frame_field_info_present_flag: %d \n", vui->frame_field_info_present_flag); 
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->default_display_window_flag = bs_read_u1(b); printf("vui->default_display_window_flag: %d \n", vui->default_display_window_flag); 
    if( vui->default_display_window_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->def_disp_win_left_offset = bs_read_ue(b); printf("vui->def_disp_win_left_offset: %d \n", vui->def_disp_win_left_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->def_disp_win_right_offset = bs_read_ue(b); printf("vui->def_disp_win_right_offset: %d \n", vui->def_disp_win_right_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->def_disp_win_top_offset = bs_read_ue(b); printf("vui->def_disp_win_top_offset: %d \n", vui->def_disp_win_top_offset); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->def_disp_win_bottom_offset = bs_read_ue(b); printf("vui->def_disp_win_bottom_offset: %d \n", vui->def_disp_win_bottom_offset); 
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_timing_info_present_flag = bs_read_u1(b); printf("vui->vui_timing_info_present_flag: %d \n", vui->vui_timing_info_present_flag); 
    if( vui->vui_timing_info_present_flag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_num_units_in_tick = bs_read_u(b, 32); printf("vui->vui_num_units_in_tick: %d \n", vui->vui_num_units_in_tick); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_time_scale = bs_read_u(b, 32); printf("vui->vui_time_scale: %d \n", vui->vui_time_scale); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_poc_proportional_to_timing_flag = bs_read_u1(b); printf("vui->vui_poc_proportional_to_timing_flag: %d \n", vui->vui_poc_proportional_to_timing_flag); 
        if( vui->vui_poc_proportional_to_timing_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_num_ticks_poc_diff_one_minus1 = bs_read_ue(b); printf("vui->vui_num_ticks_poc_diff_one_minus1: %d \n", vui->vui_num_ticks_poc_diff_one_minus1); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->vui_hrd_parameters_present_flag = bs_read_u1(b); printf("vui->vui_hrd_parameters_present_flag: %d \n", vui->vui_hrd_parameters_present_flag); 
        if( vui->vui_hrd_parameters_present_flag ) {
            read_debug_hevc_hrd_parameters( &vui->hrd, b, 1, sps->sps_max_sub_layers_minus1 );
        }
    }
    printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->bitstream_restriction_flag = bs_read_u1(b); printf("vui->bitstream_restriction_flag: %d \n", vui->bitstream_restriction_flag); 
    if( vui->bitstream_restriction_flag )
    {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->tiles_fixed_structure_flag = bs_read_u1(b); printf("vui->tiles_fixed_structure_flag: %d \n", vui->tiles_fixed_structure_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->motion_vectors_over_pic_boundaries_flag = bs_read_u1(b); printf("vui->motion_vectors_over_pic_boundaries_flag: %d \n", vui->motion_vectors_over_pic_boundaries_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->restricted_ref_pic_lists_flag = bs_read_u1(b); printf("vui->restricted_ref_pic_lists_flag: %d \n", vui->restricted_ref_pic_lists_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->min_spatial_segmentation_idc = bs_read_ue(b); printf("vui->min_spatial_segmentation_idc: %d \n", vui->min_spatial_segmentation_idc); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->max_bytes_per_pic_denom = bs_read_ue(b); printf("vui->max_bytes_per_pic_denom: %d \n", vui->max_bytes_per_pic_denom); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->max_bits_per_min_cu_denom = bs_read_ue(b); printf("vui->max_bits_per_min_cu_denom: %d \n", vui->max_bits_per_min_cu_denom); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->log2_max_mv_length_horizontal = bs_read_ue(b); printf("vui->log2_max_mv_length_horizontal: %d \n", vui->log2_max_mv_length_horizontal); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); vui->log2_max_mv_length_vertical = bs_read_ue(b); printf("vui->log2_max_mv_length_vertical: %d \n", vui->log2_max_mv_length_vertical); 
    }
}

//Appendix E.2.2 HRD parameters syntax
void read_debug_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1)
{
    if( commonInfPresentFlag ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->nal_hrd_parameters_present_flag = bs_read_u1(b); printf("hrd->nal_hrd_parameters_present_flag: %d \n", hrd->nal_hrd_parameters_present_flag); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->vcl_hrd_parameters_present_flag = bs_read_u1(b); printf("hrd->vcl_hrd_parameters_present_flag: %d \n", hrd->vcl_hrd_parameters_present_flag); 
        if( hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag ){
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->sub_pic_hrd_params_present_flag = bs_read_u1(b); printf("hrd->sub_pic_hrd_params_present_flag: %d \n", hrd->sub_pic_hrd_params_present_flag); 
            if( hrd->sub_pic_hrd_params_present_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->tick_divisor_minus2 = bs_read_u8(b); printf("hrd->tick_divisor_minus2: %d \n", hrd->tick_divisor_minus2); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->du_cpb_removal_delay_increment_length_minus1 = bs_read_u(b, 5); printf("hrd->du_cpb_removal_delay_increment_length_minus1: %d \n", hrd->du_cpb_removal_delay_increment_length_minus1); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->sub_pic_cpb_params_in_pic_timing_sei_flag = bs_read_u1(b); printf("hrd->sub_pic_cpb_params_in_pic_timing_sei_flag: %d \n", hrd->sub_pic_cpb_params_in_pic_timing_sei_flag); 
                printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->dpb_output_delay_du_length_minus1 = bs_read_u(b, 5); printf("hrd->dpb_output_delay_du_length_minus1: %d \n", hrd->dpb_output_delay_du_length_minus1); 
            }
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->bit_rate_scale = bs_read_u(b, 4); printf("hrd->bit_rate_scale: %d \n", hrd->bit_rate_scale); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->cpb_size_scale = bs_read_u(b, 4); printf("hrd->cpb_size_scale: %d \n", hrd->cpb_size_scale); 
            if( hrd->sub_pic_hrd_params_present_flag ) {
                printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->cpb_size_du_scale = bs_read_u(b, 4); printf("hrd->cpb_size_du_scale: %d \n", hrd->cpb_size_du_scale); 
            }
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->initial_cpb_removal_delay_length_minus1 = bs_read_u(b, 5); printf("hrd->initial_cpb_removal_delay_length_minus1: %d \n", hrd->initial_cpb_removal_delay_length_minus1); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->au_cpb_removal_delay_length_minus1 = bs_read_u(b, 5); printf("hrd->au_cpb_removal_delay_length_minus1: %d \n", hrd->au_cpb_removal_delay_length_minus1); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->dpb_output_delay_length_minus1 = bs_read_u(b, 5); printf("hrd->dpb_output_delay_length_minus1: %d \n", hrd->dpb_output_delay_length_minus1); 
        }
    }
    
    for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->fixed_pic_rate_general_flag[ i ] = bs_read_u1(b); printf("hrd->fixed_pic_rate_general_flag[ i ]: %d \n", hrd->fixed_pic_rate_general_flag[ i ]); 
        if( !hrd->fixed_pic_rate_general_flag[ i ] ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->fixed_pic_rate_within_cvs_flag[ i ] = bs_read_u1(b); printf("hrd->fixed_pic_rate_within_cvs_flag[ i ]: %d \n", hrd->fixed_pic_rate_within_cvs_flag[ i ]); 
        }
        if( hrd->fixed_pic_rate_within_cvs_flag[ i ] ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->elemental_duration_in_tc_minus1[ i ] = bs_read_ue(b); printf("hrd->elemental_duration_in_tc_minus1[ i ]: %d \n", hrd->elemental_duration_in_tc_minus1[ i ]); 
        } else {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->low_delay_hrd_flag[ i ] = bs_read_u1(b); printf("hrd->low_delay_hrd_flag[ i ]: %d \n", hrd->low_delay_hrd_flag[ i ]); 
        }
        if( hrd->low_delay_hrd_flag[ i ] ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); hrd->cpb_cnt_minus1[ i ] = bs_read_ue(b); printf("hrd->cpb_cnt_minus1[ i ]: %d \n", hrd->cpb_cnt_minus1[ i ]); 
        }
        if( hrd->nal_hrd_parameters_present_flag ) {
            read_debug_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_nal[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
        if( hrd->vcl_hrd_parameters_present_flag ) {
            read_debug_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_vcl[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
    }
}

//Appendix E.2.3 Sub-layer HRD parameters syntax
void read_debug_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag)
{
    for( int i = 0; i <= CpbCnt; i++ ) {
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sub_layer_hrd->bit_rate_value_minus1[i] = bs_read_ue(b); printf("sub_layer_hrd->bit_rate_value_minus1[i]: %d \n", sub_layer_hrd->bit_rate_value_minus1[i]); 
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sub_layer_hrd->cpb_size_value_minus1[i] = bs_read_ue(b); printf("sub_layer_hrd->cpb_size_value_minus1[i]: %d \n", sub_layer_hrd->cpb_size_value_minus1[i]); 
        if( sub_pic_hrd_params_present_flag ) {
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sub_layer_hrd->cpb_size_du_value_minus1[i] = bs_read_ue(b); printf("sub_layer_hrd->cpb_size_du_value_minus1[i]: %d \n", sub_layer_hrd->cpb_size_du_value_minus1[i]); 
            printf("%ld.%d: ", b->p - b->start, b->bits_left); sub_layer_hrd->bit_rate_du_value_minus1[i] = bs_read_ue(b); printf("sub_layer_hrd->bit_rate_du_value_minus1[i]: %d \n", sub_layer_hrd->bit_rate_du_value_minus1[i]); 
        }
        printf("%ld.%d: ", b->p - b->start, b->bits_left); sub_layer_hrd->cbr_flag[i] = bs_read_u1(b); printf("sub_layer_hrd->cbr_flag[i]: %d \n", sub_layer_hrd->cbr_flag[i]); 
    }
}


void write_debug_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b);
void write_debug_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b);
void write_debug_sei_rbsp(hevc_stream_t* h, bs_t* b);
void write_debug_sei_message(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_end_of_seq_rbsp();
void write_debug_end_of_bitstream_rbsp();
void write_debug_filler_data_rbsp(bs_t* b);
void write_debug_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b);
void write_debug_hevc_rbsp_slice_trailing_bits(bs_t* b);
void write_debug_hevc_rbsp_trailing_bits(bs_t* b);
void write_debug_hevc_byte_alignment(bs_t* b);
void write_debug_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1);
void write_debug_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b );
void write_debug_hevc_slice_header(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b);
void write_debug_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets );
void write_debug_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b);
void write_debug_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1);
void write_debug_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag);



//7.3.1 NAL unit syntax
int write_debug_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size)
{
    hevc_nal_t* nal = h->nal;

    int nal_size = size;
    int rbsp_size = size;
    uint8_t* rbsp_buf = (uint8_t*)calloc(1, rbsp_size);

    if( 0 )
    {
        int rc = nal_to_rbsp(buf, &nal_size, rbsp_buf, &rbsp_size);

        if (rc < 0) { free(rbsp_buf); return -1; } // handle conversion error
    }

    if( 1 )
    {
        rbsp_size = size*3/4; // NOTE this may have to be slightly smaller (3/4 smaller, worst case) in order to be guaranteed to fit
    }

    bs_t* b = bs_new(rbsp_buf, rbsp_size);
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    int forbidden_zero_bit = 1;  bs_write_u(b, forbidden_zero_bit, 0);
    printf("forbidden_zero_bit: %d ( %ld )\n", forbidden_zero_bit, decimal_to_binary( forbidden_zero_bit )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 6, nal->nal_unit_type);
    printf("nal->nal_unit_type: %d ( %ld )\n", nal->nal_unit_type, decimal_to_binary( nal->nal_unit_type )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 6, nal->nal_layer_id);
    printf("nal->nal_layer_id: %d ( %ld )\n", nal->nal_layer_id, decimal_to_binary( nal->nal_layer_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 3, nal->nal_temporal_id_plus1);
    printf("nal->nal_temporal_id_plus1: %d ( %ld )\n", nal->nal_temporal_id_plus1, decimal_to_binary( nal->nal_temporal_id_plus1 )); 

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
            
            write_debug_hevc_slice_layer_rbsp(h, b);
            break;

#ifdef HAVE_SEI
        case NAL_UNIT_TYPE_SEI:
            write_debug_hevc_sei_rbsp(h, b);
            break;
#endif

        case HEVC_NAL_UNIT_TYPE_VPS_NUT: 
            write_debug_hevc_video_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_SPS_NUT: 
            write_debug_hevc_seq_parameter_set_rbsp(h, b); 
            break;

        case HEVC_NAL_UNIT_TYPE_PPS_NUT:   
            write_debug_hevc_pic_parameter_set_rbsp(h, b);
            break;

        default:
            return -1;
    }

    if (bs_overrun(b)) { bs_free(b); free(rbsp_buf); return -1; }

    if( 1 )
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
void write_debug_hevc_video_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i, j;

    hevc_vps_t* vps = h->vps;
    if( 0 )
    {
        memset(vps, 0, sizeof(hevc_vps_t));
    }
 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
 
    bs_write_u(b, 4, vps->vps_video_parameter_set_id);
 
    printf("vps->vps_video_parameter_set_id: %d ( %ld )\n", vps->vps_video_parameter_set_id, decimal_to_binary( vps->vps_video_parameter_set_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vps->vps_base_layer_internal_flag);
    printf("vps->vps_base_layer_internal_flag: %d ( %ld )\n", vps->vps_base_layer_internal_flag, decimal_to_binary( vps->vps_base_layer_internal_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vps->vps_base_layer_available_flag);
    printf("vps->vps_base_layer_available_flag: %d ( %ld )\n", vps->vps_base_layer_available_flag, decimal_to_binary( vps->vps_base_layer_available_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 6, vps->vps_max_layers_minus1);
    printf("vps->vps_max_layers_minus1: %d ( %ld )\n", vps->vps_max_layers_minus1, decimal_to_binary( vps->vps_max_layers_minus1 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 3, vps->vps_max_sub_layers_minus1);
    printf("vps->vps_max_sub_layers_minus1: %d ( %ld )\n", vps->vps_max_sub_layers_minus1, decimal_to_binary( vps->vps_max_sub_layers_minus1 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vps->vps_temporal_id_nesting_flag);
    printf("vps->vps_temporal_id_nesting_flag: %d ( %ld )\n", vps->vps_temporal_id_nesting_flag, decimal_to_binary( vps->vps_temporal_id_nesting_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    int vps_reserved_0xffff_16bits = 16;  bs_write_u(b, vps_reserved_0xffff_16bits, 0xffff);
    printf("vps_reserved_0xffff_16bits: %d ( %ld )\n", vps_reserved_0xffff_16bits, decimal_to_binary( vps_reserved_0xffff_16bits )); 
    
    write_debug_hevc_profile_tier_level(&vps->ptl, b, 1, vps->vps_max_sub_layers_minus1); 
    
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    
    bs_write_u1(b, vps->vps_sub_layer_ordering_info_present_flag);
    
    printf("vps->vps_sub_layer_ordering_info_present_flag: %d ( %ld )\n", vps->vps_sub_layer_ordering_info_present_flag, decimal_to_binary( vps->vps_sub_layer_ordering_info_present_flag )); 
    for( i = ( vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1 ); 
            i <= vps->vps_max_sub_layers_minus1; i++ ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vps->vps_max_dec_pic_buffering_minus1[ i ]);
        printf("vps->vps_max_dec_pic_buffering_minus1[ i ]: %d ( %ld )\n", vps->vps_max_dec_pic_buffering_minus1[ i ], decimal_to_binary( vps->vps_max_dec_pic_buffering_minus1[ i ] ));         
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vps->vps_max_num_reorder_pics[ i ]);
        printf("vps->vps_max_num_reorder_pics[ i ]: %d ( %ld )\n", vps->vps_max_num_reorder_pics[ i ], decimal_to_binary( vps->vps_max_num_reorder_pics[ i ] ));         
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vps->vps_max_latency_increase_plus1[ i ]);
        printf("vps->vps_max_latency_increase_plus1[ i ]: %d ( %ld )\n", vps->vps_max_latency_increase_plus1[ i ], decimal_to_binary( vps->vps_max_latency_increase_plus1[ i ] ));         
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 6, vps->vps_max_layer_id);
    printf("vps->vps_max_layer_id: %d ( %ld )\n", vps->vps_max_layer_id, decimal_to_binary( vps->vps_max_layer_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, vps->vps_num_layer_sets_minus1);
    printf("vps->vps_num_layer_sets_minus1: %d ( %ld )\n", vps->vps_num_layer_sets_minus1, decimal_to_binary( vps->vps_num_layer_sets_minus1 )); 
    for( i = 1; i <= vps->vps_num_layer_sets_minus1; i++ )
        for( j = 0; j <= vps->vps_max_layer_id; j++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, vps->layer_id_included_flag[ i ][ j ]);
            printf("vps->layer_id_included_flag[ i ][ j ]: %d ( %ld )\n", vps->layer_id_included_flag[ i ][ j ], decimal_to_binary( vps->layer_id_included_flag[ i ][ j ] )); 
        }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vps->vps_timing_info_present_flag);
    printf("vps->vps_timing_info_present_flag: %d ( %ld )\n", vps->vps_timing_info_present_flag, decimal_to_binary( vps->vps_timing_info_present_flag )); 
    if( vps->vps_timing_info_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 32, vps->vps_num_units_in_tick);
        printf("vps->vps_num_units_in_tick: %d ( %ld )\n", vps->vps_num_units_in_tick, decimal_to_binary( vps->vps_num_units_in_tick )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 32, vps->vps_time_scale);
        printf("vps->vps_time_scale: %d ( %ld )\n", vps->vps_time_scale, decimal_to_binary( vps->vps_time_scale )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vps->vps_poc_proportional_to_timing_flag);
        printf("vps->vps_poc_proportional_to_timing_flag: %d ( %ld )\n", vps->vps_poc_proportional_to_timing_flag, decimal_to_binary( vps->vps_poc_proportional_to_timing_flag )); 
        if( vps->vps_poc_proportional_to_timing_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, vps->vps_num_ticks_poc_diff_one_minus1);
            printf("vps->vps_num_ticks_poc_diff_one_minus1: %d ( %ld )\n", vps->vps_num_ticks_poc_diff_one_minus1, decimal_to_binary( vps->vps_num_ticks_poc_diff_one_minus1 )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vps->vps_num_hrd_parameters);
        printf("vps->vps_num_hrd_parameters: %d ( %ld )\n", vps->vps_num_hrd_parameters, decimal_to_binary( vps->vps_num_hrd_parameters )); 
        for( i = 0; i < vps->vps_num_hrd_parameters; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, vps->hrd_layer_set_idx[ i ]);
            printf("vps->hrd_layer_set_idx[ i ]: %d ( %ld )\n", vps->hrd_layer_set_idx[ i ], decimal_to_binary( vps->hrd_layer_set_idx[ i ] )); 
            if (i > 0) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, vps->cprms_present_flag[ i ]);
                printf("vps->cprms_present_flag[ i ]: %d ( %ld )\n", vps->cprms_present_flag[ i ], decimal_to_binary( vps->cprms_present_flag[ i ] )); 
            }
            write_debug_hevc_hrd_parameters(&vps->hrd[i], b,
                                           vps->cprms_present_flag[ i ],
                                           vps->vps_max_sub_layers_minus1);
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vps->vps_extension_flag);
    printf("vps->vps_extension_flag: %d ( %ld )\n", vps->vps_extension_flag, decimal_to_binary( vps->vps_extension_flag )); 
    //TODO: support extension data
    //if (vps->vps_extension_flag)    

    write_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.2 Sequence parameter set RBSP syntax
void write_debug_hevc_seq_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;

    hevc_sps_t* sps = h->sps;
    if( 0 )
    {
        memset(sps, 0, sizeof(hevc_sps_t));
    }
 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
 
    bs_write_u(b, 4, sps->sps_video_parameter_set_id);
 
    printf("sps->sps_video_parameter_set_id: %d ( %ld )\n", sps->sps_video_parameter_set_id, decimal_to_binary( sps->sps_video_parameter_set_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 3, sps->sps_max_sub_layers_minus1);
    printf("sps->sps_max_sub_layers_minus1: %d ( %ld )\n", sps->sps_max_sub_layers_minus1, decimal_to_binary( sps->sps_max_sub_layers_minus1 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->sps_temporal_id_nesting_flag);
    printf("sps->sps_temporal_id_nesting_flag: %d ( %ld )\n", sps->sps_temporal_id_nesting_flag, decimal_to_binary( sps->sps_temporal_id_nesting_flag )); 
    write_debug_hevc_profile_tier_level(&sps->ptl, b, 1, sps->sps_max_sub_layers_minus1); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->sps_seq_parameter_set_id);
    printf("sps->sps_seq_parameter_set_id: %d ( %ld )\n", sps->sps_seq_parameter_set_id, decimal_to_binary( sps->sps_seq_parameter_set_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->chroma_format_idc);
    printf("sps->chroma_format_idc: %d ( %ld )\n", sps->chroma_format_idc, decimal_to_binary( sps->chroma_format_idc )); 
    if( sps->chroma_format_idc == 3 ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->separate_colour_plane_flag);
        printf("sps->separate_colour_plane_flag: %d ( %ld )\n", sps->separate_colour_plane_flag, decimal_to_binary( sps->separate_colour_plane_flag )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->pic_width_in_luma_samples);
    printf("sps->pic_width_in_luma_samples: %d ( %ld )\n", sps->pic_width_in_luma_samples, decimal_to_binary( sps->pic_width_in_luma_samples )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->pic_height_in_luma_samples);
    printf("sps->pic_height_in_luma_samples: %d ( %ld )\n", sps->pic_height_in_luma_samples, decimal_to_binary( sps->pic_height_in_luma_samples )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->conformance_window_flag);
    printf("sps->conformance_window_flag: %d ( %ld )\n", sps->conformance_window_flag, decimal_to_binary( sps->conformance_window_flag )); 
    if( sps->conformance_window_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->conf_win_left_offset);
        printf("sps->conf_win_left_offset: %d ( %ld )\n", sps->conf_win_left_offset, decimal_to_binary( sps->conf_win_left_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->conf_win_right_offset);
        printf("sps->conf_win_right_offset: %d ( %ld )\n", sps->conf_win_right_offset, decimal_to_binary( sps->conf_win_right_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->conf_win_top_offset);
        printf("sps->conf_win_top_offset: %d ( %ld )\n", sps->conf_win_top_offset, decimal_to_binary( sps->conf_win_top_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->conf_win_bottom_offset);
        printf("sps->conf_win_bottom_offset: %d ( %ld )\n", sps->conf_win_bottom_offset, decimal_to_binary( sps->conf_win_bottom_offset )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->bit_depth_luma_minus8);
    printf("sps->bit_depth_luma_minus8: %d ( %ld )\n", sps->bit_depth_luma_minus8, decimal_to_binary( sps->bit_depth_luma_minus8 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->bit_depth_chroma_minus8);
    printf("sps->bit_depth_chroma_minus8: %d ( %ld )\n", sps->bit_depth_chroma_minus8, decimal_to_binary( sps->bit_depth_chroma_minus8 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->log2_max_pic_order_cnt_lsb_minus4);
    printf("sps->log2_max_pic_order_cnt_lsb_minus4: %d ( %ld )\n", sps->log2_max_pic_order_cnt_lsb_minus4, decimal_to_binary( sps->log2_max_pic_order_cnt_lsb_minus4 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->sps_sub_layer_ordering_info_present_flag);
    printf("sps->sps_sub_layer_ordering_info_present_flag: %d ( %ld )\n", sps->sps_sub_layer_ordering_info_present_flag, decimal_to_binary( sps->sps_sub_layer_ordering_info_present_flag )); 
    for( i = ( sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1 ); 
            i <= sps->sps_max_sub_layers_minus1; i++ ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->sps_max_dec_pic_buffering_minus1 [ i ]);
        printf("sps->sps_max_dec_pic_buffering_minus1 [ i ]: %d ( %ld )\n", sps->sps_max_dec_pic_buffering_minus1 [ i ], decimal_to_binary( sps->sps_max_dec_pic_buffering_minus1 [ i ] )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->sps_max_num_reorder_pics [ i ]);
        printf("sps->sps_max_num_reorder_pics [ i ]: %d ( %ld )\n", sps->sps_max_num_reorder_pics [ i ], decimal_to_binary( sps->sps_max_num_reorder_pics [ i ] )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->sps_max_latency_increase_plus1 [ i ]);
        printf("sps->sps_max_latency_increase_plus1 [ i ]: %d ( %ld )\n", sps->sps_max_latency_increase_plus1 [ i ], decimal_to_binary( sps->sps_max_latency_increase_plus1 [ i ] )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->log2_min_luma_coding_block_size_minus3);
    printf("sps->log2_min_luma_coding_block_size_minus3: %d ( %ld )\n", sps->log2_min_luma_coding_block_size_minus3, decimal_to_binary( sps->log2_min_luma_coding_block_size_minus3 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->log2_diff_max_min_luma_coding_block_size);
    printf("sps->log2_diff_max_min_luma_coding_block_size: %d ( %ld )\n", sps->log2_diff_max_min_luma_coding_block_size, decimal_to_binary( sps->log2_diff_max_min_luma_coding_block_size )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->log2_min_luma_transform_block_size_minus2);
    printf("sps->log2_min_luma_transform_block_size_minus2: %d ( %ld )\n", sps->log2_min_luma_transform_block_size_minus2, decimal_to_binary( sps->log2_min_luma_transform_block_size_minus2 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->log2_diff_max_min_luma_transform_block_size);
    printf("sps->log2_diff_max_min_luma_transform_block_size: %d ( %ld )\n", sps->log2_diff_max_min_luma_transform_block_size, decimal_to_binary( sps->log2_diff_max_min_luma_transform_block_size )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->max_transform_hierarchy_depth_inter);
    printf("sps->max_transform_hierarchy_depth_inter: %d ( %ld )\n", sps->max_transform_hierarchy_depth_inter, decimal_to_binary( sps->max_transform_hierarchy_depth_inter )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->max_transform_hierarchy_depth_intra);
    printf("sps->max_transform_hierarchy_depth_intra: %d ( %ld )\n", sps->max_transform_hierarchy_depth_intra, decimal_to_binary( sps->max_transform_hierarchy_depth_intra )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->scaling_list_enabled_flag);
    printf("sps->scaling_list_enabled_flag: %d ( %ld )\n", sps->scaling_list_enabled_flag, decimal_to_binary( sps->scaling_list_enabled_flag )); 
    
    if( sps->scaling_list_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->sps_scaling_list_data_present_flag);
        printf("sps->sps_scaling_list_data_present_flag: %d ( %ld )\n", sps->sps_scaling_list_data_present_flag, decimal_to_binary( sps->sps_scaling_list_data_present_flag )); 
        if( sps->sps_scaling_list_data_present_flag ) {
            write_debug_hevc_scaling_list_data(&sps->scaling_list_data, b); 
        }
    }
    
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    
    bs_write_u1(b, sps->amp_enabled_flag);
    
    printf("sps->amp_enabled_flag: %d ( %ld )\n", sps->amp_enabled_flag, decimal_to_binary( sps->amp_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->sample_adaptive_offset_enabled_flag);
    printf("sps->sample_adaptive_offset_enabled_flag: %d ( %ld )\n", sps->sample_adaptive_offset_enabled_flag, decimal_to_binary( sps->sample_adaptive_offset_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->pcm_enabled_flag);
    printf("sps->pcm_enabled_flag: %d ( %ld )\n", sps->pcm_enabled_flag, decimal_to_binary( sps->pcm_enabled_flag )); 
    if( sps->pcm_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 4, sps->pcm_sample_bit_depth_luma_minus1);
        printf("sps->pcm_sample_bit_depth_luma_minus1: %d ( %ld )\n", sps->pcm_sample_bit_depth_luma_minus1, decimal_to_binary( sps->pcm_sample_bit_depth_luma_minus1 )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 4, sps->pcm_sample_bit_depth_chroma_minus1);
        printf("sps->pcm_sample_bit_depth_chroma_minus1: %d ( %ld )\n", sps->pcm_sample_bit_depth_chroma_minus1, decimal_to_binary( sps->pcm_sample_bit_depth_chroma_minus1 )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->log2_min_pcm_luma_coding_block_size_minus3);
        printf("sps->log2_min_pcm_luma_coding_block_size_minus3: %d ( %ld )\n", sps->log2_min_pcm_luma_coding_block_size_minus3, decimal_to_binary( sps->log2_min_pcm_luma_coding_block_size_minus3 )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->log2_diff_max_min_pcm_luma_coding_block_size);
        printf("sps->log2_diff_max_min_pcm_luma_coding_block_size: %d ( %ld )\n", sps->log2_diff_max_min_pcm_luma_coding_block_size, decimal_to_binary( sps->log2_diff_max_min_pcm_luma_coding_block_size )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->pcm_loop_filter_disabled_flag);
        printf("sps->pcm_loop_filter_disabled_flag: %d ( %ld )\n", sps->pcm_loop_filter_disabled_flag, decimal_to_binary( sps->pcm_loop_filter_disabled_flag )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sps->num_short_term_ref_pic_sets);
    printf("sps->num_short_term_ref_pic_sets: %d ( %ld )\n", sps->num_short_term_ref_pic_sets, decimal_to_binary( sps->num_short_term_ref_pic_sets )); 
    for( i = 0; i < sps->num_short_term_ref_pic_sets; i++) {
        write_debug_hevc_st_ref_pic_set(&sps->st_ref_pic_set[i], b, i, sps->num_short_term_ref_pic_sets);
    }
    
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    
    bs_write_u1(b, sps->long_term_ref_pics_present_flag);
    
    printf("sps->long_term_ref_pics_present_flag: %d ( %ld )\n", sps->long_term_ref_pics_present_flag, decimal_to_binary( sps->long_term_ref_pics_present_flag )); 
    if( sps->long_term_ref_pics_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sps->num_long_term_ref_pics_sps);
        printf("sps->num_long_term_ref_pics_sps: %d ( %ld )\n", sps->num_long_term_ref_pics_sps, decimal_to_binary( sps->num_long_term_ref_pics_sps )); 
        for( i = 0; i < sps->num_long_term_ref_pics_sps; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sps->lt_ref_pic_poc_lsb_sps[ i ]);
            printf("sps->lt_ref_pic_poc_lsb_sps[ i ]: %d ( %ld )\n", sps->lt_ref_pic_poc_lsb_sps[ i ], decimal_to_binary( sps->lt_ref_pic_poc_lsb_sps[ i ] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sps->used_by_curr_pic_lt_sps_flag[ i ]);
            printf("sps->used_by_curr_pic_lt_sps_flag[ i ]: %d ( %ld )\n", sps->used_by_curr_pic_lt_sps_flag[ i ], decimal_to_binary( sps->used_by_curr_pic_lt_sps_flag[ i ] )); 
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->sps_temporal_mvp_enabled_flag);
    printf("sps->sps_temporal_mvp_enabled_flag: %d ( %ld )\n", sps->sps_temporal_mvp_enabled_flag, decimal_to_binary( sps->sps_temporal_mvp_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->strong_intra_smoothing_enabled_flag);
    printf("sps->strong_intra_smoothing_enabled_flag: %d ( %ld )\n", sps->strong_intra_smoothing_enabled_flag, decimal_to_binary( sps->strong_intra_smoothing_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->vui_parameters_present_flag);
    printf("sps->vui_parameters_present_flag: %d ( %ld )\n", sps->vui_parameters_present_flag, decimal_to_binary( sps->vui_parameters_present_flag )); 
    if( sps->vui_parameters_present_flag ) {
        write_debug_hevc_vui_parameters(sps, b);
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps->sps_extension_present_flag);
    printf("sps->sps_extension_present_flag: %d ( %ld )\n", sps->sps_extension_present_flag, decimal_to_binary( sps->sps_extension_present_flag )); 
    
    if( sps->sps_extension_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->sps_range_extension_flag);
        printf("sps->sps_range_extension_flag: %d ( %ld )\n", sps->sps_range_extension_flag, decimal_to_binary( sps->sps_range_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->sps_multilayer_extension_flag);
        printf("sps->sps_multilayer_extension_flag: %d ( %ld )\n", sps->sps_multilayer_extension_flag, decimal_to_binary( sps->sps_multilayer_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sps->sps_3d_extension_flag);
        printf("sps->sps_3d_extension_flag: %d ( %ld )\n", sps->sps_3d_extension_flag, decimal_to_binary( sps->sps_3d_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 5, sps->sps_extension_5bits);
        printf("sps->sps_extension_5bits: %d ( %ld )\n", sps->sps_extension_5bits, decimal_to_binary( sps->sps_extension_5bits )); 
    }
    if( sps->sps_range_extension_flag ) {
        write_debug_hevc_sps_range_extension( &sps->sps_range_ext, b);
    }
    
    if( 0 )
    {
        memcpy(h->sps_table[sps->sps_seq_parameter_set_id], h->sps, sizeof(hevc_sps_t));
    }
}

//7.3.2.2.2 Sequence parameter set range extension syntax
void write_debug_hevc_sps_range_extension(hevc_sps_range_ext_t* sps_range_ext, bs_t* b)
{
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->transform_skip_rotation_enabled_flag);
    printf("sps_range_ext->transform_skip_rotation_enabled_flag: %d ( %ld )\n", sps_range_ext->transform_skip_rotation_enabled_flag, decimal_to_binary( sps_range_ext->transform_skip_rotation_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->transform_skip_context_enabled_flag);
    printf("sps_range_ext->transform_skip_context_enabled_flag: %d ( %ld )\n", sps_range_ext->transform_skip_context_enabled_flag, decimal_to_binary( sps_range_ext->transform_skip_context_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->implicit_rdpcm_enabled_flag);
    printf("sps_range_ext->implicit_rdpcm_enabled_flag: %d ( %ld )\n", sps_range_ext->implicit_rdpcm_enabled_flag, decimal_to_binary( sps_range_ext->implicit_rdpcm_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->explicit_rdpcm_enabled_flag);
    printf("sps_range_ext->explicit_rdpcm_enabled_flag: %d ( %ld )\n", sps_range_ext->explicit_rdpcm_enabled_flag, decimal_to_binary( sps_range_ext->explicit_rdpcm_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->extended_precision_processing_flag);
    printf("sps_range_ext->extended_precision_processing_flag: %d ( %ld )\n", sps_range_ext->extended_precision_processing_flag, decimal_to_binary( sps_range_ext->extended_precision_processing_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->intra_smoothing_disabled_flag);
    printf("sps_range_ext->intra_smoothing_disabled_flag: %d ( %ld )\n", sps_range_ext->intra_smoothing_disabled_flag, decimal_to_binary( sps_range_ext->intra_smoothing_disabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->high_precision_offsets_enabled_flag);
    printf("sps_range_ext->high_precision_offsets_enabled_flag: %d ( %ld )\n", sps_range_ext->high_precision_offsets_enabled_flag, decimal_to_binary( sps_range_ext->high_precision_offsets_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->persistent_rice_adaptation_enabled_flag);
    printf("sps_range_ext->persistent_rice_adaptation_enabled_flag: %d ( %ld )\n", sps_range_ext->persistent_rice_adaptation_enabled_flag, decimal_to_binary( sps_range_ext->persistent_rice_adaptation_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, sps_range_ext->cabac_bypass_alignment_enabled_flag);
    printf("sps_range_ext->cabac_bypass_alignment_enabled_flag: %d ( %ld )\n", sps_range_ext->cabac_bypass_alignment_enabled_flag, decimal_to_binary( sps_range_ext->cabac_bypass_alignment_enabled_flag )); 
}


//7.3.2.3 Picture parameter set RBSP syntax
void write_debug_hevc_pic_parameter_set_rbsp(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_pps_t* pps = h->pps;
    if( 0 )
    {
        memset(pps, 0, sizeof(hevc_pps_t));
    }

    printf("%d.%d: ", b->p - b->start, b->bits_left); 

    bs_write_ue(b, pps->pic_parameter_set_id);

    printf("pps->pic_parameter_set_id: %d ( %ld )\n", pps->pic_parameter_set_id, decimal_to_binary( pps->pic_parameter_set_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps->seq_parameter_set_id);
    printf("pps->seq_parameter_set_id: %d ( %ld )\n", pps->seq_parameter_set_id, decimal_to_binary( pps->seq_parameter_set_id )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->dependent_slice_segments_enabled_flag);
    printf("pps->dependent_slice_segments_enabled_flag: %d ( %ld )\n", pps->dependent_slice_segments_enabled_flag, decimal_to_binary( pps->dependent_slice_segments_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->output_flag_present_flag);
    printf("pps->output_flag_present_flag: %d ( %ld )\n", pps->output_flag_present_flag, decimal_to_binary( pps->output_flag_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b,  3 , pps->num_extra_slice_header_bits);
    printf("pps->num_extra_slice_header_bits: %d ( %ld )\n", pps->num_extra_slice_header_bits, decimal_to_binary( pps->num_extra_slice_header_bits )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->sign_data_hiding_enabled_flag);
    printf("pps->sign_data_hiding_enabled_flag: %d ( %ld )\n", pps->sign_data_hiding_enabled_flag, decimal_to_binary( pps->sign_data_hiding_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->cabac_init_present_flag);
    printf("pps->cabac_init_present_flag: %d ( %ld )\n", pps->cabac_init_present_flag, decimal_to_binary( pps->cabac_init_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps->num_ref_idx_l0_default_active_minus1);
    printf("pps->num_ref_idx_l0_default_active_minus1: %d ( %ld )\n", pps->num_ref_idx_l0_default_active_minus1, decimal_to_binary( pps->num_ref_idx_l0_default_active_minus1 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps->num_ref_idx_l1_default_active_minus1);
    printf("pps->num_ref_idx_l1_default_active_minus1: %d ( %ld )\n", pps->num_ref_idx_l1_default_active_minus1, decimal_to_binary( pps->num_ref_idx_l1_default_active_minus1 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_se(b, pps->init_qp_minus26);
    printf("pps->init_qp_minus26: %d ( %ld )\n", pps->init_qp_minus26, decimal_to_binary( pps->init_qp_minus26 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->constrained_intra_pred_flag);
    printf("pps->constrained_intra_pred_flag: %d ( %ld )\n", pps->constrained_intra_pred_flag, decimal_to_binary( pps->constrained_intra_pred_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->transform_skip_enabled_flag);
    printf("pps->transform_skip_enabled_flag: %d ( %ld )\n", pps->transform_skip_enabled_flag, decimal_to_binary( pps->transform_skip_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->cu_qp_delta_enabled_flag);
    printf("pps->cu_qp_delta_enabled_flag: %d ( %ld )\n", pps->cu_qp_delta_enabled_flag, decimal_to_binary( pps->cu_qp_delta_enabled_flag )); 
    if( pps->cu_qp_delta_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps->diff_cu_qp_delta_depth);
        printf("pps->diff_cu_qp_delta_depth: %d ( %ld )\n", pps->diff_cu_qp_delta_depth, decimal_to_binary( pps->diff_cu_qp_delta_depth )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_se(b, pps->pps_cb_qp_offset);
    printf("pps->pps_cb_qp_offset: %d ( %ld )\n", pps->pps_cb_qp_offset, decimal_to_binary( pps->pps_cb_qp_offset )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_se(b, pps->pps_cr_qp_offset);
    printf("pps->pps_cr_qp_offset: %d ( %ld )\n", pps->pps_cr_qp_offset, decimal_to_binary( pps->pps_cr_qp_offset )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->pps_slice_chroma_qp_offsets_present_flag);
    printf("pps->pps_slice_chroma_qp_offsets_present_flag: %d ( %ld )\n", pps->pps_slice_chroma_qp_offsets_present_flag, decimal_to_binary( pps->pps_slice_chroma_qp_offsets_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->weighted_pred_flag);
    printf("pps->weighted_pred_flag: %d ( %ld )\n", pps->weighted_pred_flag, decimal_to_binary( pps->weighted_pred_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->weighted_bipred_flag);
    printf("pps->weighted_bipred_flag: %d ( %ld )\n", pps->weighted_bipred_flag, decimal_to_binary( pps->weighted_bipred_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->transquant_bypass_enabled_flag);
    printf("pps->transquant_bypass_enabled_flag: %d ( %ld )\n", pps->transquant_bypass_enabled_flag, decimal_to_binary( pps->transquant_bypass_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->tiles_enabled_flag);
    printf("pps->tiles_enabled_flag: %d ( %ld )\n", pps->tiles_enabled_flag, decimal_to_binary( pps->tiles_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->entropy_coding_sync_enabled_flag);
    printf("pps->entropy_coding_sync_enabled_flag: %d ( %ld )\n", pps->entropy_coding_sync_enabled_flag, decimal_to_binary( pps->entropy_coding_sync_enabled_flag )); 
    if( pps->tiles_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps->num_tile_columns_minus1);
        printf("pps->num_tile_columns_minus1: %d ( %ld )\n", pps->num_tile_columns_minus1, decimal_to_binary( pps->num_tile_columns_minus1 )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps->num_tile_rows_minus1);
        printf("pps->num_tile_rows_minus1: %d ( %ld )\n", pps->num_tile_rows_minus1, decimal_to_binary( pps->num_tile_rows_minus1 )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->uniform_spacing_flag);
        printf("pps->uniform_spacing_flag: %d ( %ld )\n", pps->uniform_spacing_flag, decimal_to_binary( pps->uniform_spacing_flag )); 
        if( !pps->uniform_spacing_flag ) {
            for( i = 0; i < pps->num_tile_columns_minus1; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_ue(b, pps->column_width_minus1[ i ]);
                printf("pps->column_width_minus1[ i ]: %d ( %ld )\n", pps->column_width_minus1[ i ], decimal_to_binary( pps->column_width_minus1[ i ] )); 
            }
            for( i = 0; i < pps->num_tile_rows_minus1; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_ue(b, pps->row_height_minus1[ i ]);
                printf("pps->row_height_minus1[ i ]: %d ( %ld )\n", pps->row_height_minus1[ i ], decimal_to_binary( pps->row_height_minus1[ i ] )); 
            }
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->loop_filter_across_tiles_enabled_flag);
        printf("pps->loop_filter_across_tiles_enabled_flag: %d ( %ld )\n", pps->loop_filter_across_tiles_enabled_flag, decimal_to_binary( pps->loop_filter_across_tiles_enabled_flag )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->pps_loop_filter_across_slices_enabled_flag);
    printf("pps->pps_loop_filter_across_slices_enabled_flag: %d ( %ld )\n", pps->pps_loop_filter_across_slices_enabled_flag, decimal_to_binary( pps->pps_loop_filter_across_slices_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->deblocking_filter_control_present_flag);
    printf("pps->deblocking_filter_control_present_flag: %d ( %ld )\n", pps->deblocking_filter_control_present_flag, decimal_to_binary( pps->deblocking_filter_control_present_flag )); 
    if( pps->deblocking_filter_control_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->deblocking_filter_override_enabled_flag);
        printf("pps->deblocking_filter_override_enabled_flag: %d ( %ld )\n", pps->deblocking_filter_override_enabled_flag, decimal_to_binary( pps->deblocking_filter_override_enabled_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->pps_deblocking_filter_disabled_flag);
        printf("pps->pps_deblocking_filter_disabled_flag: %d ( %ld )\n", pps->pps_deblocking_filter_disabled_flag, decimal_to_binary( pps->pps_deblocking_filter_disabled_flag )); 
        if( pps->pps_deblocking_filter_disabled_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pps->pps_beta_offset_div2);
            printf("pps->pps_beta_offset_div2: %d ( %ld )\n", pps->pps_beta_offset_div2, decimal_to_binary( pps->pps_beta_offset_div2 )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pps->pps_tc_offset_div2);
            printf("pps->pps_tc_offset_div2: %d ( %ld )\n", pps->pps_tc_offset_div2, decimal_to_binary( pps->pps_tc_offset_div2 )); 
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->pps_scaling_list_data_present_flag);
    printf("pps->pps_scaling_list_data_present_flag: %d ( %ld )\n", pps->pps_scaling_list_data_present_flag, decimal_to_binary( pps->pps_scaling_list_data_present_flag )); 
    if( pps->pps_scaling_list_data_present_flag ) {
        write_debug_hevc_scaling_list_data(&pps->scaling_list_data, b);
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->lists_modification_present_flag);
    printf("pps->lists_modification_present_flag: %d ( %ld )\n", pps->lists_modification_present_flag, decimal_to_binary( pps->lists_modification_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps->log2_parallel_merge_level_minus2);
    printf("pps->log2_parallel_merge_level_minus2: %d ( %ld )\n", pps->log2_parallel_merge_level_minus2, decimal_to_binary( pps->log2_parallel_merge_level_minus2 )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->slice_segment_header_extension_present_flag);
    printf("pps->slice_segment_header_extension_present_flag: %d ( %ld )\n", pps->slice_segment_header_extension_present_flag, decimal_to_binary( pps->slice_segment_header_extension_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps->pps_extension_present_flag);
    printf("pps->pps_extension_present_flag: %d ( %ld )\n", pps->pps_extension_present_flag, decimal_to_binary( pps->pps_extension_present_flag )); 
    if( pps->pps_extension_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->pps_range_extension_flag);
        printf("pps->pps_range_extension_flag: %d ( %ld )\n", pps->pps_range_extension_flag, decimal_to_binary( pps->pps_range_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->pps_multilayer_extension_flag);
        printf("pps->pps_multilayer_extension_flag: %d ( %ld )\n", pps->pps_multilayer_extension_flag, decimal_to_binary( pps->pps_multilayer_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->pps_3d_extension_flag);
        printf("pps->pps_3d_extension_flag: %d ( %ld )\n", pps->pps_3d_extension_flag, decimal_to_binary( pps->pps_3d_extension_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pps->pps_extension_5bits);
        printf("pps->pps_extension_5bits: %d ( %ld )\n", pps->pps_extension_5bits, decimal_to_binary( pps->pps_extension_5bits )); 
    }
    if( pps->pps_range_extension_flag ) {
        write_debug_hevc_pps_range_extension( pps, b);
    }

    write_debug_hevc_rbsp_trailing_bits(b);

    if( 0 )
    {
        memcpy(h->pps_table[pps->pic_parameter_set_id], h->pps, sizeof(hevc_pps_t));
    }
}

//7.3.2.3.2 Picture parameter set range extension syntax
void write_debug_hevc_pps_range_extension(hevc_pps_t* pps, bs_t* b)
{
    hevc_pps_range_ext_t *pps_range_ext = &pps->pps_range_ext;;
    if( pps->transform_skip_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps_range_ext->log2_max_transform_skip_block_size_minus2);
        printf("pps_range_ext->log2_max_transform_skip_block_size_minus2: %d ( %ld )\n", pps_range_ext->log2_max_transform_skip_block_size_minus2, decimal_to_binary( pps_range_ext->log2_max_transform_skip_block_size_minus2 )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps_range_ext->cross_component_prediction_enabled_flag);
    printf("pps_range_ext->cross_component_prediction_enabled_flag: %d ( %ld )\n", pps_range_ext->cross_component_prediction_enabled_flag, decimal_to_binary( pps_range_ext->cross_component_prediction_enabled_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, pps_range_ext->chroma_qp_offset_list_enabled_flag);
    printf("pps_range_ext->chroma_qp_offset_list_enabled_flag: %d ( %ld )\n", pps_range_ext->chroma_qp_offset_list_enabled_flag, decimal_to_binary( pps_range_ext->chroma_qp_offset_list_enabled_flag )); 
    if( pps_range_ext->chroma_qp_offset_list_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps_range_ext->diff_cu_chroma_qp_offset_depth);
        printf("pps_range_ext->diff_cu_chroma_qp_offset_depth: %d ( %ld )\n", pps_range_ext->diff_cu_chroma_qp_offset_depth, decimal_to_binary( pps_range_ext->diff_cu_chroma_qp_offset_depth )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, pps_range_ext->chroma_qp_offset_list_len_minus1);
        printf("pps_range_ext->chroma_qp_offset_list_len_minus1: %d ( %ld )\n", pps_range_ext->chroma_qp_offset_list_len_minus1, decimal_to_binary( pps_range_ext->chroma_qp_offset_list_len_minus1 )); 
        for( int i = 0; i <= pps_range_ext->chroma_qp_offset_list_len_minus1; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pps_range_ext->cb_qp_offset_list[ i ]);
            printf("pps_range_ext->cb_qp_offset_list[ i ]: %d ( %ld )\n", pps_range_ext->cb_qp_offset_list[ i ], decimal_to_binary( pps_range_ext->cb_qp_offset_list[ i ] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pps_range_ext->cr_qp_offset_list[ i ]);
            printf("pps_range_ext->cr_qp_offset_list[ i ]: %d ( %ld )\n", pps_range_ext->cr_qp_offset_list[ i ], decimal_to_binary( pps_range_ext->cr_qp_offset_list[ i ] )); 
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps_range_ext->log2_sao_offset_scale_luma);
    printf("pps_range_ext->log2_sao_offset_scale_luma: %d ( %ld )\n", pps_range_ext->log2_sao_offset_scale_luma, decimal_to_binary( pps_range_ext->log2_sao_offset_scale_luma )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, pps_range_ext->log2_sao_offset_scale_chroma);
    printf("pps_range_ext->log2_sao_offset_scale_chroma: %d ( %ld )\n", pps_range_ext->log2_sao_offset_scale_chroma, decimal_to_binary( pps_range_ext->log2_sao_offset_scale_chroma )); 
}

#ifdef HAVE_SEI
//7.3.2.4 Supplemental enhancement information RBSP syntax
void write_debug_sei_rbsp(hevc_stream_t* h, bs_t* b)
{
    if( 0 )
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
            write_debug_sei_message(h, b);
        } while( more_rbsp_data(h, b) );
    }

    if( 1 )
    {
        for (int i = 0; i < h->num_seis; i++) {
            h->sei = h->seis[i];
            write_debug_sei_message(h, b);
        }
        h->sei = NULL;
    }

    write_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.5 Supplemental enhancement information message syntax
void write_debug_sei_message(hevc_stream_t* h, bs_t* b)
{
    if( 1 )
    {
        _write_ff_coded_number(b, h->sei->payloadType);
        _write_ff_coded_number(b, h->sei->payloadSize);
    }
    if( 0 )
    {
        h->sei->payloadType = _read_ff_coded_number(b);
        h->sei->payloadSize = _read_ff_coded_number(b);
    }
    write_debug_sei_payload( h, b, h->sei->payloadType, h->sei->payloadSize );
}
#endif

//7.3.2.5 Access unit delimiter RBSP syntax
void write_debug_hevc_access_unit_delimiter_rbsp(hevc_stream_t* h, bs_t* b)
{
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u(b, 3, h->aud->primary_pic_type);
    printf("h->aud->primary_pic_type: %d ( %ld )\n", h->aud->primary_pic_type, decimal_to_binary( h->aud->primary_pic_type )); 
    write_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.6 End of sequence RBSP syntax
void write_debug_hevc_end_of_seq_rbsp()
{
}

//7.3.2.7 End of bitstream RBSP syntax
void write_debug_end_of_bitstream_rbsp()
{
}

//7.3.2.8 Filler data RBSP syntax
void write_debug_filler_data_rbsp(bs_t* b)
{
    while( bs_next_bits(b, 8) == 0xFF )
    {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        int ff_byte = 8;  bs_write_u(b, ff_byte, 0xFF);
        printf("ff_byte: %d ( %ld )\n", ff_byte, decimal_to_binary( ff_byte )); 
    }
    write_debug_hevc_rbsp_trailing_bits(b);
}

//7.3.2.9 Slice segment layer RBSP syntax
void write_debug_hevc_slice_layer_rbsp(hevc_stream_t* h,  bs_t* b)
{
    write_debug_hevc_slice_header(h, b);
    hevc_slice_data_rbsp_t* slice_data = h->slice_data;

    if ( slice_data != NULL )
    {
        if ( slice_data->rbsp_buf != NULL ) free( slice_data->rbsp_buf ); 
        uint8_t *sptr = b->p + (!!b->bits_left); // CABAC-specific: skip alignment bits, if there are any
        slice_data->rbsp_size = b->end - sptr;
        
        slice_data->rbsp_buf = (uint8_t*)malloc(slice_data->rbsp_size);
        memcpy( slice_data->rbsp_buf, sptr, slice_data->rbsp_size );
    }

    //write_debug_hevc_slice_data(h, b); /* all categories of slice_data( ) syntax */
    write_debug_hevc_rbsp_slice_trailing_bits( b );
}

//7.3.2.10 RBSP slice trailing bits syntax
void write_debug_hevc_rbsp_slice_trailing_bits(bs_t* b)
{
    write_debug_hevc_rbsp_trailing_bits(b);
    //while( more_rbsp_trailing_data(b) )
    //{
    //    value( cabac_zero_word, f(16, 0x0000) );
    //}
}

//7.3.2.11 RBSP trailing bits syntax
void write_debug_hevc_rbsp_trailing_bits(bs_t* b)
{
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    int rbsp_stop_one_bit = 1;  bs_write_u(b, rbsp_stop_one_bit, 1);
    printf("rbsp_stop_one_bit: %d ( %ld )\n", rbsp_stop_one_bit, decimal_to_binary( rbsp_stop_one_bit )); 

    while( !bs_byte_aligned(b) )
    {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        int rbsp_alignment_zero_bit = 1;  bs_write_u(b, rbsp_alignment_zero_bit, 0);
        printf("rbsp_alignment_zero_bit: %d ( %ld )\n", rbsp_alignment_zero_bit, decimal_to_binary( rbsp_alignment_zero_bit )); 
    }
}

//7.3.2.12 Byte alignment syntax
void write_debug_hevc_byte_alignment(bs_t* b)
{
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    int alignment_bit_equal_to_one = 1;  bs_write_u(b, alignment_bit_equal_to_one, 1);
    printf("alignment_bit_equal_to_one: %d ( %ld )\n", alignment_bit_equal_to_one, decimal_to_binary( alignment_bit_equal_to_one )); 

    while( !bs_byte_aligned(b) )
    {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        int alignment_bit_equal_to_zero = 1;  bs_write_u(b, alignment_bit_equal_to_zero, 0);
        printf("alignment_bit_equal_to_zero: %d ( %ld )\n", alignment_bit_equal_to_zero, decimal_to_binary( alignment_bit_equal_to_zero )); 
    }
}

//7.3.3 Profile, tier and level syntax
void write_debug_hevc_profile_tier_level(hevc_profile_tier_level_t* ptl, bs_t* b, int profilePresentFlag, int maxNumSubLayersMinus1)
{
    int i, j;
    if( profilePresentFlag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 2, ptl->general_profile_space);
        printf("ptl->general_profile_space: %d ( %ld )\n", ptl->general_profile_space, decimal_to_binary( ptl->general_profile_space )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, ptl->general_tier_flag);
        printf("ptl->general_tier_flag: %d ( %ld )\n", ptl->general_tier_flag, decimal_to_binary( ptl->general_tier_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 5, ptl->general_profile_idc);
        printf("ptl->general_profile_idc: %d ( %ld )\n", ptl->general_profile_idc, decimal_to_binary( ptl->general_profile_idc )); 
        for( i = 0; i < 32; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_profile_compatibility_flag[ i ]);
            printf("ptl->general_profile_compatibility_flag[ i ]: %d ( %ld )\n", ptl->general_profile_compatibility_flag[ i ], decimal_to_binary( ptl->general_profile_compatibility_flag[ i ] )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, ptl->general_progressive_source_flag);
        printf("ptl->general_progressive_source_flag: %d ( %ld )\n", ptl->general_progressive_source_flag, decimal_to_binary( ptl->general_progressive_source_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, ptl->general_interlaced_source_flag);
        printf("ptl->general_interlaced_source_flag: %d ( %ld )\n", ptl->general_interlaced_source_flag, decimal_to_binary( ptl->general_interlaced_source_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, ptl->general_non_packed_constraint_flag);
        printf("ptl->general_non_packed_constraint_flag: %d ( %ld )\n", ptl->general_non_packed_constraint_flag, decimal_to_binary( ptl->general_non_packed_constraint_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, ptl->general_frame_only_constraint_flag);
        printf("ptl->general_frame_only_constraint_flag: %d ( %ld )\n", ptl->general_frame_only_constraint_flag, decimal_to_binary( ptl->general_frame_only_constraint_flag )); 
        if( ptl->general_profile_idc == 4 || ptl->general_profile_compatibility_flag[ 4 ] || 
            ptl->general_profile_idc == 5 || ptl->general_profile_compatibility_flag[ 5 ] || 
            ptl->general_profile_idc == 6 || ptl->general_profile_compatibility_flag[ 6 ] || 
            ptl->general_profile_idc == 7 || ptl->general_profile_compatibility_flag[ 7 ] ) {
                
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
                
            bs_write_u1(b, ptl->general_max_12bit_constraint_flag);
                
            printf("ptl->general_max_12bit_constraint_flag: %d ( %ld )\n", ptl->general_max_12bit_constraint_flag, decimal_to_binary( ptl->general_max_12bit_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_max_10bit_constraint_flag);
            printf("ptl->general_max_10bit_constraint_flag: %d ( %ld )\n", ptl->general_max_10bit_constraint_flag, decimal_to_binary( ptl->general_max_10bit_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_max_8bit_constraint_flag);
            printf("ptl->general_max_8bit_constraint_flag: %d ( %ld )\n", ptl->general_max_8bit_constraint_flag, decimal_to_binary( ptl->general_max_8bit_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_max_422chroma_constraint_flag);
            printf("ptl->general_max_422chroma_constraint_flag: %d ( %ld )\n", ptl->general_max_422chroma_constraint_flag, decimal_to_binary( ptl->general_max_422chroma_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_max_420chroma_constraint_flag);
            printf("ptl->general_max_420chroma_constraint_flag: %d ( %ld )\n", ptl->general_max_420chroma_constraint_flag, decimal_to_binary( ptl->general_max_420chroma_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_max_monochrome_constraint_flag);
            printf("ptl->general_max_monochrome_constraint_flag: %d ( %ld )\n", ptl->general_max_monochrome_constraint_flag, decimal_to_binary( ptl->general_max_monochrome_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_intra_constraint_flag);
            printf("ptl->general_intra_constraint_flag: %d ( %ld )\n", ptl->general_intra_constraint_flag, decimal_to_binary( ptl->general_intra_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_one_picture_only_constraint_flag);
            printf("ptl->general_one_picture_only_constraint_flag: %d ( %ld )\n", ptl->general_one_picture_only_constraint_flag, decimal_to_binary( ptl->general_one_picture_only_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->general_lower_bit_rate_constraint_flag);
            printf("ptl->general_lower_bit_rate_constraint_flag: %d ( %ld )\n", ptl->general_lower_bit_rate_constraint_flag, decimal_to_binary( ptl->general_lower_bit_rate_constraint_flag )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            int general_reserved_zero_34bits = 34;  bs_write_u(b, general_reserved_zero_34bits, 0);
            printf("general_reserved_zero_34bits: %d ( %ld )\n", general_reserved_zero_34bits, decimal_to_binary( general_reserved_zero_34bits )); 
        } else {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            int general_reserved_zero_43bits = 43;  bs_write_u(b, general_reserved_zero_43bits, 0);
            printf("general_reserved_zero_43bits: %d ( %ld )\n", general_reserved_zero_43bits, decimal_to_binary( general_reserved_zero_43bits )); 
        }
        if( ( ptl->general_profile_idc >= 1 && ptl->general_profile_idc <= 5 ) ||
              ptl->general_profile_compatibility_flag[ 1 ] ||
              ptl->general_profile_compatibility_flag[ 2 ] ||
              ptl->general_profile_compatibility_flag[ 3 ] ||
              ptl->general_profile_compatibility_flag[ 4 ] ||
              ptl->general_profile_compatibility_flag[ 5 ] ) {

            printf("%d.%d: ", b->p - b->start, b->bits_left); 

            bs_write_u1(b, ptl->general_inbld_flag);

            printf("ptl->general_inbld_flag: %d ( %ld )\n", ptl->general_inbld_flag, decimal_to_binary( ptl->general_inbld_flag )); 
        } else {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            int general_reserved_zero_bit = 1;  bs_write_u(b, general_reserved_zero_bit, 0);
            printf("general_reserved_zero_bit: %d ( %ld )\n", general_reserved_zero_bit, decimal_to_binary( general_reserved_zero_bit )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u8(b, ptl->general_level_idc);
        printf("ptl->general_level_idc: %d ( %ld )\n", ptl->general_level_idc, decimal_to_binary( ptl->general_level_idc )); 
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->sub_layer_profile_present_flag[ i ]);
            printf("ptl->sub_layer_profile_present_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_profile_present_flag[ i ], decimal_to_binary( ptl->sub_layer_profile_present_flag[ i ] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, ptl->sub_layer_level_present_flag[ i ]);
            printf("ptl->sub_layer_level_present_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_level_present_flag[ i ], decimal_to_binary( ptl->sub_layer_level_present_flag[ i ] )); 
        }
        if( maxNumSubLayersMinus1 > 0 ) {
            for( i = maxNumSubLayersMinus1; i < 8; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                int reserved_zero_xxbits = 2;  bs_write_u(b, reserved_zero_xxbits, 0);
                printf("reserved_zero_xxbits: %d ( %ld )\n", reserved_zero_xxbits, decimal_to_binary( reserved_zero_xxbits )); 
            }
        }
        for( i = 0; i < maxNumSubLayersMinus1; i++ ) { 
            if( ptl->sub_layer_profile_present_flag[ i ] ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b, 2, ptl->sub_layer_profile_space[ i ]);
                printf("ptl->sub_layer_profile_space[ i ]: %d ( %ld )\n", ptl->sub_layer_profile_space[ i ], decimal_to_binary( ptl->sub_layer_profile_space[ i ] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_tier_flag[ i ]);
                printf("ptl->sub_layer_tier_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_tier_flag[ i ], decimal_to_binary( ptl->sub_layer_tier_flag[ i ] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b, 5, ptl->sub_layer_profile_idc[ i ]);
                printf("ptl->sub_layer_profile_idc[ i ]: %d ( %ld )\n", ptl->sub_layer_profile_idc[ i ], decimal_to_binary( ptl->sub_layer_profile_idc[ i ] )); 
                for( j = 0; j < 32; j++ ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u(b, 1, ptl->sub_layer_profile_compatibility_flag[ i ][ j ]);
                    printf("ptl->sub_layer_profile_compatibility_flag[ i ][ j ]: %d ( %ld )\n", ptl->sub_layer_profile_compatibility_flag[ i ][ j ], decimal_to_binary( ptl->sub_layer_profile_compatibility_flag[ i ][ j ] )); 
                }
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_progressive_source_flag[ i ]);
                printf("ptl->sub_layer_progressive_source_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_progressive_source_flag[ i ], decimal_to_binary( ptl->sub_layer_progressive_source_flag[ i ] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_interlaced_source_flag[ i ]);
                printf("ptl->sub_layer_interlaced_source_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_interlaced_source_flag[ i ], decimal_to_binary( ptl->sub_layer_interlaced_source_flag[ i ] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_non_packed_constraint_flag[ i ]);
                printf("ptl->sub_layer_non_packed_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_non_packed_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_non_packed_constraint_flag[ i ] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_frame_only_constraint_flag[ i ]);
                printf("ptl->sub_layer_frame_only_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_frame_only_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_frame_only_constraint_flag[ i ] )); 
                if( ptl->sub_layer_profile_idc[ i ] == 4 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 4 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 5 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 5 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 6 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 6 ] ||
                    ptl->sub_layer_profile_idc[ i ] == 7 ||
                    ptl->sub_layer_profile_compatibility_flag[ i ][ 7 ] ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_12bit_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_12bit_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_12bit_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_12bit_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_10bit_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_10bit_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_10bit_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_10bit_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_8bit_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_8bit_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_8bit_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_8bit_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_422chroma_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_422chroma_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_422chroma_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_422chroma_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_420chroma_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_420chroma_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_420chroma_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_420chroma_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_max_monochrome_constraint_flag[ i ]);
                    printf("ptl->sub_layer_max_monochrome_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_max_monochrome_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_max_monochrome_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_intra_constraint_flag[ i ]);
                    printf("ptl->sub_layer_intra_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_intra_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_intra_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_one_picture_only_constraint_flag[ i ]);
                    printf("ptl->sub_layer_one_picture_only_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_one_picture_only_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_one_picture_only_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_lower_bit_rate_constraint_flag[ i ]);
                    printf("ptl->sub_layer_lower_bit_rate_constraint_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_lower_bit_rate_constraint_flag[ i ], decimal_to_binary( ptl->sub_layer_lower_bit_rate_constraint_flag[ i ] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    int sub_layer_reserved_zero_34bits = 34;  bs_write_u(b, sub_layer_reserved_zero_34bits, 0);
                    printf("sub_layer_reserved_zero_34bits: %d ( %ld )\n", sub_layer_reserved_zero_34bits, decimal_to_binary( sub_layer_reserved_zero_34bits )); 
                } else {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    int sub_layer_reserved_zero_43bits = 43;  bs_write_u(b, sub_layer_reserved_zero_43bits, 0);
                    printf("sub_layer_reserved_zero_43bits: %d ( %ld )\n", sub_layer_reserved_zero_43bits, decimal_to_binary( sub_layer_reserved_zero_43bits )); 
                }
            
                if( ( ptl->sub_layer_profile_idc[ i ] >= 1 && ptl->sub_layer_profile_idc[ i ] <= 5 ) ||
                   ptl->sub_layer_profile_compatibility_flag[ 1 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 2 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 3 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 4 ] ||
                   ptl->sub_layer_profile_compatibility_flag[ 5 ] ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, ptl->sub_layer_inbld_flag[ i ]);
                    printf("ptl->sub_layer_inbld_flag[ i ]: %d ( %ld )\n", ptl->sub_layer_inbld_flag[ i ], decimal_to_binary( ptl->sub_layer_inbld_flag[ i ] )); 
                } else {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    int sub_layer_reserved_zero_bit = 1;  bs_write_u(b, sub_layer_reserved_zero_bit, 0);
                    printf("sub_layer_reserved_zero_bit: %d ( %ld )\n", sub_layer_reserved_zero_bit, decimal_to_binary( sub_layer_reserved_zero_bit )); 
                }
            }
            if( ptl->sub_layer_level_present_flag[ i ] ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, ptl->sub_layer_level_idc[ i ]);
                printf("ptl->sub_layer_level_idc[ i ]: %d ( %ld )\n", ptl->sub_layer_level_idc[ i ], decimal_to_binary( ptl->sub_layer_level_idc[ i ] )); 
            }
        }
    }
}

//7.3.4 Scaling list data syntax
void write_debug_hevc_scaling_list_data( hevc_scaling_list_data_t *sld, bs_t* b )
{
    int nextCoef, coefNum;
    for( int sizeId = 0; sizeId < 4; sizeId++ )
        for( int matrixId = 0; matrixId < 6; matrixId += ( sizeId == 3 ) ? 3 : 1 ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ]);
            printf("sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ]: %d ( %ld )\n", sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ], decimal_to_binary( sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] )); 
            if( !sld->scaling_list_pred_mode_flag[ sizeId ][ matrixId ] ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_ue(b, sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ]);
                printf("sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ]: %d ( %ld )\n", sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ], decimal_to_binary( sld->scaling_list_pred_matrix_id_delta[ sizeId ][ matrixId ] )); 
            } else {
                nextCoef = 8;
                coefNum=MIN(64, (1 << (4+(sizeId << 1))));
                if( sizeId > 1 ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_se(b, sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ]);
                    printf("sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ]: %d ( %ld )\n", sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ], decimal_to_binary( sld->scaling_list_dc_coef_minus8[ sizeId - 2 ][ matrixId ] )); 
                }
 
                for( int i = 0; i < coefNum; i++) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_se(b, sld->scaling_list_delta_coef[ sizeId ][ matrixId ]);
                    printf("sld->scaling_list_delta_coef[ sizeId ][ matrixId ]: %d ( %ld )\n", sld->scaling_list_delta_coef[ sizeId ][ matrixId ], decimal_to_binary( sld->scaling_list_delta_coef[ sizeId ][ matrixId ] )); 
                }
            }
 
        }
}

//7.3.6 Slice header syntax
void write_debug_hevc_slice_header(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    if( 0 )
    {
        init_slice_hevc(h);
    }

    hevc_nal_t* nal = h->nal;

    printf("%d.%d: ", b->p - b->start, b->bits_left); 

    bs_write_u1(b, sh->first_slice_segment_in_pic_flag);

    printf("sh->first_slice_segment_in_pic_flag: %d ( %ld )\n", sh->first_slice_segment_in_pic_flag, decimal_to_binary( sh->first_slice_segment_in_pic_flag )); 
    if( nal->nal_unit_type >= HEVC_NAL_UNIT_TYPE_BLA_W_LP && 
       nal->nal_unit_type <= HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->no_output_of_prior_pics_flag);
            printf("sh->no_output_of_prior_pics_flag: %d ( %ld )\n", sh->no_output_of_prior_pics_flag, decimal_to_binary( sh->no_output_of_prior_pics_flag )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_ue(b, sh->pic_parameter_set_id);
    printf("sh->pic_parameter_set_id: %d ( %ld )\n", sh->pic_parameter_set_id, decimal_to_binary( sh->pic_parameter_set_id )); 
    
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];

    //set default value
    sh->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
    sh->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;

    if( !sh->first_slice_segment_in_pic_flag ) {
        if( pps->dependent_slice_segments_enabled_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->dependent_slice_segment_flag);
            printf("sh->dependent_slice_segment_flag: %d ( %ld )\n", sh->dependent_slice_segment_flag, decimal_to_binary( sh->dependent_slice_segment_flag )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b,  getSliceSegmentAddressBitLength( sps ) , sh->slice_segment_address);
        printf("sh->slice_segment_address: %d ( %ld )\n", sh->slice_segment_address, decimal_to_binary( sh->slice_segment_address )); 
    }
    
    if( !sh->dependent_slice_segment_flag ) {
        for( i = 0; i < pps->num_extra_slice_header_bits; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            int slice_reserved_flag = 1;  bs_write_u(b, slice_reserved_flag, 1);
            printf("slice_reserved_flag: %d ( %ld )\n", slice_reserved_flag, decimal_to_binary( slice_reserved_flag )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sh->slice_type);
        printf("sh->slice_type: %d ( %ld )\n", sh->slice_type, decimal_to_binary( sh->slice_type )); 
        if( pps->output_flag_present_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->pic_output_flag);
            printf("sh->pic_output_flag: %d ( %ld )\n", sh->pic_output_flag, decimal_to_binary( sh->pic_output_flag )); 
        }
        if( sps->separate_colour_plane_flag == 1 ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 2, sh->colour_plane_id);
            printf("sh->colour_plane_id: %d ( %ld )\n", sh->colour_plane_id, decimal_to_binary( sh->colour_plane_id )); 
        }
        if( nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_W_RADL &&
            nal->nal_unit_type != HEVC_NAL_UNIT_TYPE_IDR_N_LP) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sh->slice_pic_order_cnt_lsb);
            printf("sh->slice_pic_order_cnt_lsb: %d ( %ld )\n", sh->slice_pic_order_cnt_lsb, decimal_to_binary( sh->slice_pic_order_cnt_lsb )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->short_term_ref_pic_set_sps_flag);
            printf("sh->short_term_ref_pic_set_sps_flag: %d ( %ld )\n", sh->short_term_ref_pic_set_sps_flag, decimal_to_binary( sh->short_term_ref_pic_set_sps_flag )); 
            if( !sh->short_term_ref_pic_set_sps_flag ) {
                write_debug_hevc_st_ref_pic_set( &sh->st_ref_pic_set, b, sps->num_short_term_ref_pic_sets, sps->num_short_term_ref_pic_sets );
            } else if( sps->num_short_term_ref_pic_sets > 1 ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b,  ceil( log2( sps->num_short_term_ref_pic_sets ) ) , sh->short_term_ref_pic_set_idx);
                printf("sh->short_term_ref_pic_set_idx: %d ( %ld )\n", sh->short_term_ref_pic_set_idx, decimal_to_binary( sh->short_term_ref_pic_set_idx )); 
            }
            if( sps->long_term_ref_pics_present_flag ) {
                if( sps->num_long_term_ref_pics_sps > 0 ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_ue(b, sh->num_long_term_sps);
                    printf("sh->num_long_term_sps: %d ( %ld )\n", sh->num_long_term_sps, decimal_to_binary( sh->num_long_term_sps )); 
                }
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_ue(b, sh->num_long_term_pics);
                printf("sh->num_long_term_pics: %d ( %ld )\n", sh->num_long_term_pics, decimal_to_binary( sh->num_long_term_pics )); 
                for( i = 0; i < sh->num_long_term_sps + sh->num_long_term_pics; i++ ) {
                    if( i < sh->num_long_term_sps ) {
                        if( sps->num_long_term_ref_pics_sps > 1 ) {
                            printf("%d.%d: ", b->p - b->start, b->bits_left); 
                            bs_write_u(b,  ceil( log2( sps->num_long_term_ref_pics_sps ) ) , sh->lt_idx_sps[ i ]);
                            printf("sh->lt_idx_sps[ i ]: %d ( %ld )\n", sh->lt_idx_sps[ i ], decimal_to_binary( sh->lt_idx_sps[ i ] )); 
                        }
                    } else {
                        printf("%d.%d: ", b->p - b->start, b->bits_left); 
                        bs_write_u(b,  sps->log2_max_pic_order_cnt_lsb_minus4 + 4 , sh->poc_lsb_lt[ i ]);
                        printf("sh->poc_lsb_lt[ i ]: %d ( %ld )\n", sh->poc_lsb_lt[ i ], decimal_to_binary( sh->poc_lsb_lt[ i ] )); 
                        printf("%d.%d: ", b->p - b->start, b->bits_left); 
                        bs_write_u1(b, sh->used_by_curr_pic_lt_flag[ i ]);
                        printf("sh->used_by_curr_pic_lt_flag[ i ]: %d ( %ld )\n", sh->used_by_curr_pic_lt_flag[ i ], decimal_to_binary( sh->used_by_curr_pic_lt_flag[ i ] )); 
                    }
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, sh->delta_poc_msb_present_flag[ i ]);
                    printf("sh->delta_poc_msb_present_flag[ i ]: %d ( %ld )\n", sh->delta_poc_msb_present_flag[ i ], decimal_to_binary( sh->delta_poc_msb_present_flag[ i ] )); 
                    if( sh->delta_poc_msb_present_flag[ i ]) {
                        printf("%d.%d: ", b->p - b->start, b->bits_left); 
                        bs_write_ue(b, sh->delta_poc_msb_cycle_lt[ i ]);
                        printf("sh->delta_poc_msb_cycle_lt[ i ]: %d ( %ld )\n", sh->delta_poc_msb_cycle_lt[ i ], decimal_to_binary( sh->delta_poc_msb_cycle_lt[ i ] )); 
                    }
                }
            }
            if( sps->sps_temporal_mvp_enabled_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, sh->slice_temporal_mvp_enabled_flag);
                printf("sh->slice_temporal_mvp_enabled_flag: %d ( %ld )\n", sh->slice_temporal_mvp_enabled_flag, decimal_to_binary( sh->slice_temporal_mvp_enabled_flag )); 
            }
        }
        if( sps->sample_adaptive_offset_enabled_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->slice_sao_luma_flag);
            printf("sh->slice_sao_luma_flag: %d ( %ld )\n", sh->slice_sao_luma_flag, decimal_to_binary( sh->slice_sao_luma_flag )); 
            int ChromaArrayType = 0;
            if( sps->separate_colour_plane_flag == 0) {
                ChromaArrayType = sps->chroma_format_idc;
            }
            if( ChromaArrayType != 0 ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, sh->slice_sao_chroma_flag);
                printf("sh->slice_sao_chroma_flag: %d ( %ld )\n", sh->slice_sao_chroma_flag, decimal_to_binary( sh->slice_sao_chroma_flag )); 
            }
        }
        if( sh->slice_type == HEVC_SLICE_TYPE_P || sh->slice_type == HEVC_SLICE_TYPE_B ){
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->num_ref_idx_active_override_flag);
            printf("sh->num_ref_idx_active_override_flag: %d ( %ld )\n", sh->num_ref_idx_active_override_flag, decimal_to_binary( sh->num_ref_idx_active_override_flag )); 
            if( sh->num_ref_idx_active_override_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_ue(b, sh->num_ref_idx_l0_active_minus1);
                printf("sh->num_ref_idx_l0_active_minus1: %d ( %ld )\n", sh->num_ref_idx_l0_active_minus1, decimal_to_binary( sh->num_ref_idx_l0_active_minus1 )); 
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_ue(b, sh->num_ref_idx_l1_active_minus1);
                    printf("sh->num_ref_idx_l1_active_minus1: %d ( %ld )\n", sh->num_ref_idx_l1_active_minus1, decimal_to_binary( sh->num_ref_idx_l1_active_minus1 )); 
                }
            }
            if( pps->lists_modification_present_flag && getNumPicTotalCurr( sps, sh ) > 1 ) {
                write_debug_hevc_ref_pic_lists_modification( h, b );
            }
            
            if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, sh->mvd_l1_zero_flag);
                printf("sh->mvd_l1_zero_flag: %d ( %ld )\n", sh->mvd_l1_zero_flag, decimal_to_binary( sh->mvd_l1_zero_flag )); 
            }
            if( pps->cabac_init_present_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, sh->cabac_init_flag);
                printf("sh->cabac_init_flag: %d ( %ld )\n", sh->cabac_init_flag, decimal_to_binary( sh->cabac_init_flag )); 
            }
            if( sh->slice_temporal_mvp_enabled_flag ) {
                if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_u1(b, sh->collocated_from_l0_flag);
                    printf("sh->collocated_from_l0_flag: %d ( %ld )\n", sh->collocated_from_l0_flag, decimal_to_binary( sh->collocated_from_l0_flag )); 
                }
                if( ( sh->collocated_from_l0_flag && sh->num_ref_idx_l0_active_minus1 > 0 ) ||
                   ( !sh->collocated_from_l0_flag && sh->num_ref_idx_l1_active_minus1 > 0 ) ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_ue(b, sh->collocated_ref_idx);
                    printf("sh->collocated_ref_idx: %d ( %ld )\n", sh->collocated_ref_idx, decimal_to_binary( sh->collocated_ref_idx )); 
                }
            }
            if( ( pps->weighted_pred_flag && sh->slice_type == HEVC_SLICE_TYPE_P ) || 
                ( pps->weighted_bipred_flag && sh->slice_type == HEVC_SLICE_TYPE_B ) ) {
                write_debug_hevc_pred_weight_table( h, b );
            }
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, sh->five_minus_max_num_merge_cand);
            printf("sh->five_minus_max_num_merge_cand: %d ( %ld )\n", sh->five_minus_max_num_merge_cand, decimal_to_binary( sh->five_minus_max_num_merge_cand )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_se(b, sh->slice_qp_delta);
        printf("sh->slice_qp_delta: %d ( %ld )\n", sh->slice_qp_delta, decimal_to_binary( sh->slice_qp_delta )); 
        if( pps->pps_slice_chroma_qp_offsets_present_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, sh->slice_cb_qp_offset);
            printf("sh->slice_cb_qp_offset: %d ( %ld )\n", sh->slice_cb_qp_offset, decimal_to_binary( sh->slice_cb_qp_offset )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, sh->slice_cr_qp_offset);
            printf("sh->slice_cr_qp_offset: %d ( %ld )\n", sh->slice_cr_qp_offset, decimal_to_binary( sh->slice_cr_qp_offset )); 
        }
        if( pps->pps_range_ext.chroma_qp_offset_list_enabled_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->cu_chroma_qp_offset_enabled_flag);
            printf("sh->cu_chroma_qp_offset_enabled_flag: %d ( %ld )\n", sh->cu_chroma_qp_offset_enabled_flag, decimal_to_binary( sh->cu_chroma_qp_offset_enabled_flag )); 
        }
        if( pps->deblocking_filter_override_enabled_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->deblocking_filter_override_flag);
            printf("sh->deblocking_filter_override_flag: %d ( %ld )\n", sh->deblocking_filter_override_flag, decimal_to_binary( sh->deblocking_filter_override_flag )); 
        }
        if( sh->deblocking_filter_override_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->slice_deblocking_filter_disabled_flag);
            printf("sh->slice_deblocking_filter_disabled_flag: %d ( %ld )\n", sh->slice_deblocking_filter_disabled_flag, decimal_to_binary( sh->slice_deblocking_filter_disabled_flag )); 
            if( !sh->slice_deblocking_filter_disabled_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, sh->slice_beta_offset_div2);
                printf("sh->slice_beta_offset_div2: %d ( %ld )\n", sh->slice_beta_offset_div2, decimal_to_binary( sh->slice_beta_offset_div2 )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, sh->slice_tc_offset_div2);
                printf("sh->slice_tc_offset_div2: %d ( %ld )\n", sh->slice_tc_offset_div2, decimal_to_binary( sh->slice_tc_offset_div2 )); 
            }
        }
        if( pps->pps_loop_filter_across_slices_enabled_flag &&
           ( sh->slice_sao_luma_flag || sh->slice_sao_chroma_flag || !sh->slice_deblocking_filter_disabled_flag ) ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, sh->slice_loop_filter_across_slices_enabled_flag);
            printf("sh->slice_loop_filter_across_slices_enabled_flag: %d ( %ld )\n", sh->slice_loop_filter_across_slices_enabled_flag, decimal_to_binary( sh->slice_loop_filter_across_slices_enabled_flag )); 
        }
    }
    if( pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sh->num_entry_point_offsets);
        printf("sh->num_entry_point_offsets: %d ( %ld )\n", sh->num_entry_point_offsets, decimal_to_binary( sh->num_entry_point_offsets )); 
        if( sh->num_entry_point_offsets > 0 ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, sh->offset_len_minus1);
            printf("sh->offset_len_minus1: %d ( %ld )\n", sh->offset_len_minus1, decimal_to_binary( sh->offset_len_minus1 )); 
            for( i = 0; i < sh->num_entry_point_offsets; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b,  sh->offset_len_minus1 + 1 , sh->entry_point_offset_minus1[ i ]);
                printf("sh->entry_point_offset_minus1[ i ]: %d ( %ld )\n", sh->entry_point_offset_minus1[ i ], decimal_to_binary( sh->entry_point_offset_minus1[ i ] )); 
            }
        }
    }
    if( pps->slice_segment_header_extension_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sh->slice_segment_header_extension_length);
        printf("sh->slice_segment_header_extension_length: %d ( %ld )\n", sh->slice_segment_header_extension_length, decimal_to_binary( sh->slice_segment_header_extension_length )); 
        //TODO: support header extension,
        for( i = 0; i < sh->slice_segment_header_extension_length; i++) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            int slice_segment_header_extension_data_byte = 8;  bs_write_u(b, slice_segment_header_extension_data_byte, 0);
            printf("slice_segment_header_extension_data_byte: %d ( %ld )\n", slice_segment_header_extension_data_byte, decimal_to_binary( slice_segment_header_extension_data_byte )); 
        }
    }
    write_debug_hevc_byte_alignment( b );
}

//7.3.6.2 Reference picture list reordering syntax
void write_debug_hevc_ref_pic_lists_modification(hevc_stream_t* h, bs_t* b)
{
    int i;
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    
    bs_write_u1(b, sh->rpld.ref_pic_list_modification_flag_l0);
    
    printf("sh->rpld.ref_pic_list_modification_flag_l0: %d ( %ld )\n", sh->rpld.ref_pic_list_modification_flag_l0, decimal_to_binary( sh->rpld.ref_pic_list_modification_flag_l0 )); 
    if( sh->rpld.ref_pic_list_modification_flag_l0 ) {
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) , sh->rpld.list_entry_l0[ i ]);
            printf("sh->rpld.list_entry_l0[ i ]: %d ( %ld )\n", sh->rpld.list_entry_l0[ i ], decimal_to_binary( sh->rpld.list_entry_l0[ i ] )); 
        }
    }
    
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        // ERROR: value( sh->rpld.ref_pic_list_modification_flag_l1, 1 );
        printf("sh->rpld.ref_pic_list_modification_flag_l1: %d ( %ld )\n", sh->rpld.ref_pic_list_modification_flag_l1, decimal_to_binary( sh->rpld.ref_pic_list_modification_flag_l1 )); 
        if( sh->rpld.ref_pic_list_modification_flag_l1 ) {
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b,  ceil( log2( getNumPicTotalCurr( sps, sh ) ) ) , sh->rpld.list_entry_l1[ i ]);
                printf("sh->rpld.list_entry_l1[ i ]: %d ( %ld )\n", sh->rpld.list_entry_l1[ i ], decimal_to_binary( sh->rpld.list_entry_l1[ i ] )); 
            }
        }
    }
}

//7.3.6.3 Prediction weight table syntax
void write_debug_hevc_pred_weight_table(hevc_stream_t* h, bs_t* b)
{
    hevc_slice_header_t* sh = h->sh;
    hevc_pps_t* pps = &h->pps[sh->pic_parameter_set_id];
    hevc_sps_t* sps = &h->sps[pps->seq_parameter_set_id];
    hevc_pred_weight_table_t *pwt = &sh->pwt;

    int i, j;

    printf("%d.%d: ", b->p - b->start, b->bits_left); 

    bs_write_ue(b, pwt->luma_log2_weight_denom);

    printf("pwt->luma_log2_weight_denom: %d ( %ld )\n", pwt->luma_log2_weight_denom, decimal_to_binary( pwt->luma_log2_weight_denom )); 
    
    int ChromaArrayType = 0;
    if( sps->separate_colour_plane_flag == 0) {
        ChromaArrayType = sps->chroma_format_idc;
    }

    if( ChromaArrayType != 0 ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_se(b, pwt->delta_chroma_log2_weight_denom);
        printf("pwt->delta_chroma_log2_weight_denom: %d ( %ld )\n", pwt->delta_chroma_log2_weight_denom, decimal_to_binary( pwt->delta_chroma_log2_weight_denom )); 
    }
    
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, pwt->luma_weight_l0_flag[i]);
        printf("pwt->luma_weight_l0_flag[i]: %d ( %ld )\n", pwt->luma_weight_l0_flag[i], decimal_to_binary( pwt->luma_weight_l0_flag[i] )); 
    }
    if( ChromaArrayType != 0 )
        for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, pwt->chroma_weight_l0_flag[i]);
            printf("pwt->chroma_weight_l0_flag[i]: %d ( %ld )\n", pwt->chroma_weight_l0_flag[i], decimal_to_binary( pwt->chroma_weight_l0_flag[i] )); 
        }
    for( i = 0; i <= sh->num_ref_idx_l0_active_minus1; i++ ) {
        if( pwt->luma_weight_l0_flag[i] ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pwt->delta_luma_weight_l0[i]);
            printf("pwt->delta_luma_weight_l0[i]: %d ( %ld )\n", pwt->delta_luma_weight_l0[i], decimal_to_binary( pwt->delta_luma_weight_l0[i] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_se(b, pwt->luma_offset_l0[i]);
            printf("pwt->luma_offset_l0[i]: %d ( %ld )\n", pwt->luma_offset_l0[i], decimal_to_binary( pwt->luma_offset_l0[i] )); 
        }
        if( pwt->chroma_weight_l0_flag[i] ) {
            for( j =0; j < 2; j++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, pwt->delta_chroma_weight_l0[i][j]);
                printf("pwt->delta_chroma_weight_l0[i][j]: %d ( %ld )\n", pwt->delta_chroma_weight_l0[i][j], decimal_to_binary( pwt->delta_chroma_weight_l0[i][j] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, pwt->delta_chroma_offset_l0[i][j]);
                printf("pwt->delta_chroma_offset_l0[i][j]: %d ( %ld )\n", pwt->delta_chroma_offset_l0[i][j], decimal_to_binary( pwt->delta_chroma_offset_l0[i][j] )); 
            }
        }
    }
    if( sh->slice_type == HEVC_SLICE_TYPE_B ) {
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, pwt->luma_weight_l1_flag[i]);
            printf("pwt->luma_weight_l1_flag[i]: %d ( %ld )\n", pwt->luma_weight_l1_flag[i], decimal_to_binary( pwt->luma_weight_l1_flag[i] )); 
        }
        if( ChromaArrayType != 0 )
            for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, pwt->chroma_weight_l1_flag[i]);
                printf("pwt->chroma_weight_l1_flag[i]: %d ( %ld )\n", pwt->chroma_weight_l1_flag[i], decimal_to_binary( pwt->chroma_weight_l1_flag[i] )); 
            }
        for( i = 0; i <= sh->num_ref_idx_l1_active_minus1; i++ ) {
            if( pwt->luma_weight_l1_flag[i] ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, pwt->delta_luma_weight_l1[i]);
                printf("pwt->delta_luma_weight_l1[i]: %d ( %ld )\n", pwt->delta_luma_weight_l1[i], decimal_to_binary( pwt->delta_luma_weight_l1[i] )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_se(b, pwt->luma_offset_l1[i]);
                printf("pwt->luma_offset_l1[i]: %d ( %ld )\n", pwt->luma_offset_l1[i], decimal_to_binary( pwt->luma_offset_l1[i] )); 
            }
            if( pwt->chroma_weight_l1_flag[i] ) {
                for( j =0; j < 2; j++ ) {
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_se(b, pwt->delta_chroma_weight_l1[i][j]);
                    printf("pwt->delta_chroma_weight_l1[i][j]: %d ( %ld )\n", pwt->delta_chroma_weight_l1[i][j], decimal_to_binary( pwt->delta_chroma_weight_l1[i][j] )); 
                    printf("%d.%d: ", b->p - b->start, b->bits_left); 
                    bs_write_se(b, pwt->delta_chroma_offset_l1[i][j]);
                    printf("pwt->delta_chroma_offset_l1[i][j]: %d ( %ld )\n", pwt->delta_chroma_offset_l1[i][j], decimal_to_binary( pwt->delta_chroma_offset_l1[i][j] )); 
                }
            }
        }        
    }
}

//7.3.7 Short-term reference picture set syntax
void write_debug_hevc_st_ref_pic_set( hevc_st_ref_pic_set_t *st_ref_pic_set, bs_t* b, int stRpsIdx, int num_short_term_ref_pic_sets )
{
    int i, j;
    
    if( stRpsIdx != 0 ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, st_ref_pic_set->inter_ref_pic_set_prediction_flag);
        printf("st_ref_pic_set->inter_ref_pic_set_prediction_flag: %d ( %ld )\n", st_ref_pic_set->inter_ref_pic_set_prediction_flag, decimal_to_binary( st_ref_pic_set->inter_ref_pic_set_prediction_flag )); 
    }
    if( st_ref_pic_set->inter_ref_pic_set_prediction_flag ) {
        if( stRpsIdx == num_short_term_ref_pic_sets ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, st_ref_pic_set->delta_idx_minus1);
            printf("st_ref_pic_set->delta_idx_minus1: %d ( %ld )\n", st_ref_pic_set->delta_idx_minus1, decimal_to_binary( st_ref_pic_set->delta_idx_minus1 )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, st_ref_pic_set->delta_rps_sign);
        printf("st_ref_pic_set->delta_rps_sign: %d ( %ld )\n", st_ref_pic_set->delta_rps_sign, decimal_to_binary( st_ref_pic_set->delta_rps_sign )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, st_ref_pic_set->abs_delta_rps_minus1);
        printf("st_ref_pic_set->abs_delta_rps_minus1: %d ( %ld )\n", st_ref_pic_set->abs_delta_rps_minus1, decimal_to_binary( st_ref_pic_set->abs_delta_rps_minus1 )); 
        
        int RefRpsIdx = stRpsIdx - ( st_ref_pic_set->delta_idx_minus1 + 1 );
        
        for( j = 0; j <= NumDeltaPocs[ RefRpsIdx ]; j++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_flag[ j ]);
            printf("st_ref_pic_set->used_by_curr_pic_flag[ j ]: %d ( %ld )\n", st_ref_pic_set->used_by_curr_pic_flag[ j ], decimal_to_binary( st_ref_pic_set->used_by_curr_pic_flag[ j ] )); 
            if( !st_ref_pic_set->used_by_curr_pic_flag[ j ] ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, st_ref_pic_set->use_delta_flag[ j ]);
                printf("st_ref_pic_set->use_delta_flag[ j ]: %d ( %ld )\n", st_ref_pic_set->use_delta_flag[ j ], decimal_to_binary( st_ref_pic_set->use_delta_flag[ j ] )); 
            }
        }
    } else {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, st_ref_pic_set->num_negative_pics);
        printf("st_ref_pic_set->num_negative_pics: %d ( %ld )\n", st_ref_pic_set->num_negative_pics, decimal_to_binary( st_ref_pic_set->num_negative_pics )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, st_ref_pic_set->num_positive_pics);
        printf("st_ref_pic_set->num_positive_pics: %d ( %ld )\n", st_ref_pic_set->num_positive_pics, decimal_to_binary( st_ref_pic_set->num_positive_pics )); 
        for( i = 0; i < st_ref_pic_set->num_negative_pics; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, st_ref_pic_set->delta_poc_s0_minus1[ i ]);
            printf("st_ref_pic_set->delta_poc_s0_minus1[ i ]: %d ( %ld )\n", st_ref_pic_set->delta_poc_s0_minus1[ i ], decimal_to_binary( st_ref_pic_set->delta_poc_s0_minus1[ i ] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_s0_flag[ i ]);
            printf("st_ref_pic_set->used_by_curr_pic_s0_flag[ i ]: %d ( %ld )\n", st_ref_pic_set->used_by_curr_pic_s0_flag[ i ], decimal_to_binary( st_ref_pic_set->used_by_curr_pic_s0_flag[ i ] )); 
            
            //update derived field
            UsedByCurrPicS0[ stRpsIdx ][ i ] = st_ref_pic_set->used_by_curr_pic_s0_flag[ i ];
            
            if( i == 0 ) {
                DeltaPocS0[ stRpsIdx ][ i ] = -1 * ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            } else {
                DeltaPocS0[ stRpsIdx ][ i ] = DeltaPocS0[ stRpsIdx ][ i - 1 ] - ( st_ref_pic_set->delta_poc_s0_minus1[ i ] + 1 );
            }
        }
        for( i = 0; i < st_ref_pic_set->num_positive_pics; i++ ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, st_ref_pic_set->delta_poc_s1_minus1[ i ]);
            printf("st_ref_pic_set->delta_poc_s1_minus1[ i ]: %d ( %ld )\n", st_ref_pic_set->delta_poc_s1_minus1[ i ], decimal_to_binary( st_ref_pic_set->delta_poc_s1_minus1[ i ] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, st_ref_pic_set->used_by_curr_pic_s1_flag[ i ]);
            printf("st_ref_pic_set->used_by_curr_pic_s1_flag[ i ]: %d ( %ld )\n", st_ref_pic_set->used_by_curr_pic_s1_flag[ i ], decimal_to_binary( st_ref_pic_set->used_by_curr_pic_s1_flag[ i ] )); 
        
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
void write_debug_hevc_vui_parameters(hevc_sps_t* sps, bs_t* b)
{
    hevc_vui_t* vui = &sps->vui;
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->aspect_ratio_info_present_flag);
    printf("vui->aspect_ratio_info_present_flag: %d ( %ld )\n", vui->aspect_ratio_info_present_flag, decimal_to_binary( vui->aspect_ratio_info_present_flag )); 
    if( vui->aspect_ratio_info_present_flag )
    {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u8(b, vui->aspect_ratio_idc);
        printf("vui->aspect_ratio_idc: %d ( %ld )\n", vui->aspect_ratio_idc, decimal_to_binary( vui->aspect_ratio_idc )); 
        if( vui->aspect_ratio_idc == SAR_Extended )
        {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 16, vui->sar_width);
            printf("vui->sar_width: %d ( %ld )\n", vui->sar_width, decimal_to_binary( vui->sar_width )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 16, vui->sar_height);
            printf("vui->sar_height: %d ( %ld )\n", vui->sar_height, decimal_to_binary( vui->sar_height )); 
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->overscan_info_present_flag);
    printf("vui->overscan_info_present_flag: %d ( %ld )\n", vui->overscan_info_present_flag, decimal_to_binary( vui->overscan_info_present_flag )); 
    if( vui->overscan_info_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->overscan_appropriate_flag);
        printf("vui->overscan_appropriate_flag: %d ( %ld )\n", vui->overscan_appropriate_flag, decimal_to_binary( vui->overscan_appropriate_flag )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->video_signal_type_present_flag);
    printf("vui->video_signal_type_present_flag: %d ( %ld )\n", vui->video_signal_type_present_flag, decimal_to_binary( vui->video_signal_type_present_flag )); 
    if( vui->video_signal_type_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 3, vui->video_format);
        printf("vui->video_format: %d ( %ld )\n", vui->video_format, decimal_to_binary( vui->video_format )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->video_full_range_flag);
        printf("vui->video_full_range_flag: %d ( %ld )\n", vui->video_full_range_flag, decimal_to_binary( vui->video_full_range_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->colour_description_present_flag);
        printf("vui->colour_description_present_flag: %d ( %ld )\n", vui->colour_description_present_flag, decimal_to_binary( vui->colour_description_present_flag )); 
        if( vui->colour_description_present_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u8(b, vui->colour_primaries);
            printf("vui->colour_primaries: %d ( %ld )\n", vui->colour_primaries, decimal_to_binary( vui->colour_primaries )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u8(b, vui->transfer_characteristics);
            printf("vui->transfer_characteristics: %d ( %ld )\n", vui->transfer_characteristics, decimal_to_binary( vui->transfer_characteristics )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u8(b, vui->matrix_coefficients);
            printf("vui->matrix_coefficients: %d ( %ld )\n", vui->matrix_coefficients, decimal_to_binary( vui->matrix_coefficients )); 
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->chroma_loc_info_present_flag);
    printf("vui->chroma_loc_info_present_flag: %d ( %ld )\n", vui->chroma_loc_info_present_flag, decimal_to_binary( vui->chroma_loc_info_present_flag )); 
    if( vui->chroma_loc_info_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->chroma_sample_loc_type_top_field);
        printf("vui->chroma_sample_loc_type_top_field: %d ( %ld )\n", vui->chroma_sample_loc_type_top_field, decimal_to_binary( vui->chroma_sample_loc_type_top_field )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->chroma_sample_loc_type_bottom_field);
        printf("vui->chroma_sample_loc_type_bottom_field: %d ( %ld )\n", vui->chroma_sample_loc_type_bottom_field, decimal_to_binary( vui->chroma_sample_loc_type_bottom_field )); 
    }
    
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    
    bs_write_u1(b, vui->neutral_chroma_indication_flag);
    
    printf("vui->neutral_chroma_indication_flag: %d ( %ld )\n", vui->neutral_chroma_indication_flag, decimal_to_binary( vui->neutral_chroma_indication_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->field_seq_flag);
    printf("vui->field_seq_flag: %d ( %ld )\n", vui->field_seq_flag, decimal_to_binary( vui->field_seq_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->frame_field_info_present_flag);
    printf("vui->frame_field_info_present_flag: %d ( %ld )\n", vui->frame_field_info_present_flag, decimal_to_binary( vui->frame_field_info_present_flag )); 
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->default_display_window_flag);
    printf("vui->default_display_window_flag: %d ( %ld )\n", vui->default_display_window_flag, decimal_to_binary( vui->default_display_window_flag )); 
    if( vui->default_display_window_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->def_disp_win_left_offset);
        printf("vui->def_disp_win_left_offset: %d ( %ld )\n", vui->def_disp_win_left_offset, decimal_to_binary( vui->def_disp_win_left_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->def_disp_win_right_offset);
        printf("vui->def_disp_win_right_offset: %d ( %ld )\n", vui->def_disp_win_right_offset, decimal_to_binary( vui->def_disp_win_right_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->def_disp_win_top_offset);
        printf("vui->def_disp_win_top_offset: %d ( %ld )\n", vui->def_disp_win_top_offset, decimal_to_binary( vui->def_disp_win_top_offset )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->def_disp_win_bottom_offset);
        printf("vui->def_disp_win_bottom_offset: %d ( %ld )\n", vui->def_disp_win_bottom_offset, decimal_to_binary( vui->def_disp_win_bottom_offset )); 
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->vui_timing_info_present_flag);
    printf("vui->vui_timing_info_present_flag: %d ( %ld )\n", vui->vui_timing_info_present_flag, decimal_to_binary( vui->vui_timing_info_present_flag )); 
    if( vui->vui_timing_info_present_flag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 32, vui->vui_num_units_in_tick);
        printf("vui->vui_num_units_in_tick: %d ( %ld )\n", vui->vui_num_units_in_tick, decimal_to_binary( vui->vui_num_units_in_tick )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u(b, 32, vui->vui_time_scale);
        printf("vui->vui_time_scale: %d ( %ld )\n", vui->vui_time_scale, decimal_to_binary( vui->vui_time_scale )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->vui_poc_proportional_to_timing_flag);
        printf("vui->vui_poc_proportional_to_timing_flag: %d ( %ld )\n", vui->vui_poc_proportional_to_timing_flag, decimal_to_binary( vui->vui_poc_proportional_to_timing_flag )); 
        if( vui->vui_poc_proportional_to_timing_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, vui->vui_num_ticks_poc_diff_one_minus1);
            printf("vui->vui_num_ticks_poc_diff_one_minus1: %d ( %ld )\n", vui->vui_num_ticks_poc_diff_one_minus1, decimal_to_binary( vui->vui_num_ticks_poc_diff_one_minus1 )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->vui_hrd_parameters_present_flag);
        printf("vui->vui_hrd_parameters_present_flag: %d ( %ld )\n", vui->vui_hrd_parameters_present_flag, decimal_to_binary( vui->vui_hrd_parameters_present_flag )); 
        if( vui->vui_hrd_parameters_present_flag ) {
            write_debug_hevc_hrd_parameters( &vui->hrd, b, 1, sps->sps_max_sub_layers_minus1 );
        }
    }
    printf("%d.%d: ", b->p - b->start, b->bits_left); 
    bs_write_u1(b, vui->bitstream_restriction_flag);
    printf("vui->bitstream_restriction_flag: %d ( %ld )\n", vui->bitstream_restriction_flag, decimal_to_binary( vui->bitstream_restriction_flag )); 
    if( vui->bitstream_restriction_flag )
    {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->tiles_fixed_structure_flag);
        printf("vui->tiles_fixed_structure_flag: %d ( %ld )\n", vui->tiles_fixed_structure_flag, decimal_to_binary( vui->tiles_fixed_structure_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->motion_vectors_over_pic_boundaries_flag);
        printf("vui->motion_vectors_over_pic_boundaries_flag: %d ( %ld )\n", vui->motion_vectors_over_pic_boundaries_flag, decimal_to_binary( vui->motion_vectors_over_pic_boundaries_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, vui->restricted_ref_pic_lists_flag);
        printf("vui->restricted_ref_pic_lists_flag: %d ( %ld )\n", vui->restricted_ref_pic_lists_flag, decimal_to_binary( vui->restricted_ref_pic_lists_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->min_spatial_segmentation_idc);
        printf("vui->min_spatial_segmentation_idc: %d ( %ld )\n", vui->min_spatial_segmentation_idc, decimal_to_binary( vui->min_spatial_segmentation_idc )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->max_bytes_per_pic_denom);
        printf("vui->max_bytes_per_pic_denom: %d ( %ld )\n", vui->max_bytes_per_pic_denom, decimal_to_binary( vui->max_bytes_per_pic_denom )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->max_bits_per_min_cu_denom);
        printf("vui->max_bits_per_min_cu_denom: %d ( %ld )\n", vui->max_bits_per_min_cu_denom, decimal_to_binary( vui->max_bits_per_min_cu_denom )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->log2_max_mv_length_horizontal);
        printf("vui->log2_max_mv_length_horizontal: %d ( %ld )\n", vui->log2_max_mv_length_horizontal, decimal_to_binary( vui->log2_max_mv_length_horizontal )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, vui->log2_max_mv_length_vertical);
        printf("vui->log2_max_mv_length_vertical: %d ( %ld )\n", vui->log2_max_mv_length_vertical, decimal_to_binary( vui->log2_max_mv_length_vertical )); 
    }
}

//Appendix E.2.2 HRD parameters syntax
void write_debug_hevc_hrd_parameters(hevc_hrd_t* hrd, bs_t* b, int commonInfPresentFlag, int maxNumSubLayersMinus1)
{
    if( commonInfPresentFlag ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, hrd->nal_hrd_parameters_present_flag);
        printf("hrd->nal_hrd_parameters_present_flag: %d ( %ld )\n", hrd->nal_hrd_parameters_present_flag, decimal_to_binary( hrd->nal_hrd_parameters_present_flag )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, hrd->vcl_hrd_parameters_present_flag);
        printf("hrd->vcl_hrd_parameters_present_flag: %d ( %ld )\n", hrd->vcl_hrd_parameters_present_flag, decimal_to_binary( hrd->vcl_hrd_parameters_present_flag )); 
        if( hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag ){
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, hrd->sub_pic_hrd_params_present_flag);
            printf("hrd->sub_pic_hrd_params_present_flag: %d ( %ld )\n", hrd->sub_pic_hrd_params_present_flag, decimal_to_binary( hrd->sub_pic_hrd_params_present_flag )); 
            if( hrd->sub_pic_hrd_params_present_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u8(b, hrd->tick_divisor_minus2);
                printf("hrd->tick_divisor_minus2: %d ( %ld )\n", hrd->tick_divisor_minus2, decimal_to_binary( hrd->tick_divisor_minus2 )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b, 5, hrd->du_cpb_removal_delay_increment_length_minus1);
                printf("hrd->du_cpb_removal_delay_increment_length_minus1: %d ( %ld )\n", hrd->du_cpb_removal_delay_increment_length_minus1, decimal_to_binary( hrd->du_cpb_removal_delay_increment_length_minus1 )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u1(b, hrd->sub_pic_cpb_params_in_pic_timing_sei_flag);
                printf("hrd->sub_pic_cpb_params_in_pic_timing_sei_flag: %d ( %ld )\n", hrd->sub_pic_cpb_params_in_pic_timing_sei_flag, decimal_to_binary( hrd->sub_pic_cpb_params_in_pic_timing_sei_flag )); 
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b, 5, hrd->dpb_output_delay_du_length_minus1);
                printf("hrd->dpb_output_delay_du_length_minus1: %d ( %ld )\n", hrd->dpb_output_delay_du_length_minus1, decimal_to_binary( hrd->dpb_output_delay_du_length_minus1 )); 
            }
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 4, hrd->bit_rate_scale);
            printf("hrd->bit_rate_scale: %d ( %ld )\n", hrd->bit_rate_scale, decimal_to_binary( hrd->bit_rate_scale )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 4, hrd->cpb_size_scale);
            printf("hrd->cpb_size_scale: %d ( %ld )\n", hrd->cpb_size_scale, decimal_to_binary( hrd->cpb_size_scale )); 
            if( hrd->sub_pic_hrd_params_present_flag ) {
                printf("%d.%d: ", b->p - b->start, b->bits_left); 
                bs_write_u(b, 4, hrd->cpb_size_du_scale);
                printf("hrd->cpb_size_du_scale: %d ( %ld )\n", hrd->cpb_size_du_scale, decimal_to_binary( hrd->cpb_size_du_scale )); 
            }
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 5, hrd->initial_cpb_removal_delay_length_minus1);
            printf("hrd->initial_cpb_removal_delay_length_minus1: %d ( %ld )\n", hrd->initial_cpb_removal_delay_length_minus1, decimal_to_binary( hrd->initial_cpb_removal_delay_length_minus1 )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 5, hrd->au_cpb_removal_delay_length_minus1);
            printf("hrd->au_cpb_removal_delay_length_minus1: %d ( %ld )\n", hrd->au_cpb_removal_delay_length_minus1, decimal_to_binary( hrd->au_cpb_removal_delay_length_minus1 )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u(b, 5, hrd->dpb_output_delay_length_minus1);
            printf("hrd->dpb_output_delay_length_minus1: %d ( %ld )\n", hrd->dpb_output_delay_length_minus1, decimal_to_binary( hrd->dpb_output_delay_length_minus1 )); 
        }
    }
    
    for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, hrd->fixed_pic_rate_general_flag[ i ]);
        printf("hrd->fixed_pic_rate_general_flag[ i ]: %d ( %ld )\n", hrd->fixed_pic_rate_general_flag[ i ], decimal_to_binary( hrd->fixed_pic_rate_general_flag[ i ] )); 
        if( !hrd->fixed_pic_rate_general_flag[ i ] ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, hrd->fixed_pic_rate_within_cvs_flag[ i ]);
            printf("hrd->fixed_pic_rate_within_cvs_flag[ i ]: %d ( %ld )\n", hrd->fixed_pic_rate_within_cvs_flag[ i ], decimal_to_binary( hrd->fixed_pic_rate_within_cvs_flag[ i ] )); 
        }
        if( hrd->fixed_pic_rate_within_cvs_flag[ i ] ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, hrd->elemental_duration_in_tc_minus1[ i ]);
            printf("hrd->elemental_duration_in_tc_minus1[ i ]: %d ( %ld )\n", hrd->elemental_duration_in_tc_minus1[ i ], decimal_to_binary( hrd->elemental_duration_in_tc_minus1[ i ] )); 
        } else {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_u1(b, hrd->low_delay_hrd_flag[ i ]);
            printf("hrd->low_delay_hrd_flag[ i ]: %d ( %ld )\n", hrd->low_delay_hrd_flag[ i ], decimal_to_binary( hrd->low_delay_hrd_flag[ i ] )); 
        }
        if( hrd->low_delay_hrd_flag[ i ] ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, hrd->cpb_cnt_minus1[ i ]);
            printf("hrd->cpb_cnt_minus1[ i ]: %d ( %ld )\n", hrd->cpb_cnt_minus1[ i ], decimal_to_binary( hrd->cpb_cnt_minus1[ i ] )); 
        }
        if( hrd->nal_hrd_parameters_present_flag ) {
            write_debug_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_nal[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
        if( hrd->vcl_hrd_parameters_present_flag ) {
            write_debug_hevc_sub_layer_hrd_parameters(&hrd->sub_layer_hrd_vcl[i], b, hrd->cpb_cnt_minus1[ i ] + 1, hrd->sub_pic_hrd_params_present_flag);
        }
    }
}

//Appendix E.2.3 Sub-layer HRD parameters syntax
void write_debug_hevc_sub_layer_hrd_parameters(hevc_sub_layer_hrd_t* sub_layer_hrd, bs_t* b, int CpbCnt, int sub_pic_hrd_params_present_flag)
{
    for( int i = 0; i <= CpbCnt; i++ ) {
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sub_layer_hrd->bit_rate_value_minus1[i]);
        printf("sub_layer_hrd->bit_rate_value_minus1[i]: %d ( %ld )\n", sub_layer_hrd->bit_rate_value_minus1[i], decimal_to_binary( sub_layer_hrd->bit_rate_value_minus1[i] )); 
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_ue(b, sub_layer_hrd->cpb_size_value_minus1[i]);
        printf("sub_layer_hrd->cpb_size_value_minus1[i]: %d ( %ld )\n", sub_layer_hrd->cpb_size_value_minus1[i], decimal_to_binary( sub_layer_hrd->cpb_size_value_minus1[i] )); 
        if( sub_pic_hrd_params_present_flag ) {
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, sub_layer_hrd->cpb_size_du_value_minus1[i]);
            printf("sub_layer_hrd->cpb_size_du_value_minus1[i]: %d ( %ld )\n", sub_layer_hrd->cpb_size_du_value_minus1[i], decimal_to_binary( sub_layer_hrd->cpb_size_du_value_minus1[i] )); 
            printf("%d.%d: ", b->p - b->start, b->bits_left); 
            bs_write_ue(b, sub_layer_hrd->bit_rate_du_value_minus1[i]);
            printf("sub_layer_hrd->bit_rate_du_value_minus1[i]: %d ( %ld )\n", sub_layer_hrd->bit_rate_du_value_minus1[i], decimal_to_binary( sub_layer_hrd->bit_rate_du_value_minus1[i] )); 
        }
        printf("%d.%d: ", b->p - b->start, b->bits_left); 
        bs_write_u1(b, sub_layer_hrd->cbr_flag[i]);
        printf("sub_layer_hrd->cbr_flag[i]: %d ( %ld )\n", sub_layer_hrd->cbr_flag[i], decimal_to_binary( sub_layer_hrd->cbr_flag[i] )); 
    }
}
