/*
  hevc_stream.h
  Created by leslie_qiwa@gmail.com on 6/8/17
 */


#ifndef _HEVC_STREAM_H
#define _HEVC_STREAM_H        1

#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "bs.h"
#include "h264_sei.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_NUM_SUBLAYERS            32
#define MAX_NUM_HRD_PARAM            10
#define MAX_CPB_CNT                  32
#define MAX_NUM_NEGATIVE_PICS        32
#define MAX_NUM_POSITIVE_PICS        32
#define MAX_NUM_REF_PICS_L0          32
#define MAX_NUM_REF_PICS_L1          32
#define MAX_NUM_SHORT_TERM_REF_PICS  32
#define MAX_NUM_LONG_TERM_REF_PICS   32
#define MAX_NUM_LONG_TERM_REF_PICS   32
#define MAX_NUM_PALLETTE_PREDICTOR   32
#define MAX_NUM_CHROMA_QP_OFFSET_LST 32
#define MAX_NUM_ENTRY_POINT_OFFSET   32
#define MAX_NUM_TILE_COLUMN          32
#define MAX_NUM_TILE_ROW             32

/**
   hrd_parameters
   @see E.2.2 HRD parameters syntax
*/
typedef struct
{
    int bit_rate_value_minus1[MAX_CPB_CNT];
    int cpb_size_value_minus1[MAX_CPB_CNT];
    int cpb_size_du_value_minus1[MAX_CPB_CNT];
    int bit_rate_du_value_minus1[MAX_CPB_CNT];
    int cbr_flag[MAX_CPB_CNT];
} hevc_sub_layer_hrd_t;
/**
   hrd_parameters
   @see E.2.2 HRD parameters syntax
*/
typedef struct
{
    int nal_hrd_parameters_present_flag;
    int vcl_hrd_parameters_present_flag;
    int sub_pic_hrd_params_present_flag;
    int tick_divisor_minus2;
    int du_cpb_removal_delay_increment_length_minus1;
    int sub_pic_cpb_params_in_pic_timing_sei_flag;
    int dpb_output_delay_du_length_minus1;
    int bit_rate_scale;
    int cpb_size_scale;
    int cpb_size_du_scale;
    int initial_cpb_removal_delay_length_minus1;
    int au_cpb_removal_delay_length_minus1;
    int dpb_output_delay_length_minus1;
    int fixed_pic_rate_general_flag[MAX_NUM_SUBLAYERS];
    int fixed_pic_rate_within_cvs_flag[MAX_NUM_SUBLAYERS];
    int elemental_duration_in_tc_minus1[MAX_NUM_SUBLAYERS];
    int low_delay_hrd_flag[MAX_NUM_SUBLAYERS];
    int cpb_cnt_minus1[MAX_NUM_SUBLAYERS];
    hevc_sub_layer_hrd_t sub_layer_hrd_nal[MAX_NUM_SUBLAYERS];
    hevc_sub_layer_hrd_t sub_layer_hrd_vcl[MAX_NUM_SUBLAYERS];
} hevc_hrd_t;

/**
   Profile, tier and level
   @see 7.3 Profile, tier and level syntax
*/
typedef struct
{
    //profile parameters
    int general_profile_space;
    int general_tier_flag;
    int general_profile_idc;
    int general_profile_compatibility_flag[32];
    int general_progressive_source_flag;
    int general_interlaced_source_flag;
    int general_non_packed_constraint_flag;
    int general_frame_only_constraint_flag;
    int general_max_12bit_constraint_flag;
    int general_max_10bit_constraint_flag;
    int general_max_8bit_constraint_flag;
    int general_max_422chroma_constraint_flag;
    int general_max_420chroma_constraint_flag;
    int general_max_monochrome_constraint_flag;
    int general_intra_constraint_flag;
    int general_one_picture_only_constraint_flag;
    int general_lower_bit_rate_constraint_flag;
    int general_max_14bit_constraint_flag;
    int general_inbld_flag;
    // level parameters
    int general_level_idc;
    int sub_layer_profile_present_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_level_present_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_profile_space[MAX_NUM_SUBLAYERS];
    int sub_layer_tier_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_profile_idc[MAX_NUM_SUBLAYERS];
    int sub_layer_profile_compatibility_flag[MAX_NUM_SUBLAYERS][32];
    int sub_layer_progressive_source_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_interlaced_source_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_non_packed_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_frame_only_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_12bit_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_10bit_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_8bit_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_422chroma_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_420chroma_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_monochrome_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_intra_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_one_picture_only_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_lower_bit_rate_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_max_14bit_constraint_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_inbld_flag[MAX_NUM_SUBLAYERS];
    int sub_layer_level_idc[MAX_NUM_SUBLAYERS];
} hevc_profile_tier_level_t;

/**
   Scaling List Data
   @see 7.3.4 scaling list data syntax
*/
typedef struct
{
    int scaling_list_pred_mode_flag[4][6];
    int scaling_list_pred_matrix_id_delta[4][6];
    int scaling_list_dc_coef_minus8[2][6];
    int scaling_list_delta_coef[4][64];
} hevc_scaling_list_data_t;

/**
   Video Parameter Set
   @see 7.3.2.1 Video parameter set RBSP syntax
   @see read_hevc_video_parameter_set_rbsp
   @see write_hevc_video_parameter_set_rbsp
   @see debug_vps
*/
typedef struct
{
    int vps_video_parameter_set_id;
    int vps_base_layer_internal_flag;
    int vps_base_layer_available_flag;
    int vps_max_layers_minus1;
    int vps_max_sub_layers_minus1;
    int vps_temporal_id_nesting_flag;

    hevc_profile_tier_level_t ptl;

    int vps_sub_layer_ordering_info_present_flag;
    int vps_max_dec_pic_buffering_minus1[MAX_NUM_SUBLAYERS];
    int vps_max_num_reorder_pics[MAX_NUM_SUBLAYERS];
    int vps_max_latency_increase_plus1[MAX_NUM_SUBLAYERS];
    int vps_max_layer_id;
    int vps_num_layer_sets_minus1;
    int layer_id_included_flag[MAX_NUM_SUBLAYERS][MAX_NUM_SUBLAYERS]; /*use same fixed value*/
    int vps_timing_info_present_flag;
    int vps_num_units_in_tick;
    int vps_time_scale;
    int vps_poc_proportional_to_timing_flag;
    int vps_num_ticks_poc_diff_one_minus1;
    int vps_num_hrd_parameters;
    int hrd_layer_set_idx[MAX_NUM_HRD_PARAM];
    int cprms_present_flag[MAX_NUM_HRD_PARAM];
    hevc_hrd_t hrd[MAX_NUM_HRD_PARAM];
    int vps_extension_flag;
    int vps_extension_data_flag;
} hevc_vps_t;

/**
   Sequence Parameter Set
   @see 7.3.2.2 Sequence parameter set RBSP syntax
*/
typedef struct
{
    int inter_ref_pic_set_prediction_flag;
    int delta_idx_minus1;
    int delta_rps_sign;
    int abs_delta_rps_minus1;
    int used_by_curr_pic_flag[MAX_NUM_SHORT_TERM_REF_PICS];
    int use_delta_flag[MAX_NUM_SHORT_TERM_REF_PICS];
    int num_negative_pics;
    int num_positive_pics;
    int delta_poc_s0_minus1[MAX_NUM_NEGATIVE_PICS];
    int used_by_curr_pic_s0_flag[MAX_NUM_NEGATIVE_PICS];
    int delta_poc_s1_minus1[MAX_NUM_POSITIVE_PICS];
    int used_by_curr_pic_s1_flag[MAX_NUM_NEGATIVE_PICS];
} hevc_st_ref_pic_set_t;

/**
   Sequence Parameter Set
   @see 7.3.2.2 Sequence parameter set RBSP syntax
*/
typedef struct
{
    int aspect_ratio_info_present_flag;
    int aspect_ratio_idc;
    int sar_width;
    int sar_height;
    int overscan_info_present_flag;
    int overscan_appropriate_flag;
    int video_signal_type_present_flag;
    int video_format;
    int video_full_range_flag;
    int colour_description_present_flag;
    int colour_primaries;
    int transfer_characteristics;
    int matrix_coefficients;
    int chroma_loc_info_present_flag;
    int chroma_sample_loc_type_top_field;
    int chroma_sample_loc_type_bottom_field;
    int neutral_chroma_indication_flag;
    int field_seq_flag;
    int frame_field_info_present_flag;
    int default_display_window_flag; 
    int def_disp_win_left_offset;
    int def_disp_win_right_offset;
    int def_disp_win_top_offset;
    int def_disp_win_bottom_offset;
    int vui_timing_info_present_flag;
    int vui_num_units_in_tick;
    int vui_time_scale;
    int vui_poc_proportional_to_timing_flag;
    int vui_num_ticks_poc_diff_one_minus1;
    int vui_hrd_parameters_present_flag;
    hevc_hrd_t hrd;
    int bitstream_restriction_flag;
    int tiles_fixed_structure_flag;
    int motion_vectors_over_pic_boundaries_flag;
    int restricted_ref_pic_lists_flag;
    int min_spatial_segmentation_idc;
    int max_bytes_per_pic_denom;
    int max_bits_per_min_cu_denom;
    int log2_max_mv_length_horizontal;
    int log2_max_mv_length_vertical;
} hevc_vui_t;

/**
   Sequence Parameter Set range extension syntax
   @see 7.3.2.2.2 Sequence parameter set range extension syntax
*/
typedef struct
{
    int transform_skip_rotation_enabled_flag;
    int transform_skip_context_enabled_flag;
    int implicit_rdpcm_enabled_flag;
    int explicit_rdpcm_enabled_flag;
    int extended_precision_processing_flag;
    int intra_smoothing_disabled_flag;
    int high_precision_offsets_enabled_flag;
    int persistent_rice_adaptation_enabled_flag;
    int cabac_bypass_alignment_enabled_flag;
} hevc_sps_range_ext_t;

/**
   Sequence Parameter Set screen content coding extension syntax
   @see 7.3.2.2.3 Sequence parameter set screen content coding extension syntax
*/
typedef struct
{
    int sps_curr_pic_ref_enabled_flag;
    int palette_mode_enabled_flag;
    int palette_max_size;
    int delta_palette_max_predictor_size;
    int sps_palette_predictor_initializer_present_flag;
    int sps_num_palette_predictor_initializer_minus1;
    int sps_palette_predictor_initializers[3][MAX_NUM_PALLETTE_PREDICTOR];
    int motion_vector_resolution_control_idc;
    int intra_boundary_filtering_disabled_flag;
} hevc_sps_scc_ext_t;

/**
   Sequence Parameter Set
   @see 7.3.2.2 Sequence parameter set RBSP syntax
   @see read_hevc_seq_parameter_set_rbsp
   @see write_hevc_seq_parameter_set_rbsp
   @see debug_sps
*/
typedef struct
{
    int sps_video_parameter_set_id;
    int sps_max_sub_layers_minus1;
    int sps_temporal_id_nesting_flag;
    hevc_profile_tier_level_t ptl;
    int sps_seq_parameter_set_id;
    int chroma_format_idc;
    int separate_colour_plane_flag;
    int pic_width_in_luma_samples;
    int pic_height_in_luma_samples;
    int conformance_window_flag;
    int conf_win_left_offset;
    int conf_win_right_offset;
    int conf_win_top_offset;
    int conf_win_bottom_offset;
    int bit_depth_luma_minus8;
    int bit_depth_chroma_minus8;
    int log2_max_pic_order_cnt_lsb_minus4;
    int sps_sub_layer_ordering_info_present_flag;
    int sps_max_dec_pic_buffering_minus1[MAX_NUM_SUBLAYERS];
    int sps_max_num_reorder_pics[MAX_NUM_SUBLAYERS];
    int sps_max_latency_increase_plus1[MAX_NUM_SUBLAYERS];
    int log2_min_luma_coding_block_size_minus3;
    int log2_diff_max_min_luma_coding_block_size;
    int log2_min_luma_transform_block_size_minus2;
    int log2_diff_max_min_luma_transform_block_size;
    int max_transform_hierarchy_depth_inter;
    int max_transform_hierarchy_depth_intra;
    int scaling_list_enabled_flag;
    int sps_scaling_list_data_present_flag;
    hevc_scaling_list_data_t scaling_list_data;
    int amp_enabled_flag;
    int sample_adaptive_offset_enabled_flag;
    int pcm_enabled_flag;
    int pcm_sample_bit_depth_luma_minus1;
    int pcm_sample_bit_depth_chroma_minus1;
    int log2_min_pcm_luma_coding_block_size_minus3;
    int log2_diff_max_min_pcm_luma_coding_block_size;
    int pcm_loop_filter_disabled_flag;
    int num_short_term_ref_pic_sets;
    hevc_st_ref_pic_set_t st_ref_pic_set[MAX_NUM_SHORT_TERM_REF_PICS];
    int long_term_ref_pics_present_flag;
    int num_long_term_ref_pics_sps;
    int lt_ref_pic_poc_lsb_sps[MAX_NUM_LONG_TERM_REF_PICS];
    int used_by_curr_pic_lt_sps_flag[MAX_NUM_LONG_TERM_REF_PICS];
    int sps_temporal_mvp_enabled_flag;
    int strong_intra_smoothing_enabled_flag;
    int vui_parameters_present_flag;
    hevc_vui_t vui;
    int sps_extension_present_flag;
    int sps_range_extension_flag;
    int sps_multilayer_extension_flag;
    int sps_3d_extension_flag;
    int sps_extension_5bits;
    hevc_sps_range_ext_t sps_range_ext;
    //TODO: support sps_extension_data_flag;
    //TODO: support SVC/MVC extensions
} hevc_sps_t;

/**
   Picture Parameter Set range extension syntax
   @see 7.3.2.3.2 Picture parameter set range extension syntax
*/
typedef struct
{
    int log2_max_transform_skip_block_size_minus2;
    int cross_component_prediction_enabled_flag;
    int chroma_qp_offset_list_enabled_flag;
    int diff_cu_chroma_qp_offset_depth;
    int chroma_qp_offset_list_len_minus1;
    int cb_qp_offset_list[MAX_NUM_CHROMA_QP_OFFSET_LST];
    int cr_qp_offset_list[MAX_NUM_CHROMA_QP_OFFSET_LST];
    int log2_sao_offset_scale_luma;
    int log2_sao_offset_scale_chroma;
} hevc_pps_range_ext_t;

/**
   Picture Parameter Set
   @see 7.3.2.3 Picture parameter set RBSP syntax
   @see read_hevc_pic_parameter_set_rbsp
   @see write_hevc_pic_parameter_set_rbsp
   @see debug_pps
*/
typedef struct 
{
    int pic_parameter_set_id;
    int seq_parameter_set_id;
    int dependent_slice_segments_enabled_flag;
    int output_flag_present_flag;
    int num_extra_slice_header_bits;
    int sign_data_hiding_enabled_flag;
    int cabac_init_present_flag;
    int num_ref_idx_l0_default_active_minus1;
    int num_ref_idx_l1_default_active_minus1;
    int init_qp_minus26;
    int constrained_intra_pred_flag;
    int transform_skip_enabled_flag;
    int cu_qp_delta_enabled_flag;
    int diff_cu_qp_delta_depth;
    int pps_cb_qp_offset;
    int pps_cr_qp_offset;
    int pps_slice_chroma_qp_offsets_present_flag;
    int weighted_pred_flag;
    int weighted_bipred_flag;
    int transquant_bypass_enabled_flag;
    int tiles_enabled_flag;
    int entropy_coding_sync_enabled_flag;
    int num_tile_columns_minus1;
    int num_tile_rows_minus1;
    int uniform_spacing_flag;
    int column_width_minus1[MAX_NUM_TILE_COLUMN];
    int row_height_minus1[MAX_NUM_TILE_ROW];
    int loop_filter_across_tiles_enabled_flag;
    int pps_loop_filter_across_slices_enabled_flag;
    int deblocking_filter_control_present_flag;
    int deblocking_filter_override_enabled_flag;
    int pps_deblocking_filter_disabled_flag;
    int pps_beta_offset_div2;
    int pps_tc_offset_div2;
    int pps_scaling_list_data_present_flag;
    hevc_scaling_list_data_t scaling_list_data;
    int lists_modification_present_flag;
    int log2_parallel_merge_level_minus2;
    int slice_segment_header_extension_present_flag; 
    int pps_extension_present_flag;
    int pps_range_extension_flag;
    int pps_multilayer_extension_flag;
    int pps_3d_extension_flag;
    int pps_extension_5bits;
    hevc_pps_range_ext_t pps_range_ext;
    //TODO: support pps_extension_data_flag;
    //TODO: support SVC/MVC extensions
} hevc_pps_t;

/**
  Reference Picture List Modification
  @see 7.3.6.2 Reference picture list modification syntax
*/
typedef struct
{
    int ref_pic_list_modification_flag_l0;
    int list_entry_l0[MAX_NUM_REF_PICS_L0];
    int ref_pic_list_modification_flag_l1;
    int list_entry_l1[MAX_NUM_REF_PICS_L1];
} hevc_ref_pics_lists_mod_t;

/**
  Weidht Prediction Table
  @see 7.3.6.3 Weighted prediction parameters syntax
*/
typedef struct
{
    int luma_log2_weight_denom;
    int delta_chroma_log2_weight_denom;
    int luma_weight_l0_flag[MAX_NUM_REF_PICS_L0];
    int chroma_weight_l0_flag[MAX_NUM_REF_PICS_L0];
    int delta_luma_weight_l0[MAX_NUM_REF_PICS_L0];
    int luma_offset_l0[MAX_NUM_REF_PICS_L0];
    int delta_chroma_weight_l0[MAX_NUM_REF_PICS_L0][2];
    int delta_chroma_offset_l0[MAX_NUM_REF_PICS_L0][2];

    int luma_weight_l1_flag[MAX_NUM_REF_PICS_L1];
    int chroma_weight_l1_flag[MAX_NUM_REF_PICS_L1];
    int delta_luma_weight_l1[MAX_NUM_REF_PICS_L1];
    int luma_offset_l1[MAX_NUM_REF_PICS_L1];
    int delta_chroma_weight_l1[MAX_NUM_REF_PICS_L1][2];
    int delta_chroma_offset_l1[MAX_NUM_REF_PICS_L1][2];
} hevc_pred_weight_table_t;

/**
  Slice Segment Header
  @see 7.3.6 Slice segment header syntax
  @see read_hevc_slice_header_rbsp
  @see write_hevc_slice_header_rbsp
  @see debug_hevc_slice_header_rbsp
*/
typedef struct
{
    int first_slice_segment_in_pic_flag;
    int no_output_of_prior_pics_flag;
    int pic_parameter_set_id;
    int dependent_slice_segment_flag;
    int slice_segment_address;
    int slice_type;
    int pic_output_flag;
    int colour_plane_id;
    int slice_pic_order_cnt_lsb;
    int short_term_ref_pic_set_sps_flag;
    hevc_st_ref_pic_set_t st_ref_pic_set;
    int short_term_ref_pic_set_idx;
    int num_long_term_sps;
    int num_long_term_pics;
    //TODO: is MAX_NUM_LONG_TERM_REF_PICS long enough?
    int lt_idx_sps[MAX_NUM_LONG_TERM_REF_PICS];
    int poc_lsb_lt[MAX_NUM_LONG_TERM_REF_PICS];
    int used_by_curr_pic_lt_flag[MAX_NUM_LONG_TERM_REF_PICS];
    int delta_poc_msb_present_flag[MAX_NUM_LONG_TERM_REF_PICS];
    int delta_poc_msb_cycle_lt[MAX_NUM_LONG_TERM_REF_PICS];
    int slice_temporal_mvp_enabled_flag;
    int slice_sao_luma_flag;
    int slice_sao_chroma_flag;
    int num_ref_idx_active_override_flag;
    int num_ref_idx_l0_active_minus1;
    int num_ref_idx_l1_active_minus1;
    hevc_ref_pics_lists_mod_t rpld;
    int mvd_l1_zero_flag;
    int cabac_init_flag;
    int collocated_from_l0_flag;
    int collocated_ref_idx;
    hevc_pred_weight_table_t pwt;
    int five_minus_max_num_merge_cand;
    int slice_qp_delta;
    int slice_cb_qp_offset;
    int slice_cr_qp_offset;
    int cu_chroma_qp_offset_enabled_flag;
    int deblocking_filter_override_flag;
    int slice_deblocking_filter_disabled_flag;
    int slice_beta_offset_div2;
    int slice_tc_offset_div2;
    int slice_loop_filter_across_slices_enabled_flag;
    int num_entry_point_offsets;
    int offset_len_minus1;
    int entry_point_offset_minus1[MAX_NUM_ENTRY_POINT_OFFSET];
    int slice_segment_header_extension_length;
    //TODO: slice_segment_header_extension_data_byte

} hevc_slice_header_t;

/**
   Network Abstraction Layer (NAL) unit
   @see 7.3.1 NAL unit syntax
   @see read_nal_unit
   @see write_nal_unit
   @see debug_nal
*/
typedef struct
{
    int forbidden_zero_bit;
    int nal_unit_type;
    int nal_layer_id;
    int nal_temporal_id_plus1;
} hevc_nal_t;

typedef struct
{
    int rbsp_size;
    uint8_t* rbsp_buf;
} hevc_slice_data_rbsp_t;

/**
   Access unit delimiter
   @see 7.3.5 Access unit delimiter syntax
   @see read_hevc_access_unit_delimiter
   @see write_hevc_access_unit_delimiter
*/
typedef struct
{
    int primary_pic_type;
} hevc_aud_t;

/**
   HEVC stream
   Contains data structures for all NAL types that can be handled by this library.  
   When reading, data is read into those, and when writing it is written from those.  
   The reason why they are all contained in one place is that some of them depend on others, we need to 
   have all of them available to read or write correctly.
 */
typedef struct
{
    hevc_nal_t* nal;
    hevc_vps_t* vps;
    hevc_sps_t* sps;
    hevc_pps_t* pps;
    hevc_aud_t* aud;
    hevc_slice_header_t* sh;
    
    hevc_slice_data_rbsp_t* slice_data;
    
    hevc_sps_t* sps_table[32];
    hevc_pps_t* pps_table[256];
} hevc_stream_t;

hevc_stream_t* hevc_new();
void hevc_free(hevc_stream_t* h);
    
    int read_debug_hevc_nal_unit(hevc_stream_t* h, uint8_t* buf, int size);

//Table 7-1 NAL unit type codes
#define HEVC_NAL_UNIT_TYPE_TRAIL_N                    0    // Coded slice segment of a non-TSA, non-STSA trailing picture
#define HEVC_NAL_UNIT_TYPE_TRAIL_R                    1    // Coded slice segment of a non-TSA, non-STSA trailing picture
#define HEVC_NAL_UNIT_TYPE_TSA_N                      2    // Coded slice segment of a TSA picture
#define HEVC_NAL_UNIT_TYPE_TSA_R                      3    // Coded slice segment of a TSA picture
#define HEVC_NAL_UNIT_TYPE_STSA_N                     4    // Coded slice segment of an STSA picture
#define HEVC_NAL_UNIT_TYPE_STSA_R                     5    // Coded slice segment of an STSA picture
#define HEVC_NAL_UNIT_TYPE_RADL_N                     6    // Coded slice segment of a RADL picture
#define HEVC_NAL_UNIT_TYPE_RADL_R                     7    // Coded slice segment of a RADL picture
#define HEVC_NAL_UNIT_TYPE_RASL_N                     8    // Coded slice segment of a RASL picture
#define HEVC_NAL_UNIT_TYPE_RASL_R                     9    // Coded slice segment of a RASL picture
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_N10               10    // Reserved non-IRAP SLNR VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_R11               11    // Reserved non-IRAP sub-layer reference VCL NAL unit
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_N12               12    // Reserved non-IRAP SLNR VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_R13               13    // Reserved non-IRAP sub-layer reference VCL NAL unit
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_N14               14    // Reserved non-IRAP SLNR VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL_R15               15    // Reserved non-IRAP sub-layer reference VCL NAL unit
#define HEVC_NAL_UNIT_TYPE_BLA_W_LP                  16    // Coded slice segment of a BLA picture
#define HEVC_NAL_UNIT_TYPE_BLA_W_RADL                17    // Coded slice segment of a BLA picture
#define HEVC_NAL_UNIT_TYPE_BLA_N_LP                  18    // Coded slice segment of a BLA picture
#define HEVC_NAL_UNIT_TYPE_IDR_W_RADL                19    // Coded slice segment of an IDR picture
#define HEVC_NAL_UNIT_TYPE_IDR_N_LP                  20    // Coded slice segment of an IDR picture
#define HEVC_NAL_UNIT_TYPE_CRA_NUT                   21    // Coded slice segment of a CRA picture
#define HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL22            22    // Reserved IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_IRAP_VCL23            23    // Reserved IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL24                 24    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL25                 25    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL26                 26    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL27                 27    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL28                 28    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL29                 29    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL30                 30    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_RSV_VCL31                 31    // Reserved non-IRAP VCL NAL unit types
#define HEVC_NAL_UNIT_TYPE_VPS_NUT                   32    // Video parameter set
#define HEVC_NAL_UNIT_TYPE_SPS_NUT                   33    // Sequence parameter set
#define HEVC_NAL_UNIT_TYPE_PPS_NUT                   34    // Picture parameter set
#define HEVC_NAL_UNIT_TYPE_AUD_NUT                   35    // Access unit delimiter
#define HEVC_NAL_UNIT_TYPE_EOS_NUT                   36    // End of sequence
#define HEVC_NAL_UNIT_TYPE_EOB_NUT                   37    // End of bitstream
#define HEVC_NAL_UNIT_TYPE_FD_NUT                    38    // Filler data
#define HEVC_NAL_UNIT_TYPE_PREFIX_SEI_NUT            39    // Supplemental enhancement information
#define HEVC_NAL_UNIT_TYPE_SUFFIX_SEI_NUT            40    // Supplemental enhancement information

#define MAX_HEVC_VAL_UNIT_TYPE                       40
 

    
    
    //7.4.3 Table 7-6. Name association to slice_type
#define HEVC_SLICE_TYPE_B        0        // P (P slice)
#define HEVC_SLICE_TYPE_P        1        // B (B slice)
#define HEVC_SLICE_TYPE_I        2        // I (I slice)
    
    
#define HEVC_PROFILE_BASELINE  66
#define HEVC_PROFILE_MAIN      77
#define HEVC_PROFILE_EXTENDED  88
#define HEVC_PROFILE_HIGH     100

// file handle for debug output
extern FILE* h264_dbgfile;

static inline long decimal_to_binary(int n) {
    int remainder;
    long binary = 0, i = 1;
        
    while(n != 0) {
        remainder = n%2;
        n = n/2;
        binary= binary + (remainder*i);
        i = i*10;
    }
    return binary;
}
#ifdef __cplusplus
}
#endif

#endif
