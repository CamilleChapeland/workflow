/*"Author: Shan Qu, Delft University of Technology",*/
/*"Acknowledgement: The main part of this code was written during an internship",*/
/*"                 supervision of Dr. Yimin Sun at Aramco oversea company",*/
/*"First created: September 2018; latest update: April 2019",*/
/*"This code is propriatary under the Delphi Research Consortium",*/
/*" ",*/
/*"product: the header of 2D Joint Migration Inversion and FWMod",*/
/*" "*/

#ifndef SUBFUNCTIONS_H_   /* Include guard */
#define SUBFUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <su.h>
#include <segy.h>
#include <mpi.h>
#include <time.h>
#include "mkl_dfti.h"
#include "mkl_cblas.h"

//conditional debugging
//#define DEBUG

//according to the segy definition, every trace has a 240 byte trace header
#define HDRSIZE 240

#define EPS 1.175494351e-38F

//type definition
//To avoid automatic error check in vim Complex type (single precision)
typedef MKL_Complex8 fcomp;
typedef struct struct_globalconsts
{
    //parameters related to output file name
    char *pc_outfolder;
    char *pc_outlabel;                                  // add label to the output files

    //parameters related to input file names
    // for modeling only
    char *pc_vel;
    char *pc_den;

    // for inversion and modeling
    char *pc_initvel;                                     // input .su initial velocity model
    char *pc_data;                                        // input .su surface data in the time domain
    char *pc_src;                                         // input .su source wavefield in the time domain
    char *pc_model_mask;                                  // input .su mask on the model, if if_model_mask == 1
    char *pc_dip;                                         // input .su dip field used in the directional total variation, if directionalTV_lambda != 0

    //input parameters from IO input
    // for modeling only
    unsigned long stc_multiple_order;                     // the order of internal multiples considered, =0: only primaries

    // for inversion and modeling
    unsigned long stc_outpershot;                         // !=0: output per shot information: image gradient + illumination + slowness gradient + misfit per shot; (WARNING: slower than setting it as 0)
    unsigned long stc_outres;                             // !=0: output residual information: residual before image update + residual before velocity update
    /*unsigned long stc_multiple_order_start;*/
    /*unsigned long stc_multiple_order_end;*/
    unsigned long stc_velupdate_iterjump;                 // =0: angle-independent Full wavefield Migration, >0: Joint Migration Inversion, update velocity every 'velocity_update_iter_jump' iterations
    float fc_reflconstr_w;                                // =0: no reflectivity constraint, >0: weight added to the slowness update
    float fc_DTV_lam;                                     // =0: noTV, else: parameter of TV
    float fc_DTV_mu;                                      // =0: noTV, else: parameter of TV
    float fc_DTV_alp;                                     // =0: noTV, else: parameter of TV
    float fc_fmin;                                        // starting frequency
    float fc_fmax_lower;                                  // the lower ending frequency, for multi-freq strategy
    float fc_fmax_upper;                                  // the upper ending frequency, for multi-freq strategy
    float fc_fmax_tmp;                                    // the current upper frequency
    float fc_oper_vmin;                                   // possible minimum of velocity
    float fc_oper_vmax;                                   // possible maximum of velocity
    float fc_oper_dv;                                     // velocity increment in terms of calculating the propagator table
    float fc_dt;                                          // time sampling
    float fc_dx;                                          // horizontal grid spacing of models
    float fc_dz;                                          // vertical grid spacing of models
    unsigned long stc_nx;                                 // horizontal grid size of models
    unsigned long stc_nxtap;                              // horizontal boundary extension of models
    unsigned long stc_nz;                                 // vertical grid size of models
    unsigned long stc_nt;                                 // size of time samples
    unsigned long stc_nsrc;                               // size of shot numbers
    float fc_angle;                                       // parameter for tapering

    unsigned long stc_multifreqrange;                     // total multi-freq range = 4 for size_iter=15,15,10,10
    unsigned long *pstc_i_iter;                           // iteration index for each multi-freq range = 15,30,40,50 for size_iter=15,15,10,10
    /*unsigned long stc_iter_inner_image;*/
    unsigned long stc_data_masktype;                      // mask applied to the data residual, 1->in the frequency domain, 2->in the time domain
    unsigned long stc_if_model_mask;                      // if add model mask to the update, 0=no model mask, other numbers=apply model mask as multiplication to gradients

    //calculate from previous input parameters
    unsigned long stc_nx2;                                // nx2 = nx + 2 * nxtap
    unsigned long stc_nx2pow2;                            // the closest value, which is the power of 2 and larger than nx2
    unsigned long stc_nzplus1;                            // nzplus1 = nz + 1
    unsigned long stc_nf;                                 // nf = nt / 2 + 1
    float fc_df;                                          // df = 1 / (nt * dt)
    unsigned long stc_i_fmin;                             // the index of fmin
    unsigned long stc_i_fmax_lower;                       // the index of fmax_lower
    unsigned long stc_i_fmax_upper;                       // the index of fmax_upper
    unsigned long stc_i_fmax_tmp;                         // the index of fmax_tmp
    unsigned long stc_nsrcnx;                             // nsrcnx = nsrc * nx
    unsigned long stc_nsrcnx2;                            // nsrcnx2 = nsrc * nx2
    unsigned long stc_i_oper_vmin;                        // the index of oper_vmin
    unsigned long stc_i_oper_vmax;                        // the index of oper_vmax
    unsigned long stc_valid_nf;                           // the valid range of frequency points
    unsigned long stc_valid_nf_tmp;                       // the current valid range of frequency points
    unsigned long stc_valid_nv;                           // the valide range of velocity points for operator
    unsigned long stc_itertotal;                          // total iteration number
    unsigned long stc_vel_smoothnz;                       // smoothing length of velocity in the z direction
    unsigned long stc_vel_smoothnx;                       // smoothing length of velocity in the x direction
    unsigned long stc_vel_smoothN;                        // smoothing times of velocity
    unsigned long stc_illum_smoothnz;                     // smoothing length of illummination matrix in the z direction
    unsigned long stc_illum_smoothnx;                     // smoothing length of illummination matrix in the x direction
    unsigned long stc_illum_smoothN;                      // smoothing times of illummination matrix
} globalconsts;

typedef struct struct_globalarrays
{
    float *pf_inivel;
    float *pf_model_mask;
    fcomp *pfcomp_src;
    fcomp *pfcomp_data;
    float *pf_data_mask;
    float *pf_WxGxtaper;

    fcomp *pfcomp_Wx_Nslice;
    fcomp *pfcomp_Gx_Nslice;
    fcomp *pfcomp_Pmin_Nslice;
    fcomp *pfcomp_Pplus_Nslice;
    fcomp *pfcomp_Qplus_Nslice;
    fcomp *pfcomp_Res_Nslice;

    float *pf_image;
    float *pf_image_tmp;
    float *pf_dimage;

    unsigned long *pst_i_vel;
    float *pf_slow;
    float *pf_dslow;
    float *pf_vel;
    float *pf_den;
    float *pf_data_tmp;
    float *pf_dip_field;
} globalarrays;

typedef struct struct_tmparrays
{
    // temporary arrays
    fcomp *pfcomp_Wx_tmp;
    fcomp *pfcomp_Gx_tmp;
    fcomp *pfcomp_Qmin_slice_tmp;
    fcomp *pfcomp_dPmin0_slice_tmp_tmp;
    fcomp *pfcomp_tmp1;
    fcomp *pfcomp_tmp2;
    fcomp *pfcomp_tmp3;
    fcomp *pfcomp_tmp3_tmp;
    fcomp *pfcomp_tmp4;
    fcomp *pfcomp_tmp;
    fcomp *pfcomp_tmp11;
    fcomp *pfcomp_tmp33;
} tmparrays;

// for modeling
void Get_globalparameters_mod(int argc, char *argv[], globalconsts *p_gbparams);
void Print_globalparameters_mod(const globalconsts *p_gbparams);
void Create_memory_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long MPI_ID, const unsigned long Nslice);
void Create_memory_fortmpvars_mod(const globalconsts *p_gbparams, tmparrays *p_tmpvars);
void Read_input_extended_model_data_fromsu_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars);
void Modeling_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE);
void Free_memory_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars
        , const unsigned long MPI_ID);
void Free_memory_fortmpvars_mod(const globalconsts *p_gbparams, tmparrays *p_tmpvars);
void Write_results_data_2su_mod(const globalconsts *p_gbparams, const fcomp *pfcomp_Pmin_Nslice, const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp);


// for JMI
void Get_globalparameters(int argc, char *argv[], globalconsts *p_gbparams);
void Print_globalparameters(const globalconsts *p_gbparams);
void Create_memory(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long MPI_ID, const unsigned long Nslice);
void Create_memory_fortmpvars(const globalconsts *p_gbparams, tmparrays *p_tmpvars);
void Read_input_extended_model_data_fromsu(const globalconsts *p_gbparams, globalarrays *p_gbvars);
void Get_WxGxtables(const globalconsts *p_gbparams, const unsigned long ifreq, fcomp *pfcomp_Wx_slice, fcomp *pfcomp_Gx_slice);
void Get_WxGxtaper(const globalconsts *p_gbparams, float *pf_WxGxtaper);
void Get_current_operatorindex(const globalconsts *p_gbparams, unsigned long *pst_i_vel, const float *pf_vel);
void Prepare_current_freq_range(globalconsts *p_gbparams, const unsigned long iter);
void Apply_freqtaper(const globalconsts *p_gbparams, fcomp *pfcomp_data, fcomp *pfcomp_src);
void Get_current_data_intime(const globalconsts *p_gbparams, const fcomp *pfcomp_data, float *pf_data_tmp);
void Get_Pplus_Pmin_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, fcomp *pfcomp_Pplus_slice, fcomp *pfcomp_Pmin_slice, fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , const fcomp *pfcomp_src_slice, const unsigned long *pst_i_vel, const float *pf_image);
void Get_Res_freqmask(const globalconsts *p_gbparams, const float *pf_data_mask_slice, const fcomp *pfcomp_data_slice
        , const fcomp *pfcomp_Pmin_slice, fcomp *pfcomp_Res_slice);
void Get_Res_timemask(const globalconsts *p_gbparams, const float *pf_data_mask, const float *pf_data_tmp, const fcomp *pfcomp_Pmin0, fcomp *pfcomp_Res, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp);
void Modeling(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp);
void Get_image_gradient_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , float *pf_dimage_one, float *pf_illummatrix_image_one
        , const unsigned long *pst_i_vel, const float *pf_image);
void Get_image_gradient_sub_output_pershot(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , float *pf_dimage_one_allshot, float *pf_illummatrix_image_one_allshot, float *pf_illummatrix_image_one
        , const unsigned long *pst_i_vel, const float *pf_image);
void Get_image_gradient(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp);
void Get_image_gradient_output_pershot(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter);
void Combine_image_gradient(const globalconsts *p_gbparams, float *pf_dimage, float *pf_illummatrix_image, const float *pf_dimage_Nslice, const float *pf_illummatrix_image_Nslice, const unsigned long Nslice_tmp);
void Combine_image_gradient_output_pershot(const globalconsts *p_gbparams, float *pf_dimage_allshot, float *pf_illummatrix_image_allshot, float *pf_illummatrix_image, const float *pf_dimage_Nslice_allshot, const float *pf_illummatrix_image_Nslice_allshot, const float *pf_illummatrix_image_Nslice, const unsigned long Nslice_tmp);
void Smooth_model_2d(float *pf_model, const unsigned long size_z, const unsigned long size_x, const unsigned long smooth_length_z, const unsigned long smooth_length_x);
void Apply_kxkzfilter(const globalconsts *p_gbparams, const float f_fmax, const float f_dx, const float f_dz, float *pf_tmp, const unsigned long st_size_z, const float *pf_vel);
void Scale_gradient(const globalconsts *p_gbparams, float *pf_dmodel, const float *pf_model, const unsigned long st_size_z);
void Apply_model_mask(const globalconsts *p_gbparams, float *pf_model, const unsigned long st_size_z, const float *pf_model_mask);
void Apply_model_mask_extend(const globalconsts *p_gbparams, float *pf_model, const unsigned long st_size_z, const float *pf_model_mask);
void Apply_illummatrix2model(const globalconsts *p_gbparams, float *pf_dmodel, const float *pf_illummatrix, const unsigned long st_size_z);
void Get_image_tmp_withalpha(const globalconsts *p_gbparams, float *pf_image_tmp, const float *pf_image, const float *pf_dimage, const float f_alpha);
void Get_dPmin0_image_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , fcomp *pfcomp_dPmin0_slice_tmp
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , const unsigned long *pst_i_vel, const float *pf_image, const float *pf_dimage);
void Get_dPmin0_image(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE);
void Get_dPmin0(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE);
float Get_J(const globalconsts *p_gbparams, const fcomp *pfcomp_Res_Nslice
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE);
float Get_J_pershot(const globalconsts *p_gbparams, const fcomp *pfcomp_Res_Nslice
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long isrc);
void Get_slowness_gradient_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , float *pf_dslow_one, float *pf_illummatrix_velocity_one
        , const unsigned long *pst_i_vel, const float *pf_image);
void Get_slowness_gradient_sub_output_pershot(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , float *pf_dslow_one_allshot, float *pf_illummatrix_velocity_one
        , const unsigned long *pst_i_vel, const float *pf_image);
void Get_slowness_gradient(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter);
void Get_extraslownessgradient_reflectivityconst(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long iter);
void Apply_directionalTV(const globalconsts *p_gbparams, float fc_DTV_lam, float fc_DTV_mu, float fc_DTV_alp, const float *pf_dip_field, float *pf_slow, float *pf_dtv_dx_est, float *pf_dtv_dz_est, float *pf_dtv_bx_est, float *pf_dtv_bz_est);
void Get_slowness_gradient_output_pershot(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter);
void Get_dPmin0_velocity_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , fcomp *pfcomp_dPmin0_slice_tmp
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , const unsigned long *pst_i_vel, const float *pf_dslow, const float *pf_image);
void Get_dPmin0_velocity(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE);
void Combine_slowness_gradient(const globalconsts *p_gbparams, float *pf_dslow, float *pf_illummatrix_velocity, const float *pf_dslow_Nslice, const float *pf_illummatrix_velocity_Nslice, const unsigned long Nslice_tmp);
void Combine_slowness_gradient_output_pershot(const globalconsts *p_gbparams, float *pf_dslow_allshot, float *pf_illummatrix_velocity, const float *pf_dslow_Nslice_allshot, const float *pf_illummatrix_velocity_Nslice, const unsigned long Nslice_tmp);
void Get_index_velocity_tmp(const globalconsts *p_gbparams, unsigned long *pst_i_vel, float *pf_vel, const float *pf_slow, const float *pf_dslow, const float *pf_model_mask, const float f_alpha);
void Write_results_model_2su(const globalconsts *p_gbparams, const float *pf_vel, const float *pf_image, const float *pf_J, const float *pf_J_image_allshot, const unsigned long iter);
void Write_results_model_gradient_2su(const globalconsts *p_gbparams, const float *pf_dslow, const float *pf_dimage, const unsigned long iter);
void Write_final_results_2su(const globalconsts *p_gbparams, const float *pf_vel, const float *pf_image);
void Write_results_residual_2su(const globalconsts *p_gbparams, fcomp *pfcomp_Res, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter, const unsigned long if_forimage);
void Free_memory(const globalconsts *p_gbparams, globalarrays *p_gbvars
        , const unsigned long MPI_ID);
void Free_memory_fortmpvars(const globalconsts *p_gbparams, tmparrays *p_tmpvars);

fcomp *Allocate_1d_floatcomplex(unsigned long array_size);
float *Allocate_1d_float(unsigned long array_size);
int *Allocate_1d_int(unsigned long array_size);
unsigned long *Allocate_1d_unsignedlong(unsigned long array_size);
unsigned long Convert2nextpow2(unsigned long in);
float mean(unsigned long N, float *a);
void conv_jmi(const float *a, const float *b, float *c, const unsigned long size_a, const unsigned long size_b);
void DTV_DDfilter(float *pf_tvx, float *pf_tvz, const float *pf_slow, const float *pf_dip_field, const float alpha, const unsigned long nx2, const unsigned long nz);
void DTV_DDTfilter(float *pf_dtv, float *pf_tvx, float *pf_tvz, const float *pf_slow, const float *pf_dip_field, const float alpha, const unsigned long nx2, const unsigned long nz);

#endif // SUBFUNCTIONS_H_
