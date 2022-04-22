#include "subfunctions.h"

//self document
char *sdoc[] =
{
    "2D Joint Migration Inversion code",
    "Input files:",
    "initial_velocity=NULL                 input .su initial velocity model",
    "data=NULL                             input .su surface data in the time domain",
    "source=NULL                           input .su source wavefield in the time domain",
    "model_mask=NULL                       input .su mask on the model, if_model_mask >= 1",
    "dip_field=NULL                        input .su dip field used in the directional total variation, if directionalTV_lambda != 0",

    "Output files:",
    ""
    "outfile_label=jmi                     add label to the output files",

    "Parameters:",
    "output_pershot_info=0                 !=0: output per shot information: image gradient + illumination + slowness gradient + misfit per shot; (WARNING: slower than setting it as 0)",
    "output_residual_info=0                !=0: output residual information: residual before image update + residual before velocity update",
    "velocity_update_iter_jump=0           =0: angle-independent Full wavefield Migration, >0: Joint Migration Inversion, update velocity every 'velocity_update_iter_jump' iterations",
    "reflectivity_const_weight=0.0         =0: no reflectivity constraint, >0: weight added to the slowness update",
    "directionalTV_lambda=0                =0: noTV, else: parameter of TV",
    "directionalTV_mu=0                    =0: noTV, else: parameter of TV",
    "directionalTV_alpha=0                 =0: noTV, else: parameter of TV",
    "fmin=5                                starting frequency",
    "fmax_lower=10                         the lower ending frequency, for multi-freq strategy",
    "fmax_upper=40                         the upper ending frequency, for multi-freq strategy",
    "operator_vmin=500                     possible minimum of velocity",
    "operator_vmax=8000                    possible maximum of velocity",
    "operator_dv=1                         velocity increment in terms of calculating the propagator table",
    "dt=0.004                              time sampling",
    "dx=20                                 horizontal grid spacing of models",
    "dz=10                                 vertical grid spacing of models",
    "size_x=0                              horizontal grid size of models",
    "size_xtap=20                          horizontal boundary extension of models",
    "size_z=0                              vertical grid size of models",
    "size_t=0                              size of time samples",
    "size_src=0                            size of shot numbers",
    "angle=80                              parameter for tapering",
    "size_iter=10,10,10,10                 iteration number for each multi-freq range",
    "data_mask_type=1                      mask applied to the data residual, 1->in the frequency domain, 2->in the time domain",
    "if_model_mask=0                       if add model mask to the update, 0=no model mask, >0=apply model mask as multiplication to gradients",
    "                                      1=mask image gradient only; 2=mask velocity gradient only; 3=mask both gradients",
    "size_velocity_smooth_z=1              smoothing length of velocity in the z direction",
    "size_velocity_smooth_x=1              smoothing length of velocity in the x direction",
    "size_velocity_smooth_N=1              smoothing times of velocity",
    "size_illummatrix_smooth_z=1           smoothing length of illummination matrix in the z direction",
    "size_illummatrix_smooth_x=1           smoothing length of illummination matrix in the x direction",
    "size_illummatrix_smooth_N=1           smoothing times of illummination matrix",

    "",
    "Note: The input data must be in .su format",
    "",
    "Author: Shan Qu, Delft University of Technology",
    "Acknowledgement: The main part of this code was written during an internship",
    "                 supervision of Dr. Yimin Sun at Aramco oversea company",
    "First created: September 2018; latest update: June 2019",
    "This code is propriatary under the Delphi Research Consortium",
    " ",
    "product: JMI - 2D Joint Migration Inversion",
    " ",
   NULL
};

int main(int argc, char *argv[])
{
    FILE *fp = NULL;
    char filename[64];

    //MPI related
    MPI_Init(&argc, &argv);
    int MPI_SIZE_int = 0;
    int MPI_ID_int = 0;
    unsigned long MPI_SIZE = 0;
    unsigned long MPI_ID = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_SIZE_int);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_ID_int);
    MPI_SIZE = (unsigned long)MPI_SIZE_int;
    MPI_ID = (unsigned long)MPI_ID_int;

    double start = MPI_Wtime();

    //get global parameters from input
    unsigned long ierror;
    globalconsts gb_params;
    Get_globalparameters(argc, argv, &gb_params);


    unsigned long Nslice = (unsigned long)ceil((float)gb_params.stc_valid_nf / (float)MPI_SIZE);
    unsigned long Nslice_tmp;

    if (MPI_ID == 0)
    {
        Print_globalparameters(&gb_params);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    globalarrays gb_vars;
    Create_memory(&gb_params, &gb_vars, MPI_ID, Nslice);
    tmparrays tmp_vars;
    Create_memory_fortmpvars(&gb_params, &tmp_vars);

    float *pf_J = NULL;
    float *pf_J_norm = NULL;
    float *pf_J_image_allshot = NULL;
    if (MPI_ID == 0)
    {
        pf_J = Allocate_1d_float(gb_params.stc_itertotal);
        pf_J_norm = Allocate_1d_float(gb_params.stc_itertotal);
        if (gb_params.stc_outpershot != 0)
        {
            pf_J_image_allshot = Allocate_1d_float(gb_params.stc_nsrc);
        }
    }

    float *pf_dtv_dx1_est = NULL;
    float *pf_dtv_dz1_est = NULL;
    float *pf_dtv_bx1_est = NULL;
    float *pf_dtv_bz1_est = NULL;
    if (gb_params.fc_DTV_lam >= EPS)
    {
        if (MPI_ID == 0)
        {
            pf_dtv_dx1_est = Allocate_1d_float(gb_params.stc_nz * gb_params.stc_nx2);
            pf_dtv_dz1_est = Allocate_1d_float(gb_params.stc_nz * gb_params.stc_nx2);
            pf_dtv_bx1_est = Allocate_1d_float(gb_params.stc_nz * gb_params.stc_nx2);
            pf_dtv_bz1_est = Allocate_1d_float(gb_params.stc_nz * gb_params.stc_nx2);
        }
    }

    unsigned long index_WxGx = 0;
    unsigned long index_source = 0;
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_m = 0;
    unsigned long index_slice_o = 0;
    unsigned long iter;
    unsigned long izx;
    unsigned long isrcx;
    unsigned long isrc;
    unsigned long ismooth;
    unsigned long iorder;
    unsigned long iorder_inner;
    unsigned long ifreq;
    unsigned long istep;
    unsigned long i_restart;
    unsigned long stc_i_fmax_tmp_old;
    unsigned long i_alphaprotection;
    unsigned long i_alphaprotection_velocity;
    unsigned long st_if_break_alphaprotection;
    unsigned long iz;
    unsigned long ix;
    float f_numerator = 0.0;
    float f_numerator_full = 0.0;
    float f_denominator = 0.0;
    float f_denominator_full = 0.0;
    float f_alpha = 0.0;
    float f_alpha_velocity = 0.0;
    float f_J_image_old = 0.0;
    float f_J_image_old_full = 0.0;
    float f_J_image = 0.0;
    float f_J_image_pershot = 0.0;
    float f_J_image_full = 0.0;
    float f_J_velocity_old = 0.0;
    float f_J_velocity = 0.0;
    float f_J0 = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;

    fcomp *pfcomp_Res = NULL;
    if (gb_params.stc_outres != 0)
    {
        pfcomp_Res = Allocate_1d_floatcomplex(MPI_SIZE * Nslice * gb_params.stc_nsrcnx2);
    }

    //**********************************************************************************************
    //****joint migration inversion starts from now on, parallel with frequency*********************
    //**********************************************************************************************
    if (MPI_ID == 0)
    {
        printf("\nJoint migration inversion starts from now on:\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_ID == 0)
    {
        //get models and data from input
        Read_input_extended_model_data_fromsu(&gb_params, &gb_vars);

        cblas_scopy(gb_params.stc_nz * gb_params.stc_nx2, gb_vars.pf_inivel, 1, gb_vars.pf_vel, 1);
        for (izx=0; izx<gb_params.stc_nz * gb_params.stc_nx2; izx++)
        {
            gb_vars.pf_slow[izx] = 1.0 / gb_vars.pf_vel[izx];
        }

        Get_current_operatorindex(&gb_params, gb_vars.pst_i_vel, gb_vars.pf_vel);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(gb_vars.pst_i_vel, gb_params.stc_nz * gb_params.stc_nx2
            , MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (gb_params.stc_data_masktype == 1)
    {
        MPI_Bcast(gb_vars.pf_data_mask, gb_params.stc_valid_nf * gb_params.stc_nsrcnx2
                , MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(gb_vars.pfcomp_data, 2 * gb_params.stc_valid_nf * gb_params.stc_nsrcnx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(gb_vars.pfcomp_src, 2 * gb_params.stc_valid_nf * gb_params.stc_nsrcnx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //get operators Wx and Gx
    for (ifreq=MPI_ID; ifreq<gb_params.stc_valid_nf; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * gb_params.stc_valid_nv * gb_params.stc_nx2pow2;
        Get_WxGxtables(&gb_params, ifreq, gb_vars.pfcomp_Wx_Nslice + index_slice, gb_vars.pfcomp_Gx_Nslice + index_slice);
    }

    //get boundary tapering
    Get_WxGxtaper(&gb_params, gb_vars.pf_WxGxtaper);
    MPI_Barrier(MPI_COMM_WORLD);

    //main iteration loop
    for (iter=0; iter<gb_params.stc_itertotal; iter++)
    {
        //prepare frequency range for multi-frequency strategy
        Prepare_current_freq_range(&gb_params, iter);
        Nslice_tmp = ceil((float)gb_params.stc_valid_nf_tmp / (float)MPI_SIZE);
        f_J_image = 0.0;

        //set Pmin and image to 0 once in a while to avoid accumulated error
        i_restart=0;
        if (iter == 0)
        {
            memset(gb_vars.pfcomp_Pmin_Nslice, 0
                    , sizeof(float) * 2 * Nslice_tmp * gb_params.stc_nzplus1 * gb_params.stc_nsrcnx2);
            memset(gb_vars.pf_image, 0, sizeof(float) * gb_params.stc_nzplus1 * gb_params.stc_nx2);
            /*if (gb_params.stc_i_fmax_tmp < gb_params.stc_i_fmax_upper)*/
            /*{*/
                /*Apply_freqtaper(&gb_params, gb_vars.pfcomp_data, gb_vars.pfcomp_src);*/
            /*}*/
            // define that new iteration cycle is started
            i_restart=1;
        }
        else if (gb_params.stc_i_fmax_tmp > stc_i_fmax_tmp_old)
        {
            memset(gb_vars.pfcomp_Pmin_Nslice, 0
                    , sizeof(float) * 2 * Nslice_tmp * gb_params.stc_nzplus1 * gb_params.stc_nsrcnx2);
            memset(gb_vars.pf_image, 0, sizeof(float) * gb_params.stc_nzplus1 * gb_params.stc_nx2);
            /*if (gb_params.stc_i_fmax_tmp < gb_params.stc_i_fmax_upper)*/
            /*{*/
                /*Apply_freqtaper(&gb_params, gb_vars.pfcomp_data, gb_vars.pfcomp_src);*/
            /*}*/
            // define that new iteration cycle is started
            i_restart=1;
        }

        //get current frequency range data in time
        if (MPI_ID == 0)
        {
            printf("..Nslice=%zu, Nslice_tmp=%zu\n", Nslice, Nslice_tmp);
            printf("..iter=%zu, current_fmax=%f\n", iter, gb_params.fc_fmax_tmp);
            if (gb_params.stc_data_masktype == 2)
            {
                if (iter == 0)
                {
                    Get_current_data_intime(&gb_params, gb_vars.pfcomp_data, gb_vars.pf_data_tmp);
                }
                else if (gb_params.stc_i_fmax_tmp > stc_i_fmax_tmp_old)
                {
                    Get_current_data_intime(&gb_params, gb_vars.pfcomp_data, gb_vars.pf_data_tmp);
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        memset(gb_vars.pfcomp_Pplus_Nslice, 0
                , sizeof(float) * 2 * Nslice_tmp * gb_params.stc_nzplus1 * gb_params.stc_nsrcnx2);
        memset(gb_vars.pfcomp_Qplus_Nslice, 0
                , sizeof(float) * 2 * Nslice_tmp * gb_params.stc_nzplus1 * gb_params.stc_nsrcnx2);
        memset(gb_vars.pfcomp_Res_Nslice, 0
                , sizeof(float) * 2 * Nslice_tmp * gb_params.stc_nsrcnx2);

        //**************************************************************************************
        //****imaging starts from now on********************************************************
        //**************************************************************************************
        if (MPI_ID == 0)
        {
            printf("....IMAGING:\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //modeling and residual
        if (MPI_ID == 0)
        {
            printf("..........0. modeling\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        Modeling(&gb_params, &gb_vars, &tmp_vars
                , 0, MPI_ID, MPI_SIZE, Nslice_tmp);

        if (gb_params.stc_outres != 0)
        {
            MPI_Gather(gb_vars.pfcomp_Res_Nslice, 2 * Nslice_tmp * gb_params.stc_nsrcnx2, MPI_FLOAT
                    , pfcomp_Res, 2 * Nslice_tmp * gb_params.stc_nsrcnx2, MPI_FLOAT, 0, MPI_COMM_WORLD);
            if (MPI_ID == 0)
            {
                Write_results_residual_2su(&gb_params, pfcomp_Res, MPI_SIZE, Nslice_tmp, iter, 1);
            }
        }

        //calculate image gradient
        if (MPI_ID == 0)
        {
            printf("..........1. calculate image gradient\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (gb_params.stc_outpershot != 0)
        {
            Get_image_gradient_output_pershot(&gb_params, &gb_vars, &tmp_vars
                    , MPI_ID, MPI_SIZE, Nslice_tmp, iter);
        }
        else
        {
            Get_image_gradient(&gb_params, &gb_vars, &tmp_vars
                    , MPI_ID, MPI_SIZE, Nslice_tmp);
        }

        //calculate wavefield perturbation due to image gradient
        if (MPI_ID == 0)
        {
            printf("..........2. calculate wavefield perturbation due to image gradient\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        Get_dPmin0_image(&gb_params, &gb_vars, &tmp_vars
                , &f_numerator_full, &f_denominator_full
                , MPI_ID, MPI_SIZE);

        if (gb_params.stc_outpershot != 0)
        {
            for (isrc=0; isrc<gb_params.stc_nsrc; isrc++)
            {
                f_J_image_pershot = Get_J_pershot(&gb_params, gb_vars.pfcomp_Res_Nslice, MPI_ID, MPI_SIZE, isrc);
                if (MPI_ID == 0)
                {
                    f_J_image += f_J_image_pershot;
                    pf_J_image_allshot[isrc] = f_J_image_pershot;
                }
            }
        }
        else
        {
            f_J_image = Get_J(&gb_params, gb_vars.pfcomp_Res_Nslice, MPI_ID, MPI_SIZE);
        }
        if (MPI_ID == 0)
        {
            pf_J[iter] = f_J_image;
            // Normalize the energy of the misfit from first inner iteration
            if (i_restart > 0)
            {
               f_J0 = f_J_image;
            }
            pf_J_norm[iter] = f_J_image/f_J0;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //get initial alpha and current J
        if (MPI_ID == 0)
        {
            printf("..........3. calculate initial alpha\n");
            if (f_denominator_full < EPS)
            {
                f_alpha = 1.0;
            }
            else
            {
                f_alpha = f_numerator_full / f_denominator_full;
            }
            //printf("f_alpha=%e, f_numerator_full=%e, f_denominator_full=%e\n", f_alpha, f_numerator_full, f_denominator_full);

            Get_image_tmp_withalpha(&gb_params, gb_vars.pf_image_tmp, gb_vars.pf_image, gb_vars.pf_dimage, f_alpha);
            cblas_scopy(gb_params.stc_nzplus1 * gb_params.stc_nx2, gb_vars.pf_image_tmp, 1, gb_vars.pf_image, 1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(gb_vars.pf_image, gb_params.stc_nzplus1 * gb_params.stc_nx2
                , MPI_FLOAT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (MPI_ID == 0) printf("....IMAGING IS DONE\n");

        if (gb_params.stc_velupdate_iterjump != 0)
        {
            if (((iter+1) % gb_params.stc_velupdate_iterjump) == 0)
            {
                //************************************************************************************
                //***velocity inversion starts from now on********************************************
                //************************************************************************************
                if (MPI_ID == 0)
                {
                    printf("....VELOCITY INVERSION:\n");
                    printf("......velocity inversion (only do it once per iteration):\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);

                //modeling and residual
                if (MPI_ID == 0)
                {
                    printf("..........0. modeling\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);
                Modeling(&gb_params, &gb_vars, &tmp_vars
                        , 0, MPI_ID, MPI_SIZE, Nslice_tmp);

                if (gb_params.stc_outres != 0)
                {
                    MPI_Gather(gb_vars.pfcomp_Res_Nslice, 2 * Nslice_tmp * gb_params.stc_nsrcnx2, MPI_FLOAT
                            , pfcomp_Res, 2 * Nslice_tmp * gb_params.stc_nsrcnx2, MPI_FLOAT, 0, MPI_COMM_WORLD);
                    if (MPI_ID == 0)
                    {
                        Write_results_residual_2su(&gb_params, pfcomp_Res, MPI_SIZE, Nslice_tmp, iter, 0);
                    }
                }


                //calculate velocity gradient
                if (MPI_ID == 0)
                {
                    printf("..........1. calculate velocity gradient\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);

                if (gb_params.stc_outpershot != 0)
                {
                    Get_slowness_gradient_output_pershot(&gb_params, &gb_vars, &tmp_vars
                            , MPI_ID, MPI_SIZE, Nslice_tmp, iter);
                }
                else
                {
                    Get_slowness_gradient(&gb_params, &gb_vars, &tmp_vars
                            , MPI_ID, MPI_SIZE, Nslice_tmp, iter);
                }

                //calculate wavefield perturbation due to velocity gradient
                if (MPI_ID == 0)
                {
                    printf("..........2. calculate wavefield perturbation due to velocity gradient\n");
                }
                MPI_Barrier(MPI_COMM_WORLD);
                Get_dPmin0_velocity(&gb_params, &gb_vars, &tmp_vars
                        , &f_numerator_full, &f_denominator_full
                        , MPI_ID, MPI_SIZE);

                //get initial alpha and current J
                if (MPI_ID == 0)
                {
                    printf("..........3. calculate initial alpha\n");
                    if (f_denominator_full < EPS)
                    {
                        f_alpha_velocity = 1.0;
                    }
                    else
                    {
                        f_alpha_velocity = f_numerator_full / f_denominator_full;
                    }
                    //printf("f_alpha_velocity=%e, numerator=%e, denominator=%e\n", f_alpha_velocity, f_numerator_full, f_denominator_full);

                    Get_index_velocity_tmp(&gb_params, gb_vars.pst_i_vel, gb_vars.pf_vel
                            , gb_vars.pf_slow, gb_vars.pf_dslow, gb_vars.pf_model_mask, f_alpha_velocity);
                    cblas_saxpy(gb_params.stc_nz * gb_params.stc_nx2, f_alpha_velocity, gb_vars.pf_dslow, 1, gb_vars.pf_slow, 1);

                    /*apply DTV*/
                    if (gb_params.fc_DTV_lam >= EPS)
                    {
                        if (gb_params.fc_fmax_tmp > 20.0)
                        {
                            Apply_directionalTV(&gb_params, gb_params.fc_DTV_lam, gb_params.fc_DTV_mu, gb_params.fc_DTV_alp, gb_vars.pf_dip_field, gb_vars.pf_slow, pf_dtv_dx1_est, pf_dtv_dz1_est, pf_dtv_bx1_est, pf_dtv_bz1_est);
                        }
                        for (izx=0; izx<gb_params.stc_nz * gb_params.stc_nx2; izx++)
                        {
                            gb_vars.pf_vel[izx] = 1.0 / gb_vars.pf_slow[izx];
                        }
                        Get_current_operatorindex(&gb_params, gb_vars.pst_i_vel, gb_vars.pf_vel);
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);

                MPI_Bcast(gb_vars.pst_i_vel, gb_params.stc_nz * gb_params.stc_nx2
                        , MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                if (MPI_ID == 0)
                {
                    printf("....VELOCITY INVERSION IS DONE\n");
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        } //end of velocity inversion

        MPI_Barrier(MPI_COMM_WORLD);

        if (MPI_ID == 0)
        {
            Write_results_model_2su(&gb_params, gb_vars.pf_vel, gb_vars.pf_image, pf_J_norm, pf_J_image_allshot, iter);
            Write_results_model_gradient_2su(&gb_params, gb_vars.pf_dslow, gb_vars.pf_dimage, iter);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        stc_i_fmax_tmp_old = gb_params.stc_i_fmax_tmp;
    } //end of one iteration

    if (MPI_ID == 0)
    {
        Write_final_results_2su(&gb_params, gb_vars.pf_vel, gb_vars.pf_image);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    Free_memory(&gb_params, &gb_vars, MPI_ID);
    Free_memory_fortmpvars(&gb_params, &tmp_vars);
    if (gb_params.stc_outres != 0)
    {
        free(pfcomp_Res);
    }

    if (gb_params.fc_DTV_lam >= EPS)
    {
        if (MPI_ID == 0)
        {
            free(pf_dtv_dx1_est);
            free(pf_dtv_dz1_est);
            free(pf_dtv_bx1_est);
            free(pf_dtv_bz1_est);
        }
    }

    MPI_Finalize();
}

