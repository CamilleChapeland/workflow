/*"Author: Shan Qu, Delft University of Technology",*/
/*"Acknowledgement: The main part of this code was written during an internship",*/
/*"                 supervision of Dr. Yimin Sun at Aramco oversea company",*/
/*"First created: September 2018; latest update: April 2019",*/
/*"This code is propriatary under the Delphi Research Consortium",*/
/*" ",*/
/*"product: the subfunctions for 2D Joint Migration Inversion and 2D FWMod",*/
/*" "*/

#include "subfunctions.h"


void Get_globalparameters(int argc, char *argv[], globalconsts *p_gbparams)
{
    //this function:
    //get global parameters from input
    //calculate some other global parameters
    initargs(argc, argv);
    requestdoc(1);
    char filename[64];
    FILE *fp;

    if (!getparstring("outfolder", &(p_gbparams->pc_outfolder)))
    {
        p_gbparams->pc_outfolder = "Results/";
    }

    if (!getparstring("outfile_label", &(p_gbparams->pc_outlabel)))
    {
        p_gbparams->pc_outlabel = "jmi";
    }

    if (!getparstring("initial_velocity", &(p_gbparams->pc_initvel)))
    {
        printf("error! no initial_velocity .su file is specified!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!getparstring("data", &(p_gbparams->pc_data)))
    {
        printf("error! no data .su file is specified!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!getparstring("source", &(p_gbparams->pc_src)))
    {
        printf("error! no source .su file is specified!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!getparulong("if_model_mask", &(p_gbparams->stc_if_model_mask)))
    {
        p_gbparams->stc_if_model_mask = 0;
    }

    if (p_gbparams->stc_if_model_mask >= 1)
    {
        if (!getparstring("model_mask", &(p_gbparams->pc_model_mask)))
        {
            printf("error! no model_mask .su file is specified!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    if (!getparulong("output_pershot_info", &(p_gbparams->stc_outpershot)))
    {
        p_gbparams->stc_outpershot = 0;
    }

    if (!getparulong("output_residual_info", &(p_gbparams->stc_outres)))
    {
        p_gbparams->stc_outres = 0;
    }

    if (!getparulong("velocity_update_iter_jump", &(p_gbparams->stc_velupdate_iterjump)))
    {
        p_gbparams->stc_velupdate_iterjump = 1;
    }

    if (!getparfloat("reflectivity_const_weight", &(p_gbparams->fc_reflconstr_w)))
    {
        p_gbparams->fc_reflconstr_w = 0.0;
    }

    if (!getparfloat("directionalTV_lambda", &(p_gbparams->fc_DTV_lam)))
    {
        p_gbparams->fc_DTV_lam = 0.0;
    }

    if (p_gbparams->fc_DTV_lam >= EPS)
    {
        if (!getparfloat("directionalTV_mu", &(p_gbparams->fc_DTV_mu)))
        {
            p_gbparams->fc_DTV_mu = 0.0;
        }
        if (!getparfloat("directionalTV_alpha", &(p_gbparams->fc_DTV_alp)))
        {
            p_gbparams->fc_DTV_alp = 0.0;
        }
        if (!getparstring("dip_field", &(p_gbparams->pc_dip)))
        {
            printf("error! no dip_field .su file is specified!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    if (!getparfloat("fmin", &(p_gbparams->fc_fmin)))
    {
        p_gbparams->fc_fmin = 4.0;
    }

    if (!getparfloat("fmax_lower", &(p_gbparams->fc_fmax_lower)))
    {
        p_gbparams->fc_fmax_lower = 10.0;
    }

    if (!getparfloat("fmax_upper", &(p_gbparams->fc_fmax_upper)))
    {
        p_gbparams->fc_fmax_upper = 40.0;
    }

    if (!getparfloat("operator_vmin", &(p_gbparams->fc_oper_vmin)))
    {
        p_gbparams->fc_oper_vmin = 500.0;
    }

    if (!getparfloat("operator_vmax", &(p_gbparams->fc_oper_vmax)))
    {
        p_gbparams->fc_oper_vmax = 8000.0;
    }

    if (!getparfloat("operator_dv", &(p_gbparams->fc_oper_dv)))
    {
        p_gbparams->fc_oper_dv = 1.0;
    }

    if (!getparfloat("dx", &(p_gbparams->fc_dx)))
    {
        p_gbparams->fc_dx = 25.0;
    }

    if (!getparfloat("dz", &(p_gbparams->fc_dz)))
    {
        p_gbparams->fc_dz = 10.0;
    }

    if (!getparfloat("dt", &(p_gbparams->fc_dt)))
    {
        p_gbparams->fc_dt = 0.004;
    }

    if (!getparulong("size_x", &(p_gbparams->stc_nx)))
    {
        p_gbparams->stc_nx = 0;
    }

    if (!getparulong("size_xtap", &(p_gbparams->stc_nxtap)))
    {
        p_gbparams->stc_nxtap = 10;
    }

    if (!getparulong("size_z", &(p_gbparams->stc_nz)))
    {
        p_gbparams->stc_nz = 0;
    }

    if (!getparulong("size_t", &(p_gbparams->stc_nt)))
    {
        p_gbparams->stc_nt = 0;
    }

    if (!getparulong("size_src", &(p_gbparams->stc_nsrc)))
    {
        p_gbparams->stc_nsrc = 0;
    }

    if (!getparfloat("angle", &(p_gbparams->fc_angle)))
    {
        p_gbparams->fc_angle = 80.0;
    }

    p_gbparams->stc_multifreqrange = countparval("size_iter");
    p_gbparams->pstc_i_iter = Allocate_1d_unsignedlong(p_gbparams->stc_multifreqrange);

    if (!getnparulong(1, "size_iter", &(p_gbparams->pstc_i_iter[0])))
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (p_gbparams->stc_multifreqrange > 1)
    {
        unsigned long i;
        for (i=1; i<p_gbparams->stc_multifreqrange; i++)
        {
            p_gbparams->pstc_i_iter[i] += p_gbparams->pstc_i_iter[i-1];
        }
    }

    if (!getparulong("data_mask_type", &(p_gbparams->stc_data_masktype)))
    {
        p_gbparams->stc_data_masktype = 1;
    }

    if (!getparulong("velocity_smooth_z", &(p_gbparams->stc_vel_smoothnz)))
    {
        p_gbparams->stc_vel_smoothnz = 1;
    }

    if (!getparulong("velocity_smooth_x", &(p_gbparams->stc_vel_smoothnx)))
    {
        p_gbparams->stc_vel_smoothnx = 1;
    }

    if (!getparulong("velocity_smooth_N", &(p_gbparams->stc_vel_smoothN)))
    {
        p_gbparams->stc_vel_smoothN = 1;
    }

    if (!getparulong("illummatrix_smooth_z", &(p_gbparams->stc_illum_smoothnz)))
    {
        p_gbparams->stc_illum_smoothnz = 1;
    }

    if (!getparulong("illummatrix_smooth_x", &(p_gbparams->stc_illum_smoothnx)))
    {
        p_gbparams->stc_illum_smoothnx = 1;
    }

    if (!getparulong("illummatrix_smooth_N", &(p_gbparams->stc_illum_smoothN)))
    {
        p_gbparams->stc_illum_smoothN = 1;
    }


    //remove files that are existing before appending
    snprintf(filename, sizeof(char) * 64, "%s%s_inverted_image.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_inverted_velocity.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }

    snprintf(filename, sizeof(char) * 64, "%s%s_dimage.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_dslowness.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }

    snprintf(filename, sizeof(char) * 64, "%s%s_residual_imaging.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_residual_inversion.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }

    snprintf(filename, sizeof(char) * 64, "%s%s_dimage_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_illum_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_dslowness_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_misfit_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_extradslowness_fromRC.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    if ((fp = fopen(filename, "rb")) != NULL)
    {
        remove(filename);
    }

    //calculate other parameters
    p_gbparams->stc_nx2 = p_gbparams->stc_nx + 2 * p_gbparams->stc_nxtap;
    p_gbparams->stc_nx2pow2 = Convert2nextpow2(p_gbparams->stc_nx2);
    p_gbparams->stc_nzplus1 = p_gbparams->stc_nz + 1;
    p_gbparams->stc_nf = p_gbparams->stc_nt / 2 + 1;
    p_gbparams->fc_df = 1.0 / ((float) p_gbparams->stc_nt * p_gbparams->fc_dt);
    p_gbparams->stc_i_fmin = round((p_gbparams->fc_fmin / p_gbparams->fc_df));
    p_gbparams->stc_i_fmax_lower = round((p_gbparams->fc_fmax_lower / p_gbparams->fc_df));
    p_gbparams->stc_i_fmax_upper = round((p_gbparams->fc_fmax_upper / p_gbparams->fc_df));
    p_gbparams->stc_nsrcnx = p_gbparams->stc_nsrc * p_gbparams->stc_nx;
    p_gbparams->stc_nsrcnx2 = p_gbparams->stc_nsrc * p_gbparams->stc_nx2;
    p_gbparams->stc_i_oper_vmin = round((p_gbparams->fc_oper_vmin
            / p_gbparams->fc_oper_dv) + 1);
    p_gbparams->stc_i_oper_vmax = round((p_gbparams->fc_oper_vmax
            / p_gbparams->fc_oper_dv) + 1);
    p_gbparams->stc_valid_nf = p_gbparams->stc_i_fmax_upper
        - p_gbparams->stc_i_fmin + 1;
    p_gbparams->stc_valid_nv = p_gbparams->stc_i_oper_vmax
        - p_gbparams->stc_i_oper_vmin + 1;
    p_gbparams->stc_itertotal = p_gbparams->pstc_i_iter[p_gbparams->stc_multifreqrange - 1];
}

void Print_globalparameters(const globalconsts *p_gbparams)
{
    char filename[64];
    snprintf(filename, sizeof(char) * 64, "%s%s_inverted_image.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    printf("... pc_outfile_image = %s\n", filename);
    if (p_gbparams->stc_velupdate_iterjump != 0)
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_inverted_velocity.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        printf("... pc_outfile_velocity = %s\n", filename);
    }
    snprintf(filename, sizeof(char) * 64, "%s%s_misfit.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    printf("... pc_outfile_misfit = %s\n", filename);

    printf("... pc_initvel = %s\n", p_gbparams->pc_initvel);
    printf("... pc_data = %s\n", p_gbparams->pc_data);
    printf("... pc_source = %s\n", p_gbparams->pc_src);
    printf("... velocity_update_iter_jump = %zu\n", p_gbparams->stc_velupdate_iterjump);
    printf("... reflectivity_const_weight = %f\n", p_gbparams->fc_reflconstr_w);

    if (p_gbparams->fc_DTV_lam >= EPS)
    {
        printf("... directionalTV_lambda = %f\n", p_gbparams->fc_DTV_lam);
        printf("... directionalTV_mu = %f\n", p_gbparams->fc_DTV_mu);
        printf("... directionalTV_alpha = %f\n", p_gbparams->fc_DTV_alp);
    }

    printf("... output_pershot_info = %zu\n", p_gbparams->stc_outpershot);
    printf("... output_residual_info = %zu\n", p_gbparams->stc_outres);
    printf("... if_model_mask = %zu ", p_gbparams->stc_if_model_mask);
    if (p_gbparams->stc_if_model_mask <= 0) printf("meaning no mask applied \n");
    if (p_gbparams->stc_if_model_mask == 1) printf("meaning mask image only \n");
    if (p_gbparams->stc_if_model_mask == 2) printf("meaning mask velocity only \n");
    if (p_gbparams->stc_if_model_mask >= 3) printf("meaning mask both image and velocity \n");
    if (p_gbparams->stc_if_model_mask >= 1)
    {
        printf("... pc_model_mask = %s\n", p_gbparams->pc_model_mask);
    }

    printf("... fmin = %f\n", p_gbparams->fc_fmin);
    printf("... fmax_lower = %f\n", p_gbparams->fc_fmax_lower);
    printf("... fmax_upper = %f\n", p_gbparams->fc_fmax_upper);
    printf("... operator_vmin = %f\n", p_gbparams->fc_oper_vmin);
    printf("... operator_vmax = %f\n", p_gbparams->fc_oper_vmax);
    printf("... operator_dv = %f\n", p_gbparams->fc_oper_dv);
    printf("... dt = %f\n", p_gbparams->fc_dt);
    printf("... dx = %f\n", p_gbparams->fc_dx);
    printf("... dz = %f\n", p_gbparams->fc_dz);
    printf("... size_x = %zu\n", p_gbparams->stc_nx);
    printf("... size_xtap = %zu\n", p_gbparams->stc_nxtap);
    printf("... size_z = %zu\n", p_gbparams->stc_nz);
    printf("... size_t = %zu\n", p_gbparams->stc_nt);
    printf("... size_src = %zu\n", p_gbparams->stc_nsrc);
    printf("... angle = %f\n", p_gbparams->fc_angle);

    printf("... iter_multfreq_range = %zu\n", p_gbparams->stc_multifreqrange);
    printf("... size_iter = [ ");
    printf("%zu ", p_gbparams->pstc_i_iter[0]);
    unsigned long i;
    for (i=1; i<p_gbparams->stc_multifreqrange; i++)
    {
        printf("%zu ", p_gbparams->pstc_i_iter[i] - p_gbparams->pstc_i_iter[i-1]);
    }
    printf("]\n");

    printf("... data_mask_type = %zu\n", p_gbparams->stc_data_masktype);

    printf("... size_x2 = %zu\n", p_gbparams->stc_nx2);
    printf("... size_x2pow2 = %zu\n", p_gbparams->stc_nx2pow2);
    printf("... size_zplus1 = %zu\n", p_gbparams->stc_nzplus1);
    printf("... size_f = %zu\n", p_gbparams->stc_nf);
    printf("... df = %f\n", p_gbparams->fc_df);
    printf("... index_fmin = %zu\n", p_gbparams->stc_i_fmin);
    printf("... index_fmax_lower = %zu\n", p_gbparams->stc_i_fmax_lower);
    printf("... index_fmax_upper = %zu\n", p_gbparams->stc_i_fmax_upper);
    printf("... size_srcx = %zu\n", p_gbparams->stc_nsrcnx);
    printf("... size_srcx2 = %zu\n", p_gbparams->stc_nsrcnx2);
    printf("... index_operator_vmin = %zu\n", p_gbparams->stc_i_oper_vmin);
    printf("... index_operator_vmax = %zu\n", p_gbparams->stc_i_oper_vmax);
    printf("... validsize_f = %zu\n", p_gbparams->stc_valid_nf);
    printf("... validsize_v = %zu\n", p_gbparams->stc_valid_nv);
    printf("... itertotal = %zu\n", p_gbparams->stc_itertotal);
    printf("... velocity_smooth_z = %zu\n", p_gbparams->stc_vel_smoothnz);
    printf("... velocity_smooth_x = %zu\n", p_gbparams->stc_vel_smoothnx);
    printf("... velocity_smooth_N = %zu\n", p_gbparams->stc_vel_smoothN);
    printf("... illummatrix_smooth_z = %zu\n", p_gbparams->stc_illum_smoothnz);
    printf("... illummatrix_smooth_x = %zu\n", p_gbparams->stc_illum_smoothnx);
    printf("... illummatrix_smooth_N = %zu\n", p_gbparams->stc_illum_smoothN);

    fflush(stdout);
}

void Create_memory(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long MPI_ID, const unsigned long Nslice)
{
    (p_gbvars)->pf_inivel = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
    if (p_gbparams->stc_if_model_mask >= 1)
    {
        (p_gbvars)->pf_model_mask = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    }
    (p_gbvars)->pfcomp_src = Allocate_1d_floatcomplex(p_gbparams->stc_valid_nf * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_data = Allocate_1d_floatcomplex(p_gbparams->stc_valid_nf * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pf_WxGxtaper = Allocate_1d_float(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);

    if (p_gbparams->stc_data_masktype == 1)
    {
        (p_gbvars)->pf_data_mask = Allocate_1d_float(p_gbparams->stc_valid_nf * p_gbparams->stc_nsrcnx2);
    }
    else if (p_gbparams->stc_data_masktype == 2)
    {
        if (MPI_ID == 0)
        {
            (p_gbvars)->pf_data_mask = Allocate_1d_float(p_gbparams->stc_nt * p_gbparams->stc_nsrcnx2);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (p_gbparams->fc_DTV_lam >= EPS)
    {
        if (MPI_ID == 0)
        {
            (p_gbvars)->pf_dip_field = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        }
    }

    (p_gbvars)->pfcomp_Wx_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2);
    (p_gbvars)->pfcomp_Gx_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2);
    (p_gbvars)->pfcomp_Pmin_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_Pplus_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_Qplus_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_Res_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nsrcnx2);

    (p_gbvars)->pf_image = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    (p_gbvars)->pf_image_tmp = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    (p_gbvars)->pst_i_vel = Allocate_1d_unsignedlong(p_gbparams->stc_nz * p_gbparams->stc_nx2);
    (p_gbvars)->pf_dimage = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    (p_gbvars)->pf_dslow = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);

    if (MPI_ID == 0)
    {
        (p_gbvars)->pf_vel = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        (p_gbvars)->pf_slow = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        (p_gbvars)->pf_data_tmp = Allocate_1d_float(p_gbparams->stc_nt * p_gbparams->stc_nsrcnx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Create_memory_fortmpvars(const globalconsts *p_gbparams, tmparrays *p_tmpvars)
{
    (p_tmpvars)->pfcomp_Wx_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_Gx_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_Qmin_slice_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);
    (p_tmpvars)->pfcomp_dPmin0_slice_tmp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);
    (p_tmpvars)->pfcomp_tmp1 = Allocate_1d_floatcomplex(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_tmp2 = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);
    (p_tmpvars)->pfcomp_tmp3 = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_tmp3_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_tmp4 = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_tmpvars)->pfcomp_tmp11 = Allocate_1d_floatcomplex(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_tmpvars)->pfcomp_tmp33 = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);
}

void Create_memory_fortmpvars_mod(const globalconsts *p_gbparams, tmparrays *p_tmpvars)
{
    (p_tmpvars)->pfcomp_Wx_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    (p_tmpvars)->pfcomp_Qmin_slice_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);
}

void Read_input_extended_model_data_fromsu(const globalconsts *p_gbparams, globalarrays *p_gbvars)
{
    //this function:
    //read initial velocity/model_mask(if exist) from .su, and save into size_z * size_x2
    //read source and data from .su, and save into size_f * size_src * size_x2
    segy seg_tmp;
    FILE *fp = NULL;
    unsigned long iz;
    unsigned long ix;
    unsigned long ifreq;
    unsigned long it;
    unsigned long fldr_old;
    unsigned long source_index_tmp;

    //read initial_velocity
    fp = fopen(p_gbparams->pc_initvel, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_initvel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_initvel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                && (p_gbparams->stc_nz == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dz == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_initvel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pf_inivel[], pf_inivel[] size: size_z * size_x2
    cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
            , p_gbvars->pf_inivel + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
    while(fvgettr(fp, &seg_tmp))
    {
        cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                , p_gbvars->pf_inivel + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
    }

    fclose(fp);

    //extend the initial_velocity model with xtap
    unsigned long index_iz = 0;
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
        {
            p_gbvars->pf_inivel[index_iz + ix] = p_gbvars->pf_inivel[index_iz + p_gbparams->stc_nxtap];
        }
        for (ix=p_gbparams->stc_nx+p_gbparams->stc_nxtap
                ; ix<p_gbparams->stc_nx2; ix++)
        {
            p_gbvars->pf_inivel[index_iz + ix]
                = p_gbvars->pf_inivel[index_iz + p_gbparams->stc_nx + p_gbparams->stc_nxtap - 1];
        }
    }

    printf("input initial_velocity from file %s is read in successfully, nz = %zu, nx = %zu, dz = %f, dx = %f\n"
            ,p_gbparams->pc_initvel, p_gbparams->stc_nz, p_gbparams->stc_nx
            , p_gbparams->fc_dz, p_gbparams->fc_dx);
    fflush(stdout);

    //read mask of model, if exist
    if (p_gbparams->stc_if_model_mask >= 1)
    {
        fp = fopen(p_gbparams->pc_model_mask, "rb");
        if (fp == NULL)
        {
            printf("fail to open file %s\n", p_gbparams->pc_model_mask);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //read in first trace
        if (!fvgettr(fp, &seg_tmp))
        {
            printf("can't get the first trace in file %s\n!", p_gbparams->pc_model_mask);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //check if the input parameters match the ones in .su header
        if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                    && (p_gbparams->stc_nz == (unsigned long)seg_tmp.ns)
                    && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                    && (p_gbparams->fc_dz == (float)seg_tmp.d1)))
        {
            printf("the input parameters don't match the ones in .su header in file %s\n!"
                    , p_gbparams->pc_model_mask);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //copy seg_tmp.data[] to pf_model_mask[], pf_model_mask[] size: size_z * size_x2
        cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                , p_gbvars->pf_model_mask + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
        while(fvgettr(fp, &seg_tmp))
        {
            cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                    , p_gbvars->pf_model_mask + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
        }

        fclose(fp);

        //extend the model_mask model with xtap
        unsigned long index_iz = 0;
        for (iz=0; iz<p_gbparams->stc_nz; iz++)
        {
            index_iz = iz * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
            {
                p_gbvars->pf_model_mask[index_iz + ix] = p_gbvars->pf_model_mask[index_iz + p_gbparams->stc_nxtap];
            }
            for (ix=p_gbparams->stc_nx+p_gbparams->stc_nxtap
                    ; ix<p_gbparams->stc_nx2; ix++)
            {
                p_gbvars->pf_model_mask[index_iz + ix]
                    = p_gbvars->pf_model_mask[index_iz + p_gbparams->stc_nx + p_gbparams->stc_nxtap - 1];
            }
        }

        printf("input model_mask from file %s is read in successfully, nz = %zu, nx = %zu, dz = %f, dx = %f\n"
                ,p_gbparams->pc_model_mask, p_gbparams->stc_nz, p_gbparams->stc_nx
                , p_gbparams->fc_dz, p_gbparams->fc_dx);
        fflush(stdout);
    }

    #ifdef DEBUG
    fp = fopen("Tmp/vel0_test.bin", "wb");
    fwrite(p_gbvars->pf_inivel, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nx2
            , 1, fp);
    fclose(fp);
    if (p_gbparams->stc_if_model_mask >= 1)
    {
        fp = fopen("Tmp/modelmask_test.bin", "wb");
        fwrite(p_gbvars->pf_model_mask, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
                , 1, fp);
        fclose(fp);
    }
    #endif

    //read dip field, if exist
    if (p_gbparams->fc_DTV_lam >= EPS)
    {
        fp = fopen(p_gbparams->pc_dip, "rb");
        if (fp == NULL)
        {
            printf("fail to open file %s\n", p_gbparams->pc_dip);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //read in first trace
        if (!fvgettr(fp, &seg_tmp))
        {
            printf("can't get the first trace in file %s\n!", p_gbparams->pc_dip);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //check if the input parameters match the ones in .su header
        if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                    && (p_gbparams->stc_nz == (unsigned long)seg_tmp.ns)
                    && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                    && (p_gbparams->fc_dz == (float)seg_tmp.d1)))
        {
            printf("the input parameters don't match the ones in .su header in file %s\n!"
                    , p_gbparams->pc_dip);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //copy seg_tmp.data[] to pf_dip_field[], pf_dip_field[] size: size_z * size_x2
        cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                , p_gbvars->pf_dip_field + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
        while(fvgettr(fp, &seg_tmp))
        {
            cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                    , p_gbvars->pf_dip_field + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
        }

        fclose(fp);

        //extend the dip_field with xtap
        unsigned long index_iz = 0;
        for (iz=0; iz<p_gbparams->stc_nz; iz++)
        {
            index_iz = iz * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
            {
                p_gbvars->pf_dip_field[index_iz + ix] = p_gbvars->pf_dip_field[index_iz + p_gbparams->stc_nxtap];
            }
            for (ix=p_gbparams->stc_nx+p_gbparams->stc_nxtap
                    ; ix<p_gbparams->stc_nx2; ix++)
            {
                p_gbvars->pf_dip_field[index_iz + ix]
                    = p_gbvars->pf_dip_field[index_iz + p_gbparams->stc_nx + p_gbparams->stc_nxtap - 1];
            }
        }

        printf("input dip_field from file %s is read in successfully, nz = %zu, nx = %zu, dz = %f, dx = %f\n"
                ,p_gbparams->pc_dip, p_gbparams->stc_nz, p_gbparams->stc_nx
                , p_gbparams->fc_dz, p_gbparams->fc_dx);
        fflush(stdout);
    }

    //read source wavefield and convert it into frequency domain
    fp = fopen(p_gbparams->pc_src, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fldr_old = seg_tmp.fldr;

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                && (p_gbparams->stc_nt == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dt == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pfcomp_src[], pfcomp_src[] size: size_validsize_f * size_src * size_x2
    //preparation for fft 1d real2halfcomplex
    //size of vector after fft is half size_f = size_t/2 + 1
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
    status = DftiCommitDescriptor(my_handle);

    fcomp *pfcomp_tmp = NULL;
    pfcomp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nt);

    //fft
    for (it=0; it<p_gbparams->stc_nt; it++)
    {
        pfcomp_tmp[it].real = seg_tmp.data[it];
        pfcomp_tmp[it].imag = 0.0;
    }
    status = DftiComputeForward(my_handle, pfcomp_tmp);

    cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
            , p_gbvars->pfcomp_src + seg_tmp.tracf - 1  + p_gbparams->stc_nxtap, p_gbparams->stc_nsrcnx2);


    source_index_tmp = 0;
    while(fvgettr(fp, &seg_tmp))
    {
        //fft
        for (it=0; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_tmp[it].real = seg_tmp.data[it];
            pfcomp_tmp[it].imag = 0.0;
        }
        status = DftiComputeForward(my_handle, pfcomp_tmp);

        if (seg_tmp.fldr != fldr_old)
        {
            source_index_tmp += 1;
            fldr_old = seg_tmp.fldr;
        }

        cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
                , p_gbvars->pfcomp_src + source_index_tmp * p_gbparams->stc_nx2 + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap
                , p_gbparams->stc_nsrcnx2);
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    fclose(fp);

    printf("input source from file %s is read in successfully, nt = %zu, nfeval = %zu, nx = %zu, dt = %f, dx = %f\n"
            , p_gbparams->pc_src, p_gbparams->stc_nt, p_gbparams->stc_valid_nf
            , p_gbparams->stc_nx, p_gbparams->fc_dt, p_gbparams->fc_dx);
    fflush(stdout);

    float *pf_data = NULL;
    pf_data = Allocate_1d_float(p_gbparams->stc_nt * p_gbparams->stc_nsrcnx2);
    //read data and convert it into frequency domain
    fp = fopen(p_gbparams->pc_data, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_data);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_data);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fldr_old = seg_tmp.fldr;

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                && (p_gbparams->stc_nt == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dt == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_data);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pfcomp_data[], pfcomp_data[] size: size_validsize_f * size_src * size_x2
    //preparation for fft 1d real2halfcomplex
    //size of vector after fft is half size_f = size_t/2 + 1

    //fft
    for (it=0; it<p_gbparams->stc_nt; it++)
    {
        pfcomp_tmp[it].real = seg_tmp.data[it];
        pfcomp_tmp[it].imag = 0.0;
    }
    status = DftiComputeForward(my_handle, pfcomp_tmp);

    cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
            , p_gbvars->pfcomp_data + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap, p_gbparams->stc_nsrcnx2);
    cblas_scopy(p_gbparams->stc_nt, seg_tmp.data, 1
            , pf_data + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap, p_gbparams->stc_nsrcnx2);

    source_index_tmp = 0;
    while(fvgettr(fp, &seg_tmp))
    {
        //fft
        for (it=0; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_tmp[it].real = seg_tmp.data[it];
            pfcomp_tmp[it].imag = 0.0;
        }
        status = DftiComputeForward(my_handle, pfcomp_tmp);

        if (seg_tmp.fldr != fldr_old)
        {
            source_index_tmp += 1;
            fldr_old = seg_tmp.fldr;
        }

        cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
                , p_gbvars->pfcomp_data + source_index_tmp * p_gbparams->stc_nx2 + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap
                , p_gbparams->stc_nsrcnx2);
        cblas_scopy(p_gbparams->stc_nt, seg_tmp.data, 1
                , pf_data + source_index_tmp * p_gbparams->stc_nx2 + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap
                , p_gbparams->stc_nsrcnx2);
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }

    status = DftiFreeDescriptor(&my_handle);
    fclose(fp);

    printf("input data from file %s is read in successfully, nt = %zu, nfeval = %zu, nx = %zu, dt = %f, dx = %f\n"
            , p_gbparams->pc_data, p_gbparams->stc_nt, p_gbparams->stc_valid_nf
            , p_gbparams->stc_nx, p_gbparams->fc_dt, p_gbparams->fc_dx);
    fflush(stdout);

    if (p_gbparams->stc_data_masktype == 1)
    {
        for (it=0; it<p_gbparams->stc_valid_nf*p_gbparams->stc_nsrcnx2; it++)
        {
            if (!((p_gbvars->pfcomp_data[it].real * p_gbvars->pfcomp_data[it].real + p_gbvars->pfcomp_data[it].imag * p_gbvars->pfcomp_data[it].imag) < EPS))
            {
                p_gbvars->pf_data_mask[it] = 1.0;
            }
        }
    }
    else if (p_gbparams->stc_data_masktype == 2)
    {
        for (it=0; it<p_gbparams->stc_nt*p_gbparams->stc_nsrcnx2; it++)
        {
            if (!((pf_data[it] * pf_data[it]) < EPS))
            {
                p_gbvars->pf_data_mask[it] = 1.0;
            }
        }
    }

    #ifdef DEBUG
    fp = fopen("Tmp/source_freq_test.bin", "wb");
    fwrite(p_gbvars->pfcomp_src, sizeof(float) * 2 * p_gbparams->stc_valid_nf
            * p_gbparams->stc_nsrcnx2, 1, fp);
    fclose(fp);
    fp = fopen("Tmp/data_freq_test.bin", "wb");
    fwrite(p_gbvars->pfcomp_data, sizeof(float) * 2 * p_gbparams->stc_valid_nf
            * p_gbparams->stc_nsrcnx2, 1, fp);
    fclose(fp);
    fp = fopen("Tmp/data_time_test.bin", "wb");
    fwrite(pf_data, sizeof(float) * p_gbparams->stc_nt
            * p_gbparams->stc_nsrcnx2, 1, fp);
    fclose(fp);
    if (p_gbparams->stc_data_masktype == 1)
    {
        fp = fopen("Tmp/mask_freq_test.bin", "wb");
        fwrite(p_gbvars->pf_data_mask, sizeof(float) * p_gbparams->stc_valid_nf
                * p_gbparams->stc_nsrcnx2, 1, fp);
        fclose(fp);
    }
    else if (p_gbparams->stc_data_masktype == 2)
    {
        fp = fopen("Tmp/mask_time_test.bin", "wb");
        fwrite(p_gbvars->pf_data_mask, sizeof(float) * p_gbparams->stc_nt
                * p_gbparams->stc_nsrcnx2, 1, fp);
        fclose(fp);
    }
    #endif

    free(pfcomp_tmp);
    pfcomp_tmp = NULL;
    free(pf_data);
    pf_data = NULL;
}

void Get_WxGxtables(const globalconsts *p_gbparams, const unsigned long ifreq, fcomp *pfcomp_Wx_slice, fcomp *pfcomp_Gx_slice)
{
    //this function:
    //calculate the table of Wx and Gx operators along a bunch of velocity
    //pfcomp_Wx_slice and pfcomp_Gx_slice size: stc_valid_nv * stc_valid_nx2pow2
    //prepare parameters
    unsigned long index_iv = 0;
    float f_dkx;
    f_dkx = 2.0 * M_PI / ((float)p_gbparams->stc_nx2pow2 * p_gbparams->fc_dx);
    float f_dom;
    f_dom = 2.0 * M_PI / ((float)p_gbparams->stc_nt * p_gbparams->fc_dt);
    float f_s = 0.0;
    float f_om = 0.0;
    float f_k2 = 0.0;
    float f_tmp = 0.0;
    float f_kx = 0.0;
    float f_eps_coef = 0.01;
    float f_eps = 0.0;
    fcomp fcomp_kz_multiply_dz;
    fcomp_kz_multiply_dz.real = 0.0;
    fcomp_kz_multiply_dz.imag = 0.0;
    float f_kx2_minus_k2 = 0.0;
    unsigned long ix;
    unsigned long iv;
    unsigned long ikx;

    //kx -> for calculating kz, see Page 16 in Xander's thesis
    float *pf_kx = NULL;
    pf_kx = Allocate_1d_float(p_gbparams->stc_nx2pow2);
    for (ikx=0; ikx<p_gbparams->stc_nx2pow2; ikx++)
    {
        if (ikx < p_gbparams->stc_nx2pow2 / 2 + 1)
        {
            pf_kx[ikx] = (float)ikx * f_dkx;
        }
        else
        {
            pf_kx[ikx] = ((float)ikx - (float)p_gbparams->stc_nx2pow2) * f_dkx;
        }
    }

    //calculate taper, plus the scaling 1/p_gbparams->stc_nx2pow2 for ifft
    int i_tmp = floor(p_gbparams->fc_dz * tan(p_gbparams->fc_angle * M_PI / 180.0)
            / p_gbparams->fc_dx);
    i_tmp = floor((float)(p_gbparams->stc_nx2pow2 - 2 * i_tmp) / 2.0);
    float *pf_WxGxtaper = NULL;
    pf_WxGxtaper = Allocate_1d_float(p_gbparams->stc_nx2pow2);
    for (ix=0; ix<i_tmp; ix++)
    {
        pf_WxGxtaper[ix] = (0.5 + 0.5 * cos(- M_PI * (float)(ix + 1) / i_tmp))
            / (float)p_gbparams->stc_nx2pow2;
    }
    for (ix=p_gbparams->stc_nx2pow2/2-1; ix<p_gbparams->stc_nx2pow2; ix++)
    {
        pf_WxGxtaper[ix] = pf_WxGxtaper[p_gbparams->stc_nx2pow2 - ix - 1];
    }

    //preparation for ifft 1d complex2complex
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nx2pow2);
    status = DftiCommitDescriptor(my_handle);
    f_om = f_dom * (float)(p_gbparams->stc_i_fmin + ifreq);
    for (iv=0; iv<p_gbparams->stc_valid_nv; iv++)
    {
        index_iv = iv * p_gbparams->stc_nx2pow2;
        f_s = 1.0 / ((float)iv * p_gbparams->fc_oper_dv + p_gbparams->fc_oper_vmin);
        //f_k2:  k^2 = (w/c)^2 for calculating kz, see Page 16 in Xander's thesis
        f_k2 = f_s * f_s * f_om * f_om;
        f_tmp = f_s * f_om * f_om;
        f_eps = f_eps_coef * f_k2;
        for (ikx=0; ikx<p_gbparams->stc_nx2pow2; ikx++)
        {
            /*
             *if (k^2>kx^2)
             *  kz = sqrt(k^2 - kx^2)
             *else
             *  kz = -j * sqrt(kx^2 - k^2)
             *end
             *in the fk domain,
             *W = exp(-j*kz*dz)
             *G = - j * om * dz * (conj(kz) * k)/(kz^2 + eps) * W
             *see Page 16 and 27 in Xander's thesis
             */
            f_kx2_minus_k2 = pf_kx[ikx] * pf_kx[ikx] - f_k2;
            if (f_kx2_minus_k2 > 0)
            {
                fcomp_kz_multiply_dz.imag = - sqrtf(f_kx2_minus_k2) * p_gbparams->fc_dz;
                pfcomp_Wx_slice[index_iv + ikx].real = exp(fcomp_kz_multiply_dz.imag);
                pfcomp_Gx_slice[index_iv + ikx].real = - f_tmp * fcomp_kz_multiply_dz.imag
                    / (f_kx2_minus_k2 + f_eps) * pfcomp_Wx_slice[index_iv + ikx].real;
            }
            else
            {
                fcomp_kz_multiply_dz.real = sqrtf(- f_kx2_minus_k2) * p_gbparams->fc_dz;
                pfcomp_Wx_slice[index_iv + ikx].real = cos(- fcomp_kz_multiply_dz.real);
                pfcomp_Wx_slice[index_iv + ikx].imag = sin(- fcomp_kz_multiply_dz.real);
                pfcomp_Gx_slice[index_iv + ikx].real = f_tmp * fcomp_kz_multiply_dz.real
                / (- f_kx2_minus_k2 + f_eps) * pfcomp_Wx_slice[index_iv + ikx].imag;
                pfcomp_Gx_slice[index_iv + ikx].imag = - f_tmp * fcomp_kz_multiply_dz.real
                    / (- f_kx2_minus_k2 + f_eps) * pfcomp_Wx_slice[index_iv + ikx].real;
            }
        }
        //transform W and G from fk domain into fx domain
        status = DftiComputeBackward(my_handle, pfcomp_Wx_slice + index_iv);
        status = DftiComputeBackward(my_handle, pfcomp_Gx_slice + index_iv);
        //apply taper to the Wx Gx tables and also divided by size_x2pow2
        for (ix=0; ix<p_gbparams->stc_nx2pow2; ix++)
        {
            pfcomp_Wx_slice[index_iv + ix].real = pf_WxGxtaper[ix] * pfcomp_Wx_slice[index_iv + ix].real;
            pfcomp_Wx_slice[index_iv + ix].imag = pf_WxGxtaper[ix] * pfcomp_Wx_slice[index_iv + ix].imag;
            pfcomp_Gx_slice[index_iv + ix].real = pf_WxGxtaper[ix] * pfcomp_Gx_slice[index_iv + ix].real;
            pfcomp_Gx_slice[index_iv + ix].imag = pf_WxGxtaper[ix] * pfcomp_Gx_slice[index_iv + ix].imag;
        }
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    status = DftiFreeDescriptor(&my_handle);

    free(pf_kx);
    free(pf_WxGxtaper);
    pf_WxGxtaper = NULL;
    pf_kx = NULL;
}

void Get_WxGxtaper(const globalconsts *p_gbparams, float *pf_WxGxtaper)
{
    unsigned long ix = 0;

    float *pf_tmp = NULL;
    pf_tmp = Allocate_1d_float(p_gbparams->stc_nx2);
    for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
    {
        pf_tmp[ix + p_gbparams->stc_nx + p_gbparams->stc_nxtap] = (0.8 + 0.2 * cos(- M_PI * (float)(ix + 1) / (float)p_gbparams->stc_nxtap));
    }
    for (ix=p_gbparams->stc_nxtap; ix<p_gbparams->stc_nxtap + p_gbparams->stc_nx; ix++)
    {
        pf_tmp[ix] = 1.0;
    }
    for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
    {
        pf_tmp[ix] = pf_tmp[p_gbparams->stc_nx2 - ix - 1];
    }

    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans
            , p_gbparams->stc_nx2, p_gbparams->stc_nx2, 1
            , 1.0, pf_tmp, p_gbparams->stc_nx2
            , pf_tmp, p_gbparams->stc_nx2
            , 0.0, pf_WxGxtaper, p_gbparams->stc_nx2);

    free(pf_tmp);
    pf_tmp = NULL;

    #ifdef DEBUG
    FILE *fp = NULL;
    fp = fopen("Tmp/WxGxtaper_test.bin", "wb");
    fwrite(pf_WxGxtaper, sizeof(float) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2, 1, fp);
    fclose(fp);
    #endif
}

void Get_current_operatorindex(const globalconsts *p_gbparams, unsigned long *pst_i_vel, const float *pf_vel)
{
    unsigned long iz;
    unsigned long ix;
    unsigned long index_iz_m = 0;
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz_m = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pst_i_vel[index_iz_m + ix] = (unsigned long)round((pf_vel[index_iz_m + ix]
                        - p_gbparams->fc_oper_vmin) / p_gbparams->fc_oper_dv);
            if (pst_i_vel[index_iz_m + ix] > p_gbparams->stc_valid_nv)
            {
                printf("error! velocity value=%f exceeds the pre-set range!\n", pf_vel[index_iz_m + ix]);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }
}

void Prepare_current_freq_range(globalconsts *p_gbparams, const unsigned long iter)
{
    if ((p_gbparams->stc_i_fmax_lower == p_gbparams->stc_i_fmax_upper) || (p_gbparams->stc_multifreqrange == 1))
    {
        p_gbparams->fc_fmax_tmp = p_gbparams->fc_fmax_upper;
        p_gbparams->stc_i_fmax_tmp = p_gbparams->stc_i_fmax_upper;
    }
    else
    {
        unsigned long istep;
        for (istep=0; istep<p_gbparams->stc_multifreqrange; istep++)
        {
            if (iter < p_gbparams->pstc_i_iter[istep])
            {
                break;
            }
        }

        p_gbparams->fc_fmax_tmp = p_gbparams->fc_fmax_lower
            + ((float)istep / (float)(p_gbparams->stc_multifreqrange - 1))
            * (p_gbparams->fc_fmax_upper - p_gbparams->fc_fmax_lower);
        p_gbparams->stc_i_fmax_tmp = round((p_gbparams->fc_fmax_tmp / p_gbparams->fc_df));
    }
    p_gbparams->stc_valid_nf_tmp = p_gbparams->stc_i_fmax_tmp
        - p_gbparams->stc_i_fmin + 1;
}

void Apply_freqtaper(const globalconsts *p_gbparams, fcomp *pfcomp_data, fcomp *pfcomp_src)
{
    unsigned long ifreq;
    unsigned long nftap;
    unsigned long index_ifreq;
    fcomp fcomp_factor;
    fcomp_factor.real = 0.0;
    fcomp_factor.imag = 0.0;
    nftap = (unsigned long)((float)p_gbparams->stc_valid_nf_tmp * 0.2);

    for (ifreq=p_gbparams->stc_valid_nf_tmp-nftap; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq++)
    {
        index_ifreq = ifreq * p_gbparams->stc_nsrcnx2;
        fcomp_factor.real = powf(cos(((float)(ifreq - p_gbparams->stc_valid_nf_tmp + nftap) * M_PI * 0.5) / (float)(nftap + 1)), 2);
        cblas_cscal(p_gbparams->stc_nsrcnx2, &fcomp_factor, pfcomp_data + index_ifreq, 1);
        cblas_cscal(p_gbparams->stc_nsrcnx2, &fcomp_factor, pfcomp_src + index_ifreq, 1);
        printf("tapering data at ifreq=%zu, factor=%f", ifreq, fcomp_factor.real);
    }
}

void Get_current_data_intime(const globalconsts *p_gbparams, const fcomp *pfcomp_data, float *pf_data_tmp)
{
    //preparation for fft 1d halfcomplex2real
    //size of vector before fft is half size_f = size_t/2 + 1
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
    status = DftiCommitDescriptor(my_handle);

    fcomp *pfcomp_tmp = NULL;
    pfcomp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nt);

    unsigned long isrcx;
    unsigned long it;
    for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
    {
        memset(pfcomp_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nt);
        cblas_ccopy(p_gbparams->stc_valid_nf_tmp, pfcomp_data + isrcx, p_gbparams->stc_nsrcnx2
                , pfcomp_tmp + p_gbparams->stc_i_fmin, 1);

        //fft
        //symetric the trace in the frequency domain, for ifft
        for (it=p_gbparams->stc_nf; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_tmp[it].real = pfcomp_tmp[p_gbparams->stc_nt - it].real;
            pfcomp_tmp[it].imag = - pfcomp_tmp[p_gbparams->stc_nt - it].imag;
        }
        status = DftiComputeBackward(my_handle, pfcomp_tmp);
        for (it=0; it<p_gbparams->stc_nt; it++)
        {
            pf_data_tmp[it * p_gbparams->stc_nsrcnx2 + isrcx] = pfcomp_tmp[it].real / (float)p_gbparams->stc_nt;
        }
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    status = DftiFreeDescriptor(&my_handle);

    #ifdef DEBUG
    FILE *fp;
    fp = fopen("Tmp/current_data_time_test.bin", "wb");
    fwrite(pf_data_tmp, sizeof(float) * p_gbparams->stc_nt
            * p_gbparams->stc_nsrcnx2, 1, fp);
    fclose(fp);
    #endif

    free(pfcomp_tmp);
    pfcomp_tmp = NULL;
}

void Get_Pplus_Pmin_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, fcomp *pfcomp_Pplus_slice, fcomp *pfcomp_Pmin_slice, fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , const fcomp *pfcomp_src_slice, const unsigned long *pst_i_vel, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;

    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(pfcomp_Pplus_slice, 0, sizeof(float) * 2
            * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == 0
                 *   Qplus(0) = - R .* Pmin(0) + S;
                 *else
                 *   Qplus(iz) = (1 + R) .* Pplus(iz) - R .* Pmin(iz);
                 *end
                 */
                if (iz == 0)
                {
                    pfcomp_Qplus_slice[index_isrc + ix].real
                        = - pf_image[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].real
                        + pfcomp_src_slice[isrc * p_gbparams->stc_nx2 + ix].real;
                    pfcomp_Qplus_slice[index_isrc + ix].imag
                        = - pf_image[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].imag
                        + pfcomp_src_slice[isrc * p_gbparams->stc_nx2 + ix].imag;
                }
                else
                {
                    pfcomp_Qplus_slice[index_isrc + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * pfcomp_Pplus_slice[index_isrc + ix].real
                        - pf_image[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].real;
                    pfcomp_Qplus_slice[index_isrc + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * pfcomp_Pplus_slice[index_isrc + ix].imag
                        - pf_image[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].imag;
                }
            }
        }

        /*
         *calculate Pplus(iz+1) = Wx * Qplus(iz)
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Qplus_slice + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, pfcomp_Pplus_slice + index_izplus1, p_gbparams->stc_nx2);
    }

    memset(p_tmpvars->pfcomp_Qmin_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(pfcomp_Pmin_slice, 0, sizeof(float) * 2
            * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == nz
                 *   Qmin(nz) = R .* Pplus(nz);
                 *else
                 *   Qmin(iz) = (1 - R) .* Pmin(iz) + R .* Pplus(iz);
                 *end
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                    p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].real
                        + pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                    p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].imag
                        + pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
                }
            }
        }

        /*
         *calculate Pmin(iz-1) = Wx * Qmin(iz)
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_Qmin_slice_tmp, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, pfcomp_Pmin_slice + index_izminus1, p_gbparams->stc_nx2);
    }
}

void Get_Res_freqmask(const globalconsts *p_gbparams, const float *pf_data_mask_slice, const fcomp *pfcomp_data_slice
        , const fcomp *pfcomp_Pmin_slice, fcomp *pfcomp_Res_slice)
{
    unsigned long isrcx;
    for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
    {
        pfcomp_Res_slice[isrcx].real = (pfcomp_data_slice[isrcx].real - pfcomp_Pmin_slice[isrcx].real) * pf_data_mask_slice[isrcx];
        pfcomp_Res_slice[isrcx].imag = (pfcomp_data_slice[isrcx].imag - pfcomp_Pmin_slice[isrcx].imag) * pf_data_mask_slice[isrcx];
    }
}

void Write_results_residual_2su(const globalconsts *p_gbparams, fcomp *pfcomp_Res, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter, const unsigned long if_forimage)
{
    unsigned long ifreq;
    unsigned long it;
    unsigned long iMPI_ID;
    unsigned long index_slice0;
    unsigned long isrcx;
    unsigned long isrc;
    unsigned long ix;
    unsigned long index_isrcx;
    unsigned long index_ifreq;
    unsigned long index_isrc;
    unsigned long index_ix;

    fcomp *pfcomp_Res_full = NULL;
    pfcomp_Res_full = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2 * p_gbparams->stc_nt);
    for (iMPI_ID=0; iMPI_ID<MPI_SIZE; iMPI_ID++)
    {
        for (ifreq=iMPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
        {
            index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
            cblas_ccopy(p_gbparams->stc_nsrcnx2, pfcomp_Res + iMPI_ID * Nslice_tmp * p_gbparams->stc_nsrcnx2 + index_slice0, 1
                    , pfcomp_Res_full + p_gbparams->stc_i_fmin + ifreq, p_gbparams->stc_nt);
        }
    }

    //preparation for fft/ifft 1d complex2real
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
    status = DftiCommitDescriptor(my_handle);

    //open output file and write the residual data
    FILE *fp;
    char filename[64];
    segy seg_tmp;
    memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);

    if (if_forimage == 1)
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_residual_imaging.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    }
    else
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_residual_inversion.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    }
    //set some basic SEG-Y headers
    seg_tmp.trwf = (int)(p_gbparams->stc_nx);
    seg_tmp.ns = (int)p_gbparams->stc_nt;
    seg_tmp.d2 = p_gbparams->fc_dx;
    seg_tmp.timbas = 3; //sort on keyword fldr
    seg_tmp.d1 = p_gbparams->fc_dt;
    seg_tmp.duse = (int)(iter + 1);
    seg_tmp.trid = 1;
    seg_tmp.dt = (int)(0.1+1000000*seg_tmp.d1);
    fp= fopen(filename, "a+");
    if (fp == NULL)
    {
        syserr("error: cannot open file %s to write P file.\n", filename);
    }

    for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
    {
        index_isrc = isrc * p_gbparams->stc_nx2 * p_gbparams->stc_nt;
        seg_tmp.fldr = (int)(isrc + 1);
        for (ix=p_gbparams->stc_nxtap; ix<p_gbparams->stc_nxtap+p_gbparams->stc_nx; ix++)
        {
            index_ix = index_isrc + ix * p_gbparams->stc_nt;
            for (it=p_gbparams->stc_nf; it<p_gbparams->stc_nt; it++)
            {
                pfcomp_Res_full[index_ix + it].real = pfcomp_Res_full[index_ix + p_gbparams->stc_nt - it].real;
                pfcomp_Res_full[index_ix + it].imag = - pfcomp_Res_full[index_ix + p_gbparams->stc_nt - it].imag;
            }
            status = DftiComputeBackward(my_handle, &pfcomp_Res_full[index_ix]);
            for (it=0; it<p_gbparams->stc_nt; it++)
            {
                seg_tmp.data[it] = pfcomp_Res_full[index_ix + it].real / (float)p_gbparams->stc_nt;
            }

            seg_tmp.tracf = (int)(ix - p_gbparams->stc_nxtap + 1);
            seg_tmp.tracl = (int)(isrc * p_gbparams->stc_nx + ix - p_gbparams->stc_nxtap + 1);
            fwrite(&seg_tmp, sizeof(float)
                    , HDRSIZE / sizeof(float) + p_gbparams->stc_nt, fp);
        }
    }
    fclose(fp);

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    status = DftiFreeDescriptor(&my_handle);

    free(pfcomp_Res_full);
    pfcomp_Res_full = NULL;
}

void Get_Res_timemask(const globalconsts *p_gbparams, const float *pf_data_mask, const float *pf_data_tmp, const fcomp *pfcomp_Pmin0, fcomp *pfcomp_Res, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp)
{
    unsigned long ifreq;
    unsigned long it;
    unsigned long iMPI_ID;
    unsigned long index_slice0;
    unsigned long isrcx;
    unsigned long index_isrcx;
    unsigned long index_ifreq;

    fcomp *pfcomp_Pmin0_full = NULL;
    pfcomp_Pmin0_full = Allocate_1d_floatcomplex(p_gbparams->stc_nt * p_gbparams->stc_nsrcnx2);
    for (iMPI_ID=0; iMPI_ID<MPI_SIZE; iMPI_ID++)
    {
        for (ifreq=iMPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
        {
            index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
            cblas_ccopy(p_gbparams->stc_nsrcnx2, pfcomp_Pmin0 + iMPI_ID * Nslice_tmp * p_gbparams->stc_nsrcnx2 + index_slice0, 1
                    , pfcomp_Pmin0_full + p_gbparams->stc_i_fmin + ifreq, p_gbparams->stc_nt);
        }
    }

    //preparation for fft/ifft 1d complex2real
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
    status = DftiCommitDescriptor(my_handle);

    for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
    {
        index_isrcx = isrcx * p_gbparams->stc_nt;
        for (it=p_gbparams->stc_nf; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_Pmin0_full[index_isrcx + it].real = pfcomp_Pmin0_full[index_isrcx + p_gbparams->stc_nt - it].real;
            pfcomp_Pmin0_full[index_isrcx + it].imag = - pfcomp_Pmin0_full[index_isrcx + p_gbparams->stc_nt - it].imag;
        }
        status = DftiComputeBackward(my_handle, &pfcomp_Pmin0_full[index_isrcx]);
        for (it=0; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_Pmin0_full[index_isrcx + it].real
                = (pf_data_tmp[it * p_gbparams->stc_nsrcnx2 + isrcx] - pfcomp_Pmin0_full[index_isrcx + it].real / (float)p_gbparams->stc_nt)
                * pf_data_mask[it * p_gbparams->stc_nsrcnx2 + isrcx];
            pfcomp_Pmin0_full[index_isrcx + it].imag = 0.0;
        }
        status = DftiComputeForward(my_handle, &pfcomp_Pmin0_full[index_isrcx]);
    }

    for (iMPI_ID=0; iMPI_ID<MPI_SIZE; iMPI_ID++)
    {
        for (ifreq=iMPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
        {
            index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
            cblas_ccopy(p_gbparams->stc_nsrcnx2, pfcomp_Pmin0_full + p_gbparams->stc_i_fmin + ifreq, p_gbparams->stc_nt
                    , pfcomp_Res + iMPI_ID * Nslice_tmp * p_gbparams->stc_nsrcnx2 + index_slice0, 1);
        }
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    status = DftiFreeDescriptor(&my_handle);

    free(pfcomp_Pmin0_full);
    pfcomp_Pmin0_full = NULL;
}

void Modeling(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp)
{
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long iorder_inner;

    fcomp *pfcomp_Pmin0_Nslice = NULL;
    fcomp *pfcomp_Pmin0 = NULL;
    fcomp *pfcomp_Res = NULL;
    if (p_gbparams->stc_data_masktype == 2)
    {
        pfcomp_Pmin0_Nslice = Allocate_1d_floatcomplex(Nslice_tmp * p_gbparams->stc_nsrcnx2);
        if (MPI_ID == 0)
        {
            pfcomp_Pmin0 = Allocate_1d_floatcomplex(MPI_SIZE * Nslice_tmp * p_gbparams->stc_nsrcnx2);
            pfcomp_Res = Allocate_1d_floatcomplex(MPI_SIZE * Nslice_tmp * p_gbparams->stc_nsrcnx2);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;
        index_ifreq = ifreq * p_gbparams->stc_nsrcnx2;

        for (iorder_inner=0; iorder_inner<iorder+1; iorder_inner++)
        {
            Get_Pplus_Pmin_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice
                    , p_gbvars->pfcomp_Qplus_Nslice + index_slice
                    , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                    , p_gbvars->pfcomp_src + index_ifreq, p_gbvars->pst_i_vel, p_gbvars->pf_image);
        }

        if (p_gbparams->stc_data_masktype == 1)
        {
            Get_Res_freqmask(p_gbparams, p_gbvars->pf_data_mask + index_ifreq, p_gbvars->pfcomp_data + index_ifreq
                    , p_gbvars->pfcomp_Pmin_Nslice + index_slice, p_gbvars->pfcomp_Res_Nslice + index_slice0);
        }
        else if (p_gbparams->stc_data_masktype == 2)
        {
            cblas_ccopy(p_gbparams->stc_nsrcnx2, p_gbvars->pfcomp_Pmin_Nslice + index_slice, 1
                    , pfcomp_Pmin0_Nslice + index_slice0, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (p_gbparams->stc_data_masktype == 2)
    {
        MPI_Gather(pfcomp_Pmin0_Nslice, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT
                , pfcomp_Pmin0, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT, 0, MPI_COMM_WORLD);
        if (MPI_ID == 0)
        {
            Get_Res_timemask(p_gbparams, p_gbvars->pf_data_mask, p_gbvars->pf_data_tmp, pfcomp_Pmin0, pfcomp_Res, MPI_SIZE, Nslice_tmp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scatter(pfcomp_Res, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT
                , p_gbvars->pfcomp_Res_Nslice, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (p_gbparams->stc_data_masktype == 2)
    {
        free(pfcomp_Pmin0_Nslice);
        if (MPI_ID == 0)
        {
            free(pfcomp_Pmin0);
            free(pfcomp_Res);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Get_image_gradient_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , float *pf_dimage_one, float *pf_illummatrix_image_one
        , const unsigned long *pst_i_vel, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izplus1_m = 0;
    unsigned long index_izplus1_s = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_izminus1_s = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;


    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp1, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp2, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp3, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp3_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(pf_dimage_one, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    memset(pf_illummatrix_image_one, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_s = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izplus1_m = (iz + 1) * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == 0
                 *   Pmin_unified(0) = eye(nx);
                 *else
                 *   Pmin_unified(iz) = (1 - R) .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp1 on the left is Pmin_unified
                 *pfcomp_tmp1 on the right is Qmin_unified from previous iteration
                 */
                if (iz == 0)
                {
                    if (isrc == ix)
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 1.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                    else
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 0.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                }
                else
                {
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*calculate Qmin_unified(iz+1) = Wx' * Pmin_unified(iz)
         *pfcomp_tmp1 on the left is Qmin_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2);

        /*
         *Qmin_Res(iz+1) = Qmin_unified(iz+1) * Res;
         *pfcomp_tmp2 is Qmin_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *imag_gradient(iz+1) = corr(Qmin_Res(iz+1), Pplus(iz+1, :, :));
         *illum_matrix(iz+1) = ||Qmin_unified(iz+1)||^2 .* ||Pplus(iz+1))||^2;
         *pfcomp_tmp2 is Qmin_Res
         */
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            cblas_cdotc_sub(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_tmp2 + ix, p_gbparams->stc_nx2
                    , pfcomp_Pplus_slice + index_izplus1_s + ix, p_gbparams->stc_nx2
                    , &fcomp_tmp);
            pf_dimage_one[index_izplus1_m + ix] = fcomp_tmp.real;
            pf_illummatrix_image_one[index_izplus1_m + ix]
                = (powf(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Pplus_slice + index_izplus1_s + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp1 + index_izplus1 + ix, p_gbparams->stc_nx2), 2));
        }
    }

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1_s = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == nz
                 *   Pplus_unified(nz) = R .* Qmin_unified(nz);
                 *else
                 *   Pplus_unified(iz) = (1 + R) .* Qplus_unified(iz) + R .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp3 on the left is Pplus_unified currently
                 *pfcomp_tmp3 on the right is Qplus_unified from previous iteration
                 *pfcomp_tmp1 is Qmin_unified from previous code
                 */

                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*calculate Qplus_unified(iz-1) = Wx' * Pplus_unified(iz)
         *pfcomp_tmp3 on the left is Qplus_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp3_tmp, p_gbparams->stc_nx2);
        cblas_ccopy(p_gbparams->stc_nx2 * p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3_tmp, 1, p_tmpvars->pfcomp_tmp3, 1);

        /*
         *Qplus_Res(iz-1) = Qplus_unified(iz-1) * Res;
         *pfcomp_tmp2 is Qplus_Res
         *pfcomp_tmp3 is Qplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *imag_gradient(iz-1) = corr(Qplus_Res(iz-1), Pmin(iz-1, :, :));
         *illum_matrix(iz-1) = ||Qplus_unified(iz-1)||^2 .* ||Pmin(iz-1))||^2;
         *pfcomp_tmp2 is Qplus_Res
         *pfcomp_tmp3 is Qplus_unified
         */
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            cblas_cdotc_sub(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_tmp2 + ix, p_gbparams->stc_nx2
                    , pfcomp_Pmin_slice + index_izminus1_s + ix, p_gbparams->stc_nx2
                    , &fcomp_tmp);
            pf_dimage_one[index_izminus1_m + ix] -= fcomp_tmp.real;
            pf_illummatrix_image_one[index_izminus1_m + ix]
                += (powf(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Pmin_slice + index_izminus1_s + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3 + ix, p_gbparams->stc_nx2), 2));
        }
    }

}

void Get_image_gradient_sub_output_pershot(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , float *pf_dimage_one_allshot, float *pf_illummatrix_image_one_allshot, float *pf_illummatrix_image_one
        , const unsigned long *pst_i_vel, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izplus1_m = 0;
    unsigned long index_izplus1_m_allshot = 0;
    unsigned long index_izplus1_s = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_izminus1_m_allshot = 0;
    unsigned long index_izminus1_s = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long isrcx;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;


    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp1, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp2, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp3, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp3_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(pf_dimage_one_allshot, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    memset(pf_illummatrix_image_one_allshot, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    memset(pf_illummatrix_image_one, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_s = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izplus1_m = (iz + 1) * p_gbparams->stc_nx2;
        index_izplus1_m_allshot = (iz + 1) * p_gbparams->stc_nsrcnx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == 0
                 *   Pmin_unified(0) = eye(nx);
                 *else
                 *   Pmin_unified(iz) = (1 - R) .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp1 on the left is Pmin_unified
                 *pfcomp_tmp1 on the right is Qmin_unified from previous iteration
                 */
                if (iz == 0)
                {
                    if (isrc == ix)
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 1.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                    else
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 0.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                }
                else
                {
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*calculate Qmin_unified(iz+1) = Wx' * Pmin_unified(iz)
         *pfcomp_tmp1 on the left is Qmin_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2);

        /*
         *Qmin_Res(iz+1) = Qmin_unified(iz+1) * Res;
         *pfcomp_tmp2 is Qmin_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *imag_gradient(iz+1) = corr(Qmin_Res(iz+1), Pplus(iz+1, :, :));
         *illum_matrix(iz+1) = ||Qmin_unified(iz+1)||^2 .* ||Pplus(iz+1))||^2;
         *pfcomp_tmp2 is Qmin_Res
         */
        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            pf_dimage_one_allshot[index_izplus1_m_allshot + isrcx] = p_tmpvars->pfcomp_tmp2[isrcx].real * pfcomp_Pplus_slice[index_izplus1_s + isrcx].real
                + p_tmpvars->pfcomp_tmp2[isrcx].imag * pfcomp_Pplus_slice[index_izplus1_s + isrcx].imag;
            pf_illummatrix_image_one_allshot[index_izplus1_m_allshot + isrcx]
                = (pfcomp_Pplus_slice[index_izplus1_s + isrcx].real * pfcomp_Pplus_slice[index_izplus1_s + isrcx].real
                        + pfcomp_Pplus_slice[index_izplus1_s + isrcx].imag * pfcomp_Pplus_slice[index_izplus1_s + isrcx].imag)
                * (p_tmpvars->pfcomp_tmp1[index_izplus1 + isrcx].real * p_tmpvars->pfcomp_tmp1[index_izplus1 + isrcx].real
                        + p_tmpvars->pfcomp_tmp1[index_izplus1 + isrcx].imag * p_tmpvars->pfcomp_tmp1[index_izplus1 + isrcx].imag);
        }
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_illummatrix_image_one[index_izplus1_m + ix]
                = (pow(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Pplus_slice + index_izplus1_s + ix, p_gbparams->stc_nx2), 2)
                        * pow(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp1 + index_izplus1 + ix, p_gbparams->stc_nx2), 2));
        }
    }

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1_s = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;
        index_izminus1_m_allshot = (iz - 1) * p_gbparams->stc_nsrcnx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == nz
                 *   Pplus_unified(nz) = R .* Qmin_unified(nz);
                 *else
                 *   Pplus_unified(iz) = (1 + R) .* Qplus_unified(iz) + R .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp3 on the left is Pplus_unified currently
                 *pfcomp_tmp3 on the right is Qplus_unified from previous iteration
                 *pfcomp_tmp1 is Qmin_unified from previous code
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*calculate Qplus_unified(iz-1) = Wx' * Pplus_unified(iz)
         *pfcomp_tmp3 on the left is Qplus_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp3_tmp, p_gbparams->stc_nx2);
        cblas_ccopy(p_gbparams->stc_nx2 * p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3_tmp, 1, p_tmpvars->pfcomp_tmp3, 1);

        /*
         *Qplus_Res(iz-1) = Qplus_unified(iz-1) * Res;
         *pfcomp_tmp2 is Qplus_Res
         *pfcomp_tmp3 is Qplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *imag_gradient(iz-1) = corr(Qplus_Res(iz-1), Pmin(iz-1, :, :));
         *illum_matrix(iz-1) = ||Qplus_unified(iz-1)||^2 .* ||Pmin(iz-1))||^2;
         *pfcomp_tmp2 is Qplus_Res
         *pfcomp_tmp3 is Qplus_unified
         */
        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            pf_dimage_one_allshot[index_izminus1_m_allshot + isrcx] -= (p_tmpvars->pfcomp_tmp2[isrcx].real * pfcomp_Pmin_slice[index_izminus1_s + isrcx].real
                + p_tmpvars->pfcomp_tmp2[isrcx].imag * pfcomp_Pmin_slice[index_izminus1_s + isrcx].imag);
            pf_illummatrix_image_one_allshot[index_izminus1_m_allshot + isrcx]
                += ((pfcomp_Pmin_slice[index_izminus1_s + isrcx].real * pfcomp_Pmin_slice[index_izminus1_s + isrcx].real
                        + pfcomp_Pmin_slice[index_izminus1_s + isrcx].imag * pfcomp_Pmin_slice[index_izminus1_s + isrcx].imag)
                * (p_tmpvars->pfcomp_tmp3[isrcx].real * p_tmpvars->pfcomp_tmp3[isrcx].real
                        + p_tmpvars->pfcomp_tmp3[isrcx].imag * p_tmpvars->pfcomp_tmp3[isrcx].imag));
        }
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_illummatrix_image_one[index_izminus1_m + ix]
                += (pow(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Pmin_slice + index_izminus1_s + ix, p_gbparams->stc_nx2), 2)
                        * pow(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3 + ix, p_gbparams->stc_nx2), 2));
        }
    }

}

void Get_image_gradient(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp)
{
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_m = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long ismooth;

    float *pf_dimage_Nslice = NULL;
    float *pf_dimage_Nslice_full = NULL;
    float *pf_illummatrix_image = NULL;
    float *pf_illummatrix_image_Nslice = NULL;
    float *pf_illummatrix_image_Nslice_full = NULL;
    pf_dimage_Nslice = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    pf_illummatrix_image_Nslice = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    if (MPI_ID == 0)
    {
        pf_illummatrix_image = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
        pf_dimage_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
        pf_illummatrix_image_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_m = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;

        Get_image_gradient_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice
                , p_gbvars->pfcomp_Res_Nslice + index_slice0
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , pf_dimage_Nslice + index_slice_m, pf_illummatrix_image_Nslice + index_slice_m
                , p_gbvars->pst_i_vel, p_gbvars->pf_image);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(pf_dimage_Nslice, pf_dimage_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pf_illummatrix_image_Nslice, pf_illummatrix_image_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //combine image gradient from different frequency and apply illumination matrix to it and scale it
    if (MPI_ID == 0)
    {
        Combine_image_gradient(p_gbparams, p_gbvars->pf_dimage, pf_illummatrix_image, pf_dimage_Nslice_full, pf_illummatrix_image_Nslice_full, Nslice_tmp);

        if (p_gbparams->fc_fmax_tmp < 12.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 20.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 30.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else
        {
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }

        for (ismooth=0; ismooth<p_gbparams->stc_illum_smoothN; ismooth++)
        {
            Smooth_model_2d(pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbparams->stc_nx2
                    , p_gbparams->stc_illum_smoothnz, p_gbparams->stc_illum_smoothnx);
        }
        Apply_illummatrix2model(p_gbparams, p_gbvars->pf_dimage, pf_illummatrix_image, p_gbparams->stc_nzplus1);

	// apply the model_mask to the image gradient
        if (p_gbparams->stc_if_model_mask == 1 || p_gbparams->stc_if_model_mask >= 3)
        {
            //It also masked the illumination matrix, but this part is maybe not useful
            //Apply_model_mask(p_gbparams, pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_model_mask);
            Apply_model_mask(p_gbparams, p_gbvars->pf_dimage, p_gbparams->stc_nzplus1, p_gbvars->pf_model_mask);
        }

        Scale_gradient(p_gbparams, p_gbvars->pf_dimage, p_gbvars->pf_image, p_gbparams->stc_nzplus1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_gbvars->pf_dimage, p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pf_dimage_Nslice);
    free(pf_illummatrix_image_Nslice);
    if (MPI_ID == 0)
    {
        free(pf_illummatrix_image);
        free(pf_dimage_Nslice_full);
        free(pf_illummatrix_image_Nslice_full);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Get_image_gradient_output_pershot(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter)
{
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_m = 0;
    unsigned long index_slice_m_allshot = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long ismooth;
    unsigned long ix;
    unsigned long isrc;
    unsigned long iz;

    float *pf_dimage_pershot = NULL;
    float *pf_illummatrix_image = NULL;
    float *pf_illummatrix_image_pershot = NULL;
    float *pf_illummatrix_image_allshot = NULL;
    float *pf_dimage_allshot = NULL;
    float *pf_dimage_Nslice_allshot = NULL;
    float *pf_dimage_Nslice_allshot_full = NULL;
    float *pf_illummatrix_image_Nslice_allshot = NULL;
    float *pf_illummatrix_image_Nslice = NULL;
    float *pf_illummatrix_image_Nslice_allshot_full = NULL;
    float *pf_illummatrix_image_Nslice_full = NULL;
    pf_dimage_Nslice_allshot = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    pf_illummatrix_image_Nslice_allshot = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    pf_illummatrix_image_Nslice = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    if (MPI_ID == 0)
    {
        pf_dimage_pershot = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
        pf_dimage_allshot = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_image_allshot = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_image = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
        pf_illummatrix_image_pershot = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
        pf_dimage_Nslice_allshot_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_image_Nslice_allshot_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_image_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_m = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2;
        index_slice_m_allshot = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;

        Get_image_gradient_sub_output_pershot(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice
                , p_gbvars->pfcomp_Res_Nslice + index_slice0
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , pf_dimage_Nslice_allshot + index_slice_m_allshot, pf_illummatrix_image_Nslice_allshot + index_slice_m_allshot, pf_illummatrix_image_Nslice + index_slice_m
                , p_gbvars->pst_i_vel, p_gbvars->pf_image);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(pf_dimage_Nslice_allshot, pf_dimage_Nslice_allshot_full
            , Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pf_illummatrix_image_Nslice_allshot, pf_illummatrix_image_Nslice_allshot_full
            , Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pf_illummatrix_image_Nslice, pf_illummatrix_image_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    FILE *fp;
    char filename[64];

    if (MPI_ID == 0)
    {
        Combine_image_gradient_output_pershot(p_gbparams, pf_dimage_allshot, pf_illummatrix_image_allshot, pf_illummatrix_image, pf_dimage_Nslice_allshot_full, pf_illummatrix_image_Nslice_allshot_full, pf_illummatrix_image_Nslice_full, Nslice_tmp);

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            for (iz=0; iz<p_gbparams->stc_nzplus1; iz++)
            {
                cblas_scopy(p_gbparams->stc_nx2, pf_dimage_allshot + iz * p_gbparams->stc_nsrcnx2 + isrc * p_gbparams->stc_nx2, 1
                        , pf_dimage_pershot + iz * p_gbparams->stc_nx2, 1);
                cblas_scopy(p_gbparams->stc_nx2, pf_illummatrix_image_allshot + iz * p_gbparams->stc_nsrcnx2 + isrc * p_gbparams->stc_nx2, 1
                        , pf_illummatrix_image_pershot + iz * p_gbparams->stc_nx2, 1);
            }
            cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1.0
                    , pf_dimage_pershot, 1, p_gbvars->pf_dimage, 1);

            if (p_gbparams->fc_fmax_tmp < 12.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_illummatrix_image_pershot, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
            }
            else if (p_gbparams->fc_fmax_tmp < 20.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_illummatrix_image_pershot, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
            }
            else if (p_gbparams->fc_fmax_tmp < 30.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_illummatrix_image_pershot, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
            }
            else
            {
                Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_illummatrix_image_pershot, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
            }

            segy seg_tmp;
            memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);
            snprintf(filename, sizeof(char) * 64, "%s%s_dimage_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
            seg_tmp.trwf = (int)p_gbparams->stc_nx;
            seg_tmp.ns = (int)p_gbparams->stc_nzplus1;
            seg_tmp.d2 = p_gbparams->fc_dx;
            seg_tmp.timbas = 3; //sort on keyword fldr
            seg_tmp.d1 = p_gbparams->fc_dz;
            seg_tmp.fldr = (int)(isrc + 1);
            seg_tmp.duse = (int)(iter + 1);
            seg_tmp.trid = 30;
            fp= fopen(filename, "a+");
            if (fp == NULL)
            {
                syserr("error: cannot open file %s to write P file.\n", filename);
            }

            for (ix=0; ix<p_gbparams->stc_nx; ix++)
            {
                cblas_scopy(p_gbparams->stc_nzplus1, pf_dimage_pershot + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                        , seg_tmp.data, 1);
                seg_tmp.tracl = (int)(isrc * p_gbparams->stc_nx + ix + 1);
                seg_tmp.tracf = (int)(ix + 1);
                fwrite(&seg_tmp, sizeof(float)
                        , HDRSIZE / sizeof(float) + p_gbparams->stc_nzplus1, fp);
            }
            fclose(fp);

            snprintf(filename, sizeof(char) * 64, "%s%s_illum_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
            fp= fopen(filename, "a+");
            if (fp == NULL)
            {
                syserr("error: cannot open file %s to write P file.\n", filename);
            }

            for (ix=0; ix<p_gbparams->stc_nx; ix++)
            {
                cblas_scopy(p_gbparams->stc_nzplus1, pf_illummatrix_image_pershot + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                        , seg_tmp.data, 1);
                seg_tmp.tracl = (int)(isrc * p_gbparams->stc_nx + ix + 1);
                seg_tmp.tracf = (int)(ix + 1);
                fwrite(&seg_tmp, sizeof(float)
                        , HDRSIZE / sizeof(float) + p_gbparams->stc_nzplus1, fp);
            }
            fclose(fp);
        }

        if (p_gbparams->fc_fmax_tmp < 12.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 20.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 30.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }
        else
        {
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_vel);
        }

        for (ismooth=0; ismooth<p_gbparams->stc_illum_smoothN; ismooth++)
        {
            Smooth_model_2d(pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbparams->stc_nx2
                    , p_gbparams->stc_illum_smoothnz, p_gbparams->stc_illum_smoothnx);
        }
        Apply_illummatrix2model(p_gbparams, p_gbvars->pf_dimage, pf_illummatrix_image, p_gbparams->stc_nzplus1);

	// apply the model_mask to the image gradient
        if (p_gbparams->stc_if_model_mask == 1 || p_gbparams->stc_if_model_mask >= 3)
        {
            //It also masked the illumination matrix, but this part is maybe not useful
            //Apply_model_mask(p_gbparams, pf_illummatrix_image, p_gbparams->stc_nzplus1, p_gbvars->pf_model_mask);
            Apply_model_mask(p_gbparams, p_gbvars->pf_dimage, p_gbparams->stc_nzplus1, p_gbvars->pf_model_mask);
        }

        Scale_gradient(p_gbparams, p_gbvars->pf_dimage, p_gbvars->pf_image, p_gbparams->stc_nzplus1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_gbvars->pf_dimage, p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pf_dimage_Nslice_allshot);
    free(pf_illummatrix_image_Nslice_allshot);
    free(pf_illummatrix_image_Nslice);
    if (MPI_ID == 0)
    {
        free(pf_dimage_pershot);
        free(pf_dimage_allshot);
        free(pf_illummatrix_image);
        free(pf_illummatrix_image_allshot);
        free(pf_dimage_Nslice_allshot_full);
        free(pf_illummatrix_image_Nslice_allshot_full);
        free(pf_illummatrix_image_Nslice_full);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Combine_image_gradient(const globalconsts *p_gbparams, float *pf_dimage, float *pf_illummatrix_image, const float *pf_dimage_Nslice, const float *pf_illummatrix_image_Nslice, const unsigned long Nslice_tmp)
{
    unsigned long islice;
    unsigned long iz;
    memset(pf_dimage, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    memset(pf_illummatrix_image, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    for (islice=0; islice<Nslice_tmp; islice++)
    {
        cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1.0
                , pf_dimage_Nslice + islice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1, pf_dimage, 1);
        cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1.0
                , pf_illummatrix_image_Nslice + islice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1, pf_illummatrix_image, 1);
    }

    /*taper deepest part of image to avoid edge effects*/
    cblas_sscal(p_gbparams->stc_nx2, 0.25, pf_dimage + (p_gbparams->stc_nzplus1 - 1) * p_gbparams->stc_nx2, 1);
    cblas_sscal(p_gbparams->stc_nx2, 0.5, pf_dimage + (p_gbparams->stc_nzplus1 - 2) * p_gbparams->stc_nx2, 1);
    cblas_sscal(p_gbparams->stc_nx2, 0.75, pf_dimage + (p_gbparams->stc_nzplus1 - 3) * p_gbparams->stc_nx2, 1);
}

void Combine_image_gradient_output_pershot(const globalconsts *p_gbparams, float *pf_dimage_allshot, float *pf_illummatrix_image_allshot, float *pf_illummatrix_image, const float *pf_dimage_Nslice_allshot, const float *pf_illummatrix_image_Nslice_allshot, const float *pf_illummatrix_image_Nslice, const unsigned long Nslice_tmp)
{
    unsigned long islice;
    memset(pf_dimage_allshot, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    memset(pf_illummatrix_image_allshot, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    memset(pf_illummatrix_image, 0, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    for (islice=0; islice<Nslice_tmp; islice++)
    {
        cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2, 1.0
                , pf_dimage_Nslice_allshot + islice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2, 1, pf_dimage_allshot, 1);
        cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2, 1.0
                , pf_illummatrix_image_Nslice_allshot + islice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2, 1, pf_illummatrix_image_allshot, 1);
        cblas_saxpy(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1.0
                , pf_illummatrix_image_Nslice + islice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2, 1, pf_illummatrix_image, 1);
    }
}

void Smooth_model_2d(float *pf_model, const unsigned long size_z, const unsigned long size_x, const unsigned long smooth_length_z, const unsigned long smooth_length_x)
{
    //smooth_length_z and smooth_length_x must be odd integral
    unsigned long zhalflen = (smooth_length_z - 1) / 2;
    unsigned long zstart = zhalflen;
    unsigned long zstop = size_z - zhalflen - 1;

    unsigned long xhalflen = (smooth_length_x - 1) / 2;
    unsigned long xstart = xhalflen;
    unsigned long xstop = size_x - xhalflen - 1;
    unsigned long iz;
    unsigned long inner_iz;
    unsigned long ix;

    float *pf_model_tmp = NULL;
    pf_model_tmp = Allocate_1d_float(size_z * size_x);
    cblas_scopy(size_z * size_x, pf_model, 1, pf_model_tmp, 1);
    float *pf_tmp = NULL;
    pf_tmp = Allocate_1d_float(smooth_length_z * smooth_length_x);

    for (iz=0; iz<size_z; iz++)
    {
        if (iz < zstart)
        {
            zhalflen = iz;
        }
        else if (iz > zstop)
        {
            zhalflen = size_z - iz - 1;
        }
        else
        {
            zhalflen = (smooth_length_z - 1) / 2;
        }
        for (ix=0; ix<size_x; ix++)
        {
            if (ix < xstart)
            {
                xhalflen = ix;
            }
            else if (ix > xstop)
            {
                xhalflen = size_x - ix - 1;
            }
            else
            {
                xhalflen = (smooth_length_x - 1) / 2;
            }

            for (inner_iz=0; inner_iz<(zhalflen*2+1); inner_iz++)
            {
                cblas_scopy(xhalflen * 2 + 1, pf_model_tmp + (iz - zhalflen + inner_iz) * size_x + ix - xhalflen, 1
                        , pf_tmp + inner_iz * (xhalflen * 2 + 1), 1);
            }
            pf_model[iz * size_x + ix] = mean((zhalflen * 2 + 1) * (xhalflen * 2 + 1), pf_tmp);
        }
    }

    free(pf_model_tmp);
    free(pf_tmp);
    pf_model_tmp = NULL;
    pf_tmp = NULL;
}

void Apply_kxkzfilter(const globalconsts *p_gbparams, const float f_fmax, const float f_dx, const float f_dz, float *pf_tmp, const unsigned long st_size_z, const float *pf_vel)
{
    unsigned long st_size_z_2 = st_size_z * 2;
    unsigned long st_size_x_2 = p_gbparams->stc_nx2 * 2;
    float f_dkz = 2.0 * M_PI / ((float)st_size_z_2 * f_dz);
    float f_dkx = 2.0 * M_PI / ((float)st_size_x_2 * f_dx);
    float f_min_velocity
        = pf_vel[cblas_isamin(p_gbparams->stc_nz * p_gbparams->stc_nx2, pf_vel, 1)];
    float f_kmax = 2.0 * M_PI * f_fmax / f_min_velocity;

    unsigned long iz;
    unsigned long ix;
    unsigned long ikz;
    unsigned long ikx;
    unsigned long index_ikz;

    fcomp *pfcomp_tmp2 = NULL;
    float *pf_filter = NULL;
    pfcomp_tmp2 = Allocate_1d_floatcomplex(st_size_z_2 * st_size_x_2);
    pf_filter = Allocate_1d_float(st_size_z_2 * st_size_x_2);

    for (iz=0; iz<st_size_z; iz++)
    {
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pfcomp_tmp2[iz * st_size_x_2 + ix].real = pf_tmp[iz * p_gbparams->stc_nx2 + ix];
        }
    }

    for (ikz=0; ikz<(st_size_z + 1); ikz++)
    {
        for (ikx=0; ikx<(p_gbparams->stc_nx2 + 1); ikx++)
        {
            pf_filter[ikz * st_size_x_2 + ikx] = exp(-(powf((float)ikz * f_dkz, 2) + powf((float)ikx * f_dkx, 2)) / (powf(f_kmax, 2)));
        }
    }
    for (ikx=0; ikx<p_gbparams->stc_nx2; ikx++)
    {
        for (ikz=st_size_z+1; ikz<st_size_z_2; ikz++)
        {
            pf_filter[ikz * st_size_x_2 + ikx] = pf_filter[(st_size_z_2 - ikz) * st_size_x_2 + ikx];
        }
    }
    for (ikz=0; ikz<st_size_z_2; ikz++)
    {
        index_ikz = ikz * st_size_x_2;
        for (ikx=p_gbparams->stc_nx2+1; ikx<st_size_x_2; ikx++)
        {
            pf_filter[index_ikz + ikx] = pf_filter[index_ikz + st_size_x_2 - ikx];
        }
    }

    //preparation for fft 2d complex2complex
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    MKL_LONG l[2];
    l[0] = st_size_z_2;
    l[1] = st_size_x_2;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 2, l);
    status = DftiCommitDescriptor(my_handle);
    status = DftiComputeForward(my_handle, pfcomp_tmp2);

    for (ikz=0; ikz<st_size_z_2; ikz++)
    {
        index_ikz = ikz * st_size_x_2;
        for (ikx=0; ikx<st_size_x_2; ikx++)
        {
            pfcomp_tmp2[index_ikz + ikx].real *= pf_filter[index_ikz + ikx];
            pfcomp_tmp2[index_ikz + ikx].imag *= pf_filter[index_ikz + ikx];
        }
    }

    status = DftiComputeBackward(my_handle, pfcomp_tmp2);
    float Scale = 1.0 / (float)(st_size_z_2 * st_size_x_2);
    for (iz=0; iz<st_size_z; iz++)
    {
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_tmp[iz * p_gbparams->stc_nx2 + ix] = pfcomp_tmp2[iz * st_size_x_2 + ix].real * Scale;
        }
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    status = DftiFreeDescriptor(&my_handle);

    free(pfcomp_tmp2);
    free(pf_filter);
    pfcomp_tmp2 = NULL;
    pf_filter = NULL;
}

void Scale_gradient(const globalconsts *p_gbparams, float *pf_dmodel, const float *pf_model, const unsigned long st_size_z)
{
    float eps;
    eps = 0.01 * (fabs(pf_model[cblas_isamax(st_size_z * p_gbparams->stc_nx2, pf_model, 1)]))
        / (fabs(pf_dmodel[cblas_isamax(st_size_z * p_gbparams->stc_nx2, pf_dmodel, 1)]) + EPS);
    if (eps < EPS)
    {
        eps = 0.001;
    }
    cblas_sscal(st_size_z * p_gbparams->stc_nx2, eps, pf_dmodel, 1);
}

void Apply_model_mask(const globalconsts *p_gbparams, float *pf_model, const unsigned long st_size_z, const float *pf_model_mask)
{
    unsigned long index_model_left;
    unsigned long index_model_right;
    unsigned long iz;
    unsigned long ix;
    unsigned long index_iz;
/*    for (iz=0; iz<p_gbparams->stc_nz; iz++)*/
    //{
        //index_iz = iz * p_gbparams->stc_nx2;
        //for (ix=0; ix<p_gbparams->stc_nx2-1; ix++)
        //{
            //if ((pf_model_mask[index_iz + ix] == 0) && (pf_model_mask[index_iz + ix + 1] == 1))
            //{
                //index_model_left = ix + 1;
            //}
            //else if ((pf_model_mask[index_iz + ix] == 1) && (pf_model_mask[index_iz + ix + 1] == 0))
            //{
                //index_model_right = ix;
            //}
        //}
        //for (ix=0; ix<index_model_left; ix++)
        //{
            //pf_model[index_iz + ix] = pf_model[index_iz + index_model_left];
        //}
        //for (ix=index_model_right+1; ix<p_gbparams->stc_nx2; ix++)
        //{
            //pf_model[index_iz + ix] = pf_model[index_iz + index_model_right];
        //}
    //}
    //for (iz=p_gbparams->stc_nz; iz<st_size_z; iz++)
    //{
        //index_iz = iz * p_gbparams->stc_nx2;
        //for (ix=0; ix<index_model_left; ix++)
        //{
            //pf_model[index_iz + ix] = pf_model[index_iz + index_model_left];
        //}
        //for (ix=index_model_right+1; ix<p_gbparams->stc_nx2; ix++)
        //{
            //pf_model[index_iz + ix] = pf_model[index_iz + index_model_right];
        //}
    /*}*/
    // Apply the mask for the regular nz values
    for (iz=0; iz<p_gbparams->stc_nz * p_gbparams->stc_nx2; iz++)
    {
        pf_model[iz] *= pf_model_mask[iz];
    }
    // If model exceeds the regular nz, copy the last mask value
    for (iz=p_gbparams->stc_nz * p_gbparams->stc_nx2; iz<st_size_z * p_gbparams->stc_nx2; iz++)
    {
        pf_model[iz] *= pf_model_mask[iz - p_gbparams->stc_nx2];
    }
}

void Apply_model_mask_extend(const globalconsts *p_gbparams, float *pf_model, const unsigned long st_size_z, const float *pf_model_mask)
{
    unsigned long index_model_left;
    unsigned long index_model_right;
    unsigned long iz;
    unsigned long ix;
    unsigned long index_iz;
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nx2-1; ix++)
        {
            if ((pf_model_mask[index_iz + ix] < pf_model_mask[index_iz + ix + 1]) && (pf_model_mask[index_iz + ix] == 0))
            {
                index_model_left = ix + 1;
            }
            else if ((pf_model_mask[index_iz + ix] > pf_model_mask[index_iz + ix + 1]) && (pf_model_mask[index_iz + ix + 1] == 0))
            {
                index_model_right = ix;
                break;
            }
        }
        for (ix=0; ix<index_model_left; ix++)
        {
            pf_model[index_iz + ix] = pf_model[index_iz + index_model_left];
        }
        for (ix=index_model_right+1; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_model[index_iz + ix] = pf_model[index_iz + index_model_right];
        }
    }
    for (iz=p_gbparams->stc_nz; iz<st_size_z; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<index_model_left; ix++)
        {
            pf_model[index_iz + ix] = pf_model[index_iz + index_model_left];
        }
        for (ix=index_model_right+1; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_model[index_iz + ix] = pf_model[index_iz + index_model_right];
        }
    }
}

void Apply_illummatrix2model(const globalconsts *p_gbparams, float *pf_dmodel, const float *pf_illummatrix, const unsigned long st_size_z)
{
    unsigned long iz;
    unsigned long ix;
    unsigned long index_iz;
    float eps_illummatrix = 0.05 * fabs(pf_illummatrix[cblas_isamax(st_size_z * p_gbparams->stc_nx2, pf_illummatrix, 1)]);
    for (iz=0; iz<st_size_z; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_dmodel[index_iz + ix] *= (1.0 / (pf_illummatrix[index_iz + ix] + eps_illummatrix));
        }
    }
}

float Get_J(const globalconsts *p_gbparams, const fcomp *pfcomp_Res_Nslice
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE)
{
    float f_J = 0.0;
    float f_J_full = 0.0;
    unsigned long ifreq;
    unsigned long isrcx;
    unsigned long index_slice0;

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            f_J += (pfcomp_Res_Nslice[index_slice0 + isrcx].real * pfcomp_Res_Nslice[index_slice0 + isrcx].real
                    + pfcomp_Res_Nslice[index_slice0 + isrcx].imag * pfcomp_Res_Nslice[index_slice0 + isrcx].imag);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&f_J, &f_J_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return f_J_full;
}

float Get_J_pershot(const globalconsts *p_gbparams, const fcomp *pfcomp_Res_Nslice
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long isrc)
{
    float f_J = 0.0;
    float f_J_full = 0.0;
    unsigned long ifreq;
    unsigned long ix;
    unsigned long index_isrcx;
    unsigned long index_slice0;

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            index_isrcx = isrc * p_gbparams->stc_nx2 + ix;
            f_J += (pfcomp_Res_Nslice[index_slice0 + index_isrcx].real * pfcomp_Res_Nslice[index_slice0 + index_isrcx].real
                    + pfcomp_Res_Nslice[index_slice0 + index_isrcx].imag * pfcomp_Res_Nslice[index_slice0 + index_isrcx].imag);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&f_J, &f_J_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    return f_J_full;
}

void Get_dPmin0_image_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice
        , fcomp *pfcomp_dPmin0_slice_tmp
        , const fcomp *pfcomp_Wx_slice, const float *pf_WxGxtaper
        , const unsigned long *pst_i_vel, const float *pf_image, const float *pf_dimage)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;


    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_dPmin0_slice_tmp_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == 0
                 *   dQplus(0) = - dR .* Pmin(0);
                 *else
                 *   dQplus(iz) = (1 + R) .* dPplus(iz) - dR .* Pmin(iz);
                 *end
                 *left: pfcomp_tmp is dQplus, right: pfcomp_tmp is dPplus
                 */
                if (iz == 0)
                {
                    p_tmpvars->pfcomp_tmp[index_isrc + ix].real
                        = - pf_dimage[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp[index_isrc + ix].imag
                        = - pf_dimage[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp[index_isrc + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp[index_isrc + ix].real
                        - pf_dimage[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp[index_isrc + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp[index_isrc + ix].imag
                        - pf_dimage[index_iz_m + ix] * pfcomp_Pmin_slice[index_isrc + ix].imag;
                }
            }
        }

        /*
         *calculate dPplus(iz+1) = Wx * dQplus(iz)
         *left: pfcomp_tmp is dPplus, right: pfcomp_tmp is dQplus
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp + index_izplus1, p_gbparams->stc_nx2);
    }

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == nz
                 *   dQmin(nz) = R .* dPplus(nz) + dR .* Pplus(nz);
                 *else
                 *   dQmin(iz) = (1 - R) .* dPmin(iz) + R .* dPplus(iz) + dR .* Pplus(iz);
                 *end
                 *left: pfcomp_dPmin0_slice_tmp is dQmin. right: pfcomp_dPmin0_slice_tmp is dPmin0
                 *pfcomp_tmp is dPplus
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp[index_isrc + ix].real
                        + pf_dimage[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                    pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp[index_isrc + ix].imag
                        + pf_dimage[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
                }
                else
                {
                    pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp[index_isrc + ix].real
                        + pf_dimage[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                    pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp[index_isrc + ix].imag
                        + pf_dimage[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
                }
            }
        }

        /*
         *calculate dPmin(iz-1) = Wx * dQmin(iz)
         *left: pfcomp_dPmin0_slice_tmp is dPmin. right: pfcomp_dPmin0_slice_tmp is dQmin
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_dPmin0_slice_tmp, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_dPmin0_slice_tmp_tmp, p_gbparams->stc_nx2);
        cblas_ccopy(p_gbparams->stc_nsrcnx2, p_tmpvars->pfcomp_dPmin0_slice_tmp_tmp, 1, pfcomp_dPmin0_slice_tmp, 1);
    }
}

void Get_dPmin0_image(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE)
{
    unsigned long ifreq;
    unsigned long isrcx;

    unsigned long index_slice;
    unsigned long index_slice0;
    unsigned long index_slice_o;
    float f_numerator = 0.0;
    float f_denominator = 0.0;
    *pf_numerator_full = 0.0;
    *pf_denominator_full = 0.0;

    fcomp *pfcomp_dPmin0_slice_tmp = NULL;
    pfcomp_dPmin0_slice_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        memset(pfcomp_dPmin0_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;

        Get_dPmin0_image_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice
                , pfcomp_dPmin0_slice_tmp
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , p_gbvars->pst_i_vel, p_gbvars->pf_image, p_gbvars->pf_dimage);

        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            //numerator and denominator for alpha
            f_numerator += (pfcomp_dPmin0_slice_tmp[isrcx].real * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].real
                    + pfcomp_dPmin0_slice_tmp[isrcx].imag * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].imag);
            f_denominator += (pfcomp_dPmin0_slice_tmp[isrcx].real * pfcomp_dPmin0_slice_tmp[isrcx].real
                    + pfcomp_dPmin0_slice_tmp[isrcx].imag * pfcomp_dPmin0_slice_tmp[isrcx].imag);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&f_numerator, pf_numerator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&f_denominator, pf_denominator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pfcomp_dPmin0_slice_tmp);
}

void Get_dPmin0(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE)
{
    unsigned long ifreq;
    unsigned long iorder_inner;
    unsigned long isrcx;

    unsigned long index_slice;
    unsigned long index_slice0;
    unsigned long index_slice_o;
    unsigned long index_ifreq;
    float f_numerator = 0.0;
    float f_denominator = 0.0;
    *pf_numerator_full = 0.0;
    *pf_denominator_full = 0.0;

    fcomp *pfcomp_Pmin_slice_tmp = NULL;
    pfcomp_Pmin_slice_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        memset(pfcomp_Pmin_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1
            * p_gbparams->stc_nsrcnx2);
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;
        index_ifreq = ifreq * p_gbparams->stc_nsrcnx2;

        for (iorder_inner=0; iorder_inner<iorder+1; iorder_inner++)
        {
            Get_Pplus_Pmin_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, pfcomp_Pmin_slice_tmp
                    , p_gbvars->pfcomp_Qplus_Nslice + index_slice
                    , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                    , p_gbvars->pfcomp_src + index_ifreq, p_gbvars->pst_i_vel, p_gbvars->pf_image);
        }

        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            pfcomp_Pmin_slice_tmp[isrcx].real -= p_gbvars->pfcomp_Pmin_Nslice[index_slice + isrcx].real;
            pfcomp_Pmin_slice_tmp[isrcx].imag -= p_gbvars->pfcomp_Pmin_Nslice[index_slice + isrcx].imag;
            //numerator and denominator for alpha
            f_numerator += (pfcomp_Pmin_slice_tmp[isrcx].real * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].real
                    + pfcomp_Pmin_slice_tmp[isrcx].imag * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].imag);
            f_denominator += (pfcomp_Pmin_slice_tmp[isrcx].real * pfcomp_Pmin_slice_tmp[isrcx].real
                    + pfcomp_Pmin_slice_tmp[isrcx].imag * pfcomp_Pmin_slice_tmp[isrcx].imag);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&f_numerator, pf_numerator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&f_denominator, pf_denominator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pfcomp_Pmin_slice_tmp);
}


void Get_image_tmp_withalpha(const globalconsts *p_gbparams, float *pf_image_tmp, const float *pf_image, const float *pf_dimage, const float f_alpha)
{
    unsigned long izx;
    for (izx=0; izx<p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2; izx++)
    {
        pf_image_tmp[izx] = pf_image[izx] + f_alpha * pf_dimage[izx];
    }
}
void Get_slowness_gradient_sub_output_pershot(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , float *pf_dslow_one_allshot, float *pf_illummatrix_velocity_one
        , const unsigned long *pst_i_vel, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_iz_m_allshot = 0;
    unsigned long index_izplus1_w = 0;
    unsigned long index_izplus1_m = 0;
    unsigned long index_izplus1_s = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_izminus1_m_allshot = 0;
    unsigned long index_izminus1_s = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long isrcx;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;


    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Gx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Qmin_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp1, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp2, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp3, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp3_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp4, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_w = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_s = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_iz_m_allshot = iz * p_gbparams->stc_nsrcnx2;
        index_izplus1_m = (iz + 1) * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        /*
         *Qmin(iz+1) = (1 - R) .* Pmin(iz+1) + R .* Pplus(iz+1);
         */
        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_izplus1_w + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].real
                    = (1.0 - pf_image[index_izplus1_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].real
                    + pf_image[index_izplus1_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].imag
                    = (1.0 - pf_image[index_izplus1_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].imag
                    + pf_image[index_izplus1_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == 0
                 *   Pmin_unified(0) = eye(nx);
                 *else
                 *   Pmin_unified(iz) = (1 - R) .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp1 on the left is Pmin_unified
                 *pfcomp_tmp1 on the right is Qmin_unified from previous iteration
                 */
                if (iz == 0)
                {
                    if (isrc == ix)
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 1.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                    else
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 0.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                }
                else
                {
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*
         *vel_gradient(iz) = Gx' * corr(Pmin_Res(iz), Qmin(iz+1, :, :))
         *                   = corr(Gx' * Pmin_unified(iz) * Res, Qmin(iz+1, :, :));
         *illum_matrix(iz) = ||Gx' * Pmin_unified(iz)||^2 .* ||Qmin(iz+1))||^2;
         */

        /*
         *calculate tmp_unified(iz) = Gx' * Pmin_unified(iz)
         *pfcomp_tmp4 on the left is tmp_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2);

        /*
         *tmp_Res(iz) = tmp_unified(iz) * Res;
         *pfcomp_tmp4 is tmp_unified
         *pfcomp_tmp2 is tmp_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *vel_gradient(iz) = corr(tmp_Res(iz), Qmin(iz+1, :, :));
         *illum_matrix(iz) = ||tmp_unified(iz)||^2 .* ||Qmin(iz+1))||^2;
         *pfcomp_tmp2 is Qmin_Res
         */
        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            pf_dslow_one_allshot[index_iz_m_allshot + isrcx] = p_tmpvars->pfcomp_tmp2[isrcx].real * p_tmpvars->pfcomp_Qmin_slice_tmp[isrcx].real
                + p_tmpvars->pfcomp_tmp2[isrcx].imag * p_tmpvars->pfcomp_Qmin_slice_tmp[isrcx].imag;
        }

        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_illummatrix_velocity_one[index_iz_m + ix]
                = (powf(cblas_scnrm2(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_Qmin_slice_tmp + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp4 + ix, p_gbparams->stc_nx2), 2));
        }

        /*
         *calculate Qmin_unified(iz+1) = Wx' * Pmin_unified(iz)
         *pfcomp_tmp1 on the left is Qmin_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2);
    }

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1_s = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;
        index_izminus1_m_allshot = (iz - 1) * p_gbparams->stc_nsrcnx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == nz
                 *   Pplus_unified(nz) = R .* Qmin_unified(nz);
                 *else
                 *   Pplus_unified(iz) = (1 + R) .* Qplus_unified(iz) + R .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp3 on the left is Pplus_unified currently
                 *pfcomp_tmp3 on the right is Qplus_unified from previous iteration
                 *pfcomp_tmp1 is Qmin_unified from previous code
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*
         *vel_gradient(iz-1) = Gx' * corr(Pplus_Res(iz), Qplus(iz-1, :, :))
         *                   = corr(Gx' * Pplus_unified(iz) * Res, Qplus(iz-1, :, :));
         *illum_matrix(iz-1) = ||Gx' * Pplus_unified(iz)|| .* ||Qplus(iz-1, :, :))||;
         */

        /*
         *calculate tmp_unified(iz) = Gx' * Pplus_unified(iz)
         *pfcomp_tmp4 on the left is tmp_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2);

        /*
         *tmp_Res(iz) = tmp_unified(iz) * Res;
         *pfcomp_tmp4 is tmp_unified
         *pfcomp_tmp2 is tmp_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *vel_gradient(iz-1) = corr(tmp_Res(iz), Qplus(iz-1, :, :));
         *illum_matrix(iz-1) = ||tmp_unified(iz)|| .* ||Qplus(iz-1, :, :))||;
         *pfcomp_tmp2 is tmp_Res
         *pfcomp_tmp4 is tmp_unified
         */
        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            pf_dslow_one_allshot[index_izminus1_m_allshot + isrcx] += (p_tmpvars->pfcomp_tmp2[isrcx].real * pfcomp_Qplus_slice[index_izminus1_s + isrcx].real
                + p_tmpvars->pfcomp_tmp2[isrcx].imag * pfcomp_Qplus_slice[index_izminus1_s + isrcx].imag);
        }
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            pf_illummatrix_velocity_one[index_izminus1_m + ix]
                += (powf(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Qplus_slice + index_izminus1_s + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp4 + ix, p_gbparams->stc_nx2), 2));
        }

        /*
         *calculate Qplus_unified(iz-1) = Wx' * Pplus_unified(iz)
         *pfcomp_tmp3 on the left is Qplus_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp3_tmp, p_gbparams->stc_nx2);
        cblas_ccopy(p_gbparams->stc_nx2 * p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3_tmp, 1, p_tmpvars->pfcomp_tmp3, 1);
    }

}

void Get_slowness_gradient_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , const fcomp *pfcomp_Res_slice
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , float *pf_dslow_one, float *pf_illummatrix_velocity_one
        , const unsigned long *pst_i_vel, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izplus1_w = 0;
    unsigned long index_izplus1_m = 0;
    unsigned long index_izplus1_s = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_izminus1_s = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;


    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Gx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Qmin_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp1, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp2, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp3, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp3_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_tmp4, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);


    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_w = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izplus1_s = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izplus1_m = (iz + 1) * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        /*
         *Qmin(iz+1) = (1 - R) .* Pmin(iz+1) + R .* Pplus(iz+1);
         */
        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_izplus1_w + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].real
                    = (1.0 - pf_image[index_izplus1_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].real
                    + pf_image[index_izplus1_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].imag
                    = (1.0 - pf_image[index_izplus1_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].imag
                    + pf_image[index_izplus1_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*if iz == 0
                 *   Pmin_unified(0) = eye(nx);
                 *else
                 *   Pmin_unified(iz) = (1 - R) .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp1 on the left is Pmin_unified
                 *pfcomp_tmp1 on the right is Qmin_unified from previous iteration
                 */
                if (iz == 0)
                {
                    if (isrc == ix)
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 1.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                    else
                    {
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].real = 0.0;
                        p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag = 0.0;
                    }
                }
                else
                {
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*
         *vel_gradient(iz) = Gx' * corr(Pmin_Res(iz), Qmin(iz+1, :, :))
         *                   = corr(Gx' * Pmin_unified(iz) * Res, Qmin(iz+1, :, :));
         *illum_matrix(iz) = ||Gx' * Pmin_unified(iz)||^2 .* ||Qmin(iz+1))||^2;
         */

        /*
         *calculate tmp_unified(iz) = Gx' * Pmin_unified(iz)
         *pfcomp_tmp4 on the left is tmp_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2);

        /*
         *tmp_Res(iz) = tmp_unified(iz) * Res;
         *pfcomp_tmp4 is tmp_unified
         *pfcomp_tmp2 is tmp_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *vel_gradient(iz) = corr(tmp_Res(iz), Qmin(iz+1, :, :));
         *illum_matrix(iz) = ||tmp_unified(iz)||^2 .* ||Qmin(iz+1))||^2;
         *pfcomp_tmp2 is Qmin_Res
         */
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            cblas_cdotc_sub(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_tmp2 + ix, p_gbparams->stc_nx2
                    , p_tmpvars->pfcomp_Qmin_slice_tmp + ix, p_gbparams->stc_nx2
                    , &fcomp_tmp);
            pf_dslow_one[index_iz_m + ix] = fcomp_tmp.real;
            pf_illummatrix_velocity_one[index_iz_m + ix]
                = (powf(cblas_scnrm2(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_Qmin_slice_tmp + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp4 + ix, p_gbparams->stc_nx2), 2));
        }

        /*
         *calculate Qmin_unified(iz+1) = Wx' * Pmin_unified(iz)
         *pfcomp_tmp1 on the left is Qmin_unified
         *pfcomp_tmp1 on the right is Pmin_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp1 + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp1 + index_izplus1, p_gbparams->stc_nx2);
    }

    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nx2 * p_gbparams->stc_nx2;
        index_izminus1_s = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nx2; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == nz
                 *   Pplus_unified(nz) = R .* Qmin_unified(nz);
                 *else
                 *   Pplus_unified(iz) = (1 + R) .* Qplus_unified(iz) + R .* Qmin_unified(iz);
                 *end
                 *pfcomp_tmp3 on the left is Pplus_unified currently
                 *pfcomp_tmp3 on the right is Qplus_unified from previous iteration
                 *pfcomp_tmp1 is Qmin_unified from previous code
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp3[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp1[index_isrc + ix].imag;
                }
            }
        }

        /*
         *vel_gradient(iz-1) = Gx' * corr(Pplus_Res(iz), Qplus(iz-1, :, :))
         *                   = corr(Gx' * Pplus_unified(iz) * Res, Qplus(iz-1, :, :));
         *illum_matrix(iz-1) = ||Gx' * Pplus_unified(iz)|| .* ||Qplus(iz-1, :, :))||;
         */

        /*
         *calculate tmp_unified(iz) = Gx' * Pplus_unified(iz)
         *pfcomp_tmp4 on the left is tmp_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2);

        /*
         *tmp_Res(iz) = tmp_unified(iz) * Res;
         *pfcomp_tmp4 is tmp_unified
         *pfcomp_tmp2 is tmp_Res
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Res_slice, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_tmp4, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2);

        /*
         *vel_gradient(iz-1) = corr(tmp_Res(iz), Qplus(iz-1, :, :));
         *illum_matrix(iz-1) = ||tmp_unified(iz)|| .* ||Qplus(iz-1, :, :))||;
         *pfcomp_tmp2 is tmp_Res
         *pfcomp_tmp4 is tmp_unified
         */
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            cblas_cdotc_sub(p_gbparams->stc_nsrc, p_tmpvars->pfcomp_tmp2 + ix, p_gbparams->stc_nx2
                    , pfcomp_Qplus_slice + index_izminus1_s + ix, p_gbparams->stc_nx2
                    , &fcomp_tmp);
            pf_dslow_one[index_izminus1_m + ix] += fcomp_tmp.real;
            pf_illummatrix_velocity_one[index_izminus1_m + ix]
                += (powf(cblas_scnrm2(p_gbparams->stc_nsrc, pfcomp_Qplus_slice + index_izminus1_s + ix, p_gbparams->stc_nx2), 2)
                        * powf(cblas_scnrm2(p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp4 + ix, p_gbparams->stc_nx2), 2));
        }

        /*
         *calculate Qplus_unified(iz-1) = Wx' * Pplus_unified(iz)
         *pfcomp_tmp3 on the left is Qplus_unified
         *pfcomp_tmp3 on the right is Pplus_unified
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans
                , p_gbparams->stc_nx2, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp3, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp3_tmp, p_gbparams->stc_nx2);
        cblas_ccopy(p_gbparams->stc_nx2 * p_gbparams->stc_nx2, p_tmpvars->pfcomp_tmp3_tmp, 1, p_tmpvars->pfcomp_tmp3, 1);
    }
}

void Get_slowness_gradient_output_pershot(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter)
{
    //pf_slow doesn't touched in this function, which acts as slowness from last iteration
    //Instead, pst_i_vel is updated
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_m = 0;
    unsigned long index_slice_m_allshot = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long ismooth;
    unsigned long ix;
    unsigned long isrc;
    unsigned long iz;

    float *pf_dslow_pershot = NULL;
    float *pf_illummatrix_velocity = NULL;
    float *pf_dslow_allshot = NULL;
    float *pf_dslow_Nslice_allshot = NULL;
    float *pf_dslow_Nslice_allshot_full = NULL;
    float *pf_illummatrix_velocity_Nslice = NULL;
    float *pf_illummatrix_velocity_Nslice_full = NULL;
    pf_dslow_Nslice_allshot = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2);
    pf_illummatrix_velocity_Nslice = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    if (MPI_ID == 0)
    {
        pf_dslow_pershot = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        pf_dslow_allshot = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_velocity = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        pf_dslow_Nslice_allshot_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2);
        pf_illummatrix_velocity_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_m = (ifreq / MPI_SIZE) * p_gbparams->stc_nz * p_gbparams->stc_nx2;
        index_slice_m_allshot = (ifreq / MPI_SIZE) * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;
        Get_slowness_gradient_sub_output_pershot(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice, p_gbvars->pfcomp_Qplus_Nslice + index_slice
                , p_gbvars->pfcomp_Res_Nslice + index_slice0
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pfcomp_Gx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , pf_dslow_Nslice_allshot + index_slice_m_allshot, pf_illummatrix_velocity_Nslice + index_slice_m
                , p_gbvars->pst_i_vel, p_gbvars->pf_image);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(pf_dslow_Nslice_allshot, pf_dslow_Nslice_allshot_full
            , Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pf_illummatrix_velocity_Nslice, pf_illummatrix_velocity_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    FILE *fp;
    char filename[64];

    //combine velocity gradient from different frequency and apply illumination matrix to it and scale it
    if (MPI_ID == 0)
    {
        Combine_slowness_gradient_output_pershot(p_gbparams, pf_dslow_allshot, pf_illummatrix_velocity, pf_dslow_Nslice_allshot_full, pf_illummatrix_velocity_Nslice_full, Nslice_tmp);

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            for (iz=0; iz<p_gbparams->stc_nz; iz++)
            {
                cblas_scopy(p_gbparams->stc_nx2, pf_dslow_allshot + iz * p_gbparams->stc_nsrcnx2 + isrc * p_gbparams->stc_nx2, 1
                        , pf_dslow_pershot + iz * p_gbparams->stc_nx2, 1);
            }
            cblas_saxpy(p_gbparams->stc_nz * p_gbparams->stc_nx2, 1.0
                    , pf_dslow_pershot, 1, p_gbvars->pf_dslow, 1);

            if (p_gbparams->fc_fmax_tmp < 12.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_dslow_pershot, p_gbparams->stc_nz, p_gbvars->pf_vel);
            }
            else if (p_gbparams->fc_fmax_tmp < 20.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_dslow_pershot, p_gbparams->stc_nz, p_gbvars->pf_vel);
            }
            else if (p_gbparams->fc_fmax_tmp < 30.0)
            {
                Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_dslow_pershot, p_gbparams->stc_nz, p_gbvars->pf_vel);
            }
            else
            {
                Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                        , pf_dslow_pershot, p_gbparams->stc_nz, p_gbvars->pf_vel);
            }

            segy seg_tmp;
            memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);
            snprintf(filename, sizeof(char) * 64, "%s%s_dslowness_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
            seg_tmp.trwf = (int)p_gbparams->stc_nx;
            seg_tmp.ns = (int)p_gbparams->stc_nz;
            seg_tmp.d2 = p_gbparams->fc_dx;
            seg_tmp.timbas = 3; //sort on keyword fldr
            seg_tmp.d1 = p_gbparams->fc_dz;
            seg_tmp.fldr = (int)(isrc + 1);
            seg_tmp.duse = (int)(iter + 1);
            seg_tmp.trid = 30;
            fp= fopen(filename, "a+");
            if (fp == NULL)
            {
                syserr("error: cannot open file %s to write P file.\n", filename);
            }

            for (ix=0; ix<p_gbparams->stc_nx; ix++)
            {
                cblas_scopy(p_gbparams->stc_nz, pf_dslow_pershot + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                        , seg_tmp.data, 1);
                seg_tmp.tracl = (int)(isrc * p_gbparams->stc_nx + ix + 1);
                seg_tmp.tracf = (int)(ix + 1);
                fwrite(&seg_tmp, sizeof(float)
                        , HDRSIZE / sizeof(float) + p_gbparams->stc_nz, fp);
            }
            fclose(fp);
        }

        if (p_gbparams->fc_fmax_tmp < 12.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 20.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 30.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else
        {
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }

	// Apply the model mask to the velocity gradient only for option 2 or 3
        if (p_gbparams->stc_if_model_mask >= 2)
        {
            //It also masked the illumination matrix, but this part is maybe not useful
	    //printf("Apply model mask to illum matrix velocity\n");
            //Apply_model_mask(p_gbparams, pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_model_mask);
	    //printf("Apply model mask to velocity\n");
            Apply_model_mask(p_gbparams, p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_model_mask);
        }

        for (ismooth=0; ismooth<p_gbparams->stc_illum_smoothN; ismooth++)
        {
            Smooth_model_2d(pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbparams->stc_nx2
                    , p_gbparams->stc_illum_smoothnz, p_gbparams->stc_illum_smoothnx);
        }
        for (ismooth=0; ismooth<p_gbparams->stc_vel_smoothN; ismooth++)
        {
            Smooth_model_2d(p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbparams->stc_nx2
                    , p_gbparams->stc_vel_smoothnz, p_gbparams->stc_vel_smoothnx);
        }
        Apply_illummatrix2model(p_gbparams, p_gbvars->pf_dslow, pf_illummatrix_velocity, p_gbparams->stc_nz);

        Scale_gradient(p_gbparams, p_gbvars->pf_dslow, p_gbvars->pf_slow, p_gbparams->stc_nz);

        if (fabs(p_gbparams->fc_reflconstr_w) > EPS)
        {
            Get_extraslownessgradient_reflectivityconst(p_gbparams, p_gbvars, iter);
        }

        #ifdef DEBUG
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_gbvars->pf_dslow, p_gbparams->stc_nz * p_gbparams->stc_nx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pf_dslow_Nslice_allshot);
    free(pf_illummatrix_velocity_Nslice);
    if (MPI_ID == 0)
    {
        free(pf_dslow_pershot);
        free(pf_dslow_allshot);
        free(pf_illummatrix_velocity);
        free(pf_dslow_Nslice_allshot_full);
        free(pf_illummatrix_velocity_Nslice_full);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Get_slowness_gradient(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp, const unsigned long iter)
{
    //pf_slow doesn't touched in this function, which acts as slowness from last iteration
    //Instead, pst_i_vel is updated
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_m = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long ismooth;

    float *pf_illummatrix_velocity = NULL;
    float *pf_dslow_Nslice = NULL;
    float *pf_dslow_Nslice_full = NULL;
    float *pf_illummatrix_velocity_Nslice = NULL;
    float *pf_illummatrix_velocity_Nslice_full = NULL;
    pf_dslow_Nslice= Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    pf_illummatrix_velocity_Nslice = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    if (MPI_ID == 0)
    {
        pf_illummatrix_velocity = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        pf_dslow_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
        pf_illummatrix_velocity_Nslice_full = Allocate_1d_float(Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_m = (ifreq / MPI_SIZE) * p_gbparams->stc_nz * p_gbparams->stc_nx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;
        Get_slowness_gradient_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice, p_gbvars->pfcomp_Qplus_Nslice + index_slice
                , p_gbvars->pfcomp_Res_Nslice + index_slice0
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pfcomp_Gx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , pf_dslow_Nslice + index_slice_m, pf_illummatrix_velocity_Nslice + index_slice_m
                , p_gbvars->pst_i_vel, p_gbvars->pf_image);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(pf_dslow_Nslice, pf_dslow_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(pf_illummatrix_velocity_Nslice, pf_illummatrix_velocity_Nslice_full
            , Nslice_tmp * p_gbparams->stc_nz * p_gbparams->stc_nx2
            , MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //combine velocity gradient from different frequency and apply illumination matrix to it and scale it
    if (MPI_ID == 0)
    {
        Combine_slowness_gradient(p_gbparams, p_gbvars->pf_dslow, pf_illummatrix_velocity, pf_dslow_Nslice_full, pf_illummatrix_velocity_Nslice_full, Nslice_tmp);

        if (p_gbparams->fc_fmax_tmp < 12.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.3 * p_gbparams->fc_fmax_tmp, 0.2 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 20.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.4 * p_gbparams->fc_fmax_tmp, 0.3 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else if (p_gbparams->fc_fmax_tmp < 30.0)
        {
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.5 * p_gbparams->fc_fmax_tmp, 0.4 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }
        else
        {
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_vel);
            Apply_kxkzfilter(p_gbparams, 0.6 * p_gbparams->fc_fmax_tmp, 0.5 * p_gbparams->fc_dx, p_gbparams->fc_dz
                    , p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_vel);
        }

	// Apply the model mask to the velocity gradient only for option 2 or 3
        if (p_gbparams->stc_if_model_mask >= 2)
        {
            //It also masked the illumination matrix, but this part is maybe not useful
	    //printf("Apply model mask to illum matrix velocity\n");
            //Apply_model_mask(p_gbparams, pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbvars->pf_model_mask);
	    //printf("Apply model mask to velocity\n");
            Apply_model_mask(p_gbparams, p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbvars->pf_model_mask);
        }

        for (ismooth=0; ismooth<p_gbparams->stc_illum_smoothN; ismooth++)
        {
            Smooth_model_2d(pf_illummatrix_velocity, p_gbparams->stc_nz, p_gbparams->stc_nx2
                    , p_gbparams->stc_illum_smoothnz, p_gbparams->stc_illum_smoothnx);
        }
        for (ismooth=0; ismooth<p_gbparams->stc_vel_smoothN; ismooth++)
        {
            Smooth_model_2d(p_gbvars->pf_dslow, p_gbparams->stc_nz, p_gbparams->stc_nx2
                    , p_gbparams->stc_vel_smoothnz, p_gbparams->stc_vel_smoothnx);
        }
        Apply_illummatrix2model(p_gbparams, p_gbvars->pf_dslow, pf_illummatrix_velocity, p_gbparams->stc_nz);

        Scale_gradient(p_gbparams, p_gbvars->pf_dslow, p_gbvars->pf_slow, p_gbparams->stc_nz);

        if (fabs(p_gbparams->fc_reflconstr_w) > EPS)
        {
            Get_extraslownessgradient_reflectivityconst(p_gbparams, p_gbvars, iter);
        }

        #ifdef DEBUG
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(p_gbvars->pf_dslow, p_gbparams->stc_nz * p_gbparams->stc_nx2
            , MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pf_dslow_Nslice);
    free(pf_illummatrix_velocity_Nslice);
    if (MPI_ID == 0)
    {
        free(pf_illummatrix_velocity);
        free(pf_dslow_Nslice_full);
        free(pf_illummatrix_velocity_Nslice_full);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Get_extraslownessgradient_reflectivityconst(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long iter)
{
    /*add extra constraint to slowness gradient by*/
    /*- differentiate current velocity in depth*/
    /*- fit this to current reflectivity*/
    /*- integrate the mismatch and low cut filter*/
    /*this is a some what pargmatic step*/
    FILE *fp;
    char filename[64];

    float f_scale_factor1 = 0.0;
    float f_scale_factor2 = 0.0;
    unsigned long ix;
    unsigned long iz;
    float *pf_image_fromvelocity_transpose = NULL;
    float *pf_image_diff_transpose = NULL;
    float *pf_image_diff_sm_transpose = NULL;
    float *pf_smooth_filter = NULL;
    fcomp *pfcomp_tmp = NULL;
    pf_image_fromvelocity_transpose = Allocate_1d_float(p_gbparams->stc_nx2 * p_gbparams->stc_nzplus1);
    pf_image_diff_transpose = Allocate_1d_float(p_gbparams->stc_nx2 * p_gbparams->stc_nzplus1);
    pf_image_diff_sm_transpose = Allocate_1d_float(p_gbparams->stc_nx2 * p_gbparams->stc_nzplus1);
    pf_smooth_filter = Allocate_1d_float(3);
    pf_smooth_filter[0] = 0.33;
    pf_smooth_filter[1] = 0.34;
    pf_smooth_filter[2] = 0.33;
    pfcomp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nzplus1);

    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        cblas_scopy(p_gbparams->stc_nzplus1, p_gbvars->pf_image + ix, p_gbparams->stc_nx2
                , pf_image_diff_transpose + ix * p_gbparams->stc_nzplus1, 1);
    }

    // calculate derivative from current velocity model in depth
    // pf_image_diff_transpose = pf_image currently
    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        for (iz=1; iz<p_gbparams->stc_nz; iz++)
        {
            pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz]
                = p_gbvars->pf_vel[iz * p_gbparams->stc_nx2 + ix]
                - p_gbvars->pf_vel[(iz - 1) * p_gbparams->stc_nx2 + ix];
            f_scale_factor1 += (pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz] * pf_image_diff_transpose[ix * p_gbparams->stc_nzplus1 + iz]);
            f_scale_factor2 += (pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz] * pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz]);
        }
    }
    cblas_scopy(p_gbparams->stc_nx2, pf_image_fromvelocity_transpose + 1, p_gbparams->stc_nzplus1
            , pf_image_fromvelocity_transpose + 0, p_gbparams->stc_nzplus1);
    cblas_scopy(p_gbparams->stc_nx2, pf_image_fromvelocity_transpose + p_gbparams->stc_nz - 1, p_gbparams->stc_nzplus1
            , pf_image_fromvelocity_transpose + p_gbparams->stc_nz, p_gbparams->stc_nzplus1);
    f_scale_factor1 /= f_scale_factor2;
    cblas_sscal(p_gbparams->stc_nx2 * p_gbparams->stc_nzplus1, f_scale_factor1, pf_image_fromvelocity_transpose, 1);

    // define LS scale factor to fit pf_image_fromvelocity_transpose to pf_image
    // pf_image_diff_transpose = pf_image currently
    f_scale_factor1 = 0.0;
    f_scale_factor2 = 0.0;
    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        for (iz=0; iz<p_gbparams->stc_nzplus1; iz++)
        {
            f_scale_factor1 += (pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz] * pf_image_diff_transpose[ix * p_gbparams->stc_nzplus1 + iz]);
            f_scale_factor2 += (pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz] * pf_image_fromvelocity_transpose[ix * p_gbparams->stc_nzplus1 + iz]);
        }
    }
    f_scale_factor1 /= f_scale_factor2;
    cblas_saxpy(p_gbparams->stc_nx2 * p_gbparams->stc_nzplus1, -1.0 * f_scale_factor1, pf_image_fromvelocity_transpose, 1, pf_image_diff_transpose, 1);

    // integration, make sure to include the proper dz scaling
    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        for (iz=1; iz<p_gbparams->stc_nzplus1; iz++)
        {
            pf_image_diff_transpose[ix * p_gbparams->stc_nzplus1 + iz]
                = pf_image_diff_transpose[ix * p_gbparams->stc_nzplus1 + iz - 1]
                + p_gbparams->fc_dz * pf_image_diff_transpose[ix * p_gbparams->stc_nzplus1 + iz];
        }
    }

    // smoothing in z
    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        conv_jmi(pf_image_diff_transpose + ix * p_gbparams->stc_nzplus1, pf_smooth_filter, pf_image_diff_sm_transpose + ix * p_gbparams->stc_nzplus1, p_gbparams->stc_nzplus1, 3);
    }
    cblas_scopy(p_gbparams->stc_nx2, pf_image_diff_sm_transpose + 1, p_gbparams->stc_nzplus1
            , pf_image_diff_sm_transpose + 0, p_gbparams->stc_nzplus1);
    cblas_scopy(p_gbparams->stc_nx2, pf_image_diff_sm_transpose + p_gbparams->stc_nz - 1, p_gbparams->stc_nzplus1
            , pf_image_diff_sm_transpose + p_gbparams->stc_nz, p_gbparams->stc_nzplus1);

    // remove low freqency trend
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nzplus1);
    status = DftiCommitDescriptor(my_handle);

    for (ix=0; ix<p_gbparams->stc_nx2; ix++)
    {
        for (iz=0; iz<p_gbparams->stc_nzplus1; iz++)
        {
            pfcomp_tmp[iz].real = pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + iz];
            pfcomp_tmp[iz].imag = 0.0;
        }
        status = DftiComputeForward(my_handle, pfcomp_tmp);

        pfcomp_tmp[0].real *= 0.0;
        pfcomp_tmp[0].imag *= 0.0;
        pfcomp_tmp[1].real *= 0.25;
        pfcomp_tmp[1].imag *= 0.25;
        pfcomp_tmp[2].real *= 0.5;
        pfcomp_tmp[2].imag *= 0.5;
        pfcomp_tmp[3].real *= 0.5;
        pfcomp_tmp[3].imag *= 0.5;
        pfcomp_tmp[p_gbparams->stc_nz].real *= 0.25;
        pfcomp_tmp[p_gbparams->stc_nz].imag *= 0.25;
        pfcomp_tmp[p_gbparams->stc_nz - 1].real *= 0.5;
        pfcomp_tmp[p_gbparams->stc_nz - 1].imag *= 0.5;
        pfcomp_tmp[p_gbparams->stc_nz - 2].real *= 0.5;
        pfcomp_tmp[p_gbparams->stc_nz - 2].imag *= 0.5;

        status = DftiComputeBackward(my_handle, pfcomp_tmp);
        for (iz=0; iz<p_gbparams->stc_nzplus1; iz++)
        {
            // scale it such that it has the values of normal slowness updates
            // (usual numbers in order of 1 for this constraint)
            pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + iz] = - pfcomp_tmp[iz].real / (float)p_gbparams->stc_nzplus1 / 400000.0;
        }

        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + 2] *= 0.0;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + 3] *= 0.25;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + 4] *= 0.5;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + 5] *= 0.75;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + p_gbparams->stc_nz] *= 0.0;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + p_gbparams->stc_nz - 1] *= 0.25;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + p_gbparams->stc_nz - 2] *= 0.5;
        pf_image_diff_sm_transpose[ix * p_gbparams->stc_nzplus1 + p_gbparams->stc_nz - 3] *= 0.75;
    }

    segy seg_tmp;
    memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);
    snprintf(filename, sizeof(char) * 64, "%s%s_extradslowness_fromRC.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    seg_tmp.trwf = (int)p_gbparams->stc_nx;
    seg_tmp.ns = (int)p_gbparams->stc_nz;
    seg_tmp.d2 = p_gbparams->fc_dx;
    seg_tmp.timbas = 11; //sort on keyword duse
    seg_tmp.d1 = p_gbparams->fc_dz;
    seg_tmp.duse = (int)(iter + 1);
    seg_tmp.trid = 30;
    fp= fopen(filename, "a+");
    if (fp == NULL)
    {
        syserr("error: cannot open file %s to write P file.\n", filename);
    }

    for (ix=0; ix<p_gbparams->stc_nx; ix++)
    {
        cblas_scopy(p_gbparams->stc_nz, pf_image_diff_sm_transpose + (ix + p_gbparams->stc_nxtap) * p_gbparams->stc_nzplus1 + 1, 1
                , seg_tmp.data, 1);
        seg_tmp.tracl = (int)(iter * p_gbparams->stc_nx + ix + 1);
        seg_tmp.tracf = (int)(ix + 1);
        fwrite(&seg_tmp, sizeof(float)
                , HDRSIZE / sizeof(float) + p_gbparams->stc_nz, fp);
    }
    fclose(fp);

    // also skip first depth sample, as velocities have one sample less
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        cblas_saxpy(p_gbparams->stc_nx2, p_gbparams->fc_reflconstr_w, pf_image_diff_sm_transpose + iz + 1, p_gbparams->stc_nzplus1
                , p_gbvars->pf_dslow + iz * p_gbparams->stc_nx2, 1);
    }

    free(pf_image_fromvelocity_transpose);
    free(pf_image_diff_transpose);
    free(pf_image_diff_sm_transpose);
    free(pf_smooth_filter);
    free(pfcomp_tmp);
}

void Apply_directionalTV(const globalconsts *p_gbparams, float fc_DTV_lam, float fc_DTV_mu, float fc_DTV_alp, const float *pf_dip_field, float *pf_slow, float *pf_dtv_dx_est, float *pf_dtv_dz_est, float *pf_dtv_bx_est, float *pf_dtv_bz_est)
{
    unsigned long izx;
    float *pf_tvx = NULL;
    float *pf_tvz = NULL;
    float *pf_dtv = NULL;
    pf_tvx = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
    pf_tvz = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
    pf_dtv = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);

    DTV_DDfilter(pf_tvx, pf_tvz, pf_slow, pf_dip_field, fc_DTV_alp, p_gbparams->stc_nx2, p_gbparams->stc_nz);

    for (izx=0; izx<p_gbparams->stc_nz * p_gbparams->stc_nx2; izx++)
    {
        pf_tvx[izx] = pf_dtv_dx_est[izx] - pf_tvx[izx] - pf_dtv_bx_est[izx];
        pf_tvz[izx] = pf_dtv_dz_est[izx] - pf_tvz[izx] - pf_dtv_bz_est[izx];
    }

    DTV_DDTfilter(pf_dtv, pf_tvx, pf_tvz, pf_slow, pf_dip_field, fc_DTV_alp, p_gbparams->stc_nx2, p_gbparams->stc_nz);
    for (izx=0; izx<p_gbparams->stc_nz * p_gbparams->stc_nx2; izx++)
    {
        pf_slow[izx] += 1.0 / fc_DTV_mu * fc_DTV_lam * pf_dtv[izx];
    }

    DTV_DDfilter(pf_tvx, pf_tvz, pf_slow, pf_dip_field, fc_DTV_alp, p_gbparams->stc_nx2, p_gbparams->stc_nz);
    for (izx=0; izx<p_gbparams->stc_nz * p_gbparams->stc_nx2; izx++)
    {
        if (fabs(pf_tvx[izx] + pf_dtv_bx_est[izx]) <= (1.0 / fc_DTV_lam))
        {
            pf_dtv_dx_est[izx] = 0.0;
        }
        else
        {
            pf_dtv_dx_est[izx] = (fabs(pf_tvx[izx] + pf_dtv_bx_est[izx]) - (1.0 / fc_DTV_lam))
                / fabs(pf_tvx[izx] + pf_dtv_bx_est[izx]) * (pf_tvx[izx] + pf_dtv_bx_est[izx]);
        }
        if (fabs(pf_tvz[izx] + pf_dtv_bz_est[izx]) <= (1.0 / fc_DTV_lam))
        {
            pf_dtv_dz_est[izx] = 0.0;
        }
        else
        {
            pf_dtv_dz_est[izx] = (fabs(pf_tvz[izx] + pf_dtv_bz_est[izx]) - (1.0 / fc_DTV_lam))
                / fabs(pf_tvz[izx] + pf_dtv_bz_est[izx]) * (pf_tvz[izx] + pf_dtv_bz_est[izx]);
        }
        pf_dtv_bx_est[izx] += (pf_tvx[izx] - pf_dtv_dx_est[izx]);
        pf_dtv_bz_est[izx] += (pf_tvz[izx] - pf_dtv_dz_est[izx]);
    }

    free(pf_tvx);
    free(pf_tvz);
    free(pf_dtv);
}

void Get_dPmin0_velocity_sub(const globalconsts *p_gbparams, tmparrays *p_tmpvars, const fcomp *pfcomp_Pplus_slice, const fcomp *pfcomp_Pmin_slice, const fcomp *pfcomp_Qplus_slice
        , fcomp *pfcomp_dPmin0_slice_tmp
        , const fcomp *pfcomp_Wx_slice, const fcomp *pfcomp_Gx_slice, const float *pf_WxGxtaper
        , const unsigned long *pst_i_vel, const float *pf_dslow, const float *pf_image)
{
    unsigned long index_isrc = 0;
    unsigned long index_iz = 0;
    unsigned long index_izminus1 = 0;
    unsigned long index_izplus1 = 0;
    unsigned long index_iz_m = 0;
    unsigned long index_izminus1_m = 0;
    unsigned long index_isrc_s = 0;
    unsigned long index_ix = 0;
    unsigned long ix_shifted = 0;
    unsigned long Wx_index = 0;
    unsigned long isrc;
    unsigned long iz;
    unsigned long ix;
    unsigned long ix_izplus1;
    unsigned long ix_izminus1;
    fcomp fcomp_alpha;
    fcomp_alpha.real = 1.0;
    fcomp_alpha.imag = 0.0;
    fcomp fcomp_beta;
    fcomp_beta.real = 0.0;
    fcomp_beta.imag = 0.0;

    memset(p_tmpvars->pfcomp_Wx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Gx_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nx2 * p_gbparams->stc_nx2);
    memset(p_tmpvars->pfcomp_Qmin_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp2, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp11, 0, sizeof(float) * 2 * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    memset(p_tmpvars->pfcomp_tmp33, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);

    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izplus1 = (iz + 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;

        for (ix_izplus1=0; ix_izplus1<p_gbparams->stc_nx2; ix_izplus1++)
        {
            index_ix = ix_izplus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izplus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_iz_m + ix_izplus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc = index_iz + isrc * p_gbparams->stc_nx2;
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == 0
                 *   dQplus(0) = 0;
                 *else
                 *   dQplus(iz) = (1 + R) .* dPplus(iz);
                 *end
                 *pfcomp_tmp2 is dQplus; pfcomp_tmp11 is dPplus
                 */
                if (iz == 0)
                {

                }
                else
                {
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].real
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp11[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].imag
                        = (1.0 + pf_image[index_iz_m + ix]) * p_tmpvars->pfcomp_tmp11[index_isrc + ix].imag;
                }
            }
        }

        /*
         *calculate dPplus(iz+1) = Wx * dQplus(iz) + Gx * diag(dslow) * Qplus(iz) = dPplus1 + dPplus2
         */
        /*
         *calculate dPplus1 = Wx * dQplus(iz)
         *pfcomp_tmp11 is dPplus1; pfcomp_tmp2 is dQplus(iz)
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp11 + index_izplus1, p_gbparams->stc_nx2);

        /*
         *calculate dPplus2 = diag(dslow) * Gx * Qplus(iz)
         *pfcomp_tmp33 is dPplus2;
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, pfcomp_Qplus_slice + index_iz, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp33, p_gbparams->stc_nx2);
        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                p_tmpvars->pfcomp_tmp11[index_izplus1 + index_isrc_s + ix].real += (pf_dslow[index_iz_m + ix] * p_tmpvars->pfcomp_tmp33[index_isrc_s + ix].real);
                p_tmpvars->pfcomp_tmp11[index_izplus1 + index_isrc_s + ix].imag += (pf_dslow[index_iz_m + ix] * p_tmpvars->pfcomp_tmp33[index_isrc_s + ix].imag);
            }
        }
    }


    for (iz=p_gbparams->stc_nz; iz>0; iz--)
    {
        index_iz = iz * p_gbparams->stc_nsrcnx2;
        index_izminus1 = (iz - 1) * p_gbparams->stc_nsrcnx2;
        index_iz_m = iz * p_gbparams->stc_nx2;
        index_izminus1_m = (iz - 1) * p_gbparams->stc_nx2;

        for (ix_izminus1=0; ix_izminus1<p_gbparams->stc_nx2; ix_izminus1++)
        {
            index_ix = ix_izminus1 * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                ix_shifted = (unsigned long)(((((int)ix - (int)ix_izminus1) % (int)p_gbparams->stc_nx2pow2 + (int)p_gbparams->stc_nx2pow2)) % (int)p_gbparams->stc_nx2pow2);
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].real
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Wx_tmp[index_ix + ix].imag
                    = pfcomp_Wx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].real
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].real
                    * pf_WxGxtaper[index_ix + ix];
                p_tmpvars->pfcomp_Gx_tmp[index_ix + ix].imag
                    = pfcomp_Gx_slice[pst_i_vel[index_izminus1_m + ix_izminus1] * p_gbparams->stc_nx2pow2 + ix_shifted].imag
                    * pf_WxGxtaper[index_ix + ix];
            }
        }

        /*
         *Qmin(iz) = (1 - R) .* Pmin(iz) + R .* Pplus(iz);
         */
        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].real
                    = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].real
                    + pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].real;
                p_tmpvars->pfcomp_Qmin_slice_tmp[index_isrc_s + ix].imag
                    = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_Pmin_slice[index_isrc + ix].imag
                    + pf_image[index_iz_m + ix] * pfcomp_Pplus_slice[index_isrc + ix].imag;
            }
        }

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            index_isrc = index_iz + index_isrc_s;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                /*
                 *if iz == nz
                 *   dQmin(nz) = R .* dPplus(nz);
                 *else
                 *   dQmin(iz) = (1 - R) .* dPmin(iz) + R .* dPplus(iz);
                 *end
                 *pfcomp_tmp2 is dQmin; pfcomp_tmp11 is dPplus; pfcomp_dPmin0_slice_tmp is dPmin
                 */
                if (iz == p_gbparams->stc_nz)
                {
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].real
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp11[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].imag
                        = pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp11[index_isrc + ix].imag;
                }
                else
                {
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].real
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].real
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp11[index_isrc + ix].real;
                    p_tmpvars->pfcomp_tmp2[index_isrc_s + ix].imag
                        = (1.0 - pf_image[index_iz_m + ix]) * pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].imag
                        + pf_image[index_iz_m + ix] * p_tmpvars->pfcomp_tmp11[index_isrc + ix].imag;

                }
            }
        }

        /*
         *calculate dPmin(iz-1) = Wx * dQmin(iz) + Gx * diag(dslow) * Qmin(iz) = dPmin1 + dPmin2
         *calculate dPmin1 = Wx * dQmin(iz)
         *pfcomp_tmp2 is dQmin; pfcomp_dPmin0_slice_tmp is dPmin1
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_tmp2, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Wx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, pfcomp_dPmin0_slice_tmp, p_gbparams->stc_nx2);

        /*
         *calculate dPmin2 = diag(dslow) * Gx * Qmin(iz)
         *pfcomp_tmp33 is dPmin2; pfcomp_Qmin_slice is Qmin
         */
        cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans
                , p_gbparams->stc_nsrc, p_gbparams->stc_nx2, p_gbparams->stc_nx2
                , &fcomp_alpha, p_tmpvars->pfcomp_Qmin_slice_tmp, p_gbparams->stc_nx2
                , p_tmpvars->pfcomp_Gx_tmp, p_gbparams->stc_nx2
                , &fcomp_beta, p_tmpvars->pfcomp_tmp33, p_gbparams->stc_nx2);
        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc_s = isrc * p_gbparams->stc_nx2;
            for (ix=0; ix<p_gbparams->stc_nx2; ix++)
            {
                pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].real += (pf_dslow[index_izminus1_m + ix] * p_tmpvars->pfcomp_tmp33[index_isrc_s + ix].real);
                pfcomp_dPmin0_slice_tmp[index_isrc_s + ix].imag += (pf_dslow[index_izminus1_m + ix] * p_tmpvars->pfcomp_tmp33[index_isrc_s + ix].imag);
            }
        }
    }
}

void Get_dPmin0_velocity(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , float *pf_numerator_full, float *pf_denominator_full
        , const unsigned long MPI_ID, const unsigned long MPI_SIZE)
{
    unsigned long ifreq;
    unsigned long isrcx;

    unsigned long index_slice;
    unsigned long index_slice0;
    unsigned long index_slice_o;
    float f_numerator = 0.0;
    float f_denominator = 0.0;
    *pf_numerator_full = 0.0;
    *pf_denominator_full = 0.0;

    fcomp *pfcomp_dPmin0_slice_tmp = NULL;
    pfcomp_dPmin0_slice_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        memset(pfcomp_dPmin0_slice_tmp, 0, sizeof(float) * 2 * p_gbparams->stc_nsrcnx2);
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;

        Get_dPmin0_velocity_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice, p_gbvars->pfcomp_Qplus_Nslice + index_slice
                , pfcomp_dPmin0_slice_tmp
                , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pfcomp_Gx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                , p_gbvars->pst_i_vel, p_gbvars->pf_dslow, p_gbvars->pf_image);

        for (isrcx=0; isrcx<p_gbparams->stc_nsrcnx2; isrcx++)
        {
            //numerator and denominator for alpha
            f_numerator += (pfcomp_dPmin0_slice_tmp[isrcx].real * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].real
                    + pfcomp_dPmin0_slice_tmp[isrcx].imag * p_gbvars->pfcomp_Res_Nslice[index_slice0 + isrcx].imag);
            f_denominator += (pfcomp_dPmin0_slice_tmp[isrcx].real * pfcomp_dPmin0_slice_tmp[isrcx].real
                    + pfcomp_dPmin0_slice_tmp[isrcx].imag * pfcomp_dPmin0_slice_tmp[isrcx].imag);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&f_numerator, pf_numerator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&f_denominator, pf_denominator_full, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(pfcomp_dPmin0_slice_tmp);
}

void Combine_slowness_gradient_output_pershot(const globalconsts *p_gbparams, float *pf_dslow_allshot, float *pf_illummatrix_velocity, const float *pf_dslow_Nslice_allshot, const float *pf_illummatrix_velocity_Nslice, const unsigned long Nslice_tmp)
{
    unsigned long islice;
    memset(pf_dslow_allshot, 0, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2);
    memset(pf_illummatrix_velocity, 0, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    for (islice=0; islice<Nslice_tmp; islice++)
    {
        cblas_saxpy(p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2, 1.0
                , pf_dslow_Nslice_allshot + islice * p_gbparams->stc_nz * p_gbparams->stc_nsrcnx2, 1, pf_dslow_allshot, 1);
        cblas_saxpy(p_gbparams->stc_nz * p_gbparams->stc_nx2, 1.0
                , pf_illummatrix_velocity_Nslice + islice * p_gbparams->stc_nz * p_gbparams->stc_nx2, 1, pf_illummatrix_velocity, 1);
    }
}

void Combine_slowness_gradient(const globalconsts *p_gbparams, float *pf_dslow, float *pf_illummatrix_velocity, const float *pf_dslow_Nslice, const float *pf_illummatrix_velocity_Nslice, const unsigned long Nslice_tmp)
{
    unsigned long islice;
    memset(pf_dslow, 0, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    memset(pf_illummatrix_velocity, 0, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nx2);
    for (islice=0; islice<Nslice_tmp; islice++)
    {
        cblas_saxpy(p_gbparams->stc_nz * p_gbparams->stc_nx2, 1.0
                , pf_dslow_Nslice + islice * p_gbparams->stc_nz * p_gbparams->stc_nx2, 1, pf_dslow, 1);
        cblas_saxpy(p_gbparams->stc_nz * p_gbparams->stc_nx2, 1.0
                , pf_illummatrix_velocity_Nslice + islice * p_gbparams->stc_nz * p_gbparams->stc_nx2, 1, pf_illummatrix_velocity, 1);
    }
}

void Get_index_velocity_tmp(const globalconsts *p_gbparams, unsigned long *pst_i_vel, float *pf_vel, const float *pf_slow, const float *pf_dslow, const float *pf_model_mask, const float f_alpha)
{
    //pf_slow is intact in this function
    unsigned long izx;
    unsigned long ierror;
    for (izx=0; izx<p_gbparams->stc_nz * p_gbparams->stc_nx2; izx++)
    {
        pf_vel[izx] = 1.0 / (pf_slow[izx] + f_alpha * pf_dslow[izx]);
    }
    // Copy the velocity model in the taper range if the mask indicates this; only for option=2 or 3
    if (p_gbparams->stc_if_model_mask >= 2)
    {
        Apply_model_mask_extend(p_gbparams, pf_vel, p_gbparams->stc_nz, pf_model_mask);
    }

    Get_current_operatorindex(p_gbparams, pst_i_vel, pf_vel);
}

void Write_results_model_gradient_2su(const globalconsts *p_gbparams, const float *pf_dslow, const float *pf_dimage, const unsigned long iter)
{
    FILE *fp;
    char filename[64];
    unsigned long ix;

    segy seg_tmp;
    memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);

    snprintf(filename, sizeof(char) * 64, "%s%s_dimage.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    seg_tmp.trwf = (int)p_gbparams->stc_nx;
    seg_tmp.ns = (int)p_gbparams->stc_nzplus1;
    seg_tmp.d2 = p_gbparams->fc_dx;
    seg_tmp.timbas = 11; //sort on keyword duse
    seg_tmp.d1 = p_gbparams->fc_dz;
    seg_tmp.duse = (int)(iter + 1);
    seg_tmp.trid = 30;
    fp= fopen(filename, "a+");
    if (fp == NULL)
    {
        syserr("error: cannot open file %s to write P file.\n", filename);
    }

    for (ix=0; ix<p_gbparams->stc_nx; ix++)
    {
        cblas_scopy(p_gbparams->stc_nzplus1, pf_dimage + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                , seg_tmp.data, 1);
        seg_tmp.tracl = (int)(iter * p_gbparams->stc_nx + ix + 1);
        seg_tmp.tracf = (int)(ix + 1);
        fwrite(&seg_tmp, sizeof(float)
                , HDRSIZE / sizeof(float) + p_gbparams->stc_nzplus1, fp);
    }
    fclose(fp);

    if (p_gbparams->stc_velupdate_iterjump != 0)
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_dslowness.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        seg_tmp.ns = (int)p_gbparams->stc_nz;
        fp= fopen(filename, "a+");
        if (fp == NULL)
        {
            syserr("error: cannot open file %s to write P file.\n", filename);
        }

        for (ix=0; ix<p_gbparams->stc_nx; ix++)
        {
            cblas_scopy(p_gbparams->stc_nz, pf_dslow + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                    , seg_tmp.data, 1);
            seg_tmp.tracl = (int)(iter * p_gbparams->stc_nx + ix + 1);
            seg_tmp.tracf = (int)(ix + 1);
            fwrite(&seg_tmp, sizeof(float)
                    , HDRSIZE / sizeof(float) + p_gbparams->stc_nz, fp);
        }
        fclose(fp);
    }
}

void Write_results_model_2su(const globalconsts *p_gbparams, const float *pf_vel, const float *pf_image, const float *pf_J, const float *pf_J_image_allshot, const unsigned long iter)
{
    FILE *fp;
    char filename[64];
    unsigned long ix;

    segy seg_tmp;
    memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);
    seg_tmp.trwf = (int)p_gbparams->stc_itertotal;
    seg_tmp.ns = p_gbparams->stc_itertotal;
    seg_tmp.d2 = 1;
    seg_tmp.timbas = 11;//sort on header duse
    seg_tmp.d1 = 1;
    seg_tmp.tracf = 1;
    seg_tmp.tracl = 1;
    seg_tmp.duse = (int)(iter + 1);
    snprintf(filename, sizeof(char) * 64, "%s%s_misfit.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    fp= fopen(filename, "wb");
    if (fp == NULL)
    {
        printf("error: cannot open file %s to write P file.\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    cblas_scopy(p_gbparams->stc_itertotal, pf_J, 1
           , seg_tmp.data, 1);
    fwrite(&seg_tmp, sizeof(float), HDRSIZE / sizeof(float) + p_gbparams->stc_itertotal, fp);
    fclose(fp);

    if (p_gbparams->stc_outpershot != 0)
    {
        seg_tmp.ns = (int)p_gbparams->stc_nsrc;
        snprintf(filename, sizeof(char) * 64, "%s%s_misfit_pershot.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        fp= fopen(filename, "a+");
        if (fp == NULL)
        {
            printf("error: cannot open file %s to write P file.\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cblas_scopy(p_gbparams->stc_nsrc, pf_J_image_allshot, 1
                , seg_tmp.data, 1);
        fwrite(&seg_tmp, sizeof(float), HDRSIZE / sizeof(float) + p_gbparams->stc_nsrc, fp);
        fclose(fp);
    }

    snprintf(filename, sizeof(char) * 64, "%s%s_inverted_image.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    seg_tmp.trwf = (int)p_gbparams->stc_nx;
    seg_tmp.ns = (int)p_gbparams->stc_nzplus1;
    seg_tmp.d2 = p_gbparams->fc_dx;
    seg_tmp.timbas = 11;//sorting on header duse
    seg_tmp.d1 = p_gbparams->fc_dz;
    seg_tmp.duse = (int)(iter + 1);
    fp= fopen(filename, "a+");
    if (fp == NULL)
    {
        syserr("error: cannot open file %s to write P file.\n", filename);
    }

    for (ix=0; ix<p_gbparams->stc_nx; ix++)
    {
        cblas_scopy(p_gbparams->stc_nzplus1, pf_image + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                , seg_tmp.data, 1);
        seg_tmp.tracl = (int)(iter * p_gbparams->stc_nx + ix + 1);
        seg_tmp.tracf = (int)(ix + 1);
        fwrite(&seg_tmp, sizeof(float)
                , HDRSIZE / sizeof(float) + p_gbparams->stc_nzplus1, fp);
    }
    fclose(fp);

    if (p_gbparams->stc_velupdate_iterjump != 0)
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_inverted_velocity.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        seg_tmp.ns = (int)p_gbparams->stc_nz;
        fp= fopen(filename, "a+");
        if (fp == NULL)
        {
            syserr("error: cannot open file %s to write P file.\n", filename);
        }

        for (ix=0; ix<p_gbparams->stc_nx; ix++)
        {
            cblas_scopy(p_gbparams->stc_nz, pf_vel + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                    , seg_tmp.data, 1);
            seg_tmp.tracl = (int)(iter * p_gbparams->stc_nx + ix + 1);
            seg_tmp.tracf = (int)(ix + 1);
            fwrite(&seg_tmp, sizeof(float)
                    , HDRSIZE / sizeof(float) + p_gbparams->stc_nz, fp);
        }
        fclose(fp);
    }
}

void Write_final_results_2su(const globalconsts *p_gbparams, const float *pf_vel, const float *pf_image)
{
    FILE *fp;
    char filename[64];
    unsigned long ix;

    segy seg_tmp;
    memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);

    if (p_gbparams->stc_velupdate_iterjump != 0)
    {
        snprintf(filename, sizeof(char) * 64, "%s%s_final_inverted_velocity.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        seg_tmp.trwf = (unsigned long)(p_gbparams->stc_nx);
        seg_tmp.ns = (unsigned long)p_gbparams->stc_nz;
        seg_tmp.d2 = p_gbparams->fc_dx;
        seg_tmp.timbas = 3;//sort on header fldr
        seg_tmp.d1 = p_gbparams->fc_dz;
        fp= fopen(filename, "wb");
        if (fp == NULL)
        {
            syserr("error: cannot open file final_velocity to write P file.\n");
        }

        for (ix=0; ix<p_gbparams->stc_nx; ix++)
        {
            cblas_scopy(p_gbparams->stc_nz, pf_vel + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                    , seg_tmp.data, 1);
            seg_tmp.tracl = (unsigned long)(ix + 1);
            seg_tmp.tracf = (unsigned long)(ix + 1);
            fwrite(&seg_tmp, sizeof(float)
                    , HDRSIZE / sizeof(float) + p_gbparams->stc_nz, fp);
        }
        fclose(fp);
    }

    snprintf(filename, sizeof(char) * 64, "%s%s_final_inverted_image.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    seg_tmp.trwf = (unsigned long)(p_gbparams->stc_nx);
    seg_tmp.ns = (unsigned long)p_gbparams->stc_nzplus1;
    seg_tmp.d2 = p_gbparams->fc_dx;
    seg_tmp.timbas = 3;//sort on header fldr
    seg_tmp.d1 = p_gbparams->fc_dz;
    fp= fopen(filename, "wb");
    if (fp == NULL)
    {
        syserr("error: cannot open file final_velocity to write P file.\n");
    }

    for (ix=0; ix<p_gbparams->stc_nx; ix++)
    {
        cblas_scopy(p_gbparams->stc_nzplus1, pf_image + p_gbparams->stc_nxtap + ix, p_gbparams->stc_nx2
                , seg_tmp.data, 1);
        seg_tmp.tracl = (unsigned long)(ix + 1);
        seg_tmp.tracf = (unsigned long)(ix + 1);
        fwrite(&seg_tmp, sizeof(float)
                , HDRSIZE / sizeof(float) + p_gbparams->stc_nzplus1, fp);
    }
    fclose(fp);
}

void Free_memory(const globalconsts *p_gbparams, globalarrays *p_gbvars
        , const unsigned long MPI_ID)
{
    free((p_gbvars)->pf_inivel);
    if (p_gbparams->stc_if_model_mask >= 1)
    {
        free((p_gbvars)->pf_model_mask);
    }
    free((p_gbvars)->pfcomp_src);
    free((p_gbvars)->pfcomp_data);
    free((p_gbvars)->pf_WxGxtaper);

    if (p_gbparams->stc_data_masktype == 1)
    {
        free((p_gbvars)->pf_data_mask);
    }
    else if (p_gbparams->stc_data_masktype == 2)
    {
        if (MPI_ID == 0)
        {
            free((p_gbvars)->pf_data_mask);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (p_gbparams->fc_DTV_lam >= EPS)
    {
        if (MPI_ID == 0)
        {
            free((p_gbvars)->pf_dip_field);
        }
    }


    free((p_gbvars)->pfcomp_Wx_Nslice);
    free((p_gbvars)->pfcomp_Gx_Nslice);
    free((p_gbvars)->pfcomp_Pmin_Nslice);
    free((p_gbvars)->pfcomp_Pplus_Nslice);
    free((p_gbvars)->pfcomp_Qplus_Nslice);
    free((p_gbvars)->pfcomp_Res_Nslice);


    free((p_gbvars)->pf_image);
    free((p_gbvars)->pf_image_tmp);

    free((p_gbvars)->pst_i_vel);
    free((p_gbvars)->pf_dimage);
    free((p_gbvars)->pf_dslow);
    if (MPI_ID == 0)
    {
        free((p_gbvars)->pf_vel);
        free((p_gbvars)->pf_slow);
        free((p_gbvars)->pf_data_tmp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Free_memory_fortmpvars(const globalconsts *p_gbparams, tmparrays *p_tmpvars)
{
    free((p_tmpvars)->pfcomp_Wx_tmp);
    free((p_tmpvars)->pfcomp_Gx_tmp);
    free((p_tmpvars)->pfcomp_Qmin_slice_tmp);
    free((p_tmpvars)->pfcomp_dPmin0_slice_tmp_tmp);
    free((p_tmpvars)->pfcomp_tmp1);
    free((p_tmpvars)->pfcomp_tmp2);
    free((p_tmpvars)->pfcomp_tmp3);
    free((p_tmpvars)->pfcomp_tmp3_tmp);
    free((p_tmpvars)->pfcomp_tmp4);
    free((p_tmpvars)->pfcomp_tmp);
    free((p_tmpvars)->pfcomp_tmp11);
    free((p_tmpvars)->pfcomp_tmp33);
}

void Free_memory_fortmpvars_mod(const globalconsts *p_gbparams, tmparrays *p_tmpvars)
{
    free((p_tmpvars)->pfcomp_Wx_tmp);
    free((p_tmpvars)->pfcomp_Qmin_slice_tmp);
}
fcomp *Allocate_1d_floatcomplex(unsigned long array_size)
{
    fcomp *p = (fcomp *)malloc(array_size * sizeof(float) * 2);

    if (p == NULL)
    {
        printf("%s: malloc failed!\n", __FILE__);
        exit(-1);
    }

    //set it to a zero array
    memset(p, 0, sizeof(float) * 2 * array_size);

    return p;
}

float *Allocate_1d_float(unsigned long array_size)
{
    float *p = (float *)malloc(array_size * sizeof(float));

    if (p == NULL)
    {
        printf("%s: malloc failed!\n", __FILE__);
        exit(-1);
    }

    //set it to a zero array
    memset(p, 0, sizeof(float) * array_size);

    return p;
}

int *Allocate_1d_int(unsigned long array_size)
{
    int *p = (int *)malloc(array_size * sizeof(int));

    if (p == NULL)
    {
        printf("%s: malloc failed!\n", __FILE__);
        exit(-1);
    }

    //set it to a zero array
    memset(p, 0, sizeof(int) * array_size);

    return p;
}

unsigned long *Allocate_1d_unsignedlong(unsigned long array_size)
{
    unsigned long *p = (unsigned long *)malloc(array_size * sizeof(unsigned long));

    if (p == NULL)
    {
        printf("%s: malloc failed!\n", __FILE__);
        exit(-1);
    }

    //set it to a zero array
    memset(p, 0, sizeof(unsigned long) * array_size);

    return p;
}

unsigned long Convert2nextpow2(unsigned long in)
{
    if ((in % 2) == 0)
    {
        return in;
    }
    else
    {
        return (in+1);
    }
    //unsigned long in_tmp = 1;
    //while (1)
    //{
        //in_tmp *= 2;
        //if (in_tmp >= in)
        //{
            //return in_tmp;
        //}
    //}
}

float mean(unsigned long N, float *a)
{
    float sum = 0.0;
    unsigned long i;

    for (i=0; i<N; i++)
    {
        sum += a[i];
    }

    return(sum / (float)N);
}

void conv_jmi(const float *a, const float *b, float *c, const unsigned long size_a, const unsigned long size_b)
{
    unsigned long i_start;
    unsigned long i;
    unsigned long j;
    unsigned long jmin;
    unsigned long jmax;

    i_start = (unsigned long)floor((float)(size_b + 1) / 2.0) - 1;
    for (i=i_start; i<(i_start+size_a); i++)
    {
        c[i-i_start] = 0.0;
        if (i >= (size_b - 1))
        {
            jmin = i - size_b + 1;
        }
        else
        {
            jmin = 0;
        }
        if (i < (size_a - 1))
        {
            jmax = i;
        }
        else
        {
            jmax = size_a - 1;
        }
        for (j=jmin; j<=jmax; j++)
        {
            c[i - i_start] += (a[j] * b[i-j]);
        }
    }
}

void DTV_DDfilter(float *pf_tvx, float *pf_tvz, const float *pf_slow, const float *pf_dip_field, const float alpha, const unsigned long nx2, const unsigned long nz)
{
    float f_h = 0.5;
    unsigned long ix;
    unsigned long iz;
    float ux;
    float uz;

    for (iz=0; iz<nz-1; iz++)
    {
        for (ix=0; ix<nx2-1; ix++)
        {
            ux = pf_slow[iz * nx2 + ix] * f_h
                + pf_slow[(iz + 1) * nx2 + ix] * f_h
                + pf_slow[iz * nx2 + ix + 1] * (-f_h)
                + pf_slow[(iz + 1) * nx2 + ix + 1] * (-f_h);
            uz = pf_slow[iz * nx2 + ix] * f_h
                + pf_slow[(iz + 1) * nx2 + ix] * (-f_h)
                + pf_slow[iz * nx2 + ix + 1] * f_h
                + pf_slow[(iz + 1) * nx2 + ix + 1] * (-f_h);
            pf_tvx[iz * nx2 + ix] = alpha
                * (cos(pf_dip_field[iz * nx2 + ix]) * ux + sin(pf_dip_field[iz * nx2 + ix]) * uz);
            pf_tvz[iz * nx2 + ix] = (2.0 - alpha)
                * (cos(pf_dip_field[iz * nx2 + ix]) * uz - sin(pf_dip_field[iz * nx2 + ix]) * ux);
        }
    }
    for (ix=0; ix<nx2-1; ix++)
    {
        pf_tvx[(nz - 1) * nx2 + ix] = pf_tvx[(nz - 2) * nx2 + ix];
        pf_tvz[(nz - 1) * nx2 + ix] = pf_tvz[(nz - 2) * nx2 + ix];
    }
    for (iz=0; iz<nz-1; iz++)
    {
        pf_tvx[iz * nx2 + nx2 - 1] = pf_tvx[iz * nx2 + nx2 - 2];
        pf_tvz[iz * nx2 + nx2 - 1] = pf_tvz[iz * nx2 + nx2 - 2];
    }
    pf_tvx[(nz - 1) * nx2 + nx2 - 1] = pf_tvx[(nz - 2) * nx2 + nx2 - 2];
    pf_tvz[(nz - 1) * nx2 + nx2 - 1] = pf_tvz[(nz - 2) * nx2 + nx2 - 2];
}

void DTV_DDTfilter(float *pf_dtv, float *pf_tvx, float *pf_tvz, const float *pf_slow, const float *pf_dip_field, const float alpha, const unsigned long nx2, const unsigned long nz)
{
    float f_h = 0.5;
    unsigned long ix;
    unsigned long iz;
    float ux;
    float uz;

    float *pf_ux = NULL;
    float *pf_uz = NULL;
    pf_ux = Allocate_1d_float(nz * nx2);
    pf_uz = Allocate_1d_float(nz * nx2);

    for (ix=0; ix<nz*nx2; ix++)
    {
        pf_ux[ix] = alpha * pf_tvx[ix];
        pf_uz[ix] = (2.0 - alpha) * pf_tvz[ix];
        pf_tvx[ix] = cos(pf_dip_field[ix]) * pf_ux[ix] - sin(pf_dip_field[ix]) * pf_uz[ix];
        pf_tvz[ix] = cos(pf_dip_field[ix]) * pf_uz[ix] + sin(pf_dip_field[ix]) * pf_ux[ix];
    }

    for (iz=1; iz<nz; iz++)
    {
        for (ix=1; ix<nx2; ix++)
        {
            pf_ux[iz * nx2 + ix] = pf_tvx[(iz - 1) * nx2 + ix - 1] * (-f_h)
                + pf_tvx[iz * nx2 + ix - 1] * (-f_h)
                + pf_tvx[(iz - 1) * nx2 + ix] * f_h
                + pf_tvx[iz * nx2 + ix] * f_h;
            pf_uz[iz * nx2 + ix] = pf_tvz[(iz - 1) * nx2 + ix - 1] * (-f_h)
                + pf_tvz[iz * nx2 + ix - 1] * f_h
                + pf_tvz[(iz - 1) * nx2 + ix] * (-f_h)
                + pf_tvz[iz * nx2 + ix] * f_h;
        }
    }
    for (ix=1; ix<nx2; ix++)
    {
        pf_ux[ix] = pf_ux[nx2 + ix];
        pf_uz[ix] = pf_uz[nx2 + ix];
    }
    for (iz=1; iz<nz; iz++)
    {
        pf_ux[iz * nx2] = pf_ux[iz * nx2 + 1];
        pf_uz[iz * nx2] = pf_uz[iz * nx2 + 1];
    }
    pf_ux[0] = pf_ux[nx2 + 1];
    pf_uz[0] = pf_uz[nx2 + 1];

    for (ix=0; ix<nz*nx2; ix++)
    {
        pf_dtv[ix] = pf_ux[ix] + pf_uz[ix];
    }
    free(pf_ux);
    free(pf_uz);
}

void Get_globalparameters_mod(int argc, char *argv[], globalconsts *p_gbparams)
{
    //this function:
    //get global parameters from input
    //calculate some other global parameters
    initargs(argc, argv);
    requestdoc(1);
    char filename[64];
    FILE *fp;

    if (!getparstring("outfolder", &(p_gbparams->pc_outfolder)))
    {
        p_gbparams->pc_outfolder = "Results/";
    }

    if (!getparstring("outfile_label", &(p_gbparams->pc_outlabel)))
    {
        p_gbparams->pc_outlabel = "fwmod";
    }

    if (!getparstring("velocity", &(p_gbparams->pc_vel)))
    {
        p_gbparams->pc_vel = "vel0.su";
    }

    if (!getparstring("density", &(p_gbparams->pc_den)))
    {
        p_gbparams->pc_den = "den0.su";
    }

    if (!getparstring("source", &(p_gbparams->pc_src)))
    {
        printf("error! no source .su file is specified!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (!getparulong("multiple_order", &(p_gbparams->stc_multiple_order)))
    {
        p_gbparams->stc_multiple_order = 3;
    }

    if (!getparfloat("fmin", &(p_gbparams->fc_fmin)))
    {
        p_gbparams->fc_fmin = 4.0;
    }

    if (!getparfloat("fmax_upper", &(p_gbparams->fc_fmax_upper)))
    {
        p_gbparams->fc_fmax_upper = 40.0;
    }

    if (!getparfloat("operator_vmin", &(p_gbparams->fc_oper_vmin)))
    {
        p_gbparams->fc_oper_vmin = 500.0;
    }

    if (!getparfloat("operator_vmax", &(p_gbparams->fc_oper_vmax)))
    {
        p_gbparams->fc_oper_vmax = 8000.0;
    }

    if (!getparfloat("operator_dv", &(p_gbparams->fc_oper_dv)))
    {
        p_gbparams->fc_oper_dv = 1.0;
    }

    if (!getparfloat("dx", &(p_gbparams->fc_dx)))
    {
        p_gbparams->fc_dx = 25.0;
    }

    if (!getparfloat("dz", &(p_gbparams->fc_dz)))
    {
        p_gbparams->fc_dz = 10.0;
    }

    if (!getparfloat("dt", &(p_gbparams->fc_dt)))
    {
        p_gbparams->fc_dt = 0.004;
    }

    if (!getparulong("size_x", &(p_gbparams->stc_nx)))
    {
        p_gbparams->stc_nx = 0;
    }

    if (!getparulong("size_xtap", &(p_gbparams->stc_nxtap)))
    {
        p_gbparams->stc_nxtap = 10;
    }

    if (!getparulong("size_z", &(p_gbparams->stc_nz)))
    {
        p_gbparams->stc_nz = 0;
    }

    if (!getparulong("size_t", &(p_gbparams->stc_nt)))
    {
        p_gbparams->stc_nt = 0;
    }

    if (!getparulong("size_src", &(p_gbparams->stc_nsrc)))
    {
        p_gbparams->stc_nsrc = 0;
    }

    if (!getparfloat("angle", &(p_gbparams->fc_angle)))
    {
        p_gbparams->fc_angle = 80.0;
    }


    //calculate other parameters
    p_gbparams->stc_nx2 = p_gbparams->stc_nx + 2 * p_gbparams->stc_nxtap;
    p_gbparams->stc_nx2pow2 = Convert2nextpow2(p_gbparams->stc_nx2);
    p_gbparams->stc_nzplus1 = p_gbparams->stc_nz + 1;
    p_gbparams->stc_nf = p_gbparams->stc_nt / 2 + 1;
    p_gbparams->fc_df = 1.0 / ((float) p_gbparams->stc_nt * p_gbparams->fc_dt);
    p_gbparams->stc_i_fmin = round((p_gbparams->fc_fmin / p_gbparams->fc_df));
    p_gbparams->stc_i_fmax_upper = round((p_gbparams->fc_fmax_upper / p_gbparams->fc_df));
    p_gbparams->stc_nsrcnx = p_gbparams->stc_nsrc * p_gbparams->stc_nx;
    p_gbparams->stc_nsrcnx2 = p_gbparams->stc_nsrc * p_gbparams->stc_nx2;
    p_gbparams->stc_i_oper_vmin = round((p_gbparams->fc_oper_vmin
            / p_gbparams->fc_oper_dv) + 1);
    p_gbparams->stc_i_oper_vmax = round((p_gbparams->fc_oper_vmax
            / p_gbparams->fc_oper_dv) + 1);
    p_gbparams->stc_valid_nf = p_gbparams->stc_i_fmax_upper
        - p_gbparams->stc_i_fmin + 1;
    p_gbparams->stc_valid_nv = p_gbparams->stc_i_oper_vmax
        - p_gbparams->stc_i_oper_vmin + 1;
}

void Print_globalparameters_mod(const globalconsts *p_gbparams)
{
    char filename[64];
    snprintf(filename, sizeof(char) * 64, "%s%s_data0.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
    printf("... pc_outfile_data = %s\n", filename);

    printf("... pc_velocity = %s\n", p_gbparams->pc_vel);
    printf("... pc_density = %s\n", p_gbparams->pc_den);
    printf("... pc_source = %s\n", p_gbparams->pc_src);

    printf("... fmin = %f\n", p_gbparams->fc_fmin);
    printf("... fmax_upper = %f\n", p_gbparams->fc_fmax_upper);
    printf("... operator_vmin = %f\n", p_gbparams->fc_oper_vmin);
    printf("... operator_vmax = %f\n", p_gbparams->fc_oper_vmax);
    printf("... operator_dv = %f\n", p_gbparams->fc_oper_dv);
    printf("... dt = %f\n", p_gbparams->fc_dt);
    printf("... dx = %f\n", p_gbparams->fc_dx);
    printf("... dz = %f\n", p_gbparams->fc_dz);
    printf("... size_x = %zu\n", p_gbparams->stc_nx);
    printf("... size_xtap = %zu\n", p_gbparams->stc_nxtap);
    printf("... size_z = %zu\n", p_gbparams->stc_nz);
    printf("... size_t = %zu\n", p_gbparams->stc_nt);
    printf("... size_src = %zu\n", p_gbparams->stc_nsrc);
    printf("... angle = %f\n", p_gbparams->fc_angle);


    printf("... size_x2 = %zu\n", p_gbparams->stc_nx2);
    printf("... size_x2pow2 = %zu\n", p_gbparams->stc_nx2pow2);
    printf("... size_zplus1 = %zu\n", p_gbparams->stc_nzplus1);
    printf("... size_f = %zu\n", p_gbparams->stc_nf);
    printf("... df = %f\n", p_gbparams->fc_df);
    printf("... index_fmin = %zu\n", p_gbparams->stc_i_fmin);
    printf("... index_fmax_upper = %zu\n", p_gbparams->stc_i_fmax_upper);
    printf("... size_srcx = %zu\n", p_gbparams->stc_nsrcnx);
    printf("... size_srcx2 = %zu\n", p_gbparams->stc_nsrcnx2);
    printf("... index_operator_vmin = %zu\n", p_gbparams->stc_i_oper_vmin);
    printf("... index_operator_vmax = %zu\n", p_gbparams->stc_i_oper_vmax);
    printf("... valid_nf = %zu\n", p_gbparams->stc_valid_nf);
    printf("... valid_nv = %zu\n", p_gbparams->stc_valid_nv);

    fflush(stdout);
}

void Create_memory_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars, const unsigned long MPI_ID, const unsigned long Nslice)
{
    (p_gbvars)->pfcomp_src = Allocate_1d_floatcomplex(p_gbparams->stc_valid_nf * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pf_WxGxtaper = Allocate_1d_float(p_gbparams->stc_nx2 * p_gbparams->stc_nx2);

    (p_gbvars)->pfcomp_Wx_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2);
    (p_gbvars)->pfcomp_Gx_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2);
    (p_gbvars)->pfcomp_Pmin_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_Pplus_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);
    (p_gbvars)->pfcomp_Qplus_Nslice = Allocate_1d_floatcomplex(Nslice * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2);

    (p_gbvars)->pf_image = Allocate_1d_float(p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2);
    (p_gbvars)->pst_i_vel = Allocate_1d_unsignedlong(p_gbparams->stc_nz * p_gbparams->stc_nx2);

    if (MPI_ID == 0)
    {
        (p_gbvars)->pf_vel = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        (p_gbvars)->pf_den = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
        (p_gbvars)->pf_slow = Allocate_1d_float(p_gbparams->stc_nz * p_gbparams->stc_nx2);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Read_input_extended_model_data_fromsu_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars)
{
    //this function:
    //read velocity/density from .su, and save into size_z * size_x2
    //read source from .su, and save into size_f * size_src * size_x2
    segy seg_tmp;
    FILE *fp = NULL;
    unsigned long iz;
    unsigned long ix;
    unsigned long ifreq;
    unsigned long it;
    unsigned long fldr_old;
    unsigned long source_index_tmp;

    //read velocity
    fp = fopen(p_gbparams->pc_vel, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_vel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_vel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                && (p_gbparams->stc_nz == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dz == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_vel);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pf_vel[], pf_vel[] size: size_z * size_x2
    cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
            , p_gbvars->pf_vel + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
    while(fvgettr(fp, &seg_tmp))
    {
        cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                , p_gbvars->pf_vel + p_gbparams->stc_nxtap + seg_tmp.tracf - 1, p_gbparams->stc_nx2);
    }

    fclose(fp);

    //extend the velocity model with xtap
    unsigned long index_iz = 0;
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
        {
            p_gbvars->pf_vel[index_iz + ix] = p_gbvars->pf_vel[index_iz + p_gbparams->stc_nxtap];
        }
        for (ix=p_gbparams->stc_nx+p_gbparams->stc_nxtap
                ; ix<p_gbparams->stc_nx2; ix++)
        {
            p_gbvars->pf_vel[index_iz + ix]
                = p_gbvars->pf_vel[index_iz + p_gbparams->stc_nx + p_gbparams->stc_nxtap - 1];
        }
    }

    printf("input velocity from file %s is read in successfully, nz = %zu, nx = %zu, dz = %f, dx = %f\n"
            ,p_gbparams->pc_vel, p_gbparams->stc_nz, p_gbparams->stc_nx
            , p_gbparams->fc_dz, p_gbparams->fc_dx);
    fflush(stdout);

    //read density
    fp = fopen(p_gbparams->pc_den, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_den);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_den);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                &&  (p_gbparams->stc_nz == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dz == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_den);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pf_den[], pf_den[] size: size_z * size_x
    cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
            , p_gbvars->pf_den + p_gbparams->stc_nxtap + seg_tmp.tracl - 1, p_gbparams->stc_nx2);
    while(fvgettr(fp, &seg_tmp))
    {
        cblas_scopy(p_gbparams->stc_nz, seg_tmp.data, 1
                , p_gbvars->pf_den + p_gbparams->stc_nxtap + seg_tmp.tracl - 1, p_gbparams->stc_nx2);
    }
    fclose(fp);

    //extend the density model with xtap
    for (iz=0; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nxtap; ix++)
        {
            p_gbvars->pf_den[index_iz + ix] = p_gbvars->pf_den[index_iz + p_gbparams->stc_nxtap];
        }
        for (ix=p_gbparams->stc_nx+p_gbparams->stc_nxtap
                ; ix<p_gbparams->stc_nx2; ix++)
        {
            p_gbvars->pf_den[index_iz + ix]
                = p_gbvars->pf_den[index_iz + p_gbparams->stc_nx + p_gbparams->stc_nxtap - 1];
        }
    }

    printf("input density from file %s is read in successfully, nz = %zu, nx = %zu, dz = %f, dx = %f\n"
            , p_gbparams->pc_den, p_gbparams->stc_nz, p_gbparams->stc_nx
            , p_gbparams->fc_dz, p_gbparams->fc_dx);
    fflush(stdout);

    //calculate image, pf_image[] size: (size_z+1) * size_x
    float tmp1 = 0.0;
    float tmp2 = 0.0;
    for (iz=1; iz<p_gbparams->stc_nz; iz++)
    {
        index_iz = iz * p_gbparams->stc_nx2;
        for (ix=0; ix<p_gbparams->stc_nx2; ix++)
        {
            tmp1 = (p_gbvars->pf_den[index_iz - p_gbparams->stc_nx2 + ix])
                * (p_gbvars->pf_vel[index_iz - p_gbparams->stc_nx2 + ix]);
            tmp2 = (p_gbvars->pf_den[index_iz + ix])
                * (p_gbvars->pf_vel[index_iz + ix]);
            p_gbvars->pf_image[index_iz + ix] = (tmp2 - tmp1) / (tmp2 + tmp1);
        }
    }

    #ifdef DEBUG
    fp = fopen("Tmp/vel_test.bin", "wb");
    fwrite(p_gbvars->pf_vel, sizeof(float) * p_gbparams->stc_nz * p_gbparams->stc_nx2
            , 1, fp);
    fclose(fp);
    fp = fopen("Tmp/image_test.bin", "wb");
    fwrite(p_gbvars->pf_image, sizeof(float) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nx2
            , 1, fp);
    fclose(fp);
    #endif

    //read source wavefield and convert it into frequency domain
    fp = fopen(p_gbparams->pc_src, "rb");
    if (fp == NULL)
    {
        printf("fail to open file %s\n", p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //read in first trace
    if (!fvgettr(fp, &seg_tmp))
    {
        printf("can't get the first trace in file %s\n!", p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    fldr_old = seg_tmp.fldr;

    //check if the input parameters match the ones in .su header
    if (!((p_gbparams->stc_nx == (unsigned long)seg_tmp.trwf)
                && (p_gbparams->stc_nt == (unsigned long)seg_tmp.ns)
                && (p_gbparams->fc_dx == (float)seg_tmp.d2)
                && (p_gbparams->fc_dt == (float)seg_tmp.d1)))
    {
        printf("the input parameters don't match the ones in .su header in file %s\n!"
                , p_gbparams->pc_src);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //copy seg_tmp.data[] to pfcomp_src[], pfcomp_src[] size: size_valid_nf * size_src * size_x2
    //preparation for fft 1d real2halfcomplex
    //size of vector after fft is half size_f = size_t/2 + 1
    DFTI_DESCRIPTOR_HANDLE my_handle;
    MKL_LONG status;
    status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
    status = DftiCommitDescriptor(my_handle);

    fcomp *pfcomp_tmp = NULL;
    pfcomp_tmp = Allocate_1d_floatcomplex(p_gbparams->stc_nt);

    //fft
    for (it=0; it<p_gbparams->stc_nt; it++)
    {
        pfcomp_tmp[it].real = seg_tmp.data[it];
        pfcomp_tmp[it].imag = 0.0;
    }
    status = DftiComputeForward(my_handle, pfcomp_tmp);

    cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
            , p_gbvars->pfcomp_src + seg_tmp.tracf - 1  + p_gbparams->stc_nxtap, p_gbparams->stc_nsrcnx2);


    source_index_tmp = 0;
    while(fvgettr(fp, &seg_tmp))
    {
        //fft
        for (it=0; it<p_gbparams->stc_nt; it++)
        {
            pfcomp_tmp[it].real = seg_tmp.data[it];
            pfcomp_tmp[it].imag = 0.0;
        }
        status = DftiComputeForward(my_handle, pfcomp_tmp);

        if (seg_tmp.fldr != fldr_old)
        {
            source_index_tmp += 1;
            fldr_old = seg_tmp.fldr;
        }

        cblas_ccopy(p_gbparams->stc_valid_nf, pfcomp_tmp + p_gbparams->stc_i_fmin, 1
                , p_gbvars->pfcomp_src + source_index_tmp * p_gbparams->stc_nx2 + seg_tmp.tracf - 1 + p_gbparams->stc_nxtap
                , p_gbparams->stc_nsrcnx2);
    }

    if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
    {
       printf("Error: %s\n", DftiErrorMessage(status));
    }
    fclose(fp);

    printf("input source from file %s is read in successfully, nt = %zu, nfeval = %zu, nx = %zu, dt = %f, dx = %f\n"
            , p_gbparams->pc_src, p_gbparams->stc_nt, p_gbparams->stc_valid_nf
            , p_gbparams->stc_nx, p_gbparams->fc_dt, p_gbparams->fc_dx);
    fflush(stdout);

    #ifdef DEBUG
    fp = fopen("Tmp/source_freq_test.bin", "wb");
    fwrite(p_gbvars->pfcomp_src, sizeof(float) * 2 * p_gbparams->stc_valid_nf
            * p_gbparams->stc_nsrcnx2, 1, fp);
    fclose(fp);
    #endif

    free(pfcomp_tmp);
    pfcomp_tmp = NULL;
}






void Modeling_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars, tmparrays *p_tmpvars
        , const unsigned long iorder, const unsigned long MPI_ID, const unsigned long MPI_SIZE)
{
    FILE *fp;
    char filename[64];
    unsigned long index_ifreq = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;
    unsigned long index_slice_o = 0;
    unsigned long ifreq;
    unsigned long iorder_inner;

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        index_slice_o = (ifreq / MPI_SIZE) * p_gbparams->stc_valid_nv * p_gbparams->stc_nx2pow2;
        index_ifreq = ifreq * p_gbparams->stc_nsrcnx2;

        for (iorder_inner=0; iorder_inner<iorder+1; iorder_inner++)
        {
            Get_Pplus_Pmin_sub(p_gbparams, p_tmpvars, p_gbvars->pfcomp_Pplus_Nslice + index_slice, p_gbvars->pfcomp_Pmin_Nslice + index_slice
                    , p_gbvars->pfcomp_Qplus_Nslice + index_slice
                    , p_gbvars->pfcomp_Wx_Nslice + index_slice_o, p_gbvars->pf_WxGxtaper
                    , p_gbvars->pfcomp_src + index_ifreq, p_gbvars->pst_i_vel, p_gbvars->pf_image);
        }
        #ifdef DEBUG
        snprintf(filename, sizeof(char) * 32, "Tmp/Pmin%zu.bin", ifreq);
        fp = fopen(filename, "wb");
        fwrite(p_gbvars->pfcomp_Pmin_Nslice + index_slice, sizeof(float) * 2
                * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2, 1, fp);
        fclose(fp);
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Write_results_data_2su_mod(const globalconsts *p_gbparams, const fcomp *pfcomp_Pmin_Nslice, const unsigned long MPI_ID, const unsigned long MPI_SIZE, const unsigned long Nslice_tmp)
{
    FILE *fp;
    char filename[64];
    unsigned long ifreq;
    unsigned long isrc;
    unsigned long ix;
    unsigned long it;
    unsigned long impi_id;
    unsigned long index_isrc = 0;
    unsigned long index_ix = 0;
    unsigned long index_slice = 0;
    unsigned long index_slice0 = 0;

    fcomp *pfcomp_Pmin0_Nslice = NULL;
    pfcomp_Pmin0_Nslice = Allocate_1d_floatcomplex(Nslice_tmp * p_gbparams->stc_nsrcnx2);
    fcomp *pfcomp_Pmin0 = NULL;
    fcomp *pfcomp_Pmin0_full = NULL;
    if (MPI_ID == 0)
    {
        pfcomp_Pmin0 = Allocate_1d_floatcomplex(MPI_SIZE * Nslice_tmp * p_gbparams->stc_nsrcnx2);
        pfcomp_Pmin0_full = Allocate_1d_floatcomplex(p_gbparams->stc_nsrcnx2 * p_gbparams->stc_nt);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    for (ifreq=MPI_ID; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
    {
        index_slice = (ifreq / MPI_SIZE) * p_gbparams->stc_nzplus1 * p_gbparams->stc_nsrcnx2;
        index_slice0 = (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
        cblas_ccopy(p_gbparams->stc_nsrcnx2, pfcomp_Pmin_Nslice + index_slice, 1
                , pfcomp_Pmin0_Nslice + index_slice0, 1);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Gather(pfcomp_Pmin0_Nslice, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT
            , pfcomp_Pmin0, 2 * Nslice_tmp * p_gbparams->stc_nsrcnx2, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (MPI_ID == 0)
    {
        for (impi_id=0; impi_id<MPI_SIZE; impi_id++)
        {
            for (ifreq=impi_id; ifreq<p_gbparams->stc_valid_nf_tmp; ifreq+=MPI_SIZE)
            {
                index_slice = impi_id * Nslice_tmp * p_gbparams->stc_nsrcnx2
                    + (ifreq / MPI_SIZE) * p_gbparams->stc_nsrcnx2;
                cblas_ccopy(p_gbparams->stc_nsrcnx2, pfcomp_Pmin0 + index_slice, 1
                        , pfcomp_Pmin0_full + p_gbparams->stc_i_fmin + ifreq, p_gbparams->stc_nt);
            }
        }

        segy seg_tmp;
        memset(&seg_tmp, 0, sizeof(char) * HDRSIZE);
        seg_tmp.trwf = (unsigned long)p_gbparams->stc_nx;
        seg_tmp.ns = (unsigned long)p_gbparams->stc_nt;
        seg_tmp.d2 = p_gbparams->fc_dx;
        seg_tmp.d1 = p_gbparams->fc_dt;
        seg_tmp.timbas = 3;//sorting on header fldr
        seg_tmp.trid = 1;
        //unit is micro second
        seg_tmp.dt = (unsigned long)(p_gbparams->fc_dt * 1000000.0 + 0.1);
        snprintf(filename, sizeof(char) * 64, "%s%s_data0.su", p_gbparams->pc_outfolder, p_gbparams->pc_outlabel);
        fp = fopen(filename, "wb");
        if (fp == NULL)
        {
            printf("error! cannot open file %s to write P file.\n", filename);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //preparation for fft 1d halfcomplex2real
        //size of vector before fft is half size_f = size_t/2 + 1
        DFTI_DESCRIPTOR_HANDLE my_handle;
        MKL_LONG status;
        status = DftiCreateDescriptor(&my_handle, DFTI_SINGLE, DFTI_COMPLEX, 1, p_gbparams->stc_nt);
        status = DftiCommitDescriptor(my_handle);

        for (isrc=0; isrc<p_gbparams->stc_nsrc; isrc++)
        {
            index_isrc = isrc * p_gbparams->stc_nx2 * p_gbparams->stc_nt;
            for (ix=0; ix<p_gbparams->stc_nx; ix++)
            {
                index_ix = index_isrc + (ix + p_gbparams->stc_nxtap) * p_gbparams->stc_nt;
                seg_tmp.fldr = (unsigned long)(isrc + 1);
                seg_tmp.tracl = (unsigned long)(isrc * p_gbparams->stc_nx + ix + 1);
                seg_tmp.tracf = (unsigned long)(ix + 1);
                seg_tmp.f2 = ix * p_gbparams->fc_dx;

                //fft
                //symetric the trace in the frequency domain, for ifft
                for (it=p_gbparams->stc_nf; it<p_gbparams->stc_nt; it++)
                {
                    pfcomp_Pmin0_full[index_ix + it].real = pfcomp_Pmin0_full[index_ix + p_gbparams->stc_nt - it].real;
                    pfcomp_Pmin0_full[index_ix + it].imag = - pfcomp_Pmin0_full[index_ix + p_gbparams->stc_nt - it].imag;
                }
                status = DftiComputeBackward(my_handle, &pfcomp_Pmin0_full[index_ix]);
                for (it=0; it<p_gbparams->stc_nt; it++)
                {
                    seg_tmp.data[it] = pfcomp_Pmin0_full[index_ix + it].real / (float)p_gbparams->stc_nt;
                }

                fwrite(&seg_tmp, sizeof(float)
                        , HDRSIZE / sizeof(float) + p_gbparams->stc_nt, fp);
            }
        }

        if (status && !DftiErrorClass(status,DFTI_NO_ERROR))
        {
           printf("Error: %s\n", DftiErrorMessage(status));
        }
        status = DftiFreeDescriptor(&my_handle);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(pfcomp_Pmin0_full);
    free(pfcomp_Pmin0_Nslice);
    if (MPI_ID == 0)
    {
        free(pfcomp_Pmin0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Free_memory_mod(const globalconsts *p_gbparams, globalarrays *p_gbvars
        , const unsigned long MPI_ID)
{
    free((p_gbvars)->pfcomp_src);
    free((p_gbvars)->pf_WxGxtaper);

    free((p_gbvars)->pfcomp_Wx_Nslice);
    free((p_gbvars)->pfcomp_Gx_Nslice);
    free((p_gbvars)->pfcomp_Pmin_Nslice);
    free((p_gbvars)->pfcomp_Pplus_Nslice);
    free((p_gbvars)->pfcomp_Qplus_Nslice);


    free((p_gbvars)->pf_image);

    free((p_gbvars)->pst_i_vel);
    if (MPI_ID == 0)
    {
        free((p_gbvars)->pf_vel);
        free((p_gbvars)->pf_den);
        free((p_gbvars)->pf_slow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

