#include "subfunctions.h"

//self document
char *sdoc[] =
{
    "2D Full Wavefield Modeling code",
    "Input files:",
    "velocity=NULL                         input .su velocity model",
    "density=NULL                          input .su density model",
    "source=NULL                           input .su source wavefield in the time domain",

    "Output files:",
    ""
    "outfile_label=fwmod                   add label to the output files",

    "Parameters:",
    "multiple_order=3                      the order of internal multiples considered, =0: only primaries",
    "fmin=5                                starting frequency",
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

    "",
    "Note: The input data must be in .su format",
    "",
    "Author: Shan Qu, Delft University of Technology",
    "Acknowledgement: The main part of this code was written during an internship",
    "                 supervision of Dr. Yimin Sun at Aramco oversea company",
    "First created: September 2018; latest update: April 2019",
    "This code is propriatary under the Delphi Research Consortium",
    "",
    "product: JMI - 2D Full Wavefield Modeling",
    "                                        ",
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
    Get_globalparameters_mod(argc, argv, &gb_params);


    unsigned long Nslice = (unsigned long)ceil((float)gb_params.stc_valid_nf / (float)MPI_SIZE);
    unsigned long Nslice_tmp = Nslice;

    if (MPI_ID == 0)
    {
        Print_globalparameters_mod(&gb_params);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    globalarrays gb_vars;
    Create_memory_mod(&gb_params, &gb_vars, MPI_ID, Nslice);
    tmparrays tmp_vars;
    Create_memory_fortmpvars_mod(&gb_params, &tmp_vars);

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
    unsigned long iter_inner_image;
    unsigned long ismooth;
    unsigned long iorder;
    unsigned long iorder_inner;
    unsigned long ifreq;
    unsigned long istep;
    unsigned long iter_inner;
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
    float f_J_image_full = 0.0;
    float f_J_velocity_old = 0.0;
    float f_J_velocity = 0.0;
    fcomp fcomp_tmp;
    fcomp_tmp.real = 0.0;
    fcomp_tmp.imag = 0.0;

    //**********************************************************************************************
    //****modeling starts from now on, parallel with frequency**************************************
    //**********************************************************************************************
    if (MPI_ID == 0)
    {
        printf("\nmodeling starts from now on:\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (MPI_ID == 0)
    {
        //get models and data from input
        Read_input_extended_model_data_fromsu_mod(&gb_params, &gb_vars);

        for (izx=0; izx<gb_params.stc_nz * gb_params.stc_nx2; izx++)
        {
            gb_vars.pf_slow[izx] = 1.0 / gb_vars.pf_vel[izx];
        }

        Get_current_operatorindex(&gb_params, gb_vars.pst_i_vel, gb_vars.pf_vel);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(gb_vars.pst_i_vel, gb_params.stc_nz * gb_params.stc_nx2
            , MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(gb_vars.pf_image, gb_params.stc_nzplus1 * gb_params.stc_nx2
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

    gb_params.fc_fmax_tmp = gb_params.fc_fmax_upper;
    gb_params.stc_i_fmax_tmp = gb_params.stc_i_fmax_upper;
    gb_params.stc_valid_nf_tmp = gb_params.stc_valid_nf;

    Modeling_mod(&gb_params, &gb_vars, &tmp_vars
            , gb_params.stc_multiple_order, MPI_ID, MPI_SIZE);

    Write_results_data_2su_mod(&gb_params, gb_vars.pfcomp_Pmin_Nslice, MPI_ID, MPI_SIZE, Nslice_tmp);

    Free_memory_mod(&gb_params, &gb_vars, MPI_ID);
    Free_memory_fortmpvars_mod(&gb_params, &tmp_vars);

    if (MPI_ID == 0)
    {
        printf("\nmodeling is DONE\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}


