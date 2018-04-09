#include "../Parameter_files/INIT_PARAMS.H"
#include "../Parameter_files/ANAL_PARAMS.H"
#include "../Parameter_files/HEAT_PARAMS.H"
#include "../Parameter_files/Variables.h"

/*
#define HII_C_INDEX(x,y,z)((unsigned long long)((z)+(HII_MID+1llu)*((y)+HII_D*(x))))// for 3D complex array
#define HII_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(HII_MID+1llu)*((y)+HII_D*(x)))) // for 3D real array with the FFT padding
#define HII_R_INDEX(x,y,z)((unsigned long long)((z)+HII_D*((y)+HII_D*(x)))) // for 3D real array with no padding
*/


/* INDEXING MACROS */
#define LOS_RED_C_INDEX(x,y,z)((unsigned long long)((z)+(DIM_MOCK_OBS_CUBIC_MID+1llu)*((y)+DIM_MOCK_OBS_CUBIC*(x))))// for 3D complex array
#define LOS_RED_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(DIM_MOCK_OBS_CUBIC_MID+1llu)*((y)+DIM_MOCK_OBS_CUBIC*(x)))) // for 3D real array with the FFT padding
#define LOS_RED_R_INDEX(x,y,z)((unsigned long long)((z)+DIM_MOCK_OBS_CUBIC*((y)+DIM_MOCK_OBS_CUBIC*(x)))) // for 3D real array with no padding

#define HII_LOS_RED_C_INDEX(x,y,z)((unsigned long long)((z)+(DIM_MOCK_OBS_MID+1llu)*((y)+DIM_MOCK_OBS*(x))))// for 3D complex array
#define HII_LOS_RED_R_FFT_INDEX(x,y,z)((unsigned long long)((z)+2llu*(DIM_MOCK_OBS_MID+1llu)*((y)+DIM_MOCK_OBS*(x)))) // for 3D real array with the FFT padding
#define HII_LOS_RED_R_INDEX(x,y,z)((unsigned long long)((z)+DIM_MOCK_OBS*((y)+DIM_MOCK_OBS*(x)))) // for 3D real array with no padding

void init_21cmPS_arrays();

int main(int argc, char ** argv){
    
    CUBIC_BOX_LENGTH = atof(argv[2]);
    
    DIM_MOCK_OBS_CUBIC = atof(argv[3]);
    DIM_MOCK_OBS_CUBIC_MID = DIM_MOCK_OBS_CUBIC/2;
    
    DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS = (unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC);
    DIM_MOCK_OBS_CUBIC_TOT_FFT_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*2llu*(DIM_MOCK_OBS_CUBIC_MID+1llu)));
    DIM_MOCK_OBS_CUBIC_KSPACE_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS_CUBIC*DIM_MOCK_OBS_CUBIC*(DIM_MOCK_OBS_CUBIC_MID+1llu)));

    RED_BOX_LENGTH = atof(argv[4]);
    
    DIM_MOCK_OBS = atof(argv[5]);
    DIM_MOCK_OBS_MID = DIM_MOCK_OBS/2;
    
    DIM_MOCK_OBS_TOT_NUM_PIXELS = (unsigned long long)(DIM_MOCK_OBS*DIM_MOCK_OBS*DIM_MOCK_OBS);
    DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS*DIM_MOCK_OBS*2llu*(DIM_MOCK_OBS_MID+1llu)));
    DIM_MOCK_OBS_KSPACE_NUM_PIXELS = ((unsigned long long)(DIM_MOCK_OBS*DIM_MOCK_OBS*(DIM_MOCK_OBS_MID+1llu)));
    
    FILE *F;
    
    char filename[1000];
    int i,j,k,ii,jj;
    
    double ave;
    int ratio, counter;
    int n_x, n_y, n_z;
    float k_x, k_y, k_z, k_mag;
    unsigned long long ct;
    
    double *final_PS_ave, *final_k_ave, *final_PS_error_ave;
    
    strcpy(filename, argv[1]);
    
//    printf("filename = %s\n",filename);

    float *delta_T_LCBox,*delta_T_LCBox_Reduced;

    delta_T_LCBox = (float *) malloc(sizeof(float)*DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS);
    delta_T_LCBox_Reduced = (float *) malloc(sizeof(float)*DIM_MOCK_OBS_TOT_NUM_PIXELS);
    
//    printf("DIM_MOCK_OBS_CUBIC = %d DIM_MOCK_OBS = %d\n",DIM_MOCK_OBS_CUBIC,DIM_MOCK_OBS);
//    printf("DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS = %d DIM_MOCK_OBS_TOT_NUM_PIXELS = %d\n",DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS,DIM_MOCK_OBS_TOT_NUM_PIXELS);
//    printf("DIM_MOCK_OBS_CUBIC_TOT_FFT_NUM_PIXELS = %d DIM_MOCK_CUBIC_OBS_KSPACE_NUM_PIXELS = %d\n",DIM_MOCK_OBS_CUBIC_TOT_FFT_NUM_PIXELS,DIM_MOCK_OBS_CUBIC_KSPACE_NUM_PIXELS);
//    printf("DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS = %d DIM_MOCK_OBS_KSPACE_NUM_PIXELS = %d\n",DIM_MOCK_OBS_TOT_FFT_NUM_PIXELS,DIM_MOCK_OBS_KSPACE_NUM_PIXELS);
    
    if((DIM_MOCK_OBS_CUBIC % DIM_MOCK_OBS)!=0) {
        printf("The L.O.S dimensions of the smaller (asymmetric) box is not an integer multiple of the larger cubic box containing the mock observation!\n");
        printf("This piece of code only works when the ratio is an integer multiple");
        return -1;
    }
    
    F = fopen(filename, "rb");
    for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
        for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
            for (k=0; k<DIM_MOCK_OBS_CUBIC; k++){
                mod_fread(delta_T_LCBox, sizeof(float)*DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS, 1, F);
            }
        }
    }
    fclose(F);
    
    ratio = DIM_MOCK_OBS_CUBIC/DIM_MOCK_OBS;
    
    counter = 0;

    fftwf_complex *deldel_T_FullBox, *deldel_T_ReducedBox;
    fftwf_plan plan;
    
    deldel_T_FullBox = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*DIM_MOCK_OBS_CUBIC_KSPACE_NUM_PIXELS);
    deldel_T_ReducedBox = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex)*DIM_MOCK_OBS_KSPACE_NUM_PIXELS);
    
    init_21cmPS_arrays();
    
    final_PS_ave = calloc(NUM_BINS,sizeof(double));
    final_k_ave = calloc(NUM_BINS,sizeof(double));
    final_PS_error_ave = calloc(NUM_BINS,sizeof(double));
    
    /*
    ave = 0.0;
    for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
        for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
            for (k=0; k<DIM_MOCK_OBS_CUBIC; k++){
                ave += delta_T_LCBox[LOS_RED_R_INDEX(i,j,k)];
            }
        }
    }
    ave /= DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS;
     
    for (i=0; i<DIM_MOCK_OBS_CUBIC; i++){
        for (j=0; j<DIM_MOCK_OBS_CUBIC; j++){
            for (k=0; k<DIM_MOCK_OBS_CUBIC; k++){
                *((float *)deldel_T_FullBox + LOS_RED_R_FFT_INDEX(i,j,k)) = (delta_T_LCBox[LOS_RED_R_INDEX(i,j,k)]/ave - 1)*VOLUME/(DIM_MOCK_OBS_CUBIC_TOT_NUM_PIXELS+0.0);
                if (DIMENSIONAL_T_POWER_SPEC){
                    *((float *)deldel_T_FullBox + LOS_RED_R_FFT_INDEX(i,j,k)) *= ave;
                }
                // Note: we include the V/N factor for the scaling after the fft
            }
        }
    }
     
    plan = fftwf_plan_dft_r2c_3d(DIM_MOCK_OBS_CUBIC, DIM_MOCK_OBS_CUBIC, DIM_MOCK_OBS_CUBIC, (float *)deldel_T_FullBox, (fftwf_complex *)deldel_T_FullBox, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
//    fftwf_cleanup();
    
    for (n_x=0; n_x<DIM_MOCK_OBS_CUBIC; n_x++){
        if (n_x>DIM_MOCK_OBS_CUBIC_MID)
            k_x =(n_x-(int)DIM_MOCK_OBS_CUBIC) * DELTA_K;  // wrap around for FFT convention
        else
            k_x = n_x * DELTA_K;
        
        for (n_y=0; n_y<DIM_MOCK_OBS_CUBIC; n_y++){
            if (n_y>DIM_MOCK_OBS_CUBIC_MID) {
                k_y =(n_y-(int)DIM_MOCK_OBS_CUBIC) * DELTA_K;
            }
            else
                k_y = n_y * DELTA_K;
            
            for (n_z=0; n_z<=DIM_MOCK_OBS_CUBIC_MID; n_z++){
                k_z = n_z * DELTA_K;
                
                k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
     
                // now go through the k bins and update
                ct = 0;
                k_floor = 0;
                k_ceil = k_first_bin_ceil;
                while (k_ceil < k_max){
                    // check if we fal in this bin
                    if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                        in_bin_ct[ct]++;
                        p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_FullBox[LOS_RED_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*VOLUME);
                        // note the 1/VOLUME factor, which turns this into a power density in k-space
                        
                        k_ave[ct] += k_mag;
                        break;
                    }
                    
                    ct++;
                    k_floor=k_ceil;
                    k_ceil*=k_factor;
                }
            }
        }
    }
     
     for (ct=1; ct<NUM_BINS; ct++){
     if (in_bin_ct[ct]>0) {
     printf("%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
     }
     }
    */

    counter = 0;
    
    for(jj=0;jj<ratio;jj++) {
    
    
        for(ii=0;ii<4;ii++) {
        
            for (ct=0; ct<NUM_BINS; ct++){
                p_box[ct] = k_ave[ct] = 0;
                in_bin_ct[ct] = 0;
            }
        
            ave = 0.0;
            for (i=0; i<DIM_MOCK_OBS; i++){
                for (j=0; j<DIM_MOCK_OBS; j++){
                    for (k=0; k<DIM_MOCK_OBS; k++){
                        if(ii==0) {
                            delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)] = delta_T_LCBox[LOS_RED_R_INDEX(i,j,k+counter*DIM_MOCK_OBS)];
                            ave += delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)];
                        }
                        if(ii==1) {
                            delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)] = delta_T_LCBox[LOS_RED_R_INDEX(i+DIM_MOCK_OBS,j,k+counter*DIM_MOCK_OBS)];
                            ave += delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)];
                        }
                        if(ii==2) {
                            delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)] = delta_T_LCBox[LOS_RED_R_INDEX(i,j+DIM_MOCK_OBS,k+counter*DIM_MOCK_OBS)];
                            ave += delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)];
                        }
                        if(ii==3) {
                            delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)] = delta_T_LCBox[LOS_RED_R_INDEX(i+DIM_MOCK_OBS,j+DIM_MOCK_OBS,k+counter*DIM_MOCK_OBS)];
                            ave += delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)];
                        }
                    }
                }
            }
            ave /= DIM_MOCK_OBS_TOT_NUM_PIXELS;
        
            for (i=0; i<DIM_MOCK_OBS; i++){
                for (j=0; j<DIM_MOCK_OBS; j++){
                    for (k=0; k<DIM_MOCK_OBS; k++){
                        *((float *)deldel_T_ReducedBox + HII_LOS_RED_R_FFT_INDEX(i,j,k)) = (delta_T_LCBox_Reduced[HII_LOS_RED_R_INDEX(i,j,k)]/ave - 1)*(RED_BOX_LENGTH*RED_BOX_LENGTH*RED_BOX_LENGTH)/(DIM_MOCK_OBS_TOT_NUM_PIXELS+0.0);
                        if (DIMENSIONAL_T_POWER_SPEC){
                            *((float *)deldel_T_ReducedBox + HII_LOS_RED_R_FFT_INDEX(i,j,k)) *= ave;
                        }
                        // Note: we include the V/N factor for the scaling after the fft
                    }
                }
            }
        
            plan = fftwf_plan_dft_r2c_3d(DIM_MOCK_OBS, DIM_MOCK_OBS, DIM_MOCK_OBS, (float *)deldel_T_ReducedBox, (fftwf_complex *)deldel_T_ReducedBox, FFTW_ESTIMATE);
            fftwf_execute(plan);
            fftwf_destroy_plan(plan);
            //    fftwf_cleanup();
        
            for (n_x=0; n_x<DIM_MOCK_OBS; n_x++){
                if (n_x>DIM_MOCK_OBS_MID)
                    k_x =(n_x-(int)DIM_MOCK_OBS) * (TWOPI/RED_BOX_LENGTH);  // wrap around for FFT convention
                else
                    k_x = n_x * (TWOPI/RED_BOX_LENGTH);
            
                for (n_y=0; n_y<DIM_MOCK_OBS; n_y++){
					if(n_x != 0 && n_y != 0) {
                    if (n_y>DIM_MOCK_OBS_MID) {
                        k_y =(n_y-(int)DIM_MOCK_OBS) * (TWOPI/RED_BOX_LENGTH);
                    }
                    else
                        k_y = n_y * (TWOPI/RED_BOX_LENGTH);
                
                    for (n_z=0; n_z<=DIM_MOCK_OBS_MID; n_z++){
                        k_z = n_z * (TWOPI/RED_BOX_LENGTH);
                    
                        k_mag = sqrt(k_x*k_x + k_y*k_y + k_z*k_z);
                    
                        // now go through the k bins and update
                        ct = 0;
                        k_floor = 0;
                        k_ceil = k_first_bin_ceil;
                        while (k_ceil < k_max){
                            // check if we fal in this bin
                            if ((k_mag>=k_floor) && (k_mag < k_ceil)){
                                in_bin_ct[ct]++;
                                p_box[ct] += pow(k_mag,3)*pow(cabs(deldel_T_ReducedBox[HII_LOS_RED_C_INDEX(n_x, n_y, n_z)]), 2)/(2.0*PI*PI*(RED_BOX_LENGTH*RED_BOX_LENGTH*RED_BOX_LENGTH));
                                // note the 1/VOLUME factor, which turns this into a power density in k-space
                            
                                k_ave[ct] += k_mag;
                                break;
                            }
                        
                            ct++;
                            k_floor=k_ceil;
                            k_ceil*=k_factor;
                        }
                    }
                }
				}
            }
        
            sprintf(filename, "delTps_estimate_carvedup_section_individual%d_%d_%i_%.0fMpc_lighttravel.txt",ii,jj, HII_DIM, BOX_LEN);
            F=fopen(filename, "wt");
            for (ct=1; ct<NUM_BINS; ct++){
                if (in_bin_ct[ct]>0) {
                    fprintf(F, "%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
                }
            }
            fclose(F);
            
            
            for (ct=1; ct<NUM_BINS; ct++){
                if (in_bin_ct[ct]>0) {
//                    printf("%e\t%e\t%e\n", k_ave[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0), p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0));
                    final_k_ave[ct] += k_ave[ct]/(in_bin_ct[ct]+0.0);
                    final_PS_ave[ct] += p_box[ct]/(in_bin_ct[ct]+0.0);
                    final_PS_error_ave[ct] += p_box[ct]/(in_bin_ct[ct]+0.0)/sqrt(in_bin_ct[ct]+0.0);
                }
            }
        
            
            for (i=0; i<DIM_MOCK_OBS; i++){
                for (j=0; j<DIM_MOCK_OBS; j++){
                    for (k=0; k<(DIM_MOCK_OBS_MID+1); k++){
                        *((float *)deldel_T_ReducedBox + HII_LOS_RED_R_FFT_INDEX(i,j,k)) = 0.0;
                    }
                }
            }
        }
    
        counter += 1;
        
        sprintf(filename, "delTps_estimate_carvedup_section_averaged_%d_%i_%.0fMpc_lighttravel.txt",jj, HII_DIM, BOX_LEN);
        F=fopen(filename, "wt");
        for (ct=1; ct<NUM_BINS; ct++){
            if (in_bin_ct[ct]>0) {
                fprintf(F, "%e\t%e\t%e\n", final_k_ave[ct]/4., final_PS_ave[ct]/4., final_PS_error_ave[ct]/4.);
            }
        }
        fclose(F);
        
        for (ct=0; ct<NUM_BINS; ct++){
            final_k_ave[ct] = 0.;
            final_PS_ave[ct] = 0.;
            final_PS_error_ave[ct] = 0.;
        }
    }
    
    
  
    return 0;
}

void init_21cmPS_arrays() {
    
    //    k_factor = 1.25;
    k_factor = 1.35;
    //    k_factor = 1.5;
    //    k_factor = 2.0;
    k_first_bin_ceil = (TWOPI/RED_BOX_LENGTH);
    k_max = (TWOPI/RED_BOX_LENGTH)*(int)DIM_MOCK_OBS;
//    k_first_bin_ceil = DELTA_K;
//    k_max = DELTA_K*(int)DIM_MOCK_OBS_CUBIC;
    // initialize arrays
    // ghetto counting (lookup how to do logs of arbitrary bases in c...)
    NUM_BINS = 0;
    k_floor = 0;
    k_ceil = k_first_bin_ceil;
    while (k_ceil < k_max){
        NUM_BINS++;
        k_floor=k_ceil;
        k_ceil*=k_factor;
    }
    
    p_box = calloc(NUM_BINS,sizeof(double));
    k_ave = calloc(NUM_BINS,sizeof(double));
    in_bin_ct = (unsigned long long *)calloc(NUM_BINS,sizeof(unsigned long long));
    
}


