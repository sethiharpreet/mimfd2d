/* compute PML coefficients */

// void pml_coeff(float *a_x, float *a_x_half, float *b_x, float *b_x_half, float *K_x, float *K_x_half, float *d_x, float *d_x_half, float *alpha_prime_x, float *alpha_prime_x_half, \
//                float *a_z, float *a_z_half, float *b_z, float *b_z_half, float *K_z, float *K_z_half, float *d_z, float *d_z_half, float *alpha_prime_z, float *alpha_prime_z_half, \
//                float quasi_cp_max, float fpml, float dx, float dz, float dt, int npml, int nx, int nz);

void pml_coeff(float *a_x, float *a_x_half, float *b_x, float *b_x_half, float *K_x, float *K_x_half, \
               float *a_z, float *a_z_half, float *b_z, float *b_z_half, float *K_z, float *K_z_half, \
               float quasi_cp_max, float fpml, float dx, float dz, float dt, int npml, int nx, int nz);


/* ------------------------------------------ Acoustic ----------------------------------------------*/

void update_pressure_fv_PML(float **p_fv, float **vx_ff, float **vz_vv, float **t21, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad,\
                            float **conv_vzvv_fv, float **conv_vxff_fv, float *K_x, float *b_x, float *a_x, float *K_z_half, float *b_z_half, float *a_z_half, \
                            int npml, bool fsrf, bool cpld);


void update_pressure_vf_PML(float **p_vf, float **vz_ff, float **vx_vv, float **t12, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad,\
                            float **conv_vxvv_vf, float **conv_vzff_vf, float *K_z, float *b_z, float *a_z, float *K_x_half, float *b_x_half, float *a_x_half, \
                            int npml, bool fsrf, bool cpld);

void update_ac_velocity_ff_PML(float **vz_ff, float **vx_ff, float ** p_vf, float **p_fv, float **t22, float **ro, float dtx, float dtz, int nzpad, int nxpad,\
                               float **conv_pvf_ff, float **conv_pfv_ff, float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, float *a_z_half, \
                               int npml, bool fsrf, bool cpld);


void update_ac_velocity_vv_PML(float **vz_vv, float **vx_vv, float **p_vf, float **p_fv, float **t11, float **ro, float dtx, float dtz, int nzpad, int nxpad, \
                               float **conv_pfv_vv, float **conv_pvf_vv, float *K_x, float *K_z, float *b_x, float *b_z, float *a_x, float *a_z, \
                               float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, float *a_z_half, int npml, bool fsrf, bool cpld);


/* ------------------------------------------ TTI ----------------------------------------------*/

void update_stress_vf_tti_PML(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv, float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,   \
                              float **c11,    float **c13,    float **c33,    float **c55,   float **c15,    float **c35,    float **t12,    float dtx, float dtz,   \
                              int nzpad, int nxpad,  float **conv_vz_vv_vf,  float **conv_vx_vv_vf, float **conv_vz_ff_vf, float **conv_vx_ff_vf, float *K_z, float *b_z, float *a_z, \
                              float *K_x_half, float *b_x_half, float *a_x_half, int npml, bool fsrf);

void update_stress_fv_tti_PML(float **vz_ff, float **vx_ff, float **vz_vv,  float **vx_vv, float **sxx_fv, float **sxz_fv, float **szz_fv, float **rho, \
                              float **c11, float **c13, float **c33, float **c55, float **c15, float **c35, float **t21, float dtx, float dtz, int nzpad, int nxpad,\
                              float **conv_vz_ff_fv, float **conv_vx_ff_fv, float **conv_vz_vv_fv, float **conv_vx_vv_fv, float *K_x, float *b_x, float *a_x, \
                              float *K_z_half, float *b_z_half, float *a_z_half, int npml, bool fsrf );

void update_velocity_ff_PML(float **vx_ff, float ** vz_ff, float **sxx_fv, float **sxz_fv, float **sxz_vf, float **szz_vf, float **rho, float **t22, \
                            float dtx, float dtz, int nzpad, int nxpad, float **conv_sxx_fv_ff, float **conv_szz_vf_ff, float **conv_sxz_vf_ff, float **conv_sxz_fv_ff,\
                            float *K_x_half, float *b_x_half, float *a_x_half, float *K_z_half, float *b_z_half, float *a_z_half, int npml, bool fsrf);


void update_velocity_vv_PML(float **vx_vv, float **vz_vv, float **szz_fv, float **sxz_fv, float **sxz_vf, float **sxx_vf, float **rho, float **t11, \
                            float dtx, float dtz, int nzpad, int nxpad, float **conv_szz_fv_vv, float **conv_sxz_vf_vv, float **conv_sxx_vf_vv, float **conv_sxz_fv_vv,\
                            float *K_x, float *b_x, float *a_x, float *K_z, float *b_z, float *a_z, int npml, bool fsrf);
