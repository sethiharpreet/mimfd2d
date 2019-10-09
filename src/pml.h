/* compute PML coefficients */

void pml_coeff(float *a_x, float *a_x_half, float *b_x, float *b_x_half, float *K_x, float *K_x_half, float *d_x, float *d_x_half, float *alpha_prime_x, float *alpha_prime_x_half, \
               float *a_z, float *a_z_half, float *b_z, float *b_z_half, float *K_z, float *K_z_half, float *d_z, float *d_z_half, float *alpha_prime_z, float *alpha_prime_z_half, \
               float quasi_cp_max, float fpml, float dx, float dz, float dt, int npml, int nx, int nz);

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
