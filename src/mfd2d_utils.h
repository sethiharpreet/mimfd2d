/*< expand domain with free surface at top >*/
void expand2d_bottom(float **in, float **out, fdm2d fdm);

/*< expand domain with free surface at top >*/
void expand2d_top(float **in, float **out, fdm2d fdm );

/*< expand domain with free surface at top >*/
void expand2d(float **in, float **out, fdm2d fdm );

/*< extract domain asssuming free surface at top >*/
void window2d(float **in, float **out, int nz, int nx, int nb);

/*< extract domain asssuming free surface at top >*/
void extract_wfld_fsg(float** p, float** p_ff, float** p_vv, int nz, int nx);

/*< apply absorbing boundary condition >*/
void sponge_ex(float **u, int nxpad, int nzpad, int nb);

/*< extract parameter >*/
void cut2dfree(float** in, float** out, fdm2d fdm, sf_axis cz, sf_axis cx);

/*< Divergence on Dimension 1 >*/
void D1_CC(float **d1, float **dd, int nz, int nx);

/*< Divergence on Dimension 2 >*/
void D2_CC(float **d2, float **dd, int nz, int nx);

/*< Gradient on Dimension 1 >*/
void G1_CC(float **g1,float **gg, int nz, int nx);

/*< Gradient on Dimension 2 >*/
void G2_CC(float **g2,float **gg, int nz, int nx);

/*< Update pressure on fv grids >*/
void update_pressure_fv(float **p_fv, float **vx_ff, float **vz_vv, float **t21, float **rho, float **vel, float dtx, float dtz, int nzpad, int nxpad);

/*< Update pressure on vf grids >*/
void update_pressure_vf(float **p_vf, float **vz_ff, float **vx_vv, float **t12, float **rho, float **vel, float dtx, float dtz, int nzpad, int nxpad);

/*< update acoustic velocity  ff grids >*/
void update_ac_velocity_ff(float **vz_ff, float **vx_ff, float ** p_vf, float **p_fv, float **t22, float **rho, float dtx, float dtz, int nzpad, int nxpad);

/*< update acoustic velocity  vv grids >*/
void update_ac_velocity_vv(float **vz_vv, float **vx_vv, float **p_vf, float **p_fv, float **t11, float **rho, float dtx, float dtz, int nzpad, int nxpad);


/*<----------------------------- VTI ------------------------------------------------------------->*/
/*< update stress vf grids >*/
void update_stress_vf(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                      float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,            \
                      float **c11,    float **c13,    float **c33,    float **c55,            \
                      float **t12,    float dtx,      float dtz,      int nzpad,  int nxpad);


/*< update stress vf grids >*/
void update_stress_fv(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                      float **sxx_fv, float **sxz_fv, float **szz_fv, float **rho,            \
                      float **c11,    float **c13,    float **c33,    float **c55,            \
                      float **t12,    float dtx,      float dtz,      int nzpad, int nxpad);

/*<----------------------------- TTI ------------------------------------------------------------->*/
/*< update stress vf grids >*/
void update_stress_vf_tti(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                          float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,            \
                          float **c11,    float **c13,    float **c33,    float **c55,            \
                          float **c15,    float **c35,    float **t12,    float dtx, float dtz,   \
                          int nzpad,      int nxpad);

/*< update stress vf grids >*/
void update_stress_fv_tti(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                          float **sxx_fv, float **sxz_fv, float **szz_fv, float **rho,            \
                          float **c11,    float **c13,    float **c33,    float **c55,            \
                          float **c15,    float **c35,    float **t21,    float dtx,  float dtz,  \
                          int nzpad,      int nxpad);

/*< update velocity ff grids >*/
void update_velocity_ff(float **vx_ff, float **vz_ff, float **sxx_vf, float **sxz_vf, float **sxz_fv,  \
                        float **szz_fv,float **t22, float **rho, float dtx, float dtz, int nzpad, int nxpad);

/*< update velocity vv grids >*/
void update_velocity_vv(float **vx_vv, float **vz_vv, float **szz_fv, float **sxz_fv, float **sxz_vf,
                        float **sxx_vf, float **rho, float **t11, float dtx, float dtz, int nxpad, int nzpad);

/*< Inject pressure source >*/
void inject_ps_src_cart(float **p_fv,float **p_vf, float*** ww, lint2d cs, int it);

/*< Inject velocity source >*/
void inject_vel_src_cart(float **v_ff,float **v_vv, float*** ww, lint2d cs, int it);

/*< apply vv  vel boundary condition >*/
void apply_vel_cpl_vv(float **vz1_vv, float **vz2_vv, int nz, int nx);

/*< apply vf stress boundary condition >*/
void apply_stress_cpl_vf(float **p1_vf, float **szz_vf, float **sxz_vf, int nz, int nx);

/*< apply ff vel boundary condition >*/
void apply_vel_cpl_ff(float **vz1_ff, float **vz2_ff, float **vx1_ff, float **vx2_ff, \
                      float **vz2_vv, float **vx2_vv, float **vx1_vv, \
                      float **ro, float **vp, float **c13, float **c33, float idz, float idx, \
                      int nz, int nx);

/*< apply fv stress boundary condition >*/
void apply_stress_cpl_fv(float **p_fv, float **szz_fv, float **sxz_fv, float **sxz_vf, \
                         float **ro1, float **ro2, float idz, float idx, int nz, int nx);

/*< apply vf stress boundary condition >*/
void apply_stress_cpl_vf(float **p1_vf, float **szz_vf, float **sxz_vf, int nz, int nx);

/*< combine wavefields >*/
void combine_wfld(float **wout, float** win1, float** win2, int nz1, int nz2, int nx);
