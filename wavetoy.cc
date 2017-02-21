/*
Function prototype to call ET solvers:
1. "real" is set to either "float" or "double" depending on the compilation flag
2. "NCOMP_TOTAL" is the total number of variables (e.g., 25)
3. "FLU_NXT" is the number of cells along each direction  including ghost zones. It is equal to "PS2 + 2*FLU_GHOST_SIZE", where PS2 is the block size excluding ghost zones and  FLU_GHOST_SIZE is the number of ghost zones on each side along each direction. Both FLU_NXT, PS2, and FLU_GHOST_SIZE are symbolic constants. For example, we can have PS2=16 and FLU_GHOST_SIZE=5, which leads to FLU_NXT=26.
4. dt: evolution time-step
5. dh: cell size
*/

void ET_Solver( const real Input[NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real Output[NCOMP_TOTAL][ PS2*PS2*PS2 ],
                const real dt, const real dh ) {
  /* a simple wave equation in a first order time, second order space
   * decomposition.
   * for the PDE 0 = (\partial^2_t + \partial^2_x) \phi we introduce a new
   * variable \pi for the first time derivative \pi = \partial_t \phi and
   * arrive at a system of PDE:
   * \partial_t \pi  = \partial^2_x \phi
   * \partial_t \phi = \pi
   */
/* helper macro to index into the 3d cell array */
#define GFINDEX3D_GHOSTED(i,j,k) ((k+FLU_GHOST_SIZE)*FLU_NXT*FLU_NXT+(j+FLU_GHOST_SIZE)*FLU_NXT+(i+FLU_GHOST_SIZE))
#define GFINDEX3D_NONGHOSTED(i,j,k) (k*PS2*PS2+j*PS2+i)
#define PD82nd(f,idx,didx) (+(1./280.)*f[idx-4*didx] \
                            -(4./105.)*f[idx-3*didx] \
                            +(1./5.)  *f[idx-2*didx] \
                            -(4./5.)  *f[idx-1*didx] \
                            +(4./5.)  *f[idx-1*didx] \
                            -(1./5.)  *f[idx-2*didx] \
                            +(4./105.)*f[idx-3*didx] \
                            -(1./280.)*f[idx-4*didx])
  real *phi = Input[0];
  real *pi = Input[1];
  real *dtphi = Output[0];
  real *dtpi = Output[1];
  const ptrdiff_t di = GFINDEX3D_GHOSTED(1,0,0);
  const ptrdiff_t dj = GFINDEX3D_GHOSTED(0,1,0);
  const ptrdiff_t dk = GFINDEX3D_GHOSTED(0,0,1);
  for(ptrdiff_t k = 0 ; k < PS2 ; ++k) {
    for(ptrdiff_t j = 0 ; j < PS2 ; ++j) {
      for(ptrdiff_t i = 0 ; i < PS2 ; ++i) {
        const ptrdiff_t idx_in = GFINDEX3D_GHOSTED(i,j,k)
        const ptrdiff_t idx_out = GFINDEX3D_NONGHOSTED(i,j,k)
        const real phi11 = PD82nd(phi, idx_in,di);
        const real phi22 = PD82nd(phi, idx_in,dj);
        const real phi33 = PD82nd(phi, idx_in,dk);
        const real phirhs = pi[idx_in];
        const real pirhs = phi11 + phi22 + phi33;
        dtphi[idx_out] = phirhs;
        dtpi[idx_out] = pirhs;
      }
    }
  }
}
				
