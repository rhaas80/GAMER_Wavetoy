#if 0
g++ -g -O3 -Wall main.cc wavetoy.cc -o wavetoy
exit 0
#endif
 
#include <cstdio>
#include <cmath>

#include "grid.h"

void ET_Solver( const real Input[NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real Output[NCOMP_TOTAL][ PS2*PS2*PS2 ],
                const real dt, const real dh );

int main(void) {
  static real psi[NCOMP_TOTAL][PS2][PS2][PS2];
  static real rhs[NCOMP_TOTAL][PS2][PS2][PS2];

  static real psi_ghosted[NCOMP_TOTAL][FLU_NXT][FLU_NXT][FLU_NXT];

  // just some random values
  const real dt = 0.1;
  const real dh = 0.5;

  // some random ID
  for(int k = 0 ; k < PS2 ; k++) {
    for(int j = 0 ; j < PS2 ; j++) {
      for(int i = 0 ; i < PS2 ; i++) {
        psi[0][k][j][i] = sin(i*2*M_PI/PS2)*sin(j*2*M_PI/PS2)*sin(k*2*M_PI/PS2);
        psi[1][k][j][i] = 0.;

      }
    }
  }

  for(int it = 0 ; it < 100 ; it++) {
    // output
    if(it % 10 == 0) {
      for(int k = 0 ; k < 1 ; k++) { // 2d output for now only
        for(int j = 0 ; j < PS2 ; j++) {
          for(int i = 0 ; i < PS2 ; i++) {
            printf("%d %d %d %e %e\n", i,j,k, psi[0][k][j][i], psi[1][k][j][i]);
          }
          printf("\n");
        }
      }
      printf("\n");
    }

    // created ghosted array for RHS computation
    for(int c = 0 ; c < NCOMP_TOTAL ; c++) {
      for(int k = -FLU_GHOST_SIZE ; k < PS2+FLU_GHOST_SIZE ; k++) {
        for(int j = -FLU_GHOST_SIZE ; j < PS2+FLU_GHOST_SIZE ; j++) {
          for(int i = -FLU_GHOST_SIZE ; i < PS2+FLU_GHOST_SIZE ; i++) {
            // this enforces periodic boundary conditions
            psi_ghosted[c][k+FLU_GHOST_SIZE][j+FLU_GHOST_SIZE][i+FLU_GHOST_SIZE] =
              psi[c][(k+PS2)%PS2][(j+PS2)%PS2][(i+PS2)%PS2];

          }
        }
      }
    }
    
    ET_Solver((real (*)[FLU_NXT*FLU_NXT*FLU_NXT])&psi_ghosted[0][0][0][0],
              (real (*)[PS2*PS2*PS2])&rhs[0][0][0][0], dt, dh);

    // 1st order explicit Euler scheme. Sue me.
    for(int c = 0 ; c < NCOMP_TOTAL ; c++) {
      for(int k = 0 ; k < PS2 ; k++) {
        for(int j = 0 ; j < PS2 ; j++) {
          for(int i = 0 ; i < PS2 ; i++) {
            psi[c][k][j][i] += dt * rhs[c][k][j][i];
          }
        }
      }
    }

  }
}
