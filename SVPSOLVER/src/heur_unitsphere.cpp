#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"

using namespace std;

#define debug  0

void  SVPsolver::SVPSheurUnitsphere(
      const NODE&    node,
      const double*  vars_localub,
      const double*  vars_locallb
      )
{

   int      m = probdata.get_m();

   double   val;
   double   val_new;
   double   *solvals = nullptr;
   double   *solvals_new = nullptr;

   double   *relaxsolvals = node.get_relaxsolval();

   assert( m > 0 );
   assert( relaxsolvals != nullptr );

   solvals = new double[m];
   solvals_new = new double[m];

   for(int i=0; i<m; i++){
      solvals[i] = round( relaxsolvals[i] );
   }

   val = compute_objval( solvals );

   double   lmin;
   double   *lsol = new double[m];

   //cout << "best:" << bestval << endl;
   //cout << "val:" << val << endl;
   while(1){

      val_new = val;

      lmin = val;
      Copy_vec( solvals, lsol, m);

      for(int i=0; i<m; i++){
         lsol[i] += 1;
         if( lsol[i] <= vars_localub[i] ){
            lmin = compute_objval( lsol );
            if( val_new > lmin ){
               val_new = lmin;
               Copy_vec( lsol, solvals_new, m);
            }
         }
         lsol[i] -= 2;
         if( vars_locallb[i] <= lsol[i] ){
            lmin = compute_objval( lsol );
            if( val_new > lmin ){
               val_new = lmin;
               Copy_vec( lsol, solvals_new, m);
            }
         }
         lsol[i] += 1;
      }

      //cout << "val_new:" << val_new << endl;
      if( val > val_new ){
         val = val_new;
         Copy_vec( solvals_new, solvals, m);
      }else{
         //cout << "break" << endl;
         break;
      }
   }

   SOLUTION solution;
   solution.set_sol( m, solvals, val);
   pool.add_solution( solution );

   if( bestval > val ){
      bestval = val;
      bestsol = solution;
      Appfac = sqrt( bestval ) / _Appfac;
#if debug
      cout << "*********************************" << endl;
      cout << "*   heur_unitsphere  dpt: " << node.get_dpt() << endl;
      cout << "*********************************" << endl;
#endif
   }

   delete[] solvals;
   delete[] solvals_new;
   delete[] lsol;
}
