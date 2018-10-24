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
#define log    0

void  SVPsolver::SVPSheurQuadratic(
      const NODE&    node,
      const double*  vars_localub,
      const double*  vars_locallb
      )
{
   assert( vars_localub != nullptr );
   assert( vars_locallb != nullptr );

   // copy
   const auto m = probdata.get_m();
   const auto B_ = probdata.get_B_();
   const auto* colnorm = norm;
   const auto ep = epsilon;

   assert( m > 0 );
   assert( B_ != nullptr );

   double val;
   double val_new;
   double* solvals = nullptr;
   double* solvals_new = nullptr;

   const auto relaxsolvals = node.get_relaxsolval();

   assert( relaxsolvals != nullptr );

   double lam;
   double lam1;
   double lam2;
   double t;
   double t1;
   double t2;
   double s1;
   double s2;
   double* normv = nullptr;
   double* Bx = nullptr;

   vector<int> nofixed_index;

   nofixed_index.reserve( m );
   for ( auto i = 0; i < m; ++i )
   {
      if ( !Equal( vars_localub[i], vars_locallb[i], ep )  )
         nofixed_index.push_back( i );
   }

   nofixed_index.shrink_to_fit();

   assert( colnorm != nullptr );

   solvals = new double[m];
   solvals_new = new double[m];
   normv = new double[m];
   Bx = new double[m];

   for( auto i = 0; i < m; ++i)
      solvals[i] = round( relaxsolvals[i] );

   val = compute_objval( solvals );

#if debug
   cout << "best:" << bestval << endl;
   cout << "val:" << val << endl;
#endif

   for( auto i = 0; i < m; ++i )
      normv[i] = colnorm[i] * colnorm[i];

   while(1)
   {
      Com_mat_Ax( B_, m, m, solvals, Bx );
      val_new = val;

      for ( auto i : nofixed_index )
      {
         double buf;
         lam = val;
         s1 = Com_dot( &B_[i*m], Bx, m);
         s2 = normv[i]; //norm[i] * norm[i];
         t = - s1 / s2;

         if ( t < vars_locallb[i] - solvals[i] )
         {
            t = vars_locallb[i] - solvals[i];
            //lam += t * ( t * s2 + 2.0 * s1 );
            buf = t;
            buf *= s2;
            buf += 2.0 * s1;
            buf *= t;
            lam += buf;
         }
         else if ( vars_localub[i] - solvals[i] < t )
         {
            t = vars_localub[i] - solvals[i];
            //lam += t * ( t*s2 + 2*s1 );
            buf = t;
            buf *= s2;
            buf += 2.0 * s1;
            buf *= t;
            lam += buf;
         }
         else
         {
            t1 = ceil( t );
            lam1 = val;
            //lam1 += t1 * ( t1*s2 + 2*s1 );
            buf = t1;
            buf *= s2;
            buf += 2.0 * s1;
            buf *= t1;
            lam1 += buf;

            t2 = floor( t );
            lam2 = val;
            //lam2 += t2 * ( t2*s2 + 2*s1 );
            buf = t2;
            buf *= s2;
            buf += 2.0 * s1;
            buf *= t2;
            lam2 += buf;

            if( lam1 > lam2 )
            {
               lam = lam2;
               t = t2;
            }
            else
            {
               lam = lam1;
               t = t1;
            }
         }

         #if debug
            cout << "i:" << i << "-> " << lam << endl;
         #endif

         if( val_new > lam )
         {
            val_new = lam;
            Copy_vec( solvals, solvals_new, m );
            solvals_new[i] = solvals[i] + t;
         }
      }

      #if debug
         cout << "val_new:" << val_new << endl;
      #endif

      // val and val_new are integer.
      if( val - val_new > 0.1)
      {
         //cout << "val:" << val << " val_new:" << val_new << endl;
         val = val_new;
         Copy_vec( solvals_new, solvals, m);
      }
      else
      {
         #if debug
            cout << "break" << endl;
         #endif
         break;
      }
   }

   SOLUTION solution;
   solution.set_sol( m, solvals, val);
   SVPStrySol( solution, true, true, nullptr );
   //bool test = false;
   //SVPStrySol( solution, true, true, &test );
   //if ( test )
   //   cout << node.get_index() << endl;

   delete[] solvals;
   delete[] solvals_new;
   delete[] Bx;
   delete[] normv;
}
