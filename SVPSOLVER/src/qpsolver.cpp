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
#include <time.h>
#include <iomanip>
#include <omp.h>

#include "qpsolver.h"
#include "vector.h"

using namespace std;


#define debug_class  0
#define debug  0
#define RUNTIME 0

QPsolver::QPsolver(){   // default constructor
#if debug_class
   cout << "QPsolver: default constructor" << endl;
#endif
   Q = nullptr;
   p = nullptr;
   A = nullptr;
   b = nullptr;
   C = nullptr;
   d = nullptr;
   l = nullptr;
   u = nullptr;
   n = 0;
   m = 0;
   k = 0;
   warm = nullptr;
   bestsol = nullptr;
   bestval = 0.0;
   H = nullptr;
   r = nullptr;
   z = nullptr;
   indexes = nullptr;
   fixed = nullptr;
   comls = nullptr;
}

QPsolver::QPsolver( const QPsolver &source )
{  // copy constructor
#if debug_class
   cout << "QPsolver: copy constructor" << endl;
#endif
   Q = source.Q;
   p = source.p;
   A = source.A;
   b = source.b;
   C = source.C;
   d = source.d;
   l = source.l;
   u = source.u;
   n = source.n;
   m = source.m;
   k = source.k;
   warm = source.warm;
   bestval = source.bestval;
   if( source.bestsol != nullptr ){
      assert( n > 0 );
      bestsol = new double[n];
      for(int i=0; i<n; i++){
         bestsol[i] = source.bestsol[i];
      }
   }else{
      bestsol = nullptr;
   }
   H = nullptr;
   r = nullptr;
   z = nullptr;
   indexes = nullptr;
   fixed = nullptr;
   comls = nullptr;
}

// assignment operator
QPsolver& QPsolver::operator=( const QPsolver& source )
{
#if debug_class
   cout << "QPsolver: assignment operator" << endl;
#endif

   if( this != &source ){
      Q = source.Q;
      p = source.p;
      A = source.A;
      b = source.b;
      C = source.C;
      d = source.d;
      l = source.l;
      u = source.u;
      n = source.n;
      m = source.m;
      k = source.k;
      warm = source.warm;
      bestval = source.bestval;
      if( source.bestsol != nullptr ){
         assert( n > 0 );
         bestsol = new double[n];
         for(int i=0; i<n; i++){
            bestsol[i] = source.bestsol[i];
         }
      }else{
         bestsol = nullptr;
      }
      H = nullptr;
      r = nullptr;
      z = nullptr;
      indexes = nullptr;
      fixed = nullptr;
      comls = nullptr;
   }
   return *this;
}

// destructor
QPsolver::~QPsolver()
{
#if debug_class
   cout << "QPsolver: destructor" << endl;
#endif
   delete[] bestsol;
   delete[] H;
   delete[] r;
   delete[] indexes;
   delete[] z;
   delete[] fixed;
   bestsol = nullptr;
   H = nullptr;
   r = nullptr;
   z = nullptr;
   fixed = nullptr;
   indexes = nullptr;
}

void  QPsolver::disp_prob(){

   if( n > 0 ){
      cout << "dim: " << n << endl;

      if( Q ){
         cout << "Q:" << endl;
         printM( n, n, Q);
      }
      if( p ){
         cout << "p:" << endl;
         printv( n, p);
      }
      if( m ){
         cout << "m: " << m << endl;
         if( A ){
            cout << "A:" << endl;
            printM( m, n, A);
         }
         if( b ){
            cout << "b:" << endl;
            printv( m, b);
         }
      }
      if( k ){
         cout << "k: " << k << endl;
         if( C ){
            cout << "C:" << endl;
            printM( k, n, C);
         }
         if( d ){
            cout << "d:" << endl;
            printv( k, d);
         }
      }
      if( u ){
         cout << "u:" << endl;
         printv( n, u);
      }
      if( l ){
         cout << "l:" << endl;
         printv( n, l);
      }
      if( warm ){
         cout << "init:" << endl;
         printv( n, warm);
      }
   }

}

void  QPsolver::solve()
{
   assert( n > 0 );
   assert( warm != nullptr );
   assert( Q != nullptr );

   #if RUNTIME
   {
      clock_t start;
      clock_t end;
      start = clock();
      cout << setprecision(20);
   }
   #endif

   double   ep = get_ep();
   double   *x;
   x = new double[n];

   Copy_vec( warm, x, n );
   //for ( int i = 0; i < n; i++ ){
   // x[i] = warm[i];
   //}

   #if debug
      cout << "initial point:";
      for(int i=0; i<n; i++){
         cout << x[i] << " ";
      }
      cout << endl;
   #endif

   // find active constraint sets {{
   bool  *W_ineq = nullptr;
   if( m > 0 ){
      cout << "You must fix qpsolver.cpp" << endl;
      exit(-1);
      assert( A != nullptr );
      assert( b != nullptr );
      W_ineq   =  new bool[m];

      for(int i=0; i<m; i++){
         double buf = Com_dot( A+(i*n), x, n);
         if( Equal( b[i], buf, ep ) == true ){
            W_ineq[i] = true;
         }else if( b[i] - buf > 0 ){
            W_ineq[i] = false;
         }else{
            cout << "the given initial point is not a feasible solution." << endl;
            cout << b[i] - buf << endl;
            disp_prob();
            exit(-1);
         }
      }

   }

   #if debug
      cout << "active constraints:";
      if( W_ineq != nullptr ){
         for(int i=0; i<m; i++){
            if( W_ineq[i] == true ){
               cout << i << " ";
            }
         }
         cout << endl;
      }else{
         cout << "There are not inequality constraints" << endl;
      }
   #endif
   // }} find active constraint sets

   // check feasibility for equation constraints{
   if( k > 0 ){
      cout << "You must fix qpsolver.cpp" << endl;
      exit(-1);
      assert( C != nullptr );
      assert( d != nullptr );

      for(int i=0; i<k; i++){
         if( Equal( Com_dot( C+(i*n), x, n), d[i], ep) == false ){
            cout << "the given initial point is not a feasible solution." << endl;
            cout << "aa" << endl;
            exit(-1);
         }
      }
   }
   // }} check feasibility for equation constraints

   // find active bounds {{
   bool  *W_u = nullptr;
   bool  *W_l = nullptr;
   if( u != nullptr )
   {
      W_u = new bool[n];

      double   buf;
      double   u_i;
      double*  x_i;

      for ( int i = 0; i < n; i++ )
      {
         u_i = u[i];
         x_i = &x[i];
         buf = u_i - *x_i;

         if ( Equal( u_i, *x_i, ep) == true )
         {
            W_u[i] = true;
            if( buf != 0 )
               *x_i = u_i;
         }
         else if ( buf > 0)
         {
            W_u[i] = false;
         }
         else
         {
            cout << "the given initial point is not a feasible solution." << endl;
            cout << "bb:" << buf << endl;
            cout << "point:";
            printv( n, warm);
            cout << "l:";
            printv( n, l);
            cout << "u:";
            printv( n, u);
            exit(-1);
         }

      }
   }

   if( l != nullptr )
   {
      W_l = new bool[n];

      double   buf;
      double   l_i;
      double*  x_i;

      for( int i = 0; i < n; i++ )
      {
         l_i = l[i];
         x_i = &x[i];
         buf = *x_i - l_i;

         if ( Equal( *x_i, l_i, ep) )
         {
            W_l[i] = true;
            if( buf != 0 )
               *x_i = l_i;
         }
         else if ( buf > 0 )
         {
            W_l[i] = false;
         }
         else
         {
            cout << "the given initial point is not a feasible solution." << endl;
            cout << "cc:" << buf << endl;
            exit(-1);
         }
      }
   }

   #if debug
      cout << "active upper bounds:";
      if( W_u != nullptr ){
         for(int i=0; i<n; i++){
            if( W_u[i] == true ){
               cout << i << " ";
            }
         }
         cout << endl;
      }else{
         cout << "There are not upper bounds" << endl;
      }
      cout << "active lower bounds:";
      if( W_l != nullptr ){
         for(int i=0; i<n; i++){
            if( W_l[i] == true ){
               cout << i << " ";
            }
         }
         cout << endl;
      }else{
         cout << "There are not lower bounds" << endl;
      }
   #endif
   // }} find active bounds

   double   *x_new = nullptr;
   double   *dx = nullptr;
   bool     update;

   x_new = new double[n];
   dx = new double[n];

   #if RUNTIME
   {
      end = clock();
      cout << "while_start: \t" << (double) (end - start) / CLOCKS_PER_SEC << "sec." << endl;
   }
   #endif

   int   debug_ct = 0;
   int   argmin;
   bool  onlybounds = false;

   if( A == nullptr && C == nullptr && l != nullptr && u != nullptr )
   {
      assert( W_ineq == nullptr );
      assert( H == nullptr );
      assert( r == nullptr );
      assert( indexes == nullptr );

      onlybounds = true;
      H = new double[n*n];
      r = new double[n];
      z = new double[n];
      indexes = new int[n];
      fixed = new double[n];
      comls = Com_LS_dposv;
   }
   else
   {
      cout << "You must fix qpsolver.cpp" << endl;
      exit(-1);
   }

   while ( 1 )
   {
      if ( debug_ct > 1000 )
      {
         cout << "error: qpsolver.cpp" << endl;
         cout << "th: " << omp_get_thread_num() << endl;
         exit(1);
      }
      debug_ct++;

      #if debug
      {
         cout << "active constraints:";
         if( W_ineq != nullptr ){
            for(int i=0; i<m; i++){
               if( W_ineq[i] == true ){
                  cout << i << " ";
               }
            }
            cout << endl;
         }else{
            cout << "There are not inequality constraints" << endl;
         }
         cout << "active upper bounds:";
         if( W_u != nullptr ){
            for(int i=0; i<n; i++){
               if( W_u[i] == true ){
                  cout << i << " ";
               }
            }
            cout << endl;
         }else{
            cout << "There are not upper bounds" << endl;
         }
         cout << "active lower bounds:";
         if( W_l != nullptr ){
            for(int i=0; i<n; i++){
               if( W_l[i] == true ){
                  cout << i << " ";
               }
            }
            cout << endl;
         }else{
            cout << "There are not lower bounds" << endl;
         }
      }
      #endif

      //argmin = solve_activeset( x_new, W_ineq, W_u, W_l);
      //cout << "argmin: " << argmin << ", old_func.: ";
      //printv( n, x_new );
      //argmin = solve_activeset_onlybounds( x_new, W_u, W_l);
      //cout << "argmin: " << argmin << ", new_func.: ";
      //printv( n, x_new );
      //assert(0);
      // solve the problem with active set
      if( onlybounds == true )
      {
         assert( W_ineq == nullptr );
         argmin = solve_activeset_onlybounds( x_new, W_u, W_l);
      }
      else
      {
         cout << "You must fix qpsolver.cpp" << endl;
         exit(-1);
         argmin = solve_activeset( x_new, W_ineq, W_u, W_l);
      }

      #if debug
      {
         cout << "argmin: " << argmin << endl;
      }
      #endif

      // dx = xnew - x
      Com_linecomb( x_new, x, n, 1.0, -1.0, dx);

#if debug
   cout << "dx: ";
   printv( n, dx);
#endif
      if( Com_nrm( dx, n) > ep ){
         // compute stepsize
         double alpha = compute_stepsize( x, dx, W_ineq, W_u, W_l);
#if debug
   cout << "alpha: " << alpha << endl;
#endif
         assert( alpha >= 0 );
         assert( alpha <= 1 );
         // x_new = x + alpha dx
         Com_linecomb( x, dx, n, 1.0, alpha, x_new);
#if debug
   cout << "x_new: ";
   printv( n, x_new);
#endif
      }

      /* x_new is feasible for the original problem now */

      // update W_ineq, W_u and W_l
      update = update_activeset( x_new, W_ineq, W_u, W_l);

      if( update == true ){
         Copy_vec( x_new, x, n);
         continue;
      }

      if( argmin == -1 ){
         break;
      }else{
         int ct = 0;
         update = false;
         if( W_ineq != nullptr ){
            assert( m > 0 );
            if( count( m, W_ineq) > argmin ){
               for(int i=0; i<m; i++){
                  if( W_ineq[i] == true ){
                     if( ct == argmin ){
                        W_ineq[i] = false;
                        update = true;
                        break;
                     }else{
                        ct++;
                     }
                  }
               }
            }else{
               ct += count( m, W_ineq);
            }
         }

         if( update == true ){
            Copy_vec( x_new, x, n);
            continue;
         }

         if( W_u != nullptr ){
            if( ct + count( n, W_u) > argmin ){
               for(int i=0; i<n; i++){
                  if( W_u[i] == true ){
                     if( ct == argmin ){
                        W_u[i] = false;
                        update = true;
                        break;
                     }else{
                        ct++;
                     }
                  }
               }
            }else{
               ct += count( n, W_u);
            }
         }

         if( update == true ){
            Copy_vec( x_new, x, n);
            continue;
         }

         if( W_l != nullptr ){
            for(int i=0; i<n; i++){
               if( W_l[i] == true ){
                  if( ct == argmin ){
                     W_l[i] = false;
                     update = true;
                     break;
                  }else{
                     ct++;
                  }
               }
            }
         }

         if( update == true ){
            Copy_vec( x_new, x, n);
            continue;
         }else{
            cout << "error: QPsolver" << endl;
            cout << "th: " << omp_get_thread_num() << endl;
            exit(-1);
         }

      }
   }

   #if RUNTIME
   {
      end = clock();
      cout << "while_end: \t" << (double) (end - start) / CLOCKS_PER_SEC << "sec." << endl;
   }
   #endif

   bestsol = new double[n];

   Copy_vec( x_new, bestsol, n);

   bestval = compute_objval( bestsol );

#if debug
   cout << "bestsol: ";
   printv( n, bestsol);
   cout << "bestval: " << bestval << endl;
#endif

   delete[] x;
   delete[] W_ineq;
   delete[] W_u;
   delete[] W_l;
   delete[] x_new;
   delete[] dx;

   #if RUNTIME
   {
      end = clock();
      cout << "end: \t\t" << (double) (end - start) / CLOCKS_PER_SEC << "sec." << endl;
      cout << endl;
   }
   #endif
}

int QPsolver::solve_activeset(
   double   *x_new,
   bool     *W_ineq,
   bool     *W_u,
   bool     *W_l
   )
{
   /**********************************************
      min   x'Qx + p'x
      s.t.  F x = g

      by solving H z = r.

      This is

      | 2Q , F' || x | = | -p |
      | F  , 0  || y |   |  g |

   **********************************************/
   assert( x_new != nullptr );

   int   w_equ = k;
   int   w_ineq = count( m, W_ineq);
   int   w_u = count( n, W_u);
   int   w_l = count( n, W_l);
   int   w = w_equ + w_ineq + w_u + w_l;

   #if debug
   {
      cout << "w: " << w << endl;
   }
   #endif

   double   *F = nullptr;  // [w*n]
   double   *g = nullptr;  // [w]

   if( w > 0 ){
      F = new double[w*n];
      g = new double[w];
   }

   /**********************************************
      1. C
      2. A
      3. u
      4. l
   **********************************************/

   // 1. C
   int      ct=0;
   if( w_equ > 0 ){
      assert( k > 0 );
      assert( C != nullptr );
      assert( w >= k );
      assert( F != nullptr );

      Copy_vec( C, F, k*n);
      Copy_vec( d, g, k);

      ct += k;

      assert( w >= ct );
   }

   // 2. A
   if( w_ineq > 0 ){
      assert( m > 0 );
      assert( A != nullptr );
      assert( w >= w_ineq );
      assert( F != nullptr );

      for(int i=0; i<m; i++){
         if( W_ineq[i] == true ){
            Copy_vec( &A[i*n], &F[ct*n], n);
            g[ct] = b[i];
            ct++;
         }
      }

      assert( w >= ct );
   }

   // 3. u
   if( w_u > 0 ){
      assert( u != nullptr );
      assert( w >= w_u );
      assert( F != nullptr );

      for(int i=0; i<n; i++){
         if( W_u[i] == true ){
            for(int j=0; j<n; j++){
               if( j == i ){
                  F[(ct*n)+j] = 1.0;
               }else{
                  F[(ct*n)+j] = 0.0;
               }
            }
            g[ct] = u[i];
            ct++;
         }
      }

      assert( w >= ct );
   }

   // 4. l
   if( w_l > 0 ){
      assert( l != nullptr );
      assert( w >= w_l );
      assert( F != nullptr );

      for(int i=0; i<n; i++){
         if( W_l[i] == true ){
            for(int j=0; j<n; j++){
               if( j == i ){
                  F[(ct*n)+j] = -1.0;
               }else{
                  F[(ct*n)+j] = 0.0;
               }
            }
            g[ct] = - l[i];
            ct++;
         }
      }

      assert( w >= ct );
   }

   assert( w == ct );

#if debug
   cout << "F:" << endl;
   //printM( w, n, F);
   cout << "g: ";
   //printv( w, g);
#endif

   int      h = n+w;
   double   *H = nullptr;  // [h*h]

   H = new double[h*h];

   for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
         H[(i*h)+j] = 2 * Q[(i*n)+j];
      }
   }

   if( w > 0 ){
      for(int i=n; i<h; i++){
         for(int j=0; j<n; j++){
            H[(i*h)+j] = F[((i-n)*n)+j];
            H[(j*h)+i] = F[((i-n)*n)+j];
         }
         for(int j=n; j<h; j++){
            H[(i*h)+j] = 0.0;
         }
      }
   }

   double   *z = nullptr;  // [h]
   double   *r = nullptr;  // [h]

   z = new double[h];
   r = new double[h];

   if( p == nullptr ){
      for(int i=0; i<n; i++){
         r[i] = 0.0;
      }
   }else{
      for(int i=0; i<n; i++){
         r[i] = - p[i];
      }
   }

   if( w > 0 ){
      for(int i=n; i<h; i++){
         assert( w > (i-n) );
         r[i] = g[(i-n)];
      }
   }

   #if debug
   {
      cout << "H:" << endl;
      printM( h, h, H);
      cout << "r: ";
      printv( h, r);
   }
   #endif

   int   com;
   switch ( SOLVE_PROBLEM_ACTIVE ) {
      case 0:
      {
         com = Com_LS_dgesv( H, r, h, z);
         break;
      }
      case 1:
      {
         com = Com_LS_dsysv( H, r, h, z);
         break;
      }
      default:
      {
         cout << "error: SOLVE_PROBLEM_ACTIVE is " << SOLVE_PROBLEM_ACTIVE << endl;
         exit(-1);
         break;
      }
   }

   #if debug
   {
      cout << "com: " << com << endl;
      cout << "z: ";
      printv( h, z);
   }
   #endif

   if( com != 0 ){
      cout << "error: com = " << com << endl;
      exit(-1);
   }


   Copy_vec( z, x_new, n);

   double min = 0.0;
   int memo = -1;
   double ep = get_ep();

   for(int i=n+k; i<h; i++){
      if( z[i] < min && fabs( z[i] ) > ep ){
         min = z[i];
         memo = i - ( n + k );
      }
   }

   delete[] F;
   delete[] g;
   delete[] H;
   delete[] z;
   delete[] r;

   return memo;
}

int QPsolver::solve_activeset_onlybounds(
   double   *x_new,
   bool     *W_u,
   bool     *W_l
   )
{
   /**********************************************
      min   x'Qx + p'x
      s.t.  x_i = a_i (i \in active set)
            x in R^n

      -> W := { i in {1,...,n} : x _i = a_i }
         X := {1,...,n}\W

         c := sum_{i \in W} a_i b_i

         x'Qx = ||Bx||^2 = ||subB subx + c||^2
         = ||subB subx||^2 + 2 c'subB subx + c'c
         := z'Hz + 2 g'z + c'c

         where H = subB'subB, g = subB'c, z = subx

                               | b_1 'b_i |
         g = sum_{i \in W} a_i |  ....    |
                               | b_#X'b_i |

         Hence, we solve

         min   z'Hz + (2 g + subp)'z + c'c + sum_{i in W} a_i p_i
         s.t.  z in R^X

         z*: 2H z* = - ( 2 g + subp )

         Note: r = - ( 2 g + subp )

         Note: for all i in W,
            lambda = - p_i - 2 (Q x*)_i

   **********************************************/
   assert( x_new != nullptr );
   assert( k == 0 );
   assert( W_u != nullptr );
   assert( W_l != nullptr );

   int   wu = count( n, W_u );
   int   wl = count( n, W_l );
   int   w = wu + wl;
   int   nvars = n - w;

   int buf_int = 0;
   int buf_bool;
   double coef = 0.0;
   int ct1;
   int ct2;

   assert( wu >= 0 );
   assert( wl >= 0 );
   assert( wl >= 0 );
   assert( nvars >= 1 );

   #if debug
   {
      printf("(w, wu, wl, n_vars) = (%d, %d, %d, %d)\n", w, wu, wl,nvars);
   }
   #endif

   assert( H != nullptr );
   assert( r != nullptr );
   assert( indexes != nullptr );
   assert( z != nullptr );
   assert( fixed != nullptr );

   double   lambda = 0.0;
   double   min = 0.0;
   int      memo = -1;
   double   ep = get_ep();
   if ( nvars == 0 )
   {
      for ( int i = 0; i < n; i++ )
      {
         assert( W_l[i] == true || W_u[i] == true );
         assert( W_l[i] == false || W_u[i] == false );

         if ( W_l[i] == true )
            x_new[i] = l[i];
         else
            x_new[i] = u[i];
      }

      ct1 = 0; // lower
      ct2 = 0; // upper
      for ( int i = 0; i < n; i++ )
      {
         lambda = 0.0;
         if ( p != nullptr )
            lambda -= p[i];

         //cout << "lambda: " << lambda << endl;

         buf_int = i * n;
         for ( int t = 0; t < n; t++ )
            lambda -= 2 * Q[buf_int + t] * x_new[t];

         assert( W_l[i] == true || W_u[i] == true );
         assert( W_l[i] == false || W_u[i] == false );

         buf_bool = W_l[i];
         if ( buf_bool == true )
         {
            lambda *= - 1.0;
            buf_int = ct1 + wu;
            ct1++;
         }
         else
         {
            buf_int = ct2;
            ct2++;
         }
      }

      //cout << "lambda: " << lambda << endl;

      if ( lambda < min && fabs( lambda ) > ep )
      {
         min = lambda;
         memo = buf_int;
      }

      return memo;

   }

   ct1 = 0;
   ct2 = n - 1;

   for ( int i = 0; i < n; i++ )
   {
      if ( W_u[i] == true || W_l[i] == true )
         indexes[ct2--] = i;
      else
         indexes[ct1++] = i;
   }

   #if debug
   {
      cout << "indexes: ";
      printv( n, indexes );
   }
   #endif

   assert( ct1 == ct2 + 1 );

   ct1 = 0;

   // set H and initialize r
   for ( int i = 0; i < nvars; i++ )
   {
      buf_int = indexes[i]*n;
      for ( int j = 0; j < nvars; j++ )
      {
         H[ct1++] = 2.0 * Q[buf_int + indexes[j]];
      }

      r[i] = 0.0;
   }

   // set r
   ct1 = 0;
   for ( int i = n - 1, j; i >= nvars; i-- )
   {
      j = indexes[i];

      assert( W_l[j] == true || W_u[j] == true );
      assert( W_l[j] == false || W_u[j] == false );

      if ( W_l[j] == true )
      {
         coef = - 2.0 * l[j];
         fixed[ct1++] = l[j];
      }
      else
      {
         coef = - 2.0 * u[j];
         fixed[ct1++] = u[j];
      }

      buf_int = j*n;
      for ( int t = 0; t < nvars; t++ )
         r[t] += coef * Q[buf_int + indexes[t]];
   }

   if ( p != nullptr )
   {
      for ( int i = 0; i < nvars; i++ )
         r[i] -= p[indexes[i]];
   }

   #if debug
   {
      cout << "H:" << endl;
      printM( nvars, nvars, H);
      cout << "r: ";
      printv( nvars, r);
   }
   #endif

   int com;

   assert( comls != nullptr );
   com = comls( H, r, nvars, z);

   #if debug
   {
      cout << "com: " << com << endl;
      cout << "z: ";
      printv( nvars, z);
   }
   #endif

   if( com != 0 )
   {
      cout << "error: com = " << com << endl;
      printf("( n, nvars, w ) = ( %d, %d, %d )\n", n, nvars, w );
      exit(-1);
   }

   for ( int i = 0; i < nvars; i++ )
      x_new[indexes[i]] = z[i];

   for ( int i = n - 1, j = 0; i >= nvars; i--, j++ )
   {
      assert( j < w );
      x_new[indexes[i]] = fixed[j];
   }

   #if debug
   {
      cout << "x_new: ";
      printv( n, x_new);
   }
   #endif

   ct1 = 0; // lower
   ct2 = 0; // upper
   for ( int i = 0, j = n - 1, ind_j; i < w; i++, j-- )
   {
      lambda = 0.0;
      ind_j = indexes[j];
      if ( p != nullptr )
         lambda -= p[ind_j];

      //cout << "lambda: " << lambda << endl;
      //cout << "ind_j: " << ind_j << endl;

      buf_int = ind_j * n;
      for ( int t = 0; t < n; t++ )
      {
         lambda -= 2 * Q[buf_int + t] * x_new[t];
      }

      assert( W_l[ind_j] == true || W_u[ind_j] == true );
      assert( W_l[ind_j] == false || W_u[ind_j] == false );

      buf_bool = W_l[ind_j];
      if ( buf_bool == true )
      {
         lambda *= - 1.0;
         buf_int = ct1 + wu;
         ct1++;
      }
      else
      {
         buf_int = ct2;
         ct2++;
      }


      //cout << "lambda: " << lambda << endl;

      if ( lambda < min && fabs( lambda ) > ep )
      {
         min = lambda;
         memo = buf_int;
      }
   }

   assert( ct1 == wl );
   assert( ct2 == wu );

   return memo;
}

double   QPsolver::compute_stepsize(
   double   *x,
   double   *dx,
   bool     *W_ineq,
   bool     *W_u,
   bool     *W_l
   )
{
   assert( x != nullptr );
   assert( dx != nullptr );
   assert( Com_nrm( dx, n) >= get_ep() );

   double   alpha = 1.0;
   double ep = get_ep();

   if( W_ineq != nullptr ){
      assert( m > 0 );
      assert( A != nullptr );
      assert( b != nullptr );
      for(int i=0; i<m; i++){
         if( W_ineq[i] == false ){
            double   dot = Com_dot( &A[i*n], dx, n);
            if( dot > ep ){
               double   a = ( b[i] - Com_dot( &A[i*n], x, n) ) / dot;
               if( alpha > a ){
                  alpha = a;
               }
            }
         }
      }
   }


   if( W_u != nullptr ){
      assert( n > 0 );
      assert( u != nullptr );
      for(int i=0; i<n; i++){
         if( W_u[i] == false && dx[i] > ep ){
            double   a = ( u[i] - x[i] ) / dx[i];
            if( alpha > a ){
               alpha = a;
            }
         }
      }
   }

   if( W_l != nullptr ){
      assert( n > 0 );
      assert( l != nullptr );
      for(int i=0; i<n; i++){
         if( W_l[i] == false && dx[i] < -ep ){
            double   a = ( l[i] - x[i] ) / dx[i];
            if( alpha > a ){
               alpha = a;
            }
         }
      }
   }

   if( fabs( alpha ) < ep ){
      alpha = 0.0;
   }

   return alpha;
}

bool QPsolver::update_activeset(
   double   *x_new,
   bool     *W_ineq,
   bool     *W_u,
   bool     *W_l
   )
{
   double ep = get_ep();
   bool  update = false;

   if( W_ineq != nullptr ){
      assert( m > 0 );
      assert( A != nullptr );
      assert( b != nullptr );
      for(int i=0; i<m; i++){
         if( W_ineq[i] == false ){
            double buf = Com_dot( &A[i*n], x_new, n);
            if( Equal( b[i], buf, ep) == true ){
               W_ineq[i] = true;
               update = true;
            }
            //if( fabs( b[i] - Com_dot( &A[i*n], x_new, n)) < ep ){
            // W_ineq[i] = true;
            // update = true;
            //}
         }
      }
   }

   if( W_u != nullptr ){
      assert( n > 0 );
      assert( u != nullptr );
      for(int i=0; i<n; i++){
         if( W_u[i] == false ){
            if( Equal( u[i], x_new[i], ep) == true ){
            //if( fabs( u[i] - x_new[i] ) < ep ){
               W_u[i] = true;
               update = true;
               if( u[i] - x_new[i] != 0 ) x_new[i] = u[i];
            }
         }
      }
   }

   if( W_l != nullptr ){
      assert( n > 0 );
      assert( l != nullptr );
      for(int i=0; i<n; i++){
         if( W_l[i] == false ){
            if( Equal( x_new[i], l[i], ep) == true ){
            //if( fabs( x_new[i] - l[i] ) < ep ){
               W_l[i] = true;
               update = true;
               if( x_new[i] - l[i] != 0 ) x_new[i] = l[i];
            }
         }
      }
   }

   return update;

}

double QPsolver::compute_objval(
   double   *x
   )
{
   assert( n > 0 );
   assert( Q != nullptr );
   assert( x != nullptr );

   double   obj = 0.0;

   double   *Qx; //[n]
   Qx = new double[n];

   Gen_ZeroVec( n, Qx);
   Com_mat_Ax( Q, n, n, x, Qx);

   obj += Com_dot( x, Qx, n);

   if( p != nullptr ){
      obj += Com_dot( x, p, n);
   }

   delete[] Qx;

   return obj;
}
