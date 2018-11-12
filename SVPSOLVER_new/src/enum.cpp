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

//#include <omp.h>

#include "svpsolver.h"
#include "vector.h"
#include "node.h"
#include "probdata.h"

#define debug 0

struct EnumNode
{
   int m;
   int ub;
   int lb;
   double lower;
   double* SFV; // sum of fixed vectors

   EnumNode( const int s_ub, const int s_lb,
            const double s_lower, const double* s_SFV, const int s_m ) {
      m = s_m; ub = s_ub; lb = s_lb; lower = s_lower;
      assert( m > 0 && s_SFV != nullptr );
      SFV = new double[m];
      for ( int i = 0; i < m; ++i ) SFV[i] = s_SFV[i];
   }
   EnumNode( const EnumNode &source ) {
      m = source.m; ub = source.ub; lb = source.lb; lower = source.lower;
      assert( m > 0 && source.SFV != nullptr );
      SFV = new double[m];
      for ( int i = 0; i < m; ++i ) SFV[i] = source.SFV[i];
   }
   EnumNode& operator=( const EnumNode& source ) {
      if ( this != &source ) {
         m = source.m; ub = source.ub; lb = source.lb; lower = source.lower;
         assert( m > 0 && source.SFV != nullptr );
         SFV = new double[m];
         for ( int i = 0; i < m; ++i ) SFV[i] = source.SFV[i];
      }
      return *this;
   }
   ~EnumNode() { delete[] SFV; SFV = nullptr; }
};

//static auto checkEnum(
//      NODE&          node,
//      const double*  vars_localub,
//      const double*  vars_locallb
//      )
//{
//   assert( vars_localub != nullptr );
//   assert( vars_locallb != nullptr );
//   assert( !node.get_branchinfo().empty() );
//
//   const auto& lastbranchinfo = node.get_branchinfo().back();
//   const auto index = lastbranchinfo.get_index();
//
//   assert( node.get_dimention() > index && index >= 0 );
//   assert( vars_locallb[index] == lastbranchinfo.get_value()
//         || vars_localub[index] == lastbranchinfo.get_value() );
//
//   if ( vars_localub[index] != vars_locallb[index] )
//      return false;
//
//   const auto m = node.get_dimention();
//   int ct;
//   const int maxval = max( m - 45, 1 );
//   bool nonzero = false;
//
//   assert( m > 0 );
//
//   ct = 0;
//   for ( auto i = 0; i < m; ++i )
//   {
//      if ( vars_localub[i] == vars_locallb[i] )
//      {
//         ++ct;
//         if( (int) vars_localub[i] ) 
//            nonzero = true;
//         if ( ct >= maxval && nonzero )
//            return true;
//      }
//   }
//
//   return false;
//}
//
//static void setIndexes(
//      int            m,
//      vector<int>&   unfixedindexes,
//      vector<int>&   nonzero_fixedindexes,
//      const double*  vars_localub,
//      const double*  vars_locallb
//      )
//{
//   assert( vars_localub != nullptr );
//   assert( vars_locallb != nullptr );
//   assert( m > 0 );
//
//   for ( auto i = 0; i < m; ++i )
//   {
//      if ( vars_localub[i] != vars_locallb[i] )
//         unfixedindexes.push_back( i );
//      else if ( (int) vars_localub[i] )
//         nonzero_fixedindexes.push_back( i );
//
//      #if debug
//      printf(" %d <= x%d <= %d\n", (int)vars_locallb[i], i, (int)vars_localub[i] );
//      #endif
//   }
//}
//
//static void setGS(
//      // input
//      const int            m,
//      const PROB_DATA&     pd,
//      const vector<int>&   unfixedindexes,
//      const int            subdim,
//      const vector<int>&   nonzero_fixedindexes,
//      const int            numfixed,
//      // output
//      vector<double>&            SNOVs,
//      vector<vector<double>>&    CMGSO
//      )
//{
//   assert( !unfixedindexes.empty() );
//   assert( m > 0 );
//   assert( subdim > 0 && subdim <= m );
//   assert( subdim + numfixed == (int) SNOVs.size() );
//   assert( subdim + numfixed == (int) CMGSO.size() );
//   assert( subdim == (int) unfixedindexes.size() );
//   assert( nonzero_fixedindexes.empty() || numfixed == (int) nonzero_fixedindexes.size() );
//
//   const auto B_ = pd.get_B_();
//   const auto sdnf = subdim + numfixed;
//   const auto msdnf = m * sdnf;
//   const auto sdnf_1 = sdnf - 1;
//
//   double* OM = new double[msdnf];   // orthogonal matrix
//   vector<vector<double>> CF( sdnf );          // coefficient matrix
//
//   double coef;
//
//   // set SNOVs
//   for( int i = 0, im = 0, t, tm; i < subdim; ++i, im += m )
//   {
//      assert( unfixedindexes[i] < m && unfixedindexes[i] >= 0 );
//      assert( i * m == im );
//
//      t = unfixedindexes[i];
//      tm = t * m;
//
//      Copy_vec( &B_[tm], &OM[im], m );
//
//      CF[i].resize( i + 1 );
//      for ( int j = 0, jm = 0; j < i; ++j, jm += m )
//      {
//         assert( j * m == jm );
//
//         coef = Com_dot( &B_[tm], &OM[jm], m ) / SNOVs[j];
//         Com_linecomb( &OM[im], &OM[jm], m, 1.0, - coef, &OM[im] );
//         CF[i][j] = coef;
//      }
//      CF[i][i] = 1;
//
//      SNOVs[i] = Com_dot( &OM[im], &OM[im], m );
//   }
//   for( int i = subdim, im = subdim * m, t, tm; i < sdnf; ++i, im += m )
//   {
//      assert( nonzero_fixedindexes[i-subdim] < m && nonzero_fixedindexes[i-subdim] >= 0 );
//      assert( i * m == im );
//
//      t = nonzero_fixedindexes[i-subdim];
//      tm = t * m;
//
//      Copy_vec( &B_[tm], &OM[im], m );
//
//      CF[i].resize( i + 1 );
//      for ( int j = 0, jm = 0; j < i; ++j, jm += m )
//      {
//         assert( j * m == jm );
//
//         coef = Com_dot( &B_[tm], &OM[jm], m ) / SNOVs[j];
//         Com_linecomb( &OM[im], &OM[jm], m, 1.0, - coef, &OM[im] );
//         CF[i][j] = coef;
//      }
//      CF[i][i] = 1;
//
//      SNOVs[i] = Com_dot( &OM[im], &OM[im], m );
//
//   }
//
//   delete[] OM;
//   OM = nullptr;
//
//   //for ( auto cfs : CF )
//   //{
//   //   for ( auto cf : cfs )
//   //      cout << cf << " ";
//   //   cout << endl;
//   //}
//
//   // set CMGSO
//   for ( auto i = 0; i < sdnf; ++i )
//   {
//      CMGSO[i].resize( i + 1 );
//      for ( auto j = 0; j <= i; ++j )
//      {
//         assert( !CF[sdnf_1 - j].empty() );
//         assert( sdnf_1 - i < (int) CF[sdnf_1 - j].size() );
//
//         CMGSO[i][j] = CF[sdnf_1 - j][sdnf_1 - i]; // sdnf_1 = sdnf - 1
//      }
//   }
//
//   //for ( auto cfs : CMGSO )
//   //{
//   //   for ( auto cf : cfs )
//   //      cout << cf << " ";
//   //   cout << endl;
//   //}
//
//}
//
//static void setVsTB1(
//      // input
//      const int            m,
//      const PROB_DATA&     pd,
//      const vector<int>&   unfixedindexes,
//      const int            subdim,
//      // output
//      vector<vector<double>>&    VsTB
//      )
//{
//   assert( !unfixedindexes.empty() );
//   assert( m > 0 );
//   assert( subdim > 0 && subdim <= m );
//   assert( subdim == (int) VsTB.size() );
//   assert( subdim == (int) unfixedindexes.size() );
//
//   const auto B_ = pd.get_B_();
//   const auto Q = pd.get_Q();
//   const auto msubdim = m * subdim;
//   const auto subdimdim = subdim * subdim;
//   const auto subdim_1 = subdim - 1;
//
//   double*  subB = new double[msubdim];
//   double*  subBB = new double[subdimdim];
//   double*  e = new double[subdim];
//   double*  buf_b = new double[m];
//
//   int ct;
//
//   // set subB
//   for( int i = 0, im = 0, t, tm; i < subdim; ++i, im += m )
//   {
//      assert( unfixedindexes[i] < m && unfixedindexes[i] >= 0 );
//      assert( i * m == im );
//
//      t = unfixedindexes[i];
//      tm = t * m;
//
//      Copy_vec( &B_[tm], &subB[im], m );
//   }
//
//   // set subBB
//   ct = 0;
//   for ( int i = 0, t, tm; i < subdim; ++i )
//   {
//      t = unfixedindexes[i];
//      tm = t * m;
//      for ( int j = 0, s; j < subdim; ++j )
//      {
//         s = unfixedindexes[j];
//         subBB[ct] = Q[tm + s];
//         assert( subBB[ct] == Com_dot( &subB[i*m], &subB[j*m], m ) );
//         ++ct;
//      }
//   }
//   assert( ct == subdimdim );
//
//   Gen_ZeroVec( subdim, e );
//   e[subdim_1] = 1.0;
//
//   if ( Com_LS_dposv_nocopy( subBB, e, subdim ) != 0 )
//   {
//      cout << "error: enum" <<endl;
//      assert(0);
//      exit(-1);
//   }
//
//   Com_mat_Ax( subB, m, subdim, e, buf_b );
//
//   VsTB[0].resize( m );
//   for ( auto i = 0; i < m; ++i )
//      VsTB[0][i] = buf_b[i];
//
//   int dim = subdim_1;
//   double* LTM = subBB;  // lower triangular matrix
//   double* buf_a = e;
//   int dim_1;
//   double L_dimdim;
//   double* p;
//   ct = 1;
//   while ( dim >= 1 )
//   {
//      dim_1 = dim - 1;
//      L_dimdim = LTM[subdim * dim_1 + dim_1];
//
//      // buf_a[dim-1] = 1.0/(L_kk * L_kk)
//      p = &buf_a[dim_1];
//      *p = 1.0;
//      *p /= L_dimdim;
//      *p /= L_dimdim;
//
//      for ( int t = dim_1 - 1, tdim; t >= 0; --t )
//      {
//         tdim = t * subdim;
//         p = &buf_a[t];
//         *p = 0.0;
//         for ( auto j = t + 1; j < dim; ++j )
//         {
//            *p += LTM[tdim + j] * buf_a[j];
//         }
//         *p *= - 1.0;
//         *p /= LTM[tdim + t];
//      }
//
//      Com_mat_Ax( subB, m, dim, buf_a, buf_b );
//
//      VsTB[ct].resize( m );
//      for ( auto i = 0; i < m; ++i )
//         VsTB[ct][i] = buf_b[i];
//
//      ++ct;
//      --dim;
//   }
//
//   assert( ct == subdim );
//
//   delete[] subB;
//   delete[] subBB;
//   delete[] e;
//   delete[] buf_b;
//}
//
////static void setVsTB2(
////      // input
////      const int            m,
////      const PROB_DATA&     pd,
////      const vector<int>&   unfixedindexes,
////      const int            subdim,
////      // output
////      vector<vector<double>>&    VsTB
////      )
////{
////   assert( !unfixedindexes.empty() );
////   assert( m > 0 );
////   assert( subdim > 0 && subdim < m );
////   assert( subdim == (int) VsTB.size() );
////   assert( subdim == (int) unfixedindexes.size() );
////
////   const auto B_ = pd.get_B_();
////   const auto Q = pd.get_Q();
////
////   auto dim = subdim;
////   auto mdim = m * dim;
////   auto dimdim = dim * dim;
////   auto dim_1 = dim - 1;
////
////   double*  subB = new double[mdim];
////   double*  subBB = new double[dimdim];
////   double*  e = new double[dim];
////   double*  buf_a = new double[dim];
////   double*  buf_b = new double[m];
////
////   int ct;
////   int k = 0;
////
////   while ( dim >= 1 )
////   {
////      mdim = m * dim;
////      dimdim = dim * dim;
////      dim_1 = dim - 1;
////
////      // set subB
////      for( int i = 0, im = 0, t, tm; i < dim; ++i, im += m )
////      {
////         assert( unfixedindexes[i] < m && unfixedindexes[i] >= 0 );
////         assert( i * m == im );
////
////         t = unfixedindexes[i];
////         tm = t * m;
////
////         Copy_vec( &B_[tm], &subB[im], m );
////      }
////
////      // set subBB
////      ct = 0;
////      for ( int i = 0, t, tm; i < dim; ++i )
////      {
////         t = unfixedindexes[i];
////         tm = t * m;
////         for ( int j = 0, s; j < dim; ++j )
////         {
////            s = unfixedindexes[j];
////            subBB[ct] = Q[tm + s];
////            assert( subBB[ct] == Com_dot( &subB[i*m], &subB[j*m], m ) );
////            ++ct;
////         }
////      }
////      assert( ct == dimdim );
////
////      Gen_ZeroVec( dim, e );
////      e[dim_1] = 1.0;
////
////      if ( Com_LS_dposv( subBB, e, dim, buf_a ) != 0 )
////      {
////         cout << "error: enum" <<endl;
////         assert(0);
////         exit(-1);
////      }
////
////      Com_mat_Ax( subB, m, dim, buf_a, buf_b );
////
////      VsTB[k].resize( m );
////      for ( auto i = 0; i < m; ++i )
////         VsTB[k][i] = buf_b[i];
////
////      dim--;
////      k++;
////   }
////
////   assert( k == subdim );
////
////   delete[] subB;
////   delete[] subBB;
////   delete[] e;
////   delete[] buf_a;
////   delete[] buf_b;
////}
//static void computeBounds(
//      const int               m,
//      const double            sqrt_bestval,
//      const double            norm,
//      const vector<double>    vec,
//      const double*           sumfixed,
//      const double            localub,
//      const double            locallb,
//      int*                    ub,
//      int*                    lb
//      )
//{
//   assert( m > 0 );
//   assert( sqrt_bestval > 0 );
//   assert( norm > 0 );
//   assert( !vec.empty() );
//   assert( m == (int) vec.size() );
//   assert( sumfixed != nullptr );
//   assert( locallb < localub );
//   assert( ub != nullptr );
//   assert( lb != nullptr );
//
//   double compute_ub;
//   double compute_lb;
//   double vec_sumfixed;
//
//   compute_ub = sqrt_bestval;
//   compute_ub *= norm;
//
//   compute_lb = - compute_ub;
//
//   vec_sumfixed = Com_dot( &vec[0], sumfixed, m );
//
//   compute_ub -= vec_sumfixed;
//   compute_lb -= vec_sumfixed;
//
//   compute_ub = floor( compute_ub );
//   compute_lb = ceil( compute_lb );
//
//   //*ub = localub;
//   //*lb = locallb;
//   if ( compute_ub >= localub )
//      *ub = localub;
//   else
//   {
//      #if debug
//      printf("tightened ub! %d -> %d\n", (int)localub, (int)compute_ub );
//      #endif
//      *ub = compute_ub;
//   }
//
//   if ( locallb >= compute_lb )
//      *lb = locallb;
//   else
//   {
//      #if debug
//      printf("tightened lb! %d -> %d\n", (int)locallb, (int)compute_lb );
//      #endif
//      *lb = compute_lb;
//   }
//}

RelaxResult SVPsolver::SVPSenumerate(
      NODE&          node,
      const double*  vars_localub,
      const double*  vars_locallb
      )
{
   assert(0);
   return R_FEASIBLE;

   //assert( vars_localub != nullptr );
   //assert( vars_locallb != nullptr );

   //if ( !checkEnum( node, vars_localub, vars_locallb ) )
   //   return R_FEASIBLE;

   //const auto& pd = probdata;
   //const auto m = pd.get_m();
   //const auto B_ = pd.get_B_();

   //assert( m > 0 );

   //vector<int> unfixedindexes;
   //vector<int> nonzero_fixedindexes;

   //// set indexes
   //setIndexes( m, unfixedindexes, nonzero_fixedindexes, vars_localub, vars_locallb );

   //assert( !unfixedindexes.empty() );

   //const int subdim = (int) unfixedindexes.size();
   //const int numfixed = (int) nonzero_fixedindexes.size();
   //const int sdnf = subdim + numfixed;
   //vector<double> SNOVs( subdim + numfixed );           // square norm of orthogonal vectors
   //vector<vector<double>> CMGSO( subdim + numfixed );   // coefficient matrix for Gram-Schmidt orthogonalization
   //vector<vector<double>> VsTB( subdim );               // vectors for tightening bounds
   //vector<double> NVsTB;                                // norm of vectors for tightening bounds

   //assert( subdim > 0 && subdim <= m );

   //// set SNOVs and GSOM;
   //setGS( m, pd, unfixedindexes, subdim, nonzero_fixedindexes, numfixed, SNOVs, CMGSO );
   ////set VsTB
   //setVsTB1( m, pd, unfixedindexes, subdim, VsTB );

   //NVsTB.reserve( subdim );
   //for ( auto& vec : VsTB )
   //   NVsTB.push_back( Com_nrm( &vec[0], m ) );

   //#if debug
   //cout << "unfixedindexes: ";
   //for ( auto i : unfixedindexes )
   //   cout << i << " ";
   //cout << endl;
   //cout << "nonzero_fixedindexes: ";
   //for ( auto i : nonzero_fixedindexes )
   //   cout << i << " ";
   //cout << endl;
   //cout << "SNOVs: ";
   //for ( auto s : SNOVs )
   //   cout << s << " ";
   //cout << endl;
   //#endif

   //int dpt;
   //int dpt_nf; // dpt + numfixed
   //int maxdpt;
   //int index;
   //int origindex;
   //int ub;
   //int lb;
   //double lower;
   //double lower_new = 0.0;
   //int upper;
   //int upper_new;
   //double sqrt_bestval;

   //double* sumfixed;
   //int* solution;

   //vector<EnumNode> enumtree;
   //double y;
   //bool update = false;

   //dpt = 0;
   //dpt_nf = numfixed;
   //maxdpt = subdim - 1;
   //index = subdim - 1;
   //origindex = unfixedindexes[index];

   //solution = new int[sdnf];
   //sumfixed = new double[m];
   ////assert( node.get_sumfixed() != nullptr );
   ////Copy_vec( node.get_sumfixed(), sumfixed, m );
   //if ( node.get_sumfixed() != nullptr )
   //   Copy_vec( node.get_sumfixed(), sumfixed, m );
   //else
   //   Gen_ZeroVec( m, sumfixed );

   //reverse( SNOVs.begin(), SNOVs.end() );

   //for ( int i = 0, t, buf = numfixed - 1; i < numfixed; ++i )
   //{
   //   t = nonzero_fixedindexes[buf - i];
   //   solution[i] = vars_localub[t];

   //   assert( vars_localub[t] == vars_locallb[t] );
   //}

   //lower = 0.0;
   //for ( int i = 0; i < numfixed; ++i )
   //{
   //   y = 0.0;
   //   for ( int j = 0; j <= i; ++j )
   //   {
   //      y += solution[j] * CMGSO[i][j];
   //   }

   //   lower += y * y * SNOVs[i];
   //}

   //upper = bestval;
   //sqrt_bestval = sqrt( bestval );

   //if ( (int) lower >= upper )
   //{
   //   cout << "error: enum" << endl;
   //   assert(0);
   //   exit(-1);
   //}

   //computeBounds( m, sqrt_bestval, NVsTB[0], VsTB[0], sumfixed, vars_localub[origindex], vars_locallb[origindex],
   //      &ub, &lb );

   //if ( lb > ub )
   //{
   //   //cout << "enum:INFEASIBLE" << endl;
   //   delete[] sumfixed;
   //   delete[] solution;
   //   return R_INFEASIBLE;
   //}

   //solution[dpt_nf] = lb;
   //enumtree.reserve( subdim );
   //enumtree.emplace_back( ub, lb, lower, sumfixed, m );

   //#if debug
   //cout << "get_lowerbound: " << node.get_lowerbound() << endl;
   //cout << "init_lower: " << lower << endl;
   //cout << "init_upper: " << upper << endl;
   //cout << "init_sol: ";
   //for ( auto i = 0; i < numfixed; ++i )
   //   printf("x%d = %d, ", nonzero_fixedindexes[numfixed-1-i], solution[i] );
   //cout << endl;
   //#endif

   //while ( dpt >= 0 )
   //{
   //   assert( index == subdim - 1 - dpt );
   //   assert( dpt + numfixed == dpt_nf );
   //   assert( !enumtree.empty() );
   //   assert( dpt + 1 == (int)enumtree.size() );
   //   if ( dpt < maxdpt )
   //   {
   //      while( solution[dpt_nf] <= ub )
   //      {
   //         y = 0.0;
   //         for ( int j = 0; j <= dpt_nf; ++j )
   //         {
   //            y += solution[j] * CMGSO[dpt_nf][j];
   //         }
   //         //lower_new = lower + y * y * SNOVs[dpt_nf];
   //         lower_new = SNOVs[dpt_nf];
   //         lower_new *= y;
   //         if ( (int)lower_new <= upper ) // !! lower_new may be too large
   //            lower_new *= y;
   //         if ( (int)lower_new <= upper ) // !! lower_new may be too large
   //            lower_new += lower;

   //         #if debug
   //         printf("x%d = %d [%d,%d] -> %d \n", origindex, solution[dpt_nf], lb, ub, (int)lower_new );
   //         #endif

   //         if ( (int) lower_new >= upper )
   //            ++solution[dpt_nf];
   //         else
   //            break;
   //      }
   //   }
   //   else
   //   {
   //      assert( dpt == maxdpt );
   //      while ( solution[dpt_nf] <= ub )
   //      {
   //         y = 0.0;
   //         for ( int j = 0; j <= dpt_nf; ++j )
   //         {
   //            y += solution[j] * CMGSO[dpt_nf][j];
   //         }
   //         upper_new = lower + y * y * SNOVs[dpt_nf];
   //         #if debug
   //         printf("x%d = %d [%d,%d] -> *%d \n", origindex, solution[dpt_nf], lb, ub, (int)upper_new );
   //         #endif

   //         if ( upper_new < upper )
   //         {
   //            upper = upper_new;
   //            double* solvals = new double[m];
   //            for ( int i = sdnf - 1, j = 0; j < subdim; --i, ++j )
   //               solvals[unfixedindexes[j]] = solution[i];
   //            for ( int i = 0; i < m; ++i )
   //            {
   //               if ( vars_localub[i] == vars_locallb[i] )
   //                  solvals[i] = vars_localub[i];
   //            }
   //            SOLUTION setsol;
   //            setsol.set_sol( m, solvals, upper );
   //            SVPStrySol( setsol, true, true, nullptr );
   //            delete[] solvals;
   //         }
   //         else
   //            ++solution[dpt_nf];
   //      }
   //   }

   //   if ( solution[dpt_nf] > ub )
   //   {
   //      assert( !enumtree.empty() );
   //      enumtree.pop_back();
   //      if ( enumtree.empty() )
   //         break;
   //      auto& enumnode = enumtree.back();
   //      --dpt;
   //      --dpt_nf;
   //      ++index;
   //      origindex = unfixedindexes[index];
   //      lower = enumnode.lower;
   //      ub = enumnode.ub;
   //      lb = enumnode.lb;
   //      Copy_vec( enumnode.SFV, sumfixed, m );
   //      ++solution[dpt_nf];

   //   }
   //   else
   //   {
   //      for ( int i = 0, om = origindex * m; i < m; ++i )
   //         sumfixed[i] += solution[dpt_nf] * B_[om + i];

   //      ++dpt;
   //      ++dpt_nf;
   //      --index;
   //      origindex = unfixedindexes[index];
   //      lower = lower_new;
   //      computeBounds( m, sqrt_bestval, NVsTB[dpt], VsTB[dpt], sumfixed, vars_localub[origindex], vars_locallb[origindex],
   //            &ub, &lb );
   //      assert( dpt_nf >= 0 && dpt_nf < m );
   //      solution[dpt_nf] = lb;

   //      enumtree.emplace_back( ub, lb, lower, sumfixed, m );

   //      #if debug
   //      if ( ub < lb )
   //         printf("x%d:[%d,%d]\n", origindex, lb, ub);
   //      #endif
   //   }
   //}

   //delete[] sumfixed;
   //delete[] solution;

   //if ( update )
   //   return R_GETINTEGER;
   //else
   //   return R_INFEASIBLE;
}
