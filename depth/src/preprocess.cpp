
#include <assert.h>
#include <numeric>
#include <math.h>
#include <iostream>

#include "vector.h"
#include "preprocess.h"
#include "solution.h"
#include "solution_pool.h"

#define debug 0

struct EnumNode
{
   int ub;
   int lb;
   double lower;
   double sum;

   EnumNode( const int s_ub, const int s_lb,
            const double s_lower, const double s_sum )
   { ub = s_ub; lb = s_lb; lower = s_lower; sum = s_sum; }
   EnumNode( const EnumNode &source )
   { ub = source.ub; lb = source.lb; lower = source.lower; sum = source.sum; }
   EnumNode& operator=( const EnumNode& source ) {
      if ( this != &source ) {
         ub = source.ub; lb = source.lb; lower = source.lower; sum = source.sum;
      }
      return *this;
   }
   ~EnumNode() {}
};

static void setGS(
      // input
      const int                     n,
      const int                     m,
      const vector<vector<double>>& matrix,
      // output
      vector<double>&               snovs,
      vector<vector<double>>&       cmgso
      )
{
   assert( !matrix.empty() );
   assert( !matrix[0].empty() );
   assert( (int) matrix.size() == m );
   assert( (int) matrix[0].size() == n );
   assert( snovs.empty() );
   assert( cmgso.empty() );

   vector<vector<double>> OM( m );
   vector<vector<double>> CF( m );  // coefficient matrix
   const auto m_1 = m - 1;
   double coef;

   snovs.resize( m );
   cmgso.resize( m );

   // set SNOVs
   for( int i = 0; i < m; ++i  )
   {
      assert( !matrix[i].empty() );
      assert( (int) matrix[i].size() == n );

      OM[i] = matrix[i];

      CF[i].resize( i + 1 );
      for ( int j = 0; j < i; ++j )
      {
         assert( !OM[j].empty() );
         assert( (int) OM[j].size() == n );

         //CF[i][j] = Com_dot( &matrix[i], &OM[j], n ) / snovs[j];
         CF[i][j] = inner_product( matrix[i].begin(), matrix[i].end(), OM[j].begin(), 0.0 );
         CF[i][j] /= snovs[j];
         coef = CF[i][j];
         //Com_linecomb( &OM[i], &OM[j], n, 1.0, - CF[i][j], &OM[i] );
         for ( auto k = 0; k < n; ++k )
            OM[i][k] -= coef * OM[j][k];
      }
      CF[i][i] = 1;

      //snovs[i] = Com_dot( &OM[i], &OM[i], n );
      snovs[i] = inner_product( OM[i].begin(), OM[i].end(), OM[i].begin(), 0.0 );
   }

   // set CMGSO
   for ( auto i = 0; i < m; ++i )
   {
      cmgso[i].resize( i + 1 );
      for ( auto j = 0; j <= i; ++j )
      {
         assert( !CF[m_1 - j].empty() );
         assert( m_1 - i < (int) CF[m_1 - j].size() );

         cmgso[i][j] = CF[m_1 - j][m_1 - i]; // m_1 = m - 1
      }
   }
}

void FindNShorterVec(
      const int                     numofsol,
      const int                     timelimit,
      const vector<vector<double>>& matrix,
      SOLUTION_POOL*                solpool
      )
{
   assert( numofsol >= 1 );
   assert( timelimit >= 1 );
   assert( !matrix.empty() );

   const auto m = (int) matrix.size();
   const auto n = (int) matrix[0].size();

   assert( m >= 1 );
   //assert( n == (int) matrix[1].size() );

   #if debug
   printf("n:%d m:%d \n", n, m);
   #endif

   vector<double> SNOVs;
   vector<vector<double>> CMGSO;
   vector<EnumNode> enumtree;

   setGS( n, m, matrix, SNOVs, CMGSO );
   solpool->alloc( numofsol );

   for ( auto i = 0; i < m; ++i )
   {
      SOLUTION sol;
      vector<double> vals( m, 0 );
      vals[i] = 1;
      sol.set_sol( m, &vals[0], inner_product( matrix[i].begin(), matrix[i].end(), matrix[i].begin(), 0.0 ) );
      solpool->add_solution( sol );
   }

   #if debug || 0
   solpool->disp();
   cout << endl;
   #endif

   int dpt = 0;
   int branchindex = m - 1;
   const int maxdpt = m - 1;
   int ub;
   int lb;
   double lower = 0.0;
   double lower_new = 0.0;
   double upper = solpool->get_biggestval();
   double upper_new;
   vector<int> solvals( m );
   double sum = 0.0;
   double buf;
   const double ep = 1.0e-10;

   assert( ep < upper );

   buf = upper;
   buf /= SNOVs[branchindex];
   buf = sqrt( buf );

   ub = floor( buf + ep );
   lb = 0;

   assert( ub >= lb );

   solvals[0] = lb;
   enumtree.reserve( m );
   enumtree.emplace_back( ub, lb, lower, sum );

   while ( dpt >= 0 )
   {
      assert( branchindex == m - 1 - dpt );
      assert( !enumtree.empty() );
      assert( dpt + 1 == (int)enumtree.size() );
      if ( dpt < maxdpt )
      {
         while( solvals[dpt] <= ub )
         {
            lower_new = sum;
            lower_new += (double) solvals[dpt];
            lower_new *= lower_new;
            lower_new *= SNOVs[branchindex];
            lower_new += lower;

            #if debug
            printf("x%d = %d [%d,%d] -> %.8f \n", branchindex, solvals[dpt], lb, ub, lower_new );
            #endif

            if ( lower_new >= upper )
               ++solvals[dpt];
            else
               break;
         }
      }
      else
      {
         assert( dpt == maxdpt );
         while ( solvals[dpt] <= ub )
         {
            upper_new = sum;
            upper_new += (double) solvals[dpt];
            upper_new *= upper_new;
            upper_new *= SNOVs[branchindex];
            upper_new += lower;
            #if debug
            printf("x%d = %d -> upper: %.8f \n", branchindex, solvals[dpt], upper_new);
            #endif
            if ( upper_new < upper && upper_new > ep )
            {
               for ( auto i = 0, n_ct = 0; i < m; ++i )
               {
                  if ( solvals[i] )
                     ++n_ct;

                  if ( n_ct >= 2 )
                  {
                     SOLUTION setsol;
                     double* vals = new double[m];
                     for ( auto j = 0; j < m; ++j )
                        vals[j] = solvals[m-1-j];
                     setsol.set_sol( m, vals, upper_new );
                     solpool->add_solution( setsol );
                     delete[] vals;
                     upper = solpool->get_biggestval();
                     break;
                  }
               }
            }
            ++solvals[dpt];
         }
      }

      if ( solvals[dpt] > ub )
      {
         assert( !enumtree.empty() );
         enumtree.pop_back();
         if ( enumtree.empty() )
            break;
         auto& enumnode = enumtree.back();
         --dpt;
         ++branchindex;
         lower = enumnode.lower;
         ub = enumnode.ub;
         lb = enumnode.lb;
         sum = enumnode.sum;
         ++solvals[dpt];
      }
      else
      {
         const auto& vec = CMGSO[dpt+1];
         assert( (int) vec.size() == dpt + 2 );
         sum = 0;
         for ( auto i = 0; i < dpt+1; ++i )
         {
            sum += vec[i] * solvals[i];
         }


         ++dpt;
         --branchindex;
         lower = lower_new;

         buf = upper - lower;
         buf /= SNOVs[branchindex];
         buf = sqrt( buf );

         ub = floor( buf - sum + ep );

         if ( lower < ep )
            lb = 0;
         else
            lb = ceil( - buf - sum - ep );
         solvals[dpt] = lb;

         enumtree.emplace_back( ub, lb, lower, sum );
      }
   }
   #if debug || 0
   solpool->disp();
   cout << endl;
   #endif
}

