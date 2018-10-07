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
#include <utility>

#include "svpsolver.h"
#include "node.h"
#include "probdata.h"
#include "qpsolver.h"
#include "solution.h"
#include "vector.h"

using namespace std;

#define debug  0

static int find_branchingvariable(
      const int      m,
      const double*  relax_solval
   )
{
   double buf;
   for ( int i = m - 1; 0 <= i; --i )
   {
      buf = relax_solval[i];
      if( !Equal( buf, round( buf ), 1.0e-12 ) )
         return i;
   }

   return -1;
}

void  SVPsolver::SVPSbranch(
      NODE&    node,
      int      index
   )
{
   if( node.get_zero() )
      SVPSbranch_BIN( node, index );
   else
      SVPSbranch_INT( node, index );
}

void  SVPsolver::SVPSbranch_BIN(
      NODE&    node,
      int      index
   )
{
   NODE  C_BIN;

   int      m = probdata.get_m();
   double*  set_ub;
   double*  set_lb;
   double*  set_warm;
   double   set_relax_objval;
   int      set_dpt;
   bool     set_zero;
   int      set_index;

   auto  node_ub = node.get_ub();
   auto  node_lb = node.get_lb();
   auto  warm = node.get_warm();

   set_ub = new double[m];
   set_lb = new double[m];
   set_warm = new double[m];

   auto ep = epsilon;

   // x_i = 0
   int memo = -1;
   Copy_vec( node.get_ub(), set_ub, m );
   Copy_vec( node.get_lb(), set_lb, m );

   for ( int i = 0; i < m; ++i )
   {
      if ( set_lb[i] != 0.0 || set_ub[i] != 0.0 )
      {
         set_lb[i] = 0.0;
         set_ub[i] = 0.0;
         memo = i;
         break;
      }
   }

   for( int i = 0, w_i, l_i, u_i; i < m; ++i )
   {
      w_i = warm[i];
      l_i = set_lb[i];
      u_i = set_ub[i];

      if( w_i < l_i || w_i > u_i )
         set_warm[i] = ( u_i + l_i ) / 2.0;
      else
         set_warm[i] = w_i;
   }

   set_relax_objval = node.get_lowerbound();

   set_dpt = node.get_dpt();
   ++set_dpt;

   set_zero = true;
   set_index = index;

   tighten_bounds( memo, set_lb, set_ub );

   C_BIN.set_vals( m, set_ub, set_lb,
                   set_warm, set_relax_objval,
                   set_dpt,  set_zero, set_index);

   if( node.get_sumfixed() != nullptr )
   {
      bool r;
      r = C_BIN.alloc_sumfixed();
      assert( r == true );
      C_BIN.set_sumfixed( 1.0, node.get_sumfixed() );
   }

   assert( set_lb[memo] == 0 && set_ub[memo] == 0 );

   for ( int i = 0; i < m; ++i )
   {
      if ( i !=memo && Equal( set_lb[i], set_ub[i], ep ) )
      {
         double fixedvalue = set_lb[i];
         if ( !Equal( fixedvalue, node_lb[i], ep ) && !Equal( fixedvalue, 0.0, ep ) )
         {
            if( C_BIN.alloc_sumfixed() == true )
               C_BIN.set_sumfixed( fixedvalue, probdata.get_bvec(memo) );
            else
               C_BIN.add_sumfixed( fixedvalue, probdata.get_bvec(memo) );
         }
      }
   }

   (nodelist.*move_back)( C_BIN );
   // Do not use C_BIN after here

   assert( memo >= 0 );

   // x_i >= 1
   warm = node.get_warm();
   node_ub = node.get_ub();
   node_lb = node.get_lb();

   node_lb[memo] = 1.0;

   for ( int i = 0, l_i, u_i, w_i; i < m; ++i )
   {
      w_i = warm[i];
      u_i = node_ub[i];
      l_i = node_lb[i];

      if( w_i < l_i || w_i > u_i )
         warm[i] = l_i;
   }

   assert( node.get_zero() == true );

   node.set_dpt( node.get_dpt() + 1 );
   node.set_zero( false );
   node.set_index( index + 1 );

   if( 1.0 == node_ub[memo] )
   {
      // then node_lb[memo] == 1.0
      if( node.alloc_sumfixed() == true )
         node.set_sumfixed( 1.0, probdata.get_bvec(memo) );
      else
         node.add_sumfixed( 1.0, probdata.get_bvec(memo) );
   }

   delete[] set_ub;
   delete[] set_lb;
   delete[] set_warm;

}
void  SVPsolver::SVPSbranch_INT(
      NODE&    node,
      int      index
   )
{
   NODE  C_RIGHT;

   int      m = probdata.get_m();
   double   *set_ub;
   double   *set_lb;
   double   *set_warm;
   double   set_relax_objval;
   int      set_dpt;
   bool     set_zero;
   int      set_index;

   double   *relax_solval = node.get_relaxsolval();

   set_ub = new double[m];
   set_lb = new double[m];
   set_warm = new double[m];

   assert( relax_solval != nullptr );

   auto     ep = epsilon;

   // find a branching variable
   int memo = find_branchingvariable( m, relax_solval );

   assert( memo >= 0 );
   assert( memo < m );

   // ceil <= x_memo
   Copy_vec( node.get_ub(), set_ub, m);
   Copy_vec( node.get_lb(), set_lb, m);

   double ceil_value = ceil( relax_solval[memo] );
   set_lb[memo] = ceil_value;

   Copy_vec( relax_solval, set_warm, m );

   set_warm[memo] = ceil_value;

   set_relax_objval = node.get_lowerbound();

   set_dpt = node.get_dpt();
   ++set_dpt;

   set_zero = false;
   set_index = index;

   C_RIGHT.set_vals( m, set_ub, set_lb,
                   set_warm, set_relax_objval,
                   set_dpt,  set_zero, set_index );

   if( node.get_sumfixed() != nullptr )
   {
      bool r;
      r = C_RIGHT.alloc_sumfixed();
      assert( r == true );
      C_RIGHT.set_sumfixed( 1.0, node.get_sumfixed() );
   }

   if( ceil_value == set_ub[memo] && ceil_value != 0.0 )
   {
      if( C_RIGHT.alloc_sumfixed() == true )
         C_RIGHT.set_sumfixed( ceil_value, probdata.get_bvec(memo) );
      else
         C_RIGHT.add_sumfixed( ceil_value, probdata.get_bvec(memo) );
   }

   (nodelist.*move_back)( C_RIGHT );
   // Do not use C_RIGHT after here

   // x_memo <= floor
   auto  node_ub = node.get_ub();
   auto  node_lb = node.get_lb();
   auto  warm = node.get_warm();

   double floor_value = floor( relax_solval[memo] );
   node_ub[memo] = floor_value;

   Copy_vec( relax_solval, warm, m );
   warm[memo] = floor_value;

   node.set_dpt( node.get_dpt() + 1 );
   node.set_zero( false );
   node.set_index( index + 1 );

   if ( Equal( node_lb[memo], floor_value, ep ) && !Equal( floor_value, 0.0, ep ) )
   {
      if( node.alloc_sumfixed() == true )
         node.set_sumfixed( floor_value, probdata.get_bvec(memo) );
      else
         node.add_sumfixed( floor_value, probdata.get_bvec(memo) );
   }

   delete[] set_ub;
   delete[] set_lb;
   delete[] set_warm;
}


