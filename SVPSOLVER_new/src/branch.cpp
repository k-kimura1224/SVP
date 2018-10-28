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
      NODE&       node,
      const int   index,
      double*     vars_localub,
      double*     vars_locallb
   )
{
   assert( index >= 0 );
   assert( vars_localub != nullptr );
   assert( vars_locallb != nullptr );

   SVPSbranch_INT( node, index, vars_localub, vars_locallb );
}

void  SVPsolver::SVPSbranch_INT(
      NODE&       node,
      const int   index,
      double*     vars_localub,
      double*     vars_locallb
   )
{
   assert( index >= 0 );
   assert( vars_localub != nullptr );
   assert( vars_locallb != nullptr );

   // copy
   const auto m = probdata.get_m();
   auto& nodeindex = index;
   auto& pd = probdata;
   auto ep = epsilon;

   assert( m > 0 );
   assert( nodeindex >= 0 );

   // node
   NODE C_RIGHT = node;
   auto relax_solval_R = C_RIGHT.get_relaxsolval();
   const auto branchindex = find_branchingvariable( m, relax_solval_R );
   const int ceil_value = ceil( relax_solval_R[branchindex] );
   const auto set_dpt = node.get_dpt() + 1;

   assert( relax_solval_R != nullptr );
   assert( branchindex >= 0 );
   assert( branchindex < m );
   assert( set_dpt > 0 );

   // ceil <= x_i
   C_RIGHT.set_branchinfo( branchindex, ceil_value, 'l' );
   C_RIGHT.set_dpt( set_dpt );
   C_RIGHT.set_index( nodeindex );

   if( ceil_value && Equal( (double) ceil_value, vars_localub[branchindex], ep ) )
   {
      if( C_RIGHT.alloc_sumfixed() )
         C_RIGHT.set_sumfixed( ceil_value, pd.get_bvec( branchindex ) );
      else
         C_RIGHT.add_sumfixed( ceil_value, pd.get_bvec( branchindex ) );
   }

   relax_solval_R[branchindex] = ceil_value;

   (nodelist.*move_back)( C_RIGHT );
   // Do not use C_RIGHT after here

   // x_i <= floor
   auto relax_solval_L = node.get_relaxsolval();
   const int floor_value = floor( relax_solval_L[branchindex] );
   node.set_branchinfo( branchindex, floor_value, 'u' );
   node.set_dpt( set_dpt );
   node.set_index( index + 1 );

   assert( relax_solval_L != nullptr );

   if ( floor_value && Equal( vars_locallb[branchindex], floor_value, ep ) )
   {
      if( node.alloc_sumfixed() == true )
         node.set_sumfixed( floor_value, pd.get_bvec( branchindex ) );
      else
         node.add_sumfixed( floor_value, pd.get_bvec( branchindex ) );
   }

   relax_solval_L[branchindex] = floor_value;

   vars_localub[branchindex] = floor_value;
}


