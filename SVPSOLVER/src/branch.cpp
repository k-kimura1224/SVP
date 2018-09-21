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
#include "node.h"
#include "probdata.h"
#include "qpsolver.h"
#include "solution.h"
#include "vector.h"

using namespace std;

#define debug  0

static int find_branchingvariable1(
      int      m,
      double   *relax_solval
   )
{

   double   eval;
   double   max_eval = 0.0;
   int memo = -1;
   double buf;

   for ( int i = 0; i < m; i++ )
   {
      buf = relax_solval[i];
      eval = pow( buf - round(buf), 2);

      if ( max_eval <= eval )
      {
         max_eval = eval;
         memo = i;
      }
   }

   return memo;
}

static int find_branchingvariable2(
      int      m,
      double   *relax_solval,
      double   *ub,
      double   *lb
   )
{

   double   min = 1.0;
   int      memo = -1;
   double   r;
   double   buf;
   double   epsilon = 1.0e-12;
   for ( int i = 0; i < m; i++ )
   {
      buf = relax_solval[i];

      if ( !Equal( ub[i], lb[i], epsilon )
         && pow( buf, 2 ) < 1.0
         && !Equal( buf, 0.0, epsilon ) )
      {
         if ( buf >  0.0 )
            r = buf;
         else
            r = - buf;

         if( min > r )
         {
            min = r;
            memo = i;
         }
      }
   }

   return memo;
}

static int find_branchingvariable3(
      int      m,
      double   *relax_solval,
      double   *ub,
      double   *lb
   )
{

   for(int i=m-1; 0<=i; i--){
      if(  pow(relax_solval[i],2) < 1.0
         && relax_solval[i] != 0.0 ){
         return i;
      }
   }

   return -1;
}

static int find_branchingvariable4(
      int      m,
      double   *relax_solval
   )
{
   double buf;
   for ( int i = m-1; 0 <= i; i-- )
   {
      buf = relax_solval[i];
      if( !Equal( buf, round( buf ), 1.0e-12 ) )
         return i;
   }

   return -1;
}

static int find_branchingvariable5(
      int      m,
      double   *relax_solval,
      int      *order
   )
{

   assert( order != nullptr );

   int memo = -1;
   for(int i=0; i<m; i++){
      assert( order[i] >= 0 && order[i] < m );
      if( relax_solval[order[i]] - round(relax_solval[order[i]]) != 0){
         memo = order[i];
         //cout << memo << endl;
         break;
      }
   }

   assert( memo >= 0 );
   assert( memo < m );

   return memo;
}

static int find_branchingvariable6(
      int      m,
      double   *relax_solval,
      int      *order
   )
{

   assert( order != nullptr );

   int memo = -1;
   for(int i=m-1; 0<=i; i--){
      assert( order[i] >= 0 && order[i] < m );
      if( relax_solval[order[i]] - round(relax_solval[order[i]]) != 0){
         memo = order[i];
         //cout << memo << endl;
         break;
      }
   }

   assert( memo >= 0 );
   assert( memo < m );

   return memo;
}

static int find_branchingvariable7(
      int      m,
      double   *relax_solval
   )
{
   double buf;
   for ( int i = m-1; 0 <= i; i-- )
   {
      buf = relax_solval[i];
      if( fabs( buf ) < 1.0 && !Equal( buf, round( buf ), 1.0e-12 ) )
         return i;
   }

   return -1;
}

void  SVPsolver::branch(
   int   sel,
   int   index
   )
{
   auto it = NodeList.begin();
   advance( it, sel);

   if( it->get_zero() == true )
      branch_BIN( sel, index );
   else
      branch_INT( sel, index );

}

void  SVPsolver::branch_BIN(
   int   sel,
   int   index
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

   auto it = NodeList.begin();
   advance( it, sel );

   double*  warm = it->get_warm();

   set_ub = new double[m];
   set_lb = new double[m];
   set_warm = new double[m];

   // x_i = 0
   int memo = -1;
   Copy_vec( it->get_ub(), set_ub, m );
   Copy_vec( it->get_lb(), set_lb, m );
   for ( int i = 0; i < m; i++ )
   {
      if ( set_lb[i] != 0.0 || set_ub[i] != 0.0 )
      {
         set_lb[i] = 0.0;
         set_ub[i] = 0.0;
         memo = i;
         break;
      }
   }

   for( int i = 0; i < m; i++ )
   {
      if( warm[i] >= set_lb[i] && warm[i] <= set_ub[i] )
         set_warm[i] = warm[i];
      else
         set_warm[i] = (set_ub[i]+set_lb[i])/2.0;
   }

   set_relax_objval = it->get_lowerbound();
   set_dpt = it->get_dpt() + 1;
   set_zero = true;
   set_index = index;

   tighten_bounds( memo, set_lb, set_ub );

   C_BIN.set_vals( m, set_ub, set_lb,
                   set_warm, set_relax_objval,
                   set_dpt,  set_zero, set_index);

   if( it->get_sumfixed() != nullptr )
   {
      bool r;
      r = C_BIN.alloc_sumfixed();
      assert( r == true );
      C_BIN.set_sumfixed( 1.0, it->get_sumfixed() );
   }

   if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 )
   {
      if( C_BIN.alloc_sumfixed() == true )
         C_BIN.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
      else
         C_BIN.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
   }

   NodeList.push_back( C_BIN );
   listsize++;

   assert( memo >= 0 );

   // x_i >= 1
   auto  ub = it->get_ub();
   auto  lb = it->get_lb();
   warm = it->get_warm();

   lb[memo] = 1.0;

   for ( int i = 0; i < m; i++ )
   {
      if( warm[i] < lb[i] || warm[i] > ub[i] )
         warm[i] = lb[i];
   }

   assert( it->get_zero() == true );

   it->set_dpt( it->get_dpt() + 1 );
   it->set_zero( false );
   it->set_index( index + 1 );

   if( lb[memo] == ub[memo] && lb[memo] != 0.0 )
   {
      if( it->alloc_sumfixed() == true )
         it->set_sumfixed( lb[memo], probdata.get_bvec(memo) );
      else
         it->add_sumfixed( lb[memo], probdata.get_bvec(memo) );
   }

   delete[] set_ub;
   delete[] set_lb;

}
void  SVPsolver::branch_INT(
   int   sel,
   int   index
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

   auto it = NodeList.begin();
   advance( it, sel);
   double   *relax_solval = it->get_relaxsolval();

   set_ub = new double[m];
   set_lb = new double[m];
   set_warm = new double[m];

   assert( relax_solval != nullptr );

   // find a branching variable
   int memo = -1;
   switch( BRANCHINGRULE_INT )
   {
      case 0:
      {
         memo = find_branchingvariable1( m, relax_solval );
         break;
      }
      case 1:
      {
         memo = find_branchingvariable2( m, relax_solval, it->get_ub(), it->get_lb() );
         if( memo == -1 )
            memo = find_branchingvariable1( m, relax_solval );
         break;
      }
      case 2:
      {
         memo = find_branchingvariable3( m, relax_solval, it->get_ub(), it->get_lb() );
         if( memo == -1 )
            memo = find_branchingvariable1( m, relax_solval );
         break;
      }
      case 3:
      {
         //printv( m, relax_solval);
         memo = find_branchingvariable4( m, relax_solval );
         break;
      }
      case 4:
      {
         memo = find_branchingvariable5( m, relax_solval, order );
         break;
      }
      case 5:
      {
         memo = find_branchingvariable6( m, relax_solval, order );
         break;
      }

      case 6:
      {
         memo = find_branchingvariable7( m, relax_solval );
         if ( memo == -1 )
            memo = find_branchingvariable4( m, relax_solval );
         break;
      }

      default:
      {
         cout << "error: branch.cpp" << endl;
         exit(-1);
         break;
      }
   }

   assert( memo >= 0 );
   assert( memo < m );

   //cout << relax_solval[memo];
   //cout << ", " << memo << endl;


   // ceil <= x_memo
   Copy_vec( it->get_ub(), set_ub, m);
   Copy_vec( it->get_lb(), set_lb, m);

   set_lb[memo] = ceil( relax_solval[memo] );

   Copy_vec( relax_solval, set_warm, m );

   set_warm[memo] = set_lb[memo];

   set_relax_objval = it->get_lowerbound();
   set_dpt = it->get_dpt() + 1;
   set_zero = false;
   set_index = index;

   C_RIGHT.set_vals( m, set_ub, set_lb,
                   set_warm, set_relax_objval,
                   set_dpt,  set_zero, set_index );

   if( it->get_sumfixed() != nullptr )
   {
      bool r;
      r = C_RIGHT.alloc_sumfixed();
      assert( r == true );
      C_RIGHT.set_sumfixed( 1.0, it->get_sumfixed() );
   }

   if( set_lb[memo] == set_ub[memo] && set_lb[memo] != 0.0 )
   {
      if( C_RIGHT.alloc_sumfixed() == true )
         C_RIGHT.set_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
      else
         C_RIGHT.add_sumfixed( set_lb[memo], probdata.get_bvec(memo) );
   }

   NodeList.push_back( C_RIGHT );
   listsize++;

   // x_memo <= floor
   auto  ub = it->get_ub();
   auto  lb = it->get_lb();
   auto  warm = it->get_warm();

   ub[memo] = floor( relax_solval[memo] );

   Copy_vec( relax_solval, warm, m );
   warm[memo] = ub[memo];

   it->set_dpt( it->get_dpt() + 1 );
   it->set_zero( false );
   it->set_index( index + 1 );

   if ( Equal( lb[memo], ub[memo], epsilon ) && !Equal( lb[memo], 0.0, epsilon ) )
   {
      if( it->alloc_sumfixed() == true )
         it->set_sumfixed( lb[memo], probdata.get_bvec(memo) );
      else
         it->add_sumfixed( lb[memo], probdata.get_bvec(memo) );
   }

   delete[] set_ub;
   delete[] set_lb;
   delete[] set_warm;
}


