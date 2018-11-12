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

#include <thread>
//#include <omp.h>

#include "svpsolver.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "testwatch.h"
#include "probdata.h"
#include "vector.h"


bool SVPsolver::SVPSsolve()
{

   if ( TIMELIMIT <= 0 )
      return false;

   // output bounds
   SVPSoutputBounds();

   // generate oa_cpool
   //if ( CUT_OA == true )
   //{
   //   exit(-1);
   // gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
   //}

   // start time
   SVPSstartTime();

   // branch-and-bound algorithm
   bool result = SVPSrunBranchandBound();
   stopwatch.stop();

   nnode = (unsigned long int) index -1;

   return result;
}


bool SVPsolver::SVPSresolve()
{

   if ( TIMELIMIT <= 0 )
      return false;

   // start time
   SVPSstartTime();

   // branch-and-bound algorithm
   bool result = SVPSrunBranchandBound();
   stopwatch.stop();

   nnode = (unsigned long int) index -1;

   return result;
}

void SVPsolver::SVPSstartTime()
{
   stopwatch.set_timelimit( TIMELIMIT );
   stopwatch.start();
}

void SVPsolver::SVPSoutputBounds()
{
   assert( !bounds.empty() );
   assert( !bounds[0].empty() );

   if ( subsolver == false && !quiet )
   {
      const auto m = probdata.get_m();
      assert( m > 0 );

      cout << "Bounds: " << endl;
      for ( int i = 0; i < m; i++ )
         cout << "x_" << i << ": [ " <<  - bounds[0][i] << ", " << bounds[0][i] << "]" << endl;
   }
}

bool SVPsolver::SVPSrunBranchandBound()
{
   status = SOLVING;

   // copy
   const auto m = probdata.get_m();
   auto& NL = nodelist;
   const auto output = !quiet;
   const auto subsol = subsolver;
   auto& nodeindex = index;
   const auto enumeration = ENUM;
   //const int dpt_enum = m / 2;
   //const int dpt_enum = 0;

   bool result = false;
   int disp = index;
   int cutoff = 0;
   int checktiming;
   RelaxResult relaxresult = INFEASIBLE;
   double* vars_localub;
   double* vars_locallb;
   NODE* node = nullptr;

   vars_localub = new double[m];
   vars_locallb = new double[m];

   while ( 1 )
   {
      assert( (NL.*check_size)() );
      //testwatch.start();
      //testwatch.stop();

      // select a node from the list
      // get local bounds of variable
      if ( node == nullptr )
      {
         node = (NL.*nodeselection)( &GLB, bestval, nodeindex, disp );
         SVPSgetVarsLocalBound( *node, vars_localub, vars_locallb );
      }
      assert( node != nullptr );

      //for ( auto i = 0; i < m; ++i )
      //   printf("%d:[%f,%f]\n", i, vars_locallb[i], vars_localub[i]);

      // solve a relaxation problem
      //testwatch.start();
      relaxresult = SVPSsolveRelaxation( *node, vars_localub, vars_locallb );
      //testwatch.stop();

      if ( relaxresult == FEASIBLE && enumeration )
         relaxresult = SVPSenumerate( *node, vars_localub, vars_locallb );

      // run heuristics
      if( relaxresult == FEASIBLE && HEUR_APP < Appfac )
         SVPSheur( *node, vars_localub, vars_locallb );

      // branch
      if( relaxresult == UPDATE || relaxresult == FEASIBLE )
      {
         //cout << "branch-start";
         SVPSbranch( *node, nodeindex, vars_localub, vars_locallb );
         nodeindex += 2;
         //cout << "-end-";
      }
      else
      {
         // remove
         //cout << "remove-start";
         node = nullptr;
         (NL.*cut_off)();
         //cout << "-end-";
         cutoff++;

         if( NL.getListsize() == 0 )
         {
            result = true;
            status = SOLVED;
            break;
         }

      }

      //cout << NL.getListsize() << endl;
      // output
      if ( output )
      {
         auto buf = ( nodeindex - 1 ) % 1000;
         if ( ( buf == 0 || buf == 1 ) && disp < nodeindex )
         {
            if( subsol == true )
               cout << this_thread::get_id() << ":";

            disp_log( *node, relaxresult, nodeindex, cutoff );
            disp = nodeindex;
            cutoff = 0;
         }
      }

      // break
      checktiming = ( nodeindex - 1 ) % 1000;
      if ( checktiming == 0 || checktiming == 1 )
      {
         if ( SVPScheckLimit() == true )
            break;
      }
   }

   if( !subsol )
      cout << "End" << endl;

   //cout << "testwatch: " << testwatch.get_result() << endl;

   delete[] vars_localub;
   delete[] vars_locallb;

   return result;
}

void SVPsolver::SVPSgetVarsLocalBound(
      NODE&    node,
      double*  vars_localub,
      double*  vars_locallb
      )
{
   assert( vars_localub != nullptr );
   assert( vars_locallb != nullptr );
   assert( !bounds.empty() );

   //cout << "get local bound - start " << endl;
   const auto m = probdata.get_m();
   const auto& vars_globalbounds = bounds;

   const auto& branchinfo = node.get_branchinfo();
   const int nodetype = node.get_type();
   int branchindex;
   int branchvalue;
   char branchtype;

   //node.NODEdispInformation();
   assert( nodetype >= 0 );
   assert( nodetype < (int)vars_globalbounds.size() );
   assert( !vars_globalbounds[nodetype].empty() );

   for ( int i = 0; i < m; ++i )
   {
      vars_localub[i] = vars_globalbounds[nodetype][i];
      vars_locallb[i] = - vars_localub[i];
   }

   for ( auto& bi : branchinfo )
   {
      branchindex = bi.get_index();
      branchvalue = bi.get_value();
      branchtype = bi.get_type();
      switch ( branchtype )
      {
         case 'u':
            vars_localub[branchindex] = branchvalue;
            break;
         case 'l':
            vars_locallb[branchindex] = branchvalue;
            break;
         case 'e'://EQUAL
            vars_localub[branchindex] = branchvalue;
            vars_locallb[branchindex] = branchvalue;
            break;
         default:
            printf("error (%c) \n", branchtype );
            exit(-1);
            break;
      }
   }
   //cout << "get local bound - end " << endl;
}

bool SVPsolver::SVPScheckLimit()
{
   // copy
   auto& sw = stopwatch;
   auto& NL = nodelist;
   auto& leftnodelimit = LEFTNODELIMIT;
   auto& nodelimit = NODELIMIT;

   bool result = false;

   if( sw.check_time() == false )
   {
      result = true;
      status = TIMEOVER;
   }
   else if ( NL.getListsize() >= leftnodelimit )
   {
      result = true;
      status = FULL_OF_LEFTNODES;
   }
   else if ( index >= nodelimit )
   {
      result = true;
      status = FULL_OF_NODES;
   }

   if ( result )
   {
      GLB = (NL.*get_GLB)();
   }

   return result;
}
