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

#define debug 0

bool SVPsolver::SVPSsolve()
{

	if ( TIMELIMIT <= 0 )
		return false;

   // start time
   SVPSstartTime();

	// output bounds
   SVPSoutputBounds();

	// generate oa_cpool
	if ( CUT_OA == true )
   {
      exit(-1);
		gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
	}

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
	int m = probdata.get_m();
	assert( m > 0 );

	if ( subsolver == false && !quiet )
   {
		cout << "Bounds: " << endl;
		for ( int i = 0; i < m; i++ )
			cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl;
	}
}

bool SVPsolver::SVPSrunBranchandBound()
{
   status = SOLVING;

   bool        result = false;
   int         selnodeindex;
	RelaxResult r;

   int   disp = index;
	int	cutoff = 0;

	while ( 1 )
   {
		assert( (int)NodeList.size() > 0 );
		assert( listsize > 0 );
		assert( (int)NodeList.size() == listsize );

		// select a node from the list
      //testwatch.start();
		selnodeindex = SVPSselectNode( index, disp );
      //testwatch.stop();
      auto node = NodeList.begin();
      advance( node, selnodeindex );

		assert( selnodeindex >= 0 );
		assert( selnodeindex < listsize );

		// solve a relaxation problem
		r = SVPSsolveRelaxation( *node );

		if( r == INFEASIBLE || r == GETINTEGER ){
			cutoff++;
		}

		// run heuristics
		if( r == FEASIBLE && HEUR_APP < Appfac ){
			SVPSheur( selnodeindex );
		}

		// output
		if ( !quiet )
      {
         auto buf = ( index - 1 ) % 1000;
         if ( ( buf == 0 || buf == 1 ) && disp < index )
         {
			   //if( subsolver == true )
            //   cout << "t" << omp_get_thread_num() << ":";
			   if( subsolver == true )
               cout << this_thread::get_id() << ":";

			   disp_log(selnodeindex, r, index, cutoff);
            disp = index;
			   cutoff = 0;
         }
		}

		// branch
		if( r == UPDATE || r == FEASIBLE )
      {
			branch( selnodeindex, index);
			index += 2;
		}
      else
      {
		   // remove
		   node = NodeList.erase( node );
		   listsize--;
      }

		// break
		assert( (int)NodeList.size() == listsize );

      auto checktiming = ( index - 1 ) % 1000;
      if ( checktiming == 0 || checktiming == 1 )
      {
         if ( SVPScheckLimit() == true )
            break;
      }

		if( listsize == 0 )
      {
         result = true;
         status = SOLVED;
			if( !subsolver )
            cout << "End" << endl;
			break;
		}

	}

   //cout << "testwatch: " << testwatch.get_result() << endl;

   return result;
}

bool SVPsolver::SVPScheckLimit()
{
   bool result = false;

	if( stopwatch.check_time() == false )
   {
      result = true;
      status = TIMEOVER;
   }
   else if ( listsize >= LEFTNODELIMIT )
   {
      result = true;
      status = FULL_OF_LEFTNODES;
   }
   else if ( index >= NODELIMIT )
   {
      result = true;
      status = FULL_OF_NODES;
   }

   if ( result )
   {
      auto it = NodeList.begin();
		auto min_lb = it->get_lowerbound();
      double lb_i;
		for ( int i = 1; i < listsize; i++ )
      {
			++it;
         lb_i = it->get_lowerbound();
			if( min_lb > lb_i )
				min_lb = lb_i;
		}

		GLB = min_lb;
   }

   return result;
}
