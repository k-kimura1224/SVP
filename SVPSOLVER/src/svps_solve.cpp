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
#include <omp.h>

#include "svpsolver.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "probdata.h"
#include "vector.h"

#define debug 0

bool SVPsolver::solve(
	bool		para,
	bool 		s_zero,
	double	s_GLB
	)
{
   int tlimit = get_timelimit();

	if ( tlimit <= 0 )
		return false;

	stopwatch.set_timelimit( tlimit );
	stopwatch.start();

	int m = probdata.get_m();
	assert( m > 0 );

	// output bounds
	if ( para == false && !quiet )
   {
		cout << "Bounds: " << endl;
		for ( int i = 0; i < m; i++ )
			cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl;
	}

	// generate a root node
	NODE	root;
	int	index = 0;

	double	*init_warm = nullptr;

	if( para == false )
   {
		init_warm = bestsol.get_solval();
	}
   else
   {
		init_warm = new double[m];
		for( int i = 0; i < m; i++ )
			init_warm[i] = (ub[i] + lb[i])/2.0;
	}

	root.set_vals( m, ub, lb, init_warm, s_GLB, 0, s_zero, index);

	for( int i = 0; i < m; i++ )
   {
		if ( lb[i] - ub[i] == 0 && lb[i] != 0.0 )
      {
			if ( root.alloc_sumfixed() == true )
				root.set_sumfixed( lb[i], probdata.get_bvec(i) );
			else
				root.add_sumfixed( lb[i], probdata.get_bvec(i) );
		}
	}

	NodeList.push_back( root );
	listsize++;
	index++;

	if ( para == true )
		delete[] init_warm;

	init_warm = nullptr;

	// generate oa_cpool
	if ( CUT_OA == true )
   {
      exit(-1);
		gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
	}

	int         selnodeindex;
	RelaxResult r;

	int	disp = index;
	int	cutoff=0;
	while ( 1 )
   {


		assert( (int)NodeList.size() > 0 );
		assert( listsize > 0 );
		assert( (int)NodeList.size() == listsize );

		// select a node from the list
		selnodeindex = select_node( index, disp );


		assert( selnodeindex >= 0 );
		assert( selnodeindex < listsize );

		// solve a relaxation problem
		r = solve_relaxation( selnodeindex );

		if( r == INFEASIBLE || r == GETINTEGER ){
			cutoff++;
		}

		// run heuristics
		if( r == FEASIBLE && HEUR_APP < Appfac ){
			heur( selnodeindex, para);
		}

		// output
		if ( !quiet )
      {
         if ((index-1) % 1000 == 0 && disp < index )
         {
			   auto it = NodeList.begin();
			   double min_lb = it->get_lowerbound();
			   for ( int i = 1; i < listsize; i++ )
            {
			   	++it;
			   	if( min_lb > it->get_lowerbound() )
			   		min_lb = it->get_lowerbound();
			   }

			   GLB = min_lb;

			   if( para == true )
               cout << "t" << omp_get_thread_num() << ":";

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
		   auto it = NodeList.begin();
		   advance( it, selnodeindex);
		   NodeList.erase( it );
		   listsize--;
      }

		// break
		assert( (int)NodeList.size() == listsize );

		if( stopwatch.check_time() == false )
      {
         auto it = NodeList.begin();
			double min_lb = it->get_lowerbound();
			for ( int i = 1; i < listsize; i++ )
         {
				++it;
				if( min_lb > it->get_lowerbound() )
					min_lb = it->get_lowerbound();
			}

			GLB = min_lb;
         break;
      }

		if( listsize == 0 )
      {
			if( !para )
            cout << "End" << endl;
			break;
		}

	}

   //cout << testwatch.get_result() << endl;

	nnode = (unsigned long int)index -1;


	if( stopwatch.check_time() == false )
   {
	   stopwatch.stop();
		return false;
   }

	stopwatch.stop();
	return true;
}
