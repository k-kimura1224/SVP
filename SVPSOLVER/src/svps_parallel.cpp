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
#include <mutex>
//#include <omp.h>

#include "svpsolver.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "vector.h"
//#include "probdata.h"

#define debug 0

void SVPsolver::SVPSsolveSubprob(
      int&        n_running_threads,
      double&     sublb_i,
      const int   thread_id
      )
{
   SVPsolver   sub;
   double      sub_timelimit;

   // lock {
   {
      lock_guard<mutex> lock(mtx);

      const auto m = probdata.get_m();
      const auto B_ = probdata.get_B_();
      sub_timelimit = TIMELIMIT - stopwatch.get_time();

      sub.SVPSsetup( m, B_, 1, sub_timelimit, true,
               true, false, false, false, false, false );
      sub.SVPSsetBounds( bounds );
      sub.SVPSsetNorm( norm );

      sub.SVPSsetGlobalLowerBound( GLB );
      sub.SVPSsetBestval( bestval );
      sub.SVPSsetupNodelist();
      sub.SVPSsetAppfac( Appfac, _Appfac );

      const auto setup = (nodelist.*setup_para_selection)();
      sub.SVPSmoveNode( (nodelist.*para_selection)(setup) );
      SVPSpopNode( setup );
   }
   // } lock

   constexpr int init_sub_nodelimit = 100000;
   constexpr int max_sub_nodelimit = init_sub_nodelimit * 20;

   unsigned long int totalpop = 0;

   bool result;
   int sub_nodelimit = init_sub_nodelimit;
   const bool disp = !quiet;
   const double pop_maxrate = 1.0 - ( 1.0 / (double) nthreads );
   assert( pop_maxrate > 0.0 && pop_maxrate <= 1.0 );

   while ( 1 )
   {
      sub.SVPSsetNodelimit( sub_nodelimit );
      result = sub.SVPSresolve();

      // lock {
      {
         unique_lock<mutex> lock(mtx);
         auto subbestval = sub.SVPSgetBestval();

         if ( bestval > subbestval )
         {
            bestval = subbestval;
            bestsol = sub.SVPSgetBestsol();
            pool.add_solution( bestsol );
            Appfac = sqrt( bestval ) / _Appfac;

            if ( disp )
            {
               cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
               cout << "get new solution ";
               cout << "-- BESTVAL: " << bestval << "  NORM: " << sqrt( bestval ) << "  AP: " << Appfac << " --";
               cout << endl;
            }
         }
         else if ( subbestval > bestval )
         {
            sub.SVPSsetBestval( bestval );
            sub.SVPSsetAppfac( Appfac, _Appfac );
         }

         if ( result == true )
         {
            assert( sub.SVPSgetStatus() == SOLVED );
            nnode += sub.SVPSgetNnode();

            if ( !nodelist.getListsize() )
            {
               assert( nodelist.getListsize() == 0 );
               n_running_threads--;

               if ( n_running_threads == 0 )
               {
                  if ( disp )
                  {
                     cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
                     cout << "break (SOLVED)" << endl;
                  }

                  result = true;
                  status = SOLVED;
                  cv.notify_all();

                  break;
               }

               if ( disp )
               {
                  cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
                  cout << "sleep" << endl;
               }

               cv.wait( lock, [this, &n_running_threads]
                     { return ( nodelist.getListsize() || n_running_threads == 0 ); } );

               if ( n_running_threads == 0 )
               {
                  assert( !nodelist.getListsize() );
                  assert( status == SOLVED );

                  if ( disp )
                  {
                     cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
                     cout << "break (SOLVED)" << endl;
                  }

                  assert( nnode > totalpop );
                  nnode -= totalpop;

                  result = true;

                  break;
               }

               n_running_threads++;
            }

            sub.SVPSresetIndex();
            sub.SVPSsetupNodelist();

            const auto setup = (nodelist.*setup_para_selection)();
            const auto nodelistsize = nodelist.getListsize();
            const auto maxsize = (nodelist.*getSubsize)( setup );
            int pushsize = nodelistsize * 0.1;

            if ( maxsize <= 100 )
               pushsize = maxsize;
            else if ( pushsize <= 100 )
               pushsize = 100;

            assert( pushsize <= maxsize );

            double   min_lb = (nodelist.*para_selection)(setup).get_lowerbound();

            for ( int i = 0; i < pushsize; i++ )
            {
               assert( (nodelist.*getSubsize)( setup ) > 0 );

               auto&    node = (nodelist.*para_selection)(setup);
               double   lb_i = node.get_lowerbound();

               if ( min_lb > lb_i )
                  min_lb = lb_i;

               sub.SVPSmoveNode( node );
               SVPSpopNode( setup );
            }

            sub.SVPSsetGlobalLowerBound( min_lb );

            assert( (nodelist.*check_size)() );

            sub_nodelimit = init_sub_nodelimit;

            if ( disp )
            {
               cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
               cout << "push->" << pushsize << endl;
            }
         }
         else
         { // result == false
            Status substatus = sub.SVPSgetStatus();

            assert( substatus == FULL_OF_NODES
                  || substatus == FULL_OF_LEFTNODES
                  || substatus == TIMEOVER );

            if ( substatus == FULL_OF_NODES )
            {
               assert( sub.SVPSgetListsize() > 0 );
               auto sub_listsize = sub.SVPSgetListsize();
               if ( sub_listsize > nodelist.getListsize() )
               {
                  int popsize = 0;
                  const auto setup = sub.SVPSgetSetupParaSelection();
                  const auto maxsize = sub.SVPSgetSubsize( setup );

                  assert( pop_maxrate > 0.0 && pop_maxrate <= 1.0 );

                  if ( n_running_threads > 1 )
                     popsize = sub_listsize * 0.5;
                  else
                     popsize = sub_listsize * pop_maxrate;

                  if ( maxsize < popsize )
                     popsize = maxsize;

                  assert( popsize > 0 && popsize <= sub_listsize );

                  totalpop += popsize;

                  for ( int i = 0; i < popsize; i++ )
                  {
                     SVPSmoveNode( sub.SVPSgetNode_para_selection(setup) );
                     sub.SVPSpopNode( setup );
                  }

                  cv.notify_all();

                  if ( disp )
                  {
                     cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
                     cout << "pop->" << popsize << endl;
                  }
               }

               if ( max_sub_nodelimit > sub_nodelimit )
                  sub_nodelimit *= 2;
               else
                  sub_nodelimit += max_sub_nodelimit;
            }
            else
            {
               if ( disp )
               {
                  cout << stopwatch.get_time() << "s (left: " << nodelist.getListsize() << ") [t" << thread_id << "] ";
                  cout << "break (" << sub.SVPSgetStringStatus() << ")" << endl;
               }

               status = substatus;
               result = false;

               nnode += sub.SVPSgetNnode();
               nnode -= totalpop;

               assert( nnode > totalpop );
               break;
            }
         }
      }
      // } lock
   }

   if ( result == false )
   {
      assert( status != SOLVED && status != SOLVING );
      sublb_i = SVPSgetGlobalLowerBound();
   }
}

bool SVPsolver::SVPSparasolve()
{
   assert( nthreads > 0 );

   if ( TIMELIMIT <= 0 )
      return false;

   // start time
   SVPSstartTime();

   // output bounds
   SVPSoutputBounds();

   // generate oa_cpool
   //if( CUT_OA == true )
   //{
   //   exit(-1);
   //   gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
   //}

   assert( probdata.get_m() >= 40 );
   // branch-and-bound algorithm
   int   ORIG_LEFTNODELIMIT = LEFTNODELIMIT;
   LEFTNODELIMIT = 1000 * nthreads;
   bool  result = SVPSrunBranchandBound();
   LEFTNODELIMIT = ORIG_LEFTNODELIMIT;

   nnode = (unsigned long int)index -1;

   if( result == true )
   {
      stopwatch.stop();
      return true;
   }

   status = SOLVING;

   if ( type == LIST )
      nodelist.sort();

   if ( type == TWO_DEQUE )
      nodelist.setup( type, bestval, 0.3 );

   // parallel {{

   cout << "Parallel mode ---------------------------------------------";
   cout << "-------------------------------------------" << endl;
   assert( (nodelist.*check_size)() );

   cout << "-- BESTVAL: " << bestval << "  NORM: " << sqrt( bestval ) << "  AP: " << Appfac << " --";
   cout << endl;

   vector<thread> threads;
   int            n_running_threads = nthreads;
   double*        sublb = new double[nthreads];

   for ( int i = 1; i < nthreads; i++ )
   {
      threads.push_back( thread( [this, &n_running_threads, &sublb, i]()
               { this->SVPSsolveSubprob( n_running_threads, sublb[i], i ); } ) );
   }

   SVPSsolveSubprob( n_running_threads, sublb[0], 0 );

   for ( auto &th: threads )
      th.join();

      assert( status != SOLVING );

   if ( status == SOLVED )
      result = true;
   else
   {
      result = false;
      double   min_lb = (nodelist.*get_GLB)();
      double   node_lb;
      assert( sublb != nullptr );
      for ( int i = 0; i < nthreads; i++ )
      {
         node_lb = sublb[i];
         if ( node_lb < min_lb )
            min_lb = node_lb;
      }
   }

   delete[] sublb;

   stopwatch.stop();
   return result;

}

void SVPsolver::SVPSresetIndex()
{
   assert( nodelist.getListsize() == 0 );
   index = 0;
}
