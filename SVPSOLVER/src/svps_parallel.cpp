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

//bool SVPsolver::p_solve()
//{
//   assert( nthreads > 0 );
//
//   if ( TIMELIMIT <= 0 )
//      return false;
//
//   // start time
//   SVPSstartTime();
//
//   // output bounds
//   SVPSoutputBounds();
//
//   // generate oa_cpool
//   if( CUT_OA == true )
//   {
//      exit(-1);
//      gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
//   }
//
//   int         sel;
//   int m = probdata.get_m();
//   assert( m > 0 );
//
//   RelaxResult r;
//   int   disp = index;
//   int   cutoff=0;
//   bool  result;
//   while(1){
//      assert( (int)NodeList.size() > 0 );
//      assert( listsize > 0 );
//      assert( (int)NodeList.size() == listsize );
//
//      if( listsize >= NUM_INITNODES ){
//         result = false;
//         break;
//      }
//      // select a node from the list
//      sel = select_node(index, disp );
//
//      assert( sel >= 0 );
//      assert( sel < listsize );
//
//      // solve a relaxation problem
//      r = solve_relaxation( sel );
//
//      if( r == INFEASIBLE || r == GETINTEGER ){
//         cutoff++;
//      }
//
//      // run heuristics
//      if( r == FEASIBLE && HEUR_APP < Appfac ){
//         SVPSheur( sel );
//      }
//
//      // output
//      if( !quiet && (index-1) % 1000 == 0 && disp < index ){
//         auto it = NodeList.begin();
//         double min_lb = it->get_lowerbound();
//         for(int i=1; i<listsize; i++){
//            ++it;
//            if( min_lb > it->get_lowerbound() ){
//               min_lb = it->get_lowerbound();
//            }
//         }
//         GLB = min_lb;
//         disp_log(sel, r, index, cutoff);
//         disp = index;
//         cutoff = 0;
//      }
//
//      // branch
//      if( r == UPDATE || r == FEASIBLE ){
//         branch( sel, index);
//         index += 2;
//      }
//      else
//      {
//         // remove
//         auto it = NodeList.begin();
//         advance( it, sel);
//         NodeList.erase( it );
//         listsize--;
//      }
//
//      // break
//      assert( (int)NodeList.size() == listsize );
//
//      if( stopwatch.check_time() == false )
//      {
//         result = true;
//         break;
//      }
//      if( listsize == 0 ){
//         cout << "End" << endl;
//         result = true;
//         break;
//      }
//
//   }
//
//   nnode = (unsigned long int)index -1;
//
//   if( result == true )
//   {
//      stopwatch.stop();
//      return true;
//   }
//
//   auto it = NodeList.begin();
//   double min_lb = it->get_lowerbound();
//   for(int i=1; i<listsize; i++){
//      ++it;
//      if( min_lb > it->get_lowerbound() ){
//         min_lb = it->get_lowerbound();
//      }
//   }
//   GLB = min_lb;
//
//
//   // parallel {{
//   NodeList.sort();
//
//   cout << "Parallel mode ---------------------------------------------";
//   cout << "-------------------------------------------" << endl;
//   assert( (int)NodeList.size() == listsize );
//
//   bool        p_run;
//   double      p_min;
//
//#pragma omp parallel num_threads(nthreads)
//   {
//   #pragma omp for schedule(dynamic, PARASIZE) private(p_run,p_min)
//   for(int i=0; i<listsize; i++){
//      auto it = NodeList.begin();
//      advance( it, i);
//
//      SVPsolver sub;
//      int sub_timelimit = TIMELIMIT - stopwatch.get_time();
//
//      sub.SVPSsetup( m, probdata.get_B_(), 1, sub_timelimit, quiet,
//            true, it->get_zero(), false, false, false, false );
//      sub.SVPSsetGlobalLowerBound( it->get_lowerbound() );
//      sub.SVPSsetBestval( bestval );
//      sub.SVPSsetBounds( it->get_ub(), it->get_lb());
//      if( BRANCHINGRULE_INT == 4 || BRANCHINGRULE_INT == 5 ) sub.SVPSsetOrder( order );
//      sub.SVPSsetNorm( norm );
//      sub.SVPSsetAppfac( Appfac, _Appfac );
//      //sub.SVPSgenerateRootNode( it->get_zero() );
//      sub.SVPSsetNode( *it );
//
//      p_run = sub.SVPSsolve();
//
//      if( !p_run ){
//         #pragma omp critical
//         {
//            nnode += sub.SVPSgetNnode();
//         }
//      }else{
//
//         it->set_solved( true );
//
//         int p_ct = 0;
//         if( i%(int)PARASIZE == 0 ){
//            p_min = bestval;
//            list<NODE>::iterator itr = NodeList.begin();
//            for(int j=0; j<listsize; ++j,++itr) {
//               if( itr->get_solved() == true ){
//                  p_ct++;
//               }
//               if( itr->get_solved() == false && p_min > itr->get_lowerbound() ){
//                  p_min = itr->get_lowerbound();
//               }
//            }
//            #pragma omp critical
//            {
//               GLB = p_min;
//            }
//         }
//
//         #pragma omp critical
//         {
//            if( bestval > sub.SVPSgetBestval() ){
//               bestval = sub.SVPSgetBestval();
//               bestsol = sub.SVPSgetBestsol();
//               pool.add_solution( bestsol );
//               Appfac = sqrt( bestval ) / _Appfac;
//            }
//
//            nnode += sub.SVPSgetNnode();
//
//            if( i%(int)PARASIZE == 0 ){
//               cout << "*";
//               cout << stopwatch.get_time() << "s: ";
//               cout << "GLB=" << GLB << ", ";
//               cout << "best=" << bestval << ", ";
//               cout << "gap=" << 100*(bestval - GLB)/bestval << "%, ";
//               cout << "AF=" << Appfac << ", ";
//               cout << "#solved=" << p_ct;
//               cout << endl;
//            }
//
//         }
//      }
//   }
//   }
//   cout << "End" << endl;
//   // }} parallel
//
//   double   min = bestval;
//   list<NODE>::iterator itr = NodeList.begin();
//   for(int j=0; j<listsize; ++j,++itr) {
//      if( itr->get_solved() == false && min > itr->get_lowerbound() ){
//         min = itr->get_lowerbound();
//      }
//   }
//   GLB = min;
//
//   if( GLB == bestval ){
//      NodeList.clear();
//      listsize = 0;
//   }
//
//   stopwatch.stop();
//   return true;
//
//}

void SVPsolver::SVPSsolveSubprob(
      unsigned short int&  n_running_threads,
      double&              sublb_i,
      const int            thread_id
      )
{
   SVPsolver   sub;
   double      sub_timelimit;

   // lock {
   {
      lock_guard<mutex> lock(mtx);

      auto  m = probdata.get_m();
      auto  B_ = probdata.get_B_();
      sub_timelimit = TIMELIMIT - stopwatch.get_time();

      sub.SVPSsetup( m, B_, 1, sub_timelimit, true,
               true, true, false, false, false, false );
      sub.SVPSsetBounds( ub, lb );
      if( BRANCHINGRULE_INT == 4 || BRANCHINGRULE_INT == 5 )
         sub.SVPSsetOrder( order );
      sub.SVPSsetNorm( norm );

      sub.SVPSsetGlobalLowerBound( GLB );
      sub.SVPSsetBestval( bestval );
      sub.SVPSsetAppfac( Appfac, _Appfac );

      sub.SVPSsetNode( *(NodeList.begin()) );
      SVPSpopNode();
   }
   // } lock

   constexpr int  init_sub_nodelimit = 100000;
   constexpr int  max_sub_nodelimit = init_sub_nodelimit * 20;

   unsigned long int totalpop = 0;

   bool  result;
   int   sub_nodelimit = init_sub_nodelimit;

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

            cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
            cout << "get new solution ";
            cout << "-- BESTVAL: " << bestval << "  NORM: " << sqrt( bestval ) << "  AP: " << Appfac << " --";
            cout << endl;
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

            if ( NodeList.empty() )
            {
               assert( listsize == 0 );
               n_running_threads--;

               if ( n_running_threads == 0 )
               {
                  cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
                  cout << "break (SOLVED)" << endl;

                  result = true;
                  status = SOLVED;
                  cv.notify_all();

                  break;
               }

               cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
               cout << "sleep" << endl;
               cv.wait( lock, [this, &n_running_threads]
                     { return ( !NodeList.empty() || n_running_threads == 0 ); } );

               if ( n_running_threads == 0 )
               {
                  assert( NodeList.empty() );
                  assert( status == SOLVED );

                  cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
                  cout << "break (SOLVED)" << endl;

                  assert( nnode > totalpop );
                  nnode -= totalpop;

                  result = true;

                  break;
               }

               n_running_threads++;
            }

            sub.SVPSresetIndex();

            int pushsize = listsize * 0.1;

            if ( listsize < 100 )
               pushsize = listsize;
            else if ( pushsize < 100 )
               pushsize = 100;

            double   min_lb = NodeList.begin()->get_lowerbound();
            double   lb_i;

            for ( int i = 0; i < pushsize; i++ )
            {
               assert( !NodeList.empty() );

               lb_i = NodeList.begin()->get_lowerbound();

               if ( min_lb > lb_i )
                  min_lb = lb_i;

               sub.SVPSsetNode( *(NodeList.begin()) );
               SVPSpopNode();
            }

            sub.SVPSsetGlobalLowerBound( min_lb );

		      assert( (int)NodeList.size() == listsize );

            sub_nodelimit = init_sub_nodelimit;

            cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
            cout << "push->" << pushsize << endl;
         }
         else
         {
            Status substatus = sub.SVPSgetStatus();

            assert( substatus == FULL_OF_NODES
                  || substatus == FULL_OF_LEFTNODES
                  || substatus == TIMEOVER );

            if ( substatus == FULL_OF_NODES )
            {
               assert( sub.SVPSgetListsize() > 0 );
               auto sub_listsize = sub.SVPSgetListsize();
               if ( sub_listsize > listsize )
               {
                  int popsize = sub_listsize * 0.5;

                  totalpop += popsize;

                  for ( int i = 0; i < popsize; i++ )
                  {
                     SVPSsetNode( sub.SVPSgetNode() );
                     sub.SVPSpopNode();
                  }

                  cv.notify_all();

                  cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
                  cout << "pop->" << popsize << endl;
               }

               if ( max_sub_nodelimit > sub_nodelimit )
                  sub_nodelimit *= 2;
               else
                  sub_nodelimit += max_sub_nodelimit;
            }
            else
            {
               cout << stopwatch.get_time() << "s (left: " << listsize << ") [t" << thread_id << "] ";
               cout << "break (" << sub.SVPSgetStringStatus() << ")" << endl;

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
   if( CUT_OA == true )
   {
      exit(-1);
      gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
   }

   assert( probdata.get_m() >= 40 );
   // branch-and-bound algorithm
   int   ORIG_LEFTNODELIMIT = LEFTNODELIMIT;
   LEFTNODELIMIT = 2 * nthreads;
   bool  result = SVPSrunBranchandBound();
   LEFTNODELIMIT = ORIG_LEFTNODELIMIT;

   nnode = (unsigned long int)index -1;

   if( result == true )
   {
      stopwatch.stop();
      return true;
   }

   NodeList.sort();
   status = SOLVING;

   // parallel {{

   cout << "Parallel mode ---------------------------------------------";
   cout << "-------------------------------------------" << endl;
   assert( (int)NodeList.size() == listsize );

   cout << "-- BESTVAL: " << bestval << "  NORM: " << sqrt( bestval ) << "  AP: " << Appfac << " --";
   cout << endl;

   vector<thread>       threads;
   unsigned short int   n_running_threads = nthreads;
   double*              sublb = new double[nthreads];

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
      double   min_lb = bestval;
      double   node_lb;
      for ( auto& node : NodeList )
      {
         node_lb = node.get_lowerbound();
         if ( node_lb < min_lb )
            min_lb = node_lb;
      }
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
   assert( NodeList.empty() );
   assert( listsize == 0 );
   index = 0;
}
