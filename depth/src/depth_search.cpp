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
#include <future>
#include <unistd.h>
#include <shared_mutex>
#include <mutex>
#include <condition_variable>

#include "svpsolver.h"
#include "node.h"
#include "vector.h"

using namespace std;

#define debug  0
#define depth_debug  0

class DNode
{
   int ub;
   int lb;
   double lower;
   vector<double> rsol;

   public:
      DNode( const int s_ub, const int s_lb,
               const double s_lower, vector<double>& s_rsol )
      { ub = s_ub; lb = s_lb; lower = s_lower; rsol = move( s_rsol ); }
      DNode( const DNode &source )
      { ub = source.ub; lb = source.lb; lower = source.lower; rsol = source.rsol; }
      DNode& operator=( const DNode& source ) {
         if ( this != &source ) {
            ub = source.ub; lb = source.lb; lower = source.lower; rsol = source.rsol;
         }
         return *this;
      }
      ~DNode() { rsol.clear(); rsol.shrink_to_fit(); }

      auto DNode_ub() const { return ub; }
      auto DNode_lb() const { return lb; }
      auto DNode_lower() const { return lower; }
      auto DNode_rsol() const { return &rsol[0]; }
      auto DNode_rsol_i( const int i ) const { return rsol[i]; }
      auto DNode_rsol_last( const int i ) const { return rsol.back(); }
};

void setup(
      const auto     m,
      deque<NODE>&   NodeList
      )
{
   assert( NodeList.empty() );

   vector<int> empty_v;
   NodeList.emplace_back( 0, 0, 0.0, empty_v );
   NodeList.front().geti_rsol().assign( m, 0.0 );
}

bool SVPsolver::SVPSdepthSearch()
{
   auto bval = bestval;
   constexpr auto borde_rate = 0.6;
   auto border = borde_rate * bestval;
   auto bestsolval = bestsol.get_solval();

   auto nodeindex = index;
   bool result = false;
   const auto m = probdata.get_m();
   const auto& vsco = VsCO;
   const auto& vcov = VCOV;
   const auto q11 = *(probdata.get_Q());
   const auto nsubsolvers = nthreads - 1;

   assert( nthreads >=2 );

   NODE* node = nullptr;
   deque<NODE> NodeList;
   vector<deque<NODE>> NodeList_sub( nsubsolvers );

   vector<thread> threads;
   shared_mutex solmtx;
   vector<mutex> nodemtxs( nsubsolvers );
   vector<condition_variable> cond_vars( nsubsolvers );
   vector<bool> run( nsubsolvers, false );
   int push_id = 0;

   auto next_id = [nsubsolvers] ( int id ) {
      return ( (++id < nsubsolvers)? id : 0 ); };

   // setup
   setup( m, NodeList );

   assert( !NodeList.empty() );
   assert( (int) NodeList.size() == 1 );

   auto vcov_inv = vcov;
   for ( auto& v : vcov_inv )
      v = 1.0 / v;

   auto branch = [&NodeList,&nodeindex,&bval,&vsco,&vcov,&bestsolval,&border,&solmtx,this,m,q11,borde_rate] (NODE& node) {

      solmtx.lock_shared();
      auto currentbestval = bval;
      solmtx.unlock_shared();

      #if debug
      cout << "[Fix";
      cout << " best:" << (int)currentbestval;
      cout << " index:" << (int) node.get_index() << "]-------------------------------" << endl;
      #endif

      // copy member variables of node
      const auto dpt = node.get_dpt();
      const auto dpt1 = dpt + 1;
      const auto& fixedvalues = node.cget_fixedvalues();
      const auto lowerbound = node.get_lowerbound();

      assert( ( fixedvalues.empty() && dpt == 0 ) ||
            ( !fixedvalues.empty() && dpt == (int)fixedvalues.size() ) );

      #if debug
      cout << "fixedvalues: ";
      for ( auto val : fixedvalues )
         cout << val << " ";
      cout << endl;
      #endif

      const auto branchindex = m - 1 - dpt;

      assert( branchindex >= 0 );
      assert( branchindex < m );

      auto& rsol = node.geti_rsol();
      const auto lastvalue = rsol[branchindex];

      assert( !rsol.empty() );
      assert( m - dpt <= (int) rsol.size() );

      if ( dpt == m - 1 )
      {
         double fixvalue = round( lastvalue );
         double upperbound_new = fixvalue;
         upperbound_new -= lastvalue;
         upperbound_new *= upperbound_new;
         upperbound_new *= q11;
         upperbound_new += lowerbound;

         #if debug
         printf("x_%d = %d -> upper: %d (%d)\n", branchindex, (int)fixvalue, (int) upperbound_new, nodeindex);
         #endif
         if ( upperbound_new < currentbestval && (int) upperbound_new )
         {
            solmtx.lock();
            bestsolval[0] = fixvalue;
            assert( (int) fixedvalues.size() == dpt );
            for ( auto j = 1, k = dpt - 1; j < m; ++j, --k )
               bestsolval[j] = fixedvalues[k];

            currentbestval = compute_objval( bestsolval );
            bval = currentbestval;
            border = currentbestval * borde_rate;
            solmtx.unlock();
         }
         return true;
      }

      assert( dpt < m - 1 );

      const auto ceilval = (int) ceil( lastvalue );
      const auto floorval = (int) floor( lastvalue );
      const auto& vsco_i = vsco[branchindex -1];
      const auto vcov_i = vcov[branchindex -1];
      const auto m_dpt1 = m - dpt1;
      double init_lowerbound_new;
      int init_fixvalue;
      double sign = 1.0;

      assert( !vsco_i.empty() );
      assert( m - dpt1 == (int) vsco_i.size() );

      if ( ceilval - lastvalue < lastvalue - floorval )
      {
         sign = - 1.0;
         init_fixvalue = ceilval;
      }
      else
         init_fixvalue = floorval;

      assert( init_fixvalue == round( lastvalue ) );

      double init_coef = init_fixvalue;
      init_coef -= lastvalue;

      init_lowerbound_new = init_coef;
      init_lowerbound_new *= init_coef;
      init_lowerbound_new *= vcov_i;
      init_lowerbound_new += lowerbound;

      #if debug
      printf("x_%d = %d -> lower: %f (-)\n", branchindex, (int) init_fixvalue , init_lowerbound_new );
      #endif

      if ( init_lowerbound_new < currentbestval )
      {
         double lowerbound_new;
         double coef;
         auto counter = 2;
         auto fixvalue = init_fixvalue;
         fixvalue += sign;

         function<void(void)> update_fix;
         if ( lowerbound < 1.0e-10 )
            update_fix = [&fixvalue]{ ++fixvalue; };
         else
          update_fix = [&fixvalue, &sign, &counter]{ sign = -sign; fixvalue += sign * counter; ++counter; };

         while ( 1 )
         {
            coef = fixvalue;
            coef -= lastvalue;

            lowerbound_new = coef;
            lowerbound_new *= coef;
            lowerbound_new *= vcov_i;
            lowerbound_new += lowerbound;

            #if debug
            printf("x_%d = %d -> lower: %f (%d)\n", branchindex, (int) fixvalue , lowerbound_new, nodeindex );
            #endif
            if ( lowerbound_new < currentbestval )
            {
               NodeList.emplace_back( dpt1, nodeindex, lowerbound_new, fixedvalues );
               auto& nn = NodeList.back();
               nn.push_back_fixedvalue( fixvalue );

               auto& new_rsol = nn.geti_rsol();
               new_rsol.resize( m_dpt1 );

               for ( auto j = 0; j < m_dpt1 ; ++j )
               {
                  new_rsol[j] = rsol[j];
                  new_rsol[j] -= coef * vsco_i[j];
               }

               ++nodeindex;
            }
            else
               break;
            update_fix();
         }

         node.inc_dpt();
         node.set_lowerbound( init_lowerbound_new );
         node.push_back_fixedvalue( init_fixvalue );

         for ( auto j = 0; j < m_dpt1 ; ++j )
            rsol[j] -= init_coef * vsco_i[j];

         ++nodeindex;
      }
      else
         return true;

      assert( !fixedvalues.empty() && dpt+1 == (int)fixedvalues.size() );

      return false;
   };

   auto depthsearch = [&cond_vars,&run,&NodeList_sub,&nodemtxs,&bval,&vsco,&vcov,&vcov_inv,&bestsolval,&border,&solmtx,this,m,q11,borde_rate] ( const int id ) {

      unique_lock<mutex> lock( nodemtxs[id], defer_lock );

      while( 1 )
      {

         lock.lock();

         if ( NodeList_sub[id].empty() )
         {
            if ( run[id] == true )
               break;

            cond_vars[id].wait( lock, [&NodeList_sub,&run,id]
                  { return ( !NodeList_sub[id].empty() || run[id] ); } );

            if ( run[id] == true )
               break;
         }

         assert( !NodeList_sub[id].empty() );

         NODE node = move( NodeList_sub[id][0] );
         NodeList_sub[id].pop_front();

         lock.unlock();

         solmtx.lock_shared();
         auto upper = bval;
         solmtx.unlock_shared();

         #if depth_debug
         cout << "[Fix";
         cout << " best:" << (int)upper;
         cout << " index:" << (int) node.get_index() << "]-------------------------------" << endl;
         #endif

         // copy member variables of node
         const auto numfix = node.get_dpt();
         //const auto& fixedvalues = node.cget_fixedvalues();
         auto& rsol = node.geti_rsol();

         assert( !node.cget_fixedvalues().empty()
               && numfix == (int)node.cget_fixedvalues().size() ) ;
         assert( !rsol.empty() );
         assert( m - numfix == (int) rsol.size() );

         #if depth_debug
         cout << "fixedvalues: ";
         for ( auto val : node.cget_fixedvalues() )
            cout << val << " ";
         cout << endl;
         #endif


         auto dpt = 0;
         auto branchindex = m - 1 - numfix;
         auto branchindex_1 = branchindex - 1;
         auto m_numfix_1 = m - numfix - dpt - 1;
         const auto maxdpt = m - 1 - numfix;
         auto lastvalue = rsol[branchindex];
         int ub;
         int lb;
         double lower = node.get_lowerbound();
         double lower_new = 0.0;
         double upper_new;
         double buf;
         double coef = 0.0;
         constexpr double ep = 1.0e-10;
         vector<int> solvals( m - numfix );
         vector<DNode> DTree;

         assert( branchindex >= 0 );
         assert( branchindex < m );
         assert( ep < upper );
         assert( ep < lower );

         buf = upper;
         buf -= lower;
         if ( branchindex > 0 )
            buf *= vcov_inv[branchindex_1];
         else
            buf *= q11;
         buf = sqrt( buf );

         ub = floor( buf + lastvalue );
         lb = ceil( - buf  + lastvalue );

         solvals[0] = lb;
         DTree.reserve( m - numfix );
         DTree.emplace_back( ub, lb, lower, rsol );
         // Do not use rsol because of move
         //
         auto relaxsol = DTree[0].DNode_rsol();

         assert( m > numfix );

         for( ; ; )
         {
            assert( dpt >= 0 );
            assert( branchindex == m - 1 - numfix - dpt );
            assert( branchindex_1 == branchindex - 1 );
            assert( !DTree.empty() );
            assert( dpt + 1 == (int)DTree.size() );
            assert( m - numfix - dpt - 1 == m_numfix_1 );

            if ( dpt < maxdpt )
            {
               while( solvals[dpt] <= ub )
               {
                  coef = solvals[dpt];
                  coef -= lastvalue;

                  lower_new = coef;
                  lower_new *= coef;
                  lower_new *= vcov[branchindex_1];
                  lower_new += lower;

                  #if depth_debug
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
                  upper_new = solvals[dpt];
                  upper_new -= lastvalue;
                  upper_new *= upper_new;
                  upper_new *= q11;
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
                           solmtx.lock();
                           const auto& fixedvalues = node.cget_fixedvalues();
                           for ( int j = 0, k = (int)solvals.size() - 1; k != 0; ++j, --k )
                              bestsolval[j] = solvals[k];
                           for ( int j = (int)solvals.size(), k = (int)fixedvalues.size() - 1; k != 0; ++j, --k )
                              bestsolval[j] = fixedvalues[k];

                           upper = compute_objval( bestsolval );
                           bval = upper;
                           border = bval * borde_rate;
                           solmtx.unlock();
                           break;
                        }
                     }
                  }
                  ++solvals[dpt];
               }
            }

            if ( solvals[dpt] > ub )
            {
               assert( !DTree.empty() );
               DTree.pop_back();
               if ( DTree.empty() )
                  break;
               auto& dnode = DTree.back();
               --dpt;
               ++branchindex;
               ++branchindex_1;
               ++m_numfix_1;

               lower = dnode.DNode_lower();
               ub = dnode.DNode_ub();
               lb = dnode.DNode_lb();
               relaxsol = dnode.DNode_rsol();
               lastvalue = relaxsol[branchindex];
               ++solvals[dpt];
            }
            else
            {
               vector<double> new_relaxsol( m_numfix_1 );
               auto& vsco_i = vsco[branchindex_1];
               assert( !vsco_i.empty() );
               assert( m_numfix_1 == (int) vsco_i.size() );
               for ( auto j = 0; j < m_numfix_1 ; ++j )
               {
                  new_relaxsol[j] = relaxsol[j];
                  new_relaxsol[j] -= coef * vsco_i[j];
               }

               ++dpt;
               --branchindex;
               --branchindex_1;
               --m_numfix_1;

               lower = lower_new;
               lastvalue = new_relaxsol[branchindex];
               assert( (int)new_relaxsol.size() == branchindex + 1 );

               buf = upper;
               buf -= lower;
               if ( branchindex > 0 )
                  buf *= vcov_inv[branchindex_1];
               else
                  buf *= q11;
               buf = sqrt( buf );

               ub = floor( buf + lastvalue );
               lb = ceil( - buf  + lastvalue );

               solvals[dpt] = lb;

               DTree.emplace_back( ub, lb, lower, new_relaxsol );
               // Do not use new_relaxsol
               relaxsol = DTree.back().DNode_rsol();
            }
         }
      }
   };

   for ( auto i = 0; i < nsubsolvers; ++i )
      threads.push_back( thread( depthsearch, i ) );

   SVPSstartTime();

   while ( !NodeList.empty() )
   {
      assert( node == nullptr );
      assert( !NodeList.empty() );

      node = &NodeList.front();

      assert( node != nullptr );

      solmtx.lock_shared();
      if ( node->get_lowerbound() > border )
      {
         solmtx.unlock_shared();

         nodemtxs[push_id].lock();
         NodeList_sub[push_id].push_back( move( NodeList[0] ) );
         nodemtxs[push_id].unlock();

         cond_vars[push_id].notify_one();

         push_id = next_id( push_id );
         node = nullptr;
         NodeList.pop_front();
         continue;
      }
      solmtx.unlock_shared();

      while ( node != nullptr )
      {
         if ( branch( *node ) )
            node = nullptr;
      }

      NodeList.pop_front();
   }

   for ( auto i = 0; i < nsubsolvers; ++i )
   {
      assert( run[i] == false );
      nodemtxs[i].lock();
      run[i] = true;
      cond_vars[i].notify_one();
      nodemtxs[i].unlock();
   }

   for ( auto &th: threads )
      th.join();

   stopwatch.stop();

   status = SOLVED;

   bestval = bval;
   sq_bestval = sqrt( bestval );
   Appfac = sq_bestval / _Appfac;

   cout << "time: " << stopwatch.get_result() << "s" << endl;
   cout << "best value: " << bestval << endl;
   cout << "norm: " << sqrt( bestval ) << endl;
   cout << "AF: " << Appfac << endl;


   cout << "best solution:" << endl;
   for ( int i = 0; i < m; i++ )
   {
      if( bestsolval[i] != 0.0 )
         cout << " " << i << ": " << bestsolval[i] << endl;
   }

   return result;
}
