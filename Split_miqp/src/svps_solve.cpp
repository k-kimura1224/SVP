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

int SVPsolver::selection_k(
      int   i,
      int   selection
      )
{
   auto k = 0;
   switch ( selection )
   {
      case 0:
      {
         k = i;
         break;
      }

      case 1:
      {
         k = probdata.get_m() - 1 - i;
         break;
      }

      default:
      {
         assert(0);
      }
   }

   return k;
}

double SVPsolver::compute_miqp_cplex(
      int   k
      )
{
   auto Q = probdata.get_Q();
   auto m = probdata.get_m();
   assert( Q != NULL );

   ofstream lpfile;
   lpfile.open( "cplex.lp", std::ios::out );

   lpfile << "Minimize" << endl;
   lpfile << " obj: [";

   {
      int ct=0;
      int   n = (( 1 + m ) * m )/2;
      int*  coef = new int[n];
      int*  var_x1 = new int[n];
      int*  var_x2 = new int[n];
      for(int i=0; i<m; i++){
         for(int j=i; j<m; j++){
            var_x1[ct] = i;
            var_x2[ct] = j;
            if( i==j )  coef[ct] = (int)Q[(i*m)+j];
            else        coef[ct] = 2 * (int)Q[(i*m)+j];
            ct++;
         }
      }

      for(int i=0; i<n; i++){
         if( coef[i] >= 0 ){
            lpfile << " +";
         }else{
            lpfile << " ";
         }

         lpfile << 2*coef[i];

         if( var_x1[i] == var_x2[i] ){
            lpfile << " x" << var_x1[i] << "^2";
         }else{
            lpfile << " x" << var_x1[i] << " * x" << var_x2[i];
         }
      }

      delete[] coef;
      delete[] var_x1;
      delete[] var_x2;
   }

   lpfile << "]/2" << endl;

   lpfile << "Subject to" << endl;
   lpfile << "Bounds" << endl;

   for(int i=0; i<m; i++)
   {
      if( i == k )
      {
         assert( 0.99 <= ub[i] );
         lpfile << " " << 1 << " <= x" << i << " <= " << (int)ub[i] << endl;
      }
      else
         lpfile << " " << (int)lb[i] << " <= x" << i << " <= " << (int)ub[i] << endl;
   }

   lpfile << "Generals" << endl;

   for(int i=0; i<m; i++)
      lpfile << " x" << i;

   lpfile << endl;
   lpfile << "End" << endl;


   //if( 1 && system("~/Applications/IBM/ILOG/CPLEX_Studio1262/cplex/bin/x86-64_osx/cplex < run.txt") == -1 )
   if( system("./cplex.sh") == -1 )
   {
      assert(0);
      exit(0);
   }

   cout << endl;
   cout << endl;

   FILE  *file;
   int   status;
   char  s[50];

   /* open file */
   file = fopen( "cplex.sol", "r");
   if( file==NULL ){
      printf("Could not open file <%s>.\n", "cplex.sol");
      exit(0);
   }

   for(int i=0; i<7; i++){
      /* skip 6 lines */
      if( fgets( s, 50, file)==NULL ){
         printf("Error reading file <%s>.\n", "cplex.sol");
         exit(0);
      }
   }

   /* read */
   double val;
   status = fscanf( file, "   objectiveValue=\"%lf\"", &val);
   if( !status ){
      cout << "Reading error" << endl;
      exit(0);
   }

   if ( bestval > val )
   {
      for(int i=7; i<27; i++){
         /* skip 6 lines */
         if( fgets( s, 50, file)==NULL ){
            printf("Error reading file <%s>.\n", "cplex.sol");
            exit(0);
         }
      }

      char buf[50];
      int j = -1;
      vector<double> solvals( m, 0.0 );
      for ( int i = 0; i < m; i++ )
      {
         status = fscanf( file, "  %s %s index=\"%d\" value=\"%lf\"/>", buf, buf, &j, &solvals[i] );
         solvals[i] = round( solvals[i] );
         assert( i == j );
         if( !status ){
            cout << "Reading error" << endl;
            exit(0);
         }
      }

      SOLUTION sol;
      sol.set_sol( m, &solvals[0], val);
      bestsol = sol;

      bestval = val;
   }

   fclose( file);

   system("rm -f cplex.lp");
   system("rm -f cplex.sol");

   return val;
}


bool SVPsolver::solve_cplex(
   int      tlimit,
   int      selection
   )
{
   if( tlimit <= 0 ){
      return false;
   }

   stopwatch.set_timelimit( tlimit );
   stopwatch.start();

   int m = probdata.get_m();

   assert( m > 0 );

   vector<bool> list( m, true );

   // output bounds
   cout << "Bounds: " << endl;
   for ( int i = 0; i < m; i++ )
   {
      cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl;
   }

   ofstream txtfile;
   txtfile.open( "run.txt", std::ios::out);
   txtfile << "read cplex.lp" << endl;
   txtfile << "opt" << endl;
   txtfile << "write cplex.sol" << endl;
   txtfile << "q" << endl;
   system("rm -f cplex.sol");
   system("rm -f clone*.log");
   system("rm -f cplex.log");

   int k;
   double optval_Q;
   double lowerbound_P;

   for ( int i = 0; i < m - 1; i++ )
   {
      cout << "-[" << i << "]--------------------------------------------" << endl;
      k = selection_k( i, selection );
      assert( list[k] == true );
      cout << "selection: " << k << endl;

      optval_Q = compute_miqp_cplex( k );
      cout << "optval(Q): " << optval_Q;

      if ( i < m - 2 )
      {
	      for ( int j = 0; j < m; j++ )
			   sch.set_z_i( j, list[j]);

	      assert( sch.get_n() > 0 );

         sch.compute_GS();
         lowerbound_P = sch.get_min();
         cout << ", lowerbound(P): " << lowerbound_P << endl;
         if ( lowerbound_P > bestval )
            break;
      }

      list[k] = false;
      ub[k] = 0.0;
      lb[k] = 0.0;

      if( stopwatch.check_time() == false ) break;

   }

   system("rm -f run.txt");
   system("rm -f cplex.txt");
   system("rm -f clone.txt");

   if( stopwatch.check_time() == false )
      return false;

   stopwatch.stop();
   return true;
}

bool SVPsolver::solve(
   bool     para,
   bool     s_zero,
   double   s_GLB,
   int      tlimit
   )
{
   if( tlimit <= 0 ){
      return false;
   }

   stopwatch.set_timelimit( tlimit );
   stopwatch.start();

   int m = probdata.get_m();
   assert( m > 0 );

   // output bounds
   if( para == false ){
      cout << "Bounds: " << endl;
      for(int i=0; i<m; i++){
         cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl;
      }
   }

   // generate a root node
   NODE  root;
   int   index=0;

   double   *init_warm =NULL;
   if( para == false ){
      init_warm = bestsol.get_solval();
   }else{
      init_warm = new double[m];
      for(int i=0; i<m; i++){
         init_warm[i] = (ub[i] + lb[i])/2.0;
      }
   }
   root.set_vals( m, ub, lb, init_warm, s_GLB, 0, s_zero, index);

   for(int i=0; i<m; i++){
      if( lb[i] - ub[i] == 0 && lb[i] != 0.0 ){
         if( root.alloc_sumfixed() == true ){
            root.set_sumfixed( lb[i], probdata.get_bvec(i) );
         }else{
            root.add_sumfixed( lb[i], probdata.get_bvec(i) );
         }
      }
   }

   NodeList.push_back( root );
   listsize++;
   index++;

   if( para == true ){
      delete[] init_warm;
   }
   init_warm = NULL;

   // generate oa_cpool
   if( CUT_OA == true ){
      gene_OAcuts( ub, lb, probdata.get_Q(), bestval);
   }

   int         sel;
   RelaxResult r;
   QP_time  = 0.0;
   __time = 0.0;
   __start = clock();
   int   disp = index;
   int   cutoff=0;
   while(1){
      assert( (int)NodeList.size() > 0 );
      assert( listsize > 0 );
      assert( (int)NodeList.size() == listsize );

      // select a node from the list
      sel = select_node(index, disp);

      assert( sel >= 0 );
      assert( sel < listsize );

      // solve a relaxation problem
      r = solve_relaxation(sel);

      if( r == INFEASIBLE || r == GETINTEGER ){
         cutoff++;
      }

      // run heuristics
      if( r == FEASIBLE && HEUR_APP < Appfac ){
         heur( sel, para);
      }
      // output
      if( (index-1) % 1000 == 0 && disp < index ){
         list<NODE>::iterator it = NodeList.begin();
         double min_lb = it->get_lowerbound();
         for(int i=1; i<listsize; i++){
            ++it;
            if( min_lb > it->get_lowerbound() ){
               min_lb = it->get_lowerbound();
            }
         }
         GLB = min_lb;
         if( para == true ) cout << "t" << omp_get_thread_num() << ":";
         disp_log(sel, r, index, cutoff);
         disp = index;
         cutoff = 0;
      }

      // branch
      if( r == UPDATE || r == FEASIBLE ){
         branch( sel, index);
         index += 2;
      }

      // remove
      clock_t  start = clock();
      list<NODE>::iterator it = NodeList.begin();
      advance( it, sel);
      //NodeList.erase( NodeList.begin() + sel);
      NodeList.erase( it );
      listsize--;

      clock_t  end = clock();
      __time += (double)(end-start)/CLOCKS_PER_SEC;
      // break
      assert( (int)NodeList.size() == listsize );
      if( stopwatch.check_time() == false ) break;
      if( listsize == 0 ){
         if( !para ) cout << "End" << endl;
         break;
      }


   }

   nnode = (unsigned long int)index -1;


   if( stopwatch.check_time() == false ){
      return false;
   }

   stopwatch.stop();
   return true;
}
