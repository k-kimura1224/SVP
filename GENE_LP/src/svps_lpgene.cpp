#include <string>
#include <string.h>
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
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "probdata.h"
#include "vector.h"


bool SVPsolver::gene_lpfile(
   int         formulationmode,
	const char*	filename	/*	input filename			  */
   )
{
   assert( formulationmode == 1 || formulationmode == 2 );

   int m = probdata.get_m();

   assert( m > 0 );

   // output bounds
   cout << "Bounds: " << endl;
   for ( int i = 0; i < m; i++ )
      cout << "x_" << i << ": [ " << lb[i] << ", " << ub[i] << "]" << endl;

   auto Q = probdata.get_Q();
   assert( Q != NULL );

   system("rm -f cplex.lp");

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

   if ( formulationmode == 1 )
   {
      // int vvs[m][2*ub[i]]; i = 1,...,m
      vector<vector<string>> vvs;
      vvs.resize( m );
      for( int i = 0; i < m; i++ )
      {
         auto lb_i = lb[i];
         auto ub_i = ub[i];
         for( int j = lb_i; j <= ub_i; j++ )
         {
            string y("y(");
            y += to_string(i);
            if ( j < 0 )
            {
               y += ",m";
               y += to_string(-j);
            }
            else
            {
               y += ",";
               y += to_string(j);
            }
            y += ")";
            vvs[i].push_back(y);
         }
      }

      for ( int i = 0; i < m; i++ )
      {
         auto lb_i = lb[i];
         auto ub_i = ub[i];

         lpfile << " xcons" << i << ": - x" << i;
         for( int j = lb_i, ct = 0; j <= ub_i; j++, ct++ )
         {
            if ( j < 0 )
               lpfile << " ";
            else
               lpfile << " + ";

            lpfile << j << " " << vvs[i][ct];
         }
         lpfile << " = 0" << endl;
      }

      for ( int i = 0; i < m; i++ )
      {
         lpfile << " ycons" << i << ":";
         for ( auto v : vvs[i] )
            lpfile << " +  " << v;

         lpfile << " = 1" << endl;
      }

      lpfile << " zerocons:";
      for ( auto v : vvs )
      {
         lpfile << " + " << v[int(v.size()/2)];
      }
      lpfile << " <= " << m - 1 << endl;

      lpfile << "Bounds" << endl;

      for ( int i = 0; i < m; i++ )
         lpfile << " " << (int)lb[i] << " <= x" << i << " <= " << (int)ub[i] << endl;

      lpfile << "Generals" << endl;

      for ( int i = 0; i < m; i++ )
         lpfile << " x" << i;
      lpfile << endl;

      lpfile << "Binary" << endl;

      for ( auto v : vvs )
         for( auto vv : v )
            lpfile << " " << vv;
      lpfile << endl;

      lpfile << "End" << endl;
   }
   else if ( formulationmode == 2 )
   {
      vector<string> x;
      vector<string> xp;
      vector<string> xm;
      vector<string> zp;
      vector<string> zm;
      vector<vector<string>> yp;
      vector<vector<string>> ym;
      yp.resize( m );
      ym.resize( m );

      for ( int i = 0; i < m; i++ )
      {
         string s_x("x");
         string s_xp("xp");
         string s_xm("xm");
         string s_zp("zp");
         string s_zm("zm");

         s_x += to_string(i);
         s_xp += to_string(i);
         s_xm += to_string(i);
         s_zp += to_string(i);
         s_zm += to_string(i);

         x.push_back(s_x);
         xp.push_back(s_xp);
         xm.push_back(s_xm);
         zp.push_back(s_zp);
         zm.push_back(s_zm);
      }

      for ( int i = 0; i < m; i++ )
      {
         assert( ub[i] == - lb[i] );
         auto ub_i = ub[i];
         auto t = 1;
         int k;
         for ( k = 0; k <= 10; k++ )
         {
            if ( ub_i <= t ) break;
            t += pow( 2, k + 1);
         }

         for ( int j = 0; j <= k; j++ )
         {
            string s_yp("yp(");
            s_yp += to_string(i);
            s_yp += ",";
            s_yp += to_string(j);
            s_yp += ")";
            yp[i].push_back(s_yp);

            string s_ym("ym(");
            s_ym += to_string(i);
            s_ym += ",";
            s_ym += to_string(j);
            s_ym += ")";
            ym[i].push_back(s_ym);
         }
      }

      for ( int i = 0; i < m; i++ )
      {
         lpfile << " xpm" << i << ": ";
         lpfile << x[i] << " - " << xp[i] << " + " << xm[i] << " = 0";
         lpfile << endl;

         lpfile << " xp" << i << ": - " << xp[i];
         for( int j = 0; j < (int)yp[i].size(); j++)
         {
            lpfile << " + ";
            lpfile << (int)pow(2,j) << " " << yp[i][j];
         }
         lpfile << " = 0" << endl;

         lpfile << " xm" << i << ": - " << xm[i];
         for( int j = 0; j < (int)ym[i].size(); j++)
         {
            lpfile << " + ";
            lpfile << (int)pow(2,j) << " " << ym[i][j];
         }
         lpfile << " = 0" << endl;
      }

      for ( int i = 0; i < m; i++ )
      {
         assert( ub[i] >= 0 && lb[i] <= 0 );

         lpfile << " zp" << i << ": ";
         lpfile << (int)ub[i] << " " << zp[i] << " - " << xp[i] << " >= 0";
         lpfile << endl;

         lpfile << " zm" << i << ": ";
         lpfile << (int)-lb[i] << " " << zm[i] << " - " << xm[i] << " >= 0";
         lpfile << endl;

         lpfile << " zpm" << i << ": ";
         lpfile << zp[i] << " + " << zm[i] << " = 1" << endl;
      }

      lpfile << " zerocons: ";
      for ( int i = 0; i < m; i++ )
      {
         lpfile << " + " << xp[i] << " + " << xm[i];
      }
      lpfile << " >= 1" << endl;

      lpfile << "Bounds" << endl;

      for ( int i = 0; i < m; i++ )
      {
         lpfile << " " << (int)lb[i] << " <= " << x[i] << " <= " << (int)ub[i] << endl;
         lpfile << " " << 0 << " <= " << xp[i] << " <= " << (int)ub[i] << endl;
         lpfile << " " << 0 << " <= " << xm[i] << " <= " << - (int)lb[i] << endl;
      }

      lpfile << "Generals" << endl;

      for ( int i = 0; i < m; i++ )
      {
         lpfile << " x" << i;
         lpfile << " xp" << i;
         lpfile << " xm" << i;
      }
      lpfile << endl;

      lpfile << "Binary" << endl;

      for ( auto v : yp )
         for( auto vv : v )
            lpfile << " " << vv;

      for ( auto v : ym )
         for( auto vv : v )
            lpfile << " " << vv;

      for ( int i = 0; i < m; i++ )
      {
         lpfile << " " << zp[i];
         lpfile << " " << zm[i];
      }

      lpfile << endl;
      lpfile << "End" << endl;
   }
   else
   {
      assert(0);
   }

   return true;
}

