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
#include <list>
#include <iomanip>
#include <numeric>
//#include <omp.h>

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "stopwatch.h"
#include "node.h"
#include "nodelist.h"
#include "Schmidt_manager.h"

using namespace std;

#define debug_class  0
#define debug  0

SVPsolver::SVPsolver(){ // default constructor
#if debug_class
   cout << "SVPsolver: default constructor" << endl;
#endif
   bestval = -1.0;
   GLB = 0.0;

   //type = LIST;
   //type = TWO_DEQUE;
   type = TEN_DEQUE;

   push_back = nullptr;
   move_back = nullptr;
   nodeselection = nullptr;
   cut_off = nullptr;
   get_GLB = nullptr;
   check_size = nullptr;
   setup_para_selection = nullptr;
   setup_parapush_selection = nullptr;
   setup_parapop_selection = nullptr;
   para_selection = nullptr;
   pop_front = nullptr;
   getSubsize = nullptr;

   index = 0;
   nnode = 0;
   nthreads = -1;

   _Appfac = 0.0;
   Appfac = 0.0;

   quiet = false;

   norm = nullptr;
   TIMELIMIT = 86400;
   LEFTNODELIMIT = 2000000000;
   NODELIMIT = 2000000000;
   MEMORY = 8;
   CUTMODE = false;
   ENUM = false;
   IMPRELAX = false;

   epsilon = 1.0e-12;

   subsolver = false;
   status = SETUP;

   SVPSbranch = &SVPsolver::SVPSbranchDefault;
   SVPSrelax = &SVPsolver::SVPSrelaxDefault;
}

SVPsolver::SVPsolver( const SVPsolver &source )
{  // copy constructor
#if debug_class
   cout << "SVPsolver: copy constructor" << endl;
#endif

   probdata = source.probdata;
   GLB = source.GLB;
   bestval = source.bestval;
   sq_bestval = source.sq_bestval;
   bestsol = source.bestsol;
   pool = source.pool;
   _Appfac = source._Appfac;
   Appfac = source.Appfac;
   type = source.type;
   nodelist = source.nodelist;
   push_back = source.push_back;
   move_back = source.move_back;
   nodeselection = source.nodeselection;
   cut_off = source.cut_off;
   get_GLB = source.get_GLB;
   check_size = source.check_size;
   setup_para_selection = source.setup_para_selection;
   setup_parapush_selection = source.setup_parapush_selection;
   setup_parapop_selection = source.setup_parapop_selection;
   para_selection = source.para_selection;
   pop_front = source.pop_front;
   getSubsize = source.getSubsize;
   index = source.index;
   stopwatch = source.stopwatch;
   testwatch = source.testwatch;
   nnode = source.nnode;
   quiet = source.quiet;
   nthreads = source.nthreads;
   TIMELIMIT = source.TIMELIMIT;
   LEFTNODELIMIT = source.LEFTNODELIMIT;
   NODELIMIT = source.NODELIMIT;
   MEMORY = source.MEMORY;
   CUTMODE = source.CUTMODE;
   ENUM = source.ENUM;
   IMPRELAX = source.IMPRELAX;
   epsilon = source.epsilon;
   subsolver = source.subsolver;
   status = source.status;
   bounds = source.bounds;
   //oa_cuts = source.oa_cuts;
   SNOVs = source.SNOVs;
   CMGSO = source.CMGSO;
   VsCB = source.VsCB;
   LVsCB = source.LVsCB;
   VsCO = source.VsCO;
   VCOV = source.VCOV;
   SVPSbranch = source.SVPSbranch;
   SVPSrelax = source.SVPSrelax;

   if ( probdata.get_m() > 0 )
   {
      int m = probdata.get_m();
      if( source.norm != nullptr ){
         norm = new double[m];
         Copy_vec( source.norm, norm, m);
      }else{
         norm = nullptr;
      }
   }else{
      norm = nullptr;
   }
}

// assignment operator
SVPsolver& SVPsolver::operator=( const SVPsolver& source )
{
#if debug_class
   cout << "SVPsolver: assignment operator" << endl;
#endif

   if( this != &source ){
      probdata = source.probdata;
      GLB = source.GLB;
      bestval = source.bestval;
      sq_bestval = source.sq_bestval;
      bestsol = source.bestsol;
      pool = source.pool;
      _Appfac = source._Appfac;
      Appfac = source.Appfac;
      type = source.type;
      push_back = source.push_back;
      move_back = source.move_back;
      nodeselection = source.nodeselection;
      cut_off = source.cut_off;
      get_GLB = source.get_GLB;
      check_size = source.check_size;
      setup_para_selection = source.setup_para_selection;
      setup_parapush_selection = source.setup_parapush_selection;
      setup_parapop_selection = source.setup_parapop_selection;
      para_selection = source.para_selection;
      pop_front = source.pop_front;
      getSubsize = source.getSubsize;
      nodelist = source.nodelist;
      index = source.index;
      stopwatch = source.stopwatch;
      testwatch = source.testwatch;
      nnode = source.nnode;
      quiet = source.quiet;
      nthreads = source.nthreads;
      TIMELIMIT = source.TIMELIMIT;
      LEFTNODELIMIT = source.LEFTNODELIMIT;
      NODELIMIT = source.NODELIMIT;
      MEMORY = source.MEMORY;
      CUTMODE = source.CUTMODE;
      ENUM = source.ENUM;
      IMPRELAX = source.IMPRELAX;
      epsilon = source.epsilon;
      subsolver = source.subsolver;
      status = source.status;
      bounds = source.bounds;
      //oa_cuts = source.oa_cuts;
      SNOVs = source.SNOVs;
      CMGSO = source.CMGSO;
      VsCB = source.VsCB;
      LVsCB = source.LVsCB;
      VsCO = source.VsCO;
      VCOV = source.VCOV;
      SVPSbranch = source.SVPSbranch;
      SVPSrelax = source.SVPSrelax;

      if( probdata.get_m() > 0 ){
         int m = probdata.get_m();

         if( source.norm != nullptr ){
            delete[] norm;
            norm = new double[m];
            Copy_vec( source.norm, norm, m);
         }else{
            norm = nullptr;
         }
      }else{
         norm = nullptr;
      }
   }

   return *this;
}

// destructor
SVPsolver::~SVPsolver()
{
#if debug_class
   cout << "SVPsolver: destructor" << endl;
#endif
   delete[] norm;
   norm = nullptr;

   for ( auto& b : bounds )
   {
      b.clear();
      b.shrink_to_fit();
   }
   bounds.clear();
   bounds.shrink_to_fit();

   assert( bounds.empty() );

   //for ( auto& c : oa_cuts )
   //{
   //   c.clear();
   //   c.shrink_to_fit();
   //}
   //oa_cuts.clear();
   //oa_cuts.shrink_to_fit();

   //assert( oa_cuts.empty() );
   SNOVs.clear();
   SNOVs.shrink_to_fit();
   for ( auto& c : CMGSO )
   {
      c.clear();
      c.shrink_to_fit();
   }
   CMGSO.clear();
   CMGSO.shrink_to_fit();

   for ( auto& mat : VsCB )
   {
      for ( auto& vec : mat )
      {
         vec.clear();
         vec.shrink_to_fit();
      }
      mat.clear();
      mat.shrink_to_fit();
   }
   VsCB.clear();
   VsCB.shrink_to_fit();
   for ( auto& l : LVsCB )
   {
      l.clear();
      l.shrink_to_fit();
   }
   LVsCB.clear();
   LVsCB.shrink_to_fit();

   for ( auto& v : VsCO )
   {
      v.clear();
      v.shrink_to_fit();
   }
   VsCO.clear();
   VsCO.shrink_to_fit();

   VCOV.clear();
   VCOV.shrink_to_fit();

   push_back = nullptr;
   move_back = nullptr;
   nodeselection = nullptr;
   cut_off = nullptr;
   get_GLB = nullptr;
   check_size = nullptr;
   setup_para_selection = nullptr;
   setup_parapush_selection = nullptr;
   setup_parapop_selection = nullptr;
   para_selection = nullptr;
   pop_front = nullptr;
   getSubsize = nullptr;

   SVPSbranch = nullptr;
   SVPSrelax = nullptr;
}

void SVPsolver::SVPSsetup(
      const int      s_m,
      const double*  s_B_,
      const int      s_nthreads,
      const int      s_timelimit,
      const int      s_memory,
      const bool     s_quiet,
      const bool     w_subsolver,   // wheter this is subsolver
      const bool     w_heur,        // wheter heuristic is executed
      const bool     w_gs,          // wheter Gram-Schmidt is set
      const bool     w_app,         // wheter Appfac is initialized
      const bool     w_gn,          // wheter node is generated
      const bool     w_nl           // wheter nodelist is set
      )
{
   auto m = s_m;
   auto B_ = s_B_;

   // probdata
   SVPScreateProbdata( m, B_ );

   // nthreads
   assert( s_nthreads > 0 );
   nthreads = s_nthreads;

   // TIMELIMIT
   assert( s_timelimit > -1 );
   TIMELIMIT = s_timelimit;

   // memory
   assert( s_timelimit > 0 );
   MEMORY = s_memory;

   // quiet
   quiet = s_quiet;

   // pool
   pool.alloc( 5 );

   // subsolver
   subsolver = w_subsolver;

   if ( w_gs )
   {
      assert( w_subsolver == false );
      SVPSsetGS();
   }

   // Appfac
   if ( w_app )
   {
      assert( w_subsolver == false );

      _Appfac = tgamma( ((double)m/2.0) + 1 );
      _Appfac = pow( _Appfac, 1.0/(double)m );

      //double absdet = fabs( determinant( B_, m) );
      //absdet = pow( absdet, 1.0/(double)m );
      //_Appfac *= absdet;
      _Appfac /= sqrt( M_PI );

      assert( !SNOVs.empty() );
      assert( m == (int) SNOVs.size() );
      for ( int j = m - 1; j >= 0; --j )
         _Appfac *= pow( sqrt( SNOVs[j] ), 1.0/(double)m);
   }

   // heuristic
   if ( w_heur )
   {
      assert( w_subsolver == false );
      SVPSheurFindMinColumn();
   }

   // nodelist
   if ( w_nl )
   {
      assert( w_subsolver == false );
      assert( w_heur == true );
      SVPSsetupNodelist();
   }

   // generate nodes
   if ( w_gn )
   {
      assert( w_subsolver == false );
      SVPSgenerateRootNode();
   }
}


void SVPsolver::SVPScreateProbdata(
   const int      m,
   const double   *B_
   )
{
   assert( m > 0 );
   assert( B_ != nullptr );

   const int mm = m * m;

   double   *B;      // [m*m],
   double   *Q;      // [m*m],

   B  =  new double[mm];
   Q  =  new double[mm];

   TraMat( m, m, B_, B);

   Gen_ZeroVec( mm, Q);
   Com_mat_AtA( B_, m, m, Q);

   probdata.set_data( m, B, B_, Q );

   delete[] B;
   delete[] Q;

}

void SVPsolver::SVPSgenerateNodes()
{
   assert(0);
   //assert( subsolver == false );
   //assert( (int)bounds.size() == 1 );

   //// copy
   //const auto& pd = probdata;
   //const auto m = pd.get_m();
   //const auto B_ = pd.get_B_();
   //auto& vars_globalbounds = bounds;
   //const auto best = bestval;
   //auto& nodeindex = index;
   //auto& NL = nodelist;
   //const auto ep = epsilon;

   //assert( m > 0 );
   //assert( B_ != nullptr );
   //assert( !vars_globalbounds.empty() );
   //assert( best > 0 );
   //assert( nodeindex == 0 );
   //assert( push_back != nullptr );

   //// node
   //double* init_warm = new double [m];
   //int dpt = 0;
   //int type = 0;
   //constexpr int branchvalue_one = 1;
   //constexpr int branchvalue_zero = 0;

   ////  Schmidt_manager
   //SCHMIDT_M sch;
   //double schmin;

   ////
   //auto m1 = m - 1;
   //auto m2 = m - 2;

   //Gen_ZeroVec( m, init_warm );
   //sch.setup( m, B_ );

   //assert( sch.get_n() > 0 );

   //++index;
   //++dpt;

   //for ( int i = 0; i < m1; ++i )
   //{
   //   auto& vgb = vars_globalbounds[type];
   //   assert( (int) vgb.size() == m );

   //   if ( !(vgb[i]) )
   //      continue;

   //   for ( int j = 0; j < m; ++j )
   //      sch.set_z_i( j, (bool) vgb[j] );

   //   // compute lower bound
   //   sch.compute_GS();
   //   schmin = sch.get_min();

   //   // test
   //   if ( schmin > best )
   //      break;

   //   // branch 1 =< x_i and x_i = 0
   //   NODE node;
   //   init_warm[i] = 1.0;
   //   node.set_vals( m, init_warm, schmin, dpt, nodeindex, type );
   //   init_warm[i] = 0.0;
   //   ++dpt;
   //   nodeindex += 2;
   //   ++type;
   //   for ( auto j = 0; j < i; ++j )
   //      node.set_branchinfo( j, branchvalue_zero, 'e' );
   //   node.set_branchinfo( i, branchvalue_one, 'l' );

   //   if( Equal( branchvalue_one, vgb[i], ep ) )
   //   {
   //      if( node.alloc_sumfixed() )
   //         node.set_sumfixed( branchvalue_one, pd.get_bvec( i ) );
   //      else
   //         node.add_sumfixed( branchvalue_one, pd.get_bvec( i ) );
   //   }

   //   (NL.*move_back)( node );

   //   // compute bounds of variables
   //   if ( i < m2 )
   //   {
   //      assert( type == (int) vars_globalbounds.size() );

   //      SVPStightenBounds( i );

   //      assert( type + 1 == (int) vars_globalbounds.size() );
   //   }
   //}

   //delete[] init_warm;

}

void SVPsolver::SVPSgenerateRootNode()
{
   assert( subsolver == false );

   // copy
//   const auto& pd = probdata;
//   const auto m = pd.get_m();
//   const auto B_ = pd.get_B_();
//   auto& vars_globalbounds = bounds;
//   const auto best = bestval;
   auto& nodeindex = index;
   auto& NL = nodelist;
//   const auto ep = epsilon;
//

//   assert( m > 0 );
//   assert( B_ != nullptr );
//   assert( !vars_globalbounds.empty() );
//   assert( best > 0 );
   assert( nodeindex == 0 );
   assert( push_back != nullptr );

   vector<int> empty_v;
   NODE node( 0, nodeindex, 0.0, empty_v );
   if ( !ENUM )
      node.geti_rsol().assign( probdata.get_m(), 0 );
   (NL.*move_back)( node );
   ++nodeindex;
}

void SVPsolver::SVPSmoveNode(
      NODE&    movenode
      )
{

   movenode.set_index( index );
   (nodelist.*move_back)( movenode );

   index++;
   // Do not use movenode after here
}

NODE& SVPsolver::SVPSgetNode_para_selection( const int setup )
{
   assert( (nodelist.*getSubsize)(setup) > 0 );
   return (nodelist.*para_selection)(setup);
}

void SVPsolver::SVPSpopNode( const int setup )
{
   assert( nodelist.getListsize() > 0 );
   assert( (nodelist.*check_size)() );

   (nodelist.*pop_front)(setup);
}

string SVPsolver::SVPSgetStringStatus() const
{
   switch ( status )
   {
      case TIMEOVER:
         return "TIMEOVER";
      case FULL_OF_LEFTNODES:
         return "FULL_OF_LEFTNODES";
      case FULL_OF_NODES:
         return "FULL_OF_NODES";
      default:
         return "";
   }
}

void SVPsolver::disp_bestsol()
{
   int   m;
   double *val = nullptr;

   m = probdata.get_m();
   val = bestsol.get_solval();

   assert( val != nullptr );

   int nodelistsize = nodelist.getListsize();

   cout << endl;
   if( nodelistsize == 0 )
   {
      cout << "SVPSOLVER found an optimal solution" << endl;
   }
   else if ( nodelistsize >= LEFTNODELIMIT )
   {
      cout << "-- LEFT_NODE_LIMIT --" << endl;
   }
   else if ( index >= NODELIMIT )
   {
      cout << "-- NODE_LIMIT --" << endl;
   }
   else
   {
      cout << "-- TIMEOVER --" << endl;
   }
   cout << endl;

   cout << "time: " << stopwatch.get_result() << "s" << endl;
   cout << "best value: " << bestval << endl;
   cout << "norm: " << sqrt( bestval ) << endl;
   cout << "AF: " << Appfac << endl;

   if( nodelistsize != 0 ){
      cout << "lower bound: " << GLB << endl;
      cout << "gap: " << 100*(bestval - GLB)/bestval << "%" << endl;
   }
   cout << "Nodes: " << nnode << endl;

   cout << "best solution:" << endl;
   for(int i=0; i<m; i++){
      if( val[i] != 0.0 ){
         cout << " " << i << ": " << val[i] << endl;
      }
   }



}

void  SVPsolver::SVPScomputeBounds()
{
   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto B = pd.get_B();

   double *e;
   double *a;

   int r;
   auto  sqrt_bestval = sqrt( bestval );

   assert( bounds.empty() );
   assert( bestval > 0.0 );
   assert( m > 0 );

   e = new double[m];
   a = new double[m];

   Gen_ZeroVec( m, e );

   bounds.reserve( m );
   bounds.resize( 1 );
   bounds[0].resize( m );

   for ( int i = 0; i < m; ++i )
   {
      assert( e[i] == 0.0 );
      e[i] = 1.0;

      r = Com_LS_dgesv( B, e, m, a );
      //r = Com_LS_dposv( B, e, m, a );
      if ( r != 0 )
      {
         cout << "error: r = " << r << endl;
         exit(-1);
      }

      bounds[0][i] = floor( sqrt_bestval * Com_nrm( a, m ) );

      e[i] = 0.0;
   }

   //if ( CUTMODE )
   //{
   //   const auto Q = pd.get_Q();
   //   double* eigenvectors = new double[m*m]; // unit
   //   double* eigenvalues = new double[m];
   //   int ct = 0;
   //   double coef;
   //   double U;
   //   double L;

   //   Gen_ZeroVec( m, eigenvalues );

   //   r = Com_eigen( Q, m, eigenvalues, eigenvectors );
   //   if ( r != 0 )
   //   {
   //      cout << "error(cut): r = " << r << endl;
   //      exit(-1);
   //   }

   //   oa_cuts.reserve( m );
   //   oa_cuts.resize( 1 );
   //   oa_cuts[0].reserve( m );

   //   ct = 0;
   //   for ( auto i = 0; i < m; ++i )
   //   {
   //      coef = sqrt( eigenvalues[i] * bestval );
   //      for ( auto j = 0; j < m; ++j )
   //      {
   //         e[j] = round( coef * eigenvectors[ct++] );
   //      }

   //      assert( ct == (i+1) * m );
   //      r = Com_LS_dgesv( B, e, m, a );
   //      if ( r != 0 )
   //      {
   //         cout << "error: r = " << r << endl;
   //         exit(-1);
   //      }

   //      U = floor( sqrt_bestval * Com_nrm( a, m ) );
   //      L = - U;

   //      oa_cuts[0].emplace_back( m, e, &L, &U );
   //      //for ( auto j = 0; j < m; ++j )
   //      //{
   //      //   e[j] = eigenvectors[ct++];
   //      //}
   //      //assert( ct == (i+1) * m );
   //      //U = sqrt( bestval/eigenvalues[i] );
   //      //L = - U;
   //      //oa_cuts[0].emplace_back( m, e, &L, &U );
   //   }

   //   delete[] eigenvectors;
   //   delete[] eigenvalues;
   //}

   delete[] e;
   delete[] a;
}


void  SVPsolver::SVPStightenBounds(
      const int      memo
      )
{
   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto B_ = pd.get_B_();
   const auto Q = pd.get_Q();

   auto dim = m - ( memo + 1 );
   auto mdim = m * dim;
   auto dimdim = dim * dim;

   double*  subB = new double[mdim];
   double*  subBB = new double[dimdim];
   double*  e = new double[dim];
   double*  buf_a = new double[dim];
   double*  buf_b = new double[m];

   int ct ;

   Gen_ZeroVec( dim, e );

   auto constant = m * ( memo + 1 );
   for ( auto i = 0; i < mdim; i++ )
      subB[i] = B_[constant + i];

   ct = 0;
   for ( int i = memo + 1, im; i < m; i++ )
   {
      im = i * m;
      for ( auto j = memo + 1; j < m; j++ )
      {
         subBB[ct] = Q[im + j];
         assert( subBB[ct] == Com_dot( &subB[(i-(memo+1))*m], &subB[(j-(memo+1))*m], m));
         ct++;
      }
   }

   int r;
   const auto sqrt_bestval = sqrt( bestval );

   constant = memo + 1;

   vector<int> tightendbounds;
   tightendbounds.reserve( m );

   for ( auto i = 0; i <= memo; i++ )
      tightendbounds.push_back( 0 );

   for ( auto i = 0; i < dim; i++ )
   {
      assert( e[i] == 0.0 );
      e[i] = 1.0;

      //r = Com_LS_dgesv( subBB, e, dim, buf_a );
      r = Com_LS_dposv( subBB, e, dim, buf_a );
      if ( r != 0 )
      {
         cout << "error(tighten): r = " << r << endl;
         exit(-1);
      }

      Com_mat_Ax( subB, m, dim, buf_a, buf_b );

      tightendbounds.push_back( floor( sqrt_bestval * Com_nrm( buf_b, m ) ) );

      e[i] = 0.0;

   }

   assert( (int)tightendbounds.size() == m );

   bounds.push_back( move ( tightendbounds ) );

   delete[] subB;
   delete[] subBB;
   delete[] e;
   delete[] buf_a;
   delete[] buf_b;
}

void  SVPsolver::SVPScomputeVsCB()
{
   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto B_ = pd.get_B_();
   const auto Q = pd.get_Q();
   auto& vscb = VsCB;
   auto& lvscb = LVsCB;

   assert( m > 0 );
   assert( B_ != nullptr );
   assert( Q != nullptr );

   auto dim = 1;
   int ct;

   vector<double> subBB( m * m );   // dim * dim
   vector<double> e( m, 0 );        // dim
   vector<double> buf_a( m );       // dim
   vector<double> buf_b( m );       // m

   vscb.resize( m );
   for ( auto i = 0; i < m; ++i )
   {
      vscb[i].resize( i + 1 );
      for ( auto& vec : vscb[i] )
         vec.resize( m );
   }

   lvscb.resize( m );
   for ( auto i = 0; i < m; ++i )
      lvscb[i].resize( i + 1 );

   while ( dim <= m )
   {
      ct = 0;
      for ( int i = 0; i < dim; ++i )
      {
         for ( int j = 0; j < dim; ++j )
         {
            subBB[ct] = Q[i*m + j];
            assert( subBB[ct] == Com_dot( &B_[i*m], &B_[j*m], m ) );
            ++ct;
         }
      }
      assert( ct == dim * dim );

      for ( auto i = 0; i < dim; ++i )
      {
         assert( e[i] == 0.0 );
         e[i] = 1.0;

         int r = Com_LS_dposv( &subBB[0], &e[0], dim, &buf_a[0] );
         if ( r != 0 )
         {
            cout << "error(vscb): r = " << r << endl;
            exit(-1);
         }

         Com_mat_Ax( B_, m, dim, &buf_a[0], &buf_b[0] );

         assert( dim - 1 < m && dim - 1 >= 0 );
         Copy_vec( &buf_b[0], &vscb[dim-1][i][0], m );
         lvscb[dim-1][i] = 1.0 / Com_dot( &buf_b[0], &buf_b[0], m );
         //printf("%.8f ", lvscb[dim-1][i] );

         e[i] = 0.0;

         assert( Equal( lvscb[dim-1][i], 1.0/ Com_dot( &vscb[dim-1][i][0], &vscb[dim-1][i][0], m ), epsilon ) );
      }
      //cout << endl;
      ++dim;
   }

   //for ( auto i : vscb )
   //{
   //   for ( auto j : i )
   //   {
   //      for ( auto k : j )
   //         printf("%3.10f ", k);
   //      cout << endl;
   //   }
   //   cout << endl;
   //}
   //cout << endl;
}

void  SVPsolver::SVPScomputeVsCO()
{
   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto B_ = pd.get_B_();
   const auto& vscb = VsCB;
   auto& vsco = VsCO;

   assert( m > 0 );
   assert( B_ != nullptr );
   assert( !vscb.empty() );
   assert( m == (int) vscb.size() );
   assert( vsco.empty() );

   vsco.resize( m - 1 );

   for ( auto i = 0; i < m - 1; ++i )
   {
      vsco[i].resize( i + 1 );
      const auto& vscb_i = vscb[i];
      const auto b_i1 = &B_[m*(i+1)];
      assert( !vscb_i.empty() );
      assert( i + 1 == (int) vscb_i.size() );

      for ( auto j = 0; j < i + 1; ++j )
      {
         const auto& vscb_i_j = vscb_i[j];
         auto& vsco_i_j = vsco[i][j];
         assert( !vscb_i_j.empty() );
         assert( m == (int) vscb_i_j.size() );

         vsco_i_j = 0;
         for ( auto k = 0; k < m; ++k )
            vsco_i_j += vscb_i_j[k] * b_i1[k];
      }
   }

   //for ( auto v : vsco )
   //{
   //   for ( auto val : v )
   //      cout << val << " ";
   //   cout << endl;
   //}
}

void  SVPsolver::SVPScomputeVCOV()
{
   assert( !ENUM );

   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto Q = pd.get_Q();
   const auto& vsco = VsCO;
   auto& vcov = VCOV;

   assert( m > 0 );
   assert( Q != nullptr );
   assert( !vsco.empty() );
   assert( m  - 1 == (int) vsco.size() );
   assert( vcov.empty() );

   vcov.resize( m - 1 );

   for ( auto i = 0; i < m - 1; ++i )
   {
      const auto& vsco_i = vsco[i];
      const auto q = &Q[m*(i+1)];
      const auto bb = &Q[m*(i+1)+i+1];

      assert( !vsco_i.empty() );
      assert( (int) vsco_i.size() == i + 1 );
      assert( q != nullptr );
      assert( bb != nullptr );

      vcov[i] = *bb;
      vcov[i] -= inner_product( vsco_i.begin(), vsco_i.end(), q, 0.0 );
   }
}

void  SVPsolver::SVPScomputeDBM()
{
   constexpr int maxsize_dim = 30;
   constexpr int maxsize_numsol = 10;
   auto& dbms = DBMs;
   auto& rsndbm = RSNDBM;

   assert( probdata.get_m() - 1 >= maxsize_dim );

   // alloc
   dbms.resize( maxsize_dim );
   rsndbm.resize( maxsize_dim );

   // set
   for ( auto i = 0; i < maxsize_dim; ++i )
   {

      SOLUTION_POOL solpool;
      FindNShorterVec( maxsize_numsol, 10, VsCB[i], &solpool );
      const auto sols = solpool.cget_sols();
      const auto num = solpool.get_num();

      dbms[i].resize( num );
      rsndbm[i].resize( num );

      for ( auto j = 0; j < num; ++j )
      {
         assert( i + 1 == sols[j].get_dim() );

         const auto solvals = sols[j].get_solval();
         assert( solvals != nullptr );

         dbms[i][j].resize( i+1 );
         for ( auto k = 0; k <= i; ++k )
            dbms[i][j][k] = solvals[k];
         rsndbm[i][j] = 1.0 / sols[j].get_objval();
      }
   }

   //int ct1 = 0;
   //for( auto matrix : dbms )
   //{
   //   int ct2 = 0;
   //   for ( auto vec : matrix )
   //   {
   //      printf("%f ", rsndbm[ct1][ct2++] );
   //      for ( auto val : vec )
   //         cout << val << " ";
   //      cout << endl;
   //   }
   //   cout << endl;
   //   ct1++;
   //}
}

void  SVPsolver::disp_log(
   const NODE&       node,
   const RelaxResult r,
   const int         index,
   const int         cutoff
   )
{
   auto k = node.get_index();
   assert( k >= 0 );

   cout << stopwatch.get_time() << "s";
   cout << " [" << k << "/";
   cout << nodelist.getListsize() << "(" << index-1 << ")] ";
   cout << "dpt:" << node.get_dpt() << " ";
   cout << "[ " << node.get_lowerbound() << ", ";
   cout << GLB << ", ";
   cout << bestval << ", ";
   cout << 100*(bestval - GLB)/bestval << "%]";

   cout << " cutoff:" << cutoff;
   cout << " AF:";
   cout << fixed << setprecision(3) << Appfac;

   cout << endl;
}

double SVPsolver::compute_objval(
   double   *x
   )
{
   int      m = probdata.get_m();
   double   *Q = probdata.get_Q();

   assert( m > 0 );
   assert( Q != nullptr );
   assert( x != nullptr );

   double   obj = 0.0;

   double   *Qx; //[m]
   Qx = new double[m];

   Gen_ZeroVec( m, Qx);
   Com_mat_Ax( Q, m, m, x, Qx);

   obj += Com_dot( x, Qx, m);

   delete[] Qx;

   return obj;
}

void SVPsolver::SVPStrySol(
      SOLUTION&   sol,
      const bool  check_lb,
      const bool  check_vbs,
      bool*       result
      )
{
   assert( sol.get_solval() != nullptr );

   // copy
   auto& solpool = pool;
   auto* currentbestval = &bestval;
   auto* current_sq_bestval = &sq_bestval;

   auto objval = sol.get_objval();
   auto update = false;

   assert( objval > 0.0 );

   solpool.add_solution( sol );

   if ( *currentbestval > objval )
   {
      *currentbestval = objval;
      *current_sq_bestval = sqrt( objval );
      bestsol = sol;
      Appfac = *current_sq_bestval / _Appfac;
      update = true;
   }

   if ( result != nullptr )
      *result = update;

   if ( update )
   {
      // update global bounds on variables
      if ( check_vbs )
      {
         //
      }

      if ( check_lb || check_vbs )
      {
         //
      }
   }
}

void SVPsolver::SVPSsetGS()
{
   const auto& pd = probdata;
   const auto m = pd.get_m();
   const auto m_1 = m - 1;
   const auto B_ = pd.get_B_();
   auto& snovs = SNOVs;
   auto& cmgso = CMGSO;

   assert( m > 0 );

   double* OM = new double[m*m];   // orthogonal matrix
   vector<vector<double>> CF( m );          // coefficient matrix
   double coef;

   snovs.resize( m );
   cmgso.resize( m );


   // set SNOVs
   for( int i = 0, im = 0; i < m; ++i, im += m )
   {
      assert( i * m == im );

      Copy_vec( &B_[im], &OM[im], m );

      CF[i].resize( i + 1 );
      for ( int j = 0, jm = 0; j < i; ++j, jm += m )
      {
         assert( j * m == jm );

         coef = Com_dot( &B_[im], &OM[jm], m ) / snovs[j];
         Com_linecomb( &OM[im], &OM[jm], m, 1.0, - coef, &OM[im] );
         CF[i][j] = coef;
      }
      CF[i][i] = 1;

      snovs[i] = Com_dot( &OM[im], &OM[im], m );
   }

   delete[] OM;
   OM = nullptr;

   //for ( auto cfs : CF )
   //{
   //   for ( auto cf : cfs )
   //      cout << cf << " ";
   //   cout << endl;
   //}

   // set CMGSO
   for ( auto i = 0; i < m; ++i )
   {
      cmgso[i].resize( i + 1 );
      for ( auto j = 0; j <= i; ++j )
      {
         assert( !CF[m_1 - j].empty() );
         assert( m_1 - i < (int) CF[m_1 - j].size() );

         cmgso[i][j] = CF[m_1 - j][m_1 - i]; // m_1 = m - 1
      }
   }

   //for ( auto cfs : CMGSO )
   //{
   //   for ( auto cf : cfs )
   //      cout << cf << " ";
   //   cout << endl;
   //}

}
