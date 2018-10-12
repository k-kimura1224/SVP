#ifndef SVPSOLVER_H__
#define SVPSOLVER_H__

#include <vector>
#include <list>
#include <mutex>
#include <condition_variable>

#include "probdata.h"
#include "qpdata.h"
#include "solution.h"
#include "solution_pool.h"
#include "node.h"
#include "stopwatch.h"
#include "testwatch.h"
#include "Schmidt_manager.h"
#include "cut_pool.h"
#include "nodelist.h"

#define  CUT_OA   false

#define  HEUR_APP    0.95

// for parallel
#define  NUM_INITNODES     30000
#define  PARASIZE          100

#define  PARA_LOG          true


using namespace std;

enum RelaxResult {
   UPDATE,
   FEASIBLE,
   INFEASIBLE,
   GETINTEGER
};

enum Status {
   SETUP,
   SOLVING,
   SOLVED,
   TIMEOVER,
   FULL_OF_LEFTNODES,
   FULL_OF_NODES
};


class SVPsolver{
   PROB_DATA      probdata;

   vector<vector<int>>  bounds;

   double         GLB;
   double         bestval;
   SOLUTION       bestsol;
   SOLUTION_POOL  pool;
   double         _Appfac;
   double         Appfac;

   int                  index;
   TYPE_NODELIST        type;
   NODELIST             nodelist;
   unsigned long int    nnode;

   void     (NODELIST::*push_back)( const NODE& );
   void     (NODELIST::*move_back)( NODE& );
   NODE&    (NODELIST::*nodeselection)( double* globallowerbound, const double bestval, const int index, const int disp );
   void     (NODELIST::*cut_off)();
   double   (NODELIST::*get_GLB)() const;
   bool     (NODELIST::*check_size)() const;
   int      (NODELIST::*setup_para_selection)() const;
   NODE&    (NODELIST::*para_selection)( const int );
   void     (NODELIST::*pop_front)( const int );
   int      (NODELIST::*getSubsize)( const int ) const;

   // for heuristics
   double   *norm;

   // for relaxation
   QP_DATA  qpdata;

   // for parallel
   int   nthreads;
   bool  subsolver;
   mutex mtx;
   condition_variable cv;

   // parameter
   int   TIMELIMIT;
   int   LEFTNODELIMIT;
   int   NODELIMIT;
   bool  quiet;

   Status status;

   CUT_POOL    oa_cpool;   //a pool of cutting planes by using
                           //outer approximation

   STOPWATCH      stopwatch;
   TESTWATCH      testwatch;

   double   epsilon;

   public:

      SVPsolver();                                 // default constructor
      SVPsolver( const SVPsolver &source );        // copy constructor
      SVPsolver& operator=( const SVPsolver& );    // assignment operator
      ~SVPsolver();                                // destructor

      // setup {
      void  SVPSsetup(  const int s_m, const double* s_B_, const int s_nthreads,
                        const int s_timelimit, const bool s_quiet,
                        const bool w_subsolver, const bool w_bounds,
                        const bool w_heur, const bool w_app, const bool w_gn,
                        const bool w_nl );
      void  SVPScreateProbdata( const int m, const double *B_);

      void  SVPSheurFindMinColumn();
      void  SVPScomputeBounds();
      void  SVPSgenerateNodes();
      void  SVPStightenBounds( const int memo );
      void  SVPSsetupNodelist() {
         nodelist.setup( type, bestval );
         switch ( type )
         {
            case TWO_DEQUE:
            {
               push_back = &NODELIST::push_back_TDEQUE;
               move_back = &NODELIST::move_back_TDEQUE;
               nodeselection = &NODELIST::nodeselection_TDEQUE;
               cut_off = &NODELIST::cutoff_TDEQUE;
               get_GLB = &NODELIST::get_GLB_TDEQUE;
               check_size = &NODELIST::check_size_TDEQUE;
               setup_para_selection = &NODELIST::setup_para_selection_TDEQUE;
               para_selection = &NODELIST::para_selection_TDEQUE;
               pop_front = &NODELIST::pop_front_TDEQUE;
               getSubsize = &NODELIST::getSubsize_TDEQUE;
               break;
            }
            default:
            {
               push_back = &NODELIST::push_back_LIST;
               move_back = &NODELIST::move_back_LIST;
               nodeselection = &NODELIST::nodeselection_LIST;
               cut_off = &NODELIST::cutoff_LIST;
               get_GLB = &NODELIST::get_GLB_LIST;
               check_size = &NODELIST::check_size_LIST;
               setup_para_selection = &NODELIST::setup_para_selection_LIST;
               para_selection = &NODELIST::para_selection_LIST;
               pop_front = &NODELIST::pop_front_LIST;
               getSubsize = &NODELIST::getSubsize_LIST;
               break;
            }
         }
      }
      // } setup

      // setup for subsolver {
      void        SVPSsetGlobalLowerBound( const double s_GLB ){ GLB = s_GLB; }
      void        SVPSmoveNode( NODE& setnode );
      NODE&       SVPSgetNode_para_selection( const int setup );
      void        SVPSpopNode( const int setup );
      // } setup for subsolver

      // solve {
      bool        SVPSsolve();
      void        SVPSstartTime();
      void        SVPSoutputBounds();
      bool        SVPSrunBranchandBound();
      bool        SVPSresolve();
      // } solve

      // solve via parallel {
      bool  SVPSparasolve();
      void  SVPSsolveSubprob( int& n_running_threads, double& sublb_i,
                              const int thread_id );
      void  SVPSresetIndex();
      // } solve via parallel

      // branch-and-bound {
      bool  SVPScheckLimit();
      void  SVPSgetVarsLocalBound( NODE& node, double* vars_localub, double* vars_locallb );
      void  SVPStrySol( SOLUTION& sol, const bool check_lb,
                        const bool check_vbs, bool* result );
      // } branch-and-bound

      // relaxation {
      RelaxResult SVPSsolveRelaxation( NODE& node, const double* vars_localub, const double* vars_locallb );
      RelaxResult SVPSsolveRelaxationINT( NODE& node, const double* vars_localub, const double* vars_locallb );
      // } relaxation

      // node selection {
      int   SVPSselectNode( int index, int disp );
      // }  node selection

      // branch {
      void  SVPSbranch( NODE& node, const int index,
                        double* vars_localub, double* vars_locallb );
      void  SVPSbranch_BIN( NODE& node, const int index,
                        double* vars_localub, double* vars_locallb );
      void  SVPSbranch_INT( NODE& node, const int index,
                        double* vars_localub, double* vars_locallb );
      // } branch

      // heuristics {
      void        SVPSheur( const NODE& node, const double* vars_localub, const double* vars_locallb );
      void        SVPSheurUnitsphere( const NODE& node, const double* vars_localub, const double* vars_locallb );
      void        SVPSheurQuadratic( const NODE& node, const double* vars_localub, const double* vars_locallb );
      // } heuristics

      // nodelist {
      auto        SVPSgetListsize() const { return nodelist.getListsize(); }
      auto        SVPSgetSubsize( const int sub ) const { return (nodelist.*getSubsize)(sub); }
      auto        SVPSgetSetupParaSelection() const { return (nodelist.*setup_para_selection)(); }
      // } nodelist
      //
      // get
      double      SVPSgetGlobalLowerBound() const { return GLB; }
      string      SVPSgetStringStatus() const;

      void        SVPScreateSch( int m, double* B_ );
      auto        SVPSgetStatus(){ return status; }
      //bool        p_solve();
      void        disp_log( const NODE& node, const RelaxResult r, const int index, const int cutoff);
      void        disp_bestsol();
      double      compute_objval( double *x );
      auto        SVPSgetBestval(){ return bestval; }
      auto        SVPSgetBestsol(){ return bestsol; }
      PROB_DATA&  SVPSgetProbdata(){ return probdata; }
      unsigned long int SVPSgetNnode(){ return nnode; }
      void        gene_OAcuts( double *u, double *l, double *Q, double M);

      void        SVPSsetTimelimit( int tlimit ){ TIMELIMIT = tlimit; }
      int         SVPSgetTimelimit(){ return TIMELIMIT; }
      void        SVPSsetLeftnodelimit( const int nlimit ){
                     assert( nlimit > 0 && nlimit <= 2000000000 );
                     LEFTNODELIMIT = nlimit;
                  }
      void        SVPSsetNodelimit( const int nlimit ){
                     assert( nlimit > 0 && nlimit <= 2000000000 );
                     NODELIMIT = nlimit;
                  }

      int         get_runtime(){ return stopwatch.get_result(); }
      double      get_gap(){ return 100*(bestval - GLB)/bestval; }

      // for parallel mode
      void        set_num_thread(int n){ nthreads = n; }
      void        SVPSsetBestval( const double val ){ bestval = val; }
      void        SVPSsetBounds( const vector<vector<int>>& s_bounds )
      {
         assert( !s_bounds.empty() );
         assert( bounds.empty() );
         bounds = s_bounds;
         bounds.shrink_to_fit();
      }
      void        SVPSsetNorm( const double *s_norm)
      {
         assert( probdata.get_m() > 0 );
         assert( s_norm != nullptr );
         assert( norm == nullptr );
         int m = probdata.get_m();
         norm = new double[m];
         Copy_vec( s_norm, norm, m );
      }
      void        SVPSsetAppfac( const double a, const double _a )
      {
         Appfac = a;
         _Appfac = _a;
      }
      void        set_quiet( bool q )
      {
         quiet = q;
      }
};

#endif
