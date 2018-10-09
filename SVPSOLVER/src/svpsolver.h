#ifndef SVPSOLVER_H__
#define SVPSOLVER_H__

//#include <vector>
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

   int*  ub;
   int*  lb;

   double         GLB;
   double         bestval;
   SOLUTION       bestsol;
   SOLUTION_POOL  pool;
   double         _Appfac;
   double         Appfac;

   TYPE_NODELIST  type;
   NODELIST       nodelist;

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

   int   index;

   STOPWATCH      stopwatch;
   TESTWATCH      testwatch;

   unsigned long int    nnode;

   SCHMIDT_M      sch;

   bool quiet;

   double   *norm;

   CUT_POOL    oa_cpool;   //a pool of cutting planes by using
                           //outer approximation

   Status status;

   int TIMELIMIT;
   int LEFTNODELIMIT;
   int NODELIMIT;

   // for relaxation
   QP_DATA  qpdata;

   // for parallel
   int   nthreads;
   bool  subsolver;
   mutex mtx;
   condition_variable cv;

   double epsilon;

   public:

      SVPsolver();                                 // default constructor
      SVPsolver( const SVPsolver &source );        // copy constructor
      SVPsolver& operator=( const SVPsolver& );    // assignment operator
      ~SVPsolver();                                // destructor

      // setup {
      void        SVPSsetup(  const int s_m, const double* s_B_, const int s_nthreads,
                              const int s_timelimit, const bool s_quiet, const bool w_subsolver,
                              const bool w_sch, const bool w_bounds, const bool w_heur,
                              const bool w_app, const bool w_grn, const bool w_nl );
      void        SVPScreateProbdata( const int m, const double *B_);

      void        SVPSheurFindMinColumn();
      void        SVPScomputeBounds();
      void        SVPSgenerateRootNode( bool w_sch );
      void        SVPSsetupNodelist() {
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
      bool        SVPSparasolve();
      void        SVPSsolveSubprob( int& n_running_threads, double& sublb_i,
                                    const int thread_id );
      void        SVPSresetIndex();
      // } solve via parallel

      // branch-and-bound {
      bool        SVPScheckLimit();
      // } branch-and-bound

      // relaxation {
      RelaxResult SVPSsolveRelaxation( NODE& node );
      RelaxResult SVPSsolveRelaxationBIN( NODE& node );
      RelaxResult SVPSsolveRelaxationINT( NODE& node );
      // } relaxation

      // node selection {
      int   SVPSselectNode( int index, int disp );
      // }  node selection

      // branch {
      void        SVPSbranch( NODE& node, int index );
      void        SVPSbranch_BIN( NODE& node, int index );
      void        SVPSbranch_INT( NODE& node, int index );
      // } branch

      // heuristics {
      void        SVPSheur( const NODE& node );
      void        SVPSheurUnitsphere( const NODE& node );
      void        SVPSheurQuadratic( const NODE& node );
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
      void        tighten_bounds( int memo, int* set_lb, int* set_ub );
      void        disp_log( const NODE& node, const RelaxResult r, const int index, const int cutoff);
      void        disp_bestsol();
      double      compute_objval( double *x );
      auto        SVPSgetBestval(){ return bestval; }
      auto        SVPSgetBestsol(){ return bestsol; }
      PROB_DATA   get_probdata(){ return probdata; }
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
      void        SVPSsetBounds( const int* u, const int* l )
      {
         assert( u != nullptr );
         assert( l != nullptr );
         assert( ub != nullptr );
         assert( lb != nullptr );

         int m = probdata.get_m();
         for( int i = 0; i < m; i++ )
         {
            ub[i] = u[i];
            lb[i] = l[i];
         }
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
