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

#define  BRANCHINGRULE_INT 3
      // 0:
      // 1:
      // 2:
      // 3: ver 2.0
      // 4: ver 2.1
      // 5: ver 2.1
      // 6:

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

   double         *ub;
   double         *lb;

   double         GLB;
   double         bestval;
   SOLUTION       bestsol;
   SOLUTION_POOL  pool;
   double         _Appfac;
   double         Appfac;

   list<NODE>     NodeList;
   int   listsize;
   int   index;

   STOPWATCH      stopwatch;
   TESTWATCH      testwatch;

   unsigned long int    nnode;

   int            *order;     // order of norm or bounds

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
                              const bool w_app, const bool w_grn );
      void        SVPScreateProbdata( const int m, const double *B_);

      void        SVPSheurFindMinColumn();
      void        SVPScomputeBounds();
      void        SVPSgenerateRootNode( bool w_sch );
      // } setup

      // setup for subsolver {
      void        SVPSsetGlobalLowerBound( const double s_GLB ){ GLB = s_GLB; }
      void        SVPSsetNode( const NODE &setnode );
      NODE&       SVPSgetNode();
      void        SVPSpopNode();
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
      void        SVPSsolveSubprob( unsigned short int& n_running_threads, double& sublb_i,
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

      // get
      double      SVPSgetGlobalLowerBound() const { return GLB; }
      string      SVPSgetStringStatus() const;

      void        SVPScreateSch( int m, double* B_ );
      auto        SVPSgetStatus(){ return status; }
      //bool        p_solve();
      void        tighten_bounds( int memo, double* set_lb, double* set_ub );
      void        SVPSheur( int sel );
      void        heur_unitsphere(int sel);
      void        heur_quadratic(int sel);
      void        disp_log(int   sel, RelaxResult r, int index, int cutoff);
      void        disp_bestsol();
      double      compute_objval( double *x );
      void        branch(int  sel, int index);
      void        branch_BIN(int sel, int index);
      void        branch_INT(int sel, int index);
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
      auto        SVPSgetListsize(){ return listsize; }
      double      get_gap(){ return 100*(bestval - GLB)/bestval; }

      // for parallel mode
      void        set_num_thread(int n){ nthreads = n; }
      void        SVPSsetBestval( const double val ){ bestval = val; }
      void        SVPSsetBounds(const double *u, const double *l)
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
      void        SVPSsetOrder( const int *s_order)
      {
         assert( s_order != nullptr );
         int m = probdata.get_m();
         order = new int[m];
         for (int i = 0; i < m; i++)
            order[i] = s_order[i];
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
