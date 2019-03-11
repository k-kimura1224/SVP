#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <assert.h>

#include "svpsolver.h"
#include "probdata.h"
#include "vector.h"
#include "solution.h"
#include "solution_pool.h"

using namespace std;

// heur_unitsphere.cpp
#define  HEUR_UNITSPHERE            false
#define  HEUR_UNITSPHERE_FREQ       1
#define  HEUR_UNITSPHERE_FREQOFS    0
#define  HEUR_UNITSPHERE_DPT        40

#define  PARA_HEUR_UNITSPHERE             false
#define  PARA_HEUR_UNITSPHERE_FREQ        5
#define  PARA_HEUR_UNITSPHERE_FREQOFS     20
#define  PARA_HEUR_UNITSPHERE_DPT         20

// heur_quadratic.cpp
#define  HEUR_QUADRATIC          true
#define  HEUR_QUADRATIC_FREQ        3
#define  HEUR_QUADRATIC_FREQOFS     0
#define  HEUR_QUADRATIC_DPT         30

#define  PARA_HEUR_QUADRATIC              true
#define  PARA_HEUR_QUADRATIC_FREQ         5
#define  PARA_HEUR_QUADRATIC_FREQOFS      0
#define  PARA_HEUR_QUADRATIC_DPT       20

#define debug  0

void  SVPsolver::SVPSheur(
      const NODE&    node,
      const double*  vars_localub,
      const double*  vars_locallb
      )
{
   const auto dpt = node.get_dpt();
   int dpt2;

   assert( dpt >= 0 );
   assert( HEUR == true );

   if( subsolver == false )
   {
      // single
      // unit sphere algorithm
      dpt2 = dpt - (int)HEUR_UNITSPHERE_FREQOFS;
      if( dpt2 >= 0
         && HEUR_UNITSPHERE
         && !(dpt2 % (int)HEUR_UNITSPHERE_FREQ)
         && dpt2 <= (int)HEUR_UNITSPHERE_DPT )
      {
         SVPSheurUnitsphere( node, vars_localub, vars_locallb );
      }
      // solve quadratic optimization with one variable
      dpt2 = dpt - (int)HEUR_QUADRATIC_FREQOFS;
      if( dpt2 >= 0
         && HEUR_QUADRATIC
         && !(dpt2 % (int)HEUR_QUADRATIC_FREQ)
         && dpt2 <= (int)HEUR_QUADRATIC_DPT )
      {
         SVPSheurQuadratic( node, vars_localub, vars_locallb );
      }
   }
   else
   {
      // parallel
      // unit sphere algorithm
      dpt2 = dpt - (int)PARA_HEUR_UNITSPHERE_FREQOFS;
      if( dpt2 >= 0
         && PARA_HEUR_UNITSPHERE
         && !(dpt2 % (int)PARA_HEUR_UNITSPHERE_FREQ)
         && dpt2 <= (int)PARA_HEUR_UNITSPHERE_DPT )
      {
         SVPSheurUnitsphere( node, vars_localub, vars_locallb );
      }
      // solve quadratic optimization with one variable
      dpt2 = dpt - (int)PARA_HEUR_QUADRATIC_FREQOFS;
      if( dpt2 >= 0
         && PARA_HEUR_QUADRATIC
         && !(dpt2 % (int)PARA_HEUR_QUADRATIC_FREQ)
         && dpt2 <= (int)PARA_HEUR_QUADRATIC_DPT )
      {
         SVPSheurQuadratic( node, vars_localub, vars_locallb );
      }
   }
}

