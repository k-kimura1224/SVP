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
#include "node.h"

using namespace std;

#define debug  0


int   SVPsolver::SVPSselectNode(
   int index,
   int disp
   )
{
   int   rslt = -1;
   auto  buf = (index-1)%100000;

   if( ( buf == 0 || buf == 1 ) && disp < index )
   {
      if ( quiet )
      {
         NodeList.sort();
      }
      else
      {
         clock_t start = clock();
         NodeList.sort();
         clock_t end = clock();
         cout << "Sorting Time: ";
         cout << (double)(end-start)/CLOCKS_PER_SEC;
         cout << "s" << endl;
      }
      GLB = NodeList.begin()->get_lowerbound();
   }

   rslt = 0;

   assert( rslt >= 0 );
   assert( rslt < listsize );

   return rslt;
}
