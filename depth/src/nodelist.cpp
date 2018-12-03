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
#include <deque>

#include "nodelist.h"
#include "node.h"
#include "vector.h"

using namespace std;

#define debug_class  0
#define debug  0
#define SIZE_RESERVE 100000

NODELIST::NODELIST(){   // default constructor
#if debug_class
   cout << "NODELIST: default constructor" << endl;
#endif
   type = LIST;
   listsize = 0;
   standardvalue = -1.0;
   standardgap = -1.0;
   division = -1.0;
   maxnode = -1;
   currentdeq = -1;
   focusdeq = -1;
}

NODELIST::NODELIST( const NODELIST &source )
{  // copy constructor
#if debug_class
   cout << "NODELIST: copy constructor" << endl;
#endif

   type = source.type;
   standardvalue = source.standardvalue;
   standardgap = source.standardgap;
   listsize = source.listsize;
   division = source.division;
   maxnode = source.maxnode;
   currentdeq = source.currentdeq;
   focusdeq = source.focusdeq;

   copy( source.node_list.begin(), source.node_list.end(), back_inserter(node_list) );
   copy( source.node_deque_1.begin(), source.node_deque_1.end(), back_inserter(node_deque_1) );
   copy( source.node_deque_2.begin(), source.node_deque_2.end(), back_inserter(node_deque_2) );
   assert(0);
   copy( source.node_deque.begin(), source.node_deque.end(), back_inserter(node_deque) );
}

// assignment operator
NODELIST& NODELIST::operator=( const NODELIST& source )
{
#if debug_class
   cout << "NODELIST: assignment operator" << endl;
#endif

   if( this != &source )
   {
      type = source.type;
      standardvalue = source.standardvalue;
      standardgap = source.standardgap;
      listsize = source.listsize;
      division = source.division;
      maxnode = source.maxnode;
      currentdeq = source.currentdeq;
      focusdeq = source.focusdeq;

      copy( source.node_list.begin(), source.node_list.end(), back_inserter(node_list) );
      copy( source.node_deque_1.begin(), source.node_deque_1.end(), back_inserter(node_deque_1) );
      copy( source.node_deque_2.begin(), source.node_deque_2.end(), back_inserter(node_deque_2) );
      assert(0);
      copy( source.node_deque.begin(), source.node_deque.end(), back_inserter(node_deque) );
   }

   return *this;
}

// destructor
NODELIST::~NODELIST()
{
#if debug_class
   cout << "NODELIST: destructor" << endl;
#endif
   list<NODE>().swap( node_list );

   node_deque_1.clear();
   node_deque_1.shrink_to_fit();
   node_deque_2.clear();
   node_deque_2.shrink_to_fit();

   for ( auto& deq : node_deque )
   {
      deq.clear();
      deq.shrink_to_fit();
   }
   node_deque.clear();
   node_deque.shrink_to_fit();
}

// setup
void NODELIST::setup(
      const TYPE_NODELIST  s_type,
      const double         bestval,
      const int            s_memory,
      const int            m,
      const double         gap
      )
{
   type = s_type;

   switch ( type )
   {
      case TWO_DEQUE:
      {
         standardgap = gap;
         standardvalue = ( 1.0 - standardgap ) * bestval;

         break;
      }
      case TEN_DEQUE:
      {
         assert( bestval > 0.0 );
         assert( s_memory > 0 );
         assert( m > 0 );

         if ( node_deque.empty() )
         {
            node_deque.resize( 10 );
         }
         else
         {
            assert( (int) node_deque.size() == 10 );
            for ( auto& deq : node_deque )
            {
               deq.clear();
               deq.shrink_to_fit();
            }
         }

         division = bestval * 0.1;

         double buf = s_memory;
         buf /= 2;
         //buf /= 4;
         buf /= 8 * m + 18;
         buf *= 1024;
         buf *= 1024;
         buf *= 1024;
         maxnode = buf;

         currentdeq = 0;
         break;
      }
      default:
      {
         break;
      }
   }
}

int   NODELIST::getSubsize_TDEQUE( const int sub ) const {
   assert( sub == 1 || sub == 2 );
   if ( sub == 1 )
   {
      if ( node_deque_1.empty() )
         return 0;
      else
         return (int) node_deque_1.size();
   }
   else
   {
      if ( node_deque_2.empty() )
         return 0;
      else
         return (int) node_deque_2.size();
   }
}

// push
void NODELIST::push_back_LIST(
      const NODE&    node
      )
{
   assert( type == LIST );
   node_list.push_back( node );
   ++listsize;
}

void NODELIST::push_back_TDEQUE(
      const NODE&    node
      )
{
   assert( type == TWO_DEQUE );
   assert( standardvalue > 0.0 );

  if ( node.get_lowerbound() > standardvalue )
      node_deque_2.push_back( node );
   else if ( Equal( standardgap, 0.1, 1.0e-12 ) )
      node_deque_2.push_back( node );
   else
      node_deque_1.push_back( node );

   ++listsize;
}

// move
void NODELIST::move_back_LIST(
      NODE&    node
      )
{
   assert( type == LIST );
   node_list.push_back( move( node ) );
   ++listsize;
}

void NODELIST::move_back_TDEQUE(
      NODE&    node
      )
{
   assert( type == TWO_DEQUE );
   assert( standardvalue > 0.0 );

   if ( node.get_lowerbound() > standardvalue )
      node_deque_2.push_back( move( node ) );
   else if ( Equal( standardgap, 0.1, 1.0e-12 ) )
      node_deque_2.push_back( move( node ) );
   else
      node_deque_1.push_back( move( node ) );

   ++listsize;
}

// nodeselection
NODE* NODELIST::nodeselection_LIST(
      double*        globallowerbound,
      const double   bestval,
      const int      index,
      const int      disp
      )
{
   assert(0);
   //assert( (int)node_list.size() == listsize );
   //assert( listsize > 0 );

   //auto  buf = (index-1)%100000;

   //if( ( buf == 0 || buf == 1 ) && disp < index )
   //{
   //   node_list.sort();
   //   *globallowerbound = node_list.begin()->get_lowerbound();
   //}

   //return *(node_list.begin());
   return nullptr;
}

NODE* NODELIST::nodeselection_TDEQUE(
      double*        globallowerbound,
      const double   bestval,
      const int      index,
      const int      disp
      )
{
   assert(0);
   //assert( listsize > 0 );
   //assert( !node_deque_1.empty() || !node_deque_2.empty() );

   ////if ( listsize > 500000 )
   ////{
   ////   if ( !node_deque_2.empty() )
   ////      return node_deque_2.front();
   ////   else
   ////      return node_deque_1.front();
   ////}

   //if ( !node_deque_1.empty() )
   //   return node_deque_1.front();

   //assert( standardvalue > 0.0 );
   //assert( standardgap > 0.0 && standardgap < 1.0 );

   //if ( Equal( standardgap, 0.1, 1.0e-12 ) )
   //   return node_deque_2.front();

   //*globallowerbound = standardvalue;
   //standardgap -= 0.1;
   //standardvalue = ( 1.0 - standardgap ) * bestval;

   //if ( (int) node_deque_2.size() == 1 )
   //   return node_deque_2.front();

   ////cout << "node_selection - start" << endl;
   //do {
   //   assert( !node_deque_2.empty() );

   //   if ( Equal( standardgap, 0.1, 1.0e-12 ) )
   //      return node_deque_2.front();

   //   int size = (int) node_deque_2.size();
   //   for ( int i = 0; i < size; )
   //   {
   //      //cout << "i:" << i << ", size2:" << node_deque_2.size();
   //      //cout << ", size1:" << node_deque_1.size() << endl;
   //      if ( node_deque_2[i].get_lowerbound() < standardvalue )
   //      {
   //         node_deque_1.push_back( move( node_deque_2[i] ) );

   //         if ( (int) node_deque_2.size() > 1 )
   //            node_deque_2[i] = move( node_deque_2.back() );

   //         node_deque_2.pop_back();
   //         --size;
   //         assert( size == 0 || size == (int) node_deque_2.size() );
   //         continue;
   //      }
   //      ++i;
   //   }
   //   //cout << "for-end" << endl;
   //   standardgap -= 0.1;
   //   standardvalue = ( 1.0 - standardgap ) * bestval;
   //} while ( node_deque_1.empty() );
   ////cout << "node_selection - end" << endl;

   //assert( !node_deque_1.empty() );

   //return node_deque_1.front();
   return nullptr;
}

//NODE* NODELIST::nodeselection_TENDEQUE(
//      double*        globallowerbound,
//      const double   bestval,
//      const int      index,
//      const int      disp
//      )
//{
//   assert( listsize > 0 );
//   assert( type == TEN_DEQUE );
//   assert( division > 0.0 );
//   assert( maxnode > 0 );
//   assert( !node_deque.empty() );
//   assert( (int)node_deque.size() == 10 );
//
//   if ( maxnode > listsize )
//   {
//      for ( auto i = currentdeq; i < 10; ++i )
//      {
//         if ( !node_deque[i].empty() )
//         {
//            return &node_deque[i].front();
//         }
//         else
//         {
//            currentdeq = i + 1;
//            *globallowerbound =  division * ( 1 + i );
//            node_deque[i].shrink_to_fit();
//         }
//      }
//   }
//   else
//   {
//      for ( auto i = 9; i >= currentdeq; --i )
//      {
//         if ( !node_deque[i].empty() )
//            return &node_deque[i].front();
//      }
//   }
//
//   printf("error!\n");
//   assert(0);
//   exit(-1);
//
//   return &node_deque[0].front();
//}

double NODELIST::get_GLB_LIST() const
{
   auto min_lb = node_list.begin()->get_lowerbound();
   double lb_i;

   for ( auto node = node_list.begin(); node != node_list.end(); ++node )
   {
      lb_i = node->get_lowerbound();
      if( min_lb > lb_i )
         min_lb = lb_i;
   }

   return min_lb;
}

double NODELIST::get_GLB_TDEQUE() const
{
   assert( !node_deque_1.empty() || !node_deque_2.empty() );

   double min_lb;
   double lb_i;

   if ( !node_deque_1.empty() )
   {
      min_lb = node_deque_1.begin()->get_lowerbound();
      for ( auto node = node_deque_1.begin(); node != node_deque_1.end(); ++node )
      {
         lb_i = node->get_lowerbound();
         if( min_lb > lb_i )
            min_lb = lb_i;
      }

      return min_lb;
   }

   if ( !node_deque_2.empty() )
   {
      min_lb = node_deque_2.begin()->get_lowerbound();
      for ( auto node = node_deque_2.begin(); node != node_deque_2.end(); ++node )
      {
         lb_i = node->get_lowerbound();
         if( min_lb > lb_i )
            min_lb = lb_i;
      }

      return min_lb;
   }

   return -1;
}

double NODELIST::get_GLB_TENDEQUE() const
{
   assert( listsize > 0 );
   assert( type == TEN_DEQUE );
   assert( division > 0.0 );
   assert( maxnode > 0 );
   assert( !node_deque.empty() );
   assert( (int)node_deque.size() == 10 );
   assert( currentdeq >= 0 && currentdeq < 10 );

   auto min_lb = division * 100.0;
   double lb_i;

   for ( auto& deq : node_deque )
   {
      if ( !deq.empty() )
      {
         for ( auto node = deq.begin(); node !=deq.end(); ++node )
         {
            lb_i = node->get_lowerbound();
            if ( min_lb > lb_i )
               min_lb = lb_i;
         }
      }
   }

   return min_lb;
}

bool NODELIST::check_size_LIST() const
{
   if ( node_list.empty() )
   {
      if ( listsize == 0 )
         return true;
      else
         return false;
   }
   else
   {
      if ( (int)node_list.size() == listsize )
         return true;
      else
         return false;
   }
}

bool NODELIST::check_size_TDEQUE() const
{
   if ( node_deque_1.empty() )
   {
      if ( node_deque_2.empty() )
      {
         if ( listsize == 0 )
            return true;
         else
            return false;
      }
      else
      {
         if ( (int)node_deque_2.size() == listsize )
            return true;
         else
            return false;
      }
   }
   else
   {
      // node_deque_1 is noempty
      if ( node_deque_2.empty() )
      {
         if ( (int)node_deque_1.size() == listsize )
            return true;
         else
            return false;
      }
      else
      {
         if ( (int)node_deque_2.size() + (int)node_deque_1.size() == listsize )
            return true;
         else
            return false;
      }
   }
}

bool NODELIST::check_size_TENDEQUE() const
{
   assert( type == TEN_DEQUE );
   assert( division > 0.0 );
   assert( maxnode > 0 );
   assert( !node_deque.empty() );
   assert( (int)node_deque.size() == 10 );

   int total = 0;

   for ( auto& deq : node_deque )
   {
      if ( !deq.empty() )
         total += (int) deq.size();
   }

   if ( total == listsize )
      return true;
   else
      return false;
}
