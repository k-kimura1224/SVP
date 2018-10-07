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

   copy( source.node_list.begin(), source.node_list.end(), back_inserter(node_list) );
   copy( source.node_deque_1.begin(), source.node_deque_1.end(), back_inserter(node_deque_1) );
   copy( source.node_deque_2.begin(), source.node_deque_2.end(), back_inserter(node_deque_2) );
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

      copy( source.node_list.begin(), source.node_list.end(), back_inserter(node_list) );
      copy( source.node_deque_1.begin(), source.node_deque_1.end(), back_inserter(node_deque_1) );
      copy( source.node_deque_2.begin(), source.node_deque_2.end(), back_inserter(node_deque_2) );
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
}

// setup
void NODELIST::setup(
      const TYPE_NODELIST  s_type,
      const double         bestval
      )
{
   type = s_type;

   switch ( type )
   {
      case TWO_DEQUE:
      {
         standardgap = 0.9;
         standardvalue = ( 1.0 - standardgap ) * bestval;

         break;
      }
      default:
      {
         break;
      }
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

   if ( node.get_lowerbound() < standardvalue )
      node_deque_1.push_back( node );
   else
      node_deque_2.push_back( node );

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

   if ( node.get_lowerbound() < standardvalue )
      node_deque_1.push_back( move( node ) );
   else
      node_deque_2.push_back( move( node ) );

   ++listsize;
}

// nodeselection
NODE& NODELIST::nodeselection_LIST(
      double*        globallowerbound,
      const double   bestval,
      const int      index,
      const int      disp
      )
{
   assert( (int)node_list.size() == listsize );
   assert( listsize > 0 );

   auto  buf = (index-1)%100000;

   if( ( buf == 0 || buf == 1 ) && disp < index )
   {
      node_list.sort();
      *globallowerbound = node_list.begin()->get_lowerbound();
   }

   return *(node_list.begin());
}

NODE& NODELIST::nodeselection_TDEQUE(
      double*        globallowerbound,
      const double   bestval,
      const int      index,
      const int      disp
      )
{
   assert( listsize > 0 );
   assert( !node_deque_1.empty() || !node_deque_2.empty() );

   if ( !node_deque_1.empty() )
      return node_deque_1.front();

   assert( standardvalue > 0.0 );
   assert( standardgap > 0.0 && standardgap < 1.0 );

   *globallowerbound = standardvalue;
   standardgap -= 0.1;
   standardvalue = ( 1.0 - standardgap ) * bestval;

   //cout << "node_selection - start" << endl;
   do {
      assert( !node_deque_2.empty() );
      int size = (int) node_deque_2.size();
      for ( int i = 0; i < size; )
      {
         //cout << "i:" << i << ", size2:" << node_deque_2.size();
         //cout << ", size1:" << node_deque_1.size() << endl;
         if ( node_deque_2[i].get_lowerbound() < standardvalue )
         {
            node_deque_1.push_back( move( node_deque_2[i] ) );
            node_deque_2[i] = move( node_deque_2.back() );
            node_deque_2.pop_back();
            --size;
            assert( size == 0 || size == (int) node_deque_2.size() );
            continue;
         }
         ++i;
      }
      //cout << "for-end" << endl;
      standardgap -= 0.1;
      standardvalue = ( 1.0 - standardgap ) * bestval;
   } while ( node_deque_1.empty() );
   //cout << "node_selection - end" << endl;

   return node_deque_1.front();
}

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

