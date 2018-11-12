#ifndef NODELIST_H__
#define NODELIST_H__

#include <list>
#include <deque>
#include <vector>

#include "node.h"

using namespace std;

enum TYPE_NODELIST {
   LIST = 0,
   TWO_DEQUE = 1,
   TEN_DEQUE = 2
};

class NODELIST{
   TYPE_NODELIST  type;

   list<NODE>     node_list;

   deque<NODE>    node_deque_1;     // difficult nodes
   deque<NODE>    node_deque_2;
   double         standardvalue;
   double         standardgap;

   vector<deque<NODE>>  node_deque;
   double               division;
   int                  maxnode;
   int                  currentdeq;

   int            listsize;

	public:

		NODELIST();											// default constructor
		NODELIST( const NODELIST &source );			// copy constructor
		NODELIST& operator=( const NODELIST& );	// assignment operator
		~NODELIST();										// destructor

      void  setup( const TYPE_NODELIST s_type, const double bestval,
                      const int s_memory, const int m,
                      const double gap=0.9 );

      int   getListsize() const { return listsize; }

      int   getSubsize_LIST( const int sub ) const { return listsize; }
      int   getSubsize_TDEQUE( const int sub ) const;
      int   getSubsize_TENDEQUE( const int sub ) const {
         assert( sub >= 0 && sub < 10 );
         assert( !node_deque.empty() );
         if ( node_deque[sub].empty() )
            return 0;
         else
            return (int) node_deque[sub].size();
      }

      void  push_back_LIST( const NODE& node );
      void  push_back_TDEQUE( const NODE& node );
      void  push_back_TENDEQUE( const NODE& node )
      {
         assert( type == TEN_DEQUE );
         assert( division > 0.0 );
         assert( maxnode > 0 );
         int index = node.get_lowerbound() / division;
         assert( index >= 0 && index < 10 );
         node_deque[index].push_back( node );
         ++listsize;
      }

      void  move_back_LIST( NODE& node );
      void  move_back_TDEQUE( NODE& node );
      void  move_back_TENDEQUE( NODE& node )
      {
         assert( type == TEN_DEQUE );
         assert( division > 0.0 );
         assert( maxnode > 0 );
         int index = node.get_lowerbound() / division;
         if ( index >= 10 )
            index = 9;
         assert( index >= 0 && index < 10 );
         node_deque[index].push_back( move( node ) );
         ++listsize;
      }

      auto&  emplace_back_TENDEQUE(
            const int   dpt1,
            const int   nodeindex,
            const double lowerbound_new,
            const vector<int>& fixedvalues
      )
      {
         assert( type == TEN_DEQUE );
         assert( division > 0.0 );
         assert( maxnode > 0 );
         int index = lowerbound_new / division;
         if ( index >= 10 )
            index = 9;
         assert( index >= 0 && index < 10 );
         node_deque[index].emplace_back( dpt1, nodeindex, lowerbound_new, fixedvalues );
         ++listsize;
         return node_deque[index].back();
      }

      void  cutoff_LIST() {
         node_list.pop_front();
         --listsize;
      }
      void  cutoff_TDEQUE() {
         cout << "error?" << endl;
         assert(0);
         exit(-1);
         if ( !node_deque_1.empty() )
         {
            node_deque_1.pop_front();
         }
         else
         {
            assert( !node_deque_2.empty() );
            node_deque_2.pop_front();
         }

         --listsize;
      }
      void  cutoff_TENDEQUE() {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( division > 0.0 );
         assert( maxnode > 0 );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         assert( currentdeq >= 0 && currentdeq < 10 );
         if ( maxnode > listsize ) {
            for ( auto i = currentdeq; i < 10; ++i ) {
               if ( !node_deque[i].empty() ) {
                  node_deque[i].pop_front();
                  break;
               }
            }
         } else {
            for ( auto i = 9; i >= currentdeq; --i ) {
               if ( !node_deque[i].empty() ) {
                  node_deque[i].pop_front();
                  break;
               }
            }
         }
         --listsize;
      }

      NODE*    nodeselection_LIST( double* globallowerbound, const double bestval, const int index, const int disp );
      NODE*    nodeselection_TDEQUE( double* globallowerbound, const double bestval, const int index, const int disp );
      NODE*    nodeselection_TENDEQUE( double* globallowerbound, const double bestval, const int index, const int disp );

      double   get_GLB_LIST() const;
      double   get_GLB_TDEQUE() const;
      double   get_GLB_TENDEQUE() const;

      bool     check_size_LIST() const;
      bool     check_size_TDEQUE() const;
      bool     check_size_TENDEQUE() const;

      int      setup_para_selection_LIST() const {
         assert( listsize > 0 && !node_list.empty() );
         return 0;
      }
      int      setup_para_selection_TDEQUE() const {
         assert( listsize > 0 &&
                  ( !node_deque_1.empty() || !node_deque_2.empty() ) );
         if ( !node_deque_1.empty() && !node_deque_2.empty() )
         {
            if ( node_deque_1.size() > node_deque_2.size() )
               return 1;
            else
               return 2;
         }
         else
         {
            if ( node_deque_1.empty() )
               return 2;
            else
               return 1;
         }
      }
      int      setup_para_selection_TENDEQUE() const {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         int maxsize = 0;
         int size;
         int index = 0;
         int maxindex = -1;

         for ( auto& deq : node_deque )
         {
            if ( !deq.empty() )
            {
               size = (int) deq.size();
               if ( size > maxsize )
               {
                  maxsize = size;
                  maxindex = index;
               }
            }
            ++index;
         }
         assert( maxindex >= 0 && maxindex < 10 );
         return maxindex;
      }

      int      setup_parapush_selection_TENDEQUE() const {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         assert( maxnode > 0 );

         if ( maxnode > listsize )
         {
            for ( int i = 0; i < 10; ++i)
            {
               if ( !node_deque[i].empty() )
                  return i;
            }
         }
         else
         {
            for ( int i = 9; i >= 0; --i )
               if ( !node_deque[i].empty() )
                  return i;
         }
         return -1;
      }
      int      setup_parapop_selection_TENDEQUE( const bool sleep ) const {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         assert( maxnode > 0 );

         if ( sleep || maxnode < listsize )
         {
            for ( int i = 0; i < 10; ++i)
            {
               if ( !node_deque[i].empty() )
                  return i;
            }
         }
         else
         {
            for ( int i = 9; i >= 0; --i )
            {
               if ( !node_deque[i].empty() )
                  return i;
            }
         }
         return -1;
      }

      NODE&    para_selection_LIST( const int setup ) {
         assert( listsize > 0 && !node_list.empty() );
         return node_list.front(); }
      NODE&    para_selection_TDEQUE( const int setup ) {
         assert( listsize > 0 );
         assert( setup == 1 || setup == 2 );
         if ( setup == 1 )
         {
            assert( !node_deque_1.empty() );
            return node_deque_1.front();
         }
         else
         {
            assert( !node_deque_2.empty() );
            return node_deque_2.front();
         }
      }
      NODE&    para_selection_TENDEQUE( const int setup ) {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         assert( setup >= 0 && setup < 10 );
         assert( !node_deque[setup].empty() );
         return node_deque[setup].front();
      }

      void     pop_front_LIST( const int setup ) {
         node_list.pop_front();
         --listsize;
      }
      void     pop_front_TDEQUE( const int setup ) {
         assert( setup == 1 || setup == 2 );
         if ( setup == 1 )
            node_deque_1.pop_front();
         else
            node_deque_2.pop_front();
         --listsize;
      }
      void     pop_front_TENDEQUE( const int setup ) {
         assert( listsize > 0 );
         assert( type == TEN_DEQUE );
         assert( !node_deque.empty() );
         assert( (int)node_deque.size() == 10 );
         assert( setup >= 0 && setup < 10 );
         assert( !node_deque[setup].empty() );
         node_deque[setup].pop_front();
         --listsize;
      }
      void     sort() { node_list.sort(); }
};

#endif
