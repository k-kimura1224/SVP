#ifndef NODELIST_H__
#define NODELIST_H__

#include <list>
#include <deque>
#include <vector>

#include "node.h"

using namespace std;

enum TYPE_NODELIST {
   LIST = 0,
   TWO_DEQUE = 1
};

class NODELIST{
   TYPE_NODELIST  type;

   list<NODE>     node_list;

   deque<NODE>    node_deque_1;     // difficult nodes
   deque<NODE>    node_deque_2;
   double         standardvalue;
   double         standardgap;

   int            listsize;

	public:

		NODELIST();											// default constructor
		NODELIST( const NODELIST &source );			// copy constructor
		NODELIST& operator=( const NODELIST& );	// assignment operator
		~NODELIST();										// destructor

      void     setup( const TYPE_NODELIST s_type, const double bestval,
                      const double gap=0.9 );

      int      getListsize() const { return listsize; }

      int      getSubsize_LIST( const int sub ) const { return listsize; }
      int      getSubsize_TDEQUE( const int sub ) const {
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

      void     push_back_LIST( const NODE& node );
      void     push_back_TDEQUE( const NODE& node );

      void     move_back_LIST( NODE& node );
      void     move_back_TDEQUE( NODE& node );

      void     cutoff_LIST() {
         node_list.pop_front();
         --listsize;
      }
      void     cutoff_TDEQUE() {
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


      NODE&    nodeselection_LIST( double* globallowerbound, const double bestval, const int index, const int disp );
      NODE&    nodeselection_TDEQUE( double* globallowerbound, const double bestval, const int index, const int disp );

      double   get_GLB_LIST() const;
      double   get_GLB_TDEQUE() const;

      bool     check_size_LIST() const;
      bool     check_size_TDEQUE() const;

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

      void     sort() { node_list.sort(); }
};

#endif
