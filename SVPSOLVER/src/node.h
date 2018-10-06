#ifndef NODE_H__
#define NODE_H__

using namespace std;


class NODE{
	int		m;
	double	*ub;
	double	*lb;
	double	*warm;
	double	*sumfixed;
	double	relax_objval;
	double	*relax_solval;
	int		dpt;
	bool		zero;
	int		index;
	bool		solved;

	public:

		NODE();											// default constructor
		NODE( const NODE &source );			// copy constructor
		NODE& operator=( const NODE& );		// assignment operator
      NODE( NODE &&source ) noexcept;       // move constructor
		~NODE();											// destructor
		bool operator<(const NODE &rhs) const
		{ return relax_objval < rhs.relax_objval; }

		void	set_vals( int s_m, double *s_ub, double *s_lb,
								double *s_warm, double s_relax_objval,
								int s_dpt, bool s_zero, int s_index);

		double	get_lowerbound() const { return relax_objval; }
		int		get_index() const { return index; }
		int		get_dpt() const { return dpt; }
		bool		get_zero() const { return zero; }
		double*	get_ub() const { return ub; }
		double*	get_lb() const { return lb; }
		double*	get_warm() const { return warm; }
		double*	get_relaxsolval() const { return relax_solval; }
		bool		get_solved() const { return solved; }
		double*	get_sumfixed() const { return sumfixed; }

      void     NODEdispInformation();

		void	set_lowerbound( double s_ro ){
			if( relax_objval < s_ro ) relax_objval = s_ro;
		}
		void	set_relaxsolval( double *solval );
		void	set_solved( bool s_solved ){ solved = s_solved; }
		bool	alloc_sumfixed();
		void	set_sumfixed( double c, double *s_sumfixed );
		void	add_sumfixed( double c, double *s_sumfixed );

      void  set_lbval( int i, double s_lb )
      {
         assert( lb != nullptr );
         assert( i >= 0 && i < m );
         lb[i] = s_lb;
      }

      void  set_ubval( int i, double s_ub )
      {
         assert( ub != nullptr );
         assert( i >= 0 && i < m );
         ub[i] = s_ub;
      }

      void  set_dpt( int s_dpt )
      {
         assert( dpt >= 0 );
         dpt = s_dpt;
      }

      void  set_zero( int s_zero )
      {
         zero = s_zero;
      }

      void  set_index( int s_index )
      {
         assert( index >= 0 );
         index = s_index;
      }
};

#endif
