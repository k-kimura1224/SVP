#ifndef STOPWATCH_H__
#define STOPWATCH_H__

#include <time.h>

using namespace std;


class STOPWATCH{
	time_t	begin;
	time_t	end;

	bool		run;

   int      total;
	int		limit;

	public:

		STOPWATCH();											// default constructor
		STOPWATCH( const STOPWATCH &source );			// copy constructor
		STOPWATCH& operator=( const STOPWATCH& );		// assignment operator
		~STOPWATCH();											// destructor

		void	start();
		void	stop();
		void	set_timelimit(int l){ limit = l; }
		int	get_result()
      {
         assert( run == false );
         return total;
      }
		int	get_time() const;
		bool	check_time();
};

#endif
