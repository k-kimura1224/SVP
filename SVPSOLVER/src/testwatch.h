#ifndef TESTWATCH_H__
#define TESTWATCH_H__

#include <time.h>

using namespace std;


class TESTWATCH{
	clock_t	begin;
	clock_t	end;

	bool		run;

   double   total;

	public:

		TESTWATCH();											// default constructor
		TESTWATCH( const TESTWATCH &source );			// copy constructor
		TESTWATCH& operator=( const TESTWATCH& );		// assignment operator
		~TESTWATCH();											// destructor

		void	start();
		void	stop();
		double	get_result() const
      {
         assert( run == false );
         return total;
      }
		double	get_time() const;
};

#endif
