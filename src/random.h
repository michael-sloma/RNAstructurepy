#if !defined(RANDOM_H)
#define RANDOM_H
#define TABSIZE 32

class randomnumber {
	private:
		long IMM1,DIVISOR,/*RNMX,*/idum2,idum,iy,iv[TABSIZE];
		double IM1INV;

	public:
		/*	make a random number generator	*/
		randomnumber();
		/*	seed the random number generator	*/
		void seed(long seeddouble);
		/*	generate a random number between 0.0 and 1.0	*/
		double roll();
		/* generate a random integer in a range 
		   Seems to work correctly BUT I make no 
		   guarantees about rigorously correct behavior 
		   MS 11/2014 */
		int roll_int(int min, int max);
};


#endif
