#ifndef _TURBOFOLD_THREAD_
#define _TURBOFOLD_THREAD_

/*
Encapsulates the thread that refolds a sequence.
*/
class TurboFold;
class t_ansi_thread;

class t_turbofold_thread
{
public:
	t_turbofold_thread(TurboFold* _turbofold_obj, int _i_seq, int _i_iter);
	~t_turbofold_thread();

	TurboFold* turbofold_obj;
	int i_seq;
	int i_iter;

#ifdef COMPILE_SMP 
	t_ansi_thread* computation_thread;
#endif // TURBOFOLD_SMP 

	// Start the thread that refolds the sequence.
	int run_refolding_thread();

	// Threading function: Use number of threads specified by iterative_decoder per pfunction computation.
	static void* refold(void* _thread_arg);
};

#endif // _TURBOFOLD_THREAD_
