#include <stdio.h>
#include <stdlib.h>
#include "../RNA_class/RNA.h"

#ifdef COMPILE_SMP
#include "../src/phmm/utils/ansi_thread/ansi_thread.h"
#include "../src/phmm/utils/mutex/ansi_mutex.h"
#endif

#include "../src/phmm/structure/structure_object.h"
#include "../src/phmm/structure/folding_constraints.h"
#include "TurboFold_thread.h"
#include "TurboFold_object.h"

#include "../src/phmm/utils/xmath/matrix/matrix.h"
#include <math.h>

t_turbofold_thread::t_turbofold_thread(TurboFold* _turbofold_obj, int _i_seq, int _i_iter)
{
	this->turbofold_obj = _turbofold_obj;
	this->i_seq = _i_seq;
	this->i_iter = _i_iter;
#ifdef COMPILE_SMP 
	this->computation_thread = NULL;
#endif // COMPILE_SMP 
}

t_turbofold_thread::~t_turbofold_thread()
{
#ifdef COMPILE_SMP 
	if(this->computation_thread != NULL)
	{
		delete(this->computation_thread);
	}
#endif // COMPILE_SMP 
}

#ifdef COMPILE_SMP 
// Start the thread that refolds the sequence.
int t_turbofold_thread::run_refolding_thread()
{
	this->computation_thread = new t_ansi_thread(t_turbofold_thread::refold, this);

	// Start the thread.
	if(this->computation_thread->run_thread())
	{
		return(0);
	}
	else
	{
/*
		printf("Could not start refolding thread for %s @ %s(%d).\n", 
			this->turbofold_obj->sequences->at(this->i_seq)->ctlabel,
			__FILE__, __LINE__);
*/
		return(1);
	}
}

// Threading function: Use number of threads specified by iterative_decoder per pfunction computation.
void* t_turbofold_thread::refold(void* _thread_arg)
{
	t_turbofold_thread* turbofold_thread = (t_turbofold_thread*)_thread_arg;	

	// Setup the variables.
        TurboFold* turbofold_obj = turbofold_thread->turbofold_obj;
	int i_seq = turbofold_thread->i_seq;
	vector<t_structure*>* sequences = turbofold_obj->sequences;
	int i_iter = turbofold_thread->i_iter;
	int n_iterations = turbofold_obj->n_iterations;

	//printf("Refolding %s with thread %d\n", sequences->at(i_seq)->ctlabel, turbofold_thread);

	// Use the supplied save file name if this is the last iteartion. 
	if(i_iter == n_iterations &&
		turbofold_obj->saves != NULL)
	{
		turbofold_obj->folders->at(i_seq)->PartitionFunction(turbofold_obj->saves->at(i_seq));
	}
	else
	{
		turbofold_obj->folders->at(i_seq)->PartitionFunction();
	}

	turbofold_thread->computation_thread->thread_state = THREAD_TERMINAL;

	// Start a new child thread: This is how the computations continue by themselves. There is not infinite loop that checks if a new sequence can be
	// refolded or not. Also it is vital to call this after setting the state to THREAD_TERMINAL so that a place for a new thread is open.
	turbofold_obj->start_next_refolding_thread(false, i_iter);

	return(NULL);
}
#endif

