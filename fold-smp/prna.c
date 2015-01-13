#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "base.h"
#include "prna.h"
#include "util.h"
#include "param.h"

//#define DISABLE_COAXIAL

struct prna
{
  int n; /* number of bases */
  base_t *seq;  /* sequence */
  int_t *v; /* n x n array */
  int_t *w, *wm; /* n x n arrays */
  int_t *w5, *w3; /* n elements */
};

/* penalty for a helix terminated by a pair containing a U */
DEV static int_t terminal_U_penalty(const base_t *s, const int i, const int j, param_t p)
{
  return s[i] == U || s[j] == U ? p->terminal_AU_penalty : 0.0;
}

DEV static int_t dangle_3p_energy(const base_t *s,
				   const int i,
				   const int j,
				   const int ip1,
                                   param_t p)
{
  return p->dangle_3p[s[i]][s[j]][s[ip1]] + terminal_U_penalty(s,i,j,p);
}

DEV static int_t dangle_5p_energy(const base_t *s,
				   const int i,
				   const int j,
				   const int jm1,
                                   param_t p)
{
  return p->dangle_5p[s[i]][s[j]][s[jm1]] + terminal_U_penalty(s,i,j,p);
}

DEV static int_t terminal_stack(const base_t *s,
		                 const int i,
				 const int j,
				 const int ip1,
				 const int jm1,
                                 param_t p)
{
  return p->tstack[s[i]][s[j]][s[ip1]][s[jm1]] + terminal_U_penalty(s,i,j,p);
}

DEV static int_t terminal_stack_multibranch(const base_t *s,
					     const int i,
					     const int j,
					     const int ip1,
					     const int jm1,
					     param_t p)
{
  return p->tstackm[s[i]][s[j]][s[ip1]][s[jm1]] + terminal_U_penalty(s,i,j,p);
}


DEV static const int_t *lookup_find(const base_t *s, const int d, param_t p)
{
  int i;
  switch (d) {
  case 3:
    for (i = 0; i < p->ntriloop; i++)
      if (sequences_match(s, p->triloop[i].seq, d+2))
	return &p->triloop[i].val;
    break;
  case 4:
    for (i = 0; i < p->ntloop; i++)
      if (sequences_match(s, p->tloop[i].seq, d+2))
	return &p->tloop[i].val;
    break;
  case 6:
    for (i = 0; i < p->nhexaloop; i++)
      if (sequences_match(s, p->hexaloop[i].seq, d+2))
	return &p->hexaloop[i].val;
    break;
  }
  return 0;
}

/***
 * Energy of a hairpin loop with d unpaired bases, d = j-i-1
 * s[i] is paired with s[j]
 * s[i+1] is mismatched with s[j-1]
 ***/
DEV static int_t hairpin_loop_energy(const base_t *s, 
				      const int i, 
				      const int j, 
				      const int d,
                                      param_t p)
{
  /* Lookup tables for special hairpin loops */
  const int_t *val;
  if ((val = lookup_find(&s[i],d,p)))
    return *val;
  
  /* Hairpin loop initiation penalty */
  int_t e;
  if (d > LOOP_MAX)
    e = (int_t) (p->hairpin_loop_initiation[LOOP_MAX] + p->prelog * 
      LOG((float) d / LOOP_MAX));
//    e = (int_t) (p->hairpin_loop_initiation[LOOP_MAX] + p->Extrapolation_for_large_loops * 
//      LOG((float) d / LOOP_MAX));
  else
    e = p->hairpin_loop_initiation[d];
  
  if (d == 3) {
    if (contains_only_base(C,d,&s[i+1]))
      e += p->c_hairpin_of_3;
    e += terminal_U_penalty(s,i,j,p);
  } else {
    e += p->tstackh[s[i]][s[j]][s[i+1]][s[j-1]];
    if (contains_only_base(C,d,&s[i+1]))
      e += p->c_hairpin_slope*d + p->c_hairpin_intercept;
  }
  
  if (s[i] == G && s[j] == U && i > 1 && s[i-1] == G && s[i-2] == G)
    e += p->bonus_for_GGG_hairpin;
  
  return e;
  
}

DEV static int_t real_min(int_t a, int_t b) { return a < b ? a : b; }

/***
 * Energy of an internal/bulge loop with d1, d2 unpaired bases,
 *   d1 = ip-i-1,  d2 = j-jp-1
 * s[i] is paired with s[j]
 * s[i+1] is mismatched sith s[j-1]
 * s[ip-1] is mismatched with s[jp+1]
 * s[ip] is paired with s[jp]
 ***/

DEV static int_t alternative_bulge_loop_correction (const int n, const base_t *s, 
                                                  const int i, 
                                                  const int ip) //i<ip
{
  int count = 1;
  int k;
  //float result;
  if (i!=n-1){
    k = i;
    while (k>=0 && s[k]==s[i+1]) {
      count++;
      k--;   
    }
  
    k = ip;
    while (k<=n-1 && (s[k]==s[i+1])) {
      count++;
      k++;
    }
  }
  return (int_t) (-1.0f * RT * conversion_factor * log ((float) count));

//  printf("count %d, log count %f, internal loop correction %f\n",count,log((float)count),result);
//  return (int_t) result;
}

DEV static int_t internal_loop_energy(const base_t *s,
                                       const int n,
				       const int i,
				       const int j,
				       const int ip,
				       const int jp,
				       const int d1,
				       const int d2,
                                       param_t p)
{
  /* Bulge loops */
  if (d1 == 0 || d2 == 0) {
    int_t e = p->bulge_loop_initiation[d1+d2]; 
    if (d1 == 1 || d2 == 1) { /* single-nucleotide bulge */
      e += p->stack[s[i]][s[j]][s[ip]][s[jp]];
        if (d1==0) e += alternative_bulge_loop_correction(n,s,jp,j); //correction for multiple equivalent bulge loops
        //else e += alternative_bulge_loop_correction(s,i,jp);
        else e += alternative_bulge_loop_correction(n,s,i,ip);
      if ((d1 == 1 && s[i+1] == C && (s[i] == C || s[i+2] == C)) ||
          (d2 == 1 && s[j-1] == C && (s[j] == C || s[j-2] == C)))
        e += p->Bonus_for_Single_C_bulges_adjacent_to_C;
    } else {
      e += terminal_U_penalty(s,i,j,p);
      e += terminal_U_penalty(s,ip,jp,p);
    }
    return e;
  } 
  
  /* Small internal loops */
  if (d1 == 1 && d2 == 1)
    return p->int11[s[i]][s[i+1]][s[i+2]][s[j-2]][s[j-1]][s[j]];
  if (d1 == 2 && d2 == 2)
    return p->int22[s[i]][s[ip]][s[j]][s[jp]][s[i+1]][s[i+2]][s[j-1]][s[j-2]];
  if (d1 == 1 && d2 == 2)
    return p->int21[s[i]][s[j]][s[i+1]][s[j-1]][s[jp+1]][s[ip]][s[jp]];
  if (d1 == 2 && d2 == 1)
    return p->int21[s[jp]][s[ip]][s[jp+1]][s[ip-1]][s[i+1]][s[j]][s[i]];
  
  /* Larger internal loops */
  tab4_t *sp;
  if (d1 == 1 || d2 == 1)
    sp = &p->tstacki1n;
  else if ((d1 == 2 && d2 == 3) || (d1 == 3 && d2 == 2))
    sp = &p->tstacki23;
  else
    sp = &p->tstacki;
  return p->internal_loop_initiation[d1+d2] + 
    real_min(p->fm_array_first_element*abs(d1-d2), p->maximum_correction) +
    (*sp)[s[i]][s[j]][s[i+1]][s[j-1]] +
    (*sp)[s[jp]][s[ip]][s[jp+1]][s[ip-1]];
  
}

#ifndef disablecoax
DEV static int_t coaxial_flush(const base_t *s,
				const int i,
				const int j,
				const int ip,
				const int jp,
                                param_t p)
{
  return terminal_U_penalty(s,i,j,p) + terminal_U_penalty(s,ip,jp,p) +
    p->coaxial[s[i]][s[j]][s[ip]][s[jp]];
}

DEV static int_t coaxial_mismatch1(const base_t *s,
				    const int i,
				    const int j,
				    const int ip,
				    const int jp,
                                    param_t p)
{
  return terminal_U_penalty(s,i,j,p) + terminal_U_penalty(s,ip,jp,p) +
    p->tstackcoax[s[j]][s[i]][s[j+1]][s[i-1]] +
    p->coaxstack[s[j+1]][s[i-1]][s[ip]][s[jp]];
}

DEV static int_t coaxial_mismatch2(const base_t *s,
				    const int i,
				    const int j,
				    const int ip,
				    const int jp,
                                    param_t p)
{
  return terminal_U_penalty(s,i,j,p) + terminal_U_penalty(s,ip,jp,p) +
    p->tstackcoax[s[jp]][s[ip]][s[jp+1]][s[ip-1]] +
    p->coaxstack[s[j]][s[i]][s[j+1]][s[jp+1]];
} 
#endif//disablecoax

/* return -ln(e^-a + e^-b) */
DEV static int_t free_energy_sum(const int_t a, const int_t b)
{
  if (a < b)
    return a - LOG1P(EXP(a-b));
  else if (b < a)
    return b - LOG1P(EXP(b-a));
  else
    return a - LOG(2);
}

DEV static void free_energy_accumulate(int_t *a, const int_t b)
{
  *a = free_energy_sum(*a,b);
}

DEV static void free_energy_min(int_t *a, const int_t b)
{
 if(*a>b) *a = b; 
}

DEV HOST static int int_min(int a, int b) { return a < b ? a : b; }

DEV HOST static int ind(int i, int j, int n) 
{ 
  return i*n + j;
}

DEV static int upper_triangle_index(int i, int j)
{
  return (j*(j-1))/2 + i;
}

DEV HOST inline static int cp(int i, int j, const base_t *s)
{
  return j-i-1 >= LOOP_MIN && is_canonical_pair(s[i],s[j]);
}
				       
DEV HOST inline static int cp2(int i, int j, const base_t *s)
{
  return is_canonical_pair(s[i],s[j]);
}

DEV HOST inline static int can_pair(int i, int j, int n, const base_t *s)
{
  if (j < i) {
    const int tmp = i;
    i = j;
    j = tmp;
  }
  return cp(i,j,s) && ((i > 0 && j < n-1 && cp(i-1,j+1,s)) || cp(i+1,j-1,s));
}

DEV HOST inline static int not_isolated(int i,int j,int n, const base_t *s)
{
  if (j < i) {
    const int tmp = i;
    i = j;
    j = tmp;
  }
  return is_canonical_pair(s[i],s[j]) && ((i > 0 && j < n-1 && cp(i-1,j+1,s)) || cp(i+1,j-1,s));
}

DEV static int wrap(int i, int n)
{
  return i >= n ? i-n : i;
}

DEV static int is_exterior(int i, int j)
{
  return j < i;
}

DEV static int is_interior(int i, int j)
{
  return i < j;
}

DEV HOST static int_t *array_val(int_t *a, int i, int j, int n, const base_t *s)
{
  return can_pair(i,j,n,s) ? &a[ind(i,j,n)] : 0;
//  return &a[ind(i,j,n)];
}

#ifdef __CUDACC__
#define ISTART blockIdx.x
#define IINC gridDim.x
#else
#define ISTART 0
#define IINC 1
#endif

//MFE recursions begin
//TODO
//figure out source of differences in arrays
//integrate with rnastructure traceback


//when recursions work on the cpu:
//do the same thing with the calculation on the GPU
GLOBAL static void calc_V_hairpin_and_V_stack//_multibranch//calculate V(i,j)=Vhairpin+Vstack+Vmultibranch but NOT bulge/internal
(int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{
  int i;
  for (i = ISTART; i < n; i += IINC) { //for(i=blockId.x;i<numberofbases;i+=gridDim.x) so each thread will handle 0+blockId.x,256+blockId.x.. 
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
//if((i==0)&&(j==4)) printf("i and j are exterior: %d\ni and j can't pair: %d\n",is_exterior(i,j) && i-j<=LOOP_MIN, !can_pair(i,j,n,s));
    //if ((is_exterior(i,j) && i-j <= LOOP_MIN) || !can_pair(i,j,n,s)){
//    if(is_exterior(i,j) || !can_pair(i,j,n,s)){
    if((is_interior(i,j) && !can_pair(i,j,n,s)) || (is_exterior(i,j) && (!is_canonical_pair(s[i],s[j]) ))){
      v[ind(i,j,n)] = INF; //this probably isn't necessary(wrong, it is necessary)
      continue;
    }
//printf("looking for a hairpin or stack\n");
    int_t vij = INF; //temp variable to fold free energy sum
    if (i != n-1 && j != 0) { 
      /* hairpin loop */
      if (is_interior(i,j))
        vij = hairpin_loop_energy(s,i,j,d,p); 
      /* stack */
      //if (can_pair(i+1,j-1,n,s) && !((is_interior(i,j)) && (d <= LOOP_MIN-2)))//-2???
      if (can_pair(i+1,j-1,n,s) && !((is_interior(i,j)) && (d <= LOOP_MIN-2)))//-2???
//         printf("we are checking stacks BEFORE. i %d j %d Vij %d hairpin energy %d V(i+1,j-1)%d\n", i,j,vij,p->stack[s[i]][s[j]][s[i+1]][s[j-1]],v[ind(i+1,j-1,n)]);
        free_energy_min(&vij, p->stack[s[i]][s[j]][s[i+1]][s[j-1]] + v[ind(i+1,j-1,n)]);
//         printf("we are checking stacks AFTER. i %d j %d Vij %d hairpin energy %d V(i+1,j-1)%d\n", i,j,vij,p->stack[s[i]][s[j]][s[i+1]][s[j-1]],v[ind(i+1,j-1,n)]);
    }



/*	//old code for multibranch loops
    if (d > 2*LOOP_MIN + 3 && i != n-1 && j != 0) { //if i and j are far enough apart to close a MBL..
      free_energy_accumulate(&vij, x[ind((d-2)%5,i+1,n)] + terminal_U_penalty(s,i,j,p) + p->a + p->c); //have to change this, uses X
      if (i != n-2)
        free_energy_accumulate(&vij, x[ind((d-3)%5,i+2,n)] + dangle_3p_energy(s,i,j,i+1,p) + p->a + p->b + p->c);
      if (j != 1)
        free_energy_accumulate(&vij, x[ind((d-3)%5,i+1,n)] + dangle_5p_energy(s,i,j,j-1,p) + p->a + p->b + p->c);
      if (i != n-2 && j != 1)
        free_energy_accumulate(&vij, x[ind((d-4)%5,i+2,n)] + terminal_stack_multibranch(s,i,j,i+1,j-1,p) + p->a + 2*p->b + p->c);//this leaves out the coaxial stacking calculation which is done later
    }
*/
    v[ind(i,j,n)] = vij;
  }
}



#ifdef __CUDACC__

#define NTHREAD 256
#define SQRT_NTHREAD 16

DEV static void free_energy_min_reduce(int_t *x, int tid, int nt)
{
  __shared__ int_t buf[NTHREAD];
  buf[tid] = *x;
  for (nt /= 2, __syncthreads(); nt > 0; nt /= 2, __syncthreads())
    if (tid < nt)
      free_energy_min(&buf[tid], buf[tid+nt]);
  if (tid == 0)
    *x = buf[0];
}

#endif /* __CUDACC__ */
GLOBAL static void calc_V_bulge_internal (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wl, int_t *w5, int_t *w3, param_t p)
{
//	Vbi(i,j) = min[V(k,l)+ Ebulge/int(i,j,k,l)] where i<k<l<j, i!=i+1, and j!=j-1
  int i;
  for (i = ISTART; i < n; i += IINC) {
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    if ((is_exterior(i,j) && i-j <= LOOP_MIN) ||
        (is_interior(i,j) && d <= LOOP_MIN+2) ||
        !can_pair(i,j,n,s))
      continue;
    int_t vij = INF;
#ifdef __CUDACC__
    const int d1start = threadIdx.x;
    const int d1inc = blockDim.x;
#else
    const int d1start = 0;
    const int d1inc = 1;
#endif
    const int dmax = int_min(LOOP_MAX, d-2);
    const int d1max = int_min(dmax, n-i-2);
    int d1;
    for (d1 = d1start; d1 <= d1max; d1 += d1inc) { //d1start is threadid, d1max is max loop size
      const int ip = i+d1+1; //ip depends on thread's ID in x dimension
      const int d2max = int_min(dmax-d1, j-1);
#ifdef __CUDACC__
      const int d2start = d1 > 0 ? threadIdx.y : threadIdx.y + 1;
      const int d2inc = blockDim.y;
#else
      const int d2start = d1 > 0 ? 0 : 1;
      const int d2inc = 1;
#endif
      int d2;
      for (d2 = d2start; d2 <= d2max; d2 += d2inc) {
        const int jp = j-d2-1;//jp depends on thread's ID in the y dimension
        if (can_pair(ip,jp,n,s))
          free_energy_min(&vij, internal_loop_energy(s,n,i,j,ip,jp,d1,d2,p) + v[ind(ip,jp,n)]);
      }
    }
#ifdef __CUDACC__
    const int tid = threadIdx.x * blockDim.y + threadIdx.y;
    free_energy_min_reduce(&vij, tid, blockDim.x*blockDim.y); //after we have 1 value per thread, do parallel reduction
    if (tid != 0)
      continue;
#endif
    free_energy_min(&v[ind(i,j,n)], vij); //write vij to V
  }
		
}


GLOBAL static void calc_V_multibranch (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{

//	Vmb(i,j) = min[WM(i+1,j-1)+c+a, WM(i+2,j-1)+Edangle5'+a+b+c, WM(i+1,j-2)+Edangle3'+a+b+c, WM(i+2,j-2)+Edangleboth+a+2b+c, 
//		min_over_k[ V(i+1,k) + min[W(k+1,j-1), WM(k+1,j-1)]] + a+2c+Eflushcoax(i to j, i+1 to k) , //various coaxial stacking possibilities
//		min_over_k[ V(k,j-1) + min[W(i+1,k-1), WM(i+1,k-1)]] + a+2c+Eflushcoax(i to j, k to j-1) , 
//		min_over_k[ V(i+2,k) + min[W(k+2,j-1), WM(k+2,j-1)]] + a+2c+2b+Emismatch3'coax(i to j, i+2 to k) , 
//		min_over_k[ V(i+2,k) + min[W(k+1,j-2), WM(k+1,j-2)]] + a+2c+2b+Emismatch5'coax(i to j, i+2 to k) , 
//		min_over_k[ V(k,j-2) + min[W(i+2,k-1), WM(i+2,k-1)]] + a+2c+2b+Emismatch3'coax(i to j, k to j-2) ,
//		min_over_k[ V(k,j-2) + min[W(i+1,k-2), WM(i+1,k-2)]] + a+2c+2b+Emismatch5'coax(i to j, k to j-2) ] 

//	where i < k < j


	//V(i,j) = min(V(i,j), Vmb(i,j))

  int i;
  for (i = ISTART; i < n; i += IINC) { //for(i=blockId.x;i<numberofbases;i+=gridDim.x) so each thread will handle 0+blockId.x,256+blockId.x.. 
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    if ((is_exterior(i,j) && i-j <= LOOP_MIN) || !can_pair(i,j,n,s))
      continue;
    int_t vij=INF;
/*
//I may need to put i and j backwards in ind
    if (d > 2*LOOP_MIN + 3 && i != n-1 && j != 0) { //if i and j are far enough apart to close a MBL..
        free_energy_min(&vij, wm[ind(i+1,j-1,n)] + terminal_U_penalty(s,i,j,p) + p->a + p->c); 
        if (i != n-2)
          free_energy_min(&vij, wm[ind(i+1,j-2,n)] + dangle_3p_energy(s,i,j,i+1,p) + p->a + p->b + p->c);
        if (j != 1)
          free_energy_min(&vij, wm[ind(i+2,j-1,n)] + dangle_5p_energy(s,i,j,j-1,p) + p->a + p->b + p->c);
        if (i != n-2 && j != 1)
          free_energy_min(&vij, wm[ind(i+2,j-2,n)] + terminal_stack_multibranch(s,i,j,i+1,j-1,p) + p->a + 2*p->b + p->c);
    }//old code possibly wrong
*//*
    if (d > 2*LOOP_MIN + 3 && i != n-1 && j != 0) { //if i and j are far enough apart to close a MBL..
        free_energy_min(&vij, wm[ind(i+1,j-1,n)] + terminal_U_penalty(s,i,j,p) + p->a + p->c); 
        if (i != n-2)
          free_energy_min(&vij, wm[ind(i+2,j-1,n)] + dangle_3p_energy(s,j,i+1,i,p) + p->a + p->b + p->c);
        if (j != 1)
          free_energy_min(&vij, wm[ind(i+1,j-2,n)] + dangle_5p_energy(s,j-1,i,j,p) + p->a + p->b + p->c);
        if (i != n-1 && j != 1)
          free_energy_min(&vij, wm[ind(i+2,j-2,n)] + terminal_stack_multibranch(s,j-1,i+1,j,i,p) + p->a + 2*p->b + p->c);
*/


    if (d > 2*LOOP_MIN + 3 && i != n-1 && j != 0) { //if i and j are far enough apart to close a MBL..
        free_energy_min(&vij, wm[ind(i+1,j-1,n)] + terminal_U_penalty(s,i,j,p) + p->a + p->c); 
        if (i != n-2)
          free_energy_min(&vij, wm[ind(i+2,j-1,n)] + dangle_3p_energy(s,i,j,i+1,p) + p->a + p->b + p->c);
        if (j != 1)
          free_energy_min(&vij, wm[ind(i+1,j-2,n)] + dangle_5p_energy(s,i,j,j-1,p) + p->a + p->b + p->c);
        if (i != n-2 && j != 1)
          free_energy_min(&vij, wm[ind(i+2,j-2,n)] + terminal_stack_multibranch(s,i,j,i+1,j-1,p) + p->a + 2*p->b + p->c);
    }
  free_energy_min(&v[ind(i,j,n)], vij);
  }
}


GLOBAL static void calc_V_exterior (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wl, int_t *w5, int_t *w3, param_t p)
{
//	Vexterior(i,j) = min[ W3(i+1)+W3(j-1-N), W3(i+2)+W5(j-1-N)+E5'dangle, W3(i+1)+W5(j-2-N)+E3'dangle, W3(i+2)+W5(j-2-N)+Emismatch,
//		min_over_k[ V(i+1,k) + W3(k+1) + W5(j-1-N) + Eflushcoax ],
//		min_over_k[ V(k,j-1-N) + W3(i+1) + W5(k-1) + E ],
//		min_over_k[ V(i+2,k-2) + W3(k+1) + W5(j-1-N) + E ],
//		min_over_k[ V(i+2,k-1) + W3(k+1) + W5(j-2-N) + E ],
//		min_over_k[ V(k+1,j-2-N) + W3(i+1) + W5(k-1) + E ],
//		min_over_k[ V(k,j-2-N) + W3(i+2) + W5(k-1) + E ] ]

  int i;
  for (i = ISTART; i < n; i += IINC) { //for(i=blockId.x;i<numberofbases;i+=gridDim.x) so each thread will handle 0+blockId.x,256+blockId.x.. 
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
//    if ( !is_interior(i,j) || (i-j <= LOOP_MIN) || !can_pair(i,j,n,s))
    if ( is_interior(i,j))
    //if ( is_interior(i,j) )//wrong answer
      continue;
    int_t vij = INF; //temp variable to fold free energy sum
//    if (i != n-1 && j != 0) {
//this order of args is right
   // if(can_pair(i,j,n,s) ||( (is_canonical_pair(s[i],s[j]))&&(i-j<=LOOP_MIN) ) ){
    //if(can_pair(i,j,n,s)){
    if(is_canonical_pair(s[i],s[j])&&not_isolated(i,j,n,s)){
      free_energy_min(&vij, w3[i+1] + w5[j-1] + terminal_U_penalty(s,i,j,p));
      if (i != n-1)
        free_energy_min(&vij, w3[i+2] + w5[j-1] + dangle_3p_energy(s,i,j,i+1,p));
      if (j != 0)
        free_energy_min(&vij, w3[i+1] + w5[j-2] + dangle_5p_energy(s,i,j,j-1,p));
      if (i != n-1 && j != 0)
        free_energy_min(&vij, w3[i+2] + w5[j-2] + terminal_stack(s,i,j,i+1,j-1,p));
    }
/*
    free_energy_min(&vij, w3[i+1] + w5[j-1] + terminal_U_penalty(s,i,j,p));
    if (i != n-1)
      free_energy_min(&vij, w3[i+2] + w5[j-1] + dangle_3p_energy(s,j,i+1,i,p));
    if (j != 0)
      free_energy_min(&vij, w3[i+1] + w5[j-2] + dangle_5p_energy(s,j-1,i,j,p));
    if (i != n-1 && j != 0)
      free_energy_min(&vij, w3[i+2] + w5[j-2] + terminal_stack_multibranch(s,j-1,i+1,j,i,p));
*/
    free_energy_min(&v[ind(i,j,n)], vij);
    }
#ifndef disablecoax
//placeholder for coaxial stacking recursions
#endif //disablecoax
}

GLOBAL static void calc_W (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wl, int_t *w5, int_t *w3, param_t p)
{
	//W(i,j) = min[V(i,j)+c,V(i+1,j)+Edangle5',  
//			V(i,j+1)+Edangle3',  
//			V(i+1,j+1)+Edangleboth]

  int i;
  for (i = ISTART; i < n; i += IINC) { //for(i=blockId.x;i<numberofbases;i+=gridDim.x) so each thread will handle 0+blockId.x,256+blockId.x.. 
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    int_t wij = INF; //temp variable to fold free energy sum
    int_t* v_temp;
    //consider adding nucleotide to existing loop
    if(d>0){
      if (i!=n-1)
        free_energy_min(&wij, w[ind(i+1,j,n)] + p->b);
      if(j!=0)
        free_energy_min(&wij, w[ind(i,j-1,n)] + p->b);
    }
/*//this code is correct and more efficient but does not perfectly match fold, which leaves some things out.
    if((v_temp = array_val(v,i,j,n,s)))
      free_energy_min(&wij, *v_temp + terminal_U_penalty(s,i,j,p) + p->c);
    if((j!=0) && ((v_temp = array_val(v,i,j-1,n,s))))
//    if(v_temp = array_val(v,i,j-1,n,s))
      free_energy_min(&wij, *v_temp + dangle_3p_energy(s,j-1,i,j,p) + p->b + p->c);

//    if(v_temp = array_val(v,i+1,j,n,s)){
    if((i!=n-1) && ((v_temp = array_val(v,i+1,j,n,s))))
      free_energy_min(&wij, *v_temp + dangle_5p_energy(s,j,i+1,i,p) + p->b + p->c);

    if((i!=n-1) && (j!=0) && ((v_temp = array_val(v,i+1,j-1,n,s))))
//    if(v_temp = array_val(v,i+1,j-1,n,s))
      free_energy_min(&wij, *v_temp + terminal_stack_multibranch(s,j-1,i+1,j,i,p) + 2*p->b + p->c);
//      free_energy_min(&wij, *v_temp + terminal_stack(s,i,j,i-1,j+1,p) + 2*p->b + p->c);
*/
    if((is_interior(i,j) && (d>LOOP_MIN-1)) ){
      v_temp = array_val(v,i,j,n,s);
        free_energy_min(&wij, (v_temp? *v_temp:INF) + terminal_U_penalty(s,i,j,p) + p->c);
      if(j!=0){
        v_temp = array_val(v,i,j-1,n,s);
        free_energy_min(&wij, (v_temp? *v_temp:INF) + dangle_3p_energy(s,j-1,i,j,p) + p->b + p->c);
      }
   
      if(i!=n-1) {
        v_temp = array_val(v,i+1,j,n,s);
        free_energy_min(&wij, (v_temp? *v_temp:INF) + dangle_5p_energy(s,j,i+1,i,p) + p->b + p->c);
        }
   
      if((i!=n-1) && (j!=0)){
        v_temp = array_val(v,i+1,j-1,n,s);
        free_energy_min(&wij, (v_temp? *v_temp:INF) + terminal_stack_multibranch(s,j-1,i+1,j,i,p) + 2*p->b + p->c);
      }
    }
    if(is_exterior(i,j)){
      free_energy_min(&wij, v[ind(i,j,n)] + terminal_U_penalty(s,i,j,p) + p->c);
      if(j!=0){
        //v_temp = array_val(v,i,j-1,n,s);
        free_energy_min(&wij, v[ind(i,j-1,n)] + dangle_3p_energy(s,j-1,i,j,p) + p->b + p->c);
      }
   
      if(i!=n-1) {
        //v_temp = array_val(v,i+1,j,n,s);
        free_energy_min(&wij, v[ind(i+1,j,n)] + dangle_5p_energy(s,j,i+1,i,p) + p->b + p->c);
        }
   
      if((i!=n-1) && (j!=0)){
        //v_temp = array_val(v,i+1,j-1,n,s);
        free_energy_min(&wij, v[ind(i+1,j-1,n)] + terminal_stack_multibranch(s,j-1,i+1,j,i,p) + 2*p->b + p->c);
      }
    }/*
    if (is_exterior(i,j)){
      if((v_temp = array_val(v,i,j,n,s)))
        free_energy_min(&wij, *v_temp + terminal_U_penalty(s,i,j,p) + p->c);
      if((j!=0) && ((v_temp = array_val(v,i,j-1,n,s))))
        free_energy_min(&wij, *v_temp + dangle_3p_energy(s,j-1,i,j,p) + p->b + p->c);
   
      if((i!=n-1) && ((v_temp = array_val(v,i+1,j,n,s))))
        free_energy_min(&wij, *v_temp + dangle_5p_energy(s,j,i+1,i,p) + p->b + p->c);
   
      if((i!=n-1) && (j!=0) && ((v_temp = array_val(v,i+1,j-1,n,s))))
        free_energy_min(&wij, *v_temp + terminal_stack_multibranch(s,j-1,i+1,j,i,p) + 2*p->b + p->c);
    }*/

    w[ind(i,j,n)] = wij;
  }
}


GLOBAL static void calc_WM (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{
	//WM(i,j) = min[W(i,k)+W(k+1,j), 
//			V(i,k)+V(k+1,j)+2c+Eflushcoax, 
//			V(i,k)+V(k+2,j-1)+2c+Ecoax5'mismatch, 
//			V(i+1,k)+V(k+2,j)+2c+Ecoax3'mismatch]

  int i;
  for (i = ISTART; i < n; i += IINC) { //for(i=blockId.x;i<numberofbases;i+=gridDim.x) so each thread will handle 0+blockId.x,256+blockId.x.. 
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    int_t tmp = INF;

//don't need to calculate every WM
  if((is_interior(i,j) && (j-i-1 <= 2*LOOP_MIN+2))){//condition copied verbatim from algorithm.cpp
    wm[ind(i,j,n)]=INF;
    continue;
  }
//  if(is_interior(i,j)) printf("i %d j %d\n",i,j);


#ifdef __CUDACC__
    const int kstart = i + threadIdx.x;
//    const int kstart = i+1 + threadIdx.x;
    const int kinc = blockDim.x;
#else
    const int kstart = i;
//    const int kstart = i+1;
    const int kinc = 1;
#endif
    int ktmp;
//    for (ktmp = kstart; ktmp < jtmp-1; ktmp += kinc) {
    for (ktmp = kstart; ktmp < jtmp; ktmp += kinc) {
      if (ktmp != n-1) {
        const int k = wrap(ktmp,n);
        free_energy_min(&tmp, w[ind(i,k,n)] + w[ind(k+1,j,n)]);
#ifndef disablecoax
	//placeholder for coaxial stacking recursions
#endif//disablecoax
      }

    }

    if(d>0){
      if (i!=n-1)
        free_energy_min(&tmp, wm[ind(i+1,j,n)] + p->b);
      if(j!=0)
        free_energy_min(&tmp, wm[ind(i,j-1,n)] + p->b);
    }

#ifdef __CUDACC__
    free_energy_min_reduce(&tmp, threadIdx.x, blockDim.x);
    if (threadIdx.x != 0)
      continue;
#endif
    wm[ind(i,j,n)] = tmp;
    free_energy_min(&w[ind(i,j,n)],tmp);
  }
}

#ifndef disablecoax
GLOBAL static void calc_coaxial (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{
  int i;
  for (i = ISTART; i < n; i += IINC) {
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    if ((is_exterior(i,j) && i-j <= LOOP_MIN) || !can_pair(i,j,n,s))
      continue;
    const int_t *v1;
    int_t vij = INF;
    /* exterior */
    if (is_exterior(i,j)) {
      int k, kstart;
#ifdef __CUDACC__
      kstart = threadIdx.x;
      const int kinc = blockDim.x;
#else
      kstart = 0;
      const int kinc = 1;
#endif
      for (k = kstart; k < j - LOOP_MIN; k += kinc) {
	if ((v1 = array_val(v,k,j-1,n,s)))
	  free_energy_min(&vij, w3[i+1] + w5[k-1] + coaxial_flush(s,k,j-1,j,i,p) + (*v1));
	if (j-2 >= 0) {
	  if (i < n-1 && (v1 = array_val(v,k,j-2,n,s)))
	    free_energy_min(&vij, w3[i+2] + w5[k-1] + coaxial_mismatch2(s,k,j-2,j,i,p) + (*v1));
	  if ((v1 = array_val(v,k+1,j-2,n,s)))
	    free_energy_min(&vij, w3[i+1] + w5[k-1] + coaxial_mismatch1(s,k+1,j-2,j,i,p) + (*v1));
	}
      }
#ifdef __CUDACC__
      kstart = i+LOOP_MIN+1 + threadIdx.x;
#else
      kstart = i+LOOP_MIN+1;
#endif
      for (k = kstart; k < n; k += kinc) {
	if ((v1 = array_val(v,i+1,k,n,s)))
	  free_energy_min(&vij, w3[k+1] + w5[j-1] + coaxial_flush(s,j,i,i+1,k,p) + (*v1));
	if (j > 0 && (v1 = array_val(v,i+2,k,n,s)))
	  free_energy_min(&vij, w3[k+1] + w5[j-2] + coaxial_mismatch1(s,j,i,i+2,k,p) + (*v1));
	if ((v1 = array_val(v,i+2,k-1,n,s)))
	  free_energy_min(&vij, w3[k+1] + w5[j-1] + coaxial_mismatch2(s,j,i,i+2,k-1,p) + (*v1));
      }
    } /* end exterior */
    
    /* multibranch */
    if (d > 2*LOOP_MIN + 3 && i != n-1 && j != 0) { 
      int ktmp;
#ifdef __CUDACC__
      int ktmpstart = i+2 + threadIdx.x;
      const int ktmpinc = blockDim.x;
#else
      int ktmpstart = i+2;
      const int ktmpinc = 1;
#endif
      for (ktmp = ktmpstart; ktmp < jtmp-2; ktmp += ktmpinc) {
	const int k = wrap(ktmp,n);
	if (k != n-1) {
	  if ((v1 = array_val(v,i+1,k,n,s)))
	    free_energy_min(&vij, coaxial_flush(s,j,i,i+1,k,p) + (*v1) + p->a_2c + 
				   w[ind(k+1,j-1,n)]);
	  if (ktmp+2 < jtmp-1 && i+1 != n-1 && k+1 != n-1 && (v1 = array_val(v,i+2,k,n,s))) {
	    const int_t tmp = (*v1) + p->a_2b_2c;
	    free_energy_min(&vij, coaxial_mismatch2(s,j,i,i+2,k,p) + tmp + w[ind(k+2,j-1,n)]);
	    if (j != 1) {
	      free_energy_min(&vij, coaxial_mismatch1(s,j,i,i+2,k,p) + tmp + w[ind(k+1,j-2,n)]);          
	    }
	  }
	}
      }
#ifdef __CUDACC__
      ktmpstart = i+3 + threadIdx.x;
#else
      ktmpstart = i+3;
#endif
      for (ktmp = ktmpstart; ktmp < jtmp-1; ktmp += ktmpinc) {
	const int k = wrap(ktmp,n);
	if (k != 0) {
	  if ((v1 = array_val(v,k,j-1,n,s)))
	    free_energy_min(&vij, coaxial_flush(s,k,j-1,j,i,p) + (*v1) + p->a_2c + 
				   w[ind(i+1,k-1,n)]);
	  if (j != 1 && ktmp > i+3 && (v1 = array_val(v,k,j-2,n,s))) {
	    const int_t tmp = (*v1) + p->a_2b_2c;
	    if (k != 1)
	      free_energy_min(&vij, coaxial_mismatch1(s,k,j-2,j,i,p) + tmp + w[ind(i+1,k-2,n)]);
	    if (i != n-2)
	      free_energy_min(&vij, coaxial_mismatch2(s,k,j-2,j,i,p) + tmp + w[ind(i+2,k-1,n)]);
	  }
	}
      }
    } /* end multibranch */
#ifdef __CUDACC__
    free_energy_min_reduce(&vij, threadIdx.x, blockDim.x);
    if (threadIdx.x != 0)
      continue;
#endif
    free_energy_min(&v[ind(i,j,n)], vij);
  } /* end loop over i */
} /* end calc_coaxial */

GLOBAL static void calc_wl_coax(int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{
  int i;
  for (i = ISTART; i < n; i += IINC) {
    const int jtmp = i+d+1;
    const int j = wrap(jtmp,n);
    if ((is_exterior(i,j) && i-j <= LOOP_MIN) ||
	(is_interior(i,j) && d <= 2*LOOP_MIN+1))
      continue;
#ifdef __CUDACC__
    const int kstart = i+LOOP_MIN+1 + threadIdx.x;
    const int kinc = blockDim.x;
#else
    const int kstart = i+LOOP_MIN+1;
    const int kinc = 1;
#endif
    int ktmp;
    int_t tmp1 = INF, tmp2 = INF;
    for (ktmp = kstart; ktmp < jtmp-LOOP_MIN-1; ktmp += kinc) {
      const int k = wrap(ktmp,n);
      if (k == n-1)
	continue;
      int_t *v1, *v2;
      if ((v1 = array_val(v,i,k,n,s)) && (v2 = array_val(v,k+1,j,n,s)))
	free_energy_min(&tmp1, (*v1) + (*v2) + coaxial_flush(s,i,k,k+1,j,p));
      if (j == 0 || k+1 == n-1)
	continue;
      if (i != n-1 && (v1 = array_val(v,i+1,k,n,s)) && (v2 = array_val(v,k+2,j,n,s)))
	free_energy_min(&tmp2, (*v1) + (*v2) + coaxial_mismatch1(s,i+1,k,k+2,j,p));
      if ((v1 = array_val(v,i,k,n,s)) && (v2 = array_val(v,k+2,j-1,n,s)))
	free_energy_min(&tmp2, (*v1) + (*v2) + coaxial_mismatch2(s,i,k,k+2,j-1,p));
    }
#ifdef __CUDACC__
    free_energy_min_reduce(&tmp1, threadIdx.x, blockDim.x);
    free_energy_min_reduce(&tmp2, threadIdx.x, blockDim.x);
    if (threadIdx.x != 0)
      continue;
#endif
//    if (is_interior(i,j))
//      free_energy_accumulate(&wq[upper_triangle_index(i,j)], free_energy_sum(tmp1,tmp2));
//    tmp1+=2*p->c;
//    tmp2+=(2*p->b+2*p->c);
//    const int_t wcoax = (tmp1 + 2*p->c, tmp2 + 2*p->b + 2*p->c);
//    const int_t wcoax=tmp1>tmp2? tmp2:tmp1;
    free_energy_min(&wm[ind(i,j,n)], tmp1+2*p->c);
    free_energy_min(&wm[ind(i,j,n)], tmp2+2*p->b+2*p->c);
    free_energy_min(&w[ind(i,j,n)], wm[ind(i,j,n)]);
  } /* end loop over i */
} /* end calc_z */

#endif /* disablecoax */


GLOBAL static void calc_w5_and_w3 (int d, int n, base_t *s, int_t *v, int_t *w, int_t *wm, int_t *w5, int_t *w3, param_t p)
{
#ifdef __CUDACC__
  const int istart = threadIdx.x;
  const int iinc = blockDim.x;
#else
  const int istart = 0;
  const int iinc = 1;
#endif
//  int_t w5tmp = INF, w3tmp = INF;//this should probably be initialized to zero
  int_t w5tmp=0,w3tmp = 0;
  int i;
  int_t* v_temp;
  for (i = istart; i + LOOP_MIN <= d; i += iinc) {

    if((v_temp = array_val(v,i,d+1,n,s)))
//      free_energy_min(&w5tmp, w5[i-1] + *v_temp + terminal_U_penalty(s,i,d+1,p)); 
      free_energy_min(&w5tmp, w5[i-1] + *v_temp + terminal_U_penalty(s,d+1,i,p)); //the nucleotide thats more 3' has to go first in terminal_U_penalty call
    if(d-i>LOOP_MIN){//necessary, or we seg fault because we try to have a pair in a 4mer
//      if((d!=n-2) && (v_temp = array_val(v,i,d,n,s))) //d!=n-2 condition actually keeps you from considering a 3' dangle on the end of the sequence
      if((v_temp = array_val(v,i,d,n,s)))
        //free_energy_min(&w5tmp, w5[i-1] + *v_temp + dangle_3p_energy(s,i,d+1,d,p));
        free_energy_min(&w5tmp, w5[i-1] + *v_temp + dangle_3p_energy(s,d,i,d+1,p));
      if((v_temp = array_val(v,i+1,d+1,n,s)))
        //free_energy_min(&w5tmp, w5[i-1] + *v_temp + dangle_5p_energy(s,i+1,d+1,i,p)); 
        free_energy_min(&w5tmp, w5[i-1] + *v_temp + dangle_5p_energy(s,d+1,i+1,i,p)); 
    }
//    if ((d-i>LOOP_MIN+1) && (d!=n-2) && (v_temp = array_val(v,i+1,d,n,s)))//playing with conditions here
    if ((d-i>LOOP_MIN+1) && ((v_temp = array_val(v,i+1,d,n,s))))
      free_energy_min(&w5tmp, w5[i-1] + *v_temp + terminal_stack(s,d,i+1,d+1,i,p));
    //  free_energy_min(&w5tmp, w5[i-1] + *v_temp + terminal_stack(s,i,d+1,i+1,d,p));//wrong arg order

/*
//now calculate w3
    if(v_temp = array_val(v,n-d-2,n-i-1,n,s))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_U_penalty(s,n-d-2,n-i-1,p));
    if(v_temp = array_val(v,n-d-2,n-i-2,n,s))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + dangle_3p_energy(s,n-d-2,n-i-2,n-i-1,p));
    if((n-d-1 != 0) && (v_temp = array_val(v,n-d-1,n-i-1,n,s)))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + dangle_5p_energy(s,n-d-1,n-i-1,n-d-2,p));
    if((n-i-2 != n-1) && (n-d-1 != 0) && (v_temp = array_val(v,n-d-1,n-i-2,n,s)))
//      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_stack(s,n-d-1,n-i-2,n-d-2,n-i-1,p));
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_stack(s,n-d-1,n-i-2,n-d-2,n-i-1,p));
*/


    if((v_temp = array_val(v,n-d-2,n-i-1,n,s)))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_U_penalty(s,n-i-1,n-d-2,p)); 
    if((v_temp = array_val(v,n-d-2,n-i-2,n,s)))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + dangle_3p_energy(s,n-i-2,n-d-2,n-i-1,p));
    if((n-d-1 != 0) && ((v_temp = array_val(v,n-d-1,n-i-1,n,s))))
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + dangle_5p_energy(s,n-i-1,n-d-1,n-d-2,p));
    if((n-i-2 != n-1) && (n-d-1 != 0) && ((v_temp = array_val(v,n-d-1,n-i-2,n,s))))
//      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_stack(s,n-i-2,n-d-1,n-d-2,n-i-1,p));
//      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_stack(s,n-i-1,n-d-2,n-i-2,n-d-1,p));//changed order of args..
      free_energy_min(&w3tmp, w3[n-i] + *v_temp + terminal_stack(s,n-i-2,n-d-1,n-i-1,n-d-2,p));//changed order of args again..


  }//moved this bracket up
#ifdef __CUDACC__
  free_energy_min_reduce(&w5tmp, threadIdx.x, blockDim.x);
  free_energy_min_reduce(&w3tmp, threadIdx.x, blockDim.x);
  if (threadIdx.x != 0)
    return;
#endif
//}
  w5[d+1] = w5[d]; 
  w3[n-d-2] = w3[n-d-1];
  free_energy_min(&w5[d+1], w5tmp);
  free_energy_min(&w3[n-d-2], w3tmp);
} /* end calc_w5_and_w3 */


GLOBAL static void init_w5_and_w3 (int n,int_t *w5, int_t *w3)
{
  w5[-1]=0;
  w5[0]=w3[n-1]=w3[n]=0;
}

GLOBAL static void init_w(int n,int_t *w)
{
  int i;
  for(i=0;i<n;i++)
    w[ind(i,i,n)]=INF;
}

//MFE recursions end


prna_t prna_new(const char *s, param_t par)
{
  prna_t p = (prna_t) safe_malloc(sizeof(struct prna));
  memset(p, 0, sizeof(struct prna));

  const int n = p->n = strlen(s);
  p->seq = sequence_from_string(s);
  p->v = (int_t *) safe_malloc(n*n*sizeof(int_t));
  p->w = (int_t *) safe_malloc(n*n*sizeof(int_t));
  p->wm = (int_t *) safe_malloc(n*n*sizeof(int_t));
  p->w5 = (int_t *) safe_malloc((n+1)*sizeof(int_t)) + 1;
  p->w5[-1]=0;
  p->w3 = (int_t *) safe_malloc((n+1)*sizeof(int_t));
   

#ifdef __CUDACC__ /* do multithreaded fill on GPU */

  int_t *v,*w,*wm,*w5,*w3;
  cudaEvent_t start, stop, tstart, tstop;
  float t;
  CU(cudaEventCreate(&start));
  CU(cudaEventCreate(&stop));
  CU(cudaEventCreate(&tstart));
  CU(cudaEventCreate(&tstop));

  CU(cudaEventRecord(start,0));

//  int_t *v, *w5, *w3;
  
#define ALLOC(a,sz) CU(cudaMalloc(&a,(sz)*sizeof(int_t)))
  
  ALLOC(v,n*n); //best energy of structure closed by pair i,j. j>i: exterior fragment
  ALLOC(w,n*n); //best energy of structure from i to j
  ALLOC(wm,n*n); //best energy of structure i to j containing 2 or more branches
  ALLOC(w5,n+1); //best energy of structure from 1 to i
  w5++;//w5 is indexed from 1 -- is this a good idea?
  ALLOC(w3,n+1); //best energy of structure from i to numberofbases

/*  
  ALLOC(z,n*n);
  ALLOC(yl,n*n);
  ALLOC(y,n*n);
  ALLOC(wq,n*(n-1)/2);
  ALLOC(w,2*n);
  ALLOC(wl,2*n);
  ALLOC(xl,2*n);
  ALLOC(x,5*n);*/

  param_t dev_par;
  CU(cudaMalloc(&dev_par, sizeof(struct param)));
  CU(cudaMemcpy(dev_par, par, sizeof(struct param), cudaMemcpyHostToDevice));
  
  base_t *dev_s;
  CU(cudaMalloc(&dev_s,n*sizeof(base_t)));
  CU(cudaMemcpy(dev_s, p->seq, n*sizeof(base_t), cudaMemcpyHostToDevice));

  CU(cudaEventRecord(stop,0));
  CU(cudaEventSynchronize(stop));
  CU(cudaEventElapsedTime(&t,start,stop));
  fprintf(stderr, "Time for copy from CPU to GPU: %3.1f ms\n", t);

  int d;
  float t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0, tcoax=0, tcoax2=0;

  init_w5_and_w3<<<1,1>>>(n,w5,w3);
  init_w<<<1,1>>>(n,w);

  CU(cudaEventRecord(tstart,0));
//  for (d = 0; d < n; d++) { //for fragment lengths (1 : n) n-1 is the original way
  for (d = 0; d < n-1; d++) { //for fragment lengths (1 : n)

    #define ARGS d,n,dev_s,v,w,wm,w5,w3,dev_par
    #define TIMER_START CU(cudaEventRecord(start,0));
    #define TIMER_STOP     CU(cudaEventRecord(stop,0));\
                           CU(cudaEventSynchronize(stop));\
                           CU(cudaEventElapsedTime(&t,start,stop));
    TIMER_START
    calc_V_hairpin_and_V_stack<<<n,1>>>(ARGS);
    TIMER_STOP
    t1 += t;

    TIMER_START
    calc_V_bulge_internal<<<n,dim3(SQRT_NTHREAD,SQRT_NTHREAD,1)>>>(ARGS);
    TIMER_STOP
    t2 += t;

    TIMER_START
    calc_V_exterior<<<n,1>>>(ARGS);
  //  cudaError_t err = cudaDeviceSynchronize();
  //  printf(cudaGetErrorString(err);
    TIMER_STOP
    t3 += t;

    TIMER_START
    calc_V_multibranch<<<n,1>>>(ARGS);
    TIMER_STOP
    t4 += t;

#ifndef disablecoax 
    TIMER_START
    calc_coaxial<<<n,NTHREAD>>>(ARGS);
    TIMER_STOP
    tcoax +=t;
#endif

    TIMER_START
    calc_W<<<n,1>>>(ARGS);
//    cudaError_t err = cudaDeviceSynchronize();
  //  printf("%s",cudaGetErrorString(err));
    

    TIMER_STOP
    t5 += t;

    TIMER_START
    calc_WM<<<n,NTHREAD>>>(ARGS);
    TIMER_STOP
    t6 += t;
#ifndef disablecoax
    TIMER_START
    calc_wl_coax<<<n,NTHREAD>>>(ARGS);
    TIMER_STOP
    tcoax2 +=t;
#endif

    TIMER_START
    calc_w5_and_w3<<<1,NTHREAD>>>(ARGS);
    TIMER_STOP
    t7 += t;

  }
  CU(cudaEventRecord(tstop,0));
  CU(cudaEventSynchronize(tstop));
  CU(cudaEventElapsedTime(&t,tstart,tstop));
  
  fprintf(stderr, "Time for hairpin, stack: %3.1f ms\n", t1);
  fprintf(stderr, "Time for exterior: %3.1f ms\n", t2);
  fprintf(stderr, "Time for bulge/internal loop: %3.1f ms\n", t3);
  fprintf(stderr, "Time for multibranch: %3.1f ms\n", t4);
  fprintf(stderr, "Time for V coax: %3.1f ms\n", tcoax);
  fprintf(stderr, "Time for WM coax: %3.1f ms\n", tcoax2);
  fprintf(stderr, "Time for W: %3.1f ms\n", t5);
  fprintf(stderr, "Time for WM: %3.1f ms\n", t6);
  fprintf(stderr, "Time for w5 and w3: %3.1f ms\n", t7);
  fprintf(stderr, "Time for fill (total): %3.1f ms\n", t);
  
  CU(cudaEventRecord(start,0));
  
  CU(cudaMemcpy(p->v, v, n*n*sizeof(int_t), cudaMemcpyDeviceToHost));  
  CU(cudaMemcpy(p->w, w, n*n*sizeof(int_t), cudaMemcpyDeviceToHost));  
  CU(cudaMemcpy(p->wm, wm, n*n*sizeof(int_t), cudaMemcpyDeviceToHost));  
  CU(cudaMemcpy(p->w5 - 1, w5 - 1, (n+1)*sizeof(int_t), cudaMemcpyDeviceToHost));  
  CU(cudaMemcpy(p->w3, w3, (n+1)*sizeof(int_t), cudaMemcpyDeviceToHost));  
  
  CU(cudaEventRecord(stop,0));
  CU(cudaEventSynchronize(stop));
  CU(cudaEventElapsedTime(&t,start,stop));

  fprintf(stderr, "Time for copy from GPU to CPU: %3.1f ms\n", t);

  CU(cudaFree(v));
  CU(cudaFree(w5 - 1));
  CU(cudaFree(w3));
  CU(cudaFree(w));
  CU(cudaFree(wm));
  CU(cudaEventDestroy(start));
  CU(cudaEventDestroy(stop));
  CU(cudaEventDestroy(tstart));
  CU(cudaEventDestroy(tstop));

#else /* do serial fill on CPU */

#define ALLOC(a,sz) a = (int_t *) safe_malloc((sz)*sizeof(int_t))
#define ARGS d,p->n,p->seq,p->v,p->w,p->wm,p->w5,p->w3,par

/*  ALLOC(v,n*n); //best energy of structure closed by pair i,j. j>i: exterior fragment
  ALLOC(w,n*n); //best energy of structure from i to j
  ALLOC(wm,n*n); //best energy of structure i to j containing 2 or more branches
  ALLOC(w5,n+1); //best energy of structure from 1 to i
  w5++;//w5 is indexed from 1 -- is this a good idea?
  ALLOC(w3,n+1); //best energy of structure from i to numberofbases
*/


  init_w5_and_w3(n,p->w5,p->w3);
  init_w(n,p->w);
 
  int d;

  for (d = 0; d < n-1; d++) {

    calc_V_hairpin_and_V_stack(ARGS);

    calc_V_bulge_internal(ARGS);

    calc_V_exterior(ARGS);

    calc_V_multibranch(ARGS);
#ifndef disablecoax
    calc_coaxial(ARGS);
#endif
    calc_W(ARGS);

    calc_WM(ARGS);
#ifndef disablecoax
    calc_wl_coax(ARGS);
#endif
    calc_w5_and_w3(ARGS);

  }
printf("done with array fill\n");
  /* 
  free(v);    
  free(w);
  free(wm);
  free(w5-1);
  free(w3);
   */
#endif /* __CUDACC__ */
  /*
p->v[ind(5,3,p->n)] = -88;
p->v[ind(5,4,p->n)] = -54;
p->v[ind(2,1,p->n)] = -12;
p->v[ind(3,1,p->n)] = -2;

p->w[ind(5,3,p->n)] = -89;
p->w[ind(5,4,p->n)] = -90;
p->w[ind(4,3,p->n)] = -91;

p->wm[ind(3,0,p->n)] = 13986;
p->wm[ind(4,0,p->n)] = 13988;
p->wm[ind(4,2,p->n)] = 13925;
p->wm[ind(4,2,p->n)] = 13925;
p->wm[ind(4,3,p->n)] = 13911;
p->wm[ind(5,0,p->n)] = 13988;
p->wm[ind(5,1,p->n)] = 13959;
p->wm[ind(5,2,p->n)] = 13927;
p->wm[ind(5,3,p->n)] = 13915;
p->wm[ind(5,4,p->n)] = 13911;

 
p->wm[ind(6,0,p->n)] = 13988;
p->wm[ind(6,1,p->n)] = 13961;
p->wm[ind(6,2,p->n)] = 13957;
p->wm[ind(6,3,p->n)] = 13927;
p->wm[ind(6,4,p->n)] = 13915;
p->wm[ind(6,5,p->n)] = 13915;


p->wm[ind(7,0,p->n)] = 13988;
p->wm[ind(7,1,p->n)] = 13981;
p->wm[ind(7,2,p->n)] = 13961;
p->wm[ind(7,3,p->n)] = 13957;

////p->wm[ind(8,3,p->n)] = 13957;
p->wm[ind(8,1,p->n)] = 13988;
p->wm[ind(8,0,p->n)] = 13988;
//p->wm[ind(8,4,p->n)] = 13981;
p->wm[ind(8,2,p->n)] = 13981;
*/
  return p;
} /* end prna_new */

void prna_delete(prna_t p)
{
  if (p) {
    if (p->seq)
      free(p->seq);
    if (p->v)
      free(p->v);
    if (p->w)
      free(p->w);
    if (p->wm)
      free(p->wm);
    if (p->w5 - 1)
      free(p->w5 - 1);
    if (p->w3)
      free(p->w3);
    free(p);
  }
}

#define SHOWARR(a)	 \
  if (p->a) {		      \
    int i, j;	      \
    for (i = 0; i < n; i++) { \
      printf("%s%4d: ",#a,i+1);				\
      for (j = 0; j < n; j++) {				\
	const int_t *aij = array_val(p->a,i,j,n,s);		\
	printf(RF" ", aij ? (*aij) : INF);	\
      }								\
      printf("\n");						\
    }								\
  }

#define SHOW(a)							\
  if (p->a) {								\
    int i;								\
    printf("%s: ",#a);							\
    for (i = 0; i < n; i++)						\
      printf(RF" ", p->a[i]);				\
    printf("\n");							\
  }									\
  //used to be mul by RT

/*void SHOW_2DARRAY(int_t* arr, int n)
{
  int i,j;
  printf("i\tj\tvalue\n");
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      printf("%d\t%d\t%d\n",i,j,arr[ind(i,j,n)]);
}
*/

void prna_show(const prna_t p)
{
  int i,j, n = p->n;
  const base_t *s = p->seq;
  printf("n: %d\n", n);
  printf("seq: ");
  for (i = 0; i < n; i++)
    printf("%c", base_as_char(s[i]));
  printf("\n");
  printf("i\tj\tV:\tW:\tWM:\tV':\tW':\tWM':\n");
  for(j=0;j<n;j++)
    for(i=0;i<j;i++)
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i+1,j+1,p->v[ind(i,j,n)],p->w[ind(i,j,n)],p->wm[ind(i,j,n)],p->v[ind(j,i,n)],p->w[ind(j,i,n)],p->wm[ind(j,i,n)] );
        //printf("%d\t%d\t%d\n",i+1,j+1,p->w[ind(j,i,n)]);
  

/*
SHOWARR(w);
  SHOWARR(wm);
*/
//  SHOW(w5);
//  SHOW(w3);
  printf("\n\n\ni\tw5[i]\tw3[i]\n");
  printf("0\t0\t0\n");
  for(i=0;i<n;i++){
    printf("%d\t",i+1);
    printf("%d\t",p->w5[i]);
    printf("%d\n",p->w3[i]); 
  }

}

short base_as_num(base_t b)
{
  switch (b) {
  case A:
    return 1;
  case C:
    return 2;
  case G:
    return 3;
  case U:
    return 4;
  default:
    printf("unknown base %d\n",b);
    die("base_as_num: unknown base");
    return 0;
  }

}

base_t num_as_base(short x)
{
  switch (x) {
  case 1:
    return A;
  case 2:
    return C;
  case 3:
    return G;
  case 4:
    return U;
  default:
//    die("base_as_num: unknown base");
    return A;
  }
}
#ifdef SHORT
int int_pow(int val,int exp)//val^exp, integer type
{
  if (exp<0) die("int_t_pow:tried to take negative exponent");
  if (exp==0) return 1;
  else if (exp==1) return val;
  return val * int_pow(val,exp-1);
}
#endif
int_t int_t_pow(int_t val,int_t exp)//val^exp, integer type
{
  if (exp<0) die("int_t_pow:tried to take negative exponent");
  if (exp==0) return 1;
  else if (exp==1) return val;
  return val * int_t_pow(val,exp-1);
}

int_t int_t_min(int_t a, int_t b)
{
  return a>b? b:a;
}

#define write_int_t(x) fwrite(x,sizeof(int_t),1,savefile)
#define write_int(x) fwrite(x,sizeof(int),1,savefile)
#define write_short(x) fwrite(x,sizeof(short),1,savefile)
#define write_char(x) fwrite(x,sizeof(char),1,savefile)
#define write_float(x) fwrite(x,sizeof(float),1,savefile)

#define b(x) num_as_base(x)

void prna_write_save_file(const prna_t p, const param_t par)
{
  int i,j,k,l,m,n,o,q;
  FILE *savefile;
  int_t zero=0;
  int intzero = 0;
  int_t infinity=INF;
  char FALSE=0;
  int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
                                                {0,1,0,1,0,0},{0,0,0,0,0,0}};//"can-pair" array
  savefile = fopen("cuda.out","wb");
  if (!savefile) die("failed to open output file");
//save file version
  short vers=3;//this is safiversion from defines.h
  fwrite(&vers,sizeof(short),1,savefile);

//basic structure information
  int sequencelength = p->n;
  char intermolecular=0; //fold-cuda does not currentlty suport intermolecular folding
  fwrite(&sequencelength,sizeof(int),1,savefile);
  fwrite(&intermolecular,sizeof(char),1,savefile);

  int pairs = 0; //we do not currently output any pairs, just the energies for dotplot creation
  fwrite(&pairs,sizeof(int),1,savefile);
  int forbiddenpairs=0;//likewise
  fwrite(&forbiddenpairs,sizeof(int),1,savefile);
  
//write the nucleotide sequence. hnumber is the base position
  char nuc=0;
  short hnum;
  for(i=0;i<=p->n;i++){
    if(i>0) nuc = base_as_char(p->seq[i-1]);//i-1 because seq is 0-indexed but rnastructure expects 1-indexed
    hnum = i;
    fwrite(&hnum,sizeof(short),1,savefile);//hnumber[0] is nothing, as is nuc[0]
    fwrite(&nuc,sizeof(char),1,savefile);
  }

  short numseq[2*p->n+1];//numseq array contains the sequnece in a different format
  for(i=1;i<=p->n;i++){
    numseq[i] = numseq[i+p->n] = base_as_num(p->seq[i-1]);
  }
  for(i=0;i<=2*p->n;i++) 
    fwrite(&numseq[i],sizeof(short),1,savefile);

  int doubles=0;//forced double stranded nucs, not supported here
  fwrite(&doubles, sizeof(int),1,savefile);

  int singles=0;//forced single stranded nucleotides, we don't have those either
  fwrite(&singles,sizeof(int),1,savefile);

  int modified=0;//or modifed based
  fwrite(&modified,sizeof(int),1,savefile);

  int numberGU=0;//or nucs constrained to be in GU pairs
  fwrite(&numberGU,sizeof(int),1,savefile);

  const char* label="sequence";
  int length = (int) strlen(label);
  //write_int(&length);
  write_int(&intzero);
  //fwrite(label,sizeof(char),length,savefile);

  char templated=0;
  fwrite(&templated,sizeof(char),1,savefile);//we're not constrained according to a template
  
  char shaped=0;
  fwrite(&shaped,sizeof(char),1,savefile);//we're not constrained by SHAPE data

//int testnum=55;//for debugging
//fwrite(&testnum,sizeof(int),1,savefile);
  int jtmp;
//now, write the arrays
  int_t w5=0,w3=0,v=0,w=0,wm=0;
  char fce=0;
  for(i=0;i<=p->n;i++){
    if (i>0) {w5=p->w5[i-1]; w3 = p->w3[i-1];} 
      write_int_t(&w3);
      write_int_t(&w5);
    for(jtmp=0;jtmp<=p->n;jtmp++){
      j= jtmp+i > p->n? jtmp+i-p->n:jtmp+i;
      if((i>0) && (jtmp>0)){
        v = p->v[ind(i-1,j-1,p->n)];
        w = p->w[ind(i-1,j-1,p->n)];
        wm = p->wm[ind(i-1,j-1,p->n)];
      }
      else v=w=wm=14000;
      write_int_t(&v);
      write_int_t(&w);
      write_int_t(&wm);
      write_char(&fce);
      
    }
  }


  write_int_t(&p->w3[p->n]);
//int_t six=6;
//int_t seven=7;
  for(i=0;i<=2*p->n;i++){
    write_char(&FALSE);//fce[i]
    write_char(&FALSE);//mod[i]
  }

//write vmin
int_t vmintmp = INF;
//  int vmin = (int) p->w3[0];
  for(i=0;i<p->n-1;i++){
 //   printf("\ni=%d ",i); 
    for(j=i+1;j<p->n;j++) { 
   //   printf( "j=%d ",j);
      vmintmp = int_t_min(vmintmp, p->v[ind(i,j,p->n)]+p->v[ind(j,i,p->n)]);
    }
  }
int vmin = (int) vmintmp;
//printf("vmin=%d\n",vmin);

  write_int(&vmin);
//now, write a the thermodynamic data in a format that refold wants
  int_t poppen[5] = {0,6,6,6,6};
  for (i=0;i<5;i++) write_int_t(&poppen[i]); //data->poppen[i] //used in partition function

  write_int_t(&par->maximum_correction); //data->maxpen

  int_t eparam[11] = {0,0,0,0,0,par->a,par->b,30,30,-500,par->c};
  for (i=0;i<11;i++) write_int_t(&eparam[i]); //data->eparam[i]
  for (i=0;i<31;i++) {
    write_int_t(&par->internal_loop_initiation[i]);//data->inter[i]
    write_int_t(&par->bulge_loop_initiation[i]);//data->bulge[i]
    write_int_t(&par->hairpin_loop_initiation[i]);//data->hairpin[i]
  }
  for (i=0;i<6;i++) {
    for (j=0;j<6;j++) {
      for (k=0;k<6;k++) {
//        for (l=0;l<3;l++) {    
        write_int_t(&zero);//data->dangle[i][j][k][l] l[0] is meaningless
        if(i&&j&&k) write_int_t(&par->dangle_3p[num_as_base(i)][num_as_base(j)][num_as_base(k)]); 
        else write_int_t(&zero);
        if(i&&j&&k) write_int_t(&par->dangle_5p[num_as_base(i)][num_as_base(j)][num_as_base(k)]); 
        else write_int_t(&zero);
//        }
        for (l=0;l<6;l++) {
          if(i&&j&&k&&l)write_int_t(&par->stack[b(i)][b(j)][b(k)][b(l)]);//data->stack[i][j][k][l]
          else write_int_t(&zero);
         
          if(i&&j&&k&&l)write_int_t(&par->tstackh[b(i)][b(j)][b(k)][b(l)]);//data->tstckh[i][j][k][l]
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstacki[b(i)][b(j)][b(k)][b(l)]);//data->tstcki[i][j][k][l]
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->coaxial[b(i)][b(j)][b(k)][b(l)]);//data->coax[i][j][k][l]
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstackcoax[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->coaxstack[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstack[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstackm[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstacki23[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          if(i&&j&&k&&l)write_int_t(&par->tstacki1n[b(i)][b(j)][b(k)][b(l)]);
          else write_int_t(&zero);

          for (m=0;m<6;m++) {
            for (n=0;n<6;n++) {
              if(i&&j&&k&&l&&m&&n)write_int_t(&par->int11[b(i)][b(j)][b(k)][b(n)][b(m)][b(l)]);//iloop11
              else write_int_t(&infinity);

              for (o=0;o<6;o++) {
                if (inc[i][j]&&inc[n][o]){ 
                  if(i&&j&&k&&l&&m&&n&&o)write_int_t(&par->int21[b(i)][b(j)][b(k)][b(l)][b(m)][b(n)][b(o)]);//iloop21
                  else write_int_t(&infinity);
                }
                for (q=0;q<6;q++) {
                  if (inc[i][k]&&inc[j][l]){
                    if(i&&j&&k&&l&&m&&n&&o&&q&&(i!=5)&&(i!=5)&&(j!=5)&&(k!=5)&&(l!=5)&&(m!=5)&&(n!=5)&&(o!=5)&&(q!=5))
                      write_int_t(&par->int22[b(i)][b(j)][b(k)][b(l)][b(m)][b(n)][b(o)][b(q)]);//iloop222
                    else write_int_t(&infinity);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
//translate tetraloop to rnastructure format
  int_t ntloop = (int_t) par->ntloop;
  int_t tloop[par->ntloop+1][2];
  for(i=0;i<=par->ntloop;i++){
    tloop[i][0]=tloop[i][1]=0;
    for(j=0;j<6;j++){//6nt in a tetraloop
      if (i!=0) tloop[i][0] += (base_as_num(par->tloop[i-1].seq[j]) * int_t_pow(5,j));
    }
    if (i!=0) tloop[i][1] = par->tloop[i-1].val;
    else tloop[i][1]=0;
  }
//write tetraloop
  write_int_t(&ntloop);//data->numoftloops));
  for (i=0;i<=par->ntloop;i++) {
    for (j=0;j<2;j++) write_int_t(&tloop[i][j]);
  }

//translate triloop to rnastructure format
  int_t ntriloop = (int_t) par->ntriloop;//we have them as ints, they need to be shorts
  int_t triloop[par->ntriloop+1][2];
  for(i=0;i<=par->ntriloop;i++){
    triloop[i][0]=triloop[i][1]=0;
    for(j=0;j<5;j++){//5nt in a triloop
      if (i!=0) triloop[i][0] += (base_as_num(par->triloop[i-1].seq[j]) * int_t_pow(5,j));
    }
    if (i!=0) triloop[i][1] = par->triloop[i-1].val;
  }
//write triloop data
  write_int_t(&ntriloop);//data->numoftriloops));
  for (i=0;i<=par->ntriloop;i++) {
    for (j=0;j<2;j++) write_int_t(&triloop[i][j]);
  }

//translate hexaloop to rnastructure format
  int_t nhexaloop = (int_t) par->nhexaloop;//we have them as ints, they need to be shorts
  int hexaloop[par->nhexaloop+1][2];//has to be an int because the numbers to represent sequence get big
  for(i=0;i<=par->nhexaloop;i++){
    hexaloop[i][0]=triloop[i][1]=0;
    for(j=0;j<8;j++){//5nt in a triloop
#ifdef SHORT//if int_t is int, then we use the int_t_pow function instead of int_pow (there's probably a better way to do this)
      if (i!=0) hexaloop[i][0] += (base_as_num(par->hexaloop[i-1].seq[j]) * int_pow(5,j));
#else
      if (i!=0) hexaloop[i][0] += (int) (base_as_num(par->hexaloop[i-1].seq[j]) * int_t_pow(5,j));
#endif
    }
    if (i!=0) hexaloop[i][1] = (int) par->hexaloop[i-1].val;
  }
//write hexaloop data
  write_int_t(&nhexaloop);//data->numofhexaloops));
  for (i=0;i<=par->nhexaloop;i++) {
    for (j=0;j<2;j++) write_int(&hexaloop[i][j]);
  }

  write_int_t(&par->terminal_AU_penalty); //data->auend               
  write_int_t(&par->bonus_for_GGG_hairpin); //data->gubonus               
  write_int_t(&par->c_hairpin_intercept); //data->cint               
  write_int_t(&par->c_hairpin_slope); //data->cslope
  write_int_t(&par->c_hairpin_of_3); //data->c3               
  write_int_t(&par->a);//data->efn2a               
  write_int_t(&par->b);//data->efn2b  
  write_int_t(&par->c);//data->efn2c               
  write_int_t(&zero); //data->init (we skip)              
  write_int_t(&zero); //mlasym (we skip)
  write_int_t(&zero); //strain (we skip)
//  float prelog = (float) par->Extrapolation_for_large_loops;//rnastructure expects this as a float
  write_float(&par->prelog); //data->prelog  
  write_int_t(&par->Bonus_for_Single_C_bulges_adjacent_to_C);//data->singlecbulge               

  fclose(savefile);
}

static int_t free_energy_of_pair(const prna_t p, int i, int j)
{
  const int n = p->n;
  const base_t *s = p->seq;
  if (can_pair(i,j,n,s))
    return *array_val(p->v,i,j,n,s) + *array_val(p->v,j,i,n,s) - p->w3[0];
  else
    return INF;
}

static int_t probability_of_pair(const prna_t p, int i, int j)
{
  return exp(-free_energy_of_pair(p,i,j));
}

void prna_write_neg_log10_probabilities(const prna_t p, const char *fn)
{
  FILE *f = safe_fopen(fn,"w");
  int i, j;
  fprintf(f,"%d\n%-8s%-8s-log10(probability)\n",p->n,"i","j");
  for (i = 0; i < p->n; i++)
    for (j = i+1; j < p->n; j++)
      if (can_pair(i,j,p->n,p->seq))
	fprintf(f,"%-8d%-8d"RF"\n", i+1, j+1,
	        free_energy_of_pair(p,i,j)/LOG(10));
  fclose(f);
}

void prna_write_probability_matrix(const prna_t p, const char *fn)
{
  FILE *f = safe_fopen(fn,"w");
  const int n = p->n;
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      fprintf(f,RF" ", 
	      can_pair(i,j,n,p->seq) ? probability_of_pair(p,i,j) : 0);
    fprintf(f,"\n");
  }
  fclose(f);
}

static void write_ct_structure(FILE *f, const char *s, int n, const int *pair)
{
  char fmt[256];
  sprintf(fmt,"%d",n);
  int ns = strlen(fmt)+1;
  if (ns < 5)
    ns = 5;
  sprintf(fmt,"%%%dd",ns);
  int i;
  for (i = 0; i < n; i++) {
    fprintf(f,fmt,i+1);
    fprintf(f,"%2c   ",s[i]);
    fprintf(f,fmt,i);
    fprintf(f,fmt,i == n-1 ? 0 : i+2);
    fprintf(f,fmt,pair[i] == i ? 0 : pair[i]+1);
    fprintf(f,fmt,i+1);
    fprintf(f,"\n");
  }
}

static void unpair(int *pair, int i)
{
  const int j = pair[i];
  pair[i] = i;
  pair[j] = j;
}

static int is_paired(const int *pair, int i)
{
  return pair[i] != i;
}
       
static void remove_helices_shorter_than(int min_helix_length, int *pair, int n)
{
  int i;
  for (i = 0; i < n-2; i++) {
    int j = pair[i];
    if (j <= i)
      continue;
    int npair = 1;
    while (pair[i+1] == j-1 || pair[i+2] == j-1 || pair[i+1] == j-2) {
      if (pair[i+1] == j-1)
	;
      else if (pair[i+2] == j-1) {
	if (is_paired(pair,i+1))
	  unpair(pair,i+1);
	i++;
      } else
	j--;
      i++;
      j--;
      npair++;
    }
    if (npair < min_helix_length) {
      unpair(pair,i);
      if (i >= 2) {
	while (pair[i-1] == j+1 || pair[i-2] == j+1 || pair[i-1] == j+2) {
	  if (pair[i-1] == j+1)
	    unpair(pair,i-1);
	  else if (pair[i-2] == j+1) {
	    unpair(pair,i-2);
	    i--;
	  } else {
	    unpair(pair,i-1);
	    j++;
	  }
	  i--;
	  j++;
	}
      } else if (i == 1) {
	while (pair[i-1] == j+1 || pair[i-1] == j+2) {
	  if (pair[i-1] == j+1)
	    unpair(pair,i-1);
	  else {
	    unpair(pair,i-1);
	    j++;
	  }
	  i--;
	  j++;
	}
      }
    }
  }
} /* end remove_helices_shorter_than */

void prna_write_probknot(const prna_t p, const char *fn, const char *s, int min_helix_length)
{
  const int n = p->n;
  int *pair = (int *) safe_malloc(n*sizeof(int));
  int i;
  for (i = 0; i < n; i++) {
    pair[i] = i; /* unpaired */
    int j;
    for (j = 0; j < n; j++)
      if (free_energy_of_pair(p,i,j) < free_energy_of_pair(p,i,pair[i]))
	pair[i] = j;
  }
  for (i = 0; i < n; i++)
    if (pair[pair[i]] != i)
      pair[i] = i; /* unpaired */
  if (min_helix_length > 1)
    remove_helices_shorter_than(min_helix_length,pair,n);
  /* write the structure */
  if (fn) {
    FILE *f = safe_fopen(fn,"w");
    write_ct_structure(f,s,n,pair);
    fclose(f);
  } else {
    write_ct_structure(stdout,s,n,pair);
  }
  free(pair);
}
