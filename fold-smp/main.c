#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>

#include "prna.h"
#include "util.h"
#include "base.h"

static const char Usage[] = 
  "Usage: %s [options] <sequence, or file containing one>\n"
  "\n"
  "options:\n"
  "-h:        show this message\n"
  "-b <file>: read parameters from <file>, in native binary format\n"
  "-d:        use DNA parameters\n"
  "-t <file>: write probability matrix as text to <file>\n"
  "-l <file>: write -log10 probabilities as text to <file>\n"
  "-p <file>: write ProbKnot structure in ct format to <file>\n"
  "-m <length>: set minimum helix length for ProbKnot\n"
  "             (default: 3 base pairs)\n"
  "-v:        show arrays\n\n"
  "If none of -t, -l, -p, or -v is chosen,\n"
  "writes ProbKnot structure in ct format to stdout\n";

int main(int argc, char **argv)
{
  const char *cmd = *argv;
  int use_dna_params = 0, min_helix_length = 3, verbose = 0;
  const char *neg_log10_filename = 0;
  const char *text_matrix_filename = 0;
  const char *probknot_filename = 0; 
  const char *binary_parameter_filename = 0;
  /* process command-line arguments */
  int c;
  while ((c = getopt(argc, argv, "hb:dt:l:p:m:v")) != EOF)
    if (c == 'h')
      die(Usage,cmd);
    else if (c == 'b')
      binary_parameter_filename = optarg;
    else if (c == 'd')
      use_dna_params = 1;
    else if (c == 't')
      text_matrix_filename = optarg;
    else if (c == 'l')
      neg_log10_filename = optarg;
    else if (c == 'p')
      probknot_filename = optarg;
    else if (c == 'm')
      min_helix_length = atoi(optarg);
    else if (c == 'v')
      verbose = 1;
    else
      die(Usage,cmd);
  argc -= optind;
  argv += optind;
  if (argc == 0)
    die(Usage,cmd);
  /* get sequence */
  char *seq = sequence(*argv);
  /* read parameters */
  struct param par;
  if (binary_parameter_filename) {
    param_read_from_binary(binary_parameter_filename, &par);
    if (par.use_dna_params != use_dna_params)
      die("%s: -d option %s, but '%s' %s DNA parameters", cmd,
	  use_dna_params ? "set" : "not set", 
	  binary_parameter_filename, 
	  par.use_dna_params ? "is from" : "is not from");
  } else {
    const char *path = getenv("DATAPATH");
    if (!path)
      die("%s: need to set environment variable $DATAPATH", cmd);
    param_read_from_text(path, use_dna_params, &par);
  }
// param_show(par); 
  /* calculate partition function */
  prna_t p = prna_new(seq, &par);
  printf("RT is %f\n",RT);
#ifdef VERBOSE
  prna_show(p);
#endif
//  param_show(&par);
//printf("%.1f\n",conversion_factor);
  /* output */
 /* FILE *output;
  output = fopen("gpu_arrays_out.txt","w");
  fprintf(output, "i\tj\tv(i,j)\tw(i,j)\twm(i,j)\tw5(i,j)\tw3(i,j)\t\n");
  for(int j=0;j<=p->n;j++)
    for(int i=0;i<=p->n;i++)
      fprintf(output,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t",i,j,p->v(ind(i,j,p->n)),p->w(ind(i,j,p->n)),wm(ind(i,j,p->n)),p->w5(ind(i,j,p->n)),p->w3(ind(i,j,p->n)));
*/

/*  if (neg_log10_filename)
    prna_write_neg_log10_probabilities(p,neg_log10_filename);
  if (text_matrix_filename)
    prna_write_probability_matrix(p,text_matrix_filename);
  if (probknot_filename)
    prna_write_probknot(p,probknot_filename,seq,min_helix_length);
  if (verbose)
    prna_show(p);

  if (!(neg_log10_filename ||
	text_matrix_filename ||
	probknot_filename || 
	verbose))
    prna_write_probknot(p,0,seq,min_helix_length);
*/
  prna_write_save_file(p,&par);

  /* cleanup */
  prna_delete(p);
  free(seq);
  
  return 0;
  
}
