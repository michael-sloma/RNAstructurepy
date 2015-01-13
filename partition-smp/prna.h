#ifndef PRNA_H
#define PRNA_H

#include "param.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct prna *prna_t;
  prna_t prna_new(const char *s, param_t par, int quiet=true);
  void prna_delete(prna_t);
  void prna_show(const prna_t);
  void prna_write_neg_log10_probabilities(const prna_t, const char *fn);
  void prna_write_probability_matrix(const prna_t, const char *fn);
  void prna_write_probknot(const prna_t, const char *fn, const char *s, int min_helix_length);
  real_t probability_of_pair(const prna_t p, int i, int j);

#ifdef __cplusplus
}
#endif

#endif /* PRNA_H */
