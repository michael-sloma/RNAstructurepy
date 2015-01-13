#ifndef PRNA_H
#define PRNA_H

#include "param.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct prna *prna_t;
  prna_t prna_new(const char *s, param_t par);
  void prna_delete(prna_t);
  void prna_show(const prna_t);
  void prna_write_save_file(const prna_t, const param_t);
  void prna_write_neg_log10_probabilities(const prna_t, const char *fn);
  void prna_write_probability_matrix(const prna_t, const char *fn);
  void prna_write_probknot(const prna_t, const char *fn, const char *s, int min_helix_length);
  
#ifdef __cplusplus
}
#endif

#endif /* PRNA_H */
