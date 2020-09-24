#ifndef PROCESS_MSA_H
#define PROCESS_MSA_H

#ifdef PROCESS_MSA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int calc_sim_pairs(struct mumsa_data* m);
EXTERN int calc_overlap(struct mumsa_data* m);

#undef PROCESS_MSA_IMPORT
#undef EXTERN
#endif
