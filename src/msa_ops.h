#ifndef MSA_OPS_H
#define MSA_OPS_H

#ifdef MSA_OPS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa_info;
struct msa;
struct mumsa_data;

EXTERN int read_msa_into_msai(struct msa* msa,struct msa_info** msai);
EXTERN int msai_char_to_pos(struct msa_info* msai);
EXTERN int alistat(struct mumsa_data* m);


EXTERN int sanity_check_input(struct mumsa_data* md);
EXTERN int sort_msa_by_seqname(struct msa* msa);

#undef MSA_OPS_IMPORT
#undef EXTERN
#endif
