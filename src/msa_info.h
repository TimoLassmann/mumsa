#ifndef MSA_INFO_H
#define MSA_INFO_H

#ifdef MSA_INFO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct msa_info{
        int** s;
        int* adler_val;
        double pairs;
        int aln_len;
        int num_seq;

};

EXTERN int alloc_msa_info(struct msa_info** info, int numseq, int aln_len);
EXTERN void free_msa_info(struct msa_info* msai);

#undef MSA_INFO_IMPORT
#undef EXTERN
#endif
