#ifndef MUMSA_DATA_H
#define MUMSA_DATA_H

#include "msa.h"
#include "msa_info.h"

struct sets{
        double score_4_input_aln;
        double aln_sim;
        double id;
        double pcounts;
};


struct mumsa_data{
        struct msa** msa;
        struct msa_info** msai;
        struct sets* s;
        double** overlap;
        double** sim;

        double* poar;
        double* id;
        double* avg_len;
        double* al;
        int num_aln;
        int num_seq;
};


#ifdef MUMSA_DATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int alloc_mumsa_data(struct mumsa_data** mdat, int num_alignment);
EXTERN int free_mumsa_data(struct mumsa_data* m);


#undef MUMSA_DATA_IMPORT
#undef EXTERN
#endif
