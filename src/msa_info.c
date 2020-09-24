#include "tldevel.h"

#define MSA_INFO_IMPORT
#include "msa_info.h"


int alloc_msa_info(struct msa_info** info, int numseq, int aln_len)
{
        struct msa_info* msai = NULL;

        MMALLOC(msai, sizeof(struct msa_info));
        msai->aln_len = aln_len;
        msai->num_seq = numseq;
        msai->s = NULL;
        msai->adler_val = NULL;
        RUN(galloc(&msai->adler_val,msai->num_seq));

        RUN(galloc(&msai->s, msai->num_seq, msai->aln_len));

        *info = msai;
        return OK;
ERROR:
        free_msa_info(msai);
        return FAIL;
}


void free_msa_info(struct msa_info* msai)
{
        if(msai){
                if(msai->s){
                        gfree(msai->s);
                }
                if(msai->adler_val){
                        gfree(msai->adler_val);
                }
                MFREE(msai);
        }
}
