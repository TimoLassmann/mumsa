#include "tldevel.h"

#define MUMSA_DATA_IMPORT
#include "mumsa_data.h"


int alloc_mumsa_data(struct mumsa_data** mdat, int num_alignment)
{
        struct mumsa_data* m = NULL;
        int i;


        MMALLOC(m, sizeof(struct mumsa_data));
        m->msai = NULL;
        m->msa = NULL;
        m->id = NULL;
        m->poar = NULL;
        m->al = NULL;
        m->avg_len = NULL;
        m->sim = NULL;
        m->num_aln = num_alignment;
        m->num_seq = 0;


        MMALLOC(m->msa, sizeof(struct msa*) * m->num_aln);
        MMALLOC(m->msai, sizeof(struct msa_info*) * m->num_aln);
        for(i = 0; i < m->num_aln;i++){
                m->msa[i] = NULL;
                m->msai[i] = NULL;
        }
        MMALLOC(m->id, sizeof(double) * m->num_aln);
        MMALLOC(m->poar, sizeof(double) * m->num_aln);
        MMALLOC(m->avg_len , sizeof(double) * m->num_aln);
        MMALLOC(m->al, sizeof(double) * m->num_aln);



        *mdat = m;
        return OK;
ERROR:
        return FAIL;
}


int free_mumsa_data(struct mumsa_data* m)
{
        int i;
        if(m){
                if(m->sim){
                        gfree(m->sim);
                }
                if(m->msa){
                        for(i = 0; i < m->num_aln;i++){
                                free_msa(m->msa[i]);
                        }
                        MFREE(m->msa);
                }
                if(m->msai){
                        for(i = 0; i < m->num_aln;i++){
                                free_msa_info(m->msai[i]);
                        }
                        MFREE(m->msai);
                }

                if(m->id){
                        MFREE(m->id);
                }
                if(m->poar){
                        MFREE(m->poar);
                }
                if(m->avg_len){
                        MFREE(m->avg_len);
                }
                if(m->al){
                        MFREE(m->al);
                }
                MFREE(m);
        }
}
