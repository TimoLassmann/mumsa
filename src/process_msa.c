#include "tldevel.h"

#include "khash.h"

#include "mumsa_data.h"

#define PROCESS_MSA_IMPORT
#include "process_msa.h"


KHASH_MAP_INIT_INT64(PAIR,int);

static int fill_hash(khash_t(PAIR) *h, struct msa_info* msai,int aln,int i,int j);
static int pop(unsigned int x);


int calc_overlap(struct mumsa_data* m)
{
        double p1;
        double p2;
        int i;
        int j;
        int c;

        for(i = 0; i < m->num_aln;i++){
                for(j = 0;j < m->num_aln;j++){
                        m->overlap[i][j] = 0.0;
                }
        }
        for (i = (1 << m->num_aln);i--;){
                for (j = 0; j < m->num_aln;j++){
                        if(i&(1 << j)){
                                for (c =j+1;c < m->num_aln;c++){
                                        if(i&(1 << c)){
                                                m->overlap[j][c] += m->s[i].pcounts;
                                        }
                                }
                        }
                }
        }
        for (j = 0; j < m->num_aln;j++){
                for (c =j+1;c < m->num_aln;c++){
                        LOG_MSG("%d %d %f", j,c, m->overlap[j][c]);
                        m->overlap[c][j] = m->overlap[j][c];
                }
        }


        /*for (i = 0;i < m->num_aln;i++){
                p1 = m->msai[i]->pairs;

                for ( j = 0; j < m->num_aln;j++){
                        p2 = m->msai[j]->pairs;

                        m->overlap[i][j] = m->overlap[i][j]/ ((p1 + p2) / 2.0);
                        LOG_MSG("%d %d %f %f -> %f", i,j,p1,p2, m->overlap[i][j]);
                }
                }*/
        return OK;
}

int calc_sim_pairs(struct mumsa_data* m)
{
        khash_t(PAIR) *h = kh_init(PAIR);
        khiter_t k;
        double tmp;
        int len_a;
        int len_b;

        int i;
        int j;
        int c;
        int v;
        /* to make sure  */
        for (i = 0;i < (1 << m->num_aln);i++){

                m->s[i].pcounts = 0.0;
        }
        for(i = 0; i < m->num_aln;i++){
                m->msai[i]->pairs = 0.0;
        }
        for(i = 0; i < m->num_seq;i++){
                for (j = 0; j < m->num_seq;j++){
                        m->sim[i][j] = 0.0;
                }
        }


        for(i = 0; i < m->num_seq-1;i++){
                len_a = m->msa[0]->sequences[i]->len;
                for(j = i+1; j < m->num_seq;j++){
                        len_b =  m->msa[0]->sequences[j]->len;
                        kh_clear(PAIR, h);
                        for(c = 0; c < m->num_aln;c++){
                                fill_hash(h, m->msai[c],c, i, j);
                        }
                        for (k = kh_begin(h); k != kh_end(h); ++k){
                                if (kh_exist(h, k)){

                                        v = kh_value(h, k);

                                        m->s[v].pcounts += 1.0;
                                        m->sim[i][j] += pop(v);
                                }
                        }
                        //LOG_MSG("%d %d %f / %f ", i,i, m->sim[i][j], ((double) m->num_aln *(double)(MACRO_MIN(len_a,len_b))));
                        m->sim[i][j] = m->sim[i][j] / ((double) m->num_aln *(double)(MACRO_MIN(len_a,len_b)));
                        //LOG_MSG("%d %d %f", i,i, m->sim[i][j]);
                        kh_clear(PAIR, h);
                }
        }
        for(i = 0; i < m->num_seq-1;i++){
                for (j = i + 1; j < m->num_seq;j++){
                        m->sim[j][i] = m->sim[i][j];
                }
        }

        for(i = 0; i < m->num_seq-1;i++){
                tmp = 0.0;
                for (j = 0; j < m->num_seq;j++){
                        tmp+= m->sim[i][j];
                }
                m->sim[i][i] = tmp / (double) m->num_seq;
        }

        kh_destroy(PAIR, h);
        return OK;
}

int fill_hash(khash_t(PAIR) *h, struct msa_info* msai,int aln,int i,int j)
{
        uint64_t key;
        khiter_t k;
        int ret;
        int c;

        for(c = 0; c < msai->aln_len;c++){
                if(msai->s[i][c] != -1 && msai->s[j][c] != -1){
                        key = (uint64_t) msai->s[i][c] << 32UL | (uint64_t) msai->s[j][c];

                        k = kh_put(PAIR, h, key, &ret);
                        if (!ret){
                                kh_value(h, k) |= (1 << aln);
                        }else{
                                kh_value(h, k) = 1 << aln;
                        }
                        msai->pairs += 1.0;
                }
        }

        return OK;
}



int pop(unsigned int x)
{
        static char table[256] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
        return table[x & 0xff] + table[(x >> 8) & 0xff] + table[(x >> 16) & 0xff] + table[(x >> 24)];
}
