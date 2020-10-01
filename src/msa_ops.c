/*
  msa_ops.c

  Released under GPL - see the 'COPYING' file

  Copyright (C) 2020 Timo Lassmann <timolassmann@gmail.com>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Please send bug reports, comments etc. to:
  timolassmann@gmail.com
*/

#include "tldevel.h"

#include "global.h"
#include "msa.h"

#include <string.h>
#include <ctype.h>

#include "mumsa_data.h"
#include "msa_info.h"


#define MSA_OPS_IMPORT
#include "msa_ops.h"

#define MOD_ADLER 65521

static unsigned int adler(char* seq,int len);
static int cmp_seq_by_name(const void *a, const void *b);
static int do_alistat_single(struct msa* msa,struct msa_info* msai,double* avg_len_out, double* id_out, double*al_out);

int sanity_check_input(struct mumsa_data* md)
{
        int i;
        int j;
        int c;
        /* check sequence names  */
        for(i = 0; i < md->num_aln;i++){

                for(j = i+1; j <  md->num_aln;j++){
                        if(md->msa[i]->numseq != md->msa[j]->numseq){
                                ERROR_MSG("Alignment %d and %d contain different number of sequences (%d vs %d)",i,j, md->msa[i]->numseq, md->msa[j]->numseq);
                        }
                        for(c = 0; c < md->msai[i]->num_seq;c++){
                                if(md->msai[i]->adler_val[c] != md->msai[j]->adler_val[c]){
                                        ERROR_MSG("Sequence %d and %s differ!",
                                                  md->msa[i]->sequences[c]->name,
                                                  md->msa[j]->sequences[c]->name
                                                );
                                }
                        }
                }
        }
        return OK;
ERROR:
        return FAIL;
}

int sort_msa_by_seqname(struct msa* msa)
{
        int i;
        /* LOG_MSG("Before"); */
        /* for(i = 0; i < msa->numseq;i++){ */
        /*         fprintf(stdout,"%s\n",msa->sequences[i]->name); */
        /* } */
        qsort(msa->sequences,  msa->numseq, sizeof(struct msa_seq*), cmp_seq_by_name);
        /* LOG_MSG("Sorted"); */
        /* for(i = 0; i < msa->numseq;i++){ */
        /*         fprintf(stdout,"%s\n",msa->sequences[i]->name); */
        /* } */
        return OK;
}


int cmp_seq_by_name(const void *a, const void *b)
{

        struct msa_seq* const *ia = a;
        struct msa_seq* const *ib = b;

        return strcmp((*ia)->name, (*ib)->name);
}


int msai_char_to_pos(struct msa_info* msai)
{
        int i;
        int j;
        int c;
        for(i = 0; i < msai->num_seq;i++){
                c = 0;
                for(j = 0; j < msai->aln_len;j++){
                        if(msai->s[i][j] >=0){
                                msai->s[i][j] = c;
                                c++;
                        }
                }
        }
        return OK;
}

int alistat(struct mumsa_data* m)
{
        int i;
        for(i = 0; i < m->num_aln;i++){
                RUN(do_alistat_single(m->msa[i], m->msai[i],&m->avg_len[i],&m->id[i],&m->al[i]));
                /* LOG_MSG("Aln: %f %f %f",m->avg_len[i],m->id[i],m->al[i]); */
        }
        return OK;
ERROR:
        return FAIL;
}

int do_alistat_single(struct msa* msa,struct msa_info* msai,double* avg_len_out, double* id_out, double*al_out)
{
        double id;
        double avg_len;
        double al;
        double poar;

        double l_id;
        double aln_pos;
        int i;
        int j;
        int c;

        id = 0.0;
        al = 0.0;
        avg_len = 0.0;
        poar = 0.0;
        for(i = 0; i < msai->num_seq;i++){
                avg_len += (double) msa->sequences[i]->len;
                for(j = i+1; j < msai->num_seq;j++){
                        if(msa->sequences[i]->len > msa->sequences[j]->len){
                                al += (double) msa->sequences[j]->len;
                        }else{
                                al += (double) msa->sequences[i]->len;
                        }
                        l_id = 0.0;
                        aln_pos = 0.0;
                        for(c = 0; c < msai->aln_len;c++){
                                if(msai->s[i][c] >= 0){
                                        if(msai->s[j][c] >= 0){
                                                if(msai->s[i][c] == msai->s[j][c]){
                                                        l_id += 1.0;
                                                }
                                                aln_pos++;
                                                poar++;
                                        }
                                }
                        }
                        if(!aln_pos){
                                aln_pos = 1.0;
                        }
                        id += l_id / aln_pos;
                }

        }
        id =  id / (double) ((msa->numseq * (msa->numseq-1)) / 2.0);
        al = poar / al;
        avg_len /= (double)msa->numseq;

        *avg_len_out = avg_len;
        *id_out = id;
        *al_out = al;

        return OK;
ERROR:
        return FAIL;
}

int read_msa_into_msai(struct msa* msa,struct msa_info** msai)
{
        struct msa_info* m = NULL;
        uint8_t* tmp;
        int aln_len;
        int len_tmp;
        int i;
        int j;
        int c;
        int g;

        ASSERT(msa != NULL, "No alignment");
        ASSERT(msa->numseq >= 2, "Less than 2 sequences ??");

        /* determine aln_len; */
        aln_len = 0;
        for(j = 0;j < msa->sequences[0]->len;j++){
                aln_len += msa->sequences[0]->gaps[j];
                aln_len++;
        }
        aln_len += msa->sequences[0]->gaps[msa->sequences[0]->len];


        RUN(alloc_msa_info(&m, msa->numseq, aln_len));

        /* LOG_MSG("%d %d", msa->numseq, aln_len+1); */
        for(i = 0; i < msa->numseq;i++){
                g = 0;
                for(j = 0;j < msa->sequences[i]->len;j++){

                        for(c = 0;c < msa->sequences[i]->gaps[j];c++){
                                m->s[i][g] = -1;
                                g++;
                        }
                        //LOG_MSG("%d ",g);
                        m->s[i][g] = msa->sequences[i]->s[j];
                        g++;
                }
                for(c = 0;c < msa->sequences[i]->gaps[msa->sequences[i]->len];c++){
                        m->s[i][g] = -1;
                        g++;
                }
                m->adler_val[i] = (int) adler(msa->sequences[i]->seq, msa->sequences[i]->len);

        }
        *msai = m;
        return OK;
ERROR:

        free_msa_info(m);
        return FAIL;

}



unsigned int adler(char* seq,int len)
{
        unsigned int a = 1;
        unsigned int b = 0;
        unsigned int l = (unsigned int) len;
        while (l) {
                unsigned tlen = l > 5550 ? 5550 : l;
                l -= tlen;
                do {
                        if(isalpha ((int)*seq)){
                                a += (int)toupper(*seq);
                                b += a;
                        }
                        seq++;
                } while (--tlen);
                a = (a & 0xffff) + (a >> 16) * (65536-MOD_ADLER);
                b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);
        }
        if (a >= MOD_ADLER){
                a -= MOD_ADLER;
        }
        b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);
        if (b >= MOD_ADLER){
                b -= MOD_ADLER;
        }
        return (b << 16) | a;
}
