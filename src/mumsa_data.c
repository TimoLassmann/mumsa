/*
  mumsa_data.c

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
        m->s = NULL;
        m->overlap = NULL;

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


        MMALLOC(m->s, sizeof(struct sets) * (1 << m->num_aln));
        for(i = 0; i < (1 << m->num_aln);i++){
                m->s[i].id = 0.0;
                m->s[i].aln_sim = 0.0;
                m->s[i].pcounts = 0.0;
                m->s[i].score_4_input_aln = 0.0;
        }

        RUN(galloc(&m->overlap,m->num_aln,m->num_aln));

        *mdat = m;
        return OK;
ERROR:
        return FAIL;
}


int free_mumsa_data(struct mumsa_data* m)
{
        int i;
        if(m){
                if(m->s){
                        MFREE(m->s);
                }
                if(m->sim){
                        gfree(m->sim);
                }
                if(m->overlap){
                        gfree(m->overlap);
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
