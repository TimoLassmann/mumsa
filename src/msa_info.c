/*
  msa_info.c

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
        msai->pairs = 0.0;
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
