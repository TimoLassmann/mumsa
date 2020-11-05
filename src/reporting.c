/*
  reporting.c

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
#include "mumsa_data.h"

#define REPORTING_IMPORT
#include "reporting.h"


int print_scores(struct mumsa_data* m,struct parameters* param)
{
        int i;
        double* scores = NULL;
        double n_pairs;
        MMALLOC(scores,sizeof(double) * m->num_aln);
        scores[0] = 1.0;
        for(i = 1; i < m->num_aln;i++){
                scores[i] = m->overlap[0][i];
                n_pairs = 0.0;
                if(param->score_mode == MUMSA_SCORE_REF){
                        n_pairs = m->msai[0]->pairs;
                }else if(param->score_mode == MUMSA_SCORE_REF_TEST){
                        n_pairs = (double)(m->msai[0]->pairs + m->msai[i]->pairs) / 2.0;
                }else{
                        n_pairs = m->msai[i]->pairs;
                }
                if(n_pairs == 0.0){
                        scores[i] = 0.0;
                }else{
                        scores[i] = scores[i] / n_pairs;
                }
        }
        for(i = 1; i < m->num_aln;i++){
                fprintf(stdout,"auto %s %f\n", param->infile[i],scores[i]);
        }
        MFREE(scores);

        return OK;
ERROR:
        if(scores){
                MFREE(scores);
        }
        return FAIL;
}
