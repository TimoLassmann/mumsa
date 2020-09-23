/*
  boot.c

  Released under GPL - see the 'COPYING' file

  Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>

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


#include "cmsa.h"
#include "boot.h"

float cutoff = 0.0;

struct aln_space* bootstrap(struct aln_space* aln_space, struct alignment** alignments,struct parameters* param)
{
        float* scores = 0;
        float max = -1.0f;

        int i,j,c;

        scores = tmalloc(sizeof(float)*num_alignments);
        aln_space->support = tmalloc(sizeof(float)*num_alignments);
        for (i = 0;i < num_alignments;i++){
                scores[i] = 0.0f;
                aln_space->support[i] = 0.0f;
        }


        cutoff = 1.0f-param->bootcutoff;

        //fprintf(stderr,"CUTOFF:%f	%f\n",cutoff,1+cutoff);

        aln_space = make_sets(aln_space,alignments,param->gap_flag);

        fprintf(stderr,"\nBootstrapping:\n");
        for (i = 0; i < param->bootstrap;i++){
                fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)param->bootstrap*100);
                aln_space = make_boot_sets(aln_space,alignments, param->gap_flag);
                aln_space = calculate_overlap(aln_space,alignments,param->ref_flag);
                aln_space = score_sets(aln_space,alignments,param);
                aln_space = find_best_alignment(alignments,aln_space,param);
                max = -1.0f;
                c = -1;
                for (j = 0;j <  num_alignments;j++){
                        scores[j] += aln_space->scores[j];
                        if(aln_space->scores[j] > max){
                                max = aln_space->scores[j];
                                c = j;
                        }
                }
                aln_space->support[c] += 1;
        }
        fprintf(stderr,"\r%8.0f percent done\n",100.0f);




        for (i = 0; i <  num_alignments;i++){
                //	aln_space->scores[i] = scores[i] / (float) param->bootstrap;
                //	aln_space->scores[i] *= (1.00f + cutoff);
                aln_space->support[i] = aln_space->support[i] / (float) param->bootstrap;
        }

        aln_space = make_sets(aln_space,alignments,param->gap_flag);
        aln_space = calculate_overlap(aln_space,alignments,param->ref_flag);
        aln_space = score_sets(aln_space,alignments,param);
        aln_space = calculate_diff(aln_space,alignments,param->set_diff);
        aln_space = find_best_alignment(alignments,aln_space,param);


        /*aln_space = make_sets(aln_space,alignments,param->gap_flag);

          aln_space = calculate_overlap(aln_space,alignments,param->ref_flag);
          //aln_space = score_sets(aln_space,alignments,param);
          aln_space = calculate_diff(aln_space);
          aln_space = find_best_alignment(alignments,aln_space,param);
          //aln_space = find_best_alignment(alignments,aln_space,param);
          for (j = 0;j <  num_alignments;j++){
          fprintf(stderr,"%f ",aln_space->scores[j]);
          }
          fprintf(stderr,"\n");*/

        free(scores);
        return aln_space;
}

struct aln_space* make_boot_sets(struct aln_space* aln_space, struct alignment** alignments, int gap_flag)
{
        struct node** hash = 0;
        struct node* old_n = 0;
        struct node* n = 0;


        int i,j,c;
        int hashsize = 0;


        for ( i = 0; i < num_alignments;i++){
                alignments[i]->pairs = 0;
                if (hashsize < alignments[i]->len){
                        hashsize = alignments[i]->len;
                }
        }
        hashsize++;

        hash = tmalloc(sizeof(struct node*) * hashsize);
        for (i = hashsize;i--;){
                hash[i] = 0;
        }
        for (i = (1 << num_alignments);i--;){
                aln_space->s[i].pcounts = 0.0f;
        }


        for (i = 0; i < numseq;i++){
                for (j = i +1; j < numseq;j++){
                        for (c = 0; c < num_alignments;c++){
                                hash = feed_boot_hash(hash,alignments[c],i,j,c,gap_flag);
                        }

                        for (c = hashsize;c--;){
                                n = hash[c];
                                while (n){
                                        aln_space->s[n->group].pcounts++;
                                        old_n = n;
                                        n = n->next;
                                        free(old_n);
                                }
                                hash[c] = 0;
                        }
                }
        }

        free(hash);

        return aln_space;
}

struct node** feed_boot_hash(struct node** hash, struct alignment* ap, unsigned int a ,unsigned int b, unsigned int aln, int gap_flag)
{
        unsigned int i;
        unsigned int p = 0;
        unsigned int len = ap->len;
        int* seqa = 0;
        int* seqb = 0;

        seqa = ap->s[a];
        seqb = ap->s[b];
        if (!gap_flag){
                for (i = 0; i < len;i++){
                        if (seqa[i] < len && seqb[i] < len){
                                //if(ran1(seed) > cutoff){
                                hash[seqa[i]] = insert_pair(seqb[i],aln,hash[seqa[i]]);
                                        p++;
                                        //}
                        }
                }
        }else{
                for (i = 0; i < len;i++){
                        /*if (seqa[i] < len && seqb[i] < len){
                          hash[seqa[i]] = insert(seqb[i],aln,hash[seqa[i]]);
                          p++;
                          }else if (seqa[i] == len && seqb[i] < len){
                          hash[gap_flag] = insert(seqb[i],aln,hash[len]);
                          p++;
                          }else if (seqb[i] == len && seqa[i] < len){
                          hash[seqa[i]] = insert(gap_flag,aln,hash[seqa[i]]);
                          p++;
                          }*/
                        if (seqa[i] < len){
                                if (seqb[i] < len){
                                        hash[seqa[i]] = insert_pair(seqb[i],aln,hash[seqa[i]]);
                                }else{
                                        hash[seqa[i]] = insert_pair(gap_flag,aln,hash[seqa[i]]);
                                }
                                p++;
                        }else if (seqa[i] == len && seqb[i] < len){
                                hash[gap_flag] = insert_pair(seqb[i],aln,hash[len]);
                                p++;
                        }
                }
        }
        ap->pairs += p;
        return hash;
}
