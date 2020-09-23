/*
  main.c

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
#include <time.h>


unsigned int numseq = 0;
unsigned int numprofiles = 0;
unsigned int num_alignments = 0;

int main(int argc, char** argv)
{
        struct alignment** alignments = 0;
        struct aln_space* aln_space = 0;
        int *seed;
        float max = 0.0f;
        int i = 0;
        int c = 0;
        int nseed = (int) -1 * time(NULL);

        seed=(&nseed);

        struct parameters* param = 0;

        param = interface(param,argc, argv);

        alignments = detect_and_read_alignments(alignments,param);

        //Sanity check
        is_alignment(alignments,param);

        for (i = 0; i < num_alignments;i++){
                alignments[i] = change_character_numbers_and_seq_len(alignments[i]);
                alignments[i] = alistat(alignments[i]);
                alignments[i] = calculate_adler(alignments[i]);
        }

        if(param->sort){
                if(byg_start(param->sort,"inputinfileINPUTINFILE") != -1){
                        param->order = get_input_order(param->order,alignments[0]);
                }
        }

        for (i = 0; i < num_alignments;i++){
                quickalnSort(alignments[i],alignments[i]->numseq);
        }
        //Sanity check
        sanity_check(alignments,param);
        aln_space = aln_space_alloc(aln_space);
        if(!param->bootstrap){
                aln_space = make_sets(aln_space,alignments,param->gap_flag);
                aln_space = calculate_overlap(aln_space,alignments,param->ref_flag);

                if(!param->ref_flag){
                        aln_space = score_sets(aln_space,alignments,param);
                        aln_space = calculate_diff(aln_space,alignments,param->set_diff);
                        aln_space = find_best_alignment(alignments,aln_space,param);
                }
                //aln_space = eval_space(alignments,aln_space,param);
        }else{
                aln_space = bootstrap(aln_space,alignments,param);
        }

        if(!param->sort){
                //make tree out of sequence sim scores...
                param->order = make_tree_order(param->order,aln_space->sim,alignments[0]);
        }

        //aln_space = calculate_overlap(aln_space);

        if(param->ref_flag){
                free(aln_space->scores);
                aln_space->diff = -1.0f;

                if (byg_start("11",param->ref_flag) != -1){
                        aln_space->scores = aln_space->ref_11;
                }else if (byg_start("01",param->ref_flag) != -1){
                        aln_space->scores = aln_space->ref_01;
                }else if (byg_start("ID",param->ref_flag) != -1){
                        for (i = 0; i < num_alignments;i++){
                                aln_space->ref_01[i] = alignments[i] ->id;
                        }
                        aln_space->scores = aln_space->ref_01;
                }else if (byg_start("AL",param->ref_flag) != -1){
                        for (i = 0; i < num_alignments;i++){
                                aln_space->ref_01[i] = alignments[i] ->al;
                        }
                        aln_space->scores = aln_space->ref_01;
                }else{
                        aln_space->scores = aln_space->ref_10;

                }
                print_alignment_scores(aln_space,param,alignments);
                aln_space->scores = 0;

        }else if(param->print_alignment_flag){
                /*aln_space = find_best_alignment(alignments,aln_space,param);*/
                max = -1.0f;
                for (i = 0; i < num_alignments;i++){
                        if(aln_space->scores[i] > max){
                                max = aln_space->scores[i];
                                c = i;
                        }
                }
                alignments[c] = clean_alignment(alignments[c],param->order);
                output(alignments[c],param);
        }else if(param->best_flag){
                /*##########################################
                  Prints out (hopefully) an alignment consistning only the most reliable POAR's.
                  ###########################################*/
                // puts new alignment 'over' alignemnt 0...
                alignments = get_reliable_blocks(alignments,aln_space,param);


                alignments[0] = alistat(alignments[0]);
                fprintf(stderr,"ID:%f	AL:%f   SCORE:%f\n",alignments[0]->id,alignments[0]->al,alignments[0]->id*alignments[0]->al);
                alignments[0] = clean_alignment(alignments[0],param->order);
                output(alignments[0],param);
        }else{
                //aln_space->diff = 0.0f;
                /*c = 0;

                  for (i = 0;i < num_alignments;i++){
                  for ( j = i+1; j < num_alignments;j++){
                  //aln_space->o[j][i] = aln_space->o[i][j];
                  aln_space->diff += aln_space->o[i][j];
                  c++;
                  }
                  }
                  aln_space->diff = (float)aln_space->diff/(float)(c);

                  aln_space = find_best_alignment(alignments,aln_space,param);*/

                //Try to produce R output (Relation between alignments)
                if(!param->server_flag){
                        r_statistics_output(param->infile,aln_space->scores,aln_space->o,num_alignments);


                        float* seq_scores = tmalloc(sizeof(float)*numseq);
                        char** names = tmalloc(sizeof(char*) * numseq);

                        for (i =0; i < numseq;i++){
                                names[i] = alignments[0]->sn[i];
                                seq_scores[i] = aln_space->sim[i][i];
                        }
                        r_statistics_output(names,seq_scores,aln_space->sim,numseq);
                        free(names);
                        free(seq_scores);
                }

                //Try to produce R output (Relation between sequences );
                /*char** seq_names = 0;
                  float* seq_scores = tmalloc(sizeof(float)*numseq);
                  for (i =0; i < numseq;i++){
                  seq_scores[i] = aln_space->sim[i][i];
                  }
                  seq_names = get_seq_names(seq_names,aln_name[0]);


                  if(!param->server_flag){
                  r_statistics_output(seq_names,seq_scores,aln_space->sim,numseq);
                  }
                  free(seq_scores);
                  for (i =0; i < numseq;i++){
                  free(seq_names[i]);
                  }
                  free(seq_names);*/

                //Print out default output -> alignment scores....
                print_alignment_scores(aln_space,param,alignments);
        }

        free_aln_space(aln_space);

        for (i = 0; i < num_alignments;i++){
                free_aln(alignments[i]);
        }
        free(alignments);

        free_param(param);
        return 0;
}
