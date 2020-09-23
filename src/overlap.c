/*
  overlap.c

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

struct aln_space* calculate_overlap(struct aln_space* aln_space,struct alignment** alignments,char* ref)
{
        //int test = 0;
        int i,j,c;
        //int temp;
        float** o = aln_space->o;
        float tmp = 0.0;
        /*overlap = tmalloc(sizeof(float*)*num_alignments);
          for (i = 0;i < num_alignments;i++){
          overlap[i] = tmalloc(sizeof(float)*num_alignments);
          for ( j = 0; j < num_alignments;j++){
          overlap[i][j] = 0.0f;
          }
          }*/
        for (i = 0;i < num_alignments;i++){
                for ( j = 0; j < num_alignments;j++){
                        aln_space->o[i][j] = 0.0f;
                }
        }

        if(ref){
                for (i = 1;i < (1 << num_alignments);i++){
                        if(i&1){
                                for (c = 1;c < num_alignments;c++){
                                        if(i&(1 << c)){
                                                o[0][c] += aln_space->s[i].pcounts;
                                        }
                                }
                        }
                }
                aln_space->ref_10 = tmalloc(sizeof(float)*num_alignments);
                aln_space->ref_11 = tmalloc(sizeof(float)*num_alignments);
                aln_space->ref_01 = tmalloc(sizeof(float)*num_alignments);
                aln_space->ref_10[0] = 1.0f;
                aln_space->ref_11[0] = 1.0f;
                aln_space->ref_01[0] = 1.0f;
                for (i = 1;i < num_alignments;i++){
                        tmp = aln_space->o[0][i];
                        //fprintf(stderr,"%f\n",tmp);
                        if(alignments[0]->pairs){
                                aln_space->ref_10[i] = (float)tmp/(float)alignments[0]->pairs;//SPS
                        }else{
                                aln_space->ref_10[i] = 0.0f;
                        }

                        if(alignments[i]->pairs+alignments[0]->pairs){
                                aln_space->ref_11[i] = (float)tmp/(((float)alignments[i]->pairs+(float)alignments[0]->pairs)*0.5);//overlap
                        }else{
                                aln_space->ref_11[i] = 0.0f;
                        }

                        if(alignments[i]->pairs){
                                aln_space->ref_01[i] = (float)tmp/(float)alignments[i]->pairs;//%correrctness...
                        }else{
                                aln_space->ref_01[i] = 0.0f;
                        }
                }
        }else{
                for (i = (1 << num_alignments);i--;){
                        for (j = 0; j < num_alignments;j++){
                                if(i&(1 << j)){
                                        for (c =j+1;c < num_alignments;c++){
                                                if(i&(1 << c)){
                                                        o[j][c] += aln_space->s[i].pcounts;
                                                }
                                        }
                                }
                        }
                }
                for (j = 0; j < num_alignments;j++){
                        for (c =j+1;c < num_alignments;c++){
                                o[c][j] = o[j][c];
                        }
                }


                for (i = 0;i < num_alignments;i++){
                        for ( j = 0; j < num_alignments;j++){
                                aln_space->o[i][j] = aln_space->o[i][j]/(((float)alignments[i]->pairs+(float)alignments[j]->pairs)/2);
                        }
                }
        }
        return aln_space;
}



float* find_overlap_to_ref(float* overlap_to_ref,float* sets)
{
        int i = 0;
        for (i = (1 << num_alignments);--i;){
                //printf("%d	%d %f\n",i,pop(i),sets[i]);
                if(i & 1){
                        overlap_to_ref[pop(i)-1] += sets[i];
                }
        }
        return overlap_to_ref;
}
