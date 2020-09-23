/*
	score_sets.c
	
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

struct aln_space* score_sets(struct aln_space* aln_space, struct alignment** alignments,struct parameters* param)
{
	int i,j,c;
	int x,y;
	//float t = param->idal;
	
	for (i = 1; i < (1 << num_alignments);i++){
		aln_space->s[i].score_4_input_aln = 0.0f;
	}
	
	for (i = 0; i < num_alignments;i++){
	
		alignments[i]->score = 1+ (alignments[i]->al * alignments[i]->id*param->idal);
		//fprintf(stderr,"%f	%f:%f\n",alignments[i]->al,alignments[i]->id,alignments[i]->score);
		c = 1 << i;
		for (j = 1; j < (1 << num_alignments);j++){
			if(j & c){
				aln_space->s[j].score_4_input_aln += alignments[i]->score;
			}
		}
	}
	for (i = 1; i < (1 << num_alignments);i++){
		//fprintf(stderr,"%d	%f\n",i,aln_space->s[i].score_4_input_aln);
		aln_space->s[i].score_4_input_aln /= (float)pop(i);
		
		if(pop(i) > 1){
			j = 1 << (num_alignments-1);
			x = num_alignments-1;
			while(j){
				c = j >> 1;
				y = x - 1;
				while(c){
				//	fprintf(stderr,"comparing:%d-%d	%d-%d\n",j,c,x,y);
					if((i & j) && (i & c)){
						aln_space->s[i].aln_sim += aln_space->o[x][y];
					}
					y--;
					c = c >> 1;
				}
				x--;
				j = j >> 1;
			}
			aln_space->s[i].aln_sim /= ((float)pop(i)*(float)(pop(i)-1))/(float)2;
		}
		//exit(0);
		
		//fprintf(stderr,"%d	%f	%f\n",j,aln_space->s[j].pcounts,aln_space->s[j].score_4_input_aln);
		//if(j > 20){
		//	exit(0);
		//}
		//fprintf(stderr,"%d	%f	%f\n",j,aln_space->s[j].pcounts,aln_space->s[j].score_4_input_aln);
		//aln_space->s[j].pcounts *= aln_space->s[j].score_4_input_aln;
		//fprintf(stderr,"%d	%f	%f\n",j,aln_space->s[j].pcounts,aln_space->s[j].score_4_input_aln);
	}
	return aln_space;
}
