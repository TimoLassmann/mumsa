/*
 	best.c
	
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

/*
struct aln_space* eval_space(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param)
{
	int i,j,c;
	c = 0;
	
	for (i = 0;i < num_alignments;i++){
		for ( j = i+1; j < num_alignments;j++){
				//aln_space->o[j][i] = aln_space->o[i][j];
			aln_space->diff += aln_space->o[i][j];
			c++;
		}
	}
	aln_space->diff = (float)aln_space->diff/(float)(c);

	aln_space = find_best_alignment(alignments,aln_space,param);
	
	return aln_space;
}*/

struct aln_space* calculate_diff(struct aln_space* aln_space,struct alignment** alignments,int cutoff)
{
	float average_poar = 0.0f;
	int i,j,c;
	
	aln_space->diff = 0.0f;
	
	c = 0;
	
	
	if(cutoff){
		if(cutoff > num_alignments){
			cutoff = num_alignments;
		}
		
		for (i = 1;i < (1 << num_alignments);i++){
			if(pop(i) >= cutoff){
				aln_space->diff +=  aln_space->s[i].pcounts;
			}
		}
		for (i = 0;i < num_alignments;i++){
			average_poar += alignments[i]->pairs;
		}

		average_poar /= num_alignments;
		
		average_poar = average_poar * (num_alignments - (cutoff-1));
		if(aln_space->diff != 0.0f &&  average_poar != 0.0f){
			aln_space->diff = aln_space->diff/average_poar;
		}
		//fprintf(stderr,"%f	= %f\n",aln_space->diff,average_poar);
		
	}else{
		for (i = 0;i < num_alignments;i++){
			for ( j = i+1; j < num_alignments;j++){
				aln_space->diff += aln_space->o[i][j];
				c++;
			}
		}
		aln_space->diff = (float)aln_space->diff/(float)(c);
	}
	
	
	return aln_space;
}


struct aln_space* find_best_alignment(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param)
{
	float* scores = aln_space->scores;
	
	float* mos = 0;
	float* aos = 0;
	float* idal = 0;
	//int* num = 0;
		
	int i = 0;
	int j = 0;
	int c = 0;
	
	
	mos = tmalloc(sizeof(float)*num_alignments);
	aos = tmalloc(sizeof(float)*num_alignments);
	idal = tmalloc(sizeof(float)*num_alignments);
	
	//num = tmalloc(sizeof(int)*num_alignments);
	
	
	for (i = 0; i < num_alignments;i++){
		mos[i] = 0.0;
		aos[i] = 0.0;
		idal[i] = 0.0;
	//	num[i] = 0;
	}
	
	//score Alignments
	for (i = 1;i < (1 << num_alignments);i++){
		j = 1 << (num_alignments-1);
		c = num_alignments-1;
		while(j){
			if (i & j){
				//Shall we include aligned Pairs of residues without support? Possibly increases the score of global methods
				//if(pop(i) > 1){
					mos[c] += aln_space->s[i].pcounts*(pop(i)-1)  * (1-aln_space->s[i].aln_sim*param->x1) * aln_space->s[i].score_4_input_aln;
					//aos[c] += aln_space->s[i].aln_sim;
					//idal[c] += aln_space->s[i].score_4_input_aln;
					//num[c] += 1;
				//}
			}
			j = j >>1;
			c--;
		}
	}
	//exit(0);
	for (i = 0;i < num_alignments;i++){
		mos[i] /= ((num_alignments-1) *(float)alignments[i]->pairs);
		//aos[i] /= (1 << num_alignments);
		//idal[i] /= num_alignments;
	}
	
	for (i = 0;i < num_alignments;i++){
		scores[i] = mos[i];// + param->x2 * aos[i] + param->x3 * idal[i];
		//fprintf(stderr,"scores[i] = param->x1 * mos[i] + param->x2 * aos[i] + param->x3 * idal[i]\n%f = %f * %f + %f * %f + %f * %f\n",scores[i],param->x1,mos[i],param->x2,aos[i],param->x3,idal[i]);
	}
	
	free(mos);
	free(aos);
	free(idal);
	return aln_space;
}
