/*
	mem.c
	
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

#ifndef MEMORY 
void* tmalloc(int size){
	void* p;
	p = (void*)malloc(size);
	if (!p){
		fprintf(stderr,"Out of memory!\n");
		exit(0);
	}
	return p;
}
#endif

struct alignment* aln_alloc(struct alignment* aln)
{
	int i;
	aln = (struct alignment*) tmalloc(sizeof(struct alignment));
	
	aln->id = 0.0f;
	aln->av_len = 0.0f;
	aln->al = 0.0f;
	aln->score = 0.0f;
	
	aln->s = tmalloc(sizeof(int*) * (numseq ));
	aln->seq = tmalloc(sizeof(char*) * (numseq ));
	aln->ft =  tmalloc(sizeof(struct feature* ) * (numseq));
	aln->si  =  tmalloc(sizeof(struct sequence_information* ) * (numseq));
	aln->sl = tmalloc(sizeof(unsigned int) * (numprofiles));
	aln->sip = tmalloc(sizeof(unsigned int*)* numprofiles);
	
	aln->nsip = tmalloc(sizeof(unsigned int)* numprofiles);
	aln->sn = tmalloc(sizeof(char*) * numseq);
	aln->lsn = tmalloc(sizeof(unsigned int) * numseq);
	aln->adler = tmalloc(sizeof(unsigned int) * numseq);
	
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
		aln->sl[i] = 0;
	}
	
	for(i =0;i < numseq;i++){
		aln->adler[i] = 0;
		aln->lsn[i] = 0;
		aln->ft[i] = 0;
		aln->si[i] = 0;
		aln->sip[i] = tmalloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}
	return aln;
}

struct aln_space* aln_space_alloc(struct aln_space* aln_space)
{
	int i,j;
	aln_space = tmalloc(sizeof(struct aln_space));
	
	aln_space->s = tmalloc(sizeof(struct sets)* (1 << num_alignments));

	aln_space->sim = tmalloc(sizeof(float*)*numseq);
	aln_space->o =  tmalloc(sizeof(float*)*num_alignments);
	
	aln_space->scores = tmalloc(sizeof(float)*num_alignments);
	aln_space->support = 0;
	aln_space->ref_10 = 0;
	aln_space->ref_11 = 0;
	aln_space->ref_01 = 0;
	
	aln_space->diff = -1.0f;


	for (i = numseq;i--;){
		aln_space->sim[i] = tmalloc(sizeof(float)*numseq);
	}
	
	for (i = 0; i < num_alignments;i++){
		aln_space->scores[i] = 0.0f;
		aln_space->o[i] = tmalloc(sizeof(float)*num_alignments);
		for (j = 0; j < num_alignments;j++){
			aln_space->o[i][j] = 0.0f;
		}
	}
	
	for (i =0; i < numseq;i++){
		for (j = 0; j < numseq;j++){
			aln_space->sim[i][j]= 0.0f;
		}
	}
	
	for (i = (1 << num_alignments);i--;){
		aln_space->s[i].score_4_input_aln = 0.0f;
		aln_space->s[i].aln_sim = 0.0f;
		aln_space->s[i].id= 0.0f;
		aln_space->s[i].pcounts = 0.0f;
	}
	return aln_space;
}

void free_aln_space(struct aln_space* aln_space)
{
	int i;
	for (i = 0; i < num_alignments;i++){
		free(aln_space->o[i]);
	}
	
	for (i = numseq;i--;){
		free(aln_space->sim[i]);
		
	}
	
	if(aln_space->ref_10){
		free(aln_space->ref_10);
		free(aln_space->ref_11);
		free(aln_space->ref_01);
	}
	if( aln_space->scores){
		free(aln_space->scores);
	}
	if( aln_space->support){
		free(aln_space->support);
	}
	free(aln_space->o);
	free(aln_space->s);
	free(aln_space->sim);
	free(aln_space);
}


void free_aln(struct alignment* aln)
{
	int i;
	for (i = aln->numseq;i--;){
		free(aln->s[i]);
		free(aln->seq[i]);
		free(aln->sn[i]);
	}

	if(aln->ft){
		for(i = aln->numseq;i--;){
			free_ft(aln->ft[i]);
		}
		free(aln->ft);
	}
	if(aln->si){
		free(aln->si);
	}

	for (i = aln->numseq;i--;){
		if(aln->sip[i]){
			free(aln->sip[i]);
		}
	}
	free(aln->adler);
	free(aln->seq);
	free(aln->s);
	free(aln->sn);
	free(aln->sl);
	free(aln->lsn);
	free(aln->sip);
	free(aln->nsip);
	free(aln);
}

void free_columns(struct columns* columns)
{
	int i;
	for (i = 0;i < columns->len;i++){
		free(columns->c[i]);
	}
	free(columns->c);
	free(columns);
}

void free_ft(struct feature* n)
{
	struct feature* old_n = 0;
	if (n != NULL){
	 	old_n = n;
	 	n= n ->next;
	 	free(old_n->type);
	 	free(old_n->note);
 	 	free(old_n);
		free_ft(n);
	}
}


void free_param(struct parameters* param)
{
	if (param->order){
		free(param->order);
	}
	if(param->infile){
		free(param->infile);
	}
	free(param);
}

