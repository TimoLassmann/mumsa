/*
  misc.c

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
#include <string.h>

void is_alignment(struct alignment** alignments,struct parameters* param)
{
        int i,j,c;
        int conflict = 1;

        for (i = 0; i < num_alignments;i++){
                for (j = 0 ; j < alignments[i]->numseq;j++){
                        for (c = j+1; c < alignments[i]->numseq;c++){
                                if(strlen(alignments[i]->seq[j]) != strlen(alignments[i]->seq[c])){
                                        fprintf(stderr,"\nWARNING %d:\nSequence %d is not aligned to sequence %d in alignment '%s'.\n",conflict,j,c,param->infile[i]);
                                        conflict++;
                                }
                        }
                }
        }

        if (conflict > 1){
                for (i = 0; i < num_alignments;i++){
                        free_aln(alignments[i]);
                }
                free(alignments);

                free_param(param);
                exit(-1);
        }
}

void sanity_check(struct alignment** alignments,struct parameters* param)
{
        int i,j,c;
        int conflict = 1;
        int min_numseq = 100000000;

        for (i = 0; i < num_alignments;i++){
                if(alignments[i]->numseq < min_numseq){
                        min_numseq = alignments[i]->numseq;
                }
                for (j = i+1;j < num_alignments;j++){
                        if(alignments[i]->numseq != alignments[j]->numseq){
                                fprintf(stderr,"\nWARNING %d:\nAlignment '%s' contains %d sequences\nbut:\nAlignment '%s' contains %d sequences\n",conflict,param->infile[i],alignments[i]->numseq,param->infile[j],alignments[j]->numseq);
                                conflict++;
                        }
                }
        }

        for (i = 0; i < min_numseq;i++){
                for (j = 0;j < num_alignments;j++){
                        for (c = j + 1;c < num_alignments;c++){
                                if(alignments[j]->adler[i] != alignments[c]->adler[i]){
                                        fprintf(stderr,"\nWARNING %d:	%d	%d\nSequence %d (%s) from alignment '%s'\nis not identical to:\nSequence %d (%s) from alignment '%s'\n",conflict,alignments[j]->adler[i],alignments[c]->adler[i],i+1,alignments[j]->sn[i],param->infile[j],i+1,alignments[c]->sn[i],param->infile[c]);
                                        conflict++;
                                }

                        }
                }
        }

        if (conflict >  1){
                for (i = 0; i < num_alignments;i++){
                        free_aln(alignments[i]);
                }
                free(alignments);

                free_param(param);
                exit(-1);
        }
}

int* make_tree_order(int* order,float** sim,struct alignment* aln)
{
        float** dm = 0;
        int* tree = 0;
        int i,j,c;
        order = tmalloc(sizeof(int)*aln->numseq);
        tree = tmalloc(sizeof(int)*(numseq*3));

        for ( i = 0; i < (numseq*3);i++){
                tree[i] = 0;
        }

        dm = tmalloc(sizeof(float*)*numseq);
        for (i = 0;i < numseq; i++){
                dm[i] = tmalloc(sizeof(float)*numseq);
                for ( j = 0;j < numseq;j++){
                        dm[i][j] = sim[i][j];
                }
        }

        tree = upgma(dm,tree);

        c = 0;
        for (i = 0; i < (numseq-1)*3;i +=3){
                if(tree[i]  < numseq){
                        order[c] = aln->adler[tree[i]];
                        c++;
                }
                if(tree[i+1]  < numseq){
                        order[c] = aln->adler[tree[i+1]];
                        c++;
                }
        }
        free(tree);

        for (i = 0;i < numseq;i++){
                free(dm[i]);
        }
        free(dm);

        return order;
}

int* get_input_order(int* input_order,struct alignment* aln)
{
        int i;
        input_order = tmalloc(sizeof(int)*aln->numseq);
        for (i = 0; i < aln->numseq;i++){
                input_order[i] = aln->adler[i];
                //	fprintf(stderr,"%d	%d\n",input_order[i] ,aln->adler[i]);
        }
        return input_order;
}

void quickalnSort(struct alignment* aln, int array_size)
{
        q_sort(aln, 0, array_size - 1);
}

void q_sort(struct alignment* aln, int left, int right)
{
        int l_hold, r_hold;

        int pivot = aln->adler[left];
        //int pivot2;

        struct feature* pivot_ft = aln->ft[left];
        struct sequence_info* pivot_si = aln->si[left];
        unsigned int* pivot_sip = aln->sip[left];
        unsigned int pivot_nsip = aln->nsip[left];
        unsigned int pivot_sl = aln->sl[left];
        unsigned int pivot_lsn =  aln->lsn[left];
        int* pivot_s = aln->s[left];
        char* pivot_seq = aln->seq[left];
        char* pivot_sn = aln->sn[left];


        l_hold = left;
        r_hold = right;
        //pivot2 = aln->nsip[left];
        //pivot = aln->sip[left][0];// numbers[left];
        //pivot = aln=->adler[left];
        while (left < right){
                while ((aln->adler[right] <= pivot) && (left < right)){
                        right--;
                }
                if (left != right){
                        aln->adler[left] = aln->adler[right];

                        aln->ft[left] = aln->ft[right];
                        aln->si[left] = aln->si[right];
                        aln->sip[left] = aln->sip[right];
                        aln->nsip[left] = aln->nsip[right];
                        aln->sl[left] = aln->sl[right];
                        aln->lsn[left] = aln->lsn[right];
                        aln->s[left] = aln->s[right];
                        aln->seq[left] = aln->seq[right];
                        aln->sn[left] = aln->sn[right];

                        left++;
                }
                while ((aln->adler[left] >= pivot) && (left < right)){
                        left++;
                }
                if (left != right){
                        aln->adler[right] = aln->adler[left];

                        aln->ft[right] = aln->ft[left];
                        aln->si[right] = aln->si[left];
                        aln->sip[right] = aln->sip[left];
                        aln->nsip[right] = aln->nsip[left];
                        aln->sl[right] = aln->sl[left];
                        aln->lsn[right] = aln->lsn[left];
                        aln->s[right] = aln->s[left];
                        aln->seq[right] = aln->seq[left];
                        aln->sn[right] = aln->sn[left];

                        //aln->nsip[right] = aln->nsip[left];
                        right--;
                }
        }
        aln->adler[left] = pivot;

        aln->ft[left] = pivot_ft;
        aln->si[left] = pivot_si;
        aln->sip[left] = pivot_sip;
        aln->nsip[left] = pivot_nsip;
        aln->sl[left] = pivot_sl;
        aln->lsn[left] = pivot_lsn;
        aln->s[left] = pivot_s;
        aln->seq[left] = pivot_seq;
        aln->sn[left] = pivot_sn;


        pivot = left;
        left = l_hold;
        right = r_hold;
        if (left < pivot){
                q_sort(aln, left, pivot-1);
        }
        if (right > pivot){
                q_sort(aln, pivot+1, right);
        }
}


struct alignment* change_character_numbers_and_seq_len(struct alignment* aln)
{
        int i,j,n;
        for (i = 0; i < aln->numseq;i++){
                n = 0;
                for (j = 0; j < aln->len;j++){
                        if(aln->s[i][j] >= 0){
                                aln->s[i][j] = n;
                                n++;
                        }else{
                                aln->s[i][j] = aln->len;
                        }
                }
                aln->sl[i] = n;
        }
        return aln;
}

int pop(unsigned int x)
{
        static char table[256] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
        return table[x & 0xff] + table[(x >> 8) & 0xff] + table[(x >> 16) & 0xff] + table[(x >> 24)];
}
