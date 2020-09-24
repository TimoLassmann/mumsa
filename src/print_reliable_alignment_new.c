/*
  print_reliable_alignment.c

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
#include "print_reliable_alignment.h"

struct alignment** get_reliable_blocks(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param)
{
        struct columns* columns = 0;

        float** ann = 0;
        int criterion = 0;
        int i;

        columns =  make_initial_columns(alignments[0],columns);

        criterion = find_criterion(aln_space,param->best_flag);

        ann = make_columns(ann,columns,alignments,criterion);

        //quickSort_basic(columns->c, columns->len);

        //alignments[0] = make_alignment(columns,alignments[0]);

        alignments[0] = inclusion_sort(columns->c,columns->len,alignments[0],param->print_column_cutoff,param->print_spacer);

        //add annotation.....
        if(param->format){
                if (byg_start("macsim",param->format) != -1){
                        alignments[0] = add_annotation(alignments[0],ann);
                }
        }

        free_columns(columns);

        for (i = 0; i < numseq;i++){
                free(ann[i]);
        }
        free(ann);

        return alignments;
}

struct alignment* make_alignment(struct columns* columns,struct alignment* aln)
{
        struct column_node* start = 0;
        struct column_node* n = 0;
        int i,j;
        int aln_len = 0;
        struct column_node* tmp = 0;

        start = tmalloc(sizeof(struct column_node));
        start->column = tmalloc(sizeof(int)*numseq);

        for (i =0; i < numseq;i++){
                start->column[i] = 0;
        }
        start->next = 0;

        int * tc = tmalloc(sizeof(int)*numseq);
        for (i = 0; i < columns->len;i++){
                if(columns->c[i][0] < 0){
                        break;
                }
                aln_len++;
                for (j = 0; j < numseq;j++){
                        tc[j] = columns->c[i][j+1];
                }
                insert_column(tc,start,start->next);
        }

        aln->len = aln_len;

        for (i = 0; i < numseq;i++){
                free(aln->seq[i]);
                free(aln->s[i]);
                aln->seq[i] = tmalloc(sizeof(char)*(aln_len+1));
                aln->s[i] = tmalloc(sizeof(int)*(aln_len+1));
        }
        n = start->next;
        i = 0;
        while(n){
                tmp = n;
                for (j = 0; j < numseq;j++){
                        if ((char)(n->column[j]&0x000000ff)){
                                aln->seq[j][i] = (char)(n->column[j]&0x000000ff);
                        }else{
                                aln->seq[j][i] = '-';
                        }
                }
                i++;
                n = n->next;
                free(tmp->column);
                free(tmp);
                //free(c[i]);
        }
        //free(c);

        for (i = 0; i < numseq;i++){
                aln->seq[i][aln_len] = 0;
        }

        free(tc);
        free(start->column);
        free(start);
        return aln;
}

struct alignment* add_annotation(struct alignment* aln,float** ann)
{
        int i,j;
        int start = 0;
        int stop = 0;
        int old_ann = 0;

        /*for (i = 0; i < numseq;i++){
          for (j = 0; j <  aln->sl[i];j++){
          if(ann[i][j]){
          ann[i][j] = 1;
          }
          if(ann[i][j]  < 0.1){
          ann[i][j] = 0;
          }else if(ann[i][j]  < 0.2){
          ann[i][j] = 20;
          }else if(ann[i][j]  < 0.3){
          ann[i][j] = 30;
          }else if(ann[i][j]  < 0.4){
          ann[i][j] = 40;
          }else if(ann[i][j]  < 0.5){
          ann[i][j] = 50;
          }else if(ann[i][j]  < 0.6){
          ann[i][j] = 60;
          }else if(ann[i][j]  < 0.7){
          ann[i][j] = 70;
          }else if(ann[i][j]  < 0.8){
          ann[i][j] = 80;
          }else if(ann[i][j]  < 0.9){
          ann[i][j] = 90;
          }else{
          ann[i][j] = 100;
          }
          }
          }*/


        /*for (i = 0; i < numseq;i++){
          fprintf(stderr,"%d	%d\n",i,aln->sl[i]);
          for (j = 0; j <  aln->len;j++){
          fprintf(stderr,"%d ",aln->s[i][j]);
          }
          fprintf(stderr,"\n");
          }
          fprintf(stderr,"\n");
        */

        for (i = 0; i < numseq;i++){
                start = -1;
                stop = 0;
                old_ann = 0;
                for (j = 0; j <  aln->sl[i];j++){
                        if(ann[i][j] != old_ann && start != -1){
                                stop = j;
                                add_to_ft(aln,start,stop,old_ann,i);
                                old_ann = 0;
                                start = -1;
                        }
                        if(ann[i][j] && !old_ann){
                                start = j+1;
                                old_ann = ann[i][j];
                        }
                }
                fprintf(stderr,"Reading from:%d-%d	start:%d old_ann:%d\n",i,j,start,old_ann);
                if(ann[i][j] != old_ann && start != -1){
                        fprintf(stderr,"ADDDDDDDDDINGGG\n");
                        stop = j;
                        add_to_ft(aln,start,stop,old_ann,i);
                        old_ann = 0;
                        start = -1;
                }
        }
        return aln;
}

void add_to_ft(struct alignment* aln,int start,int stop,int ann,int seq)
{
        char buffer [5];
        struct feature *n = 0;
        struct feature *old_n = 0;
        int i,c;


        fprintf(stderr,"Adding to seq:%d	start:%d	stop:%d	ann:%d\n",seq,start,stop,num_alignments+1 - ann);
        n = tmalloc(sizeof(struct feature));
        n->next = 0;
        n->color = -1;
        n->score  = 0.0f;
        n->type = tmalloc(sizeof(char)*5);
        n->type[0] = 'C';
        n->type[1] = 'M';
        n->type[2] = 'S';
        n->type[3] = 'A';
        n->type[4] = 0;
        n->start = start;
        n->end = stop;
        c = sprintf (buffer,"%d",(num_alignments + 1) - ann);

        n->note = tmalloc(sizeof(char)*(c+1));
        for (i = 0; i < c;i++){
                n->note[i] = buffer[i];
        }
        n->note[c] = 0;

        fprintf(stderr,"Adding:%s	%s:%d-%d col:%d to sequence:%d len:%d\n",n->type,n->note,n->start,n->end,n->color,seq,aln->sl[seq]);

        if((old_n = aln->ft[seq])!= 0){
                while(old_n->next!=0){
                        old_n = old_n->next;
                }
                old_n->next = n;
        }else{
                aln->ft[seq] = n;
        }
        //old_ann = 0;
        //start = -1;
}


struct columns* make_initial_columns(struct alignment* aln,struct columns* columns)
{
        int i,j,n,c;
        columns = tmalloc(sizeof(struct columns));

        columns->len = 0;

        for ( i = 0; i < numseq;i++){
                columns->len += aln->sl[i];
        }
        columns->c = tmalloc(sizeof(int*)*columns->len);
        for (i = 0;i < columns->len;i++){
                columns->c[i] = tmalloc(sizeof(int)*(numseq+1));
                for(j = 0;j < numseq+1;j++){
                        columns->c[i][j] = 0;
                }
        }

        n = 0;
        for (i = 0; i < numseq;i++){
                c = 0;
                for (j = 0;j < aln->len;j++){
                        if(isalpha((int)aln->seq[i][j])){
                                columns->c[n][0] = 1;
                                columns->c[n][i+1] = ((int) aln->seq[i][j]) | (c << 8);
                                n++;
                                c++;
                        }
                }
        }
        return columns;
}

int find_criterion(struct aln_space* aln_space,float level)
{
        int criterion = 0;
        int i;
        float max = -1.0;
        int local_level = 0;

        if(level > num_alignments){
                level = num_alignments;
        }else if(level < 1){
                level = (float)num_alignments * level;
        }
        local_level = (int) level;

        if(!local_level){
                local_level = 1;
        }

        for (i = (1 << num_alignments);i--;){
                if (pop(i) == local_level){
                        if (aln_space->s[i].pcounts > max){
                                max = (int) aln_space->s[i].pcounts;
                                criterion = i;
                        }
                }
        }

        return criterion;
}

float** make_columns(float** ann,struct columns* col,struct alignment** alignments,int criterion)
{
        struct node** hash = 0;
        struct node* old_n = 0;
        struct node* n = 0;

        int** columns = col->c;
        int** cp1 = 0;
        int** cp2 = 0;

        int* c_len = 0;

        int i,j,c;
        int test;
        int hashsize = 0;

        for ( i = 0; i < num_alignments;i++){
                if (hashsize < alignments[i]->len){
                        hashsize = alignments[i]->len;
                }
        }
        hashsize++;

        ann = tmalloc(sizeof(float*)*numseq);
        for (i = 0; i < numseq;i++){
                ann[i] = tmalloc(sizeof(float)*alignments[0]->sl[i]);
                fprintf(stderr,"LEN:%d\n",alignments[0]->sl[i]);
                for (j = 0; j <  alignments[0]->sl[i];j++){
                        ann[i][j] = 0.0f;
                }
        }


        c_len = tmalloc(sizeof(int)*numseq);

        c_len[0] = 0;
        for (i = 1; i< numseq;i++){
                c_len[i] = alignments[0]->sl[i-1] + c_len[i-1];
        }


        hash = tmalloc(sizeof(struct node*) * hashsize);
        for (i = 0;i < hashsize;i++){
                hash[i] = 0;
        }

        for (i = 0; i < numseq;i++){
                cp1 = columns + c_len[i];
                for (j = i + 1; j < numseq;j++){
                        cp2 = columns + c_len[j];
                        for (c = 0; c < num_alignments;c++){
                                hash = feed_hash(hash,alignments[c],i,j,c,0); // - no gaps considered....
                        }
                        for (c = 0;c < hashsize;c++){
                                n = hash[c];
                                test = 0;
                                while (n){

                                        //fprintf(stderr,"%d	aligned to %d	in %d cases.\n",c,n->pos,pop(n->group));
                                        if((n->group & criterion) == criterion && !test){
                                                //if(pop(n->group) == num_alignments){
                                                if(ann[i][c]  < pop(n->group) ){
                                                        ann[i][c] = pop(n->group);
                                                }
                                                if(ann[j][n->pos] < pop(n->group)){
                                                        ann[j][n->pos]  = pop(n->group);
                                                }
                                                //}
                                                test = 1;

                                                //ann[i][c] += pop(n->group);
                                                //ann[j][n->pos] += pop(n->group);
                                                if (cp1[c][0] >= 0){
                                                        if(cp2[n->pos][0] >= 0){
                                                                cp1[c][0] += cp2[n->pos][0];
                                                                cp1[c][j+1] = cp2[n->pos][j+1];
                                                                cp2[n->pos][0] = -1;
                                                                cp2[n->pos][2] = c+c_len[i];
                                                                //	printf("++\n");
                                                        }else{
                                                                cp1[c][0] += columns[cp2[n->pos][2]][0];
                                                                cp1[c][j+1] = columns[cp2[n->pos][2]][j+1];
                                                                columns[cp2[n->pos][2]][0] = -1;
                                                                columns[cp2[n->pos][2]][2] = c+c_len[i];
                                                                //	printf("+-\n");
                                                        }
                                                }else{
                                                        if(cp2[n->pos][0] >= 0){
                                                                //	printf("-+\n");
                                                                columns[cp1[c][2]][0] = cp2[n->pos][0];
                                                                columns[cp1[c][2]][j+1] = cp2[n->pos][j+1];
                                                                cp2[n->pos][0] = -1;
                                                                cp2[n->pos][2] = cp1[c][2];
                                                        }//else{ /// all residues are already aligned correctly......
                                                        //	if(cp1[c][1] != cp2[n->pos][1]){
                                                        //		fprintf(stderr,"--	%d-%d\n",cp1[c][1],cp2[n->pos][1]);
                                                        //	}
                                                        //	columns[cp1[c][1]][j] = columns[cp2[n->pos][1]][j];
                                                        //	columns[cp2[n->pos][1]][0] = -1;
                                                        //	columns[cp2[n->pos][1]][1] = cp1[c][1];
                                                        //}
                                                }

                                        }
                                        //print_columns(columns);
                                        old_n = n;
                                        n = n->next;
                                        free(old_n);
                                }
                                hash[c] = 0;
                        }
                        //printf("\n");
                }
        }
        free(hash);
        free(c_len);





        /*for (i = 0; i < numseq;i++){
          for (j = 0; j <  alignments[0]->sl[i];j++){
          ann[i][j] /= ((numseq-1)*(num_alignments));
          }
          }
        */

        return ann;
}





void pretty_bubble_sort(int** array, int array_size)
{
        int i,j,c;
        int decision = 0;
        int* tmp = 0;
        tmp = tmalloc(sizeof(int)*numseq);

        for (i = 0; i < array_size;i++){
                for(j = 0; j < array_size-1;j++){
                        decision = 1;
                        for (c = 0; c < numseq;c++){
                                if(array[j][c] && array[j+1][c]){
                                        decision = 0;
                                        break;
                                }
                        }
                        if(decision){
                                decision = 0;
                                for (c = 0; c < numseq;c++){
                                        if(array[j][c] && !array[j+1][c]){
                                                decision = 1;
                                                break;
                                        }
                                }
                        }
                        if(decision){//switch...
                                for (c = 0; c < numseq;c++){
                                        tmp[c] = array[j+1][c];
                                }

                                for (c = 0; c < numseq;c++){
                                        array[j+1][c] = array[j][c];
                                }

                                for (c = 0; c < numseq;c++){
                                        array[j][c] = tmp[c];
                                }
                        }

                }
        }
        free(tmp);
}

void quickSort_basic(int** c, int array_size)
{
        q_sort_basic(c, 0, array_size - 1);
}



void q_sort_basic(int** c, int left, int right)
{
        int* pivot;
        int l_hold, r_hold;
        register int i;
        l_hold = left;
        r_hold = right;
        pivot  = tmalloc(sizeof(int)*(numseq+1));
        for (i = 0;i < numseq+1;i++){
                pivot[i] = c[left][i];
        }

        while (left < right){
                while ((c[right][0] <= pivot[0])  && (left < right)){
                        right--;
                }
                if (left != right){
                        for (i = 0;i < numseq+1;i++){
                                c[left][i] = c[right][i];
                        }
                        left++;
                }
                while ((c[left][0] >= pivot[0]) && (left < right)){
                        left++;
                }
                if (left != right){
                        for (i = 0;i < numseq+1;i++){
                                c[right][i] = c[left][i];
                        }
                        right--;
                }
        }

        for (i = 0; i < numseq+1;i++){
                c[left][i] = pivot[i];
        }

        free(pivot);

        i = left;
        left = l_hold;
        right = r_hold;
        if (left < i){
                q_sort_basic(c, left, i-1);
        }
        if (right > i){
                q_sort_basic(c, i+1, right);
        }
}


void quickSort_columns(int** c, int array_size)
{
        q_sort_columns(c, 0, array_size - 1);
}

int a_gt_b(int* a, int *b)
{
        int i;
        int decision = 1;
        for (i = 0;i < numseq;i++){
                if (a[i] && b[i]){
                        if ((a[i]>>8) < (b[i]>>8)){
                                decision = 0;
                                break;
                        }
                }
        }
        return decision;
}

int b_gt_a(int* a, int *b)
{
        int i;
        int decision = 1;
        for (i = 0;i < numseq;i++){
                if (a[i] && b[i]){
                        if ((b[i]>>8) < (a[i]>>8)){
                                decision = 0;
                                break;
                        }
                }
        }
        return decision;
}


void q_sort_columns(int** c, int left, int right)
{
        int* pivot;
        int l_hold, r_hold;
        register int i;
        l_hold = left;
        r_hold = right;
        pivot  = tmalloc(sizeof(int)*numseq);
        for (i = 0;i < numseq;i++){
                pivot[i] = c[left][i];
        }

        while (left < right){
                while ((a_gt_b(c[right],pivot)) && (left < right)){
                        right--;
                }
                if (left != right){
                        for (i = 0;i < numseq;i++){
                                c[left][i] = c[right][i];
                        }
                        left++;
                }
                while ((b_gt_a(c[left],pivot)) && (left < right)){
                        left++;
                }
                if (left != right){
                        for (i = 0;i < numseq;i++){
                                c[right][i] = c[left][i];
                        }
                        right--;
                }
        }

        for (i = 0; i < numseq;i++){
                c[left][i] = pivot[i];
        }

        free(pivot);

        i = left;
        left = l_hold;
        right = r_hold;
        if (left < i){
                q_sort_columns(c, left, i-1);
        }
        if (right > i){
                q_sort_columns(c, i+1, right);
        }
}

int test_order(int** c,int len)
{
        int* test = 0;
        int i,j;
        test = tmalloc(sizeof(int)*numseq);
        for (i = 0 ; i < numseq;i++){
                test[i] = -1;
        }
        for (i = 0; i < len;i++){
                for (j = 0; j < numseq;j++){
                        if(c[i][j]){
                                if(test[j]  > (c[i][j] >> 8)){
                                        free(test);
                                        return 1;
                                }
                                test[j] = c[i][j] >> 8;
                        }
                }
        }
        free(test);
        return 0;
}

struct alignment* inclusion_sort(int** numbers, int array_size,struct alignment* aln,float cutoff,int print_spacer)
{
        int ** c = 0;
        int* spacer = 0;
        int num_spacer = 0;
        int i,j,f,g;
        int aln_len = 0;

        int num_seq_cutoff = numseq * cutoff;
        if(num_seq_cutoff > numseq){
                num_seq_cutoff = numseq-1;
        }
        fprintf(stderr,"CUTOFF:%d\n",num_seq_cutoff);
        //num_seq_cutoff--;

        while(!aln_len){
                for (i = 0;i < array_size;i++){
                        if (numbers[i][0] > num_seq_cutoff){
                                aln_len++;
                        }
                }
                if(!aln_len){
                        num_seq_cutoff--;
                }
        }

        c = tmalloc(sizeof(int*)*aln_len);

        spacer = tmalloc(sizeof(int)*aln_len);
        for(i = 0; i < aln_len;i++){
                spacer[i] = 0;
        }

        f = 0;
        for (i = 0; i <array_size;i++){
                if (numbers[i][0] > num_seq_cutoff){
                        c[f] = tmalloc(sizeof(int)*numseq);
                        for (j = 1;j < numseq+1;j++){
                                c[f][j-1] = numbers[i][j];
                        }
                        f++;
                }
        }




        i = 1;
        while(test_order(c,aln_len)){
                fprintf(stderr,"Sorting alignment columns... %d\n",i);
                quickSort_columns(c,aln_len);
                i++;
        }

        if(print_spacer){
                for (i = 1; i < aln_len;i++){
                        for (j = 0; j < numseq;j++){
                                if(c[i-1][j] && c[i][j]){
                                        if(((c[i-1][j] >> 8) + 1) != (c[i][j] >> 8)){
                                                spacer[i] = 1;
                                                num_spacer++;
                                                break;//j = numseq;
                                        }
                                }
                        }
                        /*for (j = 0; j < numseq;j++){
                          fprintf(stderr,"%c",(char)c[i][j] & 0xff);
                          }
                          fprintf(stderr,"\n");*/
                }
        }
        /*fprintf(stderr,"\n");
          fprintf(stderr,"\n");
          for (i = 0; i < aln_len;i++){
          fprintf(stderr,"%d ",spacer[i]);
          }
          fprintf(stderr,"\n");*/

        aln->len = aln_len;

        for (i = 0; i < numseq;i++){
                free(aln->seq[i]);
                free(aln->s[i]);
                aln->seq[i] = tmalloc(sizeof(char)*(aln_len+1+(num_spacer*10)));
                aln->s[i] = tmalloc(sizeof(int)*(aln_len+1+(num_spacer*10)));
        }
        f = 0;
        for (i = 0; i < aln_len;i++){
                if(spacer[i]){
                        for (g = 0; g < 10;g++){
                                for (j = 0; j < numseq;j++){
                                        aln->seq[j][f] = '-';
                                        aln->s[j][f] = -1;
                                        //		fprintf(stderr,"-");
                                }
                                //	fprintf(stderr,"\n");
                                f++;
                        }
                }

                for (j = 0; j < numseq;j++){
                        if ((char)(c[i][j]&0x000000ff)){
                                aln->s[j][f]= c[i][j] >>8;
                                aln->seq[j][f] = (char)(c[i][j]&0x000000ff);
                        }else{
                                aln->s[j][f] = -1;
                                aln->seq[j][f] = '-';
                        }
                        //fprintf(stderr,"%c",aln->seq[j][f]);
                }
                //fprintf(stderr,"\n");
                f++;
                free(c[i]);
        }
        free(c);

        for (i = 0; i < numseq;i++){
                aln->seq[i][aln_len+(num_spacer*10)] = 0;
        }
        aln->len +=(num_spacer*10);
        free(spacer);
        return aln;
}

void  insert_column(int* column,struct column_node* old_n,struct column_node* n)
{
        int i;
        int decision;
        struct column_node* tmp;
        if (n == NULL){
                n = tmalloc(sizeof(struct column_node));
                n->column = tmalloc( sizeof(int)* numseq);
                for (i = 0;i < numseq;i++){
                        n->column[i] = column[i];
                }
                n->next = 0;
                old_n->next = n;
        }else if(n->column){
                decision = 0;
                for (i = 0;i < numseq;i++){
                        if (n->column[i] && column[i]){
                                if ((n->column[i]>>8) > (column[i]>>8)){
                                        decision = 1;
                                }
                        }
                }
                if (decision){
                        tmp = tmalloc(sizeof(struct column_node));
                        tmp->column = tmalloc( sizeof(int)* numseq);
                        for (i = 0;i < numseq;i++){
                                tmp->column[i] = column[i];
                        }
                        tmp->next = n;
                        old_n->next = tmp;
                        return;
                }else{
                        insert_column(column,n,n->next);
                }
        }else{
                insert_column(column,n,n->next);
        }
        return;
}
