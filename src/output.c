/*
  output.c

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
#include "output.h"
#include <string.h>

struct alignment* clean_alignment(struct alignment* aln,int* order)
{
        int i,j;
        int test = 0;

        //cos of sorting.....

        if(!order){
                for (i = 0; i <  aln->numseq;i++){
                        aln->nsip[i] = i;
                }
        }else{
                for (i = 0; i < aln->numseq;i++){
                        test = order[i];
                        for (j = 0;j < aln->numseq;j++){
                                if(aln->adler[j] == test){
                                        aln->nsip[i] = j;
                                        break;
                                }
                        }
                }
        }

        //for (i = 0; i <  aln->numseq;i++){
        //	fprintf(stderr,"%d	%d	in_aln:%d	order:%d\n",i,aln->nsip[i],aln->adler[i],order[i]);
        //}

        for (i = 0; i < aln->numseq;i++){
                //fprintf(stderr,"Final len:%d	%d	apparent len:%d	real:%d\n",i,aln->sl[i],aln->len,(int)strlen(aln->seq[i]));
                aln->sl[i] =  aln->len;
                for (j = 0; j <= aln->len;j++){
                        aln->s[i][j] = 0;
                }
        }
        return aln;
}

void output(struct alignment* aln,struct parameters* param)
{

        if(param->outfile){
                if ((fout = fopen(param->outfile, "w")) == NULL){
                        fprintf(stderr,"can't open output\n");
                        exit(0);
                }
        }else{
                fout = stdout;
        }

        if(!param->format){
                fasta_output(aln,param->outfile);
        }else{
                if (byg_start("msf",param->format) != -1){
                        msf_output(aln,param->outfile);
                }else if (byg_start(param->format,"clustalcluCLUSTALCLUCLUSTALW") != -1){
                        clustal_output(aln,param->outfile);
                }else if (byg_start("macsim",param->format) != -1){
                        macsim_output(aln,param->outfile,param->infile[0]);
                }else{
                        fasta_output(aln,param->outfile);
                }
        }
        if(param->outfile){
                fclose(fout);
        }

}

static void macsim_output(struct alignment* aln,char* outfile,char* infile)
{
        int i,j,f;
        int tmp;
        struct feature *fn = 0;

        fprintf(fout,"<?xml version=\"1.0\"?>\n<!DOCTYPE macsim SYSTEM \"http://www-bio3d-igbmc.u-strasbg.fr/macsim.dtd\">\n<macsim>\n<alignment>\n<aln-name>");
        if(infile){
                fprintf(fout,"%s.kalign",infile);
        }else{
                fprintf(fout,"kalign alignment");
        }
        fprintf(fout,"</aln-name>\n");

        for (i =0;i< numseq;i++){
                //c = aln->sl[i];
                f = aln->nsip[i];

                fprintf(fout,"<sequence seq-type=\"Protein\">\n");
                fprintf(fout,"<seq-name>");
                for (j =0; j < aln->lsn[f];j++){
                        if(!iscntrl((int)aln->sn[f][j])){
                                fprintf(fout,"%c",aln->sn[f][j]);
                        }
                }
                fprintf(fout,"</seq-name>");
                fprintf(fout,"<seq-info>\n");
                fprintf(fout,"<accession>1aab_</accession>\n");
                fprintf(fout,"<nid>1aab_</nid>\n");
                fprintf(fout,"<ec>0.0.0.0</ec>\n");
                fprintf(fout,"<group>0</group>\n");
                if(aln->ft){
                        if(aln->ft[f]){

                                fprintf(fout,"<ftable>\n");
                                fn = aln->ft[f];
                                while(fn){
                                        fprintf(fout,"<fitem><ftype>%s</ftype><fstart>%d</fstart><fstop>%d</fstop><fcolor>%d</fcolor><fscore>%f</fscore><fnote>%s</fnote></fitem>\n",fn->type,fn->start,fn->end,fn->color,fn->score,fn->note);
                                        fn = fn->next;
                                }
                                fprintf(fout,"</ftable>\n</seq-info>\n");
                        }
                }
                fprintf(fout,"<seq-data>\n");

                for (j = 0; j < aln->sl[f];j++){
                        tmp = aln->s[f][j];
                        while (tmp){
                                fprintf(fout,"-");
                                tmp--;
                        }
                        fprintf(fout,"%c",aln->seq[f][j]);
                }
                tmp = aln->s[f][aln->sl[f]];
                while (tmp){
                        fprintf(fout,"-");
                        tmp--;
                }
                fprintf(fout,"\n");
                fprintf(fout,"</seq-data>\n");
                fprintf(fout,"</sequence>\n");
        }
        fprintf(fout,"</alignment>\n");
        fprintf(fout,"</macsim>\n");
        //if(outfile){
        //	fclose(stdout);
        //}
        //free_aln(aln);
}

static void msf_output(struct alignment* aln,char* outfile)
{
        char** linear_seq = 0;
        int i,j,c,f,g;
        int max = 0;
        int aln_len = 0;
        int tmp;

        linear_seq = tmalloc(sizeof(char*)*numseq);

        aln_len = 0;
        for (j = 0; j <= aln->sl[0];j++){
                aln_len+= aln->s[0][j];
        }
        aln_len += aln->sl[0];

        for (i = 0; i < numseq;i++){
                linear_seq[i] = tmalloc(sizeof(char)*(aln_len+1));

                c = 0;
                for (j = 0; j < aln->sl[i];j++){
                        tmp = aln->s[i][j];
                        while (tmp){
                                linear_seq[i][c] ='-';
                                c++;
                                tmp--;
                        }
                        linear_seq[i][c] = aln->seq[i][j];
                        c++;
                }

                tmp =aln->s[i][aln->sl[i]];
                while (tmp){
                        linear_seq[i][c] ='-';
                        c++;
                        tmp--;
                }
                linear_seq[i][c] = 0;
        }

        //if(outfile){
        //	if ((stdout = fopen(outfile, "w")) == NULL){
        //		fprintf(stderr,"can't open output\n");
        //		exit(0);
        //	}
        //}

        fprintf(fout,"PileUp\n\n\n\n   MSF:   %d  Type: P    Check:  7038   ..\n\n",aln_len);

        for (j = 0; j< numseq;j++){
                if( aln->lsn[j] > max){
                        max = aln->lsn[j];
                }
        }

        for (i = 0; i< numseq;i++){
                f = aln->nsip[i];
                fprintf(fout," Name: ");
                for (c = 0; c < aln->lsn[f];c++){
                        if(!iscntrl((int)aln->sn[f][c])){
                                fprintf(fout,"%c",aln->sn[f][c]);
                        }
                }
                while(c < max+3){
                        fprintf(fout," ");
                        c++;
                }
                fprintf(fout,"Len:   ");
                fprintf(fout,"%d",aln_len);
                fprintf(fout,"  Check:  2349  Weight:  1.00\n");

        }
        fprintf(fout,"\n\n//\n\n");

        for (i = 0; i+60 < aln_len;i +=60){
                for (j = 0; j< numseq;j++){
                        f = aln->nsip[j];
                        for (c = 0; c < aln->lsn[f];c++){
                                if(!iscntrl((int)aln->sn[f][c])){
                                        fprintf(fout,"%c",aln->sn[f][c]);
                                }
                        }
                        while(c < max+3){
                                fprintf(fout," ");
                                c++;
                        }
                        g = 1;
                        for (c = 0; c < 60;c++){
                                fprintf(fout,"%c",linear_seq[f][c+i]);
                                if (g == 10){
                                        fprintf(fout," ");
                                        g = 0;
                                }
                                g++;
                        }
                        fprintf(fout,"\n");

                }
                fprintf(fout,"\n\n");
        }
        for (j = 0; j< numseq;j++){
                f = aln->nsip[j];
                for (c = 0; c< aln->lsn[f];c++){
                        if(!iscntrl((int)aln->sn[f][c])){
                                fprintf(fout,"%c",aln->sn[f][c]);
                        }
                }
                while(c < max+3){
                        fprintf(fout," ");
                        c++;
                }
                g = 1;
                for (c = i; c< aln_len;c++){
                        fprintf(fout,"%c",linear_seq[f][c]);
                        if (g == 10){
                                fprintf(fout," ");
                                g = 0;
                        }
                        g++;
                }
                fprintf(fout,"\n");

        }
        fprintf(fout,"\n\n");
        //if(outfile){
        //	fclose(stdout);
        //}
        for (i =0;i< numseq;i++){
                free(linear_seq[i]);
        }
        free(linear_seq);
        //free_aln(aln);
}

static void clustal_output(struct alignment* aln,char* outfile)
{
        int i,j,c,f;
        int tmp;
        int aln_len = 0;
        char** linear_seq = 0;

        linear_seq = tmalloc(sizeof(char*)*numseq);



        aln_len = 0;
        for (j = 0; j <= aln->sl[0];j++){
                aln_len+= aln->s[0][j];
        }
        aln_len += aln->sl[0];

        for (i =0;i< numseq;i++){
                linear_seq[i] = tmalloc(sizeof(char)*(aln_len+1));

                c = 0;
                for (j = 0; j < aln->sl[i];j++){
                        tmp = aln->s[i][j];
                        while (tmp){
                                linear_seq[i][c] ='-';
                                c++;
                                tmp--;
                        }
                        linear_seq[i][c] = aln->seq[i][j];
                        c++;
                }

                tmp =aln->s[i][aln->sl[i]];
                while (tmp){
                        linear_seq[i][c] ='-';
                        c++;
                        tmp--;
                }
                linear_seq[i][c] = 0;
        }


        //if(outfile){
        //	if ((stdout = fopen(outfile, "w")) == NULL){
        //		fprintf(stderr,"can't open output\n");
        //		exit(0);
        //	}
        //}

        fprintf(fout,"CLUSTAL W (1.83) multiple sequence alignment\n\n\n");

        for (i = 0; i+60 < aln_len;i +=60){
                for (j = 0; j< numseq;j++){
                        f = aln->nsip[j];
                        for (c = 0; c < aln->lsn[f];c++){
                                if(!iscntrl((int)aln->sn[f][c])){
                                        fprintf(fout,"%c",aln->sn[f][c]);
                                }
                        }
                        while(c < 18){
                                fprintf(fout," ");
                                c++;
                        }

                        for (c = 0; c< 60;c++){
                                fprintf(fout,"%c",linear_seq[f][c+i]);
                        }
                        fprintf(fout,"\n");
                }
                fprintf(fout,"\n\n");
        }
        for (j = 0; j< numseq;j++){
                f = aln->nsip[j];
                for (c = 0; c< aln->lsn[f];c++){
                        if(!iscntrl((int)aln->sn[f][c])){
                                fprintf(fout,"%c",aln->sn[f][c]);
                        }
                }
                while(c < 18){
                        fprintf(fout," ");
                        c++;
                }

                for (c = i; c< aln_len;c++){
                        fprintf(fout,"%c",linear_seq[f][c]);
                }
                fprintf(fout,"\n");
        }
        fprintf(fout,"\n\n");
        //if(outfile){
        //	fclose(stdout);
        //}
        for (i =0;i< numseq;i++){
                free(linear_seq[i]);
        }
        free(linear_seq);
        //free_aln(aln);
}

static void fasta_output(struct alignment* aln,char* outfile)
{
        int i,j,c,f;
        int tmp;
        //if(outfile){
        //	if ((stdout = fopen(outfile, "w")) == NULL){
        //		fprintf(stderr,"can't open output\n");
        //		exit(0);
        //	}
        //}
        for (i = 0; i < numseq;i++){
                f = aln->nsip[i];
                fprintf(fout,">%s\n",aln->sn[f]);
                c = 0;
                for (j = 0; j < aln->sl[f];j++){
                        tmp = aln->s[f][j];
                        while (tmp){
                                fprintf(fout,"-");
                                c++;
                                if(c == 60){
                                        fprintf(fout,"\n");
                                        c = 0;
                                }
                                tmp--;
                        }
                        fprintf(fout,"%c",aln->seq[f][j]);
                        c++;
                        if(c == 60){
                                fprintf(fout,"\n");
                                c = 0;
                        }
                }
                tmp =aln->s[f][aln->sl[f]];
                while (tmp){
                        fprintf(fout,"-");
                        c++;
                        if(c == 60){
                                fprintf(fout,"\n");
                                c = 0;
                        }
                        tmp--;
                }
                fprintf(fout,"\n");
        }
        //if(outfile){
        //	fclose(stdout);
        //}
        //free_aln(aln);
}

void print_alignment_scores(struct aln_space* aln_space, struct parameters* param,struct alignment** alignments)
{

        float diff = aln_space->diff;
        char** aln_name = param->infile;
        float* scores = aln_space->scores;
        float* support = aln_space->support;

        int i,j,c;
        float max = -1.0;
        int max_len = 0;
        int* len = 0;
        char **p = 0;

        if(param->outfile){
                if ((fout = fopen(param->outfile, "w")) == NULL){
                        fprintf(stderr,"can't open output\n");
                        exit(0);
                }
        }else{
                fout = stdout;
        }


        p = tmalloc(sizeof(char*)*num_alignments);
        len = tmalloc(sizeof(int)*num_alignments);

        for (i = 0;i < num_alignments;i++){
                p[i] = aln_name[i];
                j =0;
                c = 0;
                while(aln_name[i][j]){
                        if(aln_name[i][j] == 47){//47 = '/'
                                c = j;
                        }
                        j++;
                }
                if(c){
                        p[i] += c+1;
                }
                len[i] = j-c+1;
                if(len[i] > max_len){
                        max_len = len[i];
                }
                //	printf("%s\n",name[i]);
        }

        if (!param->server_flag){
                if(diff != -1){
                        fprintf(fout,"Alignment difficulty:%0.2f\n",diff);
                }
                fprintf(fout,"Input Alignment:");
                if(15 < max_len){
                        for (j = 0; j < max_len - 15;j++){
                                fprintf(fout," ");
                        }
                }
                if(support){
                        fprintf(fout,"	Score:	Support:\n");
                }else{
                        fprintf(stdout,"	Score:	%%ID	%%AL\n");
                }
                for (i = 0; i < num_alignments;i++){
                        max = -1.0;
                        c = 0;
                        for (j = 0;j < num_alignments;j++){
                                if (scores[j] > max){
                                        max = scores[j];
                                        c = j;
                                }
                        }
                        fprintf(fout,"%s",p[c]);
                        for (j = 0; j < max_len - len[c];j++){
                                fprintf(fout," ");
                        }
                        if(support){
                                fprintf(fout,"	%f	%d%%\n",scores[c],(int)(support[c]*100));
                        }else{
                                fprintf(fout,"	%0.2f	%0.2f	%0.2f\n",scores[c],alignments[c]->id,alignments[c]->al);
                        }
                        scores[c] = -2;
                }
        }else{
                if(diff != -1){
                        fprintf(fout,"%f	",diff);
                }

                for (i = 0; i < num_alignments;i++){
                        fprintf(fout,"%f	",scores[i]);
                }
                fprintf(fout,"\n");
        }
        free(len);
        free(p);
        if(param->outfile){
                fclose(fout);
        }
}
