/*
  cmsa.h

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


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/wait.h>

#ifdef MEMORY
#define tmalloc malloc
#endif

#define SEEK_START 0
#define SEEK_END 2

#define INFTY 936870912


double ran1(int *idum);


extern unsigned int numseq;
extern unsigned int numprofiles;
extern unsigned int num_alignments;

#ifndef MEMORY
void* tmalloc(int size);
#endif

struct parameters{
        int* order;
        char** infile;
        char* outfile;
        char* format;
        char* sort;
        char* ref_flag;

        float best_flag;
        float idal;
        float x1;
        float x2;
        float x3;
        float bootcutoff;
        float print_column_cutoff;

        int print_spacer;
        int gap_flag;
        int quiet_flag;
        int print_alignment_flag;
        int server_flag;
        int bootstrap;
        int set_diff;
};


struct columns{
        int** c;
        int len;
};

struct sets{
        float score_4_input_aln;
        float aln_sim;
        float id;
        float pcounts;
};

struct aln_space{
        struct sets* s;
        float** o;
        float** sim;
        float* scores;
        float* support;

        float* ref_10;
        float* ref_11;
        float* ref_01;

        float diff;


};

struct node{
        struct node* next;
        unsigned int pos;
        unsigned int group;
};


struct alignment{
        struct feature** ft;
        struct sequence_info** si;
        unsigned int** sip;
        unsigned int* nsip;
        unsigned int* sl;
        unsigned int* lsn;
        unsigned int* adler;
        int** s;
        char**seq;
        char** sn;
        float id;
        float av_len;
        float al;
        float score;
        unsigned int len;
        unsigned int pairs;
        unsigned int numseq;
};

struct feature{
        struct feature *next;
        char* type;
        char* note;
        float score;
        int start;
        int end;
        int color;

};

int* make_tree_order(int* order,float** sim,struct alignment* aln);
int* upgma(float **dm,int* tree);

void quickalnSort(struct alignment* aln, int array_size);
void q_sort(struct alignment* aln, int left, int right);

void is_alignment(struct alignment** alignments,struct parameters* param);
void sanity_check(struct alignment** alignments,struct parameters* param);

int* get_input_order(int* input_order,struct alignment* aln);

struct alignment* calculate_adler(struct alignment* alignment);


struct aln_space* score_sets(struct aln_space* aln_space, struct alignment** alignments,struct parameters* param);
struct alignment* alistat(struct alignment* aln);

struct alignment** get_reliable_blocks(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param);
//void print_reliable_blocks(struct alignment** alignments,struct aln_space* aln_space,int hashsize,struct parameters* param);

struct alignment* change_character_numbers_and_seq_len(struct alignment* aln);

struct alignment* clean_alignment(struct alignment* aln,int* order);
void output(struct alignment* aln,struct parameters* param);
struct alignment** detect_and_read_alignments(struct alignment** alignments,struct parameters* param);


int byg_detect(int* text,int n);

int check_identity(char* n,char*m);

int byg_count(char* pattern,char*text);
int byg_start(char* pattern,char*text);
int byg_end(char* pattern,char*text);


struct parameters* interface(struct parameters* param,int argc,char **argv);


void free_ft(struct feature* n);
void free_param(struct parameters* param);
struct alignment* aln_alloc(struct alignment* aln);
void free_aln(struct alignment* aln);
struct aln_space* aln_space_alloc(struct aln_space* aln_space);
void free_aln_space(struct aln_space* aln_space);
void free_columns(struct columns* columns);

//int getopt(int argc, char * const argv[],const char *optstring);
struct node * insert_pair(unsigned int pos,unsigned int aln,struct node* n);

struct node** feed_hash(struct node** hash,struct alignment* ap,unsigned int a, unsigned int b,unsigned int aln,int gap_flag);

//struct alignment* read_aln_matrix(struct alignment* alignment,char* infile);
//struct alignment* read_aln_matrix(struct alignment* ap, char* string);

//int** read_columns_from_alignment(int** columns,char* infile,int total_seq_len);

//char* get_alignment_into_string(char* string,char* infile);

//float* fill_sets(struct node** hash,float* sets,unsigned int hashsize,float test);
struct aln_space* make_sets(struct aln_space*,struct alignment** alignments,int gap_flag);

void r_statistics_output(char ** name,float* scores,float** overlap,int num);

//float** calculate_overlap(struct aln_space* aln_space,float** overlap);
struct aln_space* calculate_overlap(struct aln_space* aln_space,struct alignment** alignments,char* ref);
int pop(unsigned int x);

//void r_statistics_output(char ** name,float* scores,float** overlap,int num);
void print_alignment_scores(struct aln_space* aln_space, struct parameters* param,struct alignment** alignments);

//struct aln_space* calculate_diff(struct aln_space* aln_space);

struct aln_space* calculate_diff(struct aln_space* aln_space,struct alignment** alignments,int cutoff);

float* find_overlap_to_ref(float* overlap_to_ref,float* sets);

struct aln_space* find_best_alignment(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param);

//float* score_alignment(struct alignment** alignments,float* counts,int* ann,float* scores);

struct aln_space* eval_space(struct alignment** alignments,struct aln_space* aln_space,struct parameters* param);
struct aln_space* bootstrap(struct aln_space* aln_space, struct alignment** alignments,struct parameters* param);
//char** get_seq_names(char** seq_names,char* infile);
