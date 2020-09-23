/*
 	print_reliable_alignment.h
	
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


struct block{
	int start;
	int end; 
	int* sseq;
};


struct column_node{
	struct column_node* next;
	int* column;
};

struct columns* make_initial_columns(struct alignment* aln,struct columns* columns);
void  insert_column(int* column,struct column_node* old_n,struct column_node* n);
float** make_columns(float** ann, struct columns* col,struct alignment** alignments,int criterion);
struct alignment* inclusion_sort(int** numbers, int array_size,struct alignment* aln,float cutoff,int print_spacer);
int find_criterion(struct aln_space* aln_space,float level);
struct alignment* add_annotation(struct alignment* aln,float** ann);

struct alignment* make_alignment(struct columns* columns,struct alignment* alignment);


void make_alignment_from_all_poars_at_once(struct alignment** alignments,int criterion);
struct node*** feed_alignment_into_hash(struct node*** hash, struct alignment* aln,int* c_len,int aln_n);
struct node * insert_aln_pair( struct node* n,unsigned int pos, unsigned int aln);


struct aln_node* insert_pair_into_alignment(struct aln_node* n,int seq_a,int seq_b,int a,int b,int nl);
void print_aln_node(struct aln_node* n);
void free_aln_node(struct aln_node* n);

void pretty_bubble_sort(int** array, const int array_size);

void quickSort_basic(int** c, int array_size);
void q_sort_basic(int** c, int left, int right);


int check_block_swapability(int* a,int*b,int size);
struct block* block_alloc(struct block* b, int size);
void block_free(struct block* b);



void add_to_ft(struct alignment* aln,int start,int stop,int ann,int seq);
void quickSort_columns(int** c, int array_size);
inline int a_gt_b(int* a, int *b);
int b_gt_a(int* a, int *b);


int fncompare (const void * elem1, const void * elem2 );
int compare_a_b(int* a, int *b);

void q_sort_columns(int** c, int left, int right);

int test_order(int** c,int len);



