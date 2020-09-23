/*
  make_sets.c

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

struct aln_space* make_sets(struct aln_space* aln_space, struct alignment** alignments, int gap_flag)
{
  struct node** hash = 0;
  struct node* old_n = 0;
  struct node* n = 0;

  float* sim = 0;
  float tmp;
  int i,j,c;
  int hashsize = 0;


  for ( i = 0; i < num_alignments;i++){
    alignments[i]->pairs = 0;
    if (hashsize < alignments[i]->len){
      hashsize = alignments[i]->len;
    }
  }
  hashsize++;

  hash = tmalloc(sizeof(struct node*) * hashsize);
  for (i = hashsize;i--;){
    hash[i] = 0;
  }
  for (i = 0;i < (1 << num_alignments);i++){
    aln_space->s[i].pcounts = 0.0f;
  }


  for (i = 0; i < numseq;i++){
    sim = aln_space->sim[i];
    for (j = i +1; j < numseq;j++){
      for (c = 0; c < num_alignments;c++){
        hash = feed_hash(hash,alignments[c],i,j,c,gap_flag);
      }

      for (c = hashsize;c--;){
        n = hash[c];
        while (n){
          aln_space->s[n->group].pcounts++;
          if( pop(n->group) > 1){
            sim[j] += pop(n->group);
          }
          old_n = n;
          n = n->next;
          free(old_n);
        }
        hash[c] = 0;
      }
      if(gap_flag){
        c = (alignments[0]->sl[j] > alignments[0]->sl[i])?alignments[0]->sl[j]:alignments[0]->sl[i];
        sim[j] = sim[j]/(float)(num_alignments*c);
      }else{
        c = (alignments[0]->sl[j] < alignments[0]->sl[i])?alignments[0]->sl[j]:alignments[0]->sl[i];
        sim[j] = sim[j]/(float)(num_alignments*c);
      }
    }
  }
  for (i =0; i < numseq;i++){
    for (j = i + 1; j < numseq;j++){
      aln_space->sim[j][i] = aln_space->sim[i][j];
    }
  }
  for (i =0; i < numseq;i++){
    tmp = 0.0f;
    for (j = 0; j < numseq;j++){
      tmp+= aln_space->sim[i][j];
    }
    aln_space->sim[i][i] = tmp/(float)numseq;
  }

  free(hash);

  return aln_space;
}

struct node** feed_hash(struct node** hash, struct alignment* ap, unsigned int a ,unsigned int b, unsigned int aln, int gap_flag)
{
  unsigned int i;
  unsigned int p = 0;
  unsigned int len = ap->len;
  int* seqa = 0;
  int* seqb = 0;
  seqa = ap->s[a];
  seqb = ap->s[b];
  if (!gap_flag){
    for (i = 0; i < len;i++){
      if (seqa[i] < len && seqb[i] < len){
        hash[seqa[i]] = insert_pair(seqb[i],aln,hash[seqa[i]]);
        p++;
      }
    }
  }else{
    for (i = 0; i < len;i++){
      /*if (seqa[i] < len && seqb[i] < len){
        hash[seqa[i]] = insert(seqb[i],aln,hash[seqa[i]]);
        p++;
      }else if (seqa[i] == len && seqb[i] < len){
        hash[gap_flag] = insert(seqb[i],aln,hash[len]);
        p++;
      }else if (seqb[i] == len && seqa[i] < len){
        hash[seqa[i]] = insert(gap_flag,aln,hash[seqa[i]]);
        p++;
      }*/
      if (seqa[i] < len){
        if (seqb[i] < len){
          hash[seqa[i]] = insert_pair(seqb[i],aln,hash[seqa[i]]);
        }else{
          hash[seqa[i]] = insert_pair(gap_flag,aln,hash[seqa[i]]);
        }
        p++;
      }else if (seqa[i] == len && seqb[i] < len){
        hash[gap_flag] = insert_pair(seqb[i],aln,hash[len]);
        p++;
      }
    }
  }
  ap->pairs += p;
  return hash;
}

/*
float* fill_sets(struct node** hash, float* sets, unsigned int hashsize,float test)
{
  unsigned int i;
  struct node* n;
  for (i = 0;i < hashsize;i++){
    n = hash[i];
    while (n){
      test += pop( n->group);
      sets[n->group]++;
      n = n->next;
    }
  }
  return sets;
}
*/

struct node * insert_pair(unsigned int pos, unsigned int aln, struct node* n)
{
        if (n == NULL){
                n = tmalloc(sizeof(struct node));
                n->pos = pos;
                n->next = 0;
                aln = 1 << aln;
                n->group = 0;
                n->group |= aln;
        }else if(n->pos == pos){
                aln = 1 << aln;
                n->group |= aln;
                return n;
        }else{
                n->next = insert_pair(pos,aln,n->next);
        }
        return n;
}
