/*
 	upgma.c
	
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

int* upgma(float **dm,int* tree)
{
	int i,j,t;
	int *as = 0;
	float max;
	int node_a = 0;
	int node_b = 0;
	int cnode = numseq;

	as = tmalloc(sizeof(int)*numseq);
	for (i = numseq; i--;){
		as[i] = i+1;
	}

	
	t = 0;
	while (cnode != numprofiles){
		max = -INFTY;
		for (i = 0;i < numseq-1; i++){
			if (as[i]){
			for ( j = i + 1;j < numseq;j++){
				if (as[j]){
				if (dm[i][j] > max){
					max = dm[i][j];
					node_a = i;
					node_b = j;
				}
				}
			}
			}
		}
		
		tree[t] = as[node_a]-1;
		tree[t+1] = as[node_b]-1;
		tree[t+2] = cnode;
		t += 3;	
		
		/*deactivate  sequences to be joined*/
		as[node_a] = cnode+1;
		as[node_b] = 0;
		cnode++;    
		
		/*calculate new distances*/
		for (j = numseq;j--;){
			if (j != node_b){
				dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])/2;
			}
		}
		dm[node_a][node_a] = 0.0f;
		for (j = numseq;j--;){
			dm[j][node_a] = dm[node_a][j];
			dm[j][node_b] = 0.0f;
			dm[node_b][j] = 0.0f;
		}		
	}
	free(as);
	return tree;
}

