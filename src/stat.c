/*
  stat.c

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

struct alignment* alistat(struct alignment* aln)
{
        float id = 0.0f;
        float len = 0.0f;
        float poar = 0.0f;
        int i,j,c;

        aln->id = 0.0;
        aln->al = 0.0;
        aln->av_len = 0.0;
        for (i = 0; i < aln->numseq;i++){
                for(j = 0; j < aln->len;j++){
                        if (isalpha(aln->seq[i][j]) && !iscntrl(aln->seq[i][j]) ){
                                aln->av_len++;
                        }
                }
        }
        aln->av_len /= aln->numseq;

        for (i = 0; i < aln->numseq-1;i++){
                for (j = i+1;j < aln->numseq;j++){
                        if(aln->sl[i] > aln->sl[j]){
                                aln->al += aln->sl[j];
                        }else{
                                aln->al += aln->sl[i];
                        }
                }
        }

        for (i = 0; i< aln->numseq;i++){
                for (j = i+1;j < aln->numseq;j++){
                        id = 0.0f;
                        len = 0.0f;
                        for (c = 0;c < aln->len;c++){
                                if (isalpha(aln->seq[i][c])){
                                        if (isalpha(aln->seq[j][c])){
                                                if (aln->seq[i][c] == aln->seq[j][c]){
                                                        id++;
                                                }
                                                len++;
                                                poar++;
                                        }
                                }
                        }

                        if (!len){
                                len = 1;
                        }

                        aln->id += id / len;
                }
        }
        aln->id /= (aln->numseq*(aln->numseq-1))/2;
        aln->al = poar/aln->al;
        return aln;
}
