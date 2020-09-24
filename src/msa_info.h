/*
  msa_info.h

  Released under GPL - see the 'COPYING' file

  Copyright (C) 2020 Timo Lassmann <timolassmann@gmail.com>

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
#ifndef MSA_INFO_H
#define MSA_INFO_H

#ifdef MSA_INFO_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif


struct msa_info{
        int** s;
        int* adler_val;
        double pairs;
        int aln_len;
        int num_seq;

};

EXTERN int alloc_msa_info(struct msa_info** info, int numseq, int aln_len);
EXTERN void free_msa_info(struct msa_info* msai);

#undef MSA_INFO_IMPORT
#undef EXTERN
#endif
