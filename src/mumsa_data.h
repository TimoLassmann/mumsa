/*
  mumsa_data.h

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

#ifndef MUMSA_DATA_H
#define MUMSA_DATA_H

#include "msa.h"
#include "msa_info.h"

struct sets{
        double score_4_input_aln;
        double aln_sim;
        double id;
        double pcounts;
};


struct mumsa_data{
        struct msa** msa;
        struct msa_info** msai;
        struct sets* s;
        double** overlap;
        double** sim;

        double* poar;
        double* id;
        double* avg_len;
        double* al;
        int num_aln;
        int num_seq;
};


#ifdef MUMSA_DATA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int alloc_mumsa_data(struct mumsa_data** mdat, int num_alignment);
EXTERN int free_mumsa_data(struct mumsa_data* m);


#undef MUMSA_DATA_IMPORT
#undef EXTERN
#endif
