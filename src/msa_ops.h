/*
  msa_ops.h

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

#ifndef MSA_OPS_H
#define MSA_OPS_H

#ifdef MSA_OPS_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

struct msa_info;
struct msa;
struct mumsa_data;

EXTERN int read_msa_into_msai(struct msa* msa,struct msa_info** msai);
EXTERN int msai_char_to_pos(struct msa_info* msai);
EXTERN int alistat(struct mumsa_data* m);


EXTERN int sanity_check_input(struct mumsa_data* md);
EXTERN int sort_msa_by_seqname(struct msa* msa);

#undef MSA_OPS_IMPORT
#undef EXTERN
#endif
