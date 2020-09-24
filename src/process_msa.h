/*
  process_msa.h

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

#ifndef PROCESS_MSA_H
#define PROCESS_MSA_H

#ifdef PROCESS_MSA_IMPORT
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int calc_sim_pairs(struct mumsa_data* m);
EXTERN int calc_overlap(struct mumsa_data* m);

#undef PROCESS_MSA_IMPORT
#undef EXTERN
#endif
