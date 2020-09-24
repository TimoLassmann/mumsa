/*
  global.h

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

#ifndef GLOBAL_H
#define GLOBAL_H



#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>


/* #include "tldevel.h" */

struct parameters{
        char** infile;
        int score_mode;
        int num_infiles;
        int help_flag;

};

#define MUMSA_SCORE_REF 1
#define MUMSA_SCORE_REF_TEST 2
#define MUMSA_SCORE_TEST 3



#define ALN_STATUS_UNALIGNED 1   /* no gaps sequences may or may not have equal lengths  */
#define ALN_STATUS_ALIGNED 2   /* sequences have equal lengths and may or may not contain gaps*/
#define ALN_STATUS_UNKNOWN 3     /* sequences have un-equal length and contain gaps  */

#define BUFFER_LEN 128
#endif
