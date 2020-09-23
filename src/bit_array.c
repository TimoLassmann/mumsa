/*
 	bit_array.c
	
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

#include "bit_array.h"


void bit_array_set(int* a,int i)
{
	a[i >> 5] |= (1 << (i & 0x1F));
}

void bit_array_clr(int* a, int i)
{
	a[i >> 5] &= ~(1 << (i & 0x1F));
}

int bit_array_test(int* a, int i)
{
	return (a[i >> 5] & (1<< (i & 0x1F)));
}

