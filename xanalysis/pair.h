/*
Copyright 2005-2020 Kendall F. Morris

This file is part of the Xanalysis software suite.

    The Xanalysis software suite is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either
    version 3 of the License, or (at your option) any later version.

    The suite is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the suite.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef PAIRLOOP_H
#define PAIRLOOP_H

#include "edt.h"

typedef struct
{
  int at;
  int val;
} Xtrm;

typedef struct
{
  int cnt;
  Xtrm *lst;
} XtrmList;

typedef struct
{
  double time;
  int refidx;
  int taridx;
  int count;
  int bc;
  int ndcnt;
  int ref_suridx;
  int tar_suridx;
  int *straddle_at;
  int *notdone;
  signed char *sign;
  XtrmList *min;
  XtrmList *max;
  double *mean;
} Pair;

extern int sur_pair_count;
void init_pair (Pair *p, int refidx, int taridx, int bc, bool read);
bool read_pair (Pair *p, bool all);

#endif
