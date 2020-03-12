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

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <error.h>
#include <errno.h>
#include "pair.h"

bool
read_pair (Pair *p, bool all)
{
  static char *filename;
  free (filename);
  if (asprintf (&filename, "%d_%d", p->refidx, p->taridx) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  if (f == NULL) return false;
  if (fread (p, sizeof (int) * 7 + sizeof (double), 1, f) != 1) error (1, errno, "fread %s", filename);
  if (all || (p->ndcnt > 0 && p->count < sur_pair_count)) {
    if (fread (p->straddle_at, sizeof (int)   , p->bc   , f) != p->bc   )  error (1, errno, "fread %s straddle_at", filename);
    if (fread (p->notdone    , sizeof (int)   , p->ndcnt, f) != p->ndcnt)  error (1, errno, "fread %s notdone"    , filename);
    if (fread (p->sign       , sizeof *p->sign, p->bc   , f) != p->bc   )  error (1, errno, "fread %s sign"       , filename);
    if (fread (p->mean       , sizeof *p->mean, p->bc   , f) != p->bc   )  error (1, errno, "fread %s mean"       , filename);
    void read_list (XtrmList *x) {
      XtrmList m[p->bc];
      if (fread (m, sizeof *m, p->bc, f) != p->bc) error (1, errno, "fread %s min/max", filename);
      for (int i = 0; i < p->bc; i++) {
        x[i].lst = realloc (x[i].lst, (x[i].cnt = m[i].cnt) * sizeof *x[i].lst);
        if (fread (x[i].lst, sizeof *x[i].lst, m[i].cnt, f) != m[i].cnt) error (1, errno, "fread %s min[%d].lst", filename, i);
        if (0) {
          printf ("bin %d cnt %d val ", i, x[i].cnt);
          for (int j = 0; j < x[i].cnt; j++)
            printf ( " %d@%d", x[i].lst[j].val, x[i].lst[j].at);
          printf ("\n");
        }
      }
    }
    read_list (p->min);
    read_list (p->max);
  }
  fclose (f);
  return true;
}

void
init_pair (Pair *p, int refidx, int taridx, int bc, bool read)
{
  p->refidx = refidx;
  p->taridx = taridx;

  XtrmList *allocm (void) {
    XtrmList *x = calloc (bc,  sizeof x[0]);
    for (int i = 0; i < bc; i++)
      x[i].lst = malloc ((x[i].cnt = 1) * sizeof *x[i].lst);
    return x;
  }

  if (p->notdone     == NULL) p->notdone     = malloc (bc * sizeof p->notdone    [0]);
  if (p->straddle_at == NULL) p->straddle_at = calloc (bc,  sizeof p->straddle_at[0]);
  if (p->sign        == NULL) p->sign        = malloc (bc * sizeof p->sign       [0]);
  if (p->min         == NULL) p->min         = allocm (                             );
  if (p->max         == NULL) p->max         = allocm (                             );
  if (p->mean        == NULL) p->mean        = calloc (bc,  sizeof p->mean       [0]);

  void setval (XtrmList *x, int val) {
    for (int i = 0; i < bc; i++) {
      x[i].cnt = 1;
      x[i].lst[0].val = val;
    }
  }

  setval (p->min, INT_MAX);
  setval (p->max, 0      );

  if (read && read_pair (p, false))
    return;
  p->count = 0;
  p->bc = bc;
  p->ndcnt = bc;
  p->ref_suridx = 0;
  p->tar_suridx = 0;
  memset (p->straddle_at, 0, bc * sizeof p->straddle_at[0]);
  for (int i = 0; i < bc; i++)
    p->notdone[i] = i;
}
