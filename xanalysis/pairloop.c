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
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <error.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <sys/utsname.h>
#include "edt.h"
#include "tools.h"
#include "cch.h"
#include "surgen.h"
#include "pairloop.h"

int *bw;
int hc;

int sur_pair_count = 105178;

static inline uint64_t 
tick (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000000000 + ts.tv_nsec;
}

static bool
hcmp (Pair *p, int *dat, int *sur)
{
  p->count++;
  int k = 0;
  if (p->count == 1) {
    for (int i = 0; i < p->ndcnt; i++) {
      int j = p->notdone[i];
      if (sur[j] == dat[j])
        p->straddle_at[j] = p->count;
      else {
        p->sign[j] = sur[j] < dat[j] ? 1 : -1;
        p->notdone[k++] = j;
      }
    }
  }
  else {
    for (int i = 0; i < p->ndcnt; i++) {
      int j = p->notdone[i];
      if ((sur[j] - dat[j]) * p->sign[j] >= 0)
        p->straddle_at[j] = p->count;
      else
        p->notdone[k++] = j;
    }
  }
  p->ndcnt = k;
  return k == 0;
}

static void
write_pair (Pair *p)
{
  char *filename;
  if (asprintf (&filename, "%d_%d", p->refidx, p->taridx) == -1) exit (1);
  FILE *f = fopen (filename, "w");
  if (f == NULL) error (1, errno, "fopen %s write", filename);
  if (fwrite (p             , sizeof (int) * 7 + sizeof (double), 1       , f) != 1       ) error (1, errno, "fwrite %s"            , filename);
  if (fwrite (p->straddle_at, sizeof (int)                      , p->bc   , f) != p->bc   ) error (1, errno, "fwrite %s straddle_at", filename);
  if (fwrite (p->notdone    , sizeof (int)                      , p->ndcnt, f) != p->ndcnt) error (1, errno, "fwrite %s notdone"    , filename);
  if (fwrite (p->sign       , sizeof *p->sign                   , p->bc   , f) != p->bc   ) error (1, errno, "fwrite %s sign"       , filename);
  if (fwrite (p->mean       , sizeof *p->mean                   , p->bc   , f) != p->bc   ) error (1, errno, "fwrite %s mean"       , filename);

  void write_list (XtrmList *x) {
    if (fwrite (x, sizeof *x, p->bc, f) != p->bc) error (1, errno, "fwrite %s min/max", filename);
    for (int i = 0; i < p->bc; i++)
      if (fwrite (x[i].lst, sizeof *x[i].lst, x[i].cnt, f) != x[i].cnt) error (1, errno, "fwrite %s min[%d].lst", filename, i);
  }
  write_list (p->min);
  write_list (p->max);

  free (filename);
  fclose (f);
}

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
    //    printf ("read_pair %d %d: straddle_at[31] %d\n", p->refidx, p->taridx, p->straddle_at[31]);
    void read_list (XtrmList *x) {
      XtrmList m[p->bc];
      if (fread (m, sizeof *m, p->bc, f) != p->bc) error (1, errno, "fread %s min/max", filename);
      for (int i = 0; i < p->bc; i++) {
        x[i].lst = realloc (x[i].lst, (x[i].cnt = m[i].cnt) * sizeof *x[i].lst);
        if (fread (x[i].lst, sizeof *x[i].lst, m[i].cnt, f) != m[i].cnt) error (1, errno, "fread %s min[%d].lst", filename, i);
        if (0)
          if (i == 31) {
            printf ("read_pair: bin %d cnt %d val ", i, x[i].cnt);
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
  //  printf ("read_pair %d %d: straddle_at[0] %d\n", p->refidx, p->taridx, p->straddle_at[0]);
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

void
minmaxmean (Pair *p, int *hs)
{
  for (int i = 0; i < p->bc; i++) {
    void update (XtrmList *x) {
      if (p->count >= 40)
        x->lst = realloc (x->lst, ++x->cnt * sizeof *x->lst);
      x->lst[x->cnt - 1].at = p->count + 1;
      x->lst[x->cnt - 1].val = hs[i];
    }
    XtrmList *x;
    x = &p->min[i]; if (hs[i] < x->lst[x->cnt - 1].val) update (x);
    x = &p->max[i]; if (hs[i] > x->lst[x->cnt - 1].val) update (x);
    p->mean[i] += (hs[i] - p->mean[i]) / (p->count + 1);
  }
}

static char *
lock (int refidx, int taridx)
{
  static char *lockname;
  free (lockname);
  if (asprintf (&lockname, "lock_%d_%d", refidx, taridx) == -1) exit (1);
  FILE *f = fopen (lockname, "wx");
  if (f == NULL) {
    if (errno == EEXIST)
      return NULL;
    error (1, errno, "fopen %s", lockname);
  }
  struct utsname info;
  uname (&info);
  fprintf (f, "%s %d %ld\n", info.nodename, getpid (), time (0));
  fclose (f);
  return lockname;
}

void
pair_loop (Edt *edt)
{
  int bc = hc * BC;
  Pair p = {0};
  int *hset = malloc (bc * sizeof *hset);
  char *stopfile;
  if (asprintf (&stopfile, "stop_%d", getpid ()) == -1) exit (1);

  uint64_t start = tick ();
  for (int refidx = 0; refidx < edt->digital_count; refidx++)
    for (int taridx = refidx + 1; taridx < edt->digital_count; taridx++) {
      char *lockname;
      if ((lockname = lock (refidx, taridx)) == NULL)
        continue;
      Digital *ref = edt->digital[refidx];
      Digital *tar = edt->digital[taridx];
      SpikeTrain ref_st = {.N = ref->spike_count, .T = ref->spike};
      SpikeTrain tar_st = {.N = tar->spike_count, .T = tar->spike};
      init_pair (&p, refidx, taridx, bc, true);
      if (p.ndcnt > 0 && p.count < sur_pair_count) {
        int *h = cch_set (ref_st, tar_st, bw, 1, bw[hc - 1] * BC, hc);
        memcpy (hset, h, bc * sizeof *hset);
        while (p.count < sur_pair_count) {
          SpikeTrain r_st, t_st;
          get_sur_pair (edt, &p, &r_st, &t_st);
          int *hs = cch_set (r_st, t_st, bw, 1, bw[hc - 1] * BC, hc);
          minmaxmean (&p, hs);
          if (hcmp (&p, hset, hs))
            break;
        }
        uint64_t now = tick ();
        p.time = (now - start) / 1e9;
        start = now;
        write_pair (&p);
      }
      unlink (lockname);
      if (access (stopfile, F_OK) == 0) {
        unlink (stopfile);
        return;
      }
      if (access ("stop", F_OK) == 0)
        return;
    }
}

void
pair_status (Edt *edt, int threshold)
{
  int bc = hc * BC;
  Pair p = {0};
  init_pair (&p, 0, 0, bc, false);
  double time = 0;
  int paircnt = 0;
  int sigpaircnt = 0;
  int sigbincnt = 0;
  int ncnt[sur_pair_count + 1];
  memset (ncnt, 0, sizeof ncnt);
  int cellpaircount = edt->digital_count * (edt->digital_count - 1) / 2;
  int max_straddle[cellpaircount];
  memset (max_straddle, 0, sizeof max_straddle);
  int cpidx = -1;
  for (int refidx = 0; refidx < edt->digital_count; refidx++)
    for (int taridx = refidx + 1; taridx < edt->digital_count; taridx++) {
      cpidx++;
      p.refidx = refidx;
      p.taridx = taridx;
      if (!read_pair (&p, true))
        continue;
      paircnt++;
      sigpaircnt += p.ndcnt > 0;
      sigbincnt += p.ndcnt;
      time += p.time;
      for (int i = 0; i < bc; i++) {
        if (!(p.straddle_at[i] >= 0 && p.straddle_at[i] <= sur_pair_count))
          error_at_line (1, 0, __FILE__, __LINE__, "straddle = %d, sur_pair_count = %d\n", p.straddle_at[i], sur_pair_count);
        ncnt[p.straddle_at[i]]++;
        if (p.straddle_at[i] > max_straddle[cpidx])
          max_straddle[cpidx] = p.straddle_at[i];
      }
      if (p.ndcnt > 0)
        max_straddle[cpidx] = sur_pair_count;

      if (max_straddle[cpidx] >= threshold) {
        printf ("%3d %3d:", edt->digital[refidx]->id,  edt->digital[taridx]->id);
        int lastcch = 0, lastbin = 0, startbin = 0;
        for (int i = 0; i < p.ndcnt; i++) {
          int j = p.notdone[i];
          int cch = j / 101 + 1;
          int bin = j % 101 + 1;
          if (cch != lastcch) {
            if (startbin && lastbin != startbin)
              printf ("-%d", lastbin);
            printf ("\n        %d: %d", cch, bin);
            startbin = bin;
          }
          else if (startbin && bin != lastbin + 1) {
            if (lastbin != startbin)
              printf ("-%d", lastbin);
            printf (" %d", bin);
            startbin = bin;
          }
          else if (startbin == 0) {
            printf (" %d", bin);
            startbin = bin;
          }
          lastcch = cch;
          lastbin = bin;
        }
        if (startbin && lastbin != startbin)
          printf ("-%d", lastbin);
        printf ("\n");
      }
      if (0)
      printf ("refidx %3d, taridx %3d, count %6d, ndcnt %3d, suridx %3d %3d, secs %8.4f\n",
              p.refidx, p.taridx, p.count, p.ndcnt, p.ref_suridx, p.tar_suridx, p.time);
    }
  printf ("\n%d pairs, %d sig pairs, %d sig bins, %.0f seconds, %g sec/pair, %g hrs/6903pr\n",
          paircnt, sigpaircnt, sigbincnt, time, time / paircnt, 6903 * time / paircnt / 3600);
  FILE *f = fopen ("/tmp/ncnt", "w");
  fwrite (ncnt, sizeof ncnt, 1, f);
  fclose (f);
  f = fopen ("/tmp/max_straddle", "w");
  for (int i = 0; i < cellpaircount; i++)
    fprintf (f, "%d\n", max_straddle[i]);
  fclose (f);
}

void
threshold (Edt *edt)
{
  int bc = hc * BC;
  Pair p = {0};
  init_pair (&p, 0, 0, bc, false);
  int ncnt[sur_pair_count + 1];
  memset (ncnt, 0, sizeof ncnt);
  for (int refidx = 0; refidx < edt->digital_count; refidx++)
    for (int taridx = refidx + 1; taridx < edt->digital_count; taridx++) {
      p.refidx = refidx;
      p.taridx = taridx;
      if (!read_pair (&p, true))
        error (1, 0, "missing pair\n");
      for (int i = 0; i < bc; i++) {
        //        if (!(p.straddle_at[i] >= 0 && p.straddle_at[i] <= sur_pair_count))
        if (p.straddle_at[i] < 0)
          error_at_line (1, 0, __FILE__, __LINE__, "straddle = %d\n", p.straddle_at[i]);
        if (p.straddle_at[i] >= sur_pair_count)
          ncnt[0]++;
        else
          ncnt[p.straddle_at[i]]++;
      }
    }

  int N = 0;
  for (int i = 0; i < sur_pair_count; i++)
    N += ncnt[i];
  ncnt[sur_pair_count] = ncnt[0];
  
  int sum = N;
  for (int i = 2; i <= sur_pair_count; i++) {
    if (2. / i * N / sum <= .05) {
      printf ("%d\n", i);
      return;
    }
    sum -= ncnt[i];
  }
  printf ("%d\n", sur_pair_count + 1);
}

void
pair_count (Edt *edt)
{
  int bc = hc * BC;
  Pair p = {0};
  init_pair (&p, 0, 1, bc, true);
  p.ref_suridx = 0;
  p.tar_suridx = 0;
  printf ("%3d %3d\n", p.ref_suridx, p.tar_suridx);
  for (int i = 0; i < sur_pair_count; i++) {
    SpikeTrain r_st, t_st;
    get_sur_pair (edt, &p, &r_st, &t_st);
    printf ("%3d %3d\n", p.ref_suridx, p.tar_suridx);
    if (p.ref_suridx >= 325 || p.tar_suridx >= 325) {
      printf ("ref_suridx = %d, tar_suridx = %d, at i = %d\n", p.ref_suridx, p.tar_suridx, i);
      exit (0);
    }
  }
}

void
pair_time (Edt *edt)
{
  int bc = hc * BC;
  Pair p = {0};
  int *hset = malloc (bc * sizeof *hset);

  uint64_t start = tick ();
  {
    int refidx = 24;
    int taridx = 35;
    init_pair (&p, refidx, taridx, bc, false);
    Digital *ref = edt->digital[refidx];
    Digital *tar = edt->digital[taridx];
    SpikeTrain ref_st = {.N = ref->spike_count, .T = ref->spike};
    SpikeTrain tar_st = {.N = tar->spike_count, .T = tar->spike};
    {
      int *h = cch_set (ref_st, tar_st, bw, 1, bw[hc - 1] * BC, hc);
      memcpy (hset, h, bc * sizeof *hset);

      int64_t loop_start = tick ();
      int loop_count = 0;
      int N = 10000;
      int64_t sum1 = 0;
      int64_t sum2 = 0;
      int64_t sum3 = 0;
      while (p.count < N) {
        loop_count++;
        SpikeTrain r_st, t_st;
        int64_t t0 = tick ();
        get_sur_pair (edt, &p, &r_st, &t_st);
        int64_t t1 = tick ();
        sum1 += t1 - t0;
        int *hs = cch_set (r_st, t_st, bw, 1, bw[hc - 1] * BC, hc);
        int64_t t2 = tick ();
        sum2 += t2 - t1;
        if (hcmp (&p, hset, hs))
          break;
        int64_t t3 = tick ();
        sum3 += t3 - t2;
      }
      uint64_t now = tick ();
      printf ("%ld %ld %ld\n", sum1, sum2, sum3);
      printf ("%d %d (%d %d): %d cch sets, %g ms each\n",
              refidx, taridx, edt->digital[refidx]->id, edt->digital[taridx]->id, loop_count, (now - loop_start) / 1e6 / N);
        

      p.time = (now - start) / 1e9;
      start = now;
    }
  }
}


void
pair_plot (Edt *edt, int refidx, int taridx, int spcnt)
{
  int bc = hc * BC;
  Pair p = {0};
  int *hset = malloc (bc * sizeof *hset);
  int *min = malloc (bc * sizeof *min);
  int *max = calloc (bc, sizeof *max);
  double *mean = calloc (bc, sizeof *mean);
  for (int i = 0; i < bc; i++)
    min[i] = INT_MAX;
  Digital *ref = edt->digital[refidx];
  Digital *tar = edt->digital[taridx];
  SpikeTrain ref_st = {.N = ref->spike_count, .T = ref->spike};
  SpikeTrain tar_st = {.N = tar->spike_count, .T = tar->spike};
  init_pair (&p, refidx, taridx, bc, false);
  read_pair (&p, true);
  printf ("pair_loop %d %d: straddle_at[0] %d\n", p.refidx, p.taridx, p.straddle_at[31]);

  int *h = cch_set (ref_st, tar_st, bw, 1, bw[hc - 1] * BC, hc);
  memcpy (hset, h, bc * sizeof *hset);
  for (int i = 0; i < bc; i++) {
    int getval (XtrmList *m) {
      int j = m->cnt - 1;
      while (j > 0 && m->lst[j].at > spcnt)
        j--;
      return m->lst[j].val;
    }
    if (i == 31) {
      printf ("max:\n");
      XtrmList *m = &p.max[i];
      for (int j = 0; j < m->cnt; j++)
        printf ("%d at %d\n", m->lst[j].val, m->lst[j].at);
      printf ("min:\n");
      m = &p.min[i];
      for (int j = 0; j < m->cnt; j++)
        printf ("%d at %d\n", m->lst[j].val, m->lst[j].at);
    }
    min[i] = getval (&p.min[i]);
    max[i] = getval (&p.max[i]);
    mean[i] = p.mean[i];
  }
  printf ("pair_loop %d %d: straddle_at[0] %d\n", p.refidx, p.taridx, p.straddle_at[31]);
  for (int hi = 0; hi < hc; hi++) {
    char *filename;
    if (asprintf (&filename, "h%d", hi + 1) == -1) exit (1);
    FILE *f = fopen (filename, "w");
    free (filename);
    for (int bi = hi * 101; bi <  (hi + 1) * 101; bi++) {
      if (p.straddle_at[bi] > 0 && p.straddle_at[bi] <= spcnt && (max[bi] < hset[bi] || min[bi] > hset[bi]))
          error_at_line (1, 0, __FILE__, __LINE__, 
                         "pair %d %d (%d %d), bin %d, min %d, max %d, data %d, straddle_at: %d\n",
                         refidx, taridx, ref->id, tar->id, bi, min[bi], max[bi], hset[bi], p.straddle_at[bi]);
                         
      fprintf (f, "%d %d %d %g\n", hset[bi], min[bi], max[bi], mean[bi]);
    }
    fclose (f);
  }
  {
    FILE *f[9];
    for (int i = 1; i <= hc; i++) {
      f[i] = popen ("gnuplot", "w");
      fprintf (f[i], "set term x11 pos %d,%d size 392,400\n", ((i-1)%4) * 400, ((i-1)/4) * 468);
      fprintf (f[i], "set title '%d > %d, %.1f ms, %d surrogates'\n", edt->digital[refidx]->id, edt->digital[taridx]->id, bw[i - 1] / 10., (int)fmin (spcnt, p.count));
      //      printf ("plot 'h%d' u 1 w his lc 'black' t '', 'h%d' u 2 w his lc rgb 'red' t '', 'h%d' u 3 w his  lc rgb 'red' t '', 'h%d' u 4 lc rgbc 'red' w his  t ''\n", i, i, i, i);
      fprintf (f[i], "plot 'h%d' u 1 w his lc 'black' t '', 'h%d' u 2 w his lc rgb 'red' t '', 'h%d' u 3 w his  lc rgb 'red' t '', 'h%d' u 4 lc rgbc 'red' w his  t ''\n", i, i, i, i);
      fflush (f[i]);
    }
    printf ("plot done, cr when ready\n");
    getchar ();
    for (int i = 1; i <= hc; i++)
      pclose (f[i]);
  }
}

#define SGN(a) (((a) < 0) ? -1 : ((a) > 0) ? 1 : 0)

static double
bumpsize (int *data, double *mean, int startbin, int nbins)
{
  double sum = 0;
  for (int i = startbin; i < startbin + nbins; i++)
    sum += data[i] - mean[i];
  int sgn = SGN (sum);
  for (int i = startbin; i < startbin + nbins; i++)
    if (SGN (data[i] - mean[i]) == -sgn)
      return 0;
  return sum;
}

static void
maxbumploc (int *data, double *mean, int hidx, int *startbin, int *nbins)
{
  int maxstartbin = 0;
  int maxnbins = 0;
  double maxbinsum = 0;
  int curstartbin = 0;
  int curnbins = 0;
  double curbinsum = 0;
  int current_sgn = 2;

  void checkmax (void) {
    if (curbinsum > maxbinsum && curnbins > 1) {
      maxstartbin = curstartbin;
      maxnbins = curnbins;
      maxbinsum = curbinsum;
    }
  }
  for (int i = hidx * BC; i < (hidx + 1) * BC; i++) {
    int sgn = SGN (data[i] - mean[i]);
    if (sgn == 0)
      continue;
    if (sgn != current_sgn) {
      checkmax ();
      current_sgn = sgn;
      curstartbin = i;
      curnbins = 1;
      curbinsum = fabs (data[i] - mean[i]);
    }
    else {
      curbinsum += fabs (data[i] - mean[i]);
      curnbins++;
    }
  }
  checkmax ();
  *startbin = maxstartbin;
  *nbins = maxnbins;
}

void
pair_bump (Edt *edt, int refidx, int taridx)
{
  int bc = hc * BC;
  Pair p = {0};
  int *hset = malloc (bc * sizeof *hset);
  double *mean = calloc (bc, sizeof *mean);
  Digital *ref = edt->digital[refidx];
  Digital *tar = edt->digital[taridx];
  SpikeTrain ref_st = {.N = ref->spike_count, .T = ref->spike};
  SpikeTrain tar_st = {.N = tar->spike_count, .T = tar->spike};
  init_pair (&p, refidx, taridx, bc, false);
  read_pair (&p, true);
  int *h = cch_set (ref_st, tar_st, bw, 1, bw[hc - 1] * BC, hc);
  memcpy (hset, h, bc * sizeof *hset);
  memcpy (mean, p.mean, bc * sizeof *mean);
  int startbin, nbins;
  maxbumploc (hset, mean, hc - 1, &startbin, &nbins);
  double databump = bumpsize (hset, mean, startbin, nbins);
  double sgn = SGN (databump);
  double maxsurr = 0;
  int surpaircnt = 10913;
  bool toldhim = false;
  for (int i = 0; i < surpaircnt; i++) {
    SpikeTrain r_st, t_st;
    get_sur_pair (edt, &p, &r_st, &t_st);
    int *hs = cch_set (r_st, t_st, bw, 1, bw[hc - 1] * BC, hc);
    double surrbump = bumpsize (hs, mean, startbin, nbins);
    if (sgn * surrbump > maxsurr)
      maxsurr = sgn * surrbump;
    if (!toldhim && maxsurr >= databump) {
      printf ("not significant at surrogate pair %d\n", i + 1);
      toldhim = true;

      FILE *f = fopen ("h0", "w");
      for (int bi = (hc - 1) * 101; bi <  hc * 101; bi++)
        fprintf (f, "%d\n", hs[bi]);
      fclose (f);

    }
  }
  printf ("max bump: %d bins starting at bin %d\n", nbins, (startbin + 1) % BC);
  printf ("size of data bump:                                             %g\n", databump);
  printf ("size of most extreme surrogate bump in %6d surrogate pairs: %g\n", surpaircnt, maxsurr);
  printf ("data bump %s significant\n", maxsurr >= databump ? "is not" : "is");
}
