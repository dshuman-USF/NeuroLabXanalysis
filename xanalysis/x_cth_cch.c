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
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "util.h"
#include "params.h"

static int *t98;
static int t98cnt;
static int maxclen;
static int *first;
static int *last;

static long long reftime;
static long long tartime;
static long long cchtime;

static double factor;

__inline__ __attribute__((gnu_inline)) uint64_t 
rdtsc()
{
  uint32_t lo, hi;
  /* We cannot use "=A", since this would use %rax on x86_64 */
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return (uint64_t)hi << 32 | lo;
}

static void
get_first_last (int *t, int tcnt)
{
  memset (first, 0xff, t98cnt * sizeof first[0]);
  for (int n = 0; n < t98cnt; n++)
    last[n] = -2;
  int cn = 0;
  int max = 0;
  for (int tn = 0; tn < tcnt; tn++) {
    while (cn + 1 < t98cnt && abs (t98[cn + 1] - t[tn]) <  abs (t98[cn] - t[tn]))
      cn++;
    if (first[cn] == -1)
      first[cn] = tn;
    last[cn] = tn;
    int d = 0;
    if (t[tn] < t98[cn] && cn > 0     ) d = t98[cn] - t98[cn - 1];
    if (t[tn] > t98[cn] && cn < t98cnt) d = t98[cn + 1] - t98[cn];
    if (d > max) max = d;
  }
  maxclen = max;
  if (maxclen < 1565.5 * factor)
    maxclen = 1565.5 * factor;
}

static double *
get_ref_cth (int *t, int tcnt)
{
  long long start;
  long long end;
  start = rdtsc();
  int bincnt = 2 * (maxclen / 2) + 1;
  double *cth = calloc (bincnt, sizeof *cth);
  int cn = 0;
  for (int tn = 0; tn < tcnt; tn++) {
    while (cn + 1 < t98cnt && abs (t98[cn + 1] - t[tn]) <  abs (t98[cn] - t[tn]))
      cn++;
    int bin = maxclen / 2 + (t[tn] - t98[cn]);
    if (0)
      (bin >= 0 && bin < bincnt) || DIE;
    cth[bin]++;
  }
  end   = rdtsc();
  reftime += end - start;
  return cth;
}

static int tstart;

static double *
get_tar_cth2 (int *tar, int tarcnt, int *ref, double *refcth)
{
  long long start;
  long long end;
  start = rdtsc();
  int bincnt = 4 * (maxclen / 2) + 1;
  double *cth = calloc (bincnt, sizeof *cth);
  for (int cn = 0; cn < t98cnt; cn++)
    for (int rn = first[cn]; rn <= last[cn]; rn++) {
      while (tstart < tarcnt && tar[tstart] < ref[rn] - maxclen / 2)
        tstart++;
      for (int tn = tstart; tn < tarcnt && tar[tn] <= ref[rn] + maxclen / 2; tn++)
        cth[ 2 * (maxclen / 2) + (tar[tn] - t98[cn])]++;
    }
  int refbincnt = 2 * (maxclen / 2) + 1;
  double *refsum = malloc ((refbincnt + 1) * sizeof *refsum);
  refsum[0] = 0;
  for (int n = 0; n < refbincnt; n++)
    refsum[n + 1] = refsum[n] + refcth[n];
  for (int tn = 0; tn < bincnt; tn++) {
    if (cth[tn]) {
      int rnl = MAX (tn - 2 * (maxclen / 2), 0);
      int rnr = MIN (tn, 2 * (maxclen / 2));
      cth[tn] /= refsum[rnr + 1] - refsum[rnl];
    }
  }
  free (refsum);
  end   = rdtsc();
  tartime += end - start;
  return cth;
}

double *
get_tar_cth (int *t, int tcnt, int *reft, double *refcth)
{
  long long start;
  long long end;
  start = rdtsc();
  int bincnt = 4 * (maxclen / 2) + 1;
  double *cth = calloc (bincnt, sizeof *cth);
  int cn = 0;
  for (int tn = 0; tn < tcnt; tn++) {
    while (cn < t98cnt && t[tn] - t98[cn] > 2 * (maxclen / 2))
      cn++;
    for (int ci = cn; ci < t98cnt && t[tn] - t98[ci] > -2 * (maxclen / 2); ci++) {
      int bin = 2 * (maxclen / 2) + (t[tn] - t98[ci]);
      if (0)
        (bin >= 0 && bin < bincnt) || DIE;
      int refcnttmp = 0;
      for (int n = first[ci]; n <= last[ci]; n++)
        if (n >= 0 && abs (reft[n] - t[tn]) <= maxclen / 2)
          refcnttmp++;
      int reftottmp = 0;
      for (int refbin = 0; refbin < 2 * (maxclen / 2) + 1; refbin++) {
        if (refcth[refbin] == 0)
          continue;
        int refoffset = refbin - maxclen / 2;
        int taroffset = bin - 2 * (maxclen / 2);
        if (abs (refoffset - taroffset) <= maxclen / 2)
          reftottmp += refcth[refbin];
      }
      if (reftottmp) cth[bin] += (double) refcnttmp / reftottmp;
    }
  }
  end   = rdtsc();
  tartime += end - start;
  return cth;
}

double *
get_cth_cch (double *refcth, double *tarcth, double **varp)
{
  long long start;
  long long end;
  start = rdtsc();
  int bincnt    = 2 * (maxclen / 2) + 1;
  int refbincnt = 2 * (maxclen / 2) + 1;
  int tarbincnt = 4 * (maxclen / 2) + 1;
  
  double *cch = calloc (bincnt, sizeof *cch);
  double *var = calloc (bincnt, sizeof *var);
  for (int refbin = 0; refbin < refbincnt; refbin++) {
    if (refcth[refbin] == 0)
      continue;
    for (int cchoffset = -maxclen / 2; cchoffset <= maxclen / 2; cchoffset++) {
      int tarbin = maxclen / 2 + refbin + cchoffset;
      if (tarcth[tarbin] == 0)
        continue;
      int cchbin = maxclen / 2 + cchoffset;
      if (0) {
        (tarbin >= 0 && tarbin < tarbincnt) || DIE;
        (refbin >= 0 && refbin < refbincnt) || DIE;
        (cchbin >= 0 && cchbin < bincnt) || DIE;
      }
      cch[cchbin] += refcth[refbin] * tarcth[tarbin];
      var[cchbin] += refcth[refbin] * tarcth[tarbin] * (1 - tarcth[tarbin]);
    }
  }
  *varp = var;
  end   = rdtsc();
  cchtime += end - start;
  return cch;
}

static inline void
d2i (int *i, double *d, int count)
{
  for (int n = 0; n < count; n++)
    i[n] = rint (d[n] * factor);
}

static void
get_ref_tar_and_e (void *stp, int *mne, int *ital, int ref, int tar, int epl,
                   int **rp, int *rcntp, int **tp, int *tcntp, int **ep, int *ecntp)
{
  double (*spiketimes)[*mne] = stp;

#define CHK(x) (x >= 0 && x < *mne) || DIE
  CHK(ital[ref]);
  CHK(ital[tar]);
  CHK(ital[epl]);
#undef CHK
  int *r = malloc (ital[ref] * sizeof *r);
  int *t = malloc (ital[tar] * sizeof *t);
  int *e = malloc (ital[epl] * sizeof *e);
  d2i (r, spiketimes[ref], ital[ref]);
  d2i (t, spiketimes[tar], ital[tar]);
  d2i (e, spiketimes[epl], ital[epl]);

  *rp = r; *rcntp = ital[ref];
  *tp = t; *tcntp = ital[tar];
  *ep = e; *ecntp = ital[epl];
}

static double *cth_cch, *var;

void
dump_hist (char *fname, int *h)
{
  FILE *g = fopen (fname, "w");
  for (int n = 0; n < NUM_BINS; n++)
    fprintf (g, "%d\n", h[n]);
  fclose (g);
}

static void
rebin (double *d, int mid, double old_binw, double new_binw, int *ihist, int (*f)(double))
{
#if NUM_BINS % 2 != 1
#error "NUM_BINS needs to be odd"
#endif
  new_binw = nearbyint (new_binw);
  double hist[NUM_BINS];
  memset (hist, 0, sizeof hist);
  double new_max_delta = (NUM_BINS - 1) / 2 * new_binw + new_binw / 2;
  double old_first_time = -(mid + .5) * old_binw;
  int first_n = (-new_max_delta - old_first_time) / old_binw;
  int last_n =  ( new_max_delta - old_first_time) / old_binw;
  double new_first_time = -new_max_delta;
  int dcnt = mid * 2 + 1;
  int lskip = 0;
  int rskip = NUM_BINS;
  for (int n = first_n; n < last_n; n++) {
    double time = (n - mid) * old_binw;
    int bin = (time - new_first_time) / new_binw;
    (bin >= 0 && bin < NUM_BINS) || DIE;
    if (n < 0) lskip = bin + 1;
    else if (n >= dcnt) {if (rskip == NUM_BINS) rskip = bin;}
    else hist[bin] += d[n];
  }
  for (int n = 0; n < lskip; n++)
    hist[n] = 0;
  for (int n = rskip; n < NUM_BINS; n++)
    hist[n] = 0;
  for (int n = 0; n < NUM_BINS; n++)
    ihist[n] = f (hist[n]);
  return;

  printf ("%s line %d: writing histogram to \"dat\"\n", __FILE__, __LINE__);
  static int passcnt;
  char *t;
  char *c[] = {"", "black", "red"};
  if (new_binw == 5) {
    printf ("pass %d\n", passcnt);
    if (asprintf (&t, "pass %d", ++passcnt) == -1) exit (1);
    void preploti (int *y, int *cntp, char *title, char *color);
    void preplotf (double *y, int *cntp, char *title, char *color);
    int count = 101;
    if (passcnt == 1)
      preploti (ihist, &count, t, c[passcnt]) ;
    if (passcnt == 2) {
      double sd[NUM_BINS];
      union {int i; float f;} v;
      for (int n = 0; n < NUM_BINS; n++) {
        //        sd[n] = sqrt (hist[n]);
        v.i = ihist[n];
        sd[n] = v.f;
      }
      FILE *g = fopen ("dat", "w");
      for (int n = 0; n < NUM_BINS; n++)
        fprintf (g, "%d\n", ihist[n]);
      fclose (g);
      preplotf (sd, &count, t, c[passcnt]);
    }
    
  }
}

void
rebin_cth_cch (float *binwp, int *ihist)
{
  int f (double x) {return (int) rint (x);}
  rebin (cth_cch, maxclen / 2, 1, *binwp * factor, ihist, f);
}

void
rebin_var (float *binwp, int *ihist)
{
  typedef union {int i; float f;} Fi;
  int f (double x) {Fi v; v.f = sqrt (x); return v.i;}
  rebin (var, maxclen / 2, 1, *binwp * factor, ihist, f);
}

void
free_cth_cch (void)
{
  free (cth_cch);
  free (var);
}

static time_t start_time;
static int pass_count;

static void
add_cch (int *tref, int trefcnt, int *ttar, int ttarcnt)
{
  double *tmpvar;
  double *refcth = get_ref_cth (tref, trefcnt);
  double *tarcth = get_tar_cth2 (ttar, ttarcnt, tref, refcth);
  double *tmpcch = get_cth_cch (refcth, tarcth, &tmpvar);

  for (int n = 0; n < 2 * (maxclen / 2) + 1; n++)
    cth_cch[n] += tmpcch[n];
  for (int n = 0; n < 2 * (maxclen / 2) + 1; n++)
    var[n] += tmpvar[n];
  free (refcth);
  free (tarcth);
  free (tmpcch);
  free (tmpvar);
}


void
calc_cth_cch (void *stp, int *mne, int *ital, int *refp, int *tarp, int *eplp, double *factorp)
{
  double (*spiketimes)[*mne] = stp;
  
  if (start_time == 0)
    start_time = time (0);
  factor = *factorp;
  int ref = *refp - 1;
  int tar = *tarp - 1;
  int epl = *eplp - 1;
#define CHK(x) (x >= 0 && x < MAX_NUM_CHAN) || DIE
  CHK(ref);
  CHK(tar);
  CHK(epl);
#undef CHK

  int orig_trefcnt, orig_ttarcnt, orig_t98cnt, *orig_tref, *orig_ttar, *orig_t98;
  get_ref_tar_and_e (spiketimes, mne, ital, ref, tar, epl, &orig_tref, &orig_trefcnt, &orig_ttar, &orig_ttarcnt, &orig_t98, &orig_t98cnt);

  t98 = orig_t98;
  t98cnt = orig_t98cnt;
  TREALLOC (first, orig_t98cnt);
  TREALLOC (last, orig_t98cnt);
  get_first_last (orig_tref, orig_trefcnt);
  int *orig_first = first, *orig_last = last;

  int bincnt    = 2 * (maxclen / 2) + 1;
  cth_cch = calloc (bincnt, sizeof *cth_cch);
  var     = calloc (bincnt, sizeof *var);

  int orig_tstart = 0;
  for (int t98idx = 0; t98idx < orig_t98cnt; t98idx += 20) {
    int trefcnt, ttarcnt, *tref, *ttar;
    t98 = orig_t98 + t98idx;
    t98cnt = orig_t98cnt - t98idx;
    if (t98cnt > 20) t98cnt = 20;
    first = orig_first + t98idx;
    last = orig_last + t98idx;
      
    int fn, ln;
    for (fn = 0;      fn < t98cnt; fn++) if (first[fn] >= 0) break;
    for (ln = t98cnt - 1; ln >= 0; ln--) if ( last[ln] >= 0) break;
    if (ln < fn) continue;
    
    tref = orig_tref + first[fn];
    trefcnt = last[ln] - first[fn] + 1;

    int ftmp[20], ltmp[20];
    for (int n = 0; n < t98cnt; n++) {
      ftmp[n] = first[n] == -1 ? -1 : first[n] - first[fn];
      ltmp[n] =  last[n] == -2 ? -2 : last [n] - first[fn];
    }
    first = ftmp;
    last = ltmp;

    ttar = orig_ttar + orig_tstart;
    ttarcnt = orig_ttarcnt - orig_tstart;
    add_cch (tref, trefcnt, ttar, ttarcnt);
    orig_tstart += tstart;
    tstart = 0;
  }
  first = orig_first;
  last = orig_last;
  free (orig_tref);
  free (orig_ttar);
  free (orig_t98);
  time_t now = time(0);
  pass_count++;
  if (0) {
    printf ("%d pairs, %ld seconds, %g seconds per pair\n", pass_count, now - start_time, (double)(now - start_time) / pass_count);
    printf ("ref, tar, cch times: %11.9f %7.5f %8.5f\n", reftime / 2.4e9, tartime / 2.4e9, cchtime / 2.4e9);
  }
}
