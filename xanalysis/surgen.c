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
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include "tools.h"
#include "edt.h"
#include "surgen.h"
#include "pairloop.h"


typedef struct 
{
  int cellidx;
  int surcnt;
  int suridx;
  bool changed;
  SpikeTrain *sur;
} SurList;

static double *g;
static RateFunc *ri;

static inline uint64_t 
tick (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000000000 + ts.tv_nsec;
}

static char *
lock (int cellidx)
{
  static char *lockname;
  if (asprintf (&lockname, "lock_%d", cellidx) == -1) exit (1);
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

static void
load (SurList *sl, int cellidx)
{
  char *filename;
  if (sl->cellidx != -1 && sl->surcnt && sl->changed) {
    if (asprintf (&filename, "%d", sl->cellidx) == -1) exit (1);
    FILE *f = fopen (filename, "w");
    free (filename);
    if (f == NULL) error_at_line (1, errno, __FILE__, __LINE__, "fopen");
    if (fwrite (sl, sizeof (int), 3, f) != 3) error_at_line (1, errno, __FILE__, __LINE__, "fwrite");
    if (fwrite (sl->sur, sizeof (SpikeTrain), sl->surcnt, f) != sl->surcnt) error_at_line (1, errno, __FILE__, __LINE__, "fwrite");
    for (int i = 0; i < sl->surcnt; i++)
      if (fwrite (sl->sur[i].T, sizeof (int), sl->sur[i].N, f) != sl->sur[i].N) error_at_line (1, errno, __FILE__, __LINE__, "fwrite");
    fclose (f);
  }
  if (sl->surcnt) {
    for (int i = 0; i < sl->surcnt; i++)
      free (sl->sur[i].T);
    free (sl->sur);
    sl->surcnt = 0;
  }
  if (cellidx == -1) {
    sl->cellidx = -1;
    return;
  }
  if (asprintf (&filename, "%d", cellidx) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  free (filename);
  if (f == NULL) {
    sl->cellidx = cellidx;
    sl->surcnt = 0;
    sl->suridx = 0;
    sl->changed = false;
    sl->sur = NULL;
    return;
  }
  if (fread (sl, sizeof (int), 3, f) != 3) error_at_line (1, errno, __FILE__, __LINE__, "fread");
  sl->sur = malloc (sl->surcnt * sizeof (SpikeTrain));
  if (fread (sl->sur, sizeof (SpikeTrain), sl->surcnt, f) != sl->surcnt) error_at_line (1, errno, __FILE__, __LINE__, "fread");
  for (int i = 0; i < sl->surcnt; i++) {
    sl->sur[i].T = malloc (sl->sur[i].N * sizeof (int));
    if (fread (sl->sur[i].T, sizeof (int), sl->sur[i].N, f) != sl->sur[i].N) error_at_line (1, errno, __FILE__, __LINE__, "fread");
  }
  fclose (f);
  sl->changed = false;
}

RateFunc
get_ratefunc (SpikeTrain st, double g, int id)
{
  static char *filename;
  free (filename);
  if (asprintf (&filename, "rf_%03d", id) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    RateFunc rf = recip_isi (st, g);
    f = fopen (filename, "w");
    if (f == NULL) error_at_line (1, errno, __FILE__, __LINE__, "fopen %s", filename);
    if (fwrite (&rf.N, sizeof rf.N, 1, f) != 1) error_at_line (1, errno, __FILE__, __LINE__, "fwrite %s", filename);
    if (fwrite (rf.r, sizeof rf.r[0], rf.N, f) != rf.N) error_at_line (1, errno, __FILE__, __LINE__, "fwrite %s", filename);
    fclose (f);
    return rf;
  }
  int tcnt;
  if (fread (&tcnt, sizeof tcnt, 1, f) != 1) error_at_line (1, errno, __FILE__, __LINE__, "fread");
  double *y = malloc (tcnt * sizeof *y);
  if (fread (y, sizeof *y, tcnt, f) != tcnt) error_at_line (1, errno, __FILE__, __LINE__, "fread");
  fclose (f);
  return (RateFunc){.r = y, .N = tcnt};
}

static double alloc_time;

static void
addsur (Edt *edt, SurList *sl)
{
  sl->suridx = sl->surcnt++;
  sl->changed = true;
  int64_t start = tick ();
  sl->sur = realloc (sl->sur, sl->surcnt * sizeof sl->sur[0]);
  alloc_time += (tick () - start) / 1e9;
  int n = sl->cellidx;

  SpikeTrain st = {.T = edt->digital[n]->spike, .N = edt->digital[n]->spike_count};
  if (g[n] == 0)
    g[n] = miura_gamma (st);
  if (ri[n].r == 0)
    ri[n] = get_ratefunc (st, g[n], edt->digital[n]->id);

  sl->sur[sl->suridx] = spikes_from_rate (ri[n], lrint (g[n] < 1 ? 1.0 : g[n]));
}

void
generate_surrogates (Edt *edt, int surcnt)
{
  //  printf ("surcnt: %d, cell count: %d\n", surcnt, edt->digital_count);
  seed_rng (time (0));
  SurList sl = {.cellidx = -1, .surcnt = 0};
  g = calloc (edt->digital_count, sizeof *g);
  ri = calloc (edt->digital_count, sizeof *ri);
  set_ticks_per_sec (strcmp (edt->format, "%5d%10d\n") == 0 ? 10000 : 2000);
  int total = 0;
  time_t start = time (0);
  time_t now = time (0);
  for (int cellidx = 0; cellidx < edt->digital_count; cellidx++) {
    char *lockname;
    if ((lockname = lock (cellidx)) == NULL)
      continue;
    load (&sl, cellidx);
    int count = 0;
    for (int suridx = sl.surcnt; suridx < surcnt; suridx++) {
      count++;
      addsur (edt, &sl);
      //printf ("%d %d %g\r", cellidx, suridx, (double) (time (0) - now) / count);
      fflush (stdout);
    }
    if (count == 0)
      sl.cellidx = -1;
    else {
      free (ri[cellidx].r);
      ri[cellidx].r = 0;
      time_t then = now;
      now = time (0);
      time_t delta = now - then;
      total += count;
      printf ("cellidx %3d: %4ld secs, %3d surrogates, %5.3f sec/sur, %6ld total secs, %5d total surrogates, %7.5f sec/sur, g = %7.5f, %6d spikes\n",
              cellidx, delta, count, (double)delta/count, now - start, total, (double)(now - start)  / total, g[cellidx], edt->digital[cellidx]->spike_count);
      load (&sl, -1);
    }
    unlink (lockname);
    if (access ("stop", F_OK) == 0)
      break;
  }
}

void
check_surrogates (Edt *edt)
{
  for (int cellidx = 0; cellidx < edt->digital_count; cellidx++) {
    char *filename;
    if (asprintf (&filename, "%d", cellidx) == -1) exit (1);
    FILE *f = fopen (filename, "r");
    if (f == NULL) error (1, errno, "fopen");
    if (fseek (f, 4, SEEK_SET) != 0) error (1, errno, "cellidx %d fseek surcnt", cellidx);
    int surcnt;
    if (fread (&surcnt, sizeof surcnt, 1, f) != 1) error (1, errno, "cellidx %d fread surcnt", cellidx);
    printf ("cellidx %d, surcnt: %d\n", cellidx, surcnt);
    if (fseek (f, 12, SEEK_SET) != 0) error (1, errno, "cellidx %d fseek", cellidx);
    SpikeTrain st[surcnt];
    if (fread (&st, sizeof st, 1, f) != 1) error (1, errno, "cellidx %d fread st", cellidx); 
    int count = 0;
    for (int i = 0; i < surcnt; i++)
      count += st[i].N;
    if (fseek (f, count * 4, SEEK_CUR) != 0) error (1, errno, "cellidx %d fseek", cellidx);
    long end1 = ftell (f);
    if (fseek (f, 0, SEEK_END) != 0) error (1, errno, "cellidx %d fseek end", cellidx);
    long end2 = ftell (f);
    if (end1 != end2)
      error (1, 0, "cellidx %d, wrong size file\n", cellidx);
    fclose (f);
  }
}

void
get_sur_pair (Edt *edt, Pair *p, SpikeTrain *ref_st, SpikeTrain *tar_st)
{
  static SurList ref = {.cellidx = -1, .surcnt = 0};
  static SurList tar = {.cellidx = -1, .surcnt = 0};
  if (!g) {
    g = calloc (edt->digital_count, sizeof *g);
    ri = calloc (edt->digital_count, sizeof *ri);
    set_ticks_per_sec (strcmp (edt->format, "%5d%10d\n") == 0 ? 10000 : 2000);
  }
  void init (SurList *sl, int n) {
    if (n != sl->cellidx) {
      if (sl->cellidx >= 0) {
        free (ri[sl->cellidx].r);
        ri[sl->cellidx].r = 0;
      }
      load (sl, n);
    }
  }
  init (&ref, p->refidx);
  init (&tar, p->taridx);
  if (p->refidx == -1)
    return;
  ref.suridx = p->ref_suridx;
  tar.suridx = p->tar_suridx;
  if (ref.suridx == ref.surcnt)
    addsur (edt, &ref);
  else if (tar.suridx == tar.surcnt)
    addsur (edt, &tar);
  assert (ref.suridx >= 0 && ref.suridx < ref.surcnt);
  assert (tar.suridx >= 0 && tar.suridx < tar.surcnt);
  *ref_st = ref.sur[ref.suridx];
  *tar_st = tar.sur[tar.suridx];

  if (p->ref_suridx <= p->tar_suridx) {
    if (p->ref_suridx == p->tar_suridx)
      p->tar_suridx = 0;
    p->ref_suridx++;
  }
  else {
    p->tar_suridx++;
    if (p->ref_suridx == p->tar_suridx)
      p->ref_suridx = 0;
  }
}

void
edt_surrogate (Edt *edt, char *filename)
{
  Edt sur;
  memcpy (&sur, edt, sizeof sur);
  sur.digital = malloc (edt->digital_count * sizeof *sur.digital);
  for (int n = 0; n < edt->digital_count; n++) {
    sur.digital[n] = malloc (sizeof (Digital));
    sur.digital[n]->id = edt->digital[n]->id;
  }
  set_ticks_per_sec (strcmp (edt->format, "%5d%10d\n") == 0 ? 10000 : 2000);
    
  for (int n = 0; n < edt->digital_count; n++) {
    printf ("%d of %d\n", n + 1, edt->digital_count);
    SpikeTrain st = {.T = edt->digital[n]->spike, .N = edt->digital[n]->spike_count};
    double g = miura_gamma (st);
    RateFunc ri = get_ratefunc (st, g, edt->digital[n]->id);
    SpikeTrain sst = spikes_from_rate (ri, lrint (g < 1 ? 1.0 : g));
    free (ri.r);
    sur.digital[n]->spike = sst.T;
    sur.digital[n]->spike_count = sst.N;
  }
  write_edt (filename, &sur);

  for (int n = 0; n < edt->digital_count; n++) {
    free (sur.digital[n]->spike);
    free (sur.digital[n]);
  }
  free (sur.digital);
}
