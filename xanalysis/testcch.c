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

//gcc -O2 -Wall --std=c99 -o testcch `mpicc -showme:compile` testcch.c edt.c cch.c pairloop.c surgen.c tools.c -lgsl -lgslcblas -lrt `mpicc -showme:link`
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include "edt.h"
#include "tools.h"
#include "cch.h"
#include "pairloop.h"
#include "surgen.h"
#include "mpi.h"

int *
cch (int *ref, int *tar, int refcnt, int tarcnt, int bw)
{
  static int hist[BC];
  memset (hist, 0, sizeof hist);
  int it0 = 0;
  int d0 = -bw * BC / 2. + .5;
  for (int ir = 0; ir < refcnt; ir++) {
    for (int it = it0; it < tarcnt; it++) {
      int d = tar[it] - ref[ir];
      if (d < d0) {
        it0++;
        continue;
      }
      int bin = (d - d0) / bw;
      if (bin >= BC)
        break;
      hist[bin]++;
    }
  }
  return hist;
}

static inline uint64_t 
tick (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000000000 + ts.tv_nsec;
}

typedef struct
{
  double p;
  int N;
} P_N;

static void
get_pn (Edt *edt)
{
  FILE *f = fopen ("pn", "w");
  for (int i = 0; i < edt->digital_count; i++) {
    P_N pn = {0};
    Digital *d = edt->digital[i];
    pn.p = (double) d->spike_count / (d->spike[d->spike_count - 1] - d->spike[0]);
    pn.N = d->spike_count;
    if (fwrite (&pn, sizeof pn, 1, f) != 1) error (1, errno, "fwrite");
  }
  fclose (f);
}

static void
read_binwidths (Edt *edt, char *basename)
{
  char *filename;
  if (asprintf (&filename, "%s_binwidths", basename) == -1) exit (1);
  FILE *f = fopen (filename, "rb");
  if (f == NULL) error_at_line (1, errno, __FILE__, __LINE__, "fopen %s", filename);
  if (fread (&hc, sizeof hc, 1, f) != 1)
    error_at_line (1, errno, __FILE__, __LINE__, "fread %s", filename);
  free (filename);
  float *bwms = malloc (hc * sizeof *bwms);
  int nitems;
  if ((nitems = fread (bwms, sizeof *bwms, hc, f)) != hc)
    error_at_line (1, errno, __FILE__, __LINE__, "fread %s, hc %d, nitems %d", filename, hc, nitems);
  bw = malloc (hc * sizeof *bw);
  double tick = strcmp (edt->format, "%5d%10d\n") == 0 ? .1 : .5;
  for (int i = 0; i < hc; i++) {
    bw[i] = nearbyint (bwms[i] / tick);
    if (bw[i] < 1 || bwms[i] / tick > INT_MAX || bw[i] > INT_MAX / BC)
      error_at_line (1, 0, __FILE__, __LINE__, "bad bin width at %d: %g ms\n", i, bwms[i]);
  }
}

static void
finalize (void)
{
  MPI_Finalize ();
}

int
main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);
  atexit (finalize);
  bool include[1000];
  memset (include, 1, sizeof include);
  //  memset (include, 0, 100 * sizeof include[0]);
  char *edt_filename = "/home/roconnor/data/2003-12-11doneB_dig.edt";
  hc = 8;
  int default_bw[] = {1, 5, 15, 25, 55, 75, 105, 205};
  bw = default_bw;
  static char *basename;

  if (argc > 2 && strcmp (argv[1], "basename") == 0) {
    basename = argv[2];
    if (asprintf (&edt_filename, "%s_culled_spikefile.edt", basename) == -1) exit (1);
    argc -= 2;
    argv += 2;
  }

  if (argc > 2 && strcmp (argv[1], "sur_pair_count") == 0) {
    sur_pair_count = atoi (argv[2]);
    argc -= 2;
    argv += 2;
  }

  if (argc == 2 && strcmp (argv[1], "use_surrogate") == 0) {
    edt_filename = "/home/roconnor/data/2003-12-11doneB_dig_01.edt";
    argc = 1;
  }
  Edt *edt = read_edt (edt_filename, include);

  if (basename) {
    read_binwidths (edt, basename);
    char *dir;
    if (asprintf (&dir, "%s_surrogates", basename) == -1) exit (1);
    if (chdir (dir));
    free (dir);
  }

  if (argc > 1 && strcmp (argv[1], "pn") == 0) {
    get_pn (edt);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "edt_surrogate") == 0) {
    char *filename = "/home/roconnor/data/2003-12-11doneB_dig_01.edt";
    if (argc > 2)
      filename = argv[2];
    edt_surrogate (edt, filename);
    return 0;
  }

  if (argc > 3 && strcmp (argv[1], "bumpi") == 0) {
    int refidx = atoi (argv[2]);
    int taridx = atoi (argv[3]);
    pair_bump (edt, refidx, taridx);
    return 0;
  }

  if (argc > 3 && strcmp (argv[1], "bump") == 0) {
    int refid = atoi (argv[2]);
    int tarid = atoi (argv[3]);
    int refidx = -1, taridx = -1;
    for (int i = 0; i < edt->digital_count; i++) {
      if (edt->digital[i]->id == refid)
        refidx = i;
      if (edt->digital[i]->id == tarid)
        taridx = i;
    }
    if (refidx == -1) error (1, 0, "bad ref cell id: %s\n", argv[2]);
    if (taridx == -1) error (1, 0, "bad tar cell id: %s\n", argv[3]);
    pair_bump (edt, refidx, taridx);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "plotall") == 0) {
    int spcnt = argc > 2 ? atoi (argv[2]) : INT_MAX;
    for (int refidx = 0; refidx < edt->digital_count; refidx++)
      for (int taridx = refidx + 1; taridx < edt->digital_count; taridx++)
        pair_plot (edt, refidx, taridx, spcnt);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "plotlist") == 0) {
    int spcnt = argc > 2 ? atoi (argv[2]) : INT_MAX;
    FILE *f = fopen ("plotlist", "r");
    if (f == NULL) error (1, errno, "fopen plotlist");
    int refid, tarid;
    while (fscanf (f, "%d %d", &refid, &tarid) == 2) {
      int refidx = -1, taridx = -1;
      for (int i = 0; i < edt->digital_count; i++) {
        if (edt->digital[i]->id == refid)
          refidx = i;
        if (edt->digital[i]->id == tarid)
          taridx = i;
      }
      if (refidx == -1) error (1, 0, "bad ref cell id: %s\n", argv[2]);
      if (taridx == -1) error (1, 0, "bad tar cell id: %s\n", argv[3]);
      pair_plot (edt, refidx, taridx, spcnt);
    }
    return 0;
  }

  if (argc > 3 && strcmp (argv[1], "ploti") == 0) {
    int refidx = atoi (argv[2]);
    int taridx = atoi (argv[3]);
    int spcnt = argc > 4 ? atoi (argv[4]) : INT_MAX;
    pair_plot (edt, refidx, taridx, spcnt);
    return 0;
  }

  if (argc > 3 && strcmp (argv[1], "plot") == 0) {
    int refid = atoi (argv[2]);
    int tarid = atoi (argv[3]);
    int refidx = -1, taridx = -1;
    for (int i = 0; i < edt->digital_count; i++) {
      if (edt->digital[i]->id == refid)
        refidx = i;
      if (edt->digital[i]->id == tarid)
        taridx = i;
    }
    if (refidx == -1) error (1, 0, "bad ref cell id: %s\n", argv[2]);
    if (taridx == -1) error (1, 0, "bad tar cell id: %s\n", argv[3]);
    int spcnt = argc > 4 ? atoi (argv[4]) : INT_MAX;
    pair_plot (edt, refidx, taridx, spcnt);
    return 0;
  }
  
  if (argc > 1 && strcmp (argv[1], "list") == 0) {
    for (int i = 0; i < edt->digital_count; i++)
      printf ("%3d %3d %6d\n", i, edt->digital[i]->id, edt->digital[i]->spike_count);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "surrogates") == 0) {
    //    printf ("generating surrogates\n");
    generate_surrogates (edt, ceil (sqrt (sur_pair_count)));
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "check") == 0) {
    printf ("checking surrogates\n");
    check_surrogates (edt);
    return 0;
  }
  
  if (argc > 1 && strcmp (argv[1], "ratefuncs") == 0) {
    int start = 0;
    if (argc > 2)
      start = atoi (argv[2]);
    double *g = malloc (edt->digital_count * sizeof *g);
    RateFunc *ri = malloc (edt->digital_count * sizeof *ri);
    
    set_ticks_per_sec (strcmp (edt->format, "%5d%10d\n") == 0 ? 10000 : 2000);
    for (int n = start; n < edt->digital_count; n++) {
      SpikeTrain st = {.T = edt->digital[n]->spike, .N = edt->digital[n]->spike_count};
      printf ("getting firing rate for %d", edt->digital[n]->id); fflush (stdout);
      g[n] = miura_gamma (st);
      printf (" (gamma %g) (%d/%d)\n", g[n], n + 1, edt->digital_count);
      ri[n] = recip_isi (st, g[n]);
      char *filename;
      if (asprintf (&filename, "rf_%03d", edt->digital[n]->id) == -1) exit (1);
      FILE *f = fopen (filename, "w");
      if (f == NULL) error (1, errno, "fopen");
      if (fwrite (&ri[n].N, sizeof ri[n].N, 1, f) != 1) error (1, errno, "fwrite");
      if (fwrite (ri[n].r, sizeof ri[n].r[0], ri[n].N, f) != ri[n].N) error (1, errno, "fwrite 2");
      fclose (f);
      free (ri[n].r);
    }
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "status") == 0) {
    int threshold = sur_pair_count;
    if (argc > 2) threshold = atoi (argv[2]);
    pair_status (edt, threshold);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "count") == 0) {
    pair_count (edt);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "threshold") == 0) {
    threshold (edt);
    return 0;
  }

  if (argc > 1 && strcmp (argv[1], "time") == 0) {
    pair_time (edt);
    return 0;
  }

  if (argc > 1)
    error (1, 0, "unknown option: %s", argv[1]);

  pair_loop (edt);

  return 0;
}
