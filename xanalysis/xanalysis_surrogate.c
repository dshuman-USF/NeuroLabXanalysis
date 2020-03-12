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
#include <stdbool.h>
#include <math.h>
#include "gammafit_search.h"
#include "unpin.h"
#include "gen_control.h"
#include "xanalysis_surrogate.h"

static double g[MAX_NUM_CHAN];
static Pvoid_t t_est[MAX_NUM_CHAN];
static int *spiketime;
static bool use_smoothed_rate;

static int bdt_num;


static int
get_spikes (double *spiketimes, int mne, int ital)
{
    int n, i;
    if (spiketime == NULL)
      spiketime = malloc (mne * sizeof (double));
    for (i = n = 0; n < ital; n++) {
      spiketime[i] = floor ((spiketimes[n] * 10) + .5);
      if (i == 0 || spiketime[i] != spiketime[i - 1])
        i++;
    }
    return i;
}

static void
st_to_bdt (int code, int num, int *st, int sc)
{
  char *fname;
  if (asprintf (&fname, "tmp_%d_%d.edt", code, num) == -1) exit (1);
  FILE *f = fopen (fname, "w");
  free (fname);

  fprintf (f, "   33   3333333\n");
  fprintf (f, "   33   3333333\n");
  for (int n = 0; n < sc; n++)
    fprintf (f, "%5d%10d\n", code, st[n]);
  fclose (f);
}

static int codelist[1024];

void
ref_surrogate (double *spiketimes, int *ital, int *idp, void *ssp, int *mne, int *spikecount_array)
{
  double (*spikeshift)[*mne] = ssp;
  int ishift;
  int id = *idp - 1;
  
  printf ("%s line %d: doing id %d, ital %d, g %g\n", __FILE__, __LINE__, id, *ital, g[id]);
  
  int spikecount = get_spikes (spiketimes, *mne, *ital);
  printf ("%s line %d: ref %d spikecount: %d -> %d, g: %g\n", __FILE__, __LINE__, id, *ital, spikecount, g[id]);
  
  for (ishift = 0; ishift <= iTOTAL_SHIFTS; ishift++) {
    int *control, control_spikecount;
    //    printf ("%s line %d: doing spiketime[0] %d, id %d, spikecount %d, g %g\n", __FILE__, __LINE__, spiketime[0], id, spikecount, g[id]);
    if (use_smoothed_rate) put_t_est (t_est[id]);
    control = gen_control (spiketime, spikecount, g[id], &control_spikecount, 0, true);
    int n;
    for (n = 0; n < control_spikecount && n < *mne; n++)
      spikeshift[ishift][n] = control[n] / 10.;
    spikecount_array[ishift] = control_spikecount;
    printf ("control %d spikecount: %d\n", ishift, control_spikecount);
    st_to_bdt (100 + codelist[id], bdt_num++, control, control_spikecount);
    free (control);
  }
  if (use_smoothed_rate) clear_t_est ();
}

typedef union {Word_t w; double f;} Wrd_Flt_t;
typedef union {PWord_t w; double *f; Pvoid_t v;} P_Wrd_Flt_t;

void
plot_t_est (int *spiketimes, int spikecount, Pvoid_t t_est)
{
  if (!getenv("XANALYSIS_SURROGATE_PLOT"))
    return;

#define FNAME  "plot"
  FILE *f = fopen (FNAME, "w");
  int n;
  fprintf (f, "set format x '%%9.0f'\n"
           "plot \\\n"
           "'"FNAME"' index 0 using 1:2 with lines notitle,\\\n"
           "'"FNAME"' index 1 using 1:2 with lines notitle\n"
           "exit\n");
  
  for (n = 0; n < spikecount; n++)
    fprintf (f, "%d %d\n%d %d\n", spiketimes[n], n, spiketimes[n], n + 1);

  fprintf (f, "\n\n");


  Wrd_Flt_t Index;
  P_Wrd_Flt_t PValue;

  Index.w = 0;
  JLF(PValue.w, t_est, Index.w);
  while (PValue.w != 0) {
    fprintf (f, "%.17g %.17g\n", Index.f,  *PValue.f);
    JLN(PValue.w, t_est, Index.w);
  }

  fclose (f);
  if (system ("gnuplot -persist < " FNAME ));
  printf ("hit Enter when ready..."); fflush (stdout);
  getchar ();
  
}

void
tar_surrogate (void *spk_tm_ptr, int *mne, int *ital, int *ids,
               int *excluded, unsigned long *stp, int *surr_ital,
               double *E_begin, double *E_end, int *ecnt
               )
{
  double (*spiketimes)[*mne] = spk_tm_ptr;
  int code;
  printf ("tar_surrogate\n");
  static void *sr_tm_ptr;
  if (sr_tm_ptr == NULL)
    sr_tm_ptr = malloc (*mne * MAX_NUM_CHAN * sizeof (double));
  *stp = (unsigned long)sr_tm_ptr;
  double (*surr_times)[*mne] = sr_tm_ptr;

  use_smoothed_rate = getenv ("SURROGATE_USE_SMOOTHED_RATE");
  bool use_cth_rate = getenv ("SURROGATE_USE_CTH_RATE");
  if (use_cth_rate) use_smoothed_rate = true;
  set_use_smoothed_rate (use_smoothed_rate);
  if (use_smoothed_rate) set_preserve_t_est (true);
  for (code = 0; code < MAX_NUM_CODES; code++) {
    if (ids[code] == 0 || excluded[code])
      continue;
    int *lodif = 0;
    int *control, control_spikecount;
    int n;
    int id = ids[code] - 1;
    int spikecount = get_spikes (spiketimes[id], *mne, ital[id]);
    if (use_cth_rate) {
      printf ("calling cth_t_est, *ecnt = %d\n", *ecnt);
      t_est[id] = cth_t_est (spiketime, spikecount, E_begin, E_end, *ecnt);
      plot_t_est (spiketime, spikecount, t_est[id]);
      gammafit (spiketime, spikecount, t_est[id], &g[id]);
    }
    else {
      lodif = unpin (spiketime, spikecount);
      g[id] = gammafit_search (spiketime, lodif, spikecount);
      if (use_smoothed_rate) t_est[id] = get_t_est ();
    }
    printf ("%s line %d: doing code %d, spiketime[0] %d, id %d, spikecount %d, g %g\n", __FILE__, __LINE__, code + 1, spiketime[0], id, spikecount, g[id]);
    if (use_smoothed_rate) put_t_est (t_est[id]);
    control = gen_control (spiketime, spikecount, g[id], &control_spikecount, 0, true);
    if (use_smoothed_rate) zero_t_est ();
    for (n = 0; n < control_spikecount && n < *mne; n++)
      surr_times[id][n] = control[n] / 10.;
    surr_ital[id] = control_spikecount;
    printf ("tar %d spikecount %d -> %d, control %d\n", id, ital[id], spikecount, control_spikecount);

    codelist[id] = code;
    
    st_to_bdt (100 + code, bdt_num++, control, control_spikecount);
    
    free (control);
  }
}

