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

#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"
#include "cch.h"

static int *hist;
static int hist_bc;
static int hist_bw;
static int h_count;
static int *hset;

static int *
fake_hires_cch (SpikeTrain ref, SpikeTrain tar, int bw, int bc)
{
  if (hist_bc != bc)
    hist = realloc (hist, (hist_bc = bc) * sizeof *hist);
  hist_bw = bw;

  static gsl_rng *rng;
  if (rng == NULL) {
    gsl_rng_env_setup();
    rng = gsl_rng_alloc (gsl_rng_default);
  }

  double p = (double) tar.N / (tar.T[tar.N - 1] - tar.T[0]);
  for (int bin = 0; bin < bc; bin++)
    hist[bin] = gsl_ran_binomial (rng, p, ref.N);

  return hist;
}

static int *
hires_cch (SpikeTrain ref, SpikeTrain tar, int bw, int bc)
{
  if (hist_bc != bc)
    hist = realloc (hist, (hist_bc = bc) * sizeof *hist);
  hist_bw = bw;
  memset (hist, 0, bc * sizeof *hist);
  int it0 = 0;
  int d0 = -bw * bc / 2. + .5;
  for (int ir = 0; ir < ref.N; ir++) {
    for (int it = it0; it < tar.N; it++) {
      int d = tar.T[it] - ref.T[ir];
      if (d < d0) {
        it0++;
        continue;
      }
      int bin = (d - d0) / bw;
      if (bin >= bc)
        break;
      hist[bin]++;
    }
  }
  return hist;
}

static void
lumped_cch (int lump, int *lhist)
{
  assert (lump % 2 == 1);
  memset (lhist, 0, BC * sizeof *lhist);
  int bin = hist_bc / 2 - lump / 2 - lump * (BC / 2);
  for (int i = 0; i < BC; i++)
    for (int j = 0; j < lump; j++)
      lhist[i] += hist[bin++];
}

int *
cch_set (SpikeTrain ref, SpikeTrain tar, int *bw, int bw0, int bc, int hc)
{
  static bool env_checked = false;
  static bool fake = false;
  if (!env_checked) {
    char *env = getenv ("TESTCCH_CCH_TYPE");
    if (env && strcmp (env, "fake") == 0) {
      printf ("using fake cch\n");
      fake = true;
    }
    env_checked = true;
  }

  if (fake)
    fake_hires_cch (ref, tar, bw0, bc);
  else
    hires_cch (ref, tar, bw0, bc);

  if (hc != h_count)
    hset = realloc (hset, (h_count = hc) * BC * sizeof *hset);
  memset (hset, 0, hc * BC * sizeof *hset);
  int *lhist = hset;
  for (int i = 0; i < hc; i++, lhist += BC) {
    assert (bw[i] % hist_bw == 0);
    lumped_cch (bw[i] / hist_bw, lhist);
  }
  return hset;
}
