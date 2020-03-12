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

#include <stdio.h>
#include <string.h>

#include "params.h"

static void
doom_msg (int doomed[MAX_NUM_CHAN], int doomcnt)
{
  if (doomcnt <= 0)
    return;
  printf (" (NOTE: you may not include");
  for (int i = 0; i < doomcnt; i++)
    printf (" %d", doomed[i]);
  printf ("\n  because there are no spikes during acceptable cycles.)\n");
}

void
doomed_cells (void *stp, int *mne, int ital[MAX_NUM_CHAN], int ids[MAX_NUM_CODES],
              float *maxduration, int excluded[MAX_NUM_CODES], int *e_pulse, int *i_pulse, int doomed[MAX_NUM_CHAN])
{
  double (*spiketimes)[*mne] = stp;
  
  int echan = ids[*e_pulse - 1] - 1;
  int ichan = ids[*i_pulse - 1] - 1;

  memset (doomed, 0, MAX_NUM_CHAN * sizeof doomed[0]);
  int doomcnt = 0;
  for (int code = 1; code < MAX_NUM_CODES; code++) {
    if (ids[code - 1] == 0 || excluded[code - 1] == 1)
      continue;
    int chan = ids[code - 1] - 1;
    int iidx = 0;
    int cidx = 0;
    int eidx;
    for (eidx = 1; eidx < ital[echan]; eidx++) {
      int eidx0 = eidx - 1;
      if (spiketimes[echan][eidx] - spiketimes[echan][eidx0] > *maxduration)
        continue;               /* cycle too long */
      while (iidx < ital[ichan] && spiketimes[ichan][iidx] <= spiketimes[echan][eidx0])
        iidx++;
      if (iidx == ital[ichan] || spiketimes[ichan][iidx] >= spiketimes[echan][eidx])
        continue;               /* no I pulse */
      if (iidx + 1 < ital[ichan] && spiketimes[ichan][iidx + 1] < spiketimes[echan][eidx])
        continue;               /* multiple I pulses */
      while (cidx < ital[chan] && spiketimes[chan][cidx] < spiketimes[echan][eidx0])
        cidx++;
      if (cidx == ital[chan] || spiketimes[chan][cidx] >= spiketimes[echan][eidx])
        continue;               /* no spikes */
      break;
    }
    if (eidx == ital[echan])
      doomed[doomcnt++] = code;
  }
  doom_msg (doomed, doomcnt);
}

void
exclude_doomed (int excluded[MAX_NUM_CODES], int doomed[MAX_NUM_CHAN])
{
  for (int i = 0; i < MAX_NUM_CHAN; i++) {
    if (!doomed[i]) break;
    excluded[doomed[i] - 1] = 1;
  }
}
