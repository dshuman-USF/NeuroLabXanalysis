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

#define MAX_NUM_CHAN 200
#define MAX_NUM_CODES 999
#define iTOTAL_SHIFTS 20

void ref_surrogate (double *spiketimes, int *ital, int *idp, void *ssp, int *mne, int *spikecount_array);

void
tar_surrogate (void *spk_tm_ptr, int *mne, int *ital, int *ids,
               int *excluded, unsigned long *stp, int *surr_ital,
               double *E_begin, double *E_end, int *ecnt);

