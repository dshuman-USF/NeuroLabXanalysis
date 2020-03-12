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

void calc_cth_cch (double SPIKETIMES[MAX_NUM_CHAN][MAX_NUM_EVENTS], int ital[MAX_NUM_CHAN], int *refp, int *tarp, int *eplp, double *f);
void free_cth_cch ();
void rebin_cth_cch (float *binwp, int ihist[NUM_BINS]);
void rebin_var (float *binwp, int ihist[NUM_BINS]);
