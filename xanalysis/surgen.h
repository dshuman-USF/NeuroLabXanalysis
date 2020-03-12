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

#include "pairloop.h"
void get_sur_pair (Edt *edt, Pair *p, SpikeTrain *ref_st, SpikeTrain *tar_st);
void generate_surrogates (Edt *edt, int surcnt);
void check_surrogates (Edt *edt);
void edt_surrogate (Edt *edt, char *filename);
