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

void
write_sta_file (char *exp_name, int len)
{
  char *experiment = strndup (exp_name, len);
  while (len > 0 && experiment[len - 1] == ' ')
    len--;
  experiment[len] = 0;
  char *cmd;
  if (asprintf (&cmd, "sta_to_gaia %s.sta > /oberon/databases/GAIA\\ Database/sta.tab", experiment) == -1) exit (1);
  printf ("%s\n", cmd);
  if (system (cmd));
  free (cmd);
  free (experiment);
}
