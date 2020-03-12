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

char *
choose_spikefile (void)
{
  size_t alloc = 0;
  char *line;
  if (0) {
    FILE *f = popen ("cd ..; ls *.[abde]dt | wc -l", "r");
    int count;
    int items = fscanf (f, "%d", &count);
    fclose (f);
    if (items == 1 && count == 1) {
      f = popen ("cd ..; ls *.[abde]dt", "r");
      int n = getline (&line, &alloc, f);
      while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r'))
        line[n-- - 1] = 0;
      fclose (f);
      if (n > 0)
        return line;
    }
  }
  printf ("\n\nThe following input files are available:\n\n");
  if (system ("cd ..; ls *.[abde]dt"));
  int n = 0;
  {
    printf ("\n\nEnter the name of the input file used for tsfs (xassist),\n"
            "or just press Enter if you didn't run tsfs for this experiment >> ");
    fflush (stdout);
    n = getline (&line, &alloc, stdin);
    while (n > 0 && (line[n - 1] == '\n' || line[n - 1] == '\r' || line[n - 1] == ' '))
      line[n-- - 1] = 0;
  }
  return line;
}

void
copy_tsfs_files (char *bdt_file, int len)
{
  char *spikefilename;
  if (0) {
    spikefilename = strndup (bdt_file, len);
    while (len > 0 && spikefilename[len - 1] == ' ')
      len--;
    spikefilename[len] = 0;
    char *p = strrchr (spikefilename, '/');
    if (p) {
      char *s = strdup (p + 1); 
      free (spikefilename);
      spikefilename = s;
    }
  }
  else
    spikefilename = choose_spikefile ();
  char *t[] = {"-dur", "", "-ver", "-win", 0};
  for (int i = 0; t[i]; i++) {
    char *cmd;
#define F "/oberon/databases/GAIA\\ Database/tsfs%s.tab"
    if (strlen (spikefilename) > 0)
      {if (asprintf (&cmd, "> " F "; cp tsfs%s-%s.out " F, t[i], t[i], spikefilename, t[i]) == -1) exit (1);}
    else
      {if (asprintf (&cmd, "> " F , t[i]) == -1) exit (1);}
    printf ("%s\n", cmd);
    if (system (cmd));
    free (cmd);
  }
  free (spikefilename);
}
