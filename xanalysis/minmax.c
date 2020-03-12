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

//gcc -Wall -DDEBUG -o minmax minmax.c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <error.h>
#include <errno.h>

typedef struct
{
  float min_cycle;
  float max_cycle;
  float i_min;
  float i_max;
  float e_min;
  float e_max;
} MinMax;

void
read_minmax (char *db_files, MinMax *mm, int *mmok)
{
  char *filename;
  if (asprintf (&filename, "%s_minmax.txt", db_files) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    mm->max_cycle = 0;
    *mmok = 0;
  }
  else {
    *mmok = fscanf (f, "MIN_CYCLE %f MAX_CYCLE %f I_MIN %f I_MAX %f E_MIN %f E_MAX %f",
                    &mm->min_cycle,
                    &mm->max_cycle,
                    &mm->i_min,
                    &mm->i_max,
                    &mm->e_min,
                    &mm->e_max) == 6;
    fclose (f);
  }
  free (filename);
}

void
write_minmax (char *db_files, MinMax *mm)
{
  char *filename;
  if (asprintf (&filename, "%s_minmax.txt", db_files) == -1) exit (1);
  FILE *f = fopen (filename, "w");
  if (f == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't write %s", filename);
  else
    fprintf (f,
             "MIN_CYCLE %.9g\n"
             "MAX_CYCLE %.9g\n"
             "I_MIN %.9g\n"
             "I_MAX %.9g\n"
             "E_MIN %.9g\n"
             "E_MAX %.9g\n",
             mm->min_cycle,
             mm->max_cycle,
             mm->i_min,
             mm->i_max,
             mm->e_min,
             mm->e_max);
  free (filename);
  fclose (f);
}

#ifdef DEBUG
int
main (void)
{
  MinMax mm = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
  char *db_files = "test";
  write_minmax (db_files, &mm);
  int mmok;
  read_minmax  (db_files, &mm, &mmok);
  printf ("MIN_CYCLE %.9g\n"
          "MAX_CYCLE %.9g\n"
          "I_MIN %.9g\n"
          "I_MAX %.9g\n"
          "E_MIN %.9g\n"
          "E_MAX %.9g\n",
          mm.min_cycle,
          mm.max_cycle,
          mm.i_min,
          mm.i_max,
          mm.e_min,
          mm.e_max);

  return 0;
}
#endif
