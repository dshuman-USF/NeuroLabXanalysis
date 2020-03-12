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
#include <unistd.h>
#include <string.h>
#include <error.h>
#include <errno.h>

static char *
ask_for_name (void)
{
 LOOP:
  printf ("\nPlease enter the name of the response-to-perturbation file that was previously imported.\n"
          "The following response-to-perturbation data file(s) exists in this directory:\n\n");
  if (system ("ls -C *.per.db"));
  printf ("\n    PICK ONE of these ..or.. <cr> to cancel  >> ");
  fflush (stdout);
  static char *per_filename = NULL;
  static size_t alloc = 0;
  int len = getline (&per_filename, &alloc, stdin);
  if (len < 2)
    return per_filename;
  per_filename[len - 1] = 0;
  if (access (per_filename, F_OK) != 0) {
    printf ("File %s does not exist in this directory.  Please re-enter the file name.\n", per_filename);
    goto LOOP;
  }
  return per_filename;
}

static char *
get_persav_filename (char *db_files, char *per_prev_imported)
{
  if (per_prev_imported[0] != 'y')
    return strdup ("no_file");
  char *file_filename;
  if (asprintf (&file_filename, "%s_per_filename", db_files) == -1) exit (1);
  char *per_filename;
  FILE *f = fopen (file_filename, "r");

  if (f == NULL)
    per_filename = ask_for_name ();
  else {
    fseek (f, 0, SEEK_END);
    int len = ftell (f);
    rewind (f);
    per_filename = malloc (len);
    if (fread (per_filename, 1, len, f) != len) {
      error_at_line (0, errno, __FILE__, __LINE__, "Can't read %s, GAIA won't get xassist comments", per_filename);
      free (per_filename);
      fclose (f);
      return strdup ("no_file");
    }
    per_filename[len - 1] = 0;
    fclose (f);
  }
  int baselen = strlen (per_filename) - 3;
  if (baselen < 0 || strcmp (per_filename + baselen, ".db") != 0) {
      error_at_line (0, 0, __FILE__, __LINE__,
                     "\"%s\" doesn't end in .db, so I don't know how to make a per.sav filename out of it.\n"
                     "GAIA won't get xassist comments.", per_filename);
      free (per_filename);
      return strdup ("no_file");
  }
  per_filename[baselen] = 0;
  char *persav_filename;
  if (asprintf (&persav_filename, "%s.sav", per_filename) == -1) exit (1);
  free (per_filename);
  return persav_filename;
}


void
write_persav_file (char *DB_FILES, char *per_prev_imported, int len)
{
  char *dest = "/oberon/databases/GAIA\\ Database/per_comments.tab /oberon/databases/GAIA\\ Database/files_from_persav.tab";
  char *db_files = strndup (DB_FILES, len);
  while (len > 0 && db_files[len - 1] == ' ')
    len--;
  db_files[len] = 0;
  
  char *persav_filename = get_persav_filename (db_files, per_prev_imported);
  char *cmd;
  if (asprintf (&cmd, "persav_to_gaia %s %s", persav_filename, dest) == -1) exit (1);
  printf ("%s\n", cmd);
  if (system (cmd) != 0)
    error (0, 0, "Writing per comments to GAIA failed.");
  free (cmd);
  free (persav_filename);
}
  
