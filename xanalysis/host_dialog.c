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

//gcc -DSTANDALONE -Wall --std=c99 -o host_dialog host_dialog.c
#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <error.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>

static unsigned
next_in_list (char *p)
{
  static unsigned a = 1, b = 0;
  static char *bp, *bp0;
  if (p != NULL)
    bp0 = bp = p;

  if (a <= b)
    return a++;

  if (*bp != '\0' && *bp != '\n') {
    if (!isdigit(*bp))
      error_at_line (1, 0, __FILE__, __LINE__, "bad list format: %s", bp0);
    b = a = strtoul (bp, (char **)&bp, 10);
    if (*bp == '-') {
      bp++;
      if (!isdigit(*bp))
        error_at_line (1, 0, __FILE__, __LINE__, "bad list format: %s", bp0);
      b = strtoul (bp, (char **)&bp, 10);
    }
    if (!(a <= b))
      error_at_line (1, 0, __FILE__, __LINE__, "bad list format: %d, %d", a, b);
    if (b >= UINT_MAX)
      error_at_line (1, 0, __FILE__, __LINE__, "bad list format: %d, %d", a, b);
    if (*bp == ',')
      bp++;
    return a++;
  }
  return UINT_MAX;
}

static int
getcpuint (int i, char *s)
{
  char *path = "/sys/devices/system/cpu/cpu%d/topology/%s";
  char *filename;
  if (asprintf (&filename, path, i, s) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    error (1, 0, "Can't open %s, assuming 1 processor", filename);
    free (filename);
    return -1;
  }
  int val;
  if (fscanf (f, "%d", &val) != 1) {
    error (1, errno, "Can't read %s, assuming 1 processor", filename);
    val = -1;
  }
  free (filename);
  fclose (f);
  return val;
}

static int
num_cores (void)
{
  static char *line = NULL;
  static size_t len = 0;
  char *filename = "/sys/devices/system/cpu/online";
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    error (1, 0, "Can't open %s, assuming 1 processor", filename);
    return 1;
  }
  if (getline (&line, &len, f) < 1) {
    error (1, errno, "Can't read %s, assuming 1 processor", filename);
    fclose (f);
    return 1;
  }
  fclose (f);
  static int alloc = 0;
  int cnt = 0;
  static struct {unsigned coreid; unsigned pkgid;} *list = NULL;
  for (int cpu = next_in_list (line); cpu < UINT_MAX; cpu = next_in_list (0)) {
    int coreid, pkgid;
    if ((coreid = getcpuint (cpu, "core_id")) < 0) return 1;
    if ((pkgid  = getcpuint (cpu, "physical_package_id")) < 0) return 1;
    int i;
    for (i = 0; i < cnt; i++)
      if (coreid == list[i].coreid && pkgid == list[i].pkgid)
        break;
    if (i < cnt)
      continue;
    if (cnt + 1 > alloc)
      list = realloc (list, (alloc = cnt + 1) * sizeof *list);
    list[cnt].coreid = coreid;
    list[cnt++].pkgid = pkgid;
  }
  return cnt;
}

static bool
check_host (char *hostname)
{
  char *cmd;
  if (asprintf (&cmd, "ssh -o 'BatchMode yes' -o 'ConnectTimeout 1' %s true", hostname) == -1) exit (1);
  if (system (cmd) != 0) {
    fprintf(stderr,"\nCould not connect with host %s\n",hostname);
    return false;
  }
  free (cmd);
  if (asprintf (&cmd, "mpirun --mca pml ob1 -H %s -np 1 testmpi", hostname) == -1) exit (1);
  FILE *f = popen (cmd, "r");
  free (cmd);
  char *line = NULL;
  size_t len = 0;
  if (getline (&line, &len, f) != 3)
    goto OUT;
  if (strcmp (line, "OK\n") != 0)
    goto OUT;
  if (getline (&line, &len, f) != -1)
    goto OUT;
  free (line);
  if (pclose (f) != 0)
    return false;
  return true;
 OUT:
  pclose (f);
  free (line);
  return false;
}

static int
ask_slots (char *hostname, int maxslots)
{
  static bool hdr_done = false;
  if (hostname == NULL) {
    hdr_done = false;
    return 0;
  }
  if (!hdr_done) {
    printf ("\n\tThe surrogate generation for the surrogate control will be done on multiple processors.\n"
            "\tYou must specify how many processors to use on each system.\n"
            "\tHit 'Enter' for each to use the max.\n"
            "\tEnter 0 for each one to skip the surrogate control\n\n");
    hdr_done = true;
  }
  int slots = -1;
  while (slots < 0 || slots > maxslots) {
    printf ("\tHow many processors on %s (0 to %d)?  >> ", hostname, maxslots);
    static char *line = NULL;
    static size_t len = 0;
    int chars_read;
    if ((chars_read = getline (&line, &len, stdin)) < 1)
      continue;
    if (chars_read == 1 && line[0] == '\n')
      slots = maxslots;
    else
      sscanf (line, "%d", &slots);
  }
  return slots;
}

static int
get_slots (char *hostname, int maxslots)
{
  static FILE *f = NULL;
  if (hostname == NULL) {
    if (f) fclose (f);
    f = NULL;
    return 0;
  }
  if (!check_host (hostname))
    return 0;
  int slots = ask_slots (hostname, maxslots);
  if (slots == 0)
    return 0;
  if (f == NULL) {
    char *filename = "hostfile";
    f = fopen (filename, "w");
    if (f == NULL) {
      error_at_line (0, errno, __FILE__, __LINE__, "Can't open %s for write.  THERE WILL BE NO SURROGATE CONTROL.", filename);
      return 0;
    }
  }
  fprintf (f, "%s slots=%d\n", hostname, slots);
  return slots;
}

static void
setup_passwordless (FILE *f)
{
  char *line = NULL;
  size_t len = 0;
  char *cmd;
  if (asprintf (&cmd, "setup_passwordless_login") == -1) exit (1);
  bool gotone = false;
  while (getline (&line, &len, f) != -1) {
    char *hostname = NULL;
    int maxslots;
    if (sscanf (line, "%ms slots=%d", &hostname, &maxslots) == 2) {
      gotone = true;
      char *tmp;
      if (asprintf (&tmp, "%s %s", cmd, hostname) == -1) exit (1);
      free (cmd);
      cmd = tmp;
    }
    free (hostname);
  }
  free (line);
  rewind (f);
  if (gotone) {
    if (system (cmd));
    free (cmd);
  }
  else {
    free (cmd);
    if (system ("setup_passwordless_login localhost"));
  }
}

bool
host_dialog (void)
{
  ask_slots (NULL, 0);
  char *filename = "/etc/openmpi/openmpi-default-hostfile";
  FILE *f = fopen (filename, "r");
  int total = 0;
  int gotone = false;
  if (f == NULL)
    error_at_line (0, errno, __FILE__, __LINE__, "Can't open %s for read", filename);
  else {
    setup_passwordless (f);
    char *line = NULL;
    size_t len = 0;
    while (getline (&line, &len, f) != -1) {
      char *hostname = NULL;
      int maxslots;
      if (len <= 0)
        continue;
      if (sscanf (line + (line[0] == '#'),
                  "%ms slots=%d", &hostname, &maxslots) == 2) {
        printf("found: %s %d\n",hostname,maxslots);
        total += get_slots (hostname, maxslots);
        gotone = true;
      }
      free (hostname);
    }
    free (line);
  }
  if (!gotone) {
    int maxslots = num_cores ();
    printf("Using localhost\n");
    total += get_slots ("localhost", maxslots);
  }
  get_slots (NULL, 0);
  fflush(stdout);
  return total > 0;
}

#ifdef STANDALONE
int
main (void)
{
  if (host_dialog ())
    printf ("OK\n");
  else
    printf ("no good\n");
  return 0;
}
#endif
