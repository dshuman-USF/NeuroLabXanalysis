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



//gcc -Wall --std=c99 -o edt_from_surr edt_from_surr.c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <error.h>
#include <errno.h>

typedef struct
{
  int id;
  int spike_count;
  int spike_alloc;
  int *spike;
} Digital;

typedef struct
{
  int code;
  int time;
} Analog;

typedef struct
{
  int digital_count;
  Digital **digital;
  int analog_count;
  int analog_alloc;
  Analog *analog;
  char *header;
  char *format;
  int linelen;
  int gdt21time;
  int gdt22time;
  bool is_gdt;
  bool is_gdf;
  bool has_header;
} Edt;

typedef struct
{
  int code;
  int time;
} CodeTime;

typedef struct
{
  int N;
  int *T;
} SpikeTrain;

typedef struct 
{
  int cellidx;
  int surcnt;
  int suridx;
  bool changed;
  SpikeTrain *sur;
} SurList;

static int id2idx[1000];
static int cellcnt;

static Edt *
read_edt (char *filename)
{
  size_t read;
  char * line = NULL;
  size_t len = 0;
  FILE *f;
  Edt *edt;
  Digital **digital;

  if ((f = fopen (filename, "r")) == NULL) {
    fprintf (stderr, "can't open %s for read: %s\n", filename, strerror (errno));
    exit (1);
  }
  
  if (getline (&line, &len, f));

  if ((read = getline (&line, &len, f)) == -1) {
    fprintf (stderr, "Empty input file %s, aborting\n", filename);
    exit (1);
  }
  digital = calloc (1000, sizeof *digital);
  edt = calloc (1, sizeof (Edt));
  {
    int n = strlen (line);
    int c;
    while (n > 0 && ((c = line[n - 1]) == '\r' || c == '\n'))
      n--;
    edt->linelen = n;
    line[n] = 0;
    edt->header = strdup (line);
    if (asprintf (&edt->format, "%%5d%%%dd\n", n - 5) == -1) exit (1);
  }
  if (edt->linelen != 13 && edt->linelen != 15)
    error_at_line (1, 0, __FILE__, __LINE__, "%s has %d character lines, expected 13 or 15: unknown format, aborting", filename, edt->linelen);
  edt->digital_count = 0;
  while ((read = getline (&line, &len, f)) != -1) {
    int code, time, id;
    time = atoi (line + 5);
    line[5] = 0;
    code = atoi (line);
    id = code > 4095 ? 1000 + code / 4096 : code;
    if (id < 1000) {
      if (digital[id] == NULL) {
        digital[id] = calloc (1, sizeof *digital);
        digital[id]->id = id;
        edt->digital_count++;
      }
      if (digital[id]->spike_count == digital[id]->spike_alloc)
        digital[id]->spike = realloc (digital[id]->spike, (digital[id]->spike_alloc += 100000) * sizeof (int));
      digital[id]->spike[digital[id]->spike_count++] = time;
    }
  }
  edt->digital = malloc (edt->digital_count * sizeof (Digital));
  int i = 0;
  for (int n = 0; n < 1000; n++)
    if (digital[n] != NULL)
      edt->digital[i++] = digital[n];
  free (digital);
  assert (i = edt->digital_count);
  fclose (f);
  edt->analog_count = 0;
  return edt;
}

static int
compare_time (const void *a, const void *b)
{
  const CodeTime *da = (const CodeTime *) a;
  const CodeTime *db = (const CodeTime *) b;

  return (da->time > db->time) - (da->time < db->time);
}


static void
write_edt (char *filename, Edt *edt)
{
  FILE *f;
  int count = edt->analog_count, n, spike_idx, id_idx;
  CodeTime *s;

  if ((f = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Can't open %s for write, aborting: %s\n", filename, strerror (errno));
    exit (1);
  }
  for (n = 0; n < edt->digital_count; n++)
    count += edt->digital[n]->spike_count;
  s = malloc (count * sizeof *s);
  for (n = id_idx = 0; id_idx < edt->digital_count; id_idx++)
    for (spike_idx = 0; spike_idx < edt->digital[id_idx]->spike_count; spike_idx++) {
      s[n].code = edt->digital[id_idx]->id;
      s[n].time = edt->digital[id_idx]->spike[spike_idx];
      n++;
    }
  qsort (s, count, sizeof (CodeTime), compare_time);
  fprintf (f, "%s\n", edt->header);
  fprintf (f, "%s\n", edt->header);
  for (n = 0; n < count; n++)
    fprintf (f, edt->format, s[n].code, s[n].time);
  free (s);
  fclose (f);
}

static void
get_ids_from_culled_edt (char *prefix)
{
  char *filename;
  if (asprintf (&filename, "x2000/%s_culled_spikefile.edt", prefix) == -1) exit (1);
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    printf ("can't open %s for read: no surrogates\n", filename);
    free (filename);
    exit (1);
  }
  free (filename);
  
  char * line = NULL;
  size_t len = 0;
  if (getline (&line, &len, f));
  if (getline (&line, &len, f));

  int id[1000];
  memset (id, 0, sizeof id);
  while (getline (&line, &len, f) >= 5) {
    line[5] = 0;
    int code = atoi (line);
    if (code < 1000)
      id[code] = code;
  }
  fclose (f);

  memset (id2idx, 255, sizeof id2idx);
  int i = 0;
  for (int code = 0; code < 1000; code++)
    if (id[code] != 0)
      id2idx[code] = i++;
  cellcnt = i;
}

static SurList *
get_surrogates (char *prefix, int start, int surcnt)
{
  SurList *slary = malloc (cellcnt * sizeof *slary);
  SurList *sl;
  for (int cellidx = 0; cellidx < cellcnt; cellidx++) {
    sl = &slary[cellidx];
    char *filename;
    if (asprintf (&filename, "x2000/%s_surrogates/%d", prefix, cellidx) == -1) exit (1);
    FILE *f = fopen (filename, "r");
    if (f == NULL) error_at_line (1, errno, __FILE__, __LINE__, "fopen %s", filename);
    if (fread (sl, sizeof (int), 3, f) != 3) error_at_line (1, errno, __FILE__, __LINE__, "fread %s", filename);
    if (sl->surcnt < start + surcnt) error_at_line (1, 0, __FILE__, __LINE__, "There are only %d surrogates in %s, aborting.", sl->surcnt, filename);
    sl->sur = malloc (sl->surcnt * sizeof (SpikeTrain));
    if (fread (sl->sur, sizeof (SpikeTrain), sl->surcnt, f) != sl->surcnt) error_at_line (1, errno, __FILE__, __LINE__, "fread SpikeTrains from %s", filename);
    for (int i = 0; i < start; i++)
       fseek (f, sizeof (int) * sl->sur[i].N, SEEK_CUR);
    for (int i = start; i < start + surcnt; i++) {
      sl->sur[i].T = malloc (sl->sur[i].N * sizeof (int));
      if (fread (sl->sur[i].T, sizeof (int), sl->sur[i].N, f) != sl->sur[i].N)
        error_at_line (1, errno, __FILE__, __LINE__, "fread sur %d from %s", i, filename);
    }
    sl->surcnt = start + surcnt;
    free (filename);
  }
  return slary;
}

static void
gen_edt (Edt *edt, char *prefix, int suridx, SurList *slary)
{
  for (int didx = 0; didx < edt->digital_count; didx++) {
    int id = edt->digital[didx]->id;
    int cellidx = id2idx[id];
    if (cellidx >= 0) {
      edt->digital[didx]->spike_count = slary[cellidx].sur[suridx].N;
      edt->digital[didx]->spike       = slary[cellidx].sur[suridx].T;
      if (edt->linelen == 13)
        for (int i = 0; i < edt->digital[didx]->spike_count; i++)
          edt->digital[didx]->spike[i] /= 5;
    }
  }
  char *filename;
  if (asprintf (&filename, "%s_surrogate_%d.%s", prefix, suridx, edt->linelen == 13 ? "bdt" : "edt") == -1) exit (1);
  write_edt (filename, edt);
  free (filename);
}

static void
list_preserved (Edt *edt)
{
  printf ("These codes were preserved (not replaced with surrogates):\n");
  for (int didx = 0; didx < edt->digital_count; didx++) {
    int id = edt->digital[didx]->id;
    if (id2idx[id] < 0)
      printf (" %d", id);
  }
  printf ("\n");
}

int
main (int argc, char **argv)
{
  if (argc < 3) {
    printf ("\n"
            "  usage: %s SPIKEFILE PREFIX [COUNT [START]]\n"
            "  \n"
            "  Generates COUNT copies of SPIKEFILE with spike channels replaced with\n"
            "  surrogates from the x2000/PREFIX_surrogates directory, and without\n"
            "  analog channels.\n"
            "  \n"
            "  Only the channels that appear in x2000/PREFIX_culled_spikefile.edt are\n"
            "  replaced.  The other spike channels are copied.  The program will list\n"
            "  the channels that were not replaced.\n"
            "  \n"
            "  This program must be run in the directory above x2000.\n"
            "  \n"
            "  If COUNT and START are not specified, one file is generated, named\n"
            "  PREFIX_surrogate_0.bdt (or edt — it will match SPIKEFILE), in the\n"
            "  directory above x2000.\n"
            "  \n"
            "  If COUNT is specified, that number of copies are generated, each with\n"
            "  different surrogates, and with the number in the name incremented for\n"
            "  each copy.\n"
            "  \n"
            "  If START is also specifed, the first copy will have that number in the\n"
            "  name, and the rest will increment from there.  A copy with a certain\n"
            "  name will have the same surrogates in it regardless of the COUNT and\n"
            "  START arguments used to generate it.\n"
            "  \n"
            "  The surrogates only have spikes during cycles that were used in the\n"
            "  cross-correlation histograms.\n"
            "\n"
            , argv[0]);
    exit (0);
  }
  Edt *edt = read_edt (argv[1]);
  char *prefix = argv[2];
  get_ids_from_culled_edt (prefix);

  int surcnt = 1; if (argc > 3) surcnt = atoi (argv[3]);
  int start  = 0; if (argc > 4) start  = atoi (argv[4]);

  SurList *slary = get_surrogates (prefix, start, surcnt);
  for (int suridx = start; suridx < start + surcnt; suridx++)
    gen_edt (edt, prefix, suridx, slary);
  list_preserved (edt);
  
  return 0;
}
