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
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <error.h>
#include <errno.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "pair.h"
#include "params.h"
#include "host_dialog.h"
#include "ps_starbase.h"

void move2d (long *fildes, float *x, float *y);
void draw2d (long *fildes, float *x, float *y);
void text2d (long *fildes, float *x_in, float *y_in, char *string_in, int *xform, int *more, int slen);
void line_color (long *fildes, float *red, float *green, float *blue);
void rectangle (long *fildes, float *x1, float *y1, float *x2p, float *y2p);
void fill_color (long *fildes, float *red, float *green, float *blue);
void interior_style (long *fildes, int *style, int *edged);
void set_print (int *do_print, int *old_val);
int is_print (void);
void line_type (long *fildes, int *style);
void character_width (long *fildes, float *width);
void character_height (long *fildes, float *height);
void text_color (long *fildes, float *red, float *green, float *blue);


int sur_pair_count;             /* for pair.c */

static int cellidx[1000];
static int cellidxf[1000];
static Pair p;
static int qdtidx;
static int bwidx;
static double *min;
static double *max;
static double *mean;
static int spcnt;
static bool mirror;
static int subtract;
static bool show_control;
static bool show_cl;
static int max_straddle;
static bool *sig;
static char *surdir;

static char *
trim (char *string, int len)
{
  char *trimmed = strndup (string, len);
  while (len > 0 && trimmed[len - 1] == ' ')
    len--;
  trimmed[len] = 0;
  return trimmed;
}

/* called at the beginning of analyze_data */
/* cellidx should be equal to included - 1 */

static bool
get_ids_from_edt (int *included, char *DB_FILES, int slen)
{
  char *db_files = trim (DB_FILES, slen);
  char *filename;
  if (asprintf (&filename, "%s_culled_spikefile.edt", db_files) == -1) exit (1);
  free (db_files);
  FILE *f = fopen (filename, "r");
  if (f == NULL) {
    printf ("can't open %s for read: no surrogates\n", filename);
    free (filename);
    spcnt = 0;
    return false;
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

  memset (cellidxf, 255, sizeof cellidxf);
  int i = 0;
  for (int code = 0; code < 1000; code++)
    if (id[code] != 0)
      cellidxf[code] = i++;
  return true;
  
  for (int code = 1; code < 1000; code++)
    if (included[code - 1] != 0)
      if (cellidxf[code] != included[code - 1] - 1)
        error_at_line (1, 0, __FILE__, __LINE__, "BUG: code %d cellidx %d included %d\n", code, cellidx[code], included[code]);
  printf ("include matches edt\n");
}

static void
get_spcnt (void)
{
  FILE *f = fopen ("spcnt", "r");
  spcnt = 0;
  if (f == NULL)
    error_at_line (0, errno, __FILE__, __LINE__, "Can't open spcnt");
  else {
    if (fscanf (f, "%d", &spcnt));
    fclose (f);
  }
}

void
sc_file (int *included, char *DB_FILES, int len)
{
  if (!get_ids_from_edt (included, DB_FILES, len))
    return;
  char *db_files = trim (DB_FILES, len);
  if (asprintf (&surdir, "%s_surrogates", db_files) == -1) exit (1);
  free (db_files);
  if (chdir (surdir));
  get_spcnt ();
  if (chdir (".."));
  memset (cellidx, 255, sizeof cellidx);
  for (int code = 0; code < 1000; code++)
    if (included[code - 1] != 0)
      cellidx[code] = included[code - 1] - 1;
  for (int code = 0; code < 1000; code++)
    if (cellidx[code] != cellidxf[code])
      printf ("mismatch code %d: %d %d\n", code, cellidx[code], cellidxf[code]);
  
  if (memcmp (cellidxf, cellidx, sizeof cellidx) == 0)
    printf ("cellidx arrays match\n");
  else
    printf ("cellidx arrays don't match\n");
}

static void
free_pair (Pair *p)
{
  void freem (XtrmList *x) {
    for (int i = 0; i < p->bc; i++)
      free (x[i].lst);
    free (x);
  }
  free (p->notdone);
  free (p->straddle_at);
  free (p->sign);
  freem (p->min);
  freem (p->max);
  free (p->mean);
}

/* call at analyze_data.f:3269 */
void
sc_pair (int *ref, int *tar, int *total_num_qdts)
{
  if (!spcnt)
    return;
  int this_bc = NUM_BINS * *total_num_qdts * 4;
  static int bc;
  if (this_bc != bc) {
    free_pair (&p);
    bc = this_bc;
    min = realloc (min, bc * sizeof *min);
    max = realloc (max, bc * sizeof *max);
    mean = realloc (mean, bc * sizeof *mean);
    sig = realloc (sig, bc * sizeof *sig);
  }
  init_pair (&p, cellidx[*ref], cellidx[*tar], bc, false);
  if (chdir (surdir));
  bool read_ok = read_pair (&p, true);
  if (chdir (".."));
  if (!read_ok) {
    spcnt = 0;
    return;
  }

  max_straddle = 0;
  memset (sig, false, bc * sizeof *sig);
  for (int i = 0; i < bc; i++) {
    int getval (XtrmList *m) {
      int j = m->cnt - 1;
      while (j > 0 && m->lst[j].at > spcnt)
        j--;
      return m->lst[j].val;
    }
    min[i] = getval (&p.min[i]);
    max[i] = getval (&p.max[i]);
    if (p.straddle_at[i] == 0)
      max_straddle = INT_MAX;
    if (p.straddle_at[i] > max_straddle)
      max_straddle = p.straddle_at[i];
    if (p.straddle_at[i] == 0 || p.straddle_at[i] > spcnt)
      sig[i] = true;
  }
  mean = p.mean;
}

void
sc_qdt (int *qdt_idx_p)
{
  qdtidx = *qdt_idx_p - 1;
}

void
sc_izz (int *izz)
{
  bwidx = *izz - 1;
}

void
sc_mirror (int *mirror_p)
{
  mirror = *mirror_p;
}

void
sc_subtract (int *subtract_p)
{
  subtract = *subtract_p;
}

void
sc_control (int *control_p)
{
  if (*control_p == -1)
    show_control ^= 1;
  else
    show_control = *control_p;
}

void
sc_cl (int *cl_p)
{
  if (*cl_p == -1)
    show_cl ^= 1;
  else
    show_cl = *cl_p;
}

void
sc_plot (long *fildes, float *x0, float *y0, float *pixels_per_bin, float *scale_factor)
{
  if (spcnt == 0 || (!show_control && !show_cl))
    return;
  int left = (qdtidx * 4 + bwidx) * NUM_BINS;
  int right = left + NUM_BINS - 1;
  int start = mirror ? right : left;
  int   end = mirror ? left - 1 : right + 1;
  int   inc = mirror ? -1 : 1;

  void plot (double *hist) {
    double h[NUM_BINS];
    int j = 0;
    for (int i = start; i != end; i += inc)
      h[j++] = (hist[i] - subtract) * *scale_factor;
    float x = *x0;
    float y = *y0 + h[0];
    move2d (fildes, &x, &y);
    for (int i = 0; i < NUM_BINS - 1; i++) {
      x += *pixels_per_bin;
      draw2d (fildes, &x, &y);
      y = *y0 + h[i + 1]; 
      draw2d (fildes, &x, &y);
    }
    x += *pixels_per_bin;
    draw2d (fildes, &x, &y);
  }
  float r = 0;
  float g = 1;
  float b = 0;
  line_color (fildes, &r, &g, &b);

  int lt = DOT;
  if (is_print ())
    line_type (fildes, &lt);
  if (show_control)
    plot (mean);
  if (show_cl) {
    plot (min);
    plot (max);
  }
  lt = SOLID;
  line_type (fildes, &lt);

  r = 1; g = 0; b = 0;
  line_color (fildes, &r, &g, &b);
  fill_color (fildes, &r, &g, &b);
  int is = INT_SOLID;
  int zero = 0;
  interior_style (fildes, &is, &zero);
  bool s[NUM_BINS];
  int j = 0;
  for (int i = start; i != end; i += inc)
    s[j++] = sig[i];
  float x = *x0;
  float y1 = *y0 + 10;
  for (int i = 0; i < NUM_BINS; i++) {
    float x1 = x + *pixels_per_bin;
    if (s[i])
      rectangle (fildes, &x, y0, &x1, &y1);
    x = x1;
  }
  r = g = b = 0;
  line_color (fildes, &r, &g, &b);
  is = INT_HOLLOW;
  int one = 1;
  interior_style (fildes, &is, &one);
}

void
sc_plot_ps (long *fildes, float *x0, float *y0, float *pixels_per_bin, float *scale_factor)
{
  int old_val;
  int on = 1;
  set_print (&on, &old_val);
  sc_plot (fildes, x0, y0, pixels_per_bin, scale_factor);
  set_print (&old_val, NULL);
}

typedef struct
{
  int code;
  int time;
  int index;
  int spike;
} Event;

static int
compare_event (const void *avl_a, const void *avl_b)
{
  Event *a = (Event *) avl_a;
  Event *b = (Event *) avl_b;
  if (a->time < b->time) return -1;
  if (a->time > b->time) return  1;
  if (a->code < b->code) return -1;
  if (a->code > b->code) return  1;
  return 0;
}

static void
write_binwidths (float (*BWs)[MAX_NUM_QDTS], int *total_num_qdts, char *db_files)
{
  char *filename;
  if (asprintf (&filename, "%s_binwidths", db_files) == -1) exit (1);
  FILE *f = fopen (filename, "w");
  if (f == NULL) error (1, errno, "Can't open %s for write, aborting", filename);
  free (filename);
  int count = *total_num_qdts * 4;
  if (fwrite (&count, sizeof count, 1, f) != 1)
    error_at_line (1, errno, __FILE__, __LINE__, "fwrite");
  for (int qdt = 0; qdt < *total_num_qdts; qdt++)
    for (int i = 0; i < 4; i++)
      if (fwrite (&BWs[i][qdt], sizeof (float), 1, f) != 1) 
        error_at_line (1, errno, __FILE__, __LINE__, "fwrite");
  fclose (f);
}

static pid_t pid;

static int
system_start (const char *command)
{
  pid = fork ();
  if (pid == 0) {
#define SHELL "/bin/sh"
    printf ("path: %s\n", getenv ("PATH"));
    printf ("command: %s\n", command);
    execl (SHELL, SHELL, "-c", command, NULL);
    _exit (EXIT_FAILURE);
  }
  else if (pid < 0)
    return -1;
  return 0;
}

static int
system_finish (void)
{
  if (pid < 0)
    return -1;
  int status;
  if (waitpid (pid, &status, 0) != pid)
      status = -1;
  return status;
}

static void
start_surcontrol (char *db_files)
{
  char *cmd;
  if (asprintf (&cmd, "surcontrol.sh %s", db_files) == -1) exit (1);
  if (system_start (cmd) < 0)
    printf ("Couldn't start surrogate generation.  There will be no surrogate control.\n");
  free (cmd);
}

void
sc_finish (void)
{
  if (system_finish () == 0)
    printf ("Surrogate generation completed successfully for the surrogate control\n");
  else
    printf ("Surrogate generation returned an error.  There may be no surrogate control.\n");
}

void
sc_write_edt (void *stp, int *mne, int *IDs, int *ital, int *excluded, int *total_num_qdts, float (*BWs)[MAX_NUM_QDTS], char *DB_FILES, int len)
{
  double (*spiketimes)[*mne] = stp;
  if (!host_dialog ())
    return;
  int code[MAX_NUM_CHAN];
  memset (code, 0, sizeof code);
  int max_index = 0;
  for (int i = 0; i < MAX_NUM_CODES; i++) {
    if (IDs[i] == 0 || excluded[i])
      continue;
    if (IDs[i] < 0 || IDs[i] > MAX_NUM_CHAN)
      error_at_line (1, 0, __FILE__, __LINE__, "FATAL ERROR: IDs[%d] = %d (out of range)", i, IDs[i]);
    if (code[IDs[i] - 1] != 0)
      error_at_line (1, 0, __FILE__, __LINE__, "FATAL ERROR: codes %d and %d are both at code index %d", code[IDs[i] - 1], i + 1, IDs[i] - 1);
    code[IDs[i] - 1] = i + 1;
    if (IDs[i] - 1 > max_index)
      max_index = IDs[i] - 1;
  }
  Event event[MAX_NUM_CHAN];
  
  int j = 0;
  for (int i = 0; i <= max_index; i++) {
    if (code[i] == 0 || ital[i] <= 0)
      continue;
    event[j].index = i;
    event[j].code = code[i];
    event[j].time = nearbyint (spiketimes[i][0] * 10);
    event[j].spike = 1;
    j++;
  }
  int table_len = j;
  qsort (event, table_len, sizeof (Event), compare_event);
  char *db_files = trim (DB_FILES, len);
  write_binwidths (BWs, total_num_qdts, db_files);
  char *filename;
  if (asprintf (&filename, "%s_culled_spikefile.edt", db_files) == -1) exit (1);
  FILE *f = fopen (filename, "w");
  if (f == NULL) error_at_line (1, errno, __FILE__, __LINE__, "fopen %s", filename);
  free (filename);
  Event *evp = event;
  fprintf (f, "   33   3333333\n");
  fprintf (f, "   33   3333333\n");
  while (table_len > 0) {
    int i = evp[0].index;
    fprintf (f, "%5d%10d\n", code[i], evp[0].time);
    if (evp[0].spike < ital[i]) {
      Event e = evp[0];
      e.time = nearbyint (spiketimes[i][e.spike++] * 10);
      int j = 1;
      while (j < table_len && compare_event (&evp[j], &e) < 0) {
        evp[j - 1] = evp[j];
        j++;
      }
      evp[j - 1] = e;
    }
    else {
      evp++;
      table_len--;
    }
  }
  fclose (f);
  start_surcontrol (db_files);
  free (db_files);
}

void
sc_buttons (long *fildes)
{
  if (spcnt == 0)
    return;
  float left = 452;
  float right = 523;
  float lx3 = left + 13;
  float lx2 = left + 22;
  float tx3 = left + 23;
  float tx2 = left + 23;
  float offbot = 2;
  int at = ANNOTATION_TEXT;
  int more = 0;

  float bot = 90;
  float top = 140;
  float lbot = bot + offbot;
  float ltop = lbot + (top - bot) / 2;
  void __mod_new_draw_button_MOD_draw_button (long *, float *, float *, float *, float *, char *, float *, float *, char *, float *, float *, int, int);
  __mod_new_draw_button_MOD_draw_button (fildes, &left, &bot, &right, &top, "SUR", &lx3, &ltop, "CTL", &lx3, &lbot, 3, 3);
  float ony = top + offbot;
  float lx = show_control ? tx2 : tx3;
  text2d (fildes, &lx, &ony, show_control ? "on" : "off", &at, &more, show_control ? 2 : 3);

  bot = 5;
  top = 55;
  lbot = bot + offbot;
  ltop = lbot + (top - bot) / 2;
  __mod_new_draw_button_MOD_draw_button (fildes, &left, &bot, &right, &top, "SUR", &lx3, &ltop, "CL", &lx2, &lbot, 3, 2);

  ony = top + offbot;
  lx = show_cl ? tx2 : tx3;
  text2d (fildes, &lx, &ony, show_cl ? "on" : "off", &at, &more, show_cl ? 2 : 3);

}            

static void
warn (long *fildes, float x, float y, char *ss, char *sa, char *s2, char *s3)
{
  if (!(show_control || show_cl) || !(ss[0] == 'y' || sa[0] == 'y'))
    return;
  float h = .060;
  float w = .020;

  character_height (fildes, &h);
  character_width (fildes, &w);
  char *t = "WARNING: MIXED CONTROLS";
  float r = 1, g = 0, b = 0;
  text_color (fildes, &r, &g, &b);
  int at = ANNOTATION_TEXT;
  int more = 0;
  text2d (fildes, &x, &y, t, &at, &more, strlen (t));
  r = 0, g = 0, b = 0;
  text_color (fildes, &r, &g, &b);
}

void
sc_label (long *fildes, float *x, float *y, char *ss, char *sa, char *s2, char *s3)
{
  if (spcnt == 0)
    return;
  int at = ANNOTATION_TEXT;
  int more = 0;

  char *txt = NULL;
  if (!show_control && !show_cl)
    return;

  char *t2 = NULL;
  if (max_straddle > spcnt)
    {if (asprintf (&t2, "no straddle at %d ", spcnt) == -1) exit (1);}
  else
    {if (asprintf (&t2, "straddle at %d/%d ", max_straddle, spcnt) == -1) exit (1);}
    
  if (asprintf (&txt, "(%s%s%sdisplayed)", 
                show_control ? "surrogate control " : "", 
                show_control && show_cl ? "and " : "",
                show_cl ? t2 : "") == -1) exit (1);
  free (t2);
  text2d (fildes, x, y, txt, &at, &more, strlen (txt));
  free (txt);
  warn (fildes, 850, 760, ss, sa, s2, s3);
}

void
sc_label_ps (long *fildes, float *x, float *y, char *ss, char *sa, char *s2, char *s3)
{
  int old_val;
  int on = 1;
  set_print (&on, &old_val);
  sc_label (fildes, x, y, ss, sa, s2, s3);
  set_print (&old_val, NULL);
}

void
sc_label2 (long *fildes, float *x, float *y, char *ss, char *sa, char *s2, char *s3)
{
  if (spcnt == 0 || (!show_control && !show_cl))
    return;


  int at = ANNOTATION_TEXT;
  int more = 0;

  float xx = *x + 300;

  float y2 = *y - 20;
  char *t2 = NULL;
  if (!show_cl)
    t2 = strdup ("");
  else if (max_straddle > spcnt)
    {if (asprintf (&t2, "(no straddle at %d)", spcnt) == -1) exit (1);} 
  else
    {if (asprintf (&t2, "(straddle at %d/%d)", max_straddle, spcnt) == -1) exit (1);}
  text2d (fildes, &xx, &y2, t2, &at, &more, strlen (t2));
  free (t2);
    
  char *txt = "(surrogate_control_displayed)";
  if (show_control)
    text2d (fildes, &xx, y, txt, &at, &more, strlen (txt));
  warn (fildes, 850, 660, ss, sa, s2, s3);
}

void
sc_label2_ps (long *fildes, float *x, float *y, char *ss, char *sa, char *s2, char *s3)
{
  int old_val;
  int on = 1;
  set_print (&on, &old_val);
  sc_label2 (fildes, x, y, ss, sa, s2, s3);
  set_print (&old_val, NULL);
}

void
sc_hist (int *h, int *status)
{
  if (spcnt == 0) {
    *status = 1;
    return;
  }
  int left = (qdtidx * 4 + bwidx) * NUM_BINS;
  int right = left + NUM_BINS - 1;
  int start = mirror ? right : left;
  int   end = mirror ? left - 1 : right + 1;
  int   inc = mirror ? -1 : 1;
  int j = 0;
  for (int i = start; i != end; i += inc)
    h[j++] = nearbyint (mean[i]);
  *status = 0;
}
