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
#include <limits.h>
#include <regex.h>
#include <X11/Xlib.h>

void screensize (int *widthp, int *heightp, int *status)
{
  Display *display;

  char *line = NULL;
  size_t alloc = 0;
  regex_t compiled;
  regcomp (&compiled, " connected.* ([0-9]+)x([0-9]+)", REG_EXTENDED);
  regmatch_t matches[3];
  FILE *f = popen ("xrandr", "r");

  int wmin = INT_MAX;
  int wmax = 0;
  int h = INT_MAX;
  while (getline (&line, &alloc, f) > 0) {
    if (regexec (&compiled, line, 3, matches, 0) == 0) {
      int v[3];
      for (int i = 1; i < 3; i++)
        v[i] = atoi (line + matches[i].rm_so);
      if (v[1] < wmin) wmin = v[1];
      if (v[1] > wmax) wmax = v[1];
      if (v[2] < h) h = v[2];
    }
  }
  int w = wmax;
  if (w == INT_MAX) w = 1024;
  if (h == INT_MAX) w = 768;
  if (h > w * 3 / 4) h = w * 3 / 4;
  if (w > h * 2) w = h * 2;
  if (w < 1024) w = 1024;
  if (h < 600) w = 600;

  *widthp = w;
  *heightp = h;
  *status = 1;
  return;

  {
  Window root;
  int x, y;
  unsigned int width, height;
  unsigned int border_width;
  unsigned int depth;

  FILE *f = popen ("xrandr | sed  -nr 's/.* connected ([0-9]+).*/\\1/p' | sort -n | head -n 1", "r");
  int w;
  if (fscanf (f, "%d", &w));
  pclose (f);
  
  if ((display = XOpenDisplay (0)) == 0) {
    *status = 0;
    return;
  }
  
  XGetGeometry (display, DefaultRootWindow (display), &root, &x, &y, &width, &height, &border_width, &depth);

  if (width >= height * 8 / 3)
    width /= 2;
  
  *widthp = w;
  *heightp = height;
  *status = 1;
  return;
  }
  
}
