c Copyright 2005-2020 Kendall F. Morris

c This file is part of the Xanalysis software suite.

c     The Xanalysis software suite is free software: you can redistribute
c     it and/or modify it under the terms of the GNU General Public
c     License as published by the Free Software Foundation, either
c     version 3 of the License, or (at your option) any later version.

c     The suite is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.

c     You should have received a copy of the GNU General Public License
c     along with the suite.  If not, see <https://www.gnu.org/licenses/>.

      module mod_locate_region
      contains
      subroutine locate_region(x,y,region,screen)
*
*
*       date of last revision = 1-Mar-2004      lss
*
*       ***** This subroutine of showCCHs (which is a subroutine of x2000)      *****
*       *****   takes a look at the x and y coordinates returned from the       *****
*       *****   mouse when it is clicked or when timeout expires and then       *****
*       *****   determines which region of the screen has been selected.        *****
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                               *
*       *       REGION #        COORDINATES             BUTTON / CCH SELECTED           *
*       *                                                                               *
*       *       1               100<x<403, 550<y<750      1st CCH                       *
*       *       2               450<x<753,      "         2nd CCH                       *
*       *       3               800<x<1103,     "         3rd CCH                       *
*       *       4               1150<x<1453,    "         4th CCH                       *
*       *                                                                               *
*       *       41              125<x<226, 425<y<525      1st ACH                       *
*       *       42              275<x<576,      "         2nd ACH                       *
*       *       43              475<x<576,      "         3rd ACH                       *
*       *       44              625<x<726,      "         4th ACH                       *
*       *       45              825<x<926,      "         5th ACH                       *
*       *       46              975<x<1076,     "         6th ACH                       *
*       *       47              1175<x<1276,    "         7th ACH                       *
*       *       48              1325<x<1426,    "         8th ACH                       *
*       *                                                                               *
*       *       49              200<x<503, 175<y<375      1st CTH                       *
*       *       50              700<x<1003,     "         2nd CTH                       *
*       *       51              1175<x<1278,    "         3rd CTH                       *
*       *                                                                               *
*       *       5               5<x<205,    90<y<165      STOP, PLOT DIFFERENCE HIST.,  *
*       *                                                  PRINT (QSUM),                *
*       *                                                  RETURN TO MAIN MENU          *
*       *       205             205<x<275,      "         VIEW (QSUM) - Easter egg      *
*       *       20              1295<x<1495,    "         REDRAW                        *
*       *       6               275<x<475,      "         RESUME,ORIG.PLOT,WRITE(QSUM)  *
*       *                                                       AUTOVIEW                *
*       *       21              525<x<620,      "         +/- SINGLE SHIFT              *
*       *       121             630<x<725,      "         +/- AVERAGE SHIFT             *
*       *                                                 (21+121 = KEEP CURRENT VALUES)*
*       *       18              775<x<870,      "         +/- 2sd CONFIDENCE LIMITS     *
*       *       118             880<x<975,      "         +/- 3sd CONFIDENCE LIMITS     *
*       *       19              1025<x<1120,    "         PRINT (single CCH)            *
*       *       119             1130<x<1225,    "         WRITE (single CCH)            *
*       *       219             1120<x<1130,    "         VIEW (single CCH) - Easter egg*
      
*       *       65              800<x<850,  90<y<140      USER-DEFINED DISPLAY TIME     *
*       *       7               900<x<950,      "         SLOW, OK                      *
*       *       8               1000<x<1050,    "         MEDIUM, REPEAT                *
*       *       9               1100<x<1150,    "         FAST, CANCEL                  *
*       *       90              1300<x<1350,    "         START OVER (cycle selection)  *
*       *       10              5<x<205,     5<y<80       PREVIOUS PAIR, STATISTICS     *
*       *       11              275<x<475,      "         DIRECT ACCESS, QSUM, SHIFT 1x *
*       *       12              525<x<725,      "         PRINT (main CCH window),YES,  *
*       *                                                  MIRROR, SHIFT-CONTROL CCH,   *
*       *                                                  KEEP STORED VALUES           *
*       *       212             725<x<775,      "         VIEW (main CCH) - Easter egg  *
*       *       13              775<x<975,      "         WRITE (main CCH window),NO,   *
*       *                                                  PRINT (diff. hist.,          *
*       *                                                  shift-control),SCALE UP,     *
*       *                                                  RANGE OF BINS IN ORIGINAL CCH*
*       *       213             975<x<1025,     "         VIEW (DIFF) - Easter egg      *
*       *       14              1025<x<1225,    "         ENTER ANALYSIS RESULTS, TITLE,*
*       *                                                  WRITE (diff. hist.,          *
*       *                                                         shift-control)        *
*       *       15              1295<x<1495,    "         NEXT PAIR, RETURN,CLICKABLE   *
*       *       151             1295<x<1495,    "         NEXT SIG PAIR                 *
*       *       152             1295<x<1495,    "         NEXT FLAT PAIR                *
*       *                                                                               *
*       *       16              750<x<800,  300<y<350     OK (for Ti/Te range selection *
*       *       17                  "       400<y<450     SELECT RANGE                  *
*       *                                                                               *
*       *       101             825<x<925, 777<y<800       QDT file #1                  *
*       *       102             950<x<1050, 775<y<800      QDT file #2                  *
*       *       103             1075<x<1175, 775<y<800     QDT file #3                  *
*       *       104             1200<x<1300, 775<y<800     QDT file #4                  *
*       *       105             1325<x<1425, 775<y<800     QDT file #5                  *
*       *                                                                               *
*       *       106             5<x<90,  560<y<590        scale up CCHs on MAIN screen by 75%
*       *       107             3<x<90,  250<y<280        show/remove I and E pulses from CTHs
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       LINK with xanalysis
*
*
      real x,y
*
      integer region
      character*(*) screen
      include 'x2000parameter.defs'

      region = 0                !initialize

*
*
      if(INT(y).ge.550.and.INT(y).le.749)then !550 <= y <= 750  -->  top row of CCHs
         if(INT(x).ge.100.and.INT(x).le.403)then !100 <= x <= 403  -->  left-most CCH (least binwidth)
            region = 1
         else if(INT(x).ge.450.and.INT(x).le.753)then !2nd CCH
            region = 2
         else if(INT(x).ge.800.and.INT(x).le.1103)then !3rd CCH
            region = 3
         else if(INT(x).ge.1150.and.INT(x).le.1453)then !4th CCH
            region = 4
         else if((INT(x).ge.5.and.INT(x).le.90).and.
     +           (INT(y).ge.560.and.INT(y).le.590))then !scale up the CCHs in the main screen by 75% of min bin
            region = 106
         end if
*
      elseif(INT(y).ge.425.and.INT(y).le.525)then !425 <=y <= 525 --> row of ACHs
         if(INT(x).ge.125.and.INT(x).le.226)then
            region = 41
         elseif(INT(x).ge.275.and.INT(x).le.376)then
            region = 42
         elseif(INT(x).ge.475.and.INT(x).le.576)then
            region = 43
         elseif(INT(x).ge.625.and.INT(x).le.726)then
            region = 44
         elseif(INT(x).ge.825.and.INT(x).le.926)then
            region = 45
         elseif(INT(x).ge.975.and.INT(x).le.1076)then
            region = 46
         elseif(INT(x).ge.1175.and.INT(x).le.1276)then
            region = 47
         elseif(INT(x).ge.1325.and.INT(x).le.1426)then
            region = 48
         end if

      elseif(INT(y).ge.250.and.INT(y).le.280)then
         if(INT(x).ge.3.and.INT(x).le.90)region=107
*
*        case (175:375)                 !175 <=y <= 375 --> row of CTHs
*         select case (INT(x))          !now look at the x coordinate
*          case (200:503)               !1st CTH
*           region = 49
*          case (700:1003)              !2nd CTH
*           region = 50
*          case (1175:1278)             !3rd CTH
*           region = 51
*          end select
*
      else if(INT(y).ge.90.and.INT(y).le.165)then !90 <= y <= 165
         if(.false.)then
         else if(INT(x).ge.5.and.INT(x).le.205)then !STOP button
            region = 5
         else if(INT(x).gt.205.and.INT(x).lt.275)then !VIEW QSUM (Easter egg)
            region = 205
         else if(INT(x).ge.275.and.INT(x).le.450)then !RESUME, AUTOVIEW, ORIG.PLOT, WRITE (QSUM) button
            region = 6
         else if(INT(x).ge.452.and.INT(x).le.523)then !SURROGATE CONTROL
            region = 206
         else if(INT(x).ge.525.and.INT(x).le.725)then !+/- SHIFT buttons
            if(.false.)then
            else if(INT(x).ge.525.and.INT(x).le.620)then
               region = 21      !+/- SINGLE SHIFT
            else if(INT(x).ge.630.and.INT(x).le.725)then
               region = 121     !+/- AVERAGE SHIFT
            else if(INT(x).gt.620.and.INT(x).lt.630)then
               region = 221     !middle channel of KEEP CURRENT VALUES
            end if
         else if(INT(x).ge.775.and.INT(x).le.870)then !+/- 2 s.d. CONFIDENCE LIMITS button for single CCH
            region = 18
         else if(INT(x).ge.880.and.INT(x).le.975)then !+/- 3 s.d. CONFIDENCE LIMITS button for single CCH
            region = 118
         else if(INT(x).ge.1025.and.INT(x).le.1225)then !PRINT and WRITE buttons for single CCH
            if(.false.)then
            else if(INT(x).ge.1025.and.INT(x).le.1120)then
               region = 19      !PRINT
            else if(INT(x).ge.1130.and.INT(x).le.1225)then
               region = 119     !WRITE
            else if(INT(x).gt.1120.and.INT(x).lt.1130)then
               region = 219     !VIEW (Easter egg)
            end if
         else if(INT(x).ge.1295.and.INT(x).le.1495)then !REDRAW button
            region = 20
         end if
*
      else if(INT(y).ge.5.and.INT(y).le.80)then !5 <= y<= 80  --> bottom row of buttons
         if(.false.)then
         else if(INT(x).ge.5.and.INT(x).le.205)then !PREVIOUS PAIR or STATISTICS button
            region = 10
            if((screen.eq.'main').and.(INT(x).ge.5.and.INT(x).le.100)) !FIRST PAIR button
     +           region = 1010
         else if(INT(x).ge.275.and.INT(x).le.450)then !DIRECT ACCESS or QSUM button
            region = 11
         else if(INT(x).ge.452.and.INT(x).le.523)then !SURROGATE CONTROL CL
            region = 211
         else if(INT(x).ge.525.and.INT(x).le.725)then !RETURN TO MAIN MENU button
            region = 12
         else if(INT(x).ge.775.and.INT(x).le.975)then !STATISTICAL ANALYSIS button
            region = 13
            if(screen.eq.'single'.and.(INT(x).ge.5.and.INT(x).le.873))
     +           region = 131   !SCALE UP by subtracting 75% of minimum bin
         else if(INT(x).ge.1025.and.INT(x).le.1225)then !ENTER ANALYSIS RESULTS button
            region = 14
         else if(INT(x).ge.1295.and.INT(x).le.1495)then !NEXT  or RETURN or CLICKABLE button
            region = 15
            if((screen.eq.'main').and.
     +           (INT(x).ge.1377.and.INT(x).le.1435))region=151 !NEXT SIG FEATURES ONLY button
            if((screen.eq.'main').and.
     +           (INT(x).ge.1437.and.INT(x).le.1495))region=152 !NEXT FLAT ONLY button
         else if(INT(x).gt.725.and.INT(x).lt.775)then !VIEW main CCH (btw PRINT and WRITE buttons-Easter egg)
            region = 212
         else if(INT(x).gt.975.and.INT(x).lt.1025)then !VIEW DIFF (btw PRINT and WRITE buttons - Easter egg)
            region = 213
         end if
*
      else if(INT(y).ge.300.and.INT(y).le.350)then !300 <= y <= 350 --> OK button for range selection
         if(.false.)then
         else if(INT(x).ge.750.and.INT(x).le.800)then
            region = 16
         end if
*
      else if(INT(y).ge.400.and.INT(y).le.450)then !400 <= y <= 450 --> SELECT RANGE button 
         if(.false.)then
         else if(INT(x).ge.750.and.INT(x).le.800)then
            region = 17
         end if
*
      else if(INT(y).ge.775.and.INT(y).le.800)then
         if(.false.)then
         else if(INT(x).ge.825.and.INT(x).le.925)then
            region = 101        !QDT file #1
         else if(INT(x).ge.950.and.INT(x).le.1050)then
            region = 102        !QDT file #2
         else if(INT(x).ge.1075.and.INT(x).le.1175)then
            region = 103        !QDT file #3
         else if(INT(x).ge.1200.and.INT(x).le.1300)then
            region = 104        !QDT file #4
         else if(INT(x).ge.1325.and.INT(x).le.1425)then
            region = 105        !QDT file #5
         end if
*
      else
         return
*
      end if
*
*
*
*
      return
      end
        
      end module mod_locate_region
