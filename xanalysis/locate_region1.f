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

      module mod_locate_region1
      contains
      subroutine locate_region1(x,y,region)
*
*
*       date of last revision = 25-apr-1999     lss
*
*       ***** This subroutine of showCCHs (which is a subroutine of NOMAD)      *****
*       *****   takes a look at the x and y coordinates returned from the       *****
*       *****   mouse when it is clicked or when timeout expires and then       *****
*       *****   determines which region of the screen has been selected.        *****
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                               *
*       *       REGION #        COORDINATES             BUTTON / CCH SELECTED           *
*       *                                                                               *
*       *       16              75<x<90, 90<y<105       PRIMARY FEATURE: PEAK           *
*       *       17              "        70<y<85                        TROUGH          *
*       *       18              "        50<y<65                        MP & T          *
*       *       19              "        30<y<45                        SKIP            *
*       *       20              400<x<415, 90<y<105     LOCATION:       CENTRAL         *
*       *       21              "          70<y<85                      OFFSET RIGHT    *
*       *       22              "          50<y<65                      OFFSET LEFT     *
*       *       23              "          30<y<45                      SKIP            *
*       *       24              800<x<815, 90<y<105     SECONDARY FEATURE: NONE         *
*       *       25              "          70<y<85                      PEAK            *
*       *       26              "          50<y<65                      TROUGH          *
*       *       27              "          30<y<45                      MP & T          *
*       *       28              1150<x<1165, 90<y<105   LOCATION:       CENTRAL         *
*       *       29              "            70<y<85                    OFFSET RIGHT    *
*       *       30              "            50<y<65                    OFFSET LEFT     *
*       *       31              "            30<y<45                    SKIP            *
*       *       32              5<x<205,    130<y<165                   FLAT            *
*       *       33              275<x<475,      "                       SKIP ALL        *
*       *       34              500<x<600, 5<y<25                       OK              *
*       *       35              650<x<750,    "                         CLEAR           *
*       *       36              800<x<900,    "                         COMMENTS        *
*       *                                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       link with nomad
*
*
      real x,y
*
      integer region
      region = 0                !initialize
*
*
      if(.false.)then           !look at y-coordinate first
*                       !not a "valid" y-coordinate
      else if(INT(y).ge.5.and.INT(y).le.25)then
         if(.false.)then
         else if(INT(x).ge.500.and.INT(x).le.600)then !OK
            region = 34
         else if(INT(x).ge.650.and.INT(x).le.750)then !CLEAR
            region = 35
         else if(INT(x).ge.800.and.INT(x).le.900)then !COMMENT
            region = 36
         end if
* 
      else if(INT(y).ge.130.and.INT(y).le.165)then !90 <= y <= 165
         if(.false.)then
         else if(INT(x).ge.5.and.INT(x).le.205)then !FLAT button
            region = 32
         else if(INT(x).ge.275.and.INT(x).le.475)then !SKIP ALL button
            region = 33
         end if
*
      else if(INT(y).ge.90.and.INT(y).le.105)then !90 <= y<= 105  
         if(.false.)then
         else if(INT(x).ge.75.and.INT(x).le.210)then !PRIMARY FEATURE = PEAK
            region = 16
         else if(INT(x).ge.400.and.INT(x).le.535)then !LOCATION = CENTRAL
            region = 20
         else if(INT(x).ge.800.and.INT(x).le.935)then !SECONDARY FEATURE = NONE
            region = 24
         else if(INT(x).ge.1150.and.INT(x).le.1295)then !LOCATION = CENTRAL
            region = 28
         end if
*
      else if(INT(y).ge.70.and.INT(y).le.85)then ! 70 <= y<= 85  
         if(.false.)then
         else if(INT(x).ge.75.and.INT(x).le.210)then !PRIMARY FEATURE = TROUGH
            region = 17
         else if(INT(x).ge.400.and.INT(x).le.535)then !LOCATION = OFFSET RIGHT
            region = 21
         else if(INT(x).ge.800.and.INT(x).le.935)then !SECONDARY FEATURE = PEAK
            region = 25
         else if(INT(x).ge.1150.and.INT(x).le.1295)then !LOCATION = OFFSET RIGHT
            region = 29
         end if
*
      else if(INT(y).ge.50.and.INT(y).le.65)then ! 50 <= y<= 65  
         if(.false.)then
         else if(INT(x).ge.75.and.INT(x).le.210)then !PRIMARY FEATURE = MULT PEAKS & TROUGHS
            region = 18
         else if(INT(x).ge.400.and.INT(x).le.535)then !LOCATION = OFFSET LEFT
            region = 22
         else if(INT(x).ge.800.and.INT(x).le.935)then !SECONDARY FEATURE = TROUGH
            region = 26
         else if(INT(x).ge.1150.and.INT(x).le.1295)then !LOCATION = OFFSET LEFT
            region = 30
         end if
*
      else if(INT(y).ge.30.and.INT(y).le.45)then ! 30 <= y<= 45  
         if(.false.)then
         else if(INT(x).ge.75.and.INT(x).le.210)then !PRIMARY FEATURE = SKIP
            region = 19
         else if(INT(x).ge.400.and.INT(x).le.535)then !LOCATION = SKIP
            region = 23
         else if(INT(x).ge.800.and.INT(x).le.935)then !SECONDARY FEATURE = MULT PEAKS & TROUGHS
            region = 27
         else if(INT(x).ge.1150.and.INT(x).le.1295)then !LOCATION = SKIP
            region = 31
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
      end module mod_locate_region1
