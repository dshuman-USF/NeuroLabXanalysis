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

      module mod_locate_region2
      contains
      subroutine locate_region2(x,y,region)
*
*
*       date of last revision = 08-aug2003      lss
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
*       *       65              800<x<850,  90<y<140      USER-DEFINED DISPLAY TIME     *
*       *       7               900<x<950,      "         SLOW, OK                      *
*       *       8               1000<x<1050,    "         MEDIUM                        *
*       *       9               1100<x<1150,    "         FAST, SELECT                  *
*       *       90              1300<x<1350,    "         START OVER                    *
*       *       92              200<x<250,      "         VIEW ALL                      *
*       *       93              325<x<375,      "         VIEW SIGNIFICANT CCHs ONLY    *
*       *       931             450<x<500,      "         VIEW FLAT CCHs ONLY           *
*       *       94              1400<x<1450, 20<y<70      CANCEL                        *
*       *       91              1400<x<1450, 110<y<160    GO, DONE                      *
*       *                                                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       link with x200*
*
*
      real x,y
      integer region
      region = 0                !initialize

*       look at y-coordinate first
      if(INT(y).ge.90.and.INT(y).le.140)then !90 <= y <= 140
         if(INT(x).ge.800.and.INT(x).le.850)then !USER-DEFINED display time
            region = 65
         else if(INT(x).ge.900.and.INT(x).le.950)then !SLOW, OK button
            region = 7
         else if(INT(x).ge.1000.and.INT(x).le.1050)then !MEDIUM button
            region = 8
         else if(INT(x).ge.1100.and.INT(x).le.1150)then !FAST, SELECT button
            region = 9
         else if(INT(x).ge.1300.and.INT(x).le.1350)then !START OVER button
            region = 90
c          else if(INT(x).ge.1400.and.INT(x).le.1450)then               !GO
c             region = 91
         else if(INT(x).ge.200.and.INT(x).le.250)then !VIEW ALL button
            region = 92
         else if(INT(x).ge.325.and.INT(x).le.375)then !VIEW SIGNIFICANT CCHs ONLY button
            region = 93
         else if(INT(x).ge.450.and.INT(x).le.500)then !VIEW FLAT CCHs ONLY button
            region = 931
c          else if(INT(x).ge.50.and.INT(x).le.100)then                  !CANCEL
c             region = 94
         end if
      end if

      if(INT(x).ge.1400.and.INT(x).le.1450)then !1400 <= x <= 1450
         if(INT(y).ge.20.and.INT(y).le.70)then !CANCEL button
            region = 94
         elseif(INT(y).ge.110.and.INT(y).le.160)then !GO button
            region = 91
         end if
      end if


      return
      end
      end module mod_locate_region2
