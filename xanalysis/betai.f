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

      module mod_betai
      contains
      FUNCTION BETAI(A,B,X)
      use mod_betacfkm
      use mod_gammln
      IF(X.LT.0..OR.X.GT.1.)then
         PRINT '(''bad argument X in BETAI.  Hit ENTER to continue.'')'
         READ (*,'(A)')
      END IF
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *      +A*ALOG(X)+B*ALOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
      END
      end module mod_betai
