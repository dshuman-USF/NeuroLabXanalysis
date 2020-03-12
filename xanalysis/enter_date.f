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

      module mod_enter_date
      contains
*     subroutine enter_date

*     filename: enter_date.f            lss     24-Sep-2003
*     date of last revision: 10-Nov-2003

*     This code solicits a date (dd-MMM-yyyy) if none is passed, checks it for formatting,
*       requires a year between 1950 and 2999, requires an appropriate day,
*       forces the 3-character month to all upppercase and requires an
*       appropriate month.  Beyond that, you're on your own for accuracy.


      subroutine enter_date(new_date)

      character*11 new_date


 100  if(new_date.eq.' ')then
         print '(/,''ENTER experiment date ''
     +''(dd-MMM-yyyy; eg. 01-JAN-2000) >> '',$)'
         read (*,'(A11)',err=100) new_date
      end if

      if(new_date.eq.'0')then   !Easter egg default
         new_date = '01-APR-2003'
         return
      end if

      if(new_date(2:2).eq.'-')then !be sure that have 2 digits for day
         new_date(2:11) = new_date(1:10)
         new_date(1:1) = '0'
      end if

      if((INDEX(new_date,' ').ne.0).or. !check for blanks in the date entered
     +     (new_date(3:3).ne.'-').or. !dashes must be in the right places
     +     (new_date(7:7).ne.'-'))then 
         print '(/,T5,''This date has been entered using an '',
     +''incorrect format.'',/,T5
     +''You must use 2 digits for the day, '',
     +''3 characters for the month, and '',
     +''4 digits to designate the year.'',/,T5,
     +''Separate day, month, and year with ''
     +''a dash.'',//)'
         new_date = ' '
         goto 100
      end if

*     ***** check for format and validity of year: *****

      if((new_date(8:11).gt.'2999').or.                   
     +     (new_date(8:11).lt.'1950'))then
         print '(/,T5,''Invalid year ('',A4,'').  ''
     +''Please re-enter date.'')',new_date(8:11)
         new_date = ' '
         goto 100
      end if

      do i = 8,11
         if((ICHAR(new_date(i:i)).lt.48).and.
     +        (ICHAR(new_date(i:i)).gt.57))then
            print '(/,T5,''Invalid format for year.  ''
     +''Please re-enter date.'')'
            new_date = ' '
            goto 100
         end if
      end do

*     ***** check for format and validity of day: *****

      if((new_date(1:2).gt.'31').or.                   
     +     (new_date(1:2).lt.'01'))then
         print '(/,T5,''Invalid day ('',A2,'').  ''
     +''Please re-enter date.'')',new_date(1:2)
         new_date = ' '
         goto 100
      end if

      do i = 1,2
         if((ICHAR(new_date(i:i)).lt.48).and.
     +        (ICHAR(new_date(i:i)).gt.57))then
            print '(/,T5,''Invalid format for day ('',A2,'').  ''
     +''Please re-enter date.'')',new_date(1:2)
            new_date = ' '
            goto 100
         end if
      end do

*     ***** check for format and validity of month: *****

      do i = 4,6
         if((ICHAR(new_date(i:i)).lt.65).or.
     +        (ICHAR(new_date(i:i)).gt.122))then
            print '(/,T5,''Invalid format for month ('',A3,'').  ''
     +''Please re-enter date.'')',new_date(4:6)
            new_date = ' '
            goto 100
         end if
      end do

      do i = 4,6                !force month to be all caps
         if((ICHAR(new_date(i:i)).ge.97).and.
     +        (ICHAR(new_date(i:i)).le.122))then
            new_date(i:i)=CHAR(ICHAR(new_date(i:i))-32)
         end if
      end do

      if((new_date(4:6).ne.'JAN').and.
     +     (new_date(4:6).ne.'FEB').and.
     +     (new_date(4:6).ne.'MAR').and.
     +     (new_date(4:6).ne.'APR').and.
     +     (new_date(4:6).ne.'MAY').and.
     +     (new_date(4:6).ne.'JUN').and.
     +     (new_date(4:6).ne.'JUL').and.
     +     (new_date(4:6).ne.'AUG').and.
     +     (new_date(4:6).ne.'SEP').and.
     +     (new_date(4:6).ne.'OCT').and.
     +     (new_date(4:6).ne.'NOV').and.
     +     (new_date(4:6).ne.'DEC'))then
         print '(/,T5,''Invalid month ('',A3,'').  ''
     +''Please re-enter date.'')',new_date(4:6)
         new_date = ' '
         goto 100
      end if


      return
      end

      
      end module mod_enter_date
