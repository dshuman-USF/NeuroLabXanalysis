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

      module mod_get_info_IDcodes
      contains
*     filename: get_info_IDcodes.f        lss

*     date of last revision:    23-Sep-2005     lss

*     subroutine of xanalysis to retrieve the IDcodes of the cells in an info file.

*     10-Jan-2005       lss
*       rows corresponding to the coordinates of the reference electrode for each array have added to the info file;
*               These cells are given char*2 names.  The first character indicates the brain region; the second
*               character is always 0 (zero).
*       columns have also been added; they contain digital channel information:
*               dchan = the digital channel on which the cell's activity was recorded
*               ref = the digital channel of the reference electrode for the array used to record the cell's activity
*         every cell will have data in these 2 new columns

*     THIS SUBROUTINE READS IN AN INFO FILE, WHICH IS ALWAYS A CSV FILE.  
*     THE FOLLOWING DATA IS GLEANED FROM INFO FILES: 
*               date (dd-mmm-yyyy; char*11)
*               name of experiemnt (char*20)
*               recording (char*2)
*               cell name (char*4)
*               AP (char*5)
*               RL (char*5)
*               depth (char*5)
*               dchan (char*4)         [the digital channel on which the cell's activity was recorded]
*               ref_chan (char*4)      [the digital channel of the reference electrode for the array used to record
*                                               the cell's activity]
*               AA (char*5, x10)
*               STA (char*12, x10)
*               comments (char*40)

*          an ASCII file in comma-separated values (*.csv) format generated from an Excel worksheet
*                       - filename = *.info 
*                       - each line contains all the info in that row of the spreadsheet; fields are separated by commas
*                       - all info for one cell contained in one row
*                       - note that a field containing a comma is further set off by quotes
*                       - column headings are in row 3; cell data on rows 4 and up
*              1  EXP. DATE :   ,,,,,EXPERIMENT:,,,,,RECORDING:,,,,,,,,,,,,,,,,
*              2  ,,,,,,,,,,,,,,,,,,,,,,,,,,
*              3  ,,,,,AA,,,,,,,,,,STA,,,,,,,,,,,,
*              4  UNIT NAME ,MERGED CH. #,A/P,R/L,DEPTH,dchan,ref,cord,RLN,vagus,cVRG,rVRG,rPRG,lPRG,AA,AA,AA,phrenic,RLN,c.vagus,lumbar,cerv.sym.,exp.laryn,splanchnic,STA,STA,STA,COMMENTS
*              5  R65,65,6.5,6.5,6.5,8,3,aa,aa,aa,aa,aa,aa,aa,,,,,2 8,2 8,2 8,2 8,2 8,2 8,2 8,,,,"This, is, a, comment."
*              6  ,,,,,,,,,,,,,,,,,,,,,,,,



        subroutine get_info_IDcodes(INFO_FILENAME,IDs_info)


        INCLUDE 'x2000parameter.defs'

        integer IDs_info(MAX_NUM_CODES),commas(300),
     +       quotes(300),total_commas,total_quotes

        character*300 str
        character*(*) INFO_FILENAME


*     ***** initialize the arrays: *****

        do i = 1, MAX_NUM_CODES
           IDs_info(i) = 0
        end do

        if(INFO_FILENAME.eq.' ')return

        goto 10
 9      print '(//,T5,''Error opening '',A)',INFO_FILENAME
        goto 900

 10     OPEN (3,FILE=INFO_FILENAME,FORM='FORMATTED',STATUS='OLD',ERR=9)


*     ***** NOW FIND THE INFO FOR THE CELLS: *****

           num_extra = 0
           rewind 3             !start at the beginning of the file when searching for cell data
           read (3,'(A)') str
           read (3,'(A)') str   !skip records 1, 2, 3, and 4; cell data begins on row 5
           read (3,'(A)') str
           read (3,'(A)') str
 30        str = ' '
           read (3,'(A)',end=900) str
           if(((str(1:1).eq.',').and.(str(2:2).eq.',')).or.
     +          ((str(1:1).ne.',').and.(str(2:2).eq.'0')))
     +          goto 30         !record is empty or a reference value - get next one
           j = 0
           k = 0
           commas = 0
           quotes = 0
           do i = 1, LEN(str)
              if(str(i:i).eq.',')then
                 j=j+1
                 commas(j) = i
              end if
              if(str(i:i).eq.'"')then
                 k=k+1
                 quotes(k) = i
              end if
           end do
           total_commas = j
           total_quotes = k
           
*     ***** find the IDcodes: ****

           read (str(commas(1)+1:commas(2)-1),'(I4)') icell
           IDs_info(icell) = 1
c           print '(''IDs_info('',I,'') = '',I)',icell,IDs_info(icell)
           
           goto 30                                      !read another line of the info.tmp file

 900       close(3)

           return
           end
      end module mod_get_info_IDcodes
