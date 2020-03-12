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

      module mod_read_info_file
      contains
*     filename: read_info_file.f        lss

*     date of last revision:    10-Oct-2005     lss

*     subroutine of xanalysis
*     This code reads the *.info data files for importation into x2002 and thence to the database text files.

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



        subroutine read_info_file(INFO_FILENAME,
     +                            IDs,
     +                            included,
     +                            AA_import,STA_import,name_import,
     +                            comm_import,date_import,rec_import,
     +                            AP_import,RL_import,depth_import,
     +                            expname_import,
     +                            dchan,ref_chan,
     +                            ref_chan_refs,NAME_refs,AP_refs,
     +                            RL_refs,depth_refs,BDT_FILE,
     +                            recording,DB_FILES,IDs_extra)
      use mod_miscellaneous_subroutines
      use mod_read_and_write_pre_files


        INCLUDE 'x2000parameter.defs'

        integer included(MAX_NUM_CODES),IDs(MAX_NUM_CODES),
     +          commas(300),quotes(300),total_commas,total_quotes,
     +       IDs_extra(MAX_NUM_CODES),ihistnum

        character*300 str
        character*(*) INFO_FILENAME,expname_import,
     +       STA_import(MAX_NUM_CHAN,MAX_STA),date_import,
     +       AP_import(MAX_NUM_CHAN),RL_import(MAX_NUM_CHAN),
     +       depth_import(MAX_NUM_CHAN),AA_import(MAX_NUM_CHAN,MAX_AA),
     +       name_import(MAX_NUM_CHAN),comm_import(MAX_NUM_CHAN),
     +       rec_import,AP_refs(MAX_NUM_ARRAYS),RL_refs(MAX_NUM_ARRAYS),
     +       depth_refs(MAX_NUM_ARRAYS),NAME_refs(MAX_NUM_ARRAYS),
     +       ref_chan_refs(MAX_NUM_ARRAYS),
     +       dchan(MAX_NUM_CHAN),ref_chan(MAX_NUM_CHAN),DB_FILES,
     +       BDT_FILE

        character *220 pre_DBG_EXTRA
c        character*30 BDT_FILE_30
        character*40 blank40
        character*15 blank15
        character*9 blank9
c       character*4 c_cell
        character*4 dchan_refs                                     !dchan_refs is a discardable value
        character*2 blank2(5),recording
        character*1 ::  OK,extra = ' ',blank1
        
        logical*4 exist_2

        type extras_t
        integer ID
        character*4 cellname
        character*5 AP
        character*5 RL
        character*5 depth
        character*4 dchan
        character*4 ref_chan
        character*3 AA(MAX_AA)
        character*12 STA(MAX_STA)
        character*40 comment
        end type extras_t

        type(extras_t) extras(MAX_NUM_CHAN)
        integer :: num_extra = 0

        integer*4 :: mail_sent = 0
        save mail_sent

c        print '(''call read_info_file'')'

*     ***** initialize the arrays: *****

 1      do i = 1, MAX_NUM_CHAN
           dchan(i) = ' '
           ref_chan(i) = ' '
           AP_import(i)=' '
           RL_import(i)=' '
           depth_import=' '
           comm_import(i) = ' '
           name_import(i) = ' '
           EXTRAS(i)%ID = 0
           EXTRAS(i)%cellname = ' '
           EXTRAS(i)%AP = ' '
           EXTRAS(i)%RL = ' '
           EXTRAS(i)%depth = ' '
           EXTRAS(i)%dchan = ' '
           EXTRAS(i)%ref_chan = ' '
           do j = 1, MAX_AA
              AA_import(i,j) = ' '
              EXTRAS(i)%AA(j) = ' '
           end do
           do j = 1, MAX_STA
              STA_import(i,j) = ' '
              EXTRAS(i)%STA(j) = ' '
           end do
           EXTRAS(i)%comment = ' '
        end do
           blank40 = ' '
           blank15 = ' '
           blank9 = ' '
           blank1 = ' '
           do i = 1, 5
              blank2(i) = ' '
           end do
           do i = 1, MAX_NUM_ARRAYS
              AP_refs(i) = ' '
              RL_refs(i) = ' '
              depth_refs(i) = ' '
              NAME_refs(i) = ' '
              ref_chan_refs(i) = ' '
           end do
           rec_import = ' '
           expname_import = ' '
           dchan_refs = ' '

           if(INFO_FILENAME.eq.' ')return

           goto 10
 9         print '(//,T5,''Error opening '',A)',INFO_FILENAME
           goto 900

 10        OPEN (3,FILE=INFO_FILENAME,FORM='FORMATTED',STATUS='OLD',
     +          ERR=9)

           read (3,'(A)') str   !read in the 1st record
           do i = INDEX(str,',')+1, LEN(str) !look at the 1st record - it contains the date, exp name, & recording #
              if(str(i:i).ne.',')then !keep looking until find the date
                 date_import = str(i:i+10)
                 exit           
              end if
           end do
c           date_import = str(i:i+10)

           do i = 1,LEN(str)-11
              if(str(i:i+10).eq.'EXPERIMENT:')then
                 do j = i+11,LEN(str) !find the beginning letter of the experiment name
                    if((str(j:j).ne.',').and.(str(j:j).ne.' '))then
                       do k = j,LEN(str) !find the next comma
                          if(str(k:k).eq.',')then
                             expname_import = str(j:k-1)
                             name_length = k-j
                             expname_import(name_length+1:)=' '
                             goto 20
                          end if
                       end do
                    end if
                 end do
              end if
           end do

 20        do i = LEN(str),1,-1 !look backwards thru the 1st record until find the recording number
              if((str(i:i).ne.',').and.(str(i:i).ne.' '))exit
           end do
           if(str(i-1:i-1).ne.',')then !make sure that do not include commas in record number
              rec_import = str(i-1:i)
           else
              rec_import = str(i:i)
           end if

*     ***** FIND THE REFERENCE ELECTRODE INFORMATION FIRST: *****

           m = 0
           rewind 3             !start at the beginning of the file when searching for cell data
           read (3,'(A)') str
           read (3,'(A)') str   !skip records 1, 2, 3, and 4; cell data begins on row 5
           read (3,'(A)') str
           read (3,'(A)') str
 25        str = ' '
           read (3,'(A)',end=29) str
c           if((str(1:1).ne.',').and.(str(2:2).eq.'0'))then      !found ref electrode information
           if((str(1:1).ne.',').and.((str(2:2).eq.'0').or. !found ref electrode information
     +          ((str(2:2).eq.'A').and.(str(3:3).eq.'0')).or. !check for all possible ref electrodes (pA0 included)
     +          ((str(2:2).eq.'a').and.(str(3:3).eq.'0'))))then 
              m = m+1
              if(m.gt.MAX_NUM_ARRAYS)then
                 print '(''error:  too many reference electrodes!'')'
                 goto 29
              end if
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

*                         ***** load the arrays: ****

              do i = 1, total_commas
c                 if((i.eq.1).and.(commas(i).ne.1))then                !check for blank fields (denoted by ,,)
                 if(i.eq.1)then
                    if(str(1:commas(i)-1).ne.',')then
                       NAME_refs(m)=str(1:commas(i)-1)
                    else
                       NAME_refs(m)=' '
                    end if
                 end if
                 if(i.eq.3)then
                    if(str(commas(i)-1:commas(i)-1)
     +                   .ne.',')then !check for blank fields (denoted by ,,)
                       AP_refs(m) = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    else
                       AP_refs(m) = ' '
                    end if
                 end if
                 if(i.eq.4)then
                    if(str(commas(i)-1:commas(i)-1)
     +                   .ne.',')then !check for blank fields (denoted by ,,)
                       RL_refs(m) = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    else
                       RL_refs(m) = ' '
                    end if
                 end if
                 if(i.eq.5)then
                    if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                       depth_refs(m) =
     +                      str(commas(i-1)+1:commas(i)-1)
                    else
                       depth_refs(m) = ' '
                    end if
                 end if
                 if(i.eq.6)then
                    if(str(commas(i)-1:commas(i)-1)
     +                   .ne.',')then !check for blank fields (denoted by ,,)
                       dchan_refs =
     +                      str(commas(i-1)+1:commas(i)-1)
                    else
                       dchan_refs = ' '
                    end if
                 end if
                 if(i.eq.7)then
                    if(str(commas(i)-1:commas(i)-1)
     +                   .ne.',')then !check for blank fields (denoted by ,,)
                       ref_chan_refs(m) =
     +                      str(commas(i-1)+1:commas(i)-1)
                    else
                       ref_chan_refs(m) = ' '
                    end if
                 end if
              end do
           end if
           goto 25              !get the next record

 29        if(m.eq.0)then
              print '(//,T5,''***** WARNING *****   ''
     +             ''THE INFO FILE IS INCOMPLETE  *****'',//,T5,
     +             ''This info file does not contain reference ''
     +             ''electrode data.'',/,T5,
     +             ''Database-ready tables will not be written ''
     +             ''until a completed info file has been imported.''
     +             //,T5,''Please make a choice:'',//,T10,
     +             ''C -- Continue importing this info file'',/,T10,
     +             ''I -- Import a different info file'',//,T12,
     +             ''>> '',$)'
              read (*,'(A)') OK
              call upper_case(OK,LEN(OK))
              if(OK.eq.'C')then
                 print '(//,T5,''OK -- don''''t forget to import a ''
     +                ''completed info file at a later time!'',//,T8,
     +                ''<cr> to continue  >> '',$)'
                 read (*,'(A)')
              elseif (OK.eq.'I')then
                 close (unit=3)
                 isys=SYSTEM('retrieve_info.tcl info.tmp'//char(0))
                 inquire(FILE='info.tmp',EXIST=exist_2)
                 if(exist_2.eqv..TRUE.)then
                    INFO_FILENAME = 'info.tmp'
                    goto 1      !read the new info file
                 else
                    INFO_FILENAME = ' ' !don't want to import after all
                    goto 1
                 end if
              else
                 goto 29
              end if
           end if

*     ***** NOW FIND THE INFO FOR THE CELLS: *****

c           do icell = 1, MAX_NUM_CODES
c              if(included(icell).ne.0)then                      !for each cell in the QDT file:
c                 write(c_cell,'(I4)')icell                     !format c_cell same as IDimport
c                 call remove_all_blanks(c_cell,LEN(c_cell))

           num_extra = 0
           rewind 3             !start at the beginning of the file when searching for cell data
           read (3,'(A)') str
           read (3,'(A)') str   !skip records 1, 2, 3, and 4; cell data begins on row 5
           read (3,'(A)') str
           read (3,'(A)') str
 30        str = ' '
           read (3,'(A)',end=900) str



           if(((str(1:1).eq.',').and.(str(2:2).eq.',')).or.
     +          ((str(1:1).ne.',').and.((str(2:2).eq.'0').or. !found ref electrode information
     +          ((str(2:2).eq.'A').and.(str(3:3).eq.'0')).or. !check for all possible ref electrodes (pA0 included)
     +          ((str(2:2).eq.'a').and.(str(3:3).eq.'0')))))goto 30 !record is empty or a reference value - get next one



c           if(((str(1:1).eq.',').and.(str(2:2).eq.',')).or.
c     +          ((str(1:1).ne.',').and.(str(2:2).eq.'0')))
c     +          goto 30         !record is empty or a reference value - get next one
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
           
*     ***** load the arrays: ****

           read (str(commas(1)+1:commas(2)-1),'(I4)') icell
           if(included(icell).ne.0)then !this cell included in QDT file
              if(IDs_extra(icell).eq.0)extra = 'n'
           else
              if(IDs_extra(icell).ne.0)then
                 extra = 'y'    !this cell is an "extra" cell
                 num_extra = num_extra+1
                 EXTRAS(num_extra)%ID = icell
c                 print '(''extra = '',I4)',icell
c                 read (*,'(A)')
              else
                 print "('read_info: skipping icell = ',I4,
     +                ' because it is not included and not extra')",
     +                icell
                 goto 30
              end if
           end if

           if (IDs(icell).eq.0.and.extra.eq.'n') then
              print "('read_info: skipping icell = ',I4,
     +             ' because it is not in IDs and not extra')",icell
              if(mail_sent.eq.0) then
                 isys=system('echo not IDs or extra '
     +                //'| mailx -s xanalysis_bug '
     +                //'roconnor@cisc1.hscnet.hsc.usf.edu')
                 mail_sent = 1
              end if
              goto 30
           end if
           
c     if(str(commas(1)+1:commas(2)-1).eq.c_cell) then                !OK - have the right data record
           do i = 1, total_commas
              if((i.eq.1).and.(commas(i).ne.1))then !check for blank fields (denoted by ,,)
                 if(extra.eq.'n')then
                    name_import(IDs(icell))=str(1:commas(i)-1)
                 else if(extra.eq.'y')then
                    EXTRAS(num_extra)%cellname=str(1:commas(i)-1)
                 end if
              end if
              if(i.eq.3)then
                 if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                    if(extra.eq.'n')then
                       AP_import(IDs(icell)) = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%AP = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    end if
                 else
                    if(extra.eq.'n')then
                       AP_import(IDs(icell)) = ' '
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%AP = ' '
                    end if
                 end if
              end if
              if(i.eq.4)then
                 if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                    if(extra.eq.'n')then
                       RL_import(IDs(icell)) = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%RL = 
     +                      str(commas(i-1)+1:commas(i)-1)
                    end if
                 else
                    if(extra.eq.'n')then
                       RL_import(IDs(icell)) = ' '
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%RL = ' '
                    end if
                 end if
              end if
              if(i.eq.5)then
                 if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                    if(extra.eq.'n')then
                       depth_import(IDs(icell)) =
     +                      str(commas(i-1)+1:commas(i)-1)
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%depth =
     +                      str(commas(i-1)+1:commas(i)-1)
                    end if
                 else
                    if(extra.eq.'n')then
                       depth_import(IDs(icell)) = ' '
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%depth = ' '
                    end if
                 end if
              end if
              if(i.eq.6)then
                 if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                    if(extra.eq.'n')then
                       dchan(IDs(icell)) =
     +                      str(commas(i-1)+1:commas(i)-1)
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%dchan =
     +                      str(commas(i-1)+1:commas(i)-1)
                    end if
                 else
                    if(extra.eq.'n')then
                       dchan(IDs(icell)) = ' '
                    elseif(extra.eq.'n')then
                       EXTRAS(num_extra)%dchan = ' '
                    end if
                 end if
              end if
              if(i.eq.7)then
                 if(str(commas(i)-1:commas(i)-1).ne.',')then !check for blank fields (denoted by ,,)
                    if(extra.eq.'n')then
                       ref_chan(IDs(icell)) =
     +                      str(commas(i-1)+1:commas(i)-1)
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%ref_chan =
     +                      str(commas(i-1)+1:commas(i)-1)
                    end if
                 else
                    if(extra.eq.'n')then
                       ref_chan(IDs(icell)) = ' '
                    elseif(extra.eq.'y')then
                       EXTRAS(num_extra)%ref_chan = ' '
                    end if
                 end if
              end if
              if(i.eq.8)then
                 do j = 1,MAX_AA
                    if(str((commas(i+j-1)-1):(commas(i+j-1)-1))
     +                   .ne.',')then
                       if(extra.eq.'n')then
                          AA_import(IDs(icell),j) = 
     +                         str(commas(i+j-2)+1:commas(i+j-1)-1)
                       elseif(extra.eq.'y')then
                          EXTRAS(num_extra)%AA(j) = 
     +                         str(commas(i+j-2)+1:commas(i+j-1)-1)
                       end if
                    else
                       if(extra.eq.'n')then
                          AA_import(IDs(icell),j)=' '
                       elseif(extra.eq.'y')then
                          EXTRAS(num_extra)%AA(j)=' '
                       end if
                    end if
                 end do
              end if
              if(i.eq.(8+MAX_AA))then
                 do j = 1,MAX_STA
                    if(str((commas(i+j-1)-1):(commas(i+j-1)-1))
     +                   .ne.',')then
                       if(extra.eq.'n')then
                          STA_import(IDs(icell),j) = 
     +                         str(commas(i+j-2)+1:commas(i+j-1)-1)
                       elseif(extra.eq.'y')then
                          EXTRAS(num_extra)%STA(j) = 
     +                         str(commas(i+j-2)+1:commas(i+j-1)-1)
                       end if
                    else
                       if(extra.eq.'n')then
                          STA_import(IDs(icell),j)=' '     
                       elseif(extra.eq.'y')then
                          EXTRAS(num_extra)%STA(j)=' '     
                       end if
                    end if
                 end do
              end if
           end do

           if(extra.eq.'n')then
              comm_import(IDs(icell))=
     +             str(commas(8+MAX_AA+MAX_STA-1)+1:
     +             commas(8+MAX_AA+MAX_STA-1)+40)
           elseif(extra.eq.'y')then
              EXTRAS(num_extra)%comment=
     +             str(commas(8+MAX_AA+MAX_STA-1)+1:
     +             commas(8+MAX_AA+MAX_STA-1)+40)
           end if

           goto 30

c     else                    !cell name is not a match - read next record
c     goto 30
c     end if
c     end if
           
c900  end do
c999     close(3)

 900       close(3)

*     ***** write the "extra cells predbg" file now: *****

           if(num_extra.gt.0)then
              call concat(DB_FILES,'_extra_cell_preDBG', !precursor file for indiv EXTRA cells
     +             pre_DBG_EXTRA,l)
              open (29,FILE=pre_DBG_EXTRA(1:l),FORM='FORMATTED',
     +             ACCESS='DIRECT',
     +             RECL=DBG_RECL)
c              print '(''read_info: num_extra = '',I)',num_extra
c              read (*,'(A)')
c              ihistnum = 0
              ihistnum = 1
              do i = 1, num_extra
c                 call write_preDBG(29,i,ios,date_import,BDT_FILE_30,
                 call write_preDBG(29,i,ios,date_import,BDT_FILE,
     +                recording,ihistnum,
     +                EXTRAS(i)%cellname,EXTRAS(i)%ID,blank15,       
     +                EXTRAS(i)%AP,EXTRAS(i)%RL,
c     +             EXTRAS(i)%depth,blank1,EXTRAS(i)%comment,blank1,
     +                EXTRAS(i)%depth,blank1,blank40,blank1,
     +                EXTRAS(i)%AA(1),EXTRAS(i)%AA(2),EXTRAS(i)%AA(3),
     +                EXTRAS(i)%AA(4),
     +                EXTRAS(i)%AA(5),EXTRAS(i)%AA(6),EXTRAS(i)%AA(7),
     +                EXTRAS(i)%AA(8),
     +                EXTRAS(i)%AA(9),EXTRAS(i)%AA(10),
     +                EXTRAS(i)%STA(1),EXTRAS(i)%STA(2),
     +                EXTRAS(i)%STA(3),EXTRAS(i)%STA(4),
     +                EXTRAS(i)%STA(5),EXTRAS(i)%STA(6),
     +                EXTRAS(i)%STA(7),EXTRAS(i)%STA(8),
     +                EXTRAS(i)%STA(9),EXTRAS(i)%STA(10),
     +                blank9,blank9,blank9,blank9,blank9,blank9,blank9,
     +                blank9,blank9,blank9,
     +                blank9,blank9,blank9,blank9,blank9,blank9,blank9,
     +                blank9,blank9,blank9,
     +                blank9,blank9,blank9,blank9,blank9,blank9,blank9,
     +                blank9,blank9,blank9,
     +                blank9,blank9,blank9,blank9,blank9,
     +                blank2,blank2,blank2,blank2,blank2,blank2,blank2,
     +                blank2,blank2,blank2,
     +                blank2,blank2,blank2,blank2,blank2,blank2,blank2,
     +                blank2,blank2,blank2,
     +                blank2,blank2,blank2,blank2,blank2,blank2,blank2,
     +                blank2,blank2,blank2,
     +                blank2,blank2,blank2,blank2,blank2)
                 print '(''ID of extra = '',I12)',EXTRAS(i)%ID
              end do
              close (29)
           end if

           return
           end
      end module mod_read_info_file
