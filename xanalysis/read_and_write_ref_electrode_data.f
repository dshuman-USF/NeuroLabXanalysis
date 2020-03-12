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

      module mod_read_and_write_ref_electrode_data
      contains
*     filename = read_and_write_ref_electrode_data.f

*     date of last revision = 12-Jan-2005       lss

        subroutine read_ref_electrode_data(BDT_FILE,DB_FILES,
     +                     date_exp,recording,exp_name,
     +                     NAME_refs,AP_refs,RL_refs,depth_refs,
     +                     dchan,ref_chan,ref_chan_refs,
     +                     total_num_electrode_arrays)
      use mod_miscellaneous_subroutines


      include 'x2000parameter.defs'

      integer total_num_electrode_arrays

      character*(*) BDT_FILE,DB_FILES,date_exp,recording,exp_name,
     +              NAME_refs(MAX_NUM_ARRAYS),
     +              AP_refs(MAX_NUM_ARRAYS),
     +              RL_refs(MAX_NUM_ARRAYS),
     +              depth_refs(MAX_NUM_ARRAYS),
     +              ref_chan_refs(MAX_NUM_ARRAYS),
     +              dchan(MAX_NUM_CHAN),
     +              ref_chan(MAX_NUM_CHAN)

        character*250 IN_FILE
        character*40 fmt

        if(.false.)print *,BDT_FILE !suppress unused variable warning

        do i = 1, MAX_NUM_ARRAYS
           AP_refs(i) = ' '
           RL_refs(i) = ' '
           depth_refs(i) = ' '
           NAME_refs(i) = ' '
           ref_chan_refs(i) = ' '
        end do

        goto 10
 1      print '(I12,'': error opening '',A)',ios1,IN_FILE
        goto 100
 2      print '(I12,'': error reading '',A)',ios2,IN_FILE
        stop
 3      print '(I12,'': error reading '',A)',ios3,IN_FILE
        goto 100

 10     call concat(DB_FILES,'_ARRAY_REFS.sav',IN_FILE,l)
        open (UNIT=3,FILE=IN_FILE(1:l),FORM='FORMATTED',
     +       ACCESS='SEQUENTIAL',ERR=1,IOSTAT=ios1)

        call filesize (IN_FILE(1:l)//char(0), ifilesize)
        total_num_electrode_arrays = ifilesize / 212
        do i = 1, total_num_electrode_arrays
           read(3,'(8A)',err=2,IOSTAT=ios2) 
     +          date_exp,recording,exp_name,
     +          NAME_refs(i),AP_refs(i),RL_refs(i),depth_refs(i),
     +          ref_chan_refs(i)
        end do

        close (unit=3)
        

        call concat(DB_FILES,'_UNIT_REFS.sav',IN_FILE,l)
        open (3,FILE=IN_FILE(1:l),FORM='FORMATTED',
     +       ACCESS='SEQUENTIAL',ERR=1,IOSTAT=ios1)

        write(fmt,'(A1,I12,A5,I12,A5)')
     +       '(',
     +       MAX_NUM_CHAN,
     +       '(A4),',
     +       MAX_NUM_CHAN,
     +       '(A4))'
c     format string length = 35

        read(3,fmt,err=3,IOSTAT=ios3),
     +       dchan,ref_chan

 100    close (unit=3)

        return
        end


*     <<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        subroutine write_ref_electrode_data(BDT_FILE,DB_FILES,
     +                     date_exp,recording,exp_name,
     +                     NAME_refs,AP_refs,RL_refs,depth_refs,
     +                     dchan,ref_chan,ref_chan_refs)
      use mod_miscellaneous_subroutines


      include 'x2000parameter.defs'

      character*(*) BDT_FILE,DB_FILES,date_exp,recording,exp_name,
     +              NAME_refs(MAX_NUM_ARRAYS),
     +              AP_refs(MAX_NUM_ARRAYS),
     +              RL_refs(MAX_NUM_ARRAYS),
     +              depth_refs(MAX_NUM_ARRAYS),
     +              ref_chan_refs(MAX_NUM_ARRAYS),
     +              dchan(MAX_NUM_CHAN),
     +              ref_chan(MAX_NUM_CHAN)

        character*250 OUT_FILE
        character*40 fmt

        if(.false.)print *,BDT_FILE !suppress unused variable warning

        goto 10
 1      print '(I12,''; error opening '',A)',ios1,OUT_FILE
        return
 2      print '(I12,''; i='',I12,''; error writing '',A)',
     +       i,ios2,OUT_FILE
        goto 50
 3      print '(I12,''error writing '',A)',ios3,OUT_FILE
        goto 100

c 10     OUT_FILE=DB_FILES//'_ARRAY_REFS.sav'
c        call remove_all_blanks(OUT_FILE,LEN(OUT_FILE))
c        call strlength(OUT_FILE,LEN(OUT_FILE),l)

 10     call concat(DB_FILES,'_ARRAY_REFS.sav',OUT_FILE,l)
        open (UNIT=3,FILE=OUT_FILE(1:l),FORM='FORMATTED',
     +       ACCESS='SEQUENTIAL',ERR=1,IOSTAT=ios1)
        rewind 3                                               !overwrite anything already in the file

        j = 0
        do i = 1, MAX_NUM_ARRAYS
           if(NAME_refs(i).ne.' ')then
              write(3,'(8A)',err=2,IOSTAT=ios2) 
     +             date_exp,recording,exp_name,
     +             NAME_refs(i),AP_refs(i),RL_refs(i),depth_refs(i),
     +             ref_chan_refs(i)
              j = j+1
           end if
        end do
        total_num_electrode_arrays = j
 50     close (unit=3)

        call concat(DB_FILES,'_UNIT_REFS.sav',OUT_FILE,l)
        open (3,FILE=OUT_FILE(1:l),FORM='FORMATTED',
     +       ACCESS='SEQUENTIAL',ERR=1,IOSTAT=ios1)
        rewind 3                                               !overwrite anything already in the file

        write(fmt,'(A1,I12,A5,I12,A5)')
     +       '(',
     +       MAX_NUM_CHAN,
     +       '(A4),',
     +       MAX_NUM_CHAN,
     +       '(A4))'
c     format string length = 35

        write(3,fmt,err=3,IOSTAT=ios3)
     +       dchan,ref_chan

 100    close (unit=3)
        
        return
        end

      
      end module mod_read_and_write_ref_electrode_data
