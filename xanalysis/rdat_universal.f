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

      module mod_rdat_universal
      contains
*
*       ************************************************************************
*       THIS VERSION IS rdat_universal_int.f
*
*       date of last revision = 03-Nov-2004     lss
*
*       aug-2002        lss
*        rdat modified to accommodate edt files too (require I5,I10 read format)
*
*       lss 04-apr-2002 version 4
*        rdat modified to accommodate ddt files as well as adt and bdt files
*
*       lss  07-may-98
*         ID codes now loaded into IDs according to code number
*           i.e., for ID code = 23, IDs(23) = 23
*       lss  01-apr-98
*         code added to place ID codes into IDs in numerical order
*       lss  20-apr-98
*         increase in max # spike trains and # events
*       *************************************************************************
*       *                                                                       *
*       *       for adt and bdt files (0.5 msec resolution):                    *
*       *               CLOCK COUNTS (MUL BY 0.5 TO GET MSEC., F.P.)            *
*       *                       (real*8 factor = 0.5)                           *
*       *       for ddt and edt files (0.1 msec resolution):                    *
*       *               CLOCK COUNTS (MUL BY 0.1 TO GET MSEC., F.P.)            *
*       *                       (real*8 factor = 0.1)                           *
*       *                                                                       *
*       *       READS IN UP TO N SPIKE CODES                                    *
*       *       MAX_NUM_EVENTS=48000;MAX NUM EVENTS PER CODE                    *
*       *       *.adt INPUT DATA HAS THE FORMAT I2,I8: EVENT CODE,CLOCK COUNTSC *
*       *       *.bdt INPUT DATA HAS THE FORMAT I5,I8 ..FIRST 2 RECORDS         *
*       *           ARE: 11,1111111                                             *
*       *                11,LSTTIM WHERE LSTTIM IS LAST CLOCK COUNT IN FILE     *
*       *            ...THUS 1111111 IS USED TO DISTINGUISH .ADT FROM .BDT      *
*       *               FILES.                                                  *
*       *       *.ddt INPUT DATA HAS THE FORMAT I5,I8 ..FIRST 2 RECORDS         *
*       *           ARE: 22,2222222                                             *
*       *                22,LSTTIM WHERE LSTTIM IS LAST CLOCK COUNT IN FILE     *
*       *            ...THUS 2222222 IS USED TO DISTINGUISH .ddt FROM .bdt      *
*       *               AND .adt FILES.                                         *
*       *       *.edt INPUT DATA HAS THE FORMAT I5,I10 ..FIRST 2 RECORDS        *
*       *           ARE: 33,3333333                                             *
*       *                33,LSTTIM WHERE LSTTIM IS LAST CLOCK COUNT IN FILE     *
*       *            ...THUS 3333333 IS USED TO DISTINGUISH .edt FROM .ddt,.bdt *
*       *               AND .adt FILES.                                         *
*       *                                                                       *
*       *  ** BE SURE CLOCK, MAX_NUM_EVENTS ARE 32 BIT INTEGER VARIABLES        *
*       *                                                                       *
*       *************************************************************************
*       
      SUBROUTINE rdat6(SPIKETIMES,ITAL,IDs,IDT,
     +     BDT_FILE)
      use mod_miscellaneous_subroutines
*       
*
*
      include 'x2000parameter.defs'
*
      DOUBLE PRECISION, allocatable, intent(out) :: SPIKETIMES(:,:) !real*8
      DIMENSION ITOTAL_TALLY(MAX_NUM_CODES)
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),count
      double precision factor,fl_clock,endtime
      integer*4 CLOCK,CODE,CLOCK_FIRST
      integer skip, dots
      character*(*) BDT_FILE
      character*20 c_format
      character*3 file_ext,data_ext
      CLOCK_FIRST = 0
*
*
      call strlength(BDT_FILE,LEN(BDT_FILE),length)
      file_ext = BDT_FILE(length-2:length)


      goto 2

 1    print '(///,''*** DATA TYPE MISMATCH ***'',//,
     +T3,''The file header indicates that the data in file '',
     +A,'' is in '',A3,'' format,'',
     +/,T3,''but the file name indicates that the data is in '',
     +A3,'' format.'',//,
     +T3,''Please resolve the discepancy before running this ''
     +''program again.'',
     +//,T8,''<cr> to continue  >> '',$)',
     +     BDT_FILE(1:length),data_ext,file_ext
      read (*,'(A)')
      BDT_FILE = 'x'
      return

      
 2    if(file_ext.eq.'adt')then !different types of data files require different read formats
         c_format = '(I2,I8)'
         factor = 0.5D0         !precision = 0.5 msec.
      elseif(file_ext.eq.'bdt')then                  
         c_format = '(I5,I8)'
         factor = 0.5D0         !precision = 0.5 msec.
c        elseif(file_ext.eq.'ddt')then                  
c           c_format = '(I5,I8)'
c           factor = 0.1                !precision = 0.1 msec.
      elseif(file_ext.eq.'edt')then                  
         c_format = '(I5,I10)'
         factor = 0.1D0         !precision = 0.1 msec.
      else
         print '(//,T5,''invalid file type ... abort'')'
         file_ext='x'
         return
      endif


*
*       ***** check for compatibility between file name and data: *****
*
      skip=0                    !count the number of records that need to be skipped
 10   READ (1,FMT=c_format,END=2000)CODE,CLOCK !read 1st record 
      skip = skip + 1
      if (CODE.eq.0) goto 10    !this is a blank record -- read another
      if(CLOCK.eq.1111111)then
         data_ext = 'bdt'
         skip = skip + 1
c        elseif(CLOCK.eq.2222222)then
c           data_ext = 'ddt'
c           skip = skip + 1
      elseif(CLOCK.eq.3333333)then
         data_ext = 'edt'
         skip = skip + 1
      else
         data_ext = 'adt'
      end if

      if(file_ext.ne.data_ext)goto 1 !uh-oh -- data type mismatch
      print '(/,''Please wait ... I''''m initializing ''
     +''the data arrays and sacrificing the chicken.'')'
      
*       
*       ***** INITIALIZE ARRAYS BEFORE FILLING *****
*

      ITAL=0                  
      IDs=0
      IDT=0
      ITOTAL_TALLY = 0
      
*
*       *********************************************************
*       *                                                       *
*       *               READ IN THE DATA FILE                   *
*       *                                                       *
*       *********************************************************
*
*
      print '(/,''Please wait a little longer ... I''''m reading '',
     +A,'' ... This may take a few minutes.'')',
     +     BDT_FILE(1:length)
      CLOCK_FIRST = 0
      max_num_events = 0
      rewind 1                  !go back to the beginning of the input file
      do i = 1,skip             !read <skip> number of non-data records and discard them
         read (1,FMT=c_format,end=2000) CODE,CLOCK 
      end do                
*
*
      dots = 0
      count = 0
c110    READ (1,FMT=c_format,END=150) CODE,CLOCK        !read a data record
 110  READ (1,FMT=c_format,END=150) CODE !read a data record
      count=count+1
      if(MOD(count,100000).eq.0)then
         dots = dots + 1
         do i = 1, dots
            print '(''. '',$)'
         end do
         if(dots.eq.20)then
            dots = 0
            print '(//,''I''''m still working!'')'
         end if
         print '()'
      end if
*
*          ***** reject data if CODE indicates it is from an analog channel or is  *****
*          *****                a voltage data point:                           *****
*
      if(IDT.gt.MAX_NUM_CHAN) then !too many ID codes 
         print '(''MORE THAN '',I3,'' ID CODES!  FATAL ERROR'')',
     +        MAX_NUM_CHAN
         stop                   !halt the program
      end if
*
      if((CODE.eq.0).or.(CODE.gt.999)) goto 110 !reject data...read another record 
*
      if(IDs(CODE).eq.0)then    !this code has not already been "registered"
         IDs(CODE)=CODE         !load the array with this ID code
         IDT=IDT+1              !keep track of how many ID codes have been registered
      end if
      ITOTAL_TALLY(IDs(CODE))=ITOTAL_TALLY(IDs(CODE))+1
      if (ITOTAL_TALLY(IDs(CODE)).gt.max_num_events) max_num_events = ITOTAL_TALLY(IDs(CODE))
      goto 110                  !read another record
*
*
*       ***** Go through IDs(); number sequentially the locations       
*       *****  that correspond to the existence of a code. For example: 
*       *****
*       *****   IDs(3) = 3              IDs(36) = 36    IDs(401) = 401 
*       *****                   (in this case, IDT = 3)
*       *****
*       ***** Since the contents of IDs() will be used to point to
*       *****   the location of data in SPIKETIMES(), we have to change the
*       *****   contents of IDs to reflect a sequential order:
*       *****
*       *****   IDs(3) = 1              IDs(36) = 2             IDs(401) = 3
*       *****
*       *****   and data will be loaded for unit #3 in SPIKETIMES(1:ITAL(3),IDs(3))
*       
*       
 150  j=0
      do i = 1,MAX_NUM_CODES
         if(IDs(i).ne.0)then
            j=j+1
            IDs(i) = j
         end if
      end do                    !OK - contents of IDs now changed
      allocate (spiketimes (max_num_events, max_num_chan))
      spiketimes = 0.0
*
*
      rewind 1                  !go back to beginning of file
      do i = 1,skip             !read <skip> number of non-data records and discard them
         read (1,FMT=c_format,end=2000) CODE,CLOCK 
      end do                

 160  read (1,FMT=c_format,end=2000) CODE,CLOCK !read real data records
      if(CLOCK_FIRST.eq.0)CLOCK_FIRST=CLOCK !remember the 1st valid event time
      count=count+1
      if(MOD(count,100000).eq.0)then
         dots = dots + 1
         do i = 1, dots
            print '(''. '',$)'
         end do
         if(dots.eq.20)then
            dots = 0
            print '(//,''I''''m still working '')'
         end if
         print '()'
      end if

      if((CODE.eq.0).or.(CODE.gt.999))goto 160 !this is a blank record -- get another one
      if(ITAL(IDs(CODE)).lt.MAX_NUM_EVENTS)then
         ITAL(IDs(CODE))=ITAL(IDs(CODE))+1 !increase the tally of events that have been detected
*                                                for this ID code only if MAX_NUM_EVENTS has not 
*                                                been exceeded
         fl_clock = CLOCK       !change the integer CLOCK into double precision fl_clock

         SPIKETIMES(ITAL(IDs(CODE)),IDs(CODE))=fl_clock*factor !load this event into SPIKETIMES()
c         if((CODE.eq.38.or.CODE.eq.69).and.(ITAL(IDs(CODE)).le.5))
c     +         print '(I3,'': SPIKETIMES('',I1,'','',I3,'') = '',f)',
c     +         CODE,ITAL(IDs(CODE)),IDs(CODE),
c     +         SPIKETIMES(ITAL(IDs(CODE)),IDs(CODE))
      end if

      goto 160                  !read another record
*
*
*       **************************     RETURN     ************************************
*
*
 2000 do i = 1,MAX_NUM_CODES
         if(IDs(i).eq.0)cycle
         if(ITOTAL_TALLY(IDs(i)).gt.MAX_NUM_EVENTS)
     +        print '(/,T3,''***** WARNING *****''
     +''  Data array cannot hold all '',I8,
     +'' events for IDcode '',I3,
     +'' (MAX_NUM_EVENTS = '',I8,'')'')',
     +        ITOTAL_TALLY(IDs(i)),i,MAX_NUM_EVENTS

      end do

      endtime=(((CLOCK-CLOCK_FIRST)*factor)/1000.0D0)/60.0D0 !length of sample in minutes
      print '(//,T5,''Duration of '',A,'' = '',F10.2,'' minutes'',
     +//,T10,''Press <cr> to continue  >> '',$)',
     +     BDT_FILE(1:length),endtime
      read '(A)'
*
      RETURN
      END
      
      end module mod_rdat_universal
