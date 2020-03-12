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


      module mod_menus
      contains
      subroutine open_and_read(BDT_FILE,SPIKETIMES,
     +     ITAL,IDs,IDT,file_ext)
      use mod_miscellaneous_subroutines
      use mod_rdat_universal

      include 'x2000parameter.defs'
      include 'x2002array.defs'

      character*(*) BDT_FILE
      character*(*) file_ext
      logical input_exist
      double precision, allocatable, intent(out) :: spiketimes(:,:)

 5    print '(//,''The following input files are available:'',//)'
      isys=SYSTEM ('ls -C *.adt'//char(0))
      isys=SYSTEM ('ls -C *.bdt'//char(0))
*     isys=SYSTEM ('ls -C *.ddt'//char(0))
      isys=SYSTEM ('ls -C *.edt'//char(0))
      PRINT '(//,''Enter the name of the input file'',
     +     '' ...or... <cr> to exit  >> '',$)'
      READ (*,fmt='(A)',err=5) BDT_FILE
      if(BDT_FILE.eq.' ')then
         file_ext='x'
         return                 !return to main menu
      end if
      call remove_all_blanks(BDT_FILE,LEN(BDT_FILE)) !no blanks allowed in filename
      inquire (FILE=BDT_FILE,EXIST=input_exist)
      if(input_exist.eqv..TRUE.)then
         file_ext = BDT_FILE(INDEX(BDT_FILE,' ')-3:)
         OPEN(UNIT=1,file=BDT_FILE,FORM='FORMATTED',
     +        STATUS='OLD',err=5)
         CALL rdat6(SPIKETIMES,ITAL,IDs,IDT,BDT_FILE) !retrieve codes and times from input file
         close (unit=1)
         if(BDT_FILE.eq.'x')file_ext='x' !???? need this ???
      else                      !input data file does not exist, so ask again
         goto 5
      end if

      return
      end

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine enter_date_rec_protocol_name(date_exp,recording,
     +     protocol,exp_name,total_num_qdts,
     +     qdt_files,BDT_FILE)
      use mod_enter_date
      use mod_miscellaneous_subroutines

      include 'x2000parameter.defs'

      integer total_num_qdts,change_it
      character*(*) date_exp,recording,protocol,exp_name,
     +     qdt_files(MAX_NUM_QDTS),BDT_FILE
      character*11 default_date
      character*2 default_rec
      character*1 continu,offer_default
      
      if(.false.)print *,total_num_qdts !suppress unused variable warning

      date_exp=' '              !       date of experiment
      default_date=' '
      protocol=' '              !       experimental protocol
      recording=' '             !       recording#
      default_rec=' '
      exp_name=' '              !       name of experiment
      change_it=0

 9    print '(3(/))'
 15   PRINT '(/,''ENTER experiment name ''
     +     ''(175 characters; no commas) >> '',$)'
      exp_name=' '
      READ (*,fmt='(A)',err=15) exp_name
      if(exp_name.eq.' ')goto 15 !force an appropriate response
      do i = 1, 175
         if(exp_name(i:i).eq.',')exp_name(i:i)=';'
      end do
      offer_default = 'n'
      if((exp_name(1:1).ge.'0').and.(exp_name(1:1).le.'9'))then !if exp_name starts with a number...
         offer_default = 'y'
         default_date(1:2) = exp_name(9:10) ! we can use it to offer the user a
         if(exp_name(6:7).eq.'01')then ! default exp_date and recording
            default_date(4:6) = 'JAN'
         elseif(exp_name(6:7).eq.'02')then
            default_date(4:6) = 'FEB'
         elseif(exp_name(6:7).eq.'03')then
            default_date(4:6) = 'MAR'
         elseif(exp_name(6:7).eq.'04')then
            default_date(4:6) = 'APR'
         elseif(exp_name(6:7).eq.'05')then
            default_date(4:6) = 'MAY'
         elseif(exp_name(6:7).eq.'06')then
            default_date(4:6) = 'JUN'
         elseif(exp_name(6:7).eq.'07')then
            default_date(4:6) = 'JUL'
         elseif(exp_name(6:7).eq.'08')then
            default_date(4:6) = 'AUG'
         elseif(exp_name(6:7).eq.'09')then
            default_date(4:6) = 'SEP'
         elseif(exp_name(6:7).eq.'10')then
            default_date(4:6) = 'OCT'
         elseif(exp_name(6:7).eq.'11')then
            default_date(4:6) = 'NOV'
         elseif(exp_name(6:7).eq.'12')then
            default_date(4:6) = 'DEC'
         end if
         default_date(8:11) = exp_name(1:4)
         default_date(3:3) = '-'
         default_date(7:7) = '-'
         default_rec(1:2) = exp_name(13:14)
         if (default_rec(1:1) .eq. '0') then
            default_rec(1:1) = default_rec(2:2)
            default_rec(2:2) = ' '
         end if
      end if
      if(change_it.ne.0)goto 30 !user has modified a previous entry -
                                !  re-display info and ask for the OK
 20   if(offer_default.eq.'y')then
         PRINT '(/,''ENTER recording # with no leading zeros (default:
     +'',A2,'')  >> '',$)',default_rec
         READ (*,fmt='(A)',err=20) recording
         if(recording.eq.' ')recording=default_rec
      else
         PRINT '(/,''ENTER recording # with no leading zeros  >> '',$)'
         READ (*,fmt='(A)',err=20) recording
         if(recording.eq.' ')goto 20
      end if
      if(change_it.ne.0)goto 30 !user has modified a previous entry -
                                !  re-display info and ask for the OK
 10   if(offer_default.eq.'y')then
         PRINT '(/,''ENTER experiment date (default: '',A11,
     +        '')  >> '',$)',default_date
         READ (*,fmt='(A)',err=10) date_exp
         if(date_exp.eq.' ')date_exp=default_date
      else
         call enter_date(date_exp)
      end if
*     
      if(change_it.ne.0)goto 30 !user has modified a previous entry -
                                !  re-display info and ask for the OK

 25   PRINT '(/,''ENTER experimental protocol ''
     +     ''(200 characters; no commas) >> '',$)'
      READ (*,fmt='(A)',err=25) protocol
      do i = 1,200
         if(protocol(i:i).eq.',')protocol(i:i)=';'
      end do
      if(change_it.ne.0)goto 30 !user has modified a previous entry -
                                !  re-display info and ask for the OK

 30   call strlength(exp_name,LEN(exp_name),l_exp)
      call strlength(protocol,LEN(protocol),l_prot)
      call strlength(BDT_FILE,LEN(BDT_FILE),l_bdt)
      call strlength(protocol,LEN(protocol),l_prot)
      call strlength(qdt_files(1),LEN(qdt_files(1)),l_qdt1)
      call strlength(qdt_files(2),LEN(qdt_files(2)),l_qdt2)
      call strlength(qdt_files(3),LEN(qdt_files(3)),l_qdt3)
      call strlength(qdt_files(4),LEN(qdt_files(4)),l_qdt4)
      call strlength(qdt_files(5),LEN(qdt_files(5)),l_qdt5)
      PRINT '(30(/),5X,70(''_''),/,5X,''|'',T75,''|'',
     +     /,5X,''| 1. Name of experiment:  '',A,T75,''|'',
     +     /,5x,''| 2. Recording #:  '',A,T75,''|'',
     +     /,5X,''| 3. Date of experiment:  '',A,T75,''|'',
     +     /,5X,''| 4. Protocol:  '',A,T75,''|'',
     +     /,5x,''|'',T75,''|'',
     +     /,5X,''|    Input data filename:  '',T45,A,T75,''|'',
     +     /,5x,''|'',T75,''|'',
     +     /,5X,''|    Output data filename(s):  '',T45,A,T75,''|'',
     +     /,5X,''|'',T45,A,T75,''|'',
     +     /,5X,''|'',T45,A,T75,''|'',
     +     /,5X,''|'',T45,A,T75,''|'',
     +     /,5X,''|'',T45,A,T75,''|'',
     +     /,5X,70(''-''))',
     +     exp_name(1:l_exp),recording,date_exp,protocol(1:l_prot),
     +     BDT_FILE(1:l_bdt),
     +     qdt_files(1)(1:l_qdt1),
     +     qdt_files(2)(1:l_qdt2),
     +     qdt_files(3)(1:l_qdt3),
     +     qdt_files(4)(1:l_qdt4),
     +     qdt_files(5)(1:l_qdt5)

 40   print '(//,T10,''ACCEPT these data?  (y/n)  >> '',$)'
      read (*,fmt='(A1)',err=40) continu
      call upper_case(continu,LEN(continu))
      if(.false.)then
      else if(continu.eq.'Y')then
         return
      else if(continu.eq.'N')then
         change_it=0
 41      print '(/,T2,''Enter the number of the item you wish ''
     +        ''to change or <cr> to change all  >> '',$)'
         read (*,fmt='(I12)',err=41) change_it
*     
         if(.false.)then
         else if(change_it.eq.1)then !change name of experiment
            goto 15
         else if(change_it.eq.3)then !change date of experiment
            date_exp = ' '
            goto 10
         else if(change_it.eq.2)then !change recording number
            goto 20
         else if(change_it.eq.4)then !change protocol comment
            goto 25
         else if(change_it.eq.0)then !change output filenames (0--> change all data)
            date_exp = ' '
            GOTO 9
         end if
*     
      else
         goto 30                !continu <> Y or N; force a valid response
      end if

      return
      end


*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine enter_qdt_filename(number,QDT,mode)
      use mod_miscellaneous_subroutines

      integer length,number
      character*(*) QDT
      character*(*) mode
      character*1 over
      logical*4 exist

 10   if(mode.eq.'vt')then      !just looking at TQDT file(s)
         print '(//,T5,''Enter:'',T21,
     +        ''the name of an existing TQDT file to ''
     +        ''VIEW previously generated data,'',/,T11,
     +        ''.. or ..  a new filename to CREATE a new TQDT file,'',
     +        /,T11,
     +        ''.. or ..  <cr> to return to the main menu'',//,T21,
     +        ''>> '',$)'
      else
         print '(/,T5,''ENTER QDT filename #'',I1,''...or... ''
     +        ''<cr> when all filenames have been entered.''
     +        /,T5,''The QDT extension will be appended if ''
     +        ''you do not include it.'',//,T10,''>> '',$)',
     +        number
      end if
      QDT=' '
      read (*,'(A)',err=10) QDT
      if(QDT.ne.' ')then
         call remove_all_blanks(QDT,LEN(QDT))
         call strlength(QDT,LEN(QDT),length)
         if(mode.eq.'vt')then   !create / view tqdt file
            if(length.ge.6)then
               if(QDT(length-4:length).ne.'.tqdt') then
                  QDT = QDT(1:length)//'.tqdt' !be sure that all qdt filenames end in '.tqdt'
                  length = length + 5
               end if
            else
               QDT = QDT(1:length)//'.tqdt'
               length = length + 5
            end if
         else                   !create / view qdt file
            if(length.ge.5)then
               if(QDT(length-3:length).ne.'.qdt') then
                  QDT = QDT(1:length)//'.qdt' !be sure that all qdt filenames end in '.qdt'
                  length = length + 5
               end if
            else
               QDT = QDT(1:length)//'.qdt'
               length = length + 5
            end if
         end if
      else
         return
      end if
      if((mode.eq.'vt').or.(mode.eq.'cr'))return

 35   inquire(FILE='x2000/'//qdt(1:length),EXIST=exist)
      if((exist.eqv..TRUE.).and.(mode.eq.'cr'))then !warn user if about to overwrite a file
 36      over=' '
         print '(//,T10,''!!!!! WARNING !!!!!'',
     +        ''  File '',A,'' already exists!'',
     +        //,15X,''O -- write over existing QDT file'',
     +        /, 15X,''R -- rename the QDT file'',
     +        //,15X,''>> '',$)',
     +        QDT(1:length)
         READ (*,fmt='(A)',err=36) over

         call upper_case(over,LEN(over))
         if(over.eq.'R')then
            goto 10
         else if(over.eq.'O')then
            isys=SYSTEM ('rm x2000/'//QDT//char(0)) !remove old QDT file
            INQUIRE(FILE='x2000/'//QDT(1:length)//'.sav',EXIST=exist)
            if(exist.eqv..TRUE.)then
               isys=SYSTEM ('rm x2000/'//QDT(1:length)//'.sav' !remove old QDTSAV file
     +              //char(0))
            end if
            goto 35
         else                   !force an appropriate response
            goto 36
         end if
      end if
      if((exist.eqv..FALSE.).and.(mode.eq.'ed'))then !QDT file must exist if editing or viewing
         print '(''File '',A,'' does not exist in this directory.''
     +        ''  Please re-enter the qdt filename.'')',
     +        '"x2000/'//QDT(1:length)//'"'
         goto 10
      end if

      return
      end


*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine enter_binwidths(i,QDT,BWs,default,file_ext)

      include 'x2000parameter.defs'

      real BWs(MAX_NUM_QDTS,4)
      integer i
      character*(*)QDT
      character*(*) file_ext
      character*1 OK,default

      if(.false.)print *,qdt    !suppress unused variable warning

      do j = 1,4
         BWs(i,j)=0.0           !clear the array subset
      end do
      if(i.eq.1)default = 'n'

 113  if(default.eq.'n')then
         PRINT '(/,T5,''ENTER binwidths for calculation of CCHs ''
     +        ''(smallest --> largest) ...OR... '')'
         if(file_ext.eq.'bdt')then
            print '(T8,''<cr> for default binwidths of ''
     +           ''0.5, 1.5, 2.5, 5.5 ms:'')'
         else if(file_ext.eq.'edt')then
            print '(T8,''<cr> for default binwidths of ''
     +           ''0.1, 0.5, 1.5, 2.5 ms:'')'
         end if
         print '(/,T10,''Enter binwidth #1 (fmt=f5.1) ''
     +        ''(<cr>=default)  >> '',$)'
c     read (*,fmt='(F5.1)',err=114),BWs(i,1)
         read (*,fmt='(F5.1)',err=113)BWs(i,1)
         if(BWs(i,1).eq.0.0)default = 'y'
      else if(default.eq.'y')then
         PRINT '(/,T5,''ENTER binwidths for calculation of CCHs ''
     +        ''(smallest --> largest) ...OR... '')'

         if(file_ext.eq.'bdt')then
            print '(T8,''<cr> for default binwidths of ''
     +           ''7.5, 10.5, 12.5, 15.5 ms:'')'
         else if(file_ext.eq.'edt')then
            print '(T8,''<cr> for default binwidths of ''
     +           ''5.5, 7.5, 10.5, 15.5  ms:'')'
         end if
         print '(/,T10,''Enter binwidth #1 (fmt=f5.1) ''
     +        ''(<cr>=default)  >> '',$)'
         read (*,fmt='(F5.1)',err=113)BWs(i,1)
         if(BWs(i,1).eq.0.0)default = 'b'
      else
         PRINT '(/,T5,''ENTER binwidths for calculation of CCHs ''
     +        ''(smallest --> largest) '')'
         print '(/,T10,''Enter binwidth #1 (fmt=f5.1)  >> '',$)'
         read (*,fmt='(F5.1)',err=113)BWs(i,1)
      end if

*     
      if(BWs(i,1).eq.0.0)then   !go with the default binwidths
         if(file_ext.eq.'bdt')then
            if(default.eq.'y')then
               BWs(i,1)=0.5
               BWs(i,2)=1.5
               BWs(i,3)=2.5
               BWs(i,4)=5.5
            else if(default.eq.'b')then
               BWs(i,1)=7.5
               BWs(i,2)=10.5
               BWs(i,3)=12.5
               BWs(i,4)=15.5
            end if
         else if(file_ext.eq.'edt')then
            if(default.eq.'y')then
               BWs(i,1)=0.1
               BWs(i,2)=0.5
               BWs(i,3)=1.5
               BWs(i,4)=2.5
            else if(default.eq.'b')then
               BWs(i,1)=5.5
               BWs(i,2)=7.5
               BWs(i,3)=10.5
               BWs(i,4)=15.5
            end if
         end if
      else
 115     print '(T10,''Enter binwidth #2  >> ''$)'
         read (*,fmt='(F5.1)',err=115)BWs(i,2)
 116     print '(T10,''Enter binwidth #3  >> ''$)'
         read (*,fmt='(F5.1)',err=116)BWs(i,3)
 117     print '(T10,''Enter binwidth #4  >> ''$)'
         read (*,fmt='(F5.1)',err=117)BWs(i,4)
 118     OK=' '
         print '(/,T10,''Enter X to change these binwidths ''
     +        ''..or.. <cr> to continue  >> '',$)'
         read '(A)',OK
         if((OK.eq.'x').or.(OK.eq.'X'))goto 113
         if(OK.ne.' ')goto 118  !force a valid response
         if(file_ext.eq.'bdt')then
            do j=1,4
               if (BWs(i,j) < .5) BWs(i,j) = .5
               BWs(i,j) = ANINT (BWs(i,j) * 2) / 2
            end do
         end if
      end if

      return
      end

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine enter_database_filename(DB,mode)
      use mod_miscellaneous_subroutines

*     do not allow the user to:
*     over-write existing DBfiles when creating or
*     attempt to create new DBfiles when editing

      character*(*) DB
      character*(*) mode
      character*1 over
      integer length
      logical*4 exist

      print '(//,''The following ''''datasave'''' files exist'',
     +     '' in this directory:'',/)'
      icd=CHDIR('x2000'//char(0))
      isys=SYSTEM ('ls -C *.db.sav | sed s/.db.sav//g'//char(0))
      icd=CHDIR('..'//char(0))
 40   if(mode.eq.'cr')then
         print'(/,''Enter the base filename for a NEW datasave ''
     +        ''file'',
     +        ''[lower case ONLY; no extension; no quotes; ''
     +        ''20 char. max.]'',//,T10,''  >> '',$)'
      else                      !mode = 'ed' or 'jl'
         print '(/,''Enter the base filename for an EXISTING ''
     +        ''datasave file '',/,
     +        ''(no extension; lower case ONLY; no quotes; ''
     +        ''20 char. max.)'',
     +        //,T10,''  >> '',$)'
      end if
      read (*,fmt='(A)',err=40)DB
      call remove_all_blanks(DB,LEN(DB))
      call strlength(DB,LEN(DB),length)

      if(DB.eq.' ')then
         return                 !user changed mind -- return to calling code
      end if

      inquire(FILE='x2000/'//DB(1:length)
     +     //'.pre_dbg',EXIST=exist)
      if((exist.eqv..TRUE.).and.(mode.eq.'cr'))then !warn user if about to overwrite a file
         over=' '
         print '(//,T10,''!!!!! WARNING !!!!!'',
     +        ''  File '',A,'' already exists!'',
     +        ''  Enter a NEW filename!'')',
     +        DB(1:length)
         goto 40
      end if
      if((exist.eqv..FALSE.).and.(mode.ne.'cr'))then !DB file must exist if editing or viewing
         print '(''These datasave files do not exist in this ''
     +        ''directory.  ''
     +        ''Please re-enter the base filename.'',$)'
         goto 40
      end if

      return
      end


*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine enter_info_filename(FILE,mode,exp_name,
     +     recording,QDTSAV)
      use mod_miscellaneous_subroutines
      use mod_read_and_write_DBSAV

      include 'x2000parameter.defs'

      integer IDs(MAX_NUM_CODES),
     +     included_in_QDT(MAX_NUM_CODES),
     +     excluded_from_QDT(MAX_NUM_CODES),
     +     total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     total_histograms,
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     ITAL(MAX_NUM_CHAN),TOTAL_NUM_SHIFTS
      integer*4 SYSTEM,ios
      real BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,start_time,
     +     end_time,coefnum(MAX_NUM_CHAN)

      character*200 BDT_FILE
      character*200 FILE_1,dir_name,protocol
      character*50 QDT_FILENAME
      character*11 date_exp
      character*10 meanISI(MAX_NUM_CHAN),sdISI(MAX_NUM_CHAN),
     +     fiveHT(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT
      character*(*) FILE,exp_name,mode,recording,QDTSAV
      CHARACTER*5 coef(MAX_NUM_CHAN),ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN),DELTA2(MAX_NUM_CHAN)
      character*3 zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),zmodsig2_6(MAX_NUM_CHAN),
     +     card_type(MAX_NUM_CHAN)
      character*2 QDT_version,tedfactor(MAX_NUM_CHAN)
      character*1 OK,OK1,BOUNDARY
      logical*4 exist,exist_2
      external SYSTEM

      FILE = ' '
 1    print '(20(/),''Do you want to import a ''
     +     ''coordinate/AA/STA/comment info data file?  ''
     +     ''(y/n)  >> '',$)'
      read (*,'(A)') OK
      call upper_case(OK,LEN(OK))
      if(OK.eq.'Y')then
         if(mode.eq.'ed')then
            icd=CHDIR('x2000'//char(0))
            call read_QDTSAV(QDT_version,QDTSAV,ios,
     +           date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +           IDs,excluded_from_QDT,included_in_QDT,
     +           total_num_cells,
     +           I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +           BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +           BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +           end_time,icycles,
     +           total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +           ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +           zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +           zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +           zmodsig2_6,coef,coefnum,card_type,exp_name,
     +           DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,
     +           sd_rISI,
     +           num_rej_ISI,num_rej_rISI,c_MAX_INT,
     +           TOTAL_NUM_SHIFTS,num_acc_cycles,sp_per_cycle,ITAL)
            icd=CHDIR('..'//char(0))
         end if
         call strlength(exp_name,LEN(exp_name),l)
         if(l.eq.3)then         !might be an "old" experiment name (e.g., k15) -
            if((((exp_name(1:1).ge.'A').and.(exp_name(1:1).le.'Z')) !so append "m" and the recording number
     +           .or.((exp_name(1:1).ge.'a').and.
     +           (exp_name(1:1).le.'z')))
     +           .and.((exp_name(2:2).ge.'0').and.
     +           (exp_name(2:2).le.'9'))
     +           .and.((exp_name(3:3).ge.'0').and.
     +           (exp_name(3:3).le.'9')))
     +           dir_name=exp_name(1:3)//'m'//recording(1:1)
         else                   !this is a "new" experiment name (YYYY-MM-DD_xxx)
            dir_name = exp_name
         end if
         call strlength(dir_name,LEN(dir_name),l)
c     print '(''dir_name: '',A)', dir_name

         FILE_1 = '/oberon/experiments/'//dir_name(1:l)//'/'//
     +        dir_name(1:l)//'_info_curr.xls'
         inquire(FILE=FILE_1,EXIST=exist)
         if(exist.eqv..TRUE.)then
 5          print '(/,''The following coordinate/AA/STA/comment'',
     +           '' info data file exists on OBERON:  '',//,T10,A,
     +           /,T5,''Is this the one you want?  (y/n)  >> ''$)',
     +           FILE_1
            read (*,'(A)') OK1
            call upper_case(OK1,LEN(OK1))
            if(OK1.eq.'Y')then
               call strlength(FILE_1,LEN(FILE_1),m)
               isys=SYSTEM('XLS2CSV.pl '//FILE_1(1:m)//
     +              ' > x2000/info.tmp'//char(0))
               FILE = 'info.tmp'
               return
            else if(OK1.eq.'N')then
            else
               goto 5
            end if
         end if

         if((exist.eqv..FALSE.).or.(OK1.eq.'N'))then
            if(OK1.eq.'N')print '(/,''All righty, then -- ''
     +           ''go find it yourself.  Just be sure ''
     +           ''it''''s a *.info_curr.xls file!'',10(/))'
            if(exist.eqv..FALSE.)
     +           print '(/,''Fetch the *.info_curr.xls file.  Be ''
     +           ''sure it''''s a *.info_curr.xls file!'',10(/))'
            isys=SYSTEM('retrieve_info.tcl x2000/info.tmp'
     +           //char(0)//char(0))
            inquire(FILE='x2000/info.tmp',EXIST=exist_2)
            if(exist_2.eqv..TRUE.)then
               FILE = 'info.tmp'
            else
               FILE = ' '
            end if

            return
         end if
      else if(OK.eq.'N')then
         return
      else
         goto 1                 !force an appropriate answer
      end if
      end

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine enter_per_filename(FILE,mode)
      use mod_miscellaneous_subroutines

      integer length
      integer*4 SYSTEM
      character*200 FILE_1
      character*(*) FILE
      character*(*) mode
      character*1 OK
      logical*4 exist
      external SYSTEM


      if(.false.)print *,mode   !suppress unused variable warning

      FILE = ' '
 1    print '(20(/),''Do you want to import a ''
     +     ''response-to-perturbation data file?  ''
     +     ''(y/n)  >> '',$)'
      read (*,'(A)') OK
      call upper_case(OK,LEN(OK))
      if(OK.eq.'Y')then
 5       print '(/,''The following response-to-perturbation'',
     +        '' data file(s) exists in this directory:'',/)'
         icd=CHDIR('x2000'//char(0))
         isys=SYSTEM ('ls -C *.per.db'//char(0))
         icd=CHDIR('..'//char(0))
         print '(/,T5,''PICK ONE of these ..or.. ''
     +        ''<cr> to cancel  >> ''$)'
         read (*,'(A)',err=5) FILE_1
         if(FILE_1.eq.' ')then  !user entered <cr>, so CANCEL the file import
            FILE = ' '
            return
         else                   !user manually entered the name of a file in the x2000 dir.
            call remove_all_blanks(FILE_1,LEN(FILE_1))
            call strlength(FILE_1,LEN(FILE_1),length)
            if(FILE_1(length-6:length).ne.'.per.db')
     +           FILE_1 = FILE_1(1:length)//'.per.db' !be sure that all per filenames end in '.per.db'
            inquire(FILE='x2000/'//FILE,EXIST=exist) !be sure that the per.db file exists in the x2000 directory
            if(exist.eqv..FALSE.)then !PER file must exist if importing
               print '(''File '',A,'' does not exist in this ''
     +              ''directory.''
     +              ''Please re-enter the filename.'')',
     +              FILE_1(1:INDEX(FILE_1,' ')-1)
               goto 5
            else
               FILE = FILE_1
               return
            end if
         end if
      else if(OK.eq.'N')then    !user does not want to import a per.db file
         FILE = ' '
         return
      else
         goto 1                 !force an appropriate answer
      end if
      end

      end module mod_menus
      program xanalysis
      use mod_analyze_data
      use mod_calc_and_write_newshift_cccs
      use mod_database_write_new
      use mod_menus
      use mod_miscellaneous_subroutines
      use mod_qdttxt_write
*     
*     filename = menus.f
*     
*     mar-2007        lss
*     number of respiratory cycles that can be used to calculate
*     CTHs increased to 10,000 from 500; made changes to QDTSAV
*     files to accommodate the change from I3 to I5 (acutally,
*     formatted to I6 to be long-sighted)
*     apr-2002        lss
*     beginnings of handling of ddt files!  Also, allow import of coordinate data from a text file.
*     
*     23-Jan-02       lss
*     complete re-working of initial menu to simplify choices in
*     an attempt to offer a more complete analysis package;
*     autoACH can now be called from this menu
*     
*     *****  This program is the "gateway" into the pairs analysis    *****
*     *****   program (working name = x2000).  In conjunction with    *****
*     *****   the startup shell (st_x2000), it allows multiple users  *****
*     *****   to run the program -- a user-specific copy of the       *****
*     *****   is created (eg.,x2000_<YOUR LOGIN NAME HERE>) each time *****
*     *****   x2000 is run.  When the user exits the program, this    *****
*     *****   user-specific copy is deleted.  Likewise, the graphics  *****
*     *****   windows are uniquely named.                             *****
*     
*     *****  Overall, this program provides an improvement over the   *****
*     *****   XAUTOX series of analysis programs in that it reduces   *****
*     *****   "wait time".  All CCHs are calculated according to      *****
*     *****   user-entered parameters and are written to a file that  *****
*     *****   may be opened at a later time for viewing and entry     *****
*     *****   of analysis results.  Output files may be re-opened and *****
*     *****   appended.                                               *****
*     
*     
*     
*     
*     
      SAVE                      !forces arrays in this pgm unit to be stored in static memory
      include 'x2000parameter.defs'

      double precision, allocatable :: SPIKETIMES(:,:)
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)

      integer*4 isys,icd,total_num_qdts,SYSTEM,CHDIR

      real BWs(MAX_NUM_QDTS,4)
*     
      logical exist,exist_old,exist_new
*     
      character*220 DBSAV
      character*200 BDT_FILE,protocol,DB_FILES
      character*200 INFO_FILENAME,PER_FILENAME
      character*175 exp_name
      character*60 QDTTXT_old,QDTTXT_new
      character*50 qdt_files(MAX_NUM_QDTS)
      character*60 QDTSAV
      character*20 c_format,cntrl
      character*11 date_exp
      character*3 file_ext
      character*2 version,version1,mode,recording
      character*1 abort,import_info,import_per,
     +     default,choice1
*     
      EXTERNAL CHDIR,SYSTEM
*     *****************************************************************
*     *****************************************************************
*     
*     
      version = '8'             !version 7 instituted in March-2007
      call getenv ("XANALYSIS_CONTROL",cntrl)
      IF((trim(cntrl).ne.'').and.(cntrl.NE.'cth_cch'))version = '7'

*     
      if(version(1:1).lt.'3')then !format of the QDT file varies with the  program version
         c_format = '(109I6)'   !version 2 or earlier
      else
         c_format = '(2I6,107I10)' !version 3 or higher
      end if

      isys=SYSTEM('mkdir x2000 2>/dev/null'//char(0)) !create a subdirectory named "x2000" which
                                !  will hold qdt and database files
                                ! (2>/dev/null sends any error messages into
                                !   UNIX's "black hole" so the user won't see it)
      isys=SYSTEM('mkdir x2000/backups 2>/dev/null'//char(0)) !create a subdirectory of x2000 named "backups" which
                                !will hold numbered, zipped backups of database files
      goto 10                   !MAIN_MENU                                       !skip this introductory text unless user chooses option R

 1    print '(20(/),T10,''Welcome to ''$)'
      call version_message
      print '(//,T2,'' ******* F.Y.I. *******'',
     +     //,T2,''* xAnalysis may be used to perform initial ''
     +     ''data surveys or a more in-depth analysis ''
     +     /,T2,''  which will produce database-ready output files.''
     +     /,T2,''* Additional programs which can be run from ''
     +     ''xAnalysis''''s main menu include ''
     +     /,T2,''  autoACH, autoCCH, and autoSTA.''
     +     /,T2,''  These programs individually read the input file, ''
     +     ''provide limited statistical analysis and ''
     +     /,T2,''  do not produce an output file, but are useful for ''
     +     ''quick viewing of data.''
     +     /,T2,''* Menu option 5 results in the creation of ''
     +     ''ACHs, CTHs, and CCHs for all included cells;  ''
     +     /,T2,''  these histograms, as well as some statistical ''
     +     ''analyses, are saved in a TQDT and a QDTSAV file''
     +     /,T2,''  which may be reopened for later viewing.''
     +     /,T2,''* Menu option C produces a QDT file containing ''
     +     ''ACHs, CTHs, CCHs, and ''
     +     /,t2,''  shift-control CCHs for all included cells.''
     +     /,T2,''* Menu option E opens QDT and associated datasave ''
     +     ''files and ''
     +     /,T2,''  allows you to ''
     +     ''perform additional statistical analyses of CTHs, ''
     +     /,T2,''  expects additional input ''
     +     ''files (re cell location, pattern, response, etc.), ''
     +     /,T2,''  and records your descriptions of CCH features.''
     +     /,T2,''* Menu option W produces quote- and ''
     +     ''comma-delimited text files for importation ''
     +     ''into a database.''
     +     /,T2,''* QDT and datasave files are saved in subdirectory ''
     +     ''./x2000/''
     +     /,T2,''* This version of xAnalysis is limited '',
     +     ''to 200 channels and 450,000 events/channel.'',
     +     /,T2,''* IDcode numbers for spike trains must be '',
     +     ''between 1 and 999, inclusive.'',
     +     /,T2,''* BDT and EDT files may be used as input.''
     +     /,T2,''* A disorderly exit may produce incomplete output''
     +     '' files - use the menu to exit.'',
     +     /,T2,''* Suggestions?  Questions?  Complaints?  See Lauren.''
     +     //,T5,''Press <cr> to continue >>>>>  '',$)'
      read '(A)'

*     
*     ***** Main menu:        *****
*     
*     
 10   do i = 1, MAX_NUM_QDTS
         qdt_files (i)= ' '
      end do
      global_compare_and_convert_flag = ' '
      BDT_FILE = ' '
      QDTSAV = ' '
      DBSAV = ' '
      exp_name = ' '
      recording = ' '
      mode = ' '                !re-initialize the mode flag
      abort = 'n'
      import_info = 'n'
      import_per = 'n'
      INFO_FILENAME = ' '
      PER_FILENAME = ' '
      inquire(FILE='x2000/info.tmp',EXIST=exist)
      if(exist.eqv..TRUE.)isys=SYSTEM('rm x2000/info.tmp'//char(0))

      print '(20(/),''MAIN MENU for '',$)'
      call version_message
      print '(/,''------------------------------'',
     +     T48,''R   Read a brief introduction to Xanalysis'',//,
     +     T2, ''PRELIMINARY SCREENING:'',
     +     T48,''1   autoACH'',/,
     +     T6, ''NO FILES ARE CREATED by these '',
     +     T48,''2   autoCTH'',/,
     +     T7,''''''calculated-on-the-fly'''' routines:'',
     +     T48,''3   autoCCH'',/,
     +     T48,''4   autoSTA'',//,
     +     T6, ''PARTIAL HISTOGRAM FILES:'',
     +     T48,''5   CREATE new or VIEW existing TQDT file'',//,
     +     T2, ''MAIN DATA ANALYSIS and DATASAVE'',
     +     T48,''C   CREATE new set of datasave / QDT files'',/,
     +     T2, ''  FILE MAINTENANCE:'',
     +     T48,''E   EDIT an existing set of datasave '',
     +     ''/ QDT files'',/,
     +     T48,''V   VIEW an existing set of datasave '',
     +     ''/ QDT files'',/,
     +     T48,''W   WRITE a set of database-ready tables'',/,
     +     T48,''      to ./x2000, /oberon/experiments/, and '',/,
     +     T48,''      /oberon/databases/GAIA database/'',/,
     +     T2, ''OTHER PROGRAMS:'',
     +     T48,''M   matrix'',//,
     +     T48,''X   EXIT'',/,
     +     T30,        '' >> '',$)'

      read (*,fmt='(A1)',err=10) choice1
      call upper_case(choice1,LEN(choice1))


*     ************ case 1: view ACHs ************

      if(choice1.eq.'1')then    !run autoACH
c     if(file_ext.eq.'x')goto 10 !MAIN_MENU               !we have a problem - return to main menu
         isys=SYSTEM('autoACH'//char(0))
         goto 10                !MAIN_MENU                                  !return to main menu

*     ************ case 2: view CTHs ************

      else if(choice1.eq.'2')then !run autoCTH
c     if(file_ext.eq.'x')goto 10 !MAIN_MENU               !we have a problem - return to main menu
         isys=SYSTEM('autoCTH'//char(0))
         goto 10                !MAIN_MENU                                  !return to main menu

*     ************ case 3: view CCHs ************

      else if(choice1.eq.'3')then !run autoACH
c     if(file_ext.eq.'x')goto 10 !MAIN_MENU                 !we have a problem - return to main menu
         isys=SYSTEM('autoCCH'//char(0))
         goto 10                !MAIN_MENU                            !return to main menu

*     ************ case 4: view STAs ************

      else if(choice1.eq.'4')then !run autoSTA
         isys=SYSTEM('autoSTA'//char(0))
         goto 10                !MAIN_MENU

*     ************ case 5: create or view CTHS or CCHs in "partial" histogram file: (mode = 'vt') ************

      else if(choice1.eq.'5')then !might need to open an input data file (eg., bdt)

         mode = 'vt'
         print '(//,''The following .tqdt files exist in '',
     +        ''this directory:'',//)'
         icd=CHDIR('x2000'//char(0)) !change to the directory that holds tqdt files
         isys=SYSTEM ('ls -C *.tqdt'//char(0))
         icd=CHDIR('..'//char(0)) !return to the bdt file directory

         call enter_qdt_filename(1,qdt_files(1),mode)
         total_num_qdts = 1
         if(qdt_files(1).eq.' ')goto 10 !MAIN_MENU          !no tqdt filename entered

         inquire (FILE='x2000/'//qdt_files(1),EXIST=exist)
         if(exist.eqv..TRUE.)then
            goto 500            !ANALYZE                       !tqdt file exists
         else                   !tqdt file does not exist --> create it
            call open_and_read(BDT_FILE,SPIKETIMES,
     +           ITAL,IDs,IDT,
     +           file_ext)
            if(file_ext.eq.'x')goto 10 !MAIN_MENU                           !we have a problem - return to main menu
            call enter_binwidths(1,qdt_files(1),BWs,default,
     +           file_ext)
            call enter_date_rec_protocol_name(date_exp,recording,
     +           protocol,exp_name,total_num_qdts,
     +           qdt_files,BDT_FILE)
            goto 400            !CALCULATE
         end if


*     ************ case C: create new qdt and new datasave files (mode = 'cr')************

      else if(choice1.eq.'C')then
         mode='cr'              !create new datasave files
         call open_and_read(BDT_FILE,SPIKETIMES,
     +        ITAL,IDs,IDT,
     +        file_ext)
         if(file_ext.eq.'x')goto 10 !MAIN_MENU               !we have a problem - return to main menu

         print '(//,''The following .qdt files already exist in '',
     +        ''this directory:'',//)'
         icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
         isys=SYSTEM ('ls -C *.qdt'//char(0))
         icd=CHDIR('..'//char(0))
         call strlength(BDT_FILE,LEN(BDT_FILE),l)
         print '(//,''You may generate up to '',I2,
     +        '' histogram ''
     +        ''files from '',A,'' using ''
     +        ''different binwidths.'')',MAX_NUM_QDTS,
     +        BDT_FILE(1:l)
         total_num_qdts = 0

         do i = 1, MAX_NUM_QDTS
            qdt_files (i)= ' '
         end do
         do i = 1,MAX_NUM_QDTS
            call enter_qdt_filename(i,qdt_files(i),mode)
            if(qdt_files(i).ne.' ')then
               total_num_qdts=total_num_qdts+1 !keep track of how many total number of qdt files to be created
               call enter_binwidths(i,qdt_files(i),BWs,default,
     +              file_ext)
            else
               exit             !no more qdt filenames to be entered
            end if
         end do
         if(qdt_files(1).eq.' ')goto 10 !MAIN_MENU           !user changed mind - back to main menu

         call enter_database_filename(DB_FILES,mode)
         if(DB_FILES.eq.' ')goto 10 !MAIN_MENU        !user changed mind - back to main menu

         call enter_date_rec_protocol_name(date_exp,recording,
     +        protocol,exp_name,total_num_qdts,
     +        qdt_files,BDT_FILE)

         call enter_info_filename(INFO_FILENAME,mode,exp_name,
     +        recording,QDTSAV)
         if(INFO_FILENAME.ne.' ')import_info='y'

         call enter_per_filename(PER_FILENAME,mode)
         if(PER_FILENAME.ne.' ')import_per='y'

c     call enter_tsfs_filename(TSFS_FILENAME,mode)
c     if(TSFS_FILENAME.ne.' ')import_per='y'

         goto 400               !CALCULATE

*     ************ case E or V: edit or view existing datasave files (mode = 'ed' or 'jl') ************

      else if((choice1.eq.'E').or.(choice1.eq.'V'))then !choice1=E --> edit/append existing datasave files
*     !choice1=V --> view-only existing datasave files
         if(choice1.eq.'E')mode='ed' !edit data and save it in the gamesave files and database-ready tables
         if(choice1.eq.'V')mode='jl' !just look at the pretty histograms
*     
         print '(//,''The following *.qdt files exist in '',
     +        ''this directory:'',//)'

         icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
         isys=SYSTEM('ls -C *.qdt'//char(0))
         icd=CHDIR('..'//char(0))

         print '(//,''You may view up to '',I2,
     +        '' histogram files.'')',MAX_NUM_QDTS
         total_num_qdts = 0
         do i = 1, MAX_NUM_QDTS
            qdt_files (i)= ' '
         end do
         do i = 1,MAX_NUM_QDTS
            call enter_qdt_filename(i,qdt_files(i),mode)
            if(qdt_files(i).ne.' ')then
               call strlength(qdt_files(i),LEN(qdt_files(i)),l)
               QDTSAV = qdt_files(i)(1:l)//'.sav'
               total_num_qdts=total_num_qdts+1 !keep track of how many total number of qdt files to be viewed
               call strlength(qdt_files(i),LEN(qdt_files(i)),l)
               QDTTXT_old = qdt_files(i)(1:l)//'txt'
               QDTTXT_new = qdt_files(i)(1:l)//'.txt'
               icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
               inquire(FILE=QDTTXT_old,EXIST=exist_old)
               inquire(FILE=QDTTXT_new,EXIST=exist_new)
               if(exist_new.eqv..FALSE.)then
                  print '(/,'' Please wait a moment.  ''
     +                 ''I''''m writing the QDTTXT file ...'')'
                  call write_QDTTXT(qdt_files(i),0.0,0.0) !making sure that all QDTTXT files are current
               end if
               if(exist_old.eqv..TRUE.) !delete the old-style QDTTXT file
     +              isys=SYSTEM('rm '//QDTTXT_old//char(0))
               icd=CHDIR('..'//char(0))
            else
               exit             !no more qdt filenames to be entered
            end if
         end do
         if(qdt_files(1).eq.' ')goto 10 !MAIN_MENU           !user changed mind - go back to main menu

         call enter_database_filename(DB_FILES,mode)
         if(DB_FILES.eq.' ')then
            goto 10             !MAIN_MENU                                !user changed mind - go back to main menu
         else
            call concat(DB_FILES,'.db.sav',DBSAV,l)
            icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
            OPEN (14,FILE=DBSAV,FORM='FORMATTED',STATUS='OLD')
            read (14,'(A2)')version1
            close (14)
            icd=CHDIR('..'//char(0))
            if(version1(1:1).lt.'6')then !this version of xanalysis is only for versions 6 and up
               call strlength(DB_FILES,LEN(DB_FILES),l)
               print '(10(/),''You are trying to use datasave ''
     +              ''files that were created with v5 or ''
     +              ''earlier of XAnalysis.'',/,
     +              ''Run st_xa_converter to transform the '',A,
     +              '' files before using them in this ''
     +              ''version of XAnalysis.'',//,T10,
     +              ''<cr> to return to main menu  >> '',$)',
     +              QDTSAV(1:l)
               read (*,'(A)')
               goto 10          !MAIN_MENU
            end if
         end if

         if(mode.eq.'ed')then
            call strlength(qdt_files(1),LEN(qdt_files(1)),l)
            QDTSAV = qdt_files(1)(1:l)//'.sav'


            call enter_info_filename(INFO_FILENAME,mode,exp_name,
     +           recording,QDTSAV)
            if(INFO_FILENAME.ne.' ')import_info='y'
            call enter_per_filename(PER_FILENAME,mode)
            if(PER_FILENAME.ne.' ')import_per='y'
         end if
         goto 500

*     ************ case M: run matrix to make a text file of significant pairs (uses pre_dbp file) ************

      else if(choice1.eq.'M')then
         isys=SYSTEM('xmatrix'//char(0))
         print '(/,T10,''Press <cr> to continue  >> '',$)'
         read (*,'(A)')
         goto 10                !MAIN_MENU                                  !return to main menu

*     ************ case W: WRITE the database-ready tables that will imported into Access (or similar) database: ****

      else if(choice1.eq.'W')then
         mode = 'wr'
*     ***** what about including a checking routine to see if ready to write everything?!? *****
         call enter_database_filename(DB_FILES,mode)
         if(DB_FILES.eq.' ')goto 10 !MAIN_MENU     !user changed mind - go back to main menu
         icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
         call database_write(DB_FILES) !write the database-ready table files
         icd=CHDIR('..'//char(0)) !return to input file subdirectory
         goto 10                !MAIN_MENU

*     ************ case R: READ some introductory text:  ************

      else if(choice1.eq.'R')then !README
         goto 1

*     ************ case X: EXIT ************

      else if(choice1.eq.'X')then !EXIT -- end the program
         close (unit=2)         !be sure the qdt file is closed
         inquire(FILE='x2000/info.tmp',EXIST=exist)
         if(exist.eqv..TRUE.)isys=SYSTEM('rm x2000/info.tmp'//char(0))
         inquire(FILE='x2000/pert.tmp',EXIST=exist)
         if(exist.eqv..TRUE.)isys=SYSTEM('rm x2000/pert.tmp'//char(0))
         goto 999               !exit
*     
      else                      !force a valid response
         goto 10                !MAIN_MENU
      end if


*     ************ CALCULATE: ************

 400  icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files
      
      call calculate_histograms(version,mode,
     +     BDT_FILE,
     +     qdt_files,SPIKETIMES,IDs,ITAL,abort,c_format,
     +     total_num_qdts,BWs,date_exp,recording,protocol,
     +     exp_name,DB_FILES)

      icd=CHDIR('..'//char(0))  !return to input file subdirectory
      if(abort.eq.'y')then
         goto 10                !MAIN_MENU                        !abort the process - back to menu
      end if


*     ************ ANALYZE: ************

 500  icd=CHDIR('x2000'//char(0)) !change to the directory that holds qdt and datasave files

      if (mode .ne. 'jl') then
         OPEN (2,FILE=TRIM(DB_FILES)//'_per_filename',
     +        ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
         WRITE(2,'(A)')TRIM(PER_FILENAME)
         CLOSE(2)
      end if
c      print *,'Call analyze_data'
      call analyze_data(version,mode,
     +     BDT_FILE,
     +     qdt_files,
     +     DB_FILES,
     +     import_info,
     +     INFO_FILENAME,import_per,PER_FILENAME,
     +     c_format,total_num_qdts,date_exp,
     +     recording,protocol,exp_name)

      close(unit=2)
      icd=CHDIR('..'//char(0))  !return to input file subdirectory
      goto 10                   !MAIN_MENU                           !back to main menu
*     
*     
*     
 999  end



*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>

