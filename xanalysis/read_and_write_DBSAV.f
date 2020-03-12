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

      module mod_read_and_write_DBSAV
      contains
*     filename = read_and_write_DBSAV

*     date of last revision:  17-Mar-2006       lss

*     This module of Xanalysis (aka x2004) contains subroutines to read and write the DBSAV and QDTSAV files


      subroutine read_DBSAV(file_version,DBSAV,ios,included,excluded,
     +     analyzed_cells,analyzed_pairs,CELL_NAMES,perturb_applied,
     +     perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +     info_prev_imported,
     +     per_prev_imported,total_num_cells_db,qdtfilenames)

      include 'x2000parameter.defs'

      integer ios,total_num_cells_db,
     +     excluded(MAX_NUM_CODES),
     +     included(MAX_NUM_CODES),
     +     analyzed_cells(MAX_NUM_CHAN),
     +     analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +     perturb(MAX_PERTURB)

      character*(*) qdtfilenames
      character*(*) resp_type(MAX_NUM_CHAN)
      character*(*) CELL_NAMES(MAX_NUM_CHAN)
      character*(*) perturb_applied,AA_applied,STA_applied,
     +     per_prev_imported,info_prev_imported,AA(MAX_AA),
     +     STA(MAX_STA)
      character*(*) file_version
      character*200 qdtfilenames_200
      character*2 version
      character*(*) DBSAV
      character*151 fmt

      goto 10
 100  print '(''error opening DBSAV '',A,/,''error='',I5)',DBSAV,ios
      goto 999
 301  print '(''301err writing DBSAV '',A,/,''error='',I5)',DBSAV,ios
 302  print '(''302err writing DBSAV '',A,/,''error='',I5)',DBSAV,ios
 306  print '(''306err writing DBSAV '',A,/,''error='',I5)',DBSAV,ios
      goto 999
      
 10   OPEN (15,FILE=DBSAV,FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL',ERR=100,IOSTAT=ios)
      rewind 15
      read (15,'(A2)',ERR=11,IOSTAT=ios),version
 11   rewind 15
c        print '(''ios = '',I)',ios
c        print '(''version = '',A1)',version(1:1)
      if(version(1:1).eq.' ')then
c           print '(''read DBSAV the REALLY old way'')'
c           read (*,'(A)')
         MAX_NUM_CHAN_SQ = 120*120 !when MAX_NUM_CHAN = 120:

         write(fmt,'(A1,I12,A5,I12,A11,I12,A54)')
     +        '(',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),120I1,',
     +        MAX_NUM_CHAN_SQ,
     +        '(I1),120A4,A1,50I1,A1,10A1,A1,10A1,120A15,2A1,I6,A200)'

         read (15,fmt,err=306,IOSTAT=ios),
     +        included,excluded,
     +        (analyzed_cells(i),i=1,120),
     +        ((analyzed_pairs(i,j),i=1,120),j=1,120),
     +        (CELL_NAMES(i),i=1,120),perturb_applied,
     +        perturb,AA_applied,
     +        AA,STA_applied,STA,(resp_type(i),i=1,120),
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames_200
c           print '(''read DBSAV: total_num_cells_db = '',I)',
c     +          total_num_cells_db
         if(total_num_cells_db.le.0)then
            total_num_cells_db = 0
            do j = 1, MAX_NUM_CODES
               if(included(j).ne.0)
     +              total_num_cells_db = total_num_cells_db + 1
c                 print '(''read DBSAV: total_num_cells_db = '',I)',
c     +                total_num_cells_db
            end do
         end if
         read (*,'(A)')
      elseif(version(1:1).lt.'6')then
c           print '(''read DBSAV the old way'')'
c           read (*,'(A)')
         MAX_NUM_CHAN_SQ = 120*120 !when MAX_NUM_CHAN = 120:
         read (15,'(A2,999(I3),999(I3),               
     +120I1,14400(I1),
     +120A4,A1,50I1,A1,10A1,A1,10A1,
     +120A15,2A1,I6,A200)',err=301,IOSTAT=ios),
     +        file_version,included,excluded,
     +        (analyzed_cells(i),i=1,120),
     +        ((analyzed_pairs(i,j),i=1,120),j=1,120),
     +        (CELL_NAMES(i),i=1,120),perturb_applied,
     +        perturb,AA_applied,
     +        AA,STA_applied,STA,(resp_type(i),i=1,120),
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames_200
         qdtfilenames = qdtfilenames_200
      elseif(version(1:1).ge.'6')then
c           print '(''read DBSAV the new way'')'
c           read (*,'(A)')
         MAX_NUM_CHAN_SQ = MAX_NUM_CHAN**2
         write(fmt,'(A4,I12,A5,I12,A5,I12,A5,I12,A5,I12,A35,I12,A20)')
     +        '(A2,',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CHAN,
     +        '(I1),',
     +        MAX_NUM_CHAN_SQ,
     +        '(I1),',
     +        MAX_NUM_CHAN,
     +        '(A4),A1,50(I1),A1,10(A1),A1,10(A1),',
     +        MAX_NUM_CHAN,
     +        '(A15),2(A1),I6,A250)'
c     format string length = 151

         read (15,fmt,err=302,IOSTAT=ios),
     +        file_version,included,excluded,analyzed_cells,
     +        analyzed_pairs,
     +        CELL_NAMES,perturb_applied,perturb,AA_applied,
     +        AA,STA_applied,STA,resp_type,info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames
      end if

      close (unit=15)

 999  return
      end


*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      subroutine write_DBSAV(pgm_version,DBSAV,ios,included,
     +     excluded,
     +     analyzed_cells,analyzed_pairs,CELL_NAMES,perturb_applied,
     +     perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +     info_prev_imported,
     +     per_prev_imported,total_num_cells_db,qdtfilenames)

      include 'x2000parameter.defs'

      integer ios,total_num_cells_db,
     +     excluded(MAX_NUM_CODES),
     +     included(MAX_NUM_CODES),
     +     analyzed_cells(MAX_NUM_CHAN),
     +     analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +     perturb(MAX_PERTURB)

      character*(*) qdtfilenames,pgm_version
      character*(*) resp_type(MAX_NUM_CHAN)
      character*(*) CELL_NAMES(MAX_NUM_CHAN)
      character*(*) perturb_applied,AA_applied,STA_applied,
     +     per_prev_imported,info_prev_imported,AA(MAX_AA),
     +     STA(MAX_STA)
      character*(*) DBSAV

      character*333 fmt

c        print '(''write DBSAV: pgm_version = '',A)',pgm_version


      if(global_mode.eq.'jl'.and.global_compare_and_convert_flag.eq.' ')
     +     then
         print '
     +(''bug: global_compare_and_convert_flag is blank in write_preDB
     +SAV'')'
         stop
      end if
      
      if (global_mode.eq.'jl'.and.global_compare_and_convert_flag
     +     .eq.'n') return

      MAX_NUM_CHAN_SQ = MAX_NUM_CHAN**2

      goto 10
 100  print '(''error opening DBSAV '',A,/,''error='',I5)',DBSAV,ios
      goto 999
      print '(''207err writing DBSAV '',A,/,''error='',I5)',DBSAV,ios
      goto 999
      
 10   OPEN (15,FILE=DBSAV,FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL',ERR=100,IOSTAT=ios)
      rewind 15

      write(fmt,'(A4,I12,A5,I12,A5,I12,A5,I12,A5,I12,A35,I12,A20)')
     +     '(A2,',
     +     MAX_NUM_CODES,
     +     '(I3),',
     +     MAX_NUM_CODES,
     +     '(I3),',
     +     MAX_NUM_CHAN,
     +     '(I1),',
     +     MAX_NUM_CHAN_SQ,
     +     '(I1),',
     +     MAX_NUM_CHAN,
     +     '(A4),A1,50(I1),A1,10(A1),A1,10(A1),',
     +     MAX_NUM_CHAN,
     +     '(A15),2(A1),I6,A250)'
c     format string length = 151

      write (15,fmt,IOSTAT=ios)
     +     pgm_version,included,excluded,analyzed_cells,
     +     analyzed_pairs,
     +     CELL_NAMES,perturb_applied,perturb,AA_applied,
     +     AA,STA_applied,STA,resp_type,info_prev_imported,
     +     per_prev_imported,total_num_cells_db,qdtfilenames

      close (unit=15)

 999  return
      end


*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      subroutine read_QDTSAV(file_version,QDTSAV,ios,
     +     date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded,included,total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +     end_time,icycles,
     +     total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +     ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +     zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +     zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +     zmodsig2_6,coef,coefnum,card_type,exp_name,
     +     DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +     num_rej_ISI,num_rej_rISI,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +     num_acc_cycles,sp_per_cycle,ITAL)

      include 'x2000parameter.defs'
      
      integer IDs(MAX_NUM_CODES),excluded(MAX_NUM_CODES),
     +     included(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),
     +     total_num_cells,I_pulse,
     +     E_pulse,phrenic,BNDRY,cardiac_pls,icycles,
     +     total_histograms,TOTAL_NUM_SHIFTS,num_acc_cycles,
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     SPC_array

      real BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,start_time,end_time,
     +     coefnum(MAX_NUM_CHAN)

      character*(*) exp_name,protocol,BDT_FILE,QDT_FILENAME,QDTSAV
      character*30 BDT_FILE_30
      character*20 exp_name_20,protocol_20,QDT_FILENAME_20
      character*11 date_exp
      character*10 meanISI(MAX_NUM_CHAN),sdISI(MAX_NUM_CHAN),
     +     fiveHT(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT
      character*5 ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN),coef(MAX_NUM_CHAN),
     +     DELTA2(MAX_NUM_CHAN)
      character*3 zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),zmodsig2_6(MAX_NUM_CHAN),
     +     card_type(MAX_NUM_CHAN)
      character*2 file_version,recording,tedfactor(MAX_NUM_CHAN),
     +     version
      character*1 BOUNDARY
      character*333 fmt

c      print *,'begin read_qdtsav'
      goto 10
 100  print '(''error opening QDTSAV '',A,/,''error='',I5)',
     +     QDTSAV,ios
      goto 999
 408  print '(''408err reading QDTSAV '',A,/,''error='',I5)',
     +     QDTSAV,ios
 409  print '(''409err reading QDTSAV '',A,/,''error='',I5)',
     +     QDTSAV,ios
c        print '(''total # histograms = '',I10)',total_histograms
c        do i = 1, 10
c           print '(I3,'':  '',$)',i
c           print '(''ETA2_1 = '',A5)',ETA2_1(i)
c           print '(T5,''zmodsig1_6 = '',A3)',zmodsig_6(i)
c           print '(T5,''zmodsig2_6 = '',A3)',zmodsig2_6(i)
c           print '('' coef = '',A5,''; coefnum = '',f7.2)',
c     +          coef(i),coefnum(i)
c           print '(T5,''<cr>  >> '',$)'
c           read (*,'(A)')
c        end do
      goto 999

 10   OPEN (14,FILE=QDTSAV,FORM='FORMATTED',STATUS='OLD',
     +     ERR=100,IOSTAT=ios)
      rewind 14
      read (14,'(A2)',ERR=408,IOSTAT=ios),version
      file_version = version
      rewind 14
c      print *,'file_version = ', file_version
      if(version(1:1).lt.'6')then !for versions 1 thru 5: MAX_NUM_CHAN = 120
c           print '(''read QDTSAV the old way'')'                !                       char*30 BDT_FILE
         SPC_array = 120 * MAX_NUM_ACC_CYCLES

         write(fmt,'(A23,I12,A5,I12,A5,I12,A247,I12,A5)')
     +        '(A2,A11,A2,A20,A30,A20,',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),6I3,2f7.1,2f9.1,5f7.1,A1,2f5.1,I3,I10,120A5,120A5,
     +120A5,120A5,120A5,120A5,120A3,120A3,120A3,120A3,120A3,
     +120A3,120A3,120A3,120A3,120A3,120A3,120A3,120A5,120f7.2,
     +120A3,A20,120A5,120A2,120A10,120A10,120A10,120A10,
     +120A10,120A10,120A10,A10,I3,I10,',
     +        SPC_array,
     +        '(I5))'
c     format string length = 333

         write(fmt,'(A23,I12,A5,I12,A5,I12,A247,I12,A5)')
     +        '(A2,A11,A2,A20,A30,A20,',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),',
     +        MAX_NUM_CODES,
     +        '(I3),6I3,2f7.1,2f9.1,5f7.1,A1,2f5.1,I3,I10,120A5,120A5,
     +120A5,120A5,120A5,120A5,120A3,120A3,120A3,120A3,120A3,
     +120A3,120A3,120A3,120A3,120A3,120A3,120A3,120A5,120f7.2,
     +120A3,A20,120A5,120A2,120A10,120A10,120A10,120A10,
     +120A10,120A10,120A10,A10,I3,I10,',
     +        SPC_array,
     +        '(I5))'
c     format string length = 333

         read (14,fmt,ERR=409,IOSTAT=ios),
     +        version,
     +        date_exp,recording,protocol_20,BDT_FILE_30,
     +        QDT_FILENAME_20,
     +        IDs,excluded,included,total_num_cells,
     +        I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +        BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +        BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +        end_time,icycles,total_histograms,
     +        (ETA2_1(i),i=1,120),(ETA2_2(i),i=1,120),
     +        (ETA2_3(i),i=1,120),(ETA2_4(i),i=1,120),
     +        (ETA2_5(i),i=1,120),(ETA2_6(i),i=1,120),
     +        (zmodsig_1(i),i=1,120),
     +        (zmodsig_2(i),i=1,120),(zmodsig_3(i),i=1,120),
     +        (zmodsig_4(i),i=1,120),(zmodsig_5(i),i=1,120),
     +        (zmodsig_6(i),i=1,120),
     +        (zmodsig2_1(i),i=1,120),
     +        (zmodsig2_2(i),i=1,120),(zmodsig2_3(i),i=1,120),
     +        (zmodsig2_4(i),i=1,120),
     +        (zmodsig2_5(i),i=1,120),
     +        (zmodsig2_6(i),i=1,120),(coef(i),i=1,120),
     +        (coefnum(i),i=1,120),
     +        (card_type(i),i=1,120),exp_name_20,
     +        (DELTA2(i),i=1,120),(tedfactor(i),i=1,120),
     +        (meanISI(i),i=1,120),(sdISI(i),i=1,120),
     +        (fiveHT(i),i=1,120),(mean_rISI(i),i=1,120),
     +        (sd_rISI(i),i=1,120),
     +        (num_rej_ISI(i),i=1,120),(num_rej_rISI(i),i=1,120),
     +        c_MAX_INT,TOTAL_NUM_SHIFTS,num_acc_cycles,
     +        ((sp_per_cycle(i,j),i=1,120),j=1,MAX_NUM_ACC_CYCLES)

         exp_name = exp_name_20
         protocol = protocol_20
         BDT_FILE = BDT_FILE_30
         QDT_FILENAME = QDT_FILENAME_20

      else
         goto 39
 1200    print *,'0cannot read header for QDTSAV ',QDTSAV
         stop
 1201    print *,'1cannot read header for QDTSAV ',QDTSAV
         stop
c           goto 50

 39      rewind 14
c           print '(''reading newest QDTSAV'')'

         read (14,'(A2,A11,A2,A200,A200,A50)',ERR=1200,IOSTAT=ios),
     +        version,
     +        date_exp,recording,protocol,BDT_FILE,QDT_FILENAME

c           print '(''version = '',A)',version
c           print '(''date_exp = '',A)',date_exp
c           print '(''recording = '',A)',recording
c           print '(''protocol = '',A)',protocol
c           print '(''BDT_FILE = '',A)',BDT_FILE
c           print '(''QDT_FILENAME = '',A)',QDT_FILENAME

         rewind 14

         if(version(1:1).eq.'6')then
            write(fmt,'(A25,I12,A5,I12,A5,I12,A59,I12,A6)')
     +           '(A2,A11,A2,A200,A200,A50,',
     +           MAX_NUM_CODES,
     +           '(I3),',
     +           MAX_NUM_CODES,
     +           '(I3),',
     +           MAX_NUM_CODES,
     +           '(I3),6I3,2f7.1,2f9.1,5f7.1,A1,'//
     +           '2f5.1,I3,I10,A175,A10,I3,I10,',
     +           MAX_NUM_CHAN,
     +           '(I10))'
c              print '(''version 6 read'')'
c     format string length = 148
         else if(version(1:1).ge.'7')then
c            print *,'version 7 read A'
            write(fmt,'(A25,I12,A5,I12,A5,I12,A59,I12,A6)')
     +           '(A2,A11,A2,A200,A200,A50,',
     +           MAX_NUM_CODES,
     +           '(I3),',
     +           MAX_NUM_CODES,
     +           '(I3),',
     +           MAX_NUM_CODES,
     +           '(I3),6I3,2f7.1,2f9.1,5f7.1,A1,'//
     +           '2f5.1,I6,I10,A175,A10,I3,I10,',
     +           MAX_NUM_CHAN,
     +           '(I10))'
c             print '(''version 7 read'')'
c     format string length = 151
         end if
         

         read (14,fmt,ERR=1201,IOSTAT=ios),
     +        version,
     +        date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +        IDs,excluded,included,total_num_cells,
     +        I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +        BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +        BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +        end_time,icycles,
     +        total_histograms,exp_name,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +        num_acc_cycles,ITAL

         if (icycles.eq.0.and.version.eq.'6') then
            print *,'icycles is 0 in ', QDTSAV
            print *,'(historically due to a version 7'
            print *, '.qdt.sav labeled as a version 6),  aborting'
            stop
         end if
            
 40      read (14,'(I3)',ERR=50,END=999,IOSTAT=ios),IDcode
c           print '(''reading info for IDcode '',I3)',IDcode
         i = IDcode
c           print '(''IDs('',I3,'') = '',I3)',IDcode,IDs(i)

         write(fmt,'(A32,I12,A5)')
     +        '(6A5,12A3,A5,f7.2,A3,A5,A2,7A10,',
     +        num_acc_cycles,
     +        '(I5))'
c     format string length = 49

         read (14,fmt,ERR=51,IOSTAT=ios),

     +        ETA2_1(IDs(i)),ETA2_2(IDs(i)),ETA2_3(IDs(i)),
     +        ETA2_4(IDs(i)),ETA2_5(IDs(i)),ETA2_6(IDs(i)),
     +        zmodsig_1(IDs(i)),zmodsig_2(IDs(i)),zmodsig_3(IDs(i)),
     +        zmodsig_4(IDs(i)),zmodsig_5(IDs(i)),zmodsig_6(IDs(i)),
     +        zmodsig2_1(IDs(i)),zmodsig2_2(IDs(i)),zmodsig2_3(IDs(i)),
     +        zmodsig2_4(IDs(i)),zmodsig2_5(IDs(i)),zmodsig2_6(IDs(i)),
     +        coef(IDs(i)),coefnum(IDs(i)),card_type(IDs(i)),
     +        DELTA2(IDs(i)),tedfactor(IDs(i)),meanISI(IDs(i)),
     +        sdISI(IDs(i)),fiveHT(IDs(i)),mean_rISI(IDs(i)),
     +        sd_rISI(IDs(i)),
     +        num_rej_ISI(IDs(i)),num_rej_rISI(IDs(i)),
     +        (sp_per_cycle(IDs(i),j),j=1,num_acc_cycles)

         goto 40

 51      print '(''problem reading data for IDcode '',I3)',IDcode
         print *, 'from ', QDTSAV, ', (historically due to a version 7'
         print *, '.qdt.sav labeled as a version 6),  aborting'
         stop

 50      close (unit = 14)

      end if

 999  close (unit=14)

      return
      end


*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*       >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


      subroutine write_QDTSAV(pgm_version,QDTSAV,ios,
     +     date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded,included,total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +     end_time,icycles,
     +     total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +     ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +     zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +     zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +     zmodsig2_6,coef,coefnum,card_type,exp_name,
     +     DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +     num_rej_ISI,num_rej_rISI,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +     num_acc_cycles,sp_per_cycle,ITAL)

      include 'x2000parameter.defs'
      
      integer IDs(MAX_NUM_CODES),excluded(MAX_NUM_CODES),
     +     included(MAX_NUM_CODES),total_num_cells,I_pulse,
     +     E_pulse,phrenic,BNDRY,cardiac_pls,icycles,
     +     total_histograms,TOTAL_NUM_SHIFTS,num_acc_cycles,
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     SPC_array,ITAL(MAX_NUM_CHAN)

      real BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,start_time,end_time,
     +     coefnum(MAX_NUM_CHAN)

      character*(*) exp_name,protocol,BDT_FILE,QDT_FILENAME,QDTSAV
      character*11 date_exp
      character*10 meanISI(MAX_NUM_CHAN),sdISI(MAX_NUM_CHAN),
     +     fiveHT(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT
      character*5 ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN),coef(MAX_NUM_CHAN),
     +     DELTA2(MAX_NUM_CHAN)
      character*3 zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),zmodsig2_6(MAX_NUM_CHAN),
     +     card_type(MAX_NUM_CHAN)
      character*2 pgm_version,recording,tedfactor(MAX_NUM_CHAN)
      character*1 BOUNDARY
      character*333 fmt


      SPC_array = MAX_NUM_CHAN*MAX_NUM_ACC_CYCLES

      goto 10
 100  print '(''error opening QDTSAV '',A,/,''error='',I5)',
     +     QDTSAV,ios
      goto 999
 410  print '(''410err'')'
 411  print '(''411err'')'
 412  print '(''412err writing QDTSAV '',A,/,''error='',I5)',
     +     QDTSAV,ios
      stop
c        goto 999


 10   OPEN (14,FILE=QDTSAV,FORM='FORMATTED',
     +     ERR=100,IOSTAT=ios)
c        print '(''writing newest QDTSAV'')'
      rewind 14

      write(fmt,'(A25,I12,A5,I12,A5,I12,A59,I12,A6)')
     +     '(A2,A11,A2,A200,A200,A50,',
     +     MAX_NUM_CODES,
     +     '(I3),',
     +     MAX_NUM_CODES,
     +     '(I3),',
     +     MAX_NUM_CODES,
     +     '(I3),6I3,2f7.1,2f9.1,5f7.1,A1,'//
     +     '2f5.1,I6,I10,A175,A10,I3,I10,',
     +     MAX_NUM_CHAN,
     +     '(I10))'
c     format string length = 151

      write (14,fmt,ERR=410,IOSTAT=ios)
     +     pgm_version,
     +     date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded,included,total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +     end_time,icycles,
     +     total_histograms,exp_name,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +     num_acc_cycles,ITAL

      do IDcode = 1, MAX_NUM_CODES
         if((IDs(IDcode).ne.0).and.(excluded(IDcode).ne.1))then
            i = IDcode
c              print '(''IDcode = '',I)',IDcode
            write (14,'(I3)',ERR=411,IOSTAT=ios)IDcode

            write(fmt,'(A32,I12,A5)')
     +           '(6A5,12A3,A5,f7.2,A3,A5,A2,7A10,',
     +           num_acc_cycles,
     +           '(I5))'
c     format string length = 53

            write (14,fmt,ERR=412,IOSTAT=ios)
     +           ETA2_1(IDs(i)),ETA2_2(IDs(i)),ETA2_3(IDs(i)),
     +           ETA2_4(IDs(i)),ETA2_5(IDs(i)),ETA2_6(IDs(i)),
     +           zmodsig_1(IDs(i)),zmodsig_2(IDs(i)),zmodsig_3(IDs(i)),
     +           zmodsig_4(IDs(i)),zmodsig_5(IDs(i)),zmodsig_6(IDs(i)),
     +           zmodsig2_1(IDs(i)),zmodsig2_2(IDs(i)),
     +           zmodsig2_3(IDs(i)),
     +           zmodsig2_4(IDs(i)),zmodsig2_5(IDs(i)),
     +           zmodsig2_6(IDs(i)),
     +           coef(IDs(i)),coefnum(IDs(i)),card_type(IDs(i)),
     +           DELTA2(IDs(i)),tedfactor(IDs(i)),meanISI(IDs(i)),
     +           sdISI(IDs(i)),fiveHT(IDs(i)),mean_rISI(IDs(i)),
     +           sd_rISI(IDs(i)),
     +           num_rej_ISI(IDs(i)),num_rej_rISI(IDs(i)),
     +           (sp_per_cycle(IDs(i),j),j=1,num_acc_cycles)
         end if
      end do

      close (unit=14)
*     ******************

 999  return
      end
      end module mod_read_and_write_DBSAV
