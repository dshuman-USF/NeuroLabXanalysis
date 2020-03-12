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

      module mod_database_write_new
      contains
*       filename = database_write_new.f
*
*       date of last revision = 14-Sep-2006             lss
*
*       sep-2003        lss
*        increased number of perturbation choices
*
*       sep-2001        lss
*        modified to automatically update PAIR database if data
*          for a CELL has been changed
*
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*       *** INDIRECT POINTERS are now used to access unit data
*       ***  within the following arrays:
*       ***     DATA_ARRAY, ITAL, card_type, coef, coefnum, TALLY_NORM,
*       ***     ETA2_*, zmodsig_*, zmodsig2_*, analyzed_cells, analyzed_pairs,
*       ***     resp_type, card, CELL_NAMES
*       ***     [these arrays are now dimensioned to MAX_NUM_CHAN;
*       ***      the location of a cell's information (user-assigned ID code = i)
*       ***      within these arrays is array(IDs(i))]
*       ***
*       *** DIRECT POINTERS are used with these arrays: excluded, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       This subroutine of x2000 produces text files containing database information
*               that will eventually be imported into the CELLS and PAIRS databases maintained
*               with Microsoft Access.  Specifications for these text files is as follows:
*
*
*               *****************************************************************
*               *                                                               *
*               * FORMAT OF OUTPUT ASCII FILE THAT WILL BE LOADED INTO THE      *
*               *   ACCESS DATABASE 'CELLS':                                    *
*               *                                                               *
*               * -- each field in the output file corresponds to a data        *
*               *    item for a cell (eg., name, type, coordinates) or for      *
*               *    a cell pair (eg., primary feature of histogram)            *
*               *                                                               *
*               * -- each record in the output file contains all the data       *
*               *    for one cell (each record contains ____ fields) or for     *
*               *    one cell pair (each record contains ____ fields)           *
*               *                                                               *
*               * -- records must be separated by <CR><LF>                      *
*               *                                                               *
*               * -- the number of records in an output file equals the         *
*               *    number of simultaneously recorded cells in that group      *
*               *    plus the number of pairs derived from those cells          *
*               *                                                               *
*               * -- text strings containing blanks must be enclosed in         *
*               *    quotation (") marks                                        *
*               *                                                               *
*               * -- all fields may be enclosed in quotation marks              *
*               *                                                               *
*               * -- fields must be separated from one another by either a      *
*               *    comma or a blank space                                     *
*               *                                                               *
*               * -- the output filename must be lower case letters             *
*               *                                                               *
*               *****************************************************************
*
*
      subroutine database_write(DB_FILES)
*
      use mod_miscellaneous_subroutines
      use mod_read_and_write_DBSAV
      use mod_read_and_write_pre_files
      use mod_read_and_write_ref_electrode_data
      include 'x2000parameter.defs'
*
      integer ID(2),rec_num,hist_num,included(MAX_NUM_CODES),
     +     excluded(MAX_NUM_CODES),analyzed_cells(MAX_NUM_CHAN),
     +     analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +     perturb(50),
     +     total_num_cells,IDs(MAX_NUM_CODES),
     +     total_num_cells_db,total_histograms,
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     REF,TAR,ITAL(MAX_NUM_CHAN),total_num_electrode_arrays,
     +     total_num_qdts,cardiac_pls,BNDRY,phrenic,I_pulse,
     +     E_pulse

      integer*4 isys

      logical QDTexist,EXIST_extra_preDBG
*
      real NORM_OFFSET,STIM_OFFSET,Q_pos,BINW_1,
     +     BINW_2,BINW_3,BINW_4,BINW,BINW2,start_time,end_time,
     +     coefnum(MAX_NUM_CHAN),eff,det,vis,NORM_BW
*
      character*220 DBG,DBG_AA,DBG_STA,DBR_xass,DBR_tsfs,DBP,
     +     DBSAV,pre_DBG,pre_DBP,EL_REF,DBP_eff,
     +     DBG_num_events,DBP_statsbw
      character*250 qdtfilenames,qdtfilenames_1,FILE_1,pre_DBG_EXTRA,
     +     exp_name_250
      character*200 blank200,dir_name
      character*60 QDTSAV
      character*50 QDT_FILENAME,QDTfiles(MAX_NUM_QDTS+5),text50
      CHARACTER*40 comm, rescom,statcomm,text2,blank40
*
      character*30 BDT_FILE_30
      character*(*) DB_FILES
      character*200 BDT_FILE,protocol,exp_name_200
      character*175 exp_name,exp_name_x
*
      CHARACTER*20 prim,sec,text3
*
      CHARACTER*15 loc1,loc2,resp_type_1,resp_type(MAX_NUM_CHAN)
*
      character*12 STA_phrenic,STA_RLN,STA_lumbar,
     +     STA_ELN,STA_cent_vagus,STA_splanch,
     +     STA_cerv_symp,
     +     STA_new1,STA_new2,STA_new3,STA_results(MAX_STA)
*
      character*11 date
      character*10 fiveHT(MAX_NUM_CHAN),meanISI(MAX_NUM_CHAN),
     +     sdISI(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT,blank10

*
      CHARACTER*8 czk,cprobk,cdet,cvis,chalfwd,czlat,stats_bw

*
      character*9 carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,
     +     dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,carotidCO2_2_x,
     +     pertnew1_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +     sw_x,SLNsw_x,per_results(MAX_PERTURB)

      CHARACTER*5 aAP, aRL,adep,coef(MAX_NUM_CHAN),
     +     ETA2_1(MAX_NUM_CHAN),ETA2_2(MAX_NUM_CHAN),
     +     ETA2_3(MAX_NUM_CHAN),ETA2_4(MAX_NUM_CHAN),
     +     ETA2_5(MAX_NUM_CHAN),ETA2_6(MAX_NUM_CHAN),
     +     DELTA2(MAX_NUM_CHAN),AP_refs(MAX_NUM_ARRAYS),
     +     RL_refs(MAX_NUM_ARRAYS),
     +     depth_refs(MAX_NUM_ARRAYS),blank5
*
      CHARACTER*4 cell_name(2),CELL_NAMES(MAX_NUM_CHAN),
     +     NAME_refs(MAX_NUM_ARRAYS),dchan(MAX_NUM_CHAN),
     +     ref_chan(MAX_NUM_CHAN),ref_chan_refs(MAX_NUM_ARRAYS),
     +     blank4
*
      CHARACTER*3 AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +     AA_rVRG,AA_rtPRG,AA_ltPRG,
     +     AA_new1,AA_new2,AA_new3,AA_results(MAX_AA),
     +     card_type(MAX_NUM_CHAN),serotonin,blank3
      character*3 zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),
     +     zmodsig_4(MAX_NUM_CHAN),zmodsig_5(MAX_NUM_CHAN),
     +     zmodsig_6(MAX_NUM_CHAN),zmodsig2_1(MAX_NUM_CHAN),
     +     zmodsig2_2(MAX_NUM_CHAN),zmodsig2_3(MAX_NUM_CHAN),
     +     zmodsig2_4(MAX_NUM_CHAN),zmodsig2_5(MAX_NUM_CHAN),
     +     zmodsig2_6(MAX_NUM_CHAN)
*               
      character*2 carotidCO2(5),vertCO2(5),hypercap_5_O2(5),
     +     hypercap_5_air(5),hypercap_tbd(5),hypoxia_5(5),
     +     hypoxia_12(5),gasp(5),
     +     lobel(5),aorta_cath(5),pe(5),carotid_occ(5),
     +     nitrop(5),vc_cath(5),dopamine(5),vagus(5),
     +     capsaicin(5),pinch(5),deep_pain(5),
     +     codeine(5),nalox(5),methyserg(5),mucimol(5),
     +     dexameth(5),
     +     noinf(5),hyperinf(5),hypervent(5),carotidCO2_2(5),
     +     pertnew1(5),cgh(5),lcgh(5),SLNcgh(5),
     +     exp_reflex(5),
     +     sw(5),SLNsw(5),blank2,
     +     recording,recording_x,tedfactor(MAX_NUM_CHAN),version,
     +     subset,which_subset

      character*1 perturb_applied,
     +     AA_applied,STA_applied,AA(10),STA(10),
     +     info_prev_imported,per_prev_imported,
     +     sort,card,BOUNDARY,write_the_preDBP_file,blank1,dup,
     +     write_anyway

      logical exist
      INTEGER*4 SYSTEM
      EXTERNAL SYSTEM

      print '(''NEW DB_WRITE'')'

      blank200 = ' '
      blank40 = ' '
      blank10 = ' '
      blank5 = ' '
      blank4 = ' '
      blank3 = ' '
      blank2 = ' '
      blank1 = ' '
      subset = ' '
      which_subset = ' '

 11   print '(10(/),''Are these files a subset of existing data?''
     +''  (e.g., gasping, pre-vagotomy, or no-inflation''
     +'' data) (y/n)  >> '',$)'
      read (*,'(A)') subset
      call upper_case(subset,LEN(subset))
      if(subset.eq.'Y')then
         print '(/,T10,''g --> gasping'',
     +/,T10,''v --> pre-vagotomy'',
     +/,T10,''i --> no-inflation'',
c     +          /,T10,'' --> '',
     +//,T5,''Which condition applies? (<cr> to esc)  >> ''
     +,$)'
         read (*,'(A)') which_subset
         call upper_case(which_subset,LEN(which_subset))
         if(which_subset.eq.'G')then
         elseif(which_subset.eq.'V')then
         elseif(which_subset.eq.'I')then
         elseif(which_subset.eq.' ')then
            subset = 'N'
            which_subset = ' '
         else                   !force an appropriate response
            which_subset = ' '  !re-initialize the variable and ask the question again
            goto 11
         end if
      elseif(subset.eq.'N')then
      else
         goto 11
      end if

      call concat(DB_FILES,'_ARRAY_REFS.sav',EL_REF,l)
      inquire (FILE=EL_REF(1:l),EXIST=exist)
      if(exist.eqv..TRUE.)then
         call concat(DB_FILES,'_ARRAY_REFS.txt',EL_REF,l)
         OPEN(4,FILE=EL_REF(1:l)//char(0),FORM='FORMATTED',
     +        ACCESS='SEQUENTIAL')
      else
         print '(//,''Not ready to write DB-ready files. '',//,T5,
     +''... could not find data for reference electrodes ''
     +''(*_ARRAY_REFS.sav file)'',/,T5,
     +''... no files written'')'
         print '(/,T5,''<cr> to continue  >> '',$)'
         read (*,'(A)')
         return
      end if
      call concat(DB_FILES,'.db.sav',DBSAV,l)
      OPEN (15,FILE=DBSAV(1:l)//char(0),FORM='FORMATTED',
     +     STATUS='OLD',ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_UNITS.txt',DBG,l) !database info for indiv cells (during eupneic-like cond)
      OPEN (16,FILE=DBG(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_UNITS_AA.txt',DBG_AA,l) !database info for indiv cells (antidromic activation data)
      OPEN (161,FILE=DBG_AA(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_UNITS_STA.txt',DBG_STA,l) !database info for indiv cells (spike-triggered avg data)
      OPEN (162,FILE=DBG_STA(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_num_events.txt',DBG_num_events,l) !database info for number of events in each spike train
      OPEN (164,FILE=DBG_num_events(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'.pre_dbg',pre_DBG,l) !precursor file to final database file for cells
      OPEN (9,FILE=pre_DBG(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='DIRECT',RECL=DBG_RECL) !+ (6*char*9) + (29*char*2) = 112 for more perturbation responses'
      call concat(DB_FILES,'_RESPONSES_xass.txt',DBR_xass,l) !database info for perturbation/stimulus responses of indiv cells
      OPEN (20,FILE=DBR_xass(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_RESPONSES_tsfs.txt',DBR_tsfs,l) !database info for perturbation/stimulus responses of indiv cells
      OPEN (21,FILE=DBR_tsfs(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      
      call concat(DB_FILES,'_CROSSES.txt',DBP,l) !database info for pairs 
      OPEN (17,FILE=DBP(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_effectiveness.txt',DBP_eff,l) !database info for effectiveness of correlated pairs with peaks
      OPEN (171,FILE=DBP_eff(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'_stats_bw.txt',DBP_statsbw,l) !database info for binwidth used to calculate the feature stats
      OPEN (172,FILE=DBP_statsbw(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
      call concat(DB_FILES,'.pre_dbp',pre_DBP,l) !precursor file to final database file for pairs
      OPEN (10,FILE=pre_DBP(1:l)//char(0),FORM='FORMATTED',
     +     ACCESS='DIRECT',RECL=iDBP_RECL)


*       *****   read the database "gamesave" file:      *****

      call read_DBSAV(version,DBSAV,ios,included,excluded,
     +     analyzed_cells,analyzed_pairs,CELL_NAMES,
     +     perturb_applied,
     +     perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +     info_prev_imported,
     +     per_prev_imported,total_num_cells_db,qdtfilenames)

*       ***** parse qdtfilenames to see which QDT files have been associated with these database files: *****

      j = 0
      k = 1
      do i = 1,LEN(qdtfilenames)

         if(qdtfilenames(i:i).eq.';')then !filenames separated with a semicolon
            dup = 'n'
            do m = 1, MAX_NUM_QDTS+5 !be sure that QDT filenames are not duplicated
               if(qdtfilenames(k:i-1).eq.QDTfiles(m))dup = 'y'
            end do
            if(dup.eq.'n')then
               j = j+1
               QDTfiles(j) = qdtfilenames(k:i-1)
               k = i+1
            end if
         end if
      end do
      call strlength(DB_FILES,LEN(DB_FILES),l)
      print '(//,''The following QDT files have been associated ''
     +''with '',
     +A,'':''/)',DB_FILES(1:l)
      qdtfilenames_1 = ' '
      ii = 0
      total_num_qdts = 0
      do i = 1, j
         inquire(FILE=QDTfiles(i),EXIST=QDTexist)
         if(QDTexist.eqv..TRUE.)then
            call strlength(QDTfiles(i),LEN(QDTfiles(i)),l)
            qdtfilenames_1(ii+1:ii+l) = QDTfiles(i)(1:l)
            qdtfilenames_1(ii+l+1:ii+l+1) = ';'
            ii = ii+l+1
            print '(3x,A,$)',QDTfiles(i)(1:l)
            total_num_qdts = total_num_qdts + 1
         end if
      end do
c        print '(/,''qdtfilenames = '',A)',qdtfilenames
c        print '(''qdtfilenames_1 = '',A)',qdtfilenames_1
      if(qdtfilenames.ne.qdtfilenames_1)then !update the QDTfile list
         qdtfilenames = qdtfilenames_1
c           print '(''update DBSAV'')'
         call write_DBSAV(version,DBSAV,ios,included,excluded,
     +        analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)
      end if
 5    print '(//,''Choose one which contains the correct ''
     +''data for statistical evaluations of ''
     +/,''respiratory and cardiac modulation.''
     +'' (i.e., ETA2, BINARY, and DELTA2 tests)''
     +//,T15,''>> '',$)'
      read (*,'(A)') QDT_FILENAME
      if(QDT_FILENAME.eq.' ')return !esc from writing the files
      
      
      inquire (FILE=QDT_FILENAME,EXIST=QDTexist)
      if(QDTexist.eqv..TRUE.)then
         call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
         QDTSAV = QDT_FILENAME(1:l)//'.sav'
      else
         print '(//,''*** OOPS -- invalid filename.  ''
     +''Try again, please.'')'
         goto 5
      end if
      print '(//,''please wait ... writing text files for later ''
     +''import into database application ...'')'

*       ***** read the QDT.SAV file: ****

      call read_QDTSAV(version,QDTSAV,ios,
     +     date,recording,protocol,BDT_FILE,QDT_FILENAME,
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
     +     num_rej_ISI,num_rej_rISI,c_MAX_INT,jTOTAL_SHIFTS,
     +     num_acc_cycles,sp_per_cycle,ITAL)


*       ***** Check to be sure that meanISI, sdISI, and Peggy Mason statistic for serotonergic properties is valid: 
*       ***** Also, check to be sure that mean_rISI, sd_rISI are good. 
*       ***** Theoretically, this block of 'CYA' code can eventually be dropped. 
*       *****   Code now exists in serotergic_discriminant_function.f to handle these situations up front on the 
*       *****   calculation end. 

      do k = 1, MAX_NUM_CODES
         if (IDs(k).eq.0) cycle
         if(meanISI(IDs(k)).eq.'*******'.or.
     +        meanISI(IDs(k)).eq.'No data')then
            meanISI(IDs(k))='000000.00'
            sdISI(IDs(k))='0.0'
         endif
         print '(I4,''; '',I4,''; meanISI = '',A10)',
     +        k,IDs(k),meanISI(IDs(k))
         if((IDs(k).ne.0).and.(excluded(k).ne.1))then
            if((fiveHT(IDs(k)).eq.'NaN').or.
     +           (fiveHT(IDs(k)).eq.'+INF').or.
     +           (fiveHT(IDs(k)).eq.'-INF'))then
               meanISI(IDs(k))='No data'
               sdISI(IDs(k))='No data'
               fIVEHT(IDS(K))='NO DATA'
            END IF
            IF((MEAN_RISI(IDS(K)).EQ.'NAN').OR.
     +           (MEAN_RISI(IDS(K)).EQ.'+INF').OR.
     +           (MEAN_RISI(IDS(K)).EQ.'-INF'))THEN
               MEAN_RISI(IDS(K))='NO DATA'
               SD_RISI(IDS(K))='NO DATA'
            END IF
            IF(FIVEHT(IDS(K)).NE.'NO DATA')THEN
C               PRINT '(''ID: '',I3,''; 5-HT = '',A,''; MEANISI = '',   
C     +            A10,''; SDISI = '',A10)',
C     +            K,FIVEHT(IDS(K)),MEANISI(IDS(K)),SDISI(IDS(K))
               READ(MEANISI(IDS(K)),'(F7.2)') XMEAN
               read(sdISI(IDs(k)),'(f7.2)') xsd
               xfive = 146 - xmean + 0.98*xsd !re-calculate and store in scientific notation
c                print '(''xfive = '',f,''; xmean = '',f,
c     +                 ''; xsd = '',f)',xfive,xmean,xsd
               write (fiveHT(IDs(k)),'(ES10.2E2)')xfive
c                print '(''fiveHT = '',A)',fiveHT(IDs(k))
c                read (*,'(A)')
            end if
         end if
      end do
      if (version .eq. '6') version = '7'
      if (version(1:1) .lt. '6') then
         print '(10(/),''You are trying to use QDTSAV ''
     +''files that were created with v5 or ''
     +''earlier of XAnalysis.'',/,
     +''This version of xanalysis will not do that.'',/,
     +''Database files not written.'')'
         return
      end if
      call write_QDTSAV(version,QDTSAV,ios,
     +     date,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded,included,total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,
     +     start_time,end_time,icycles,
     +     total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +     ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +     zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +     zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +     zmodsig2_6,coef,coefnum,card_type,exp_name,
     +     DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +     num_rej_ISI,num_rej_rISI,c_MAX_INT,jTOTAL_SHIFTS,
     +     num_acc_cycles,sp_per_cycle,ITAL)
c             end if
c           end if
c        end do
c        read (*,'(A)')

      

      do i = 1, MAX_NUM_CHAN    !force all these to UPPER CASE LETTERS
         call upper_case(zmodsig_1(i),LEN(zmodsig_1(i)))
         call upper_case(zmodsig_2(i),LEN(zmodsig_2(i)))
         call upper_case(zmodsig_3(i),LEN(zmodsig_3(i)))
         call upper_case(zmodsig_4(i),LEN(zmodsig_4(i)))
         call upper_case(zmodsig_5(i),LEN(zmodsig_5(i)))
         call upper_case(zmodsig_6(i),LEN(zmodsig_6(i)))
         call upper_case(zmodsig2_1(i),LEN(zmodsig2_1(i)))
         call upper_case(zmodsig2_2(i),LEN(zmodsig2_2(i)))
         call upper_case(zmodsig2_3(i),LEN(zmodsig2_3(i)))
         call upper_case(zmodsig2_4(i),LEN(zmodsig2_4(i)))
         call upper_case(zmodsig2_5(i),LEN(zmodsig2_5(i)))
         call upper_case(zmodsig2_6(i),LEN(zmodsig2_6(i)))
         call upper_case(card_type(i),LEN(card_type(i)))
      end do
      do i = 1, MAX_NUM_ARRAYS
         AP_refs(i) = ' '
         RL_refs(i) = ' '
         depth_refs(i) = ' '
         NAME_refs(i) = ' '
         ref_chan_refs(i) = ' '
      end do

*       *********************************************************************************
*       *********************   create the UNITS database file:  ************************
*       *********************************************************************************

*             *****  write the header information to allow import into the Access db using a macro: *****
*                       (WARNING - these must EXACTLY match the db fieldnames!) 

      write (16,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","Cell Name","ID code","Resp Pattern",''
     +''"AP coord","RL coord","Depth coord",''
     +''"Datamax channel",''
     +''"Channel of ref elec",''
     +''"Cleanness of Sort","Delta-squared value",''
     +''"Ted Factor (TF)","Cardiac Modulation Significant?",''
     +''"Cardiac Modulation - Visual","Cell Comments",''
     +''"Coef of Variance","ETA2 - 1","ETA2 - 2","ETA2 - 3",''
     +''"ETA2 - 4","ETA2 - 5","ETA2 - 6",''
     +''"Resp mod: anova test - 1","anova test - 2",''
     +''"anova test - 3","anova test - 4","anova test - 5",''
     +''"anova test - 6","Resp Mod: binary test - 1",''
     +''"binary test - 2","binary test - 3",''
     +''"binary test - 4",''
     +''"binary test - 5","binary test - 6",''
     +''""Serotonin-Like" Stat (PMS)","Mean ISI",''
     +''"SD of ISI",''
     +''"Serotonin-Like?","Mean rISI","SD of rISI",''
     +''"# of ints cut from analysis"'',$)')
      write (16,'(A1)')char(13)

      write (161,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","Cell Name","ID code",''
     +''"Results of AA- Spinal Cord","AA-RLN","AA-Vagus",''
     +''"AA-cVRG","AA-rVRG","AA-right PRG","AA-left PRG",''
     +''"AA-extra1","AA-extra2","AA-extra3"'',$)')
      write (161,'(A1)')char(13)

      write (162,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","Cell Name","ID code",''
     +''"Results of STA-Phrenic","STA-RLN",''
     +''"STA-central Vagus","STA-Lumbar",''
     +''"STA-Cervical sympathetic","STA-ELN",''
     +''"STA-Splanchnic",''
     +''"STA-extra1","STA-extra2","STA-extra3"'',$)')
      write (162,'(A1)')char(13)

      write (164,'(''"Exp Name","Recording Number","Exp Date",''
     +''"ID code","Spike total"'',$)')
      write (164,'(A1)')char(13)

      write (20,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","Cell Name","ID code",''
     +''"Carotid CO2-general","Vertebral CO2-general",''
     +''"Hypercapnia (5%CO2 in O2)-general",''
     +''"Hypercapnia (5%CO2 in air)-general",''
     +''"Hypercapnia (tbd)-general","Hypoxia (12%)-general",''
     +''"Hypoxia (5%)-general","repeated Hypoxia (12%)-general",''
     +''"Carotid CO2 during hyperventilation-general",''
     +''"repeated Carotid CO2-general","Lobelline-general",''
     +''"EMB Baro up (Aortic balloon)-general",''
     +''"phenylephrine-general",''
     +''"Carotid Occlusion-general","Nitroprusside-general",''
     +''"EMB Baro down (vena cava balloon)-general",''
     +''"Dopamine-general","Vagotomy-general",''
     +''"Capsaicin-general","Pinch-general",''
     +''"Deep pain-general","Codeine-general",''
     +''"Naloxone-general","Methysergide-general",''
     +''"Muscimol-general","Dexamethasone-general",''
     +''"No lung inflation-general",''
     +''"Hyperinflation-general","Hyperventilation-general",''
     +''"Tracheal cough-general",''
     +''"Laryngeal Cough - general",''
     +''"ER - general","SLN Cough - general",''
     +''"Mechanical swallow - general",''
     +''"SLN swallow - general"'',$)')
      write (20,'(A1)')char(13)

      write (21,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","Cell Name","ID code","Carotid CO2-1",''
     +''"Carotid CO2-2","Carotid CO2-3","Carotid CO2-4",''
     +''"Carotid CO2-5","Vertebral CO2-1","Vertebral CO2-2",''
     +''"Vertebral CO2-3","Vertebral CO2-4",''
     +''"Vertebral CO2-5","Hypercapnia (5%CO2 in O2)-1",''
     +''"Hypercapnia (5%CO2 in O2)-2",''
     +''"Hypercapnia (5%CO2 in O2)-3",''
     +''"Hypercapnia (5%CO2 in O2)-4",''
     +''"Hypercapnia (5%CO2 in O2)-5",''
     +''"Hypercapnia (5%CO2 in air)-1",''
     +''"Hypercapnia (5%CO2 in air)-2",''
     +''"Hypercapnia (5%CO2 in air)-3",''
     +''"Hypercapnia (5%CO2 in air)-4",''
     +''"Hypercapnia (5%CO2 in air)-5",''
     +''"Hypercapnia (tbd)-1","Hypercapnia (tbd)-2",''
     +''"Hypercapnia (tbd)-3","Hypercapnia (tbd)-4",''
     +''"Hypercapnia (tbd)-5","Hypoxia (12%)-1",''
     +''"Hypoxia (12%)-2","Hypoxia (12%)-3",''
     +''"Hypoxia (12%)-4","Hypoxia (12%)-5",''
     +''"Hypoxia (5%)-1","Hypoxia (5%)-2","Hypoxia (5%)-3",''
     +''"Hypoxia (5%)-4","Hypoxia (5%)-5","Gasping-1",''
     +''"Gasping-2","Gasping-3","Gasping-4","Gasping-5",''
     +''"blank perturbation 1- 1","blank perturbation 1- 2",''
     +''"blank perturbation 1- 3","blank perturbation 1- 4",''
     +''"blank perturbation 1- 5",''
     +''"repeated Carotid CO2-1",'',$)')
      write (21,'(''"repeated Carotid CO2-2",''
     +''"repeated Carotid CO2-3",''
     +''"repeated Carotid CO2-4","repeated Carotid CO2-5",''
     +''"Lobelline-1","Lobelline-2","Lobelline-3",''
     +''"Lobelline-4","Lobelline-5",''
     +''"EMB Baro up (Aortic balloon)-1",''
     +''"EMB Baro up (Aortic balloon)-2",''
     +''"EMB Baro up (Aortic balloon)-3",''
     +''"EMB Baro up (Aortic balloon)-4",''
     +''"EMB Baro up (Aortic balloon)-5",''
     +''"Phenylephrine-1","Phenylephrine-2",''
     +''"Phenylephrine-3","Phenylephrine-4",''
     +''"Phenylephrine-5","Carotid Occlusion-1",''
     +''"Carotid Occlusion-2","Carotid Occlusion-3",''
     +''"Carotid Occlusion-4","Carotid Occlusion-5",''
     +''"Nitroprusside-1","Nitroprusside-2",''
     +''"Nitroprusside-3","Nitroprusside-4",''
     +''"Nitroprusside-5",''
     +''"EMB Baro down (vena cava balloon)-1",''
     +''"EMB Baro down (vena cava balloon)-2",''
     +''"EMB Baro down (vena cava balloon)-3",''
     +''"EMB Baro down (vena cava balloon)-4",''
     +''"EMB Baro down (vena cava balloon)-5",'',$)')
      write (21,'(''"Dopamine-1","Dopamine-2",''
     +''"Dopamine-3","Dopamine-4",''
     +''"Dopamine-5","Vagotomy-1","Vagotomy-2","Vagotomy-3",''
     +''"Vagotomy-4","Vagotomy-5","Capsaicin-1",''
     +''"Capsaicin-2","Capsaicin-3","Capsaicin-4",''
     +''"Capsaicin-5","Pinch-1","Pinch-2","Pinch-3",''
     +''"Pinch-4","Pinch-5","Deep pain-1","Deep pain-2",''
     +''"Deep pain-3","Deep pain-4","Deep pain-5",''
     +''"Codeine-1","Codeine-2","Codeine-3","Codeine-4",''
     +''"Codeine-5","Naloxone-1","Naloxone-2","Naloxone-3",''
     +''"Naloxone-4","Naloxone-5","Methysergide-1",''
     +''"Methysergide-2","Methysergide-3","Methysergide-4",''
     +''"Methysergide-5","Muscimol-1","Muscimol-2",''
     +''"Muscimol-3","Muscimol-4","Muscimol-5",''
     +''"Dexamethasone-1","Dexamethasone-2",''
     +''"Dexamethasone-3","Dexamethasone-4",''
     +''"Dexamethasone-5","No lung inflation-1",''
     +''"No lung inflation-2","No lung inflation-3",''
     +''"No lung inflation-4","No lung inflation-5",''
     +''"Hyperinflation-1","Hyperinflation-2",''
     +''"Hyperinflation-3","Hyperinflation-4",''
     +''"Hyperinflation-5","Hyperventilation-1",''
     +''"Hyperventilation-2","Hyperventilation-3",''
     +''"Hyperventilation-4","Hyperventilation-5",''
     +''"Tracheal cough-1","Tracheal cough-2",''
     +''"Tracheal cough-3","Tracheal cough-4",''
     +''"Tracheal cough-5","LAR-1","LAR-2","LAR-3","LAR-4",''
     +''"LAR-5","ER-1","ER-2","ER-3","ER-4","ER-5","SLN-1",''
     +''"SLN-2","SLN-3","SLN-4","SLN-5","SWM-1","SWM-2",''
     +''"SWM-3","SWM-4","SWM-5","SWSLN-1","SWSLN-2",''
     +''"SWSLN-3","SWSLN-4","SWSLN-5"'',$)')
      write (21,'(A1)')char(13)

      write (17,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Spike Time Filename","Histogram Filenames","TBD1",''
     +''"TBD2","TBD3","R-Cell Name","R-ID Code",''
     +''"T-Cell Name","T-ID Code","Primary feature",''
     +''"Location of primary","Secondary feature",''
     +''"Location of secondary","comments - cch (pair)",''
     +''"k-statistic value","k-significance","Detect Index",''
     +''"Visibility Index","QDT File// binwidth",''
     +''"time lag to feature","1/2 width of feature"'',$)')
      write (17,'(A1)')char(13)

      write (171,'(''"Exp Name","Recording Number","Exp Date",''
     +''"R-ID Code","T-ID Code","Effectiveness"'',$)')
      write (171,'(A1)')char(13)

      write (172,'(''"Exp Name","Recording Number","Exp Date",''
     +''"R-ID Code","T-ID Code","stats_bw"'',$)')
      write (172,'(A1)')char(13)
      
      rewind 4                  !overwrite anything already in the file
      write (4,'(''"Exp Name","Recording Number","Exp Date",''
     +''"Ref unit","AP coord","RL coord","Depth coord",''
     +''"Channel number"'',$)')
      write (4,'(A1)')char(13)

*               *****  read info from the pre_DBG file and the UNITS_REFS file and use it to create  *****
*               *****   the DBG file:                        *****

      print '(//,''Writing UNITS and RESPONSES files ....'')'

      call read_ref_electrode_data(BDT_FILE,DB_FILES,
     +     date,recording,exp_name_x,
     +     NAME_refs,AP_refs,RL_refs,depth_refs,
     +     dchan,ref_chan,ref_chan_refs,
     +     total_num_electrode_arrays)

      write_anyway = '?'
      
      do i = 1,MAX_NUM_CODES
         if(IDs(i).eq.0)cycle

         if(write_anyway.eq.'?'.and.analyzed_cells(IDs(i)).ne.2
     +        .and.included(i).gt.0)then
            print *, "Data not analyzed for cell ", IDs(i)
            print '(A,$)', "Write this cell's data and all subsequent an
     +alysed and unanalyzed data for all cells and pairs anyway? (Y/N):
     +"
            read (*,fmt='(A1)') write_anyway
            call upper_case(write_anyway,LEN(write_anyway))
         endif
         
         if(analyzed_cells(IDs(i)).eq.2.or.(write_anyway.eq.'Y'
     +        .and.included(i).gt.0))then
            rec_num = included(i) !record number of this cell's info in pre_dbg file

*               *****  read a record from the .pre_dbg file:  *****

            call read_preDBG(9,rec_num,ios,date,BDT_FILE_30,
     +           recording_x,hist_num,
     +           cell_name(1),ID(1),resp_type_1,aAP,aRL,adep,
     +           sort,comm,
     +           card,
     +           AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +           AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +           STA_phrenic,STA_RLN,STA_cent_vagus,
     +           STA_lumbar,STA_cerv_symp,STA_ELN,
     +           STA_splanch,STA_new1,STA_new2,STA_new3,
     +           carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +           hypercap_5_air_x,
     +           hypercap_tbd_x,hypoxia_12_x,hypoxia_5_x,gasp_x,
     +           lobel_x,
     +           aorta_cath_x,pe_x,carotid_occ_x,nitrop_x,
     +           vc_cath_x,dopamine_x,
     +           vagus_x,capsaicin_x,pinch_x,deep_pain_x,
     +           codeine_x,nalox_x,
     +           methyserg_x,mucimol_x,dexameth_x,noinf_x,
     +           hyperinf_x,hypervent_x,
     +           pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +           SLNcgh_x,exp_reflex_x,sw_x,
     +           SLNsw_x,
     +           carotidCO2,vertCO2,hypercap_5_O2,
     +           hypercap_5_air,
     +           hypercap_tbd,hypoxia_12,hypoxia_5,
     +           gasp,lobel,
     +           aorta_cath,pe,carotid_occ,
     +           nitrop,vc_cath,dopamine,
     +           vagus,capsaicin,pinch,
     +           deep_pain,codeine,nalox,
     +           methyserg,mucimol,dexameth,
     +           noinf,hyperinf,hypervent,
     +           pertnew1,carotidCO2_2,cgh,lcgh,
     +           SLNcgh,exp_reflex,sw,
     +           SLNsw,per_results,AA_results,STA_results)

            if(cardiac_pls.eq.0)card=' '
            
            if(cell_name(1).eq.' ')cycle !do not write info if cell has not been named
c              if(IDs(i).ne.0)
c     +          print '(''5-HT = '',A)',fiveHT(IDs(i))
*       ***************************************************************************************************
*       THIS CODE SHOULD BE SUPERFLUOUS.  THIS SITUATION SHOULD HAVE BEEN HANDLED ABOVE WHEN THE
*       QDTSAV FILE WAS RE-WRITTEN.
            if(fiveHT(IDs(i)).eq.'**********')then !uh-oh; recalculate this value
               print '(''fiveHT('',I3,'') = ***'')',IDs(i)
               read(meanISI(IDs(i)),'(f7.2)') xmean
               read(sdISI(IDs(i)),'(f7.2)') xsd
               xfive = 146 - xmean + 0.98*xsd
               print '(''xfive = '',f15.7,''; xmean = '',f15.7,
     +''; xsd = '',f15.7)',xfive,xmean,xsd
               write (fiveHT(IDs(i)),'(ES10.2E2)')xfive
               print '(''fiveHT = '',A)',fiveHT
               read (*,'(A)')
            end if
*       ***************************************************************************************************
            do ii = 1,LEN(fiveHT(IDs(i)))
               if(fiveHT(IDs(i))(ii:ii).ne.' ')exit !find the first "real" character in the string
            end do
            if(fiveHT(IDs(i))(ii:ii).eq.'-')then
               serotonin='SL'   !this cell displays serotonergic-like properties
            else
               serotonin='NSL'
            end if
            text2 = ' '
            text3 = ' '
            text3 = '(>'//c_MAX_INT//'ms.)'
            call remove_all_blanks(text3,LEN(text3))
            text50 = num_rej_ISI(IDs(i))//' ISIs rejected '
     +           //trim(text3)
            call remove_leading_blanks(text50,LEN(text50))
            text2=text50(1:40)

*               *****  write the record to the UNITS file:  *****

            write (16,'(9(''"'',A,''",''),
     +''"'',I5,''",'',37(''"'',A,''",''),
     +''"'',A,''"'',$)')
     +           exp_name,recording,date,BDT_FILE,
     +           qdtfilenames,blank200,blank200,blank200,
     +           cell_name(1),ID(1),
     +           resp_type_1,aAP,aRL,adep,dchan(IDs(i)),
     +           ref_chan(IDs(i)),
     +           sort,DELTA2(IDs(i)),tedfactor(IDs(i)),
     +           card_type(IDs(i)),
     +           card,comm,coef(IDs(i)),ETA2_1(IDs(i)),
     +           ETA2_2(IDs(i)),ETA2_3(IDs(i)),ETA2_4(IDs(i)),
     +           ETA2_5(IDs(i)),ETA2_6(IDs(i)),
     +           zmodsig_1(IDs(i)),zmodsig_2(IDs(i)),zmodsig_3(IDs(i)),
     +           zmodsig_4(IDs(i)),zmodsig_5(IDs(i)),zmodsig_6(IDs(i)),
     +           zmodsig2_1(IDs(i)),zmodsig2_2(IDs(i)),
     +           zmodsig2_3(IDs(i)),
     +           zmodsig2_4(IDs(i)),zmodsig2_5(IDs(i)),
     +           zmodsig2_6(IDs(i)),
     +           fiveHT(IDs(i)),meanISI(IDs(i)),sdISI(IDs(i)),serotonin,
     +           mean_rISI(IDs(i)),
     +           sd_rISI(IDs(i)),text2
            write (16,'(A1)')char(13)

*               *****  write the record to the NUM_EVENTS file:  *****

            write (164,'(3(''"'',A,''",''),
     +''"'',I5,''","'',I10,''"'',$)')
     +           exp_name,recording,date,ID(1),ITAL(IDs(ID(1)))
            write (164,'(A1)')char(13)
            print '(''ID(1) = '',I5,''; ITAL = '',I10)',
     +           ID(1),ITAL(IDs(ID(1)))


*               *****  write the record to the AA file:  *****

            write (161,'(9(''"'',A,''",''),
     +''"'',I5,''",'',9(''"'',A,''",''),
     +''"'',A,''"'',$)')
     +           exp_name,recording,date,BDT_FILE,
     +           qdtfilenames,blank200,blank200,blank200,
     +           cell_name(1),ID(1),
     +           AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +           AA_rtPRG,AA_ltPRG,
     +           AA_new1,AA_new2,AA_new3
            write (161,'(A1)')char(13)

*               *****  write the record to the STA file:  *****

            write (162,'(9(''"'',A,''",''),
     +''"'',I5,''",'',9(''"'',A,''",''),
     +''"'',A,''"'',$)')
     +           exp_name,recording,date,BDT_FILE,
     +           qdtfilenames,blank200,blank200,blank200,
     +           cell_name(1),ID(1),
     +           STA_phrenic,STA_RLN,STA_cent_vagus,
     +           STA_lumbar,STA_cerv_symp,STA_ELN,
     +           STA_splanch,
     +           STA_new1,STA_new2,STA_new3
            write (162,'(A1)')char(13)

*               *****  write the record to the RESPONSES_xassist.txt file:  *****

            write (20,'(9(''"'',A,''",''),
     +''"'',I5,''",'',34(''"'',A,''",''),
     +''"'',A,''"'',$)')
     +           exp_name,recording,date,BDT_FILE,
     +           qdtfilenames,blank200,blank200,blank200,
     +           cell_name(1),ID(1),
     +           carotidCO2_x,
     +           vertCO2_x,
     +           hypercap_5_O2_x,
     +           hypercap_5_air_x,
     +           hypercap_tbd_x,
     +           hypoxia_12_x,
     +           hypoxia_5_x,
     +           gasp_x,
     +           pertnew1_x,
     +           carotidCO2_2_x,
     +           lobel_x,
     +           aorta_cath_x,
     +           pe_x,
     +           carotid_occ_x,
     +           nitrop_x,
     +           vc_cath_x,
     +           dopamine_x,
     +           vagus_x,
     +           capsaicin_x,
     +           pinch_x,
     +           deep_pain_x,
     +           codeine_x,
     +           nalox_x,
     +           methyserg_x,
     +           mucimol_x,
     +           dexameth_x,
     +           noinf_x,
     +           hyperinf_x,
     +           hypervent_x,
     +           cgh_x,
     +           lcgh_x,
     +           exp_reflex_x,
     +           SLNcgh_x,
     +           sw_x,
     +           SLNsw_x
            write (20,'(A1)')char(13)

*       ***** write the record to the RESPONSES_tsfs.txt file: *****

            write (21,'(9(''"'',A,''",''),
     +''"'',I5,''",'',174(''"'',A,''",''),
     +''"'',A,''"'',$)')
     +           exp_name,recording,date,BDT_FILE,
     +           qdtfilenames,blank200,blank200,blank200,
     +           cell_name(1),ID(1),

     +           carotidCO2,vertCO2,hypercap_5_O2,
     +           hypercap_5_air,hypercap_tbd,hypoxia_12,
     +           hypoxia_5,gasp,pertnew1,carotidCO2_2,
     +           lobel,aorta_cath,pe,carotid_occ,
     +           nitrop,vc_cath,dopamine,vagus,
     +           capsaicin,pinch,deep_pain,codeine,
     +           nalox,methyserg,mucimol,dexameth,
     +           noinf,hyperinf,hypervent,cgh,
     +           lcgh,exp_reflex,SLNcgh,sw,SLNsw
            write (21,'(A1)')char(13)

         else if (included(i).le.0) then
            print *, "Cell ", IDs(i)
     +           ," not included so not written"
         else
            print *, "Data not analyzed for cell ", IDs(i)
     +           ," so not written"
         end if
        end do

*       ********************************************************************************
*       *******  write "extra" cells to the UNITS and RESPONSES_tsfs.txt files:  *******
*       ********************************************************************************

        call concat(DB_FILES,'_extra_cell_preDBG',pre_DBG_EXTRA,l)  !precursor file for indiv EXTRA cells
        inquire(FILE=pre_DBG_EXTRA(1:l)//char(0),
     +          EXIST=EXIST_extra_preDBG)
        if(EXIST_extra_preDBG.eqv..TRUE.)then
           print '(''writing info for EXTRA cells ...'')'
           open (29,FILE=pre_DBG_EXTRA(1:l)//char(0),
     +                     FORM='FORMATTED',
     +                     STATUS='OLD',ACCESS='DIRECT',
     +                     RECL=DBG_RECL)
           exp_name_200 = exp_name
           exp_name_250 = exp_name
           call filesize (pre_DBG_EXTRA, ifilesize)
           iCellsInFile = ifilesize/DBG_RECL
           do rec_num = 1, iCellsInFile

*               *****  read a record from the .pre_dbg file:  *****

              call read_preDBG(29,rec_num,ios,date,BDT_FILE_30,
     +          recording_x,hist_num,
     +          cell_name(1),ID(1),resp_type_1,aAP,aRL,adep,
     +          sort,comm,card,AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +          AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +          STA_phrenic,STA_RLN,STA_cent_vagus,
     +          STA_lumbar,STA_cerv_symp,STA_ELN,
     +          STA_splanch,STA_new1,STA_new2,STA_new3,
     +          carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +          hypercap_5_air_x,
     +          hypercap_tbd_x,hypoxia_12_x,hypoxia_5_x,gasp_x,
     +          lobel_x,aorta_cath_x,pe_x,carotid_occ_x,nitrop_x,
     +          vc_cath_x,dopamine_x,
     +          vagus_x,capsaicin_x,pinch_x,deep_pain_x,
     +          codeine_x,nalox_x,
     +          methyserg_x,mucimol_x,dexameth_x,noinf_x,
     +          hyperinf_x,hypervent_x,
     +          pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +          SLNcgh_x,exp_reflex_x,sw_x,
     +          SLNsw_x,carotidCO2,vertCO2,hypercap_5_O2,
     +          hypercap_5_air,hypercap_tbd,hypoxia_12,hypoxia_5,
     +          gasp,lobel,aorta_cath,pe,carotid_occ,
     +          nitrop,vc_cath,dopamine,vagus,capsaicin,pinch,
     +          deep_pain,codeine,nalox,methyserg,mucimol,dexameth,
     +          noinf,hyperinf,hypervent,pertnew1,carotidCO2_2,cgh,lcgh,
     +          SLNcgh,exp_reflex,sw,
     +          SLNsw,per_results,AA_results,STA_results)


              if(ios.eq.0)then                                  !write only after a successful read
                 print '(T3,''EXTRA cell ID = '',I4)',ID(1)

*       *****  write the record to the UNITS file:  *****
*       *****  INSERT "No data" for things like PMS, meanISI, etc *****

                 write (16,'(9(''"'',A,''",''),
     +              ''"'',I5,''",'',37(''"'',A,''",''),
     +              ''"'',A,''"'',$)')
c     +         exp_name,recording,date,BDT_FILE,
c     +         qdtfilenames,blank200,blank200,blank200,
     +          exp_name,recording,date,exp_name_200,
     +          exp_name_250,blank200,blank200,blank200,
     +          cell_name(1),ID(1),
     +          resp_type_1,aAP,aRL,adep,blank4,blank4,
     +          blank1,blank5,blank2,blank3,
     +          blank1,blank40,blank5,blank5,blank5,blank5,
     +          blank5,blank5,blank5,
     +          blank3,blank3,blank3,blank3,blank3,blank3,
     +          blank3,blank3,blank3,blank3,blank3,blank3,
     +          blank10,blank10,blank10,blank3,
     +          blank10,blank10,blank40
                write (16,'(A1)')char(13)

*       *****  write the record to the RESPONSES_xassist.txt file:  *****

                write (20,'(9(''"'',A,''",''),
     +              ''"'',I5,''",'',34(''"'',A,''",''),
     +              ''"'',A,''"'',$)')
c     +         exp_name,recording,date,BDT_FILE,
c     +         qdtfilenames,blank200,blank200,blank200,
     +          exp_name,recording,date,exp_name_200,
     +          exp_name_250,blank200,blank200,blank200,
     +          cell_name(1),ID(1),
     +          carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +          hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +          hypoxia_5_x,gasp_x,pertnew1_x,carotidCO2_2_x,
     +          lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +          nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +          capsaicin_x,pinch_x,deep_pain_x,codeine_x,
     +          nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +          noinf_x,hyperinf_x,hypervent_x,cgh_x,
     +          lcgh_x,exp_reflex_x,SLNcgh_x,sw_x,SLNsw_x
                write (20,'(A1)')char(13)

             end if
          end do
          close (29)
        end if





*       *********************************************************************
*       ***************  create the CROSSES database file:  *****************
*       *********************************************************************

          print '(''Writing CROSSES file ....'')'
*
*
        GET_REF: do REF = 1,(MAX_NUM_CODES-1)
           if(IDs(REF).eq.0)cycle
           GET_TAR: do TAR = REF+1,MAX_NUM_CODES
              if(IDs(TAR).eq.0)cycle
*
              WRITE_TO_PAIRS_DB: 
     +         if((write_anyway.eq.'Y'.and.included(REF).gt.0
     +             .and.included(TAR).gt.0).or.(
     +             (analyzed_pairs(IDs(REF),IDs(TAR)).eq.2).and.
     +           (analyzed_cells(IDs(REF)).eq.2).and.
     +           (analyzed_cells(IDs(TAR)).eq.2)))then           !info exists for this cell pair AND EACH OF ITS CELLS
*
              statcomm = ' '
              Q_pos = 0
              write_the_preDBP_file = 'n'
              if(included(REF).eq.1)goto 301
              do k = 1,included(REF)-1                          !calculate record location 
                 Q_pos = Q_pos + (total_num_cells-k)    
              end do
 301          Q_pos = Q_pos + (included(TAR) - included(REF))
              rec_num = Q_pos

*               ***** read the appropriate record from the .pre_dbp file:  *****

              call read_preDBP(10,rec_num,ios,recording_x,date,
     +          BDT_FILE_30,cell_name(1),ID(1),cell_name(2),ID(2),
     +          prim,loc1,sec,loc2,rescom,czk,cprobk,
     +          cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +          hist_num)

         if((prim.eq.'Flat').or.(prim.eq.'M P & T'))then        !make sure that no statistics for a non-existent feature is
            cdet=' '                                            ! included in the database files
            statcomm = ' '
            stats_bw= ' '
            cprobk = ' '
            czk = ' '
            cvis = ' '
            czlat = ' '
            chalfwd = ' '
         end if
         if((cell_name(1).eq.' ').and.
     +      (CELL_NAMES(IDs(REF)).ne.' '))then
            cell_name(1) = CELL_NAMES(IDs(REF))
            write_the_preDBP_file = 'y'
           end if
         if((cell_name(2).eq.' ').and.
     +      (CELL_NAMES(IDs(TAR)).ne.' '))then
               cell_name(2) = CELL_NAMES(IDs(TAR))
               write_the_preDBP_file = 'y'
           end if

           if(write_the_preDBP_file.eq.'y')
     +          call write_preDBP(10,rec_num,ios,recording,date,
     +          BDT_FILE_30,cell_name(1),ID(1),cell_name(2),ID(2),
     +          prim,loc1,sec,loc2,rescom,czk,cprobk,
     +          cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +          hist_num)
           if((cell_name(1).eq.' ').or.(cell_name(2).eq.' '))          !do not write info if either cell not named
     +          cycle GET_TAR                                          ! -- get the next pair

*               ***** write the record to the CROSSES file:  *****

             write (17,'(8(''"'',A,''",''),
     +          2(''"'',A,''",'',  ''"'',I5,''",''),11(''"'',A,''",''),
     +          ''"'',A,''"'',$)')
     +            exp_name,recording,date,BDT_FILE,qdtfilenames,
     +            blank200,blank200,blank200,
     +            cell_name(1),ID(1),cell_name(2),ID(2),
     +            prim,loc1,sec,loc2,rescom,
     +            czk,cprobk,cdet,cvis,statcomm,czlat,chalfwd
             write (17,'(A1)')char(13)

*     ***** write the binwidth used to calculate the statistics to stats_bw.txt: *****

             write (172,'(3(''"'',A,''",''),2(''"'',I5,''",''),
     +           ''"'',A,''"'',$)')
     +            exp_name,recording,date,ID(1),ID(2),stats_bw
             write (172,'(A1)')char(13)

*               ***** write the record to the EFFECTIVENESS file:  *****

             if(prim.eq.'Peak'.and.loc1.eq.'Offset right')then
                read (cdet,'(F6.2)') det
                read (cvis,'(F6.2)') vis
c                eff=det**2/vis/ITAL(IDs(ID(1)))
                print '(''det = '',f6.2)',det
                print '(''vis = '',f6.2)',vis
                print '(''eff = '',f7.5)',eff
                write (171,'(3(''"'',A,''",''),
     +          2(''"'',I5,''",''),''"'',F7.5,''"'',$)')
     +            exp_name,recording,date,ID(1),ID(2),eff
                write (171,'(A1)')char(13)
             end if

             else
                print *,"Data not analyzed for pair ",IDs(REF),">"
     +               ,IDs(TAR)," so not written"
           end if WRITE_TO_PAIRS_DB
          end do GET_TAR
         end do GET_REF

*       ************************************************************************************
*       ***************  create the ELECTRODE REFERENCES database file:  *******************
*       ************************************************************************************

          print '(''Writing ELECTRODE REFERENCES file ....'')'

          do i = 1, total_num_electrode_arrays
             if(NAME_refs(i).ne.' ')then
                write (4,'(7(''"'',A,''",''),''"'',A,''"'',$)')
     +                  exp_name,recording,date,
     +                  NAME_refs(i),AP_refs(i),RL_refs(i),
     +                  depth_refs(i),ref_chan_refs(i)
                write (4,'(A1)')char(13)
             end if
          end do
*                                       * = DB_FILES
        close (unit=3)          !       *.ARRAY_REFS.sav
        close (unit=4)          !       *.ARRAY_REFS.txt
        close (unit=9)          !       *.preDBG
        close (unit=10)         !       *.preDBP
        close (unit=15)         !       *.DBSAV
        close (unit=16)         !       *.DBG = "UNITS"
        close (unit=17)         !       *.DBP = "CROSSES"
        close (unit=171)        !       *.DBP_eff
        close (unit=172)        !       *.DBP_statsbw
        close (unit=20)         !       *.DBR_xass = "RESPONSES" from xassist
        close (unit=21)         !       *.DBR_tsfs = "RESPONSES" from tsfs
        close (unit=161)        !       *.DBG_AA = "AA" data for UNITS
        close (unit=162)        !       *.DBG_STA = "STA" data for UNITS
        close (unit=164)        !       *.DBG_NUM_EVENTS

*       ***** create a copy of the .txt files on Oberon for GAIA database and experiment folder: *****

        print '(/,''Copying experiment files to Oberon ....'')'

        call copy_tsfs_files (BDT_FILE)
        call write_sta_file (exp_name)
        call write_persav_file (DB_FILES, per_prev_imported)

        call strlength(exp_name,LEN(exp_name),l)
c       print '(''exp_name: '',A)', exp_name
        if(l.eq.3)then                                          !might be an "old" experiment name (e.g., k15) -
           if((((exp_name(1:1).ge.'A').and.(exp_name(1:1).le.'Z')) !so append "m" and the recording number
     +           .or.((exp_name(1:1).ge.'a').and.
     +                (exp_name(1:1).le.'z')))
     +           .and.((exp_name(2:2).ge.'0').and.
     +                 (exp_name(2:2).le.'9'))
     +           .and.((exp_name(3:3).ge.'0').and.
     +                  (exp_name(3:3).le.'9')))then

c                   dir_name=exp_name(1:3)//'m'//recording(1:1)
                    dir_name=exp_name//'m'//recording
              call remove_all_blanks(dir_name,LEN(dir_name))
              print '(''dir_name = '',A)',dir_name
           end if
        else                                                    !this is a "new" experiment name (YYYY-MM-DD_xxx)
           call strlength(exp_name,LEN(exp_name),l)
           dir_name = exp_name(1:l)
        end if
        call strlength(dir_name,LEN(dir_name),l)
        print '(''dir_name = '',A)',dir_name

        do i = 1, total_num_qdts
           if(qdtfiles(i).ne.' ')then
              call strlength(qdtfiles(i),LEN(qdtfiles(i)),m)
c             print '(''qdtfiles('',I1,'') = '',A)',i,qdtfiles(i)(1:m)
              isys=SYSTEM('cp '//qdtfiles(i)(1:m)//
     +          ' /oberon/experiments/'//                              !QDT_FILES
     +          dir_name(1:l)//'/'//qdtfiles(i)(1:m)//char(0))
              call concat(qdtfiles(i)(1:m),'.txt',FILE_1,n)
c             print '(''QDTTXT = '',A)',FILE_1(1:n)
              isys=SYSTEM('todos < '//FILE_1(1:n)
     +          //' > /oberon/experiments/'//dir_name(1:l)//'/'
     +          //FILE_1(1:n)//char(0))
              call concat(qdtfiles(i)(1:m),'.sav',FILE_1,n)
c             print '(''QDTSAV = '',A)',FILE_1(1:n)
              isys=SYSTEM('cp '//FILE_1(1:n)//' /oberon/experiments/'//!QDTSAV
     +          dir_name(1:l)//'/'//FILE_1(1:n)//char(0))
           end if
        end do

        call strlength(DBSAV,LEN(DBSAV),m)
        isys=SYSTEM('cp '//DBSAV(1:m)//' /oberon/experiments/'//       !DBSAV
     +          dir_name(1:l)//'/'//DBSAV(1:m)//char(0))
        call strlength(pre_DBG,LEN(pre_DBG),m)
        isys=SYSTEM('cp '//pre_DBG(1:m)//' /oberon/experiments/'//     !preDBG file
     +          dir_name(1:l)//'/'//pre_DBG(1:m)//char(0))
        call strlength(pre_DBP,LEN(pre_DBP),m)
        isys=SYSTEM('cp '//pre_DBP(1:m)//' /oberon/experiments/'//     !preDBP file
     +          dir_name(1:l)//'/'//pre_DBP(1:m)//char(0))

        call concat(DB_FILES,'_UNIT_REFS.sav',FILE_1,m)
        isys=SYSTEM('cp '//FILE_1(1:m)//' /oberon/experiments/'//      !*_UNIT_REFS.sav
     +          dir_name(1:l)//'/'//FILE_1(1:m)//char(0))
        call concat(DB_FILES,'_ARRAY_REFS.sav',FILE_1,m)
        isys=SYSTEM('cp '//FILE_1(1:m)//' /oberon/experiments/'//      !*_ARRAY_REFS.sav
     +          dir_name(1:l)//'/'//FILE_1(1:m)//char(0))
        
*       ******** COPY DATASAVE FILES TO OBERON ********

        print '(/,''Copying database files to Oberon ....'')'
        if(which_subset.eq.'G')then
           print '(''copying GASPING text files to oberon ...'')'
           isys=SYSTEM('cp '//EL_REF(1:m)//
     +      ' "/oberon/databases/GAIA Database/ARRAY_REFS_gasp.txt"'
     +          //char(0))
           call strlength(DBG,LEN(DBG),m)
           isys=SYSTEM('cp '//DBG(1:m)//
     +     ' "/oberon/databases/GAIA Database/UNITS_gasp.txt"'//char(0))
           call strlength(DBG_AA,LEN(DBG_AA),m)
           isys=SYSTEM('cp '//DBG_AA(1:m)//
     +      ' "/oberon/databases/GAIA Database/UNITS_AA_gasp.txt"'
     +          //char(0))
           call strlength(DBG_STA,LEN(DBG_STA),m)
           isys=SYSTEM('cp '//DBG_STA(1:m)//
     +      ' "/oberon/databases/GAIA Database/UNITS_STA_gasp.txt"'
     +          //char(0))
           call strlength(DBP,LEN(DBP),m)
           isys=SYSTEM('cp '//DBP(1:m)//
     +      ' "/oberon/databases/GAIA Database/CROSSES_gasp.txt"'
     +          //char(0))
           call strlength(DBR_xass,LEN(DBR_xass),m)
           isys=SYSTEM('cp '//DBR_xass(1:m)//
     +     ' "/oberon/databases/GAIA Database/RESPONSES_xass_gasp.txt"'
     +          //char(0))
           call strlength(DBR_tsfs,LEN(DBR_tsfs),m)
           isys=SYSTEM('cp '//DBR_tsfs(1:m)//
     +     ' "/oberon/databases/GAIA Database/RESPONSES_tsfs_gasp.txt"'
     +          //char(0))
         elseif(which_subset.eq.'I')then
           print '(''copying NO INFLATION text files to oberon ...'')'
           call strlength(DBG,LEN(DBG),m)
           isys=SYSTEM('cp '//DBG(1:m)//
     +     ' "/oberon/databases/GAIA Database/UNITS_no_inf.txt"'
     +          //char(0))
        elseif(which_subset.eq.'V')then
           print '(''copying PRE-VAGOTOMY text files to oberon ...'')'
c          isys=SYSTEM('cp '//EL_REF(1:m)//
c     +     ' "/oberon/databases/GAIA Database/ARRAY_REFS_prevx.txt"'//char(0))
           call strlength(DBG,LEN(DBG),m)
           isys=SYSTEM('cp '//DBG(1:m)//
     +      ' "/oberon/databases/GAIA Database/UNITS_prevx.txt"'
     +          //char(0))
           call strlength(DBG_AA,LEN(DBG_AA),m)
           isys=SYSTEM('cp '//DBG_AA(1:m)//
     +      ' "/oberon/databases/GAIA Database/UNITS_AA_prevx.txt"'
     +          //char(0))
           call strlength(DBG_STA,LEN(DBG_STA),m)
           isys=SYSTEM('cp '//DBG_STA(1:m)//
     +      ' "/oberon/databases/GAIA Database/UNITS_STA_prevx.txt"'
     +          //char(0))
           call strlength(DBP,LEN(DBP),m)
           isys=SYSTEM('cp '//DBP(1:m)//
     +      ' "/oberon/databases/GAIA Database/CROSSES_prevx.txt"'
     +          //char(0))
           call strlength(DBR_xass,LEN(DBR_xass),m)
           isys=SYSTEM('cp '//DBR_xass(1:m)//
     +    ' "/oberon/databases/GAIA Database/RESPONSES_xass_prevx.txt"'
     +          //char(0))
           call strlength(DBR_tsfs,LEN(DBR_tsfs),m)
           isys=SYSTEM('cp '//DBR_tsfs(1:m)//
     +    ' "/oberon/databases/GAIA Database/RESPONSES_tsfs_prevx.txt"'
     +          //char(0))
        else                                    !write "normal" files to oberon
           isys=SYSTEM('cp '//EL_REF(1:m)//
     +          ' "/oberon/databases/GAIA Database/ARRAY_REFS.txt"'
     +          //char(0))
           call strlength(DBG,LEN(DBG),m)
           isys=SYSTEM('cp '//DBG(1:m)//
     +          ' "/oberon/databases/GAIA Database/UNITS.txt"'//char(0))
           call strlength(DBG_AA,LEN(DBG_AA),m)
           isys=SYSTEM('cp '//DBG_AA(1:m)//
     +          ' "/oberon/databases/GAIA Database/UNITS_AA.txt"'
     +          //char(0))
           call strlength(DBG_STA,LEN(DBG_STA),m)
           isys=SYSTEM('cp '//DBG_STA(1:m)//
     +          ' "/oberon/databases/GAIA Database/UNITS_STA.txt"'
     +          //char(0))
           call strlength(DBG_num_events,LEN(DBG_num_events),m)
           isys=SYSTEM('cp '//DBG_num_events(1:m)//
     +          ' "/oberon/databases/GAIA Database/NUM_EVENTS.txt"'
     +          //char(0))
           call strlength(DBP,LEN(DBP),m)
           isys=SYSTEM('cp '//DBP(1:m)//
     +          ' "/oberon/databases/GAIA Database/CROSSES.txt"'
     +          //char(0))
           call strlength(DBP_eff,LEN(DBP_eff),m)
           isys=SYSTEM('cp '//DBP_eff(1:m)//
     +          ' "/oberon/databases/GAIA Database/EFFECTIVENESS.txt"'
     +          //char(0))
           call strlength(DBP_statsbw,LEN(DBP_statsbw),m)
           isys=SYSTEM('cp '//DBP_statsbw(1:m)//
     +          ' "/oberon/databases/GAIA Database/STATS_BW.txt"'
     +          //char(0))
           call strlength(DBR_xass,LEN(DBR_xass),m)
           isys=SYSTEM('cp '//DBR_xass(1:m)//
     +    ' "/oberon/databases/GAIA Database/RESPONSES_xass.txt"'
     +          //char(0))
           call strlength(DBR_tsfs,LEN(DBR_tsfs),m)
           isys=SYSTEM('cp '//DBR_tsfs(1:m)//
     +    ' "/oberon/databases/GAIA Database/RESPONSES_tsfs.txt"'
     +          //char(0))
        end if

        print '(//,''Finished writing database files!'')'

        return
        end
      end module mod_database_write_new
