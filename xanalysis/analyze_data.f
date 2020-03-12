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

      module mod_analyze_data
      contains
*       file = XAnalysis/analyze_data.f 
*
*       date of latest revision = 21-Sep-2006   lss
*
*     aug-2003          lss
*        work done on importation of *.info and *.per data files
*
*       nov-2002        lss & ro'c      
*        subroutine screensize added
*
*       aug-2001        lss
*       changes made to names of output files (eg., *_dbg.txt)
*
*       sep-2000        lss
*        x2000 now allows users to simply view CCH windows; no database files
*          are written or saved
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*       *** INDIRECT POINTERS are now used to access unit data
*       ***  within the following arrays:
*       ***     SPIKETIMES, ITAL, card_type, coef, coefnum, TALLY_NORM,
*       ***     ETA2_*, zmodsig_*, zmodsig2_*, analyzed_cells, analyzed_pairs,
*       ***     resp_type, card, CELL_NAMES
*       ***     [these arrays are now dimensioned to MAX_NUM_CHAN;
*       ***      the location of a cell's information (user-assigned ID code = i)
*       ***      within these arrays is array(IDs())]
*       ***
*       *** DIRECT POINTERS are used with these arrays: excluded, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       link with x2002 code
*
*       new_plot incorporated   01-apr-99
*
*
*
*       Format of .qdt file:    [# cells = N; # pairs = (N(N-1)/2)]
*
*               histogram #     1       phrenic overlay         total # overlays = 3
*                               2       normalized phrenic overlay
*                               3       cardiac overlay
*               -----------------------------------------------------
*                               4       CTHs for 1st cell:      respiratory CTH
*                               5                               norm. resp. CTH
*                               6                               cardiac CTH
*                               .
*                               .                       total # CTHs = 3N
*                               .
*                            3N+3       cardiac CTH for cell N
*               -------------------------------------------------------
*                            3N+4       ACHs for 1st cell:      @ binwidth 1
*                            3N+5                               @ binwidth 2
*                            3N+6                               @ binwidth 3
*                            3N+7                               @ binwidth 4
*                               .
*                               .                       total # ACHs = 4N
*                               .
*                            7N+3       ACH for cell N @ binwidth 4
*               -------------------------------------------------------
*                            7N+4       CCHs for 1st pair:      @ binwidth 1
*                            7N+5                               @ binwidth 2
*                            7N+6                               @ binwidth 3
*                            7N+7                               @ binwidth 4
*                               .
*                               .                       total # CCHs = 2(N**2)-2N
*                               .
*                    2(N**2)+5N+3       CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*                    2(N**2)+5N+4       control CCHs for 1st pair:      @ binwidth 1
*                    2(N**2)+5N+5                                       @ binwidth 2
*                    2(N**2)+5N+6                                       @ binwidth 3
*                    2(N**2)+5N+7                                       @ binwidth 4
*                               .
*                               .                       total # control CCHs = 2(N**2)-2N
*                       .
*                    4(N**2)+3N+3       control CCH for last pair @ binwidth 4
*               -------------------------------------------------------
*                                
*               total # of histograms written to .qdt = 4(N**2) + 3N + 3
*
*       ************************************************************************************
*
*               for cell n within the .qdt file:        [where N=total # cells]                          
*                                                                                
*                       location of resp CTH = (((included(n)*3)+1)      
*                       location of ACH @ binwidth 1 = (((included(n)*4)+3N      
*
*                                                                                
*          To determine the location (record number) of the CCH @ binwidth 1:  
*                                                                                
*               record number = 3 + #CTHs + #ACHs + relative position of CCH w/in all CCHs
*                                                                                
*               1. determine the position of the cell pair w/in the pair queue:  
*                  (sample queue(REF,TAR): 1,2  1,3  1,4  2,3  2,4  3,4)  N=4 in this case
*        
*                1st occurrence of REF cell in the queue=(N-1)+(N-2)+..+(N-included(REF)-1))+1
*                               Q_pos = 0                                                        
*                               if(included(REF).eq.1) goto >A<          
*                               do i=1,included(REF)-1                           
*                                 Q_pos = Q_pos + (N-i)                          
*                               end do                                                   
*        
*          relative position of TAR w/in the group of REF cells=included(TAR)-included(REF)
*        
*                       >A<             Q_pos = Q_pos + (included(TAR)-included(REF))
*                                                                                
*                       so Q_pos now equals the position of the pair within the queue
*
*                (N.B.-- could use the variable 'pair_counter' as an index of the pair's 
*                        position within the pair queue, but that would not be valid when the
*                        user has chosen to view CCHs according to direct access)
*                                                                                
*               2. now use the queue position to determine the relative position of      
*                       the 1st CCH for a pair among all CCHs:                           
*                                                                                
*                       rel_loc = (Q_pos *4)-3                                   
*                                                                                
*               SO..........                                             
*                       record number = 3 + 3N + 4N + rel_loc = 3 + 7N + rel_loc
*                                                                                
*          To determine the record number of a control CCH @ binwidth 1:         
*                                                                                
*               1. find Q_pos as above                                   
*               2. calculate rel_loc as above                    
*               3. record number=3+#CTHs+#ACHs+#CCHs+relative position of control CCH
*                  record number=3+3N+4N+(2(N**2)-2N)+rel_loc
*                               = 3 + 5N + 2(N**2) + rel_loc
*                                                                                
*
      
      SUBROUTINE analyze_data(pgm_version,mode,
     +     BDT_FILE,qdt_files,
     +     DB_FILES,
     +     import_info,INFO_FILENAME,import_per,
     +     PER_FILENAME,c_format,total_num_qdts,
     +     date_exp,recording,protocol,exp_name)
      use mod_clear_routines
      use mod_compare_and_convert
      use mod_enter_date
      use mod_get_info_IDcodes
      use mod_miscellaneous_subroutines
      use mod_picktypetest2
      use mod_plot_spikes_per_cycle
      use mod_print_and_write_routines
      use mod_read_and_write_DBSAV
      use mod_read_and_write_pre_files
      use mod_read_and_write_ref_electrode_data
      use mod_read_info_file
      use mod_showCCHs6_newshift
      use mod_showCTHs_v3

      include 'x2000parameter.defs'
      integer total_num_qdts
      integer IDs(MAX_NUM_CODES),IDs_new(MAX_NUM_CODES),
     +     ITAL(MAX_NUM_CHAN),IDs_extra(MAX_NUM_CODES),
     +     IDs_info(MAX_NUM_CODES)

      integer cardiac_hist(NUM_BINS), phrenic_hist(NUM_BINS),
     +     CTH(NUM_BINS), CCH(NUM_BINS), CONTROL(NUM_BINS), 
     +     CONTROL_AVG(NUM_BINS),
     +     CONTROL_2(NUM_BINS),
     +     excluded(MAX_NUM_CODES),
     +     excluded_from_QDT(MAX_NUM_CODES),
     +     excluded_from_QDT_new(MAX_NUM_CODES),
     +     excluded_from_DB(MAX_NUM_CODES),
     +     included(MAX_NUM_CODES),
     +     included_in_QDT(MAX_NUM_CODES),
     +     included_in_QDT_new(MAX_NUM_CODES),
     +     included_in_DB(MAX_NUM_CODES),
     +     norm_phrenic_hist(NUM_BINS),MAX_NUM_CHAN_SQ,
     +     norm_CTH(NUM_BINS),REF_ACH(NUM_BINS),TAR_ACH(NUM_BINS),
     +     cardCCH(NUM_BINS),
     +     perturb(MAX_PERTURB),perturbation,
     +     TAR_ACH_1(NUM_BINS),
     +     cardiac_REF(NUM_BINS),cardiac_TAR(NUM_BINS),
     +     CTH_REF(NUM_BINS),CTH_TAR(NUM_BINS),
     +     valid_pair_REF(MAX_NUM_PAIRS),
     +     valid_pair_TAR(MAX_NUM_PAIRS),pair_counter,
     +     analyzed_cells(MAX_NUM_CHAN),
     +     analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +     width,height1,STATUS,TOTAL_STIM,
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),ICN(6),
     +     qdt_for_ACHs

      real BINW,BINW2,STIM_OFFSET,NORM_BINW,
     +     NORM_OFFSET,coefnum(MAX_NUM_CHAN),
     +     min_bw,min_bw_1(total_num_qdts),mean_E,NORM_BW

      integer REF,TAR,E_pulse, cardiac_pls, phrenic, BNDRY, 
     +     histogram_number,histogram_number_x,
     +     qdttxt_count,Revents,Revents_p,
     +     Revents_n,Revents_c,Revents_ach1,Revents_ach2,
     +     Tevents,cell,
     +     total_histograms,I_pulse,REC_num,
     +     Rcode,Tcode,REC_num_x,REC_num_xx,total_num_cells,
     +     total_num_cells_new,Q_pos,rel_loc,
     +     REFtemp,TARtemp,IDs_1(2),REC_num_control_1,flag,
     +     IDs_per(MAX_NUM_CODES),total_num_cells_db,
     +     current_qdt,REC_num_control_avg,
     +     REC_num_control_2,SPC_array,TOTAL_NUM_SHIFTS,
     +     total_num_sig,total_num_pairs

      integer*4 mouse,isys
      INCLUDE 'gopen_type.defs'

      logical*4 exist,EXIST_extra_predbg,exist10

      character*(*) exp_name,protocol,BDT_FILE,
     +     qdt_files(MAX_NUM_QDTS),DB_FILES,
     +     INFO_FILENAME,PER_FILENAME

      character*250 qdtfilenames
      character*200 qdtfilenames_200
      character*120 text

      character*50 QDT_FILENAME,dirname,
     +     files(MAX_NUM_QDTS+10)
      CHARACTER*40 comm, rescom,statcomm,char40,
     +     comm_import(MAX_NUM_CHAN)
      character*40 AA_text(MAX_AA),
     +     STA_text(MAX_STA)
      character*30 BDT_FILE_30
      character*12 STA_text1(MAX_STA),AA_text1(MAX_AA),
     +     per_text_abbrev(MAX_PERTURB)

      character*20 c_format
      character*30  textstring,char30

      CHARACTER*20 prim,sec,char20,prim1,prim_old

      CHARACTER*15 loc1,loc2,pattern,
     +     resp_type(MAX_NUM_CHAN),button_choice,
     +     STA_res_text(122),char15,
     +     resptype

      character*12 STA_phrenic,STA_RLN,STA_lumbar,
     +     STA_ELN,STA_cent_vagus,STA_splanch,
     +     STA_cerv_symp,STA_import(MAX_NUM_CHAN,MAX_STA),
     +     STA_new1,STA_new2,STA_new3,STA_results(MAX_STA),
     +     char12,char12x(MAX_STA)

      character*9 WINDOW1
      character*5 WINDOW2
c$$$        character*10 USER

      character*11 date_exp,date_exp_new,date_import,char11 !DD-MMM-YYYY -- 4-digit year required
      character*175 expname_import

      character*220 DBG,        !"UNITS" text file for later importation into MSAccess database
     +     DBG_AA,              !"AA" text file for later importation into MSAccess database
     +     DBG_STA,             !"STA" text file for later importation into MSAccess database
     +     DBP,                 !"CROSSES"      "       "       "
     +     DBR_xass,DBR_tsfs,   !"RESPONSES"    "       "       "
     +     DBSAV,
     +     pre_DBG,pre_DBP,
     +     DBSAV_BAK,
     +     pre_DBG_BAK,pre_DBP_BAK,MEAN__E
      character*250 pre_DBG_EXTRA
      character*60 QDTSAV,qdt_sav(MAX_NUM_QDTS)
      character*10 meanISI(MAX_NUM_CHAN),sdISI(MAX_NUM_CHAN),
     +     fiveHT(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +     sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +     num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT

      CHARACTER*9 per_res(MAX_NUM_CHAN,MAX_PERT),char9

      CHARACTER*8 czk, cprobk,cdet, cvis,
     +     chalfwd, czlat,stats_bw,char8,c_num_sig

      character*9 carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,
     +     dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +     sw_x,SLNsw_x

      character*9 per_results(MAX_PERTURB),
     +     per_results_REF(MAX_PERTURB),
     +     per_results_TAR(MAX_PERTURB),char9x(MAX_PERTURB)

      character*12 bwtext1,bwtext2,bwtext3,bwtext4,bwtext5

      character*8 bwtext

      CHARACTER*5 aAP, aRL,adep,REF_coords(3),TAR_coords(3),
     +     AP_refs(MAX_NUM_ARRAYS),RL_refs(MAX_NUM_ARRAYS),
     +     depth_refs(MAX_NUM_ARRAYS),char5,
     +     coef(MAX_NUM_CHAN),ETA2_1(MAX_NUM_CHAN),
     +     ETA2_2(MAX_NUM_CHAN),ETA2_3(MAX_NUM_CHAN),
     +     ETA2_4(MAX_NUM_CHAN),ETA2_5(MAX_NUM_CHAN),
     +     ETA2_6(MAX_NUM_CHAN),AP_import(MAX_NUM_CHAN),
     +     RL_import(MAX_NUM_CHAN),depth_import(MAX_NUM_CHAN),sort1,
     +     DELTA2(MAX_NUM_CHAN)

      CHARACTER*4 CELL_NAMES(MAX_NUM_CHAN),
     +     cell_name(2),name_import(MAX_NUM_CHAN),
     +     char4,NAME_refs(MAX_NUM_ARRAYS),
     +     dchan(MAX_NUM_CHAN),
     +     ref_chan(MAX_NUM_CHAN),ref_chan_refs(MAX_NUM_ARRAYS)

      CHARACTER*3 AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +     AA_rVRG,AA_rtPRG,AA_ltPRG,
     +     AA_import(MAX_NUM_CHAN,MAX_AA),
     +     AA_new1,AA_new2,AA_new3,AA_results(MAX_AA),
     +     card1,c_i,c_j,
     +     zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),zmodsig2_6(MAX_NUM_CHAN),
     +     card_type(MAX_NUM_CHAN),char3,char3x(MAX_AA),c_percent

      character*2 carotidCO2(5),vertCO2(5),hypercap_5_O2(5),
     +     hypercap_5_air(5),hypercap_tbd(5),hypoxia_5(5),
     +     hypoxia_12(5),gasp(5),
     +     lobel(5),aorta_cath(5),pe(5),carotid_occ(5),
     +     nitrop(5),vc_cath(5),dopamine(5),vagus(5),
     +     capsaicin(5),pinch(5),deep_pain(5),
     +     codeine(5),nalox(5),methyserg(5),mucimol(5),
     +     dexameth(5),
     +     noinf(5),hyperinf(5),hypervent(5),pertnew1(5),
     +     carotidCO2_2(5),cgh(5),lcgh(5),SLNcgh(5),
     +     exp_reflex(5),
     +     sw(5),SLNsw(5),
c     +       tsfs_REF(MAX_PERTURB,5),
c     +       tsfs_TAR(MAX_PERTURB,5),
     +     recording,recording_new,rec_import,recording_x,
     +     file_version,mode,tedfactor(MAX_NUM_CHAN),
     +     mode_original,char2,char2x(5),pgm_version,QDT_version

      CHARACTER*1 OK,menu_choice,old_choice,
     +     perturb_applied,AA_applied,STA_applied,
     +     STA(MAX_STA),AA(MAX_AA),card(MAX_NUM_CHAN),
     +     auto_mode,click,cardd,
     +     cells_completed,pairs_completed,screen_CTHs,
     +     screen_CCHs,BOUNDARY,
     +     show_control,import_info,import_per,
     +     info_prev_imported,per_prev_imported,sort,
     +     single_shift,avg_shift,char1,
     +     show_single_shift,show_avg_shift,task,
     +     show_2_sd,show_3_sd,show_conf_lim,enter_info,IDmismatch,
     +     warning,warning_date,warning_rec,warning_ID,warning_inc,
     +     warning_exc,warning_num_cells,sig_only,CTHsaved,
     +     extra,flats_only,esc,scaledMAIN,changed_comment,
     +     show_pulses

c        STRUCTURE /EXTRAS/
c           character*4 ID
c           character*4 cellname
c           character*5 AP
c           character*5 RL
c           character*5 depth
c           character*4 dchan
c           character*4 refchan
c           character*3 AA(MAX_AA)
c           character*12 STA(MAX_STA)
c           character*40 comment
c        END STRUCTURE
c        RECORD /EXTRAS/ EXTRAS(MAX_NUM_CHAN)

*     ***** temp files for read and write of pre-dbg files if importing: *****
      character*5 AP_temp,RL_temp,depth_temp
      character*1 sort_temp
      character*40 comm_temp
      character*4 names_temp
      
*
      INTEGER*4 SYSTEM
      EXTERNAL SYSTEM,RENAME,CHDIR

      if(.false.)print *,fildes2 !suppress unused variable warning
*       
*
*
*       *****************************************************************
*       *****************************************************************
*
*
      id = 0
      sig_only = 'n'
      total_num_pairs = 0
      total_num_cells_db = 0
      warning = 'n'
      warning_date = 'n'
      warning_rec = 'n'
      warning_num_cells = 'n'
      warning_ID = 'n'
      warning_exc = 'n'
      warning_inc = 'n'
      qdt_for_ACHs = 0
      old_choice = ' '

      call gerr_print_control(NO_ERROR_PRINTING)
*
      WINDOW1='XAnalysis'                             
      WINDOW2='Q_SUM'
      call screensize(width,height1,STATUS)
      if(STATUS.eq.0) then
         print '(''cannot open display'')'
         return
      end if
c$$$        write (c_width,'(I6)') width

*       * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                               *
*       *       initialize variables and arrays:        *
*       *                                               *
*       * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
c      print *,'start analyze_data'
      do i = 1, MAX_NUM_CHAN
         CELL_NAMES(i)=' '
      end do
      SPC_array = MAX_NUM_CHAN*MAX_NUM_ACC_CYCLES
      BINW=0.0
      do i = 1, NUM_BINS
         CTH(i)=0               !array that holds:  resp cycle-triggered histograms
         norm_CTH(i)=0          !            normalized respiratory CTHs
         CCH(i)=0               !            cross-correlograms
         phrenic_hist(i)=0      !            phrenic overlay CTH
         norm_phrenic_hist(i)=0 !            normalized phrenic overlayCTH
         cardiac_hist(i)=0      !            cardiac CCH
         cardCCH(i)=0           !            cardiac overlay CCH
         cardiac_REF(i)=0
         cardiac_TAR(i)=0
         CONTROL(i)=0           !            shift-control cross-correlogram
         CONTROL_AVG(i)=0       !            shift-control cross-correlogram
         REF_ACH(i)=0           !            reference cell ACH
         TAR_ACH(i)=0           !            target cell ACH
         TAR_ACH_1(i)=0         !            target cell ACH (@ next-to-greatest BW)
      end do
      do i = 1, MAX_NUM_CODES
         included(i) = 0
         included_in_DB(i) = 0
         included_in_QDT(i) = 0
         excluded(i) = 0
         excluded_from_DB(i) = 0
         excluded_from_QDT(i) = 0
         IDs_info(i) = 0
         IDs_per(i) = 0
         IDs_extra(i)= 0
      end do
      NUM_STIM =0               !number of stimuli used to compute CTHs
      STIM_OFFSET=0.0           !offset applied to resp CTHs
      NORM_OFFSET=0.0
      histogram_number=0
      enter_info=' '
      qdttxt_count=0
      do i = 1, MAX_NUM_PAIRS
         valid_pair_REF(i)=0    !used together, these arrays contain ref and tar codes
         valid_pair_TAR(i)=0    ! of cells whose CCHs will be displayed
                                !  e.g., valid_pair_REF(i) and valid_pair_TAR(i) hold
                                !    the ID codes of the REF and TAR cells for the pair
                                !    that is in position i in the pair queue
      end do
      click = 'n'

      do i = 1, MAX_NUM_CHAN
         analyzed_cells(i)=0          
      end do
      do j = 1, MAX_NUM_CHAN
         do i = 1, MAX_NUM_CHAN
            analyzed_pairs(i,j)=0     
         end do                       
      end do                          

      prim = ' '
      sec = ' '
      loc1 = ' '
      loc2 = ' '
      stats_bw = ' '
      czk=' '
      cprobk=' '
      cdet=' '
      cvis=' '
      czlat=' '
      chalfwd=' '
      rescom=' '
      statcomm=' '
      MAX_NUM_CHAN_SQ=MAX_NUM_CHAN**2
      screen_CTHs='n'
      screen_CCHs='n'
      pre_DBG_EXTRA = ' '                     
      pre_DBG = ' '                   
      pre_DBP = ' '
      DBG=' '
      DBG_AA=' '
      DBG_STA=' '
      DBP=' '
      DBR_xass=' '
      DBR_tsfs=' '
      DBSAV=' '
      OK=' '
      date_exp=' '              !       date of experiment
      exp_name=' '              !       name of experiment
      protocol=' '              !       experimental protocol
      recording=' '             !       recording#
      pattern=' '               !       cell's respiratory modulation
*                               !corresponding column names in Access
*                                 database 'CELLS': 
*
      comm=' '                  !       comment, comment1, comment2
      aAP  =' '                 !       anat_AP, an_AP_1, an_AP_2
      aRL  =' '                 !       anat_RL, an_RL_1, an_RL_2
      adep  =' '                !       anat_dep, an_dep_1, an_dep_2
      sort  =' '                !       c=clean, m=messy sort
      do i = 1, MAX_NUM_CHAN
         resp_type(i)=' '       !       cell's respiratory firing pattern
         card(i)=' '            !       c = card related visual inspect.
         card_type(i)=' '       !         cell's cardiac firing pattern
         num_rej_ISI(i)= ' '
         num_rej_rISI(i) = ' '
         zmodsig_1(i)  =' '     !       anova resp mod sig values
         zmodsig_2(i)  =' '
         zmodsig_3(i)  =' '
         zmodsig_4(i)  =' '
         zmodsig_5(i)  =' '
         zmodsig_6(i)  =' '
         zmodsig2_1(i)  =' '    !binary test resp mod sig values
         zmodsig2_2(i)  =' '
         zmodsig2_3(i)  =' '
         zmodsig2_4(i)  =' '
         zmodsig2_5(i)  =' '
         zmodsig2_6(i)  =' '
         coef(i)=' '            !coef of var needed for cumsum
         coefnum(i)=0.0         !real num version of coef
         ETA2_1(i)  =' '        !
         ETA2_2(i)  =' '
         ETA2_3(i)  =' '
         ETA2_4(i)  =' '
         ETA2_5(i)  =' '
         ETA2_6(i)  =' '
      end do
      do i = 1, MAX_PERTURB
         per_results (i)= ' '
         per_results_REF(i) = ' '
         per_results_TAR(i) = ' '
      end do
      call clear_pert_fields(carotidCO2_x,vertCO2_x,
     +     hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,
     +     dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +     exp_reflex_x,
     +     sw_x,SLNsw_x,carotidCO2,vertCO2,hypercap_5_O2,
     +     hypercap_5_air,hypercap_tbd,hypoxia_5,
     +     hypoxia_12,gasp,
     +     lobel,aorta_cath,pe,carotid_occ,
     +     nitrop,vc_cath,dopamine,vagus,
     +     capsaicin,pinch,deep_pain,
     +     codeine,nalox,methyserg,mucimol,
     +     dexameth,
     +     noinf,hyperinf,hypervent,pertnew1,
     +     carotidCO2_2,cgh,lcgh,SLNcgh,
     +     exp_reflex,
     +     sw,SLNsw)

      perturbation=0
      do i = 1, MAX_PERTURB
         perturb(i)=0
      end do
      perturb_applied='n'       ! default --> no perturbations applied
*
      AA_cord=' '               !initialize AA data fields
      AA_RLN=' '
      AA_vagus=' '
      AA_cVRG=' '
      AA_rVRG=' '
      AA_rtPRG=' '
      AA_ltPRG=' '
      AA_new1=' '
      AA_new2=' '
      AA_new3=' '
      do i = 1, MAX_AA
         AA(i)='n'
      end do
      AA_applied='n'            ! default --> no AA tests applied
*
      STA_phrenic=' '           !initialize STA data fields
      STA_RLN=' '
      STA_ELN=' '
      STA_lumbar=' '
      STA_cent_vagus=' '
      STA_splanch=' '
      STA_cerv_symp=' '
      STA_new1=' '
      STA_new2=' '
      STA_new3=' '
      do i = 1, MAX_STA
         STA(i)='n'
      end do
      STA_applied='n'           ! default --> no STA tests applied

      qdtfilenames = ' '
      qdtfilenames_200 = ' '
      do i = 1, MAX_NUM_QDTS
         files (i)= ' '
      end do
      button_choice = ' '
      mode_original = ' '
      dirname = ' '
      show_pulses = 'y'         !default = show I and E pulses on resp CTHs
*
*       *****  Load a text array with descriptions of each perturbation: *****

      per_text_abbrev(1)='carotidCO2'
      per_text_abbrev(2)='vert.CO2'
      per_text_abbrev(3)='5% CO2/O2'
      per_text_abbrev(4)='5% CO2/air'
      per_text_abbrev(5)='new hypercap'
      per_text_abbrev(6)='12% O2 in N2'
      per_text_abbrev(7)='5% O2 in N2'
      per_text_abbrev(8)='r 12% O2/N2'
      per_text_abbrev(9)='caroCO2hyper'
      per_text_abbrev(10)='r carotidCO2'
      per_text_abbrev(11)='lobelline'
      per_text_abbrev(12)='aort.cath'
      per_text_abbrev(13)='pe'
      per_text_abbrev(14)='carotid occ.'
      per_text_abbrev(15)='nitropruss'
      per_text_abbrev(16)='v.c. cath'
      per_text_abbrev(17)='nitropruss'
      per_text_abbrev(18)='vagotomy'
      per_text_abbrev(19)='capsaicin'
      per_text_abbrev(20)='pinch'
      per_text_abbrev(21)='deep pain'
      per_text_abbrev(22)='codeine'
      per_text_abbrev(23)='naloxone'
      per_text_abbrev(24)='methyserg.'
      per_text_abbrev(25)='muscimol'
      per_text_abbrev(26)='dexameth.'
      per_text_abbrev(27)='no infl.'
      per_text_abbrev(28)='hyperinfl.'
      per_text_abbrev(29)='hypervent.'
      per_text_abbrev(30)='trach. cgh.'
      per_text_abbrev(31)='laryn. cgh.'
      per_text_abbrev(32)='SLN cgh'
      per_text_abbrev(33)='exp. reflex'
      per_text_abbrev(34)='mech. sw.'
      per_text_abbrev(35)='SLN sw.'

*       ***** load a text array with descriptions of AA possibilities:  *****
      AA_text(1)='1. spinal cord'
      AA_text(2)='2. recurrent laryngeal n.'
      AA_text(3)='3. vagus n.'
      AA_text(4)='4. microstimulation of caudal VRG'
      AA_text(5)='5. microstimulation of rostral VRG'
      AA_text(6)='6. microstimulation of right PRG'
      AA_text(7)='7. microstimulation of left PRG'
      AA_text(8)=' '
      AA_text(9)=' '
      AA_text(10)=' '
      AA_text1(1)='cord'
      AA_text1(2)='RLN'
      AA_text1(3)='vagus'
      AA_text1(4)='micro cVRG'
      AA_text1(5)='micro rVRG'
      AA_text1(6)='micro rtPRG'
      AA_text1(7)='micro ltPRG'
      AA_text1(8)=' '
      AA_text1(9)=' '
      AA_text1(10)=' '

*       ***** load a text array with descriptions of STA possibilities:  *****
      STA_text(1)='1. phrenic n.'
      STA_text(2)='2. recurrent laryngeal n.'
      STA_text(3)='3. central vagus n.'
      STA_text(4)='4. lumbar n.'
      STA_text(5)='5. cervical sympathetic n.'
      STA_text(6)='6. expiratory laryngeal n.'
      STA_text(7)='7. splanchnic n.'
      STA_text(8)=' '
      STA_text(9)=' '
      STA_text(10)=' '
      STA_text1(1)='phrenic'
      STA_text1(2)='RLN'
      STA_text1(3)='cent.vagus'
      STA_text1(4)='lumbar n.'
      STA_text1(5)='cerv.symp.'
      STA_text1(6)='ELN'
      STA_text1(7)='splanchnic'
      STA_text1(8)=' '
      STA_text1(9)=' '
      STA_text1(10)=' '
      do i = 1, 115
         STA_res_text(i)=' '
      end do
      STA_res_text(48)='nothing' ! 0             48
      STA_res_text(49)='sh pk rt' ! 1             49
      STA_res_text(77)='sh pk rt- motor' ! M             77
      STA_res_text(109)='sh pk rt->motor' ! m             109
      STA_res_text(80)='sh pk rt->pre' ! P             80
      STA_res_text(112)='sh pk rt->pre' ! p             112
      STA_res_text(50)='br pk rt' ! 2             50
      STA_res_text(51)='pk ctr' ! 3             51
      STA_res_text(76)='pk left' ! L             76   
      STA_res_text(108)='pk left' ! l             108
      STA_res_text(68)='pk ctr/dip' ! D             68
      STA_res_text(100)='pk ctr/dip' ! d             100
      STA_res_text(52)='tr rt'  ! 4             52
      STA_res_text(53)='tr ctr' ! 5             53
      STA_res_text(84)='tr left' ! T             84
      STA_res_text(116)='tr left' ! t             116
      STA_res_text(54)='MP&T'   ! 6             54
      STA_res_text(55)='other'  ! 7             55
      STA_res_text(56)='nt'     ! 8             56
      STA_res_text(57)='EKG'    ! 9             57

      call concat(DB_FILES,'_mean_E.sav',MEAN__E,l)
      inquire(FILE=MEAN__E,EXIST=exist)

      if(exist.eqv..TRUE.)then
         OPEN (27,FILE=MEAN__E,STATUS='OLD',
     +        FORM='UNFORMATTED')
         read (27) mean_E
         close (27)
      else
         mean_E = 0.0
      end if
*
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                       *
*       *       Construct filenames for output files.           *
*       *       Allow user to overwrite existing files          *
*       *       if desired.                                     *
*       *                                                       *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*               ***** open data output files to hold histograms & text *****
*
*
      do i = 1, MAX_NUM_QDTS
         qdt_sav(i) = ' '
      end do
c      print *,'total_num_qdts = ',total_num_qdts
      do i = 1, total_num_qdts
         call concat(qdt_files(i),'.sav',qdt_sav(i),l)
c         print *,'qdt_sav(',i,') = ', qdt_sav(i)
         call read_QDTSAV(file_version,qdt_sav(i),ios,
     +        date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +        IDs,excluded_from_QDT,included_in_QDT,total_num_cells,
     +        I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +        BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +        min_bw_1(i),BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +        end_time,icycles,
     +        total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +        ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +        zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +        zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +        zmodsig2_6,coef,coefnum,card_type,exp_name,
     +        DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +        num_rej_ISI,num_rej_rISI,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +        num_acc_cycles,sp_per_cycle,ITAL)
        
      end do
      qdt_file_version = file_version
c      print *,'qdt_file_version: ',qdt_file_version
      min_bw = min_bw_1(1)
      qdt_for_ACHs = 1
      do i = 2, total_num_qdts  !if more than one qdt file involved, figure out
         if(min_bw_1(i).lt.min_bw)qdt_for_ACHs = i ! which qdt file has the ACHs at the smallest binwidth and
      end do                    ! use that one to view CTHs
      QDT_FILENAME = qdt_files(qdt_for_ACHs)
      current_qdt = qdt_for_ACHs ! show the data from the file with the lowest binwidths first
      QDTSAV = qdt_sav(qdt_for_ACHs)
      call remove_all_blanks(QDTSAV,LEN(QDTSAV))
*
*       ***** check to be sure that the version of the program is correct: *****
*
      OPEN (14,FILE=QDTSAV,FORM='FORMATTED',
     +     STATUS='OLD')
      rewind 14
      read (14,'(A2)')QDT_version
      close (unit=14)
      call sc_qdt (current_qdt)
      if(c_format.eq.'(109I6)')then !version 2 or lower
         OPEN (2,FILE=QDT_FILENAME,FORM='FORMATTED',
     +        STATUS='OLD',ACCESS='DIRECT',RECL=654)
      else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
         OPEN (2,FILE=QDT_FILENAME,FORM='FORMATTED',
     +        STATUS='OLD',ACCESS='DIRECT',RECL=1082)
      else
      end if
*
*
      call concat(DB_FILES,'_UNITS.txt',DBG,l)

      call concat(DB_FILES,'_UNITS_AA.txt',DBG_AA,l) !database info for indiv cells (antidromic activation data)
      call concat(DB_FILES,'_UNITS_STA.txt',DBG_STA,l) !database info for indiv cells (spike-triggered avg data)
      call concat(DB_FILES,'.pre_dbg',pre_DBG,l) !precursor file to final database file for cells
      call concat(pre_DBG,'.bak',pre_DBG_BAK,l)
      call concat(DB_FILES,'_RESPONSES_xass.txt',DBR_xass,l) !database info for perturbation/stimulus responses of indiv cells
      call concat(DB_FILES,'_RESPONSES_tsfs.txt',DBR_tsfs,l) !database info for perturbation/stimulus responses of indiv cells
      call concat(DB_FILES,'_CROSSES.txt',DBP,l) !database info for pairs of cells
      call concat(DB_FILES,'.pre_dbp',pre_DBP,l) !precursor file to final database file for pairs
      call concat(pre_DBP,'.bak',pre_DBP_BAK,l)
      call concat(DB_FILES,'.db.sav',DBSAV,l) !database "gamesave" file
      call concat(DBSAV,'.bak',DBSAV_BAK,l)
      call concat(DB_FILES,'_extra_cell_preDBG',pre_DBG_EXTRA,l) !precursor file for indiv EXTRA cells


*       ********** open the database-related files: **********

      global_mode = mode
      OPEN_FILES: if((mode.eq.'ed').or.(mode.eq.'cr').or.
     +     (mode.eq.'jl'))then  !open files for read/write
         if(mode.eq.'ed')then
c              isys=SYSTEM('pwd'//char(0))
            isys = SYSTEM ('cp '//pre_DBG//' '//pre_DBG_BAK
     +           //char(0))     !make backup copies of all 
            isys = 0
            isys = SYSTEM ('savelog -q -c 1000 -r backups '
     +           //pre_DBG_BAK//char(0))+isys
            isys = SYSTEM ('cp '//pre_DBP//' '//pre_DBP_BAK
     +           //char(0))+isys !  database-related files
            isys = SYSTEM ('savelog -q -c 1000 -r backups '
     +           //pre_DBP_BAK//char(0))+isys
            isys = SYSTEM ('cp '//DBSAV//' '//DBSAV_BAK
     +           //char(0))+isys
            isys = SYSTEM ('savelog -q -c 1000 -r backups '
     +           //DBSAV_BAK//char(0))+isys
            if(isys.ne.0.and.mode.eq.'ed')then
 2139          print '(//,''*** FATAL ERROR ***''
     +              //,''Backup operation failed. ''
     +              ''Possible directory permissions issue. ''
     +              /,''Recommend abort.  Abort? (y/n)  >> '',$)'
               read (*,'(A)') OK
               call lower_case(OK,LEN(OK))
               if(OK.eq.'y')then
                  stop 'Xanalysis aborted.' !stop the program
               elseif(OK.eq.'n')then
 2140             print '(//''Be aware that data may not be saved''
     +                 '' appropriately if you continue.''
     +                 ''Still want to proceed? (y/n)  >> '',$)'
                  read (*,'(A)') OK
                  call lower_case(OK,LEN(OK))
                  if(OK.eq.'y')then
                  elseif(OK.eq.'n')then
                     stop 'Xanalysis aborted.' !stop the program
                  else
                     goto 2140
                  end if
               else
                  goto 2139
               end if
            end if
            OPEN (15,FILE=DBSAV,FORM='FORMATTED',
     +           STATUS='OLD',ACCESS='SEQUENTIAL')
            rewind 15
            char2 = ' '
            read (15,'(A2)')char2
            close (15)
         end if
         
         OPEN (9,FILE=pre_DBG,FORM='FORMATTED',
     +        ACCESS='DIRECT',RECL=DBG_RECL) !+ (6*char*9) + (29*char*2) = 112 for more perturbation responses
         OPEN (10,FILE=pre_DBP,FORM='FORMATTED',
     +        ACCESS='DIRECT',RECL=iDBP_RECL)
      endif OPEN_FILES

*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       *                                                                         *
*       *  Open the *.QDTSAV file and read in the information entered             * 
*       *    previously when the histograms were generated                        *
*       *                                                                         *
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
      call read_QDTSAV(file_version,QDTSAV,ios,
     +     date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded_from_QDT,included_in_QDT,total_num_cells,
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


*     ***** IF MODE = 'cr':  *****
*     ***** We want every included cell to have a space in the preDBG file.  So ... we have to write a 'blank'
*     *****  record to the preDBG, containing only the experiment date, the recording number, the name of the 
*     *****  spike times file (BDT/EDT), and the IDcode of the cell - a kind of place-holder record. 

      if(mode.eq.'cr')then
         do cell = 1, MAX_NUM_CODES
            if(included_in_QDT(cell).eq.0)cycle
            histogram_number=(included_in_QDT(cell)*3)+1 !calculate the histogram number of resp CTH
            call write_preDBG(9,included_in_QDT(cell),ios,
     +           date_exp,BDT_FILE, !write the cell's info to the preDBG file
     +           recording,histogram_number,CELL_NAMES(IDs(cell)),
     +           cell,resp_type(IDs(cell)),aAP,aRL,adep,sort,
     +           comm,
     +           card(IDs(cell)),AA_cord,                          
     +           AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +           AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,
     +           AA_new3,STA_phrenic,STA_RLN,
     +           STA_cent_vagus,STA_lumbar,STA_cerv_symp,
     +           STA_ELN,STA_splanch,STA_new1,STA_new2,
     +           STA_new3,carotidCO2_x,vertCO2_x,
     +           hypercap_5_O2_x,hypercap_5_air_x,
     +           hypercap_tbd_x,hypoxia_5_x,hypoxia_12_x,
     +           gasp_x,lobel_x,aorta_cath_x,pe_x,
     +           carotid_occ_x,nitrop_x,vc_cath_x,
     +           dopamine_x,vagus_x,capsaicin_x,pinch_x,
     +           deep_pain_x,codeine_x,nalox_x,
     +           methyserg_x,mucimol_x,dexameth_x,
     +           noinf_x,hyperinf_x,hypervent_x,
     +           pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +           SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +           carotidCO2,vertCO2,
     +           hypercap_5_O2,hypercap_5_air,
     +           hypercap_tbd,hypoxia_12,
     +           hypoxia_5,gasp,lobel,
     +           aorta_cath,pe,carotid_occ,
     +           nitrop,vc_cath,dopamine,
     +           vagus,capsaicin,pinch,
     +           deep_pain,codeine,nalox,
     +           methyserg,mucimol,dexameth,
     +           noinf,hyperinf,hypervent,
     +           pertnew1,carotidCO2_2,cgh,
     +           lcgh,SLNcgh,exp_reflex,
     +           sw,SLNsw)
         end do

*     ***** There will be no data in the preDBP if mode = 'cr'.  We want every pair to have a space in the preDBP file.   *****
*     ***** So ... we have to write a 'blank' record to the preDBP, containing only the experiment date, the recording    *****
*     ***** number, the name of the spike times file (BDT), and the IDcodes of the cells - a kind of place-holder record. *****
         
         do REF = 1, MAX_NUM_CODES !GET_REF_cr
            if(included_in_QDT(REF).eq.0)cycle !GET_REF_cr !this cell not valid - get another one
            do TAR = REF+1, MAX_NUM_CODES !GET_TAR_cr
               if(included_in_QDT(TAR).eq.0)cycle !GET_TAR_cr !this cell not valid - get another one
               Q_pos = 0
               if(included_in_QDT(REF).gt.1)then
                  do i2 = 1,included_in_QDT(REF)-1                    
                     Q_pos = Q_pos + (total_num_cells-i2)             
                  end do
               end if
               Q_pos=Q_pos+(included_in_QDT(TAR)-included_in_QDT(REF))
               rel_loc = (Q_pos*4)-3  
               REC_num_x = 3 + (7*total_num_cells) + rel_loc !calculate the # of the CCH at first binwidth=hist.num.
               call write_preDBP(10,Q_pos,ios,recording,date_exp,
     +              BDT_FILE,CELL_NAMES(IDs(REF)),REF,
     +              CELL_NAMES(IDs(TAR)),TAR,
     +              prim,loc1,sec,loc2,rescom,czk,cprobk,
     +              cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +              REC_num_x)

            end do              !GET_TAR_cr
         end do                 !GET_REF_cr
      end if

*     ***** Check to be sure that the chosen QDT files are compatible (same experiment and recording,   *****
*     *****  same # of cells, same IDcodes, etc.)  Warn the user if not compatible and do not allow     *****
*     *****  the user to continue -- go back to the main menu.                                          *****

      min_bw_1(1) = BINW_1      !want to know which set of QDTs begins with
*                                                               !with the smallest BW so know which ACHs to display
*                                                               !when looking at CTHs
      if(total_num_qdts.gt.1)then
         if(mode.eq.'cr')then
            print '(//,T2,''Please wait while I sort out a few ''
     +           ''things ...'')'
         else
            print '(//,T2,''Please wait while I check to be ''
     +           ''sure that these QDT files are ''
     +           ''compatible ...'')'
         end if
         do i = 2, total_num_qdts
            warning = 'n'
            warning_date = 'n'
            warning_rec = 'n'
            warning_num_cells = 'n'
            warning_ID = 'n'
            warning_exc = 'n'
            warning_inc = 'n'
            date_exp_new = ' '
            recording_new = ' '
            do j = 1, MAX_NUM_CODES
               IDs_new(j) = 0
               included_in_QDT_new(j) = 0
               excluded_from_QDT_new(j) = 0
            end do
            call read_QDTSAV(file_version,qdt_sav(i),ios,
     +           date_exp_new,recording_new,protocol,BDT_FILE,
     +           QDT_FILENAME,IDs_new,excluded_from_QDT_new,
     +           included_in_QDT_new,total_num_cells_new,
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
     +           sd_rISI,num_rej_ISI,num_rej_rISI,c_MAX_INT,
     +           TOTAL_NUM_SHIFTS,num_acc_cycles,sp_per_cycle,ITAL)

            min_bw_1(i) = BINW_1
            if(date_exp.ne.date_exp_new)then
               warning='y'
               warning_date='y'
            end if
            if(recording.ne.recording_new)then
               warning='y'
               warning_rec='y'
            end if
            if(total_num_cells.ne.total_num_cells_new)then
               warning='y'
               warning_num_cells='y'
            end if
            do j = 1, MAX_NUM_CODES
               if(IDs(j).ne.IDs_new(j))then
                  warning='y'
                  warning_ID='y'
               end if
               if(excluded_from_QDT(j).ne.
     +              excluded_from_QDT_new(j))then
                  warning='y'
                  warning_exc='y'
               end if
               if(included_in_QDT(j).ne.included_in_QDT_new(j))then
                  warning='y'
                  warning_inc='y'
               end if
               if((warning_ID.eq.'y').and.(warning_exc.eq.'y') !no need to look at all if have already
     +              .and.(warning_inc.eq.'y'))exit ! found differences in IDcodes,included,
            end do              !  and excluded cells
            if(warning.eq.'y')exit
         end do

         
         if(warning.eq.'y')then
            print '(//,T2,''*** WARNING! ***  The QDT files are ''
     +           ''not compatible in the following ''
     +           ''respects:''/)'
            if(warning_date.eq.'y')
     +           print '(T10,''-- date of experiment: '',A11,
     +           '' vs. '',A11)',date_exp,date_exp_new
            if(warning_rec.eq.'y')
     +           print '(T10,''-- recording #: '',A2,'' vs. '',
     +           A2)',recording,recording_new
            if(warning_num_cells.eq.'y')
     +           print '(T10,''-- total number of cells used to ''
     +           ''calculate the QDT file: '',I10,
     +           '' vs. '',I10)',total_num_cells,total_num_cells_new
            if(warning_ID.eq.'y')
     +           print '(T10,''-- IDcodes of cells in the original ''
     +           ''spike times file (BDT, EDT, etc.)'')'
            if(warning_inc.eq.'y')
     +           print '(T10,''-- cells included when calculating '' 
     +           ''the QDT file'')'
            if(warning_exc.eq.'y')
     +           print '(T10,''-- cells excluded when calculating '' 
     +           ''the QDT file)'')'
 1141       print '(//,T2,''You may not look at mismatched QDT ''
     +           ''files!  Press <cr> to continue  >> '',$)'
            read (*,'(A1)')OK
            call upper_case(OK,LEN(OK))
            if(OK.eq.' ')then
               close (unit=9)   !close: *.pre_dbg
               close (unit=10)  !       *.pre_dbp
               close (unit=2)   !       *.qdt
               return
            else
               goto 1141        !force an appropriate response
            end if
         else                   !warning='n', so QDT files match and we can proceed
            if(mode.eq.'cr')then
               print '(/,T10,''OK - I''''m ready now'',//)'
            else
               print '(/,T10,''OK - they''''re compatible!  ''
     +              ''Let''''s move on ...'',//)'
            end if
         end if
      end if


*
*       *****  figure out which cells and pairs are valid for analysis and which ones are to be ignored:  *****
*
      if(mode.eq.'cr')then      !initialize the analyzed_cells and analyzed_pairs arrays
         do i = 1,MAX_NUM_CODES
            if((included_in_QDT(i).ne.0))then !if this cell was included in the .qdt file but
               if(excluded_from_QDT(i).eq.1)then ! was excluded from CTH and CCH calculations,
                  analyzed_cells(IDs(i))=0 ! then tag this cell as DO NOT ANALYZE - INVALID CELL
                  do m = 1,MAX_NUM_CODES
                     if((included_in_QDT(m).ne.0).and.
     +                    (excluded_from_QDT(m).ne.1))then
                        analyzed_pairs(IDs(i),IDs(m))=0 ! and tag its possible pairs as DO NOT ANALYZE - INVALID PAIR
                        analyzed_pairs(IDs(m),IDs(i))=0
                     end if
                  end do
               else             ! if CTHs and CCHs exist in the qdt file for this cell, then
                  analyzed_cells(IDs(i))=1 ! tag it as "NOT YET ANALYZED"
                  do j = i+1,MAX_NUM_CODES
                     if((included_in_QDT(j).ne.0).and. !if a possible TARGET cell was included in the .qdt file
     +                    (excluded_from_QDT(j).ne.1))then ! AND has not been excluded THEN
                        analyzed_pairs(IDs(i),IDs(j))=1 ! tag this pair as "NOT YET ANALYZED"
                     end if
                  end do
               end if
            end if
         end do
      end if

      if(mode.ne.'vt')call enter_date(date_exp) !be sure that experiment date is in correct format

      if((exp_name.eq.' ').and.
     +     (mode.ne.'vt'))then
 1142    print '(/,''Enter name of experiment  >> '',$)'
         read (*,'(A175)',err=1142) exp_name
         if(exp_name.eq.' ')goto 1142 !force an appropriate response
      end if
*
*
*
*               ***** read the database "gamesave" file: *****
*
*
      if((mode.eq.'ed').or.(mode.eq.'jl'))then
         print *,'goto read DBSAV'
         call read_DBSAV(file_version,DBSAV,ios,included_in_DB,
     +        excluded_from_DB,
     +        analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)
         qdtfilenames = ' '

*       ***** compare the QDTSAV and DBSAV files to see if user wants to associate a "new" QDT file with *****
*       *****  a set of existing database files; if the new QDT file contains more or fewer cells than   *****
*       *****  the original QDT file, then will have to convert the preDBG and preDBP files in order to  *****
*       *****  avoid re-entering all that data!                                                          *****

         call compare_and_convert_if_add_or_delete_units
     +        (included_in_QDT,included_in_DB,excluded_from_QDT,
     +        excluded_from_DB,IDs,pre_DBG,
     +        pre_DBP,total_num_cells,analyzed_cells,
     +        analyzed_pairs,BDT_FILE,CELL_NAMES,resp_type,
     +        date_exp,recording,qdt_files,qdt_sav,DBSAV,esc,mode,
     +        total_num_qdts,dirname)
         if(esc.eq.'y')return   !do not force QDT and DB files to match cells
         
         do i = 1, MAX_NUM_CODES
            included(i) = included_in_QDT(i) !restore current values of included()
            excluded(i) = excluded_from_QDT(i)
         end do
         
         call write_DBSAV(pgm_version,DBSAV,ios,included,
     +        excluded,analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)

*     ***** Be sure that every cell and pair has a record in the preDBG or preDBP file: *****

         do cell = 1, MAX_NUM_CODES
            if(included(cell).eq.0)cycle
            histogram_number=(included(cell)*3)+1 !calculate the histogram number of resp CTH
            call read_preDBG(9,included(cell),ios,char11,char30, !read preDBG record - fill dummy variables so don't
     +           char2,integer,char4,integer,char15, ! take a chance on corrupting data
     +           char5,char5,char5,char1,char40,char1,
     +           char3,char3,char3,char3,char3,char3,char3,char3,
     +           char3,char3,
     +           char12,char12,char12,char12,char12,char12,char12,
     +           char12,char12,char12,
     +           char9,char9,char9,char9,char9,char9,char9,char9,
     +           char9,char9,char9,char9,char9,char9,char9,
     +           char9,char9,char9,char9,char9,char9,char9,
     +           char9,char9,char9,char9,char9,char9,char9,char9,
     +           char9,char9,char9,char9,char9,
     +           char2x,char2x,char2x,char2x,char2x,char2x,char2x,
     +           char2x,char2x,char2x,
     +           char2x,char2x,char2x,char2x,char2x,char2x,char2x,
     +           char2x,char2x,char2x,
     +           char2x,char2x,char2x,char2x,char2x,char2x,char2x,
     +           char2x,char2x,char2x,char2x,char2x,char2x,char2x,
     +           char2x,
     +           char9x,char3x,char12x)
            if(ios.gt.0)then
               call write_preDBG(9,included(cell),ios, !if cannot read the record, write a place-holder record
     +              date_exp,BDT_FILE, !write the cell's info to the preDBG file
     +              recording,histogram_number,
     +              CELL_NAMES(IDs(cell)),
     +              cell,resp_type(IDs(cell)),aAP,aRL,adep,sort,
     +              comm,
     +              card(IDs(cell)),AA_cord,                          
     +              AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +              AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,
     +              AA_new3,STA_phrenic,STA_RLN,
     +              STA_cent_vagus,STA_lumbar,STA_cerv_symp,
     +              STA_ELN,STA_splanch,STA_new1,STA_new2,
     +              STA_new3,carotidCO2_x,vertCO2_x,
     +              hypercap_5_O2_x,hypercap_5_air_x,
     +              hypercap_tbd_x,hypoxia_5_x,hypoxia_12_x,
     +              gasp_x,lobel_x,aorta_cath_x,pe_x,
     +              carotid_occ_x,nitrop_x,vc_cath_x,
     +              dopamine_x,vagus_x,capsaicin_x,pinch_x,
     +              deep_pain_x,codeine_x,nalox_x,
     +              methyserg_x,mucimol_x,dexameth_x,
     +              noinf_x,hyperinf_x,hypervent_x,
     +              pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +              SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +              carotidCO2,vertCO2,
     +              hypercap_5_O2,hypercap_5_air,
     +              hypercap_tbd,hypoxia_12,
     +              hypoxia_5,gasp,lobel,
     +              aorta_cath,pe,carotid_occ,
     +              nitrop,vc_cath,dopamine,
     +              vagus,capsaicin,pinch,
     +              deep_pain,codeine,nalox,
     +              methyserg,mucimol,dexameth,
     +              noinf,hyperinf,hypervent,
     +              pertnew1,carotidCO2_2,cgh,
     +              lcgh,SLNcgh,exp_reflex,
     +              sw,SLNsw)
            end if
         end do
         
         do REF = 1, MAX_NUM_CODES !GET_REF_ed
            if(included(REF).eq.0)cycle !GET_REF_ed !this cell not valid - get another one
            do TAR = REF+1, MAX_NUM_CODES !GET_TAR_ed
               if(included(TAR).eq.0)cycle !GET_TAR_ed !this cell not valid - get another one
               Q_pos = 0
               if(included(REF).gt.1)then
                  do i2 = 1,included(REF)-1
                     Q_pos = Q_pos + (total_num_cells-i2)
                  end do
               end if
               Q_pos=Q_pos+(included(TAR)-included(REF))
               rel_loc = (Q_pos*4)-3  
               REC_num_x = 3 + (7*total_num_cells) + rel_loc !calculate the # of the CCH at first binwidth=hist.num.

               call read_preDBP(10,Q_pos,ios,char2,char11,
     +              char30,char4,integer,char4,integer,
     +              char20,char15,char20,char15,char40,
     +              char8,char8,char8,char8,char8,char8,char40,
     +              char8,integer)

               if(ios.gt.0)call write_preDBP(10,Q_pos,ios,
     +              recording,date_exp,
     +              BDT_FILE,CELL_NAMES(IDs(REF)),REF,
     +              CELL_NAMES(IDs(TAR)),TAR,
     +              prim,loc1,sec,loc2,rescom,czk,cprobk,
     +              cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +              REC_num_x)

            end do              !GET_TAR_ed
         end do                 !GET_REF_ed

*       ***** OK - comparison and conversion completed - move ahead *****

*       ***** load names of current qdt files assoc w/ the database files into the "qdtfilenames" string: *****
         j = 0
         k = 1
         qdtfilenames=' '
         A: do i = 1, total_num_qdts
             call strlength(qdt_files(i),LEN(qdt_files(i)),l) 
               call strlength(qdtfilenames,LEN(qdtfilenames),k)
               if (k.eq.1)k = 0
               qdtfilenames(k+1:k+l+1) = qdt_files(i)(1:l)//';'
c               print *,'qdt_files(',i,') = ',qdt_files(i)
c               print *,'qdtfilenames = ',qdtfilenames
         end do A
      end if

      do i = 1, MAX_NUM_CODES
         included(i) = included_in_QDT(i) !restore current values of included()
         if((excluded_from_QDT(i).eq.1).or. !combine excluded codes
     +        (excluded_from_DB(i).eq.1))then
            excluded(i) = 1
            if(IDs(i).eq.0) cycle
            if((analyzed_cells(IDs(i)).eq.1).or.
     +           (analyzed_cells(IDs(i)).eq.2))
     +           analyzed_cells(IDs(i))=analyzed_cells(IDs(i))+2

            do j = 1, MAX_NUM_CODES
               if(IDs(j).ne.0)then
                  if((analyzed_pairs(IDs(i),IDs(j)).eq.1).or.
     +                 (analyzed_pairs(IDs(i),IDs(j)).eq.2))
     +                 analyzed_pairs(IDs(i),IDs(j))=
     +                 analyzed_pairs(IDs(i),IDs(j))+2
                  if((analyzed_pairs(IDs(j),IDs(i)).eq.1).or.
     +                 (analyzed_pairs(IDs(j),IDs(i)).eq.2))
     +                 analyzed_pairs(IDs(j),IDs(i))=
     +                 analyzed_pairs(IDs(j),IDs(i))+2
               end if
            end do
            
         end if
      end do


      if((mode.eq.'cr').or.(qdtfilenames.eq.' '))then !create the qdtfilenames string - holds names of all
                                !  QDT files associated with these database files
         j = 1
         do i = 1, MAX_NUM_QDTS
            if(qdt_files(i).ne.' ')then !add this qdt filename to qdtfilenames - separate with semicolon
               call strlength(qdt_files(i),LEN(qdt_files(i)),l)
               qdtfilenames(j:j+l+1) = qdt_files(i)(1:l)//';'
               j = INDEX(qdtfilenames,' ')
            end if
         end do
      end if

      if((mode.eq.'cr'))then
         info_prev_imported = 'n'
         per_prev_imported = 'n'
         call write_DBSAV(pgm_version,DBSAV,ios,included, excluded
     +        ,analyzed_cells,analyzed_pairs,CELL_NAMES,perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported, per_prev_imported,total_num_cells_db
     +        ,qdtfilenames)
      end if

*
      if((import_info.eq.'y').or.(import_per.eq.'y'))then
         if(import_info.eq.'y')then
            call get_info_IDcodes(INFO_FILENAME,IDs_info)
c            print *,'opening ',trim(PER_FILENAME),' in directory'
            isys = system ('echo -n opening `pwd`/'//char(0))
            print '(A)',trim(PER_FILENAME)
            if(import_per.eq.'y')then !check for "extra" cells (cells in xassist pert file
               OPEN (3,FILE=PER_FILENAME, ! but not in xanalysis qdt file)
     +              STATUS='OLD',ACCESS='SEQUENTIAL',
     +              FORM='UNFORMATTED')
               goto 2167
 2168          goto 12677
 2167          read (3,ERR=2168,IOSTAT=ios) date_import,rec_import,
     +              IDs_per,TOTAL_STIM,per_res,perturb
12677          if(ios.ne.0)then
                  print '(''ios = '',I12)', ios
                  goto 2170
 2169             print '(''cannot import per data - ''
     +                 ''array mismatch'')'
                  goto 2171 
 2170             rewind 3
                  read(3, ERR=2169,IOSTAT=ios) date_import,
     +                 rec_import,IDs_per,TOTAL_STIM,
     +                 ((per_res(i,j),i=1,120),j=1,MAX_PERT),
     +                 perturb
               end if
 2171          close (3)
               do i = 1, MAX_NUM_CODES
                  IDs_extra(i) = 0
                  if(IDs_per(i).ne.0.and.IDs(i).eq.0.and.
     +                 IDs_info(i).ne.0)then
                     IDs_extra(i) = IDs_per(i)
                  end if
               end do
            elseif(import_per.eq.'n')then
               inquire(FILE=pre_DBG_EXTRA,EXIST=EXIST_extra_predbg)
               if(EXIST_extra_predbg.eqv..TRUE.)then
                  open (29,FILE=pre_DBG_EXTRA,
     +                 FORM='FORMATTED',
     +                 ACCESS='DIRECT',
     +                 RECL=DBG_RECL)
                  do i = 1, MAX_NUM_CHAN
                     iID = 0    !initialize the variable
                     call read_preDBG(29,i,ios,char11,char30, 
     +                    char2,iID,char4,integer,char15,                  
     +                    char5,char5,char5,char1,char40,char1,
     +                    char3,char3,char3,char3,char3,char3,
     +                    char3,char3,char3,char3,
     +                    char12,char12,char12,char12,char12,
     +                    char12,char12,char12,char12,char12,
     +                    char9,char9,char9,char9,char9,char9,
     +                    char9,char9,
     +                    char9,char9,char9,char9,char9,char9,
     +                    char9,
     +                    char9,char9,char9,char9,char9,char9,
     +                    char9,
     +                    char9,char9,char9,char9,char9,char9,
     +                    char9,char9,
     +                    char9,char9,char9,char9,char9,
     +                    char2x,char2x,char2x,char2x,char2x,
     +                    char2x,char2x,
     +                    char2x,char2x,char2x,
     +                    char2x,char2x,char2x,char2x,char2x,
     +                    char2x,char2x,
     +                    char2x,char2x,char2x,
     +                    char2x,char2x,char2x,char2x,char2x,
     +                    char2x,char2x,
     +                    char2x,char2x,char2x,char2x,char2x,
     +                    char2x,char2x,
     +                    char2x,char9x,char3x,char12x)
                     if(ios.eq.0)then !read was successful
                        IDs_info(iID) = 1
                        if((included(iID).eq.0).and. !if this cell is in per file but not
     +                       (IDs_per(iID).ne.0))then ! in the qdt file, it is an "extra" cell
                           IDs_extra(iID) = iID
                        end if
                     end if
                  end do
                  close (29)
               end if
            end if
            print '(/,''Importing cell information from '',A)', 
     +           INFO_FILENAME
            call read_info_file(INFO_FILENAME,
     +           IDs,
     +           included,
     +           AA_import,STA_import,name_import,
     +           comm_import,date_import,
     +           rec_import,AP_import,RL_import,
     +           depth_import,expname_import,
     +           dchan,ref_chan,
     +           ref_chan_refs,NAME_refs,AP_refs,
     +           RL_refs,depth_refs,BDT_FILE,
     +           recording,DB_FILES,IDs_extra)
            if(INFO_FILENAME.eq.' ')import_info='n'

 1155       if(exp_name.ne.expname_import)then !compare the name of the experiment in the QDT and the info files
               call strlength(exp_name,LEN(exp_name),jj)
               call strlength(expname_import,LEN(expname_import),jk)
               print '(///,T5,''***** WARNING *****  '',
     +              ''NAMES OF EXPERIMENT DO NOT MATCH ***** '',
     +              '' ('',A,'' vs '',A,'')'',
     +              /,T5,''Proceed with cell coordinate/AA/STA ''
     +              ''data import? (y/n)  >> '',$)',exp_name(1:jj),
     +              expname_import(1:jk)
               OK = ' '
               read (*,'(A1)') OK
               if((OK.eq.'n').or.(OK.eq.'N'))then
                  import_info = 'n'
               else if ((OK.eq.'y').or.(OK.eq.'Y'))then
 1156             print '(/,T5,''Please enter the correct ''
     +                 ''experiment name  >> '',$)'
                  read (*,'(A175)',err=1156) exp_name
               else
                  goto 1155     !force an appropriate response
               end if
            end if


            call enter_date(date_import) !be sure that date_import is in the correct format
            
 1162       if(date_exp.ne.date_import)then !compare the dates of the QDT and the info files
               print '(/,T5,''***** WARNING *****  '',
     +              ''DATES OF EXPERIMENT DO NOT MATCH ***** '',
     +              '' ('',A,'' vs '',A,'')'',
     +              /,T5,''Proceed with cell coordinate/AA/STA ''
     +              ''data import? (y/n)  >> '',$)',
     +              date_exp,date_import
               OK = ' '
               read (*,'(A1)') OK
               if((OK.eq.'n').or.(OK.eq.'N'))then
                  import_info = 'n'
               else if ((OK.eq.'y').or.(OK.eq.'Y'))then
 1163             print '(/,T5,''Please enter the correct date ''
     +                 ''(dd-mmm-yyyy)  >> '',$)'
                  read (*,'(A11)',err=1163) date_exp
                  call enter_date(date_exp) !be sure that date_exp is in the correct format
               else
                  goto 1162     !force an appropriate response
               end if
            end if
 1165       if((recording.ne.rec_import).and.
     +           (import_info.eq.'y'))then !compare the recording numbers of the QDT and the info files
               if(rec_import(1:1).ne.' '.and.rec_import(2:2).eq.' ')then
                  rec_import(2:2) = rec_import(1:1)
                  rec_import(1:1) =' '
               end if
               if(recording.ne.rec_import)then
                  print '(/,T5,''***** WARNING *****  '',
     +              ''RECORDING NUMBERS DO NOT MATCH *****  '',
     +              '' ('',A,'' vs '',A,'')'',
     +              /,T5,''Proceed with cell coordinate/AA/STA ''
     +              ''data import? (y/n)  >> '',$)',
     +              recording,rec_import
                  OK = ' '
                  read (*,'(A1)') OK
                  if((OK.eq.'n').or.(OK.eq.'N'))then
                     import_info = 'n'
                  else if ((OK.eq.'y').or.(OK.eq.'Y'))then
 1166                print '(/,T5,''Please enter the correct ''
     +                 ''recording number  >> '',$)'
                     read (*,'(A2)',err=1166) recording
                  else
                     goto 1165  !force an appropriate response
                  end if
               end if
            end if
         end if
         
         if(import_per.eq.'y')then
            print '(/,''Importing perturbation data from '',A)', 
     +           PER_FILENAME
            OPEN (3,FILE=PER_FILENAME,
     +           STATUS='OLD',ACCESS='SEQUENTIAL',
     +           FORM='UNFORMATTED')
            goto 1167
 1168       goto 11677
 1167       read (3,ERR=1168,IOSTAT=ios) date_import,rec_import,
     +           IDs_per,TOTAL_STIM,per_res,perturb
11677       if(ios.ne.0)then
               print '(''ios = '',I12)', ios
               goto 1170
 1169          print '(''cannot import per data - array mismatch'')'
               goto 1171
 1170          rewind 3
               read(3, ERR=1169,IOSTAT=ios) date_import,
     +              rec_import,IDs_per,TOTAL_STIM,
     +              ((per_res(i,j),i=1,120),j=1,MAX_PERT),
     +              perturb
            end if
 1171       close (3)
            
            call enter_date(date_import) !be sure that date_import is in the correct format
 1172       if(date_exp.ne.date_import)then !compare the dates of the QDT and the per files
c                   call strlength(exp_name,LEN(exp_name),jj)
c                   call strlength(expname_import,LEN(expname_import),jk)
               print '(/,T5,''***** WARNING *****  '',
     +              ''DATES OF EXPERIMENT DO NOT MATCH ***** '',
     +              '' ('',A,'' vs '',A,'')'',
     +              /,T5,''Proceed with cell perturbation ''
     +              ''data import? ''
     +              ''(y/n)  >> '',$)',date_exp,date_import
               OK = ' '
               read (*,'(A1)') OK
               if((OK.eq.'n').or.(OK.eq.'N'))then
                  import_per = 'n'
               else if ((OK.eq.'y').or.(OK.eq.'Y'))then
 1173             print '(/,T5,''Please enter the correct date ''
     +                 ''(dd-mmm-yyyy)  >> '',$)'
                  read (*,'(A11)',err=1173) date_exp
                  call enter_date(date_exp) !be sure that date_exp is in the correct format
               else
                  goto 1172     !force an appropriate response
               end if
            end if
 1175       if((recording.ne.rec_import).and.
     +           (import_per.eq.'y'))then !compare the recording numbers of the QDT and the per files
               if(rec_import(1:1).ne.' '.and.rec_import(2:2).eq.' ')then
                  rec_import(2:2) = rec_import(1:1)
                  rec_import(1:1) =' '
               end if
               if(recording.ne.rec_import)then
                  print '(/,T5,''***** WARNING *****  '',
     +              ''RECORDING NUMBERS DO NOT MATCH *****  '',
     +              '' ('',A,'' vs '',A,'')'',
     +              /,T5,''Proceed with cell perturbation ''
     +              ''data import? (y/n)  >> '',$)',
     +              recording,rec_import
                  OK = ' '
                  read (*,'(A1)') OK
                  if((OK.eq.'n').or.(OK.eq.'N'))then
                     import_per = 'n'
                  else if ((OK.eq.'y').or.(OK.eq.'Y'))then
 1176                print '(/,T5,''Please enter the correct ''
     +                 ''recording number  >> '',$)'
                     read (*,'(A2)',err=1176) recording
                  else
                     goto 1175  !force an appropriate response
                  end if
               end if
            end if
         end if
         
         call write_ref_electrode_data(BDT_FILE,DB_FILES,
     +        date_exp,recording,exp_name,
     +        NAME_refs,AP_refs,RL_refs,depth_refs,
     +        dchan,ref_chan,ref_chan_refs)
         do i = 1, MAX_NUM_CODES
            inquire(FILE=pre_DBG_EXTRA,
     +           EXIST=EXIST_extra_predbg)
c                print '(''IDs_per('',I3,'') = '',I)',i,IDs_per(i)
c                print '(''EXIST = '',L)',EXIST_extra_predbg
c                if(IDs_per(i).ne.0.and.IDs_info(i).ne.0)read (*,'(A)')
            extra = ' '
            iextra = 0
            if(included(i).ne.0)extra = 'n'
            if((included(i).eq.0).and. !if the cell does not exist in the QDT file
     +           ((EXIST_extra_predbg.eqv..TRUE.).and. ! but is in the PER file (from xassist), then
     +           IDs_per(i).ne.0))then ! look for it in the EXTRAS preDBG file
               do nn = 1, MAX_NUM_CHAN
                  open (29,FILE=pre_DBG_EXTRA,
     +                 FORM='FORMATTED',
     +                 ACCESS='DIRECT',
     +                 RECL=DBG_RECL)
                  call read_preDBG(29,nn,ios,date_exp,
     +                 BDT_FILE_30,recording_x,histogram_number_x,
     +                 names_temp,icellx,resptype,
     +                 AP_temp,RL_temp,
     +                 depth_temp,sort_temp,
     +                 comm_temp,cardd,AA_cord,AA_RLN,
     +                 AA_vagus,AA_cVRG,AA_rVRG,AA_rtPRG,AA_ltPRG,
     +                 AA_new1,AA_new2,AA_new3,STA_phrenic,STA_RLN,
     +                 STA_cent_vagus,STA_lumbar,STA_cerv_symp,
     +                 STA_ELN, STA_splanch,STA_new1,STA_new2,
     +                 STA_new3,carotidCO2_x,vertCO2_x,
     +                 hypercap_5_O2_x,
     +                 hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +                 hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +                 carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,
     +                 vagus_x,capsaicin_x,pinch_x,deep_pain_x,
     +                 codeine_x,nalox_x,methyserg_x,mucimol_x,
     +                 dexameth_x,noinf_x,hyperinf_x,hypervent_x,
     +                 pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +                 SLNcgh_x,
     +                 exp_reflex_x,sw_x,SLNsw_x,carotidCO2,vertCO2,
     +                 hypercap_5_O2,hypercap_5_air,hypercap_tbd,
     +                 hypoxia_12,hypoxia_5,gasp,lobel,aorta_cath,pe,
     +                 carotid_occ,nitrop,vc_cath,dopamine,vagus,
     +                 capsaicin, pinch,deep_pain, codeine,nalox,
     +                 methyserg,mucimol,dexameth,noinf,hyperinf,
     +                 hypervent,pertnew1,carotidCO2_2,cgh,lcgh,
     +                 SLNcgh,exp_reflex,sw,SLNsw,per_results,
     +                 AA_results,STA_results)
                  close (29)
                  if(ios.ne.0)exit
c                      print '(''icellx = '',I)',icellx
c                      read (*,'(A)')
                  if(icellx.eq.i)then !this is an "extra" cell
c                         print '(''found extra cell = '',I)',i
c                         read (*,'(A)')
                     extra = 'y'
                     iextra = nn
                     exit
                  end if
               end do
            end if
            if(included(i).ne.0)then !info for this cell already exists - read it
               icell=i
               histogram_number=(included(icell)*3)+1 !histogram number of resp CTH
               call read_preDBG(9,included(icell),ios,date_exp,
     +              BDT_FILE_30,recording_x,histogram_number_x,
     +              names_temp,icellx,resptype,
     +              AP_temp,RL_temp,
     +              depth_temp,sort_temp,
     +              comm_temp,cardd,AA_cord,AA_RLN,
     +              AA_vagus,AA_cVRG,AA_rVRG,AA_rtPRG,AA_ltPRG,
     +              AA_new1,AA_new2,AA_new3,STA_phrenic,STA_RLN,
     +              STA_cent_vagus,STA_lumbar,STA_cerv_symp,
     +              STA_ELN, STA_splanch,STA_new1,STA_new2,
     +              STA_new3,carotidCO2_x,vertCO2_x,
     +              hypercap_5_O2_x,
     +              hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +              hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +              carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,
     +              vagus_x,capsaicin_x,pinch_x,deep_pain_x,
     +              codeine_x,nalox_x,methyserg_x,mucimol_x,
     +              dexameth_x,noinf_x,hyperinf_x,hypervent_x,
     +              pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +              SLNcgh_x,
     +              exp_reflex_x,sw_x,SLNsw_x,carotidCO2,vertCO2,
     +              hypercap_5_O2,hypercap_5_air,hypercap_tbd,
     +              hypoxia_12,hypoxia_5,gasp,lobel,aorta_cath,pe,
     +              carotid_occ,nitrop,vc_cath,dopamine,vagus,
     +              capsaicin, pinch,deep_pain, codeine,nalox,
     +              methyserg,mucimol,dexameth,noinf,hyperinf,
     +              hypervent,pertnew1,carotidCO2_2,cgh,lcgh,
     +              SLNcgh,exp_reflex,sw,SLNsw,per_results,
     +              AA_results,STA_results)
            end if
            
*         ***** load the imported data: *****

            if(import_info.eq.'y'.and.extra.eq.'n')then
               info_prev_imported = 'y'
               names_temp=name_import(IDs(i))
               AP_temp=AP_import(IDs(i))
               RL_temp=RL_import(IDs(i))
               depth_temp=depth_import(IDs(i))
               do j = 1, MAX_AA
                  if(j.eq.1)AA_cord=AA_import(IDs(i),j)
                  if(j.eq.2)AA_RLN=AA_import(IDs(i),j)
                  if(j.eq.3)AA_vagus=AA_import(IDs(i),j)
                  if(j.eq.4)AA_cVRG=AA_import(IDs(i),j)
                  if(j.eq.5)AA_rVRG=AA_import(IDs(i),j)
                  if(j.eq.6)AA_rtPRG=AA_import(IDs(i),j)
                  if(j.eq.7)AA_ltPRG=AA_import(IDs(i),j)
                  if(j.eq.8)AA_new1=AA_import(IDs(i),j)
                  if(j.eq.9)AA_new2=AA_import(IDs(i),j)
                  if(j.eq.10)AA_new3=AA_import(IDs(i),j)
                  if(AA_import(IDs(i),j).ne.' ')then
                     AA_applied='y'
                     AA(j) = 'y'
                  end if
               end do
               do j = 1, MAX_STA
                  if(j.eq.1)STA_phrenic=STA_import(IDs(i),j)
                  if(j.eq.2)STA_RLN=STA_import(IDs(i),j)
                  if(j.eq.3)STA_cent_vagus=
     +                 STA_import(IDs(i),j)
                  if(j.eq.4)STA_lumbar=STA_import(IDs(i),j)
                  if(j.eq.5)STA_cerv_symp=STA_import(IDs(i),j)
                  if(j.eq.6)STA_ELN=STA_import(IDs(i),j)
                  if(j.eq.7)STA_splanch=STA_import(IDs(i),j)
                  if(j.eq.8)STA_new1=STA_import(IDs(i),j)
                  if(j.eq.9)STA_new2=STA_import(IDs(i),j)
                  if(j.eq.10)STA_new3=STA_import(IDs(i),j)
                  if(STA_import(IDs(i),j).ne.' ')then
                     STA_applied='y'
                     STA(j) = 'y'
                  end if
               end do
               do j = 1, MAX_AA
                  if(AA_import(IDs(i),j).ne.' ')then
                     AA_applied='y'
                     AA(j) = 'y'
                  end if
               end do
               do j = 1, MAX_STA
                  if(STA_import(IDs(i),j).ne.' ')then
                     STA_applied='y'
                     STA(j) = 'y'
                  end if
               end do
            end if

            if(import_per.eq.'y')then !perturbation data was read-in before - wipe out obsolete data
               call clear_pert_fields(carotidCO2_x,vertCO2_x,
     +              hypercap_5_O2_x,
     +              hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +              hypoxia_12_x,gasp_x,
     +              lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +              nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +              capsaicin_x,pinch_x,deep_pain_x,
     +              codeine_x,nalox_x,methyserg_x,mucimol_x,
     +              dexameth_x,
     +              noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +              carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +              exp_reflex_x,
     +              sw_x,SLNsw_x,carotidCO2,vertCO2,hypercap_5_O2,
     +              hypercap_5_air,hypercap_tbd,hypoxia_5,
     +              hypoxia_12,gasp,
     +              lobel,aorta_cath,pe,carotid_occ,
     +              nitrop,vc_cath,dopamine,vagus,
     +              capsaicin,pinch,deep_pain,
     +              codeine,nalox,methyserg,mucimol,
     +              dexameth,
     +              noinf,hyperinf,hypervent,pertnew1,
     +              carotidCO2_2,cgh,lcgh,SLNcgh,
     +              exp_reflex,
     +              sw,SLNsw)
               per_prev_imported = 'y'


               if((IDs_per(i).ne.0.and.(included(i).ne.0)).or. !if the cell is in the QDT file ..or..
     +              (IDs_extra(i).ne.0))then ! the cell is a legitimate extra ...
                  do k = 1, MAX_PERT
                     do m = 1, MAX_PERTURB
                        if(perturb(m).eq.k)then
                           perturb_applied = 'y'
                           if(m.eq.1)then
                              carotidCO2_x=per_res(IDs_per(i),k)
                           elseif(m.eq.2)then
                              vertCO2_x=per_res(IDs_per(i),k)
                           elseif(m.eq.3)then
                              hypercap_5_O2_x=per_res(IDs_per(i),k)
                           elseif(m.eq.4)then
                              hypercap_5_air_x=per_res(IDs_per(i),k)
                           elseif(m.eq.5)then
                              hypercap_tbd_x=per_res(IDs_per(i),k)
                           elseif(m.eq.6)then
                              hypoxia_12_x=per_res(IDs_per(i),k)
                           elseif(m.eq.7)then
                              hypoxia_5_x=per_res(IDs_per(i),k)
                           elseif(m.eq.8)then
                              gasp_x=per_res(IDs_per(i),k)
                           elseif(m.eq.9)then
                              pertnew1_x=per_res(IDs_per(i),k)
                           elseif(m.eq.10)then
                              carotidCO2_2_x=per_res(IDs_per(i),k)
                           elseif(m.eq.11)then
                              lobel_x=per_res(IDs_per(i),k)
                           elseif(m.eq.12)then
                              aorta_cath_x=per_res(IDs_per(i),k)
                           elseif(m.eq.13)then
                              pe_x=per_res(IDs_per(i),k)
                           elseif(m.eq.14)then
                              carotid_occ_x=per_res(IDs_per(i),k)
                           elseif(m.eq.15)then
                              nitrop_x=per_res(IDs_per(i),k)
                           elseif(m.eq.16)then
                              vc_cath_x=per_res(IDs_per(i),k)
                           elseif(m.eq.17)then
                              dopamine_x=per_res(IDs_per(i),k)
                           elseif(m.eq.18)then
                              vagus_x=per_res(IDs_per(i),k)
                           elseif(m.eq.19)then
                              capsaicin_x=per_res(IDs_per(i),k)
                           elseif(m.eq.20)then
                              pinch_x=per_res(IDs_per(i),k)
                           elseif(m.eq.21)then
                              deep_pain_x=per_res(IDs_per(i),k)
                           elseif(m.eq.22)then
                              codeine_x=per_res(IDs_per(i),k)
                           elseif(m.eq.23)then
                              nalox_x=per_res(IDs_per(i),k)
                           elseif(m.eq.24)then
                              methyserg_x=per_res(IDs_per(i),k)
                           elseif(m.eq.25)then
                              mucimol_x=per_res(IDs_per(i),k)
                           elseif(m.eq.26)then
                              dexameth_x=per_res(IDs_per(i),k)
                           elseif(m.eq.27)then
                              noinf_x=per_res(IDs_per(i),k)
                           elseif(m.eq.28)then
                              hyperinf_x=per_res(IDs_per(i),k)
                           elseif(m.eq.29)then
                              hypervent_x=per_res(IDs_per(i),k)
                           elseif(m.eq.30)then
                              cgh_x=per_res(IDs_per(i),k)
                           elseif(m.eq.31)then
                              lcgh_x=per_res(IDs_per(i),k)
                           elseif(m.eq.32)then
                              SLNcgh_x=per_res(IDs_per(i),k)
                           elseif(m.eq.33)then
                              exp_reflex_x=per_res(IDs_per(i),k)
                           elseif(m.eq.34)then
                              sw_x=per_res(IDs_per(i),k)
                           elseif(m.eq.35)then
                              SLNsw_x=per_res(IDs_per(i),k)
                           end if
                        end if
                     end do
                  end do
               end if
            end if

            if(extra.ne.' ')then
               if(extra.eq.'n')then
                  iunit = 9
                  irecnum = included(icell)
                  ihistnum = ((included(icell)*3)+1)
               elseif(extra.eq.'y')then
                  iunit = 29
                  irecnum = iextra
                  ihistnum = 0
                  open (29,FILE=pre_DBG_EXTRA,
     +                 FORM='FORMATTED',
     +                 ACCESS='DIRECT',
     +                 RECL=DBG_RECL)
               else
                  print '(I4,''is neither fish nor fowl!'')',i
               endif
               
               call write_preDBG(iunit,irecnum,ios,
     +              date_exp,BDT_FILE,
     +              recording,ihistnum,
     +              names_temp,
     +              i,resptype,       
     +              AP_temp,RL_temp,
     +              depth_temp,
     +              sort_temp,comm_temp,
     +              cardd,AA_cord,
     +              AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +              AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,
     +              AA_new3,STA_phrenic,STA_RLN,
     +              STA_cent_vagus,STA_lumbar,
     +              STA_cerv_symp,
     +              STA_ELN,STA_splanch,STA_new1,STA_new2,
     +              STA_new3,carotidCO2_x,vertCO2_x,
     +              hypercap_5_O2_x,hypercap_5_air_x,
     +              hypercap_tbd_x,hypoxia_5_x,
     +              hypoxia_12_x,
     +              gasp_x,lobel_x,aorta_cath_x,pe_x,
     +              carotid_occ_x,nitrop_x,vc_cath_x,
     +              dopamine_x,vagus_x,capsaicin_x,
     +              pinch_x,
     +              deep_pain_x,codeine_x,nalox_x,
     +              methyserg_x,mucimol_x,dexameth_x,
     +              noinf_x,hyperinf_x,hypervent_x,
     +              pertnew1_x,carotidCO2_2_x,cgh_x,
     +              lcgh_x,
     +              SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +              carotidCO2,vertCO2,
     +              hypercap_5_O2,hypercap_5_air,
     +              hypercap_tbd,hypoxia_12,
     +              hypoxia_5,gasp,lobel,
     +              aorta_cath,pe,carotid_occ,
     +              nitrop,vc_cath,dopamine,
     +              vagus,capsaicin,pinch,
     +              deep_pain,codeine,nalox,
     +              methyserg,mucimol,dexameth,
     +              noinf,hyperinf,hypervent,
     +              pertnew1,carotidCO2_2,cgh,
     +              lcgh,SLNcgh,exp_reflex,
     +              sw,SLNsw)
               if(extra.eq.'y')close (29)
            end if
            
         end do
      end if

      if((mode.eq.'ed').or.(mode.eq.'cr'))then !re-write the DBSAV file now that have finished importing
         call write_DBSAV(pgm_version,DBSAV,ios,included,
     +        excluded,analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)
      end if

      if(mode.eq.'cr')RETURN    !return to the main menu if in 'create' mode

 1150 screen_CTHs = 'n'         !default values
      screen_CCHs = 'n'
      if(mode_original.ne.' ')then
         mode=mode_original     !restore the mode (may have changed from edit to create)
         mode_original = ' '
      end if
      if(mode.eq.'cr')
     +     print '(5(/),10x,''ENTER DATA for:'',/)'
      if(mode.eq.'ed')
     +     print '(5(/),10x,''ENTER or EDIT DATA for:'',/)'
      if(mode.eq.'vt')
     +     print '(5(/),10x,''VIEW:'',/)'
      if(mode.eq.'jl')
     +     print '(5(/),10x,''VIEW DATA for:'',/)'
      print'(15x,''C -- individual Cells'',/,
     +     15x,''P -- Pairs'',//,
     +     15x,''M -- return to main Menu  '',//)'
      print '(20x,''>> '',$)'
      read (*,'(A)',err=1150)OK
      call upper_case(OK,LEN(OK))
      if(OK.eq.'C')then         !set the flags
         screen_CTHs = 'y'
      else if(OK.eq.'P')then
         screen_CCHs = 'y'
      else if(OK.eq.'M')then
         if((mode.eq.'cr').or.(mode.eq.'ed'))then 
            call write_DBSAV(pgm_version,DBSAV,ios,included,
     +           excluded,      !be sure that DBSAV is current
     +           analyzed_cells,analyzed_pairs,CELL_NAMES,
     +           perturb_applied,
     +           perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +           info_prev_imported,
     +           per_prev_imported,total_num_cells_db,
     +           qdtfilenames)
         end if   
         if(mode.eq.'jl'.and.dirname.ne.' ')then
            call strlength(DBSAV,LEN(DBSAV),l_db)
            call strlength(dirname,LEN(dirname),l_d)
            inquire(FILE=dirname(1:l_d)//'/'//DBSAV(1:l_db),
     +           EXIST=exist10)
            if(exist10.eqv..TRUE.)then
               print '(''restore the data files'')'
               call strlength(pre_DBG,LEN(pre_DBG),l_dbg)
               call strlength(pre_DBP,LEN(pre_DBP),l_dbp)
               isys = SYSTEM('rm '//DBSAV(1:l_db)//char(0))
               isys = SYSTEM('rm '//pre_DBG(1:l_dbg)//char(0))
               isys = SYSTEM('rm '//pre_DBP(1:l_dbp)//char(0))

               icd = CHDIR(dirname(1:l_d)//char(0))
               isys=SYSTEM('mv '//DBSAV(1:l_db)//' ..'//char(0))
               isys=SYSTEM('mv '//pre_DBG(1:l_dbg)//' ..'
     +              //char(0))
               isys=SYSTEM('mv '//pre_DBP(1:l_dbp)//' ..'
     +              //char(0))
               icd = CHDIR('..'//char(0))

               print '(A)','rmdir '//dirname(1:l_d)
               isys=SYSTEM('rmdir '//dirname(1:l_d)//char(0))

            end if
         end if
         close (unit=9)         !close: *.pre_dbg
         close (unit=10)        !       *.pre_dbp
         close (unit=2)         !       *.qdt
         return
      else
         goto 1150              !C_P_M_menu            !force an appropriate answer
      end if

*
*               *****  examine analyzed_cells and analyzed_pairs to determine if  *****
*               *****   all eligible cells / pairs have been analyzed:  *****
*
      cells_completed = 'y'     !assume that all eligible cells / pairs have
      pairs_completed = 'y'     !  been analyzed UNLESS can find a cell / pair
                                !  flagged as "valid -- not yet analyzed"
                                !0 --> invalid / ineligible cell
                                !1 --> valid cell not yet analyzed
                                !2 --> analysis of cell is completed
      do i = 1, MAX_NUM_CHAN
         if(analyzed_cells(i).eq.1)cells_completed = 'n'
         do j = i+1,MAX_NUM_CHAN
            if(analyzed_pairs(i,j).eq.1)pairs_completed = 'n'
         end do
      end do


*
*
*       ********************************************************
*       ********************************************************
*       ********************************************************
*         
*       ***** begin showing CTHs:  *****
*
*               ***** read in phrenic CTH data  *****
*
      read (2,c_format,REC=1,err=191) itype,NHW,ICN,phrenic_hist
 191  BINW=(FLOAT(NHW)/50.)
      write (bwtext,'(f7.1)') BINW
      bwtext1=bwtext//' ms.'
      Rcode=ICN(1)              !ID code of E pulse
      Tcode=ICN(3)              !ID code of phrenic
      Revents=ICN(2)            !# of E pulses used
      Tevents=ICN(4)            !# of phrenic spikes used
*
*
*               ***** read normalized phrenic CTH  *****
*
      read (2,c_format,REC=2,err=192) itype,NHW,ICN,norm_phrenic_hist
 192  NORM_BINW=(FLOAT(NHW)/50.)
      write (bwtext,'(f7.1)') NORM_BINW
      bwtext2=bwtext//' ms.'
*
*
*               *****   read cardiac overlay histogram  *****
*
      read (2,c_format,REC=3,err=193) itype,NHW,ICN,cardiac_hist
 193  BINW2=(FLOAT(NHW)/50.)
      write (bwtext,'(f7.1)') BINW2
      bwtext3=bwtext//' ms.'

      if(screen_CTHs.eq.'y')goto 198  
      if(screen_CCHs.eq.'y')goto 298
*
*
*               ***** read in cell CTH data  (check for excluded cells) *****
*
*

 198  if((mode.eq.'ed').and.(cells_completed.eq.'n'))then
         print '(//,''In case you''''re wondering, ''
     +        ''CELL ANALYSIS IS INCOMPLETE.  '',//,T5,
     +        ''Would you like to see a list of unanalyzed cells?''
     +        '' (y/n)  >> '',$)'
         read (*,'(A)') OK
         call upper_case(OK,LEN(OK))
         if(OK.eq.'Y')then
            j = 0
            print '(//,''IDcodes of UNANALYZED CELLS:'',/)'
            do i = 1, MAX_NUM_CODES
               if(IDs(i).eq.0)cycle
               if((excluded(i).eq.0).and.
     +              (analyzed_cells(IDs(i)).eq.1))then       
                  print '(I3,''  '',$)',i ! print codes of unanalyzed cells one by one
                  j = j + 1
                  if(j.eq.20)then ! print 20 IDcodes per line
                     print '(/)'
                     j = 0
                  end if
               end if   
            end do
         else if(OK.eq.'N')then
         else
            goto 198            !force an appropriate response
         end if
 1198    print '(//,T5,''Would you like to enter data for ''
     +        ''the unanalyzed cells only? (y/n) >> '',$)'
         read (*,'(A1)') OK
         call upper_case(OK,LEN(OK))
         if(OK.eq.'Y')then
            print '(/,T5,''OK.  Just so you know -- ''
     +           ''You will be shown ''
     +           ''only the unanalyzed cells ''
     +           /,T30,''until you ''
     +           ''return to the CELLS/PAIRS menu.''
     +           //,T10,''Press <cr> to continue  >> '',$)'
            read (*,'(A)')
            mode_original = mode
            mode = 'cr'         !temporarily switch to "create" mode
         else if(OK.eq.'N')then !leave well enough alone
         else
            goto 1198           !force an appropriate response
         end if
      end if

*     ***** open the window for CTH display and set its parameters: *****

c$$$         integer height
c$$$         character*40 geometry
c$$$        character*6 c_width,c_height
c$$$         height=height1*0.5                             !prepare to open the CTH viewing window
c$$$     write (c_height,'(I6)') height
c$$$     geometry=' '
c$$$     geometry=c_width//'x'//c_height//'-700-5'
c$$$         call remove_all_blanks(geometry,LEN(geometry))     
c$$$         call strlength(geometry,LEN(geometry),l)
c$$$         call strlength(USER,LEN(USER),m)
c$$$             isys=SYSTEM ('xwcreate -wmdir /dev/screen/'//USER(1:m)//
c$$$     +          ' -title '//WINDOW1//
c$$$     +          ' -geometry ='//geometry(1:l)//char(0))
c$$$             fildes=gopen ('/dev/screen/'//USER(1:m)//'/'//
c$$$     +              WINDOW1//char(0),
c$$$     +      OUTINDEV,'sox11'// char(0),INIT)
      fildes=gopen (width,height1/2,-700,-5,'XAnalysis'//char(0))

      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call shade_mode(fildes,IOR(INIT,CMAP_NORMAL),0)
      call mapping_mode(fildes,1)
      call view_port(fildes,0.,0.,1.,1.)
      call view_window(fildes,1.,1.,1700.,800.)
      call text_font_index(fildes,6)
      call background_color(fildes,1.,1.,1.) !white background
      call text_color(fildes,0.,0.,0.)
      call line_color(fildes,0.,0.,0.)
      call perimeter_color(fildes,0.,0.,0.)
*
      call echo_type(fildes,0,1,.625,.1,0.0)
      call set_locator(fildes,1,.625,.1,0.0)
      call track(fildes,fildes,1)
*
      call character_width(fildes,0.010)
      call character_height(fildes,0.050)

      call background_color(fildes,1.,1.,1.) !white background
      call clear(fildes)        !display the new background color immediately
*

      id = 0
      cell = 0

*
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do cell = 1, MAX_NUM_CODES !find the first included cell
         if(IDs(cell).eq.0)then
            cycle
         else
         end if
         if((mode.eq.'cr').and. !if creating, only want to look at unanalyzed cells
     +        (analyzed_cells(IDs(cell)).eq.2))cycle
         if(((mode.eq.'cr').or.(mode.eq.'ed').or.(mode.eq.'jl'))
     +        .and.((analyzed_cells(IDs(cell)).eq.1).or.
     +        (analyzed_cells(IDs(cell)).eq.2)))goto 210 !a valid cell - use it!
         if((mode.eq.'vt').and.(included(cell).ne.0).and.
     +        (excluded(cell).ne.1))goto 210 !a valid cell - use it!
      end do               

*     ***** if have reached the next lines of code, then have "fallen through" the do loop to an undesired cell, *****
*     *****   usually (hopefully) the phrenic code, boundary code, or cardiac pulse code -- go back to previous menu    *****

      print '(///,T5,''Data has already been entered for ''
     +     ''all cells.  <cr> to return to previous menu  >> '',$)'
      read (*,'(A)')                                       
      call clear(fildes)
      close = gclose(fildes)
      goto 1150                 !C_P_M_menu                                            

*
 210  if(mode.eq.'cr')then
         CTHsaved='n'
      elseif(mode.eq.'ed')then
         if(excluded(cell).eq.1)then !excluded cells shown by direct access must be re-saved
            CTHsaved = 'n'
         else
            CTHsaved = 'y'
         end if
         if(analyzed_cells(IDs(cell)).ne.2)CTHsaved='n'
      else
         CTHsaved= ' '
      end if

      enter_info = ' '          !reset an editing flag (default='A' --> enter all info for cell)
      aAP = ' '
      aRL = ' '
      adep = ' '
      sort  =' '
      comm = ' '
      changed_comment = 'n'
      AA_cord=' '               !initialize AA data fields
      AA_RLN=' '
      AA_vagus=' '
      AA_cVRG=' '
      AA_rVRG=' '
      AA_rtPRG=' '
      AA_ltPRG=' '
      AA_new1=' '
      AA_new2=' '
      AA_new3=' '
      do i = 1, MAX_AA
         AA_results=' '
      end do
      STA_phrenic=' '           !initialize STA data fields
      STA_RLN=' '
      STA_ELN=' '
      STA_lumbar=' '
      STA_cent_vagus=' '
      STA_splanch=' '
      STA_cerv_symp=' '
      STA_new1=' '
      STA_new2=' '
      STA_new3=' '
      do i = 1, MAX_STA
         STA_results=' '
      end do
      call clear_pert_fields(carotidCO2_x,vertCO2_x,
     +     hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,
     +     dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +     exp_reflex_x,
     +     sw_x,SLNsw_x,carotidCO2,vertCO2,hypercap_5_O2,
     +     hypercap_5_air,hypercap_tbd,hypoxia_5,
     +     hypoxia_12,gasp,
     +     lobel,aorta_cath,pe,carotid_occ,
     +     nitrop,vc_cath,dopamine,vagus,
     +     capsaicin,pinch,deep_pain,
     +     codeine,nalox,methyserg,mucimol,
     +     dexameth,
     +     noinf,hyperinf,hypervent,pertnew1,
     +     carotidCO2_2,cgh,lcgh,SLNcgh,
     +     exp_reflex,
     +     sw,SLNsw)
      do i = 1, MAX_PERTURB
         per_results=' '
      end do

      call clear(fildes)        !clear the display screen
*
      REC_num=(included(cell)*3)+1 !calculate the location of the resp CTH 
*
      do i = 1, NUM_BINS
         CTH(i)=0
      end do
      read (2,c_format,REC=REC_num,err=211) itype,NHW,ICN,CTH
 211  if(itype.ne.2)print '(''NOT A CTH'')'

      Rcode=ICN(1)              !ID code of E pulse
      Tcode=ICN(3)              !ID code of target cell
      Revents=ICN(2)            !# of E pulses used
      Tevents=ICN(4)            !# of cell spikes used
      Revents_p = ICN(2)
*
      if(cell.ne.Tcode)then
         IDmismatch='y'
         print '(''IDcode mismatch!!'')'
         read (*,'(A)')
      else
         IDmismatch='n'
      end if

      do i = 1, NUM_BINS
         norm_CTH(i)=0
      end do
      read (2,c_format,REC=(REC_num+1),err=212) 
     +     itype,NHW,ICN,norm_CTH !read the norm. resp. CTH
      Revents_n = ICN(2)
 212  do i = 1, NUM_BINS
         cardCCH(i)=0
      end do
      read (2,c_format,REC=(REC_num+2),err=213) 
     +     itype,NHW,ICN,cardCCH !read the cardiac CTH
      Revents_c = ICN(2)
*
*               ***** read in the cell's ACH calculated at the smallest binwidth: *****
*
 213  REC_num=((included(cell)*4)+(3*total_num_cells)) !location of 1st ACH for this cell
      do i = 1, NUM_BINS
         TAR_ACH(i)=0
      end do
      
      read(2,c_format,REC=REC_num,err=214) 
     +     itype,NHW,ICN,TAR_ACH
 214  bw_ACH_1 = FLOAT(NHW)/50.
      write (bwtext,'(f7.1)') bw_ACH_1
      bwtext4=bwtext
      Revents_ach1 = ICN(2)
*
*               ***** read in the cell's ACH calculated at the next-to-greatest binwidth: *****
*
      REC_num=((included(cell)*4)+(3*total_num_cells))+2 !location of 3rd ACH for this cell
      do i = 1, NUM_BINS
         TAR_ACH_1(i)=0
      end do
      read(2,c_format,REC=REC_num,err=215) 
     +     itype,NHW,ICN,TAR_ACH_1
 215  bw_ACH_3 = FLOAT(NHW)/50.
      write (bwtext,'(f7.1)') bw_ACH_3
      bwtext5=bwtext//' ms.'
      Revents_ach2 = ICN(2)
*
*               ***** if this cell has already been analyzed, read in the information,  *****
*               *****  display it along with the CTHs and ask user to verify  *****
*
      if(((mode.eq.'ed').or.(mode.eq.'jl')).or.
     +     ((mode.eq.'cr').and.(info_prev_imported.eq.'y').or.
     +     (per_prev_imported.eq.'y')))then          
*
         histogram_number = (included(cell)*3)+1
         call read_preDBG(9,included(cell),ios,date_exp,BDT_FILE_30,
     +        recording_x,histogram_number_x,
     +        CELL_NAMES(IDs(cell)),cell,resp_type(IDs(cell)),
     +        aAP,aRL,adep,sort,comm,
     +        card(IDs(cell)),
     +        AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +        AA_rVRG,AA_rtPRG,AA_ltPRG,
     +        AA_new1,AA_new2,AA_new3,
     +        STA_phrenic,STA_RLN,STA_cent_vagus,
     +        STA_lumbar,STA_cerv_symp,STA_ELN,
     +        STA_splanch,STA_new1,STA_new2,STA_new3,

     +        carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +        hypercap_5_air_x,
     +        hypercap_tbd_x,hypoxia_12_x,
     +        hypoxia_5_x,gasp_x,lobel_x,
     +        aorta_cath_x,pe_x,carotid_occ_x,
     +        nitrop_x,vc_cath_x,dopamine_x,
     +        vagus_x,capsaicin_x,pinch_x,
     +        deep_pain_x,codeine_x,nalox_x,
     +        methyserg_x,mucimol_x,dexameth_x,
     +        noinf_x,hyperinf_x,hypervent_x,
     +        pertnew1_x,carotidCO2_2_x,
     +        cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +        sw_x,SLNsw_x,
     +        carotidCO2,vertCO2,
     +        hypercap_5_O2,
     +        hypercap_5_air,
     +        hypercap_tbd,
     +        hypoxia_12,hypoxia_5,
     +        gasp,lobel,
     +        aorta_cath,pe,
     +        carotid_occ,nitrop,
     +        vc_cath,dopamine,
     +        vagus,capsaicin,
     +        pinch,deep_pain,
     +        codeine,nalox,
     +        methyserg,mucimol,
     +        dexameth,noinf,
     +        hyperinf,hypervent,
     +        pertnew1,carotidCO2_2,
     +        cgh,lcgh,
     +        SLNcgh,exp_reflex,
     +        sw,SLNsw,per_results,AA_results,STA_results)
         if(histogram_number.ne.histogram_number_x)then
            call strlength(pre_DBG,LEN(pre_DBG),l)
            print '(//,T5,
     +           ''***** WARNING ***** DATA MIS-MATCH *****''
     +           //,T5,
     +           ''Datasave file '',A,'' may be corrupted.''
     +           ''  You need to address this issue ''
     +           ''before continuing.''
     +           /,T5,''Try re-importing the *.info file.''
     +           //,T10,
     +           ''Press <cr> to return to main menu''
     +           ''  >> '',$)',pre_DBG(1:l)
            read (*,'(A1)') OK
            return
         end if

      end if

      if(ios.gt.0)then          !if ios>0, an error occurred when reading the record
c              import_info = 'n'
c           else
c              import_info = 'y'
      end if
*
*               ***** display ACH and all CTHs for this cell;   *****
*
*
*
*               *** GET REMAINDER OF INFORMATION FOR CURRENT CELL ***
*
      enter_info = 'A'          !ask questions ... get answers

 218  ENTER_DATA: if(((mode.eq.'cr').or.(mode.eq.'ed')).and. !creating new database files, or
*                                                                      ! editing a cell with no previous data, so....
*                                                                      !ask user to specify the respiratory pattern,
*                                                                      !  sort, cardiac modulation, etc.
     +     ((analyzed_cells(IDs(cell)).eq.1).or. !incomplete data exists for this cell
     +     analyzed_cells(IDs(cell)).eq.3))then 
      call showCTHs(1,c_format,fildes,CTH,phrenic_hist,
     +     norm_CTH,norm_phrenic_hist,cardCCH,
     +     cardiac_hist,BINW,NORM_BINW,BINW2,
     +     bwtext1,bwtext2,bwtext3,bwtext4,bwtext5,
     +     Rcode,Revents_p,Revents_n,Revents_c,
     +     Tcode,Tevents,
     +     BDT_FILE,TAR_ACH,TAR_ACH_1,
     +     CELL_NAMES(IDs(cell)),QDT_FILENAME,
     +     sort,
     +     mode,resp_type(IDs(cell)),aAP,aRL,adep,
     +     AA,AA_results,AA_text1,AA_applied,
     +     STA,STA_results,STA_text1,STA_applied,
     +     STA_res_text,per_results,per_text_abbrev,
     +     zmodsig_6(IDs(cell)),zmodsig2_6(IDs(cell)),
     +     ETA2_6(IDs(cell)),coef(IDs(cell)),
     +     card_type(IDs(cell)),card(IDs(cell)),
     +     DELTA2(IDs(cell)),tedfactor(IDs(cell)),
     +     fiveHT(IDs(cell)),mean_rISI(IDs(cell)),
     +     sd_rISI(IDs(cell)),c_MAX_INT,
     +     num_rej_ISI(IDs(cell)),BINW_1,BINW_3,
     +     Revents_ach1,Revents_ach2,IDs,ITAL,cardiac_pls,
     +     STIM_OFFSET,NORM_OFFSET,mean_E,comm,show_pulses)

 220  if(CELL_NAMES(IDs(cell)).eq.' ')then 
         print '(//,T15,''!!! WARNING !!!''
     +        //,T5,''Either there is no data for IDcode '',I3,
     +        '' in the *.info file or you have not yet ''
     +        ''imported the file.''
     +        //,T5,''You will continue to get this error ''
     +        ''message until you include IDcode '',I3,
     +        '' in the *.info file''
     +        /,T5,''and re-import the file ...or... '',
     +        ''exclude '',I3,'' from further analysis.'',
     +        //,T10,''<cr> to continue   >> '',$)',cell,cell,cell
         read (*,'(A1)')OK
         if(OK.ne.' ')goto 220
      end if

      if((enter_info.eq.'A').or.(enter_info.eq.'R'))then
         call picktypetest1 (pattern,CELL_NAMES(IDs(cell)),cell) !enter respiratory activity pattern
         if(pattern.eq.'exclude')then
            menu_choice = 'X'
            goto 269
         end if
         resp_type(IDs(cell)) = pattern !  for this cell
      end if
*

 234  if((enter_info.eq.'A').or.(enter_info.eq.'S'))then
         PRINT '(/,T5''Based on visual inspection of the ACH, ''
     +        ''the sort is:''
     +        //,T10,''c  -->  clean''
     +        /,T10,''m  -->  messy''
     +        //,T15,''>> '',$)'
         read (*,fmt='(A1)',err=234) sort
         call lower_case(sort,LEN(sort))
         if((sort.ne.'c').and.(sort.ne.'m'))goto 234 !force an appropriate response
      end if
*
 236  if((enter_info.eq.'A').or.(enter_info.eq.'C'))then
         card(IDs(cell))=' '
         if(cardiac_pls.ne.0)then
            print '(/,T5,''Upon visual inspection, does this cell''
     +           '' exhibit cardiac-modulated activity? (Y/N)  >> '',$)'
            read (*,fmt='(A)',err=236) card(IDs(cell))
            call upper_case(card(IDs(cell)),LEN(card(IDs(cell))))
            if(card(IDs(cell)).eq.'Y')then
               card(IDs(cell))='c'
            else if(card(IDs(cell)).eq.'N')then
               card(IDs(cell))='n'
            else
               goto 236
            end if
         end if
      end if
      end if ENTER_DATA         !end of loop

*       ***** Display data for this cell: *****
*

 250  call showCTHs(1,c_format,fildes,CTH,phrenic_hist,
     +     norm_CTH,norm_phrenic_hist,cardCCH,
     +     cardiac_hist,BINW,NORM_BINW,BINW2,
     +     bwtext1,bwtext2,bwtext3,bwtext4,bwtext5,
     +     Rcode,Revents_p,Revents_n,Revents_c,
     +     Tcode,Tevents,
     +     BDT_FILE,TAR_ACH,TAR_ACH_1,
     +     CELL_NAMES(IDs(cell)),QDT_FILENAME,
     +     sort,
     +     mode,resp_type(IDs(cell)),aAP,aRL,adep,
     +     AA,AA_results,AA_text1,AA_applied,
     +     STA,STA_results,STA_text1,STA_applied,
     +     STA_res_text,per_results,per_text_abbrev,
     +     zmodsig_6(IDs(cell)),zmodsig2_6(IDs(cell)),
     +     ETA2_6(IDs(cell)),coef(IDs(cell)),
     +     card_type(IDs(cell)),card(IDs(cell)),
     +     DELTA2(IDs(cell)),tedfactor(IDs(cell)),
     +     fiveHT(IDs(cell)),mean_rISI(IDs(cell)),
     +     sd_rISI(IDs(cell)),c_MAX_INT,
     +     num_rej_ISI(IDs(cell)),BINW_1,BINW_3,
     +     Revents_ach1,Revents_ach2,IDs,ITAL,cardiac_pls,
     +     STIM_OFFSET,NORM_OFFSET,mean_E,comm,show_pulses)

      if(card(IDs(cell)).eq.'c')then
         card1='yes'
      else if(card(IDs(cell)).eq.'n')then
         card1='no'
      else
         card1=' '
      end if

      if(sort.eq.'m')then
         sort1='messy'
      else if(sort.eq.'c')then
         sort1='clean'
      else
         sort1=' '
      end if

      if(IDmismatch.eq.'y')then
         if((mode.ne.'vt').and.(import_info.ne.'y'))
     +        PRINT '(10(/),T50,''Coordinates (anatomic):'',    
     +        /,T5,''1. Cell name:  '',A,
     +        ''  (ID code = '',I3,'')'',
     +        T53,''4. AP'',T65,A,
     +        /,T5,''2. Cell type:  '',A,
     +        T53,''5. RL'',T65,A,
     +        /,T5,''3. Cardiac-related (visual inspection):  '',A,
     +        T53,''6. depth'',T65,A,
     +        /,T5,''7. Sort:  '',A,
     +        T22,''8. Comments:  '',A)',
     +        CELL_NAMES(IDs(cell)),cell,aAP,resp_type(IDs(cell)),
     +        aRL,card(IDs(cell)),adep,sort1,comm
         if((mode.ne.'vt').and.(import_info.eq.'y'))
     +        PRINT '(10(/),
     +        T40,''Coordinates (anatomic):'',      
     +        /,T5,''   IDcode: '',I4,
     +        T43,''   AP'',T55,A,
     +        T65,''3. Cardiac-related (vis. insp.):  '',A3,
     +        /,T5,''   Name:  '',A,
     +        T43,''   RL'',T55,A,
     +        T65,''7. Sort:  '',A5,
     +        /,T5,''2. Type:  '',A,
     +        T43,''   depth'',T55,A,
     +        /, T5,''   Comments:  '',A)',
     +        cell,aAP,card1,CELL_NAMES(IDs(cell)),
     +        aRL,sort1,resp_type(IDs(cell)),adep,comm
      end if

 251  print '(20(/))'
      print '(T5,''This is cell # '',I3,'' of '',I3,
     +     '' total cells.''/)',included(cell),total_num_cells

      if(CTHsaved.eq.'y')print '(T5,''Data for cell '',
     +     A,'' has been saved.'',/)',CELL_NAMES(IDs(cell))
      if(CTHsaved.eq.'n')print '(T5,''Data for cell '',
     +     A,'' has NOT been saved.'',/)',CELL_NAMES(IDs(cell))
      
      print '(T5,''CHOOSE ONE:'')'
      if(CTHsaved.eq.'n')print 
     +     '(T21,''s --> Save and go to the next cell'',/
     +     T21,''l --> save and return to the Last menu'',/)'

      if(mode.eq.'jl'.and.comm.ne.' ')then
         call strlength(comm,LEN(comm),l)
         print '(T10,''Comment: '',A,/)',comm(1:l)
      end if
      if(((CTHsaved.eq.'y').and.(mode.eq.'ed')).or.
     +     (mode.eq.'jl').or.(mode.eq.'vt'))print
     +     '(T18,''<cr> --> go to next cell'',/)'

      if((mode.eq.'cr').or.(mode.eq.'ed'))then
         print'(T21,''e --> Edit'')'
         if(analyzed_cells(IDs(cell)).lt.3)print
     +        '(T21,''x --> eXclude '',A4,'' from further analysis'')',
     +        CELL_NAMES(IDs(cell))

      end if

      if(mode.ne.'cr')print
     +     '(T21,''d --> Directly access a cell''''s CTHs''
     +     /,T21,''b --> Backup to previous cell'')'

      if(show_pulses.ne.'n')then
         print '(T21,''t --> remove E and I pulses'')'
      else
         print '(T21,''t --> display E and I pulses'')'
      end if
      print '(T21,''p --> Print the window'',
     +     /,T21,''w --> Write the window to disk in *.ps format'',
     +     /,T21,''v --> print/write preView'',
     +     /,T21,''r --> Redraw''
     +     /,T21,''q --> view plot of spikes/cycle'')'
      call strlength(comm,LEN(comm),l)
      if((mode.eq.'cr').or.(mode.eq.'ed'))
     +     print '(T21,''c --> edit Comment ('',A,'')'')',comm(1:l)

      if((mode.eq.'cr').or.(mode.eq.'ed'))then
         if(CTHsaved.eq.'y')then
            print '(T21,''m --> return to last Menu'')'
         else
            print '(T21,''m --> return to last Menu without ''
     +           ''saving'')'
         end if
      else
         print '(T21,''m --> return to previous Menu'')'
      end if

      print '(/,T28,''>> '',$)'
      
      old_choice = ' '
      menu_choice= ' '
      read (*,fmt='(A)',err=250) menu_choice
      call upper_case(menu_choice,LEN(menu_choice))

*     ***** the following choice/mode combinations are not allowed: *****

c          if((CTHsaved.eq.'n').and.(menu_choice.eq.' '))goto 251
      if((mode.eq.'cr').and.(menu_choice.eq.' '))goto 251
      if((mode.eq.'cr').and.(menu_choice.eq.'B '))goto 251
      if((CTHsaved.eq.'n').and.(mode.eq.'ed').and.
     +     (menu_choice.eq.' '))goto 251
      if((menu_choice.eq.' ').and.(CTHsaved.eq.'n'))goto 251
      if((mode.eq.'jl').or.(mode.eq.'vt'))then
         if((menu_choice.eq.'S').or.(menu_choice.eq.'E').or.
     +        (menu_choice.eq.'X').or.(menu_choice.eq.'L').or.
     +        (menu_choice.eq.'C'))goto 251
      end if
      if((mode.eq.'cr').and.(menu_choice.eq.'D'))goto 251
*     ***********

      if((menu_choice.eq.'X').and.(resp_type(IDs(cell)).ne.' '))then
         menu_choice = 'Z'
         old_choice = 'X'
      end if
      
 269  if(((menu_choice.eq.'S').or.(menu_choice.eq.'L').or.
     +     (menu_choice.eq.'Z')).or.
     +     ((menu_choice.eq.' ').and.(mode.eq.'ed')).or.
     +     ((menu_choice.eq.'M').and.(mode.eq.'ed').and.
     +     (changed_comment.eq.'y')))then !SAVE data 
 270     if(excluded(cell).eq.1)then      
            if(old_choice.ne.'X')then
               print '(//,5x,''Data for cell '',A4,
     +              '' (ID = '',I3,'') has been previously ''
     +              ''excluded.''
     +              //,5x,''Do you now wish to include this cell ''
     +              ''AND ALL ITS PAIRS ''
     +              ''in the database files? (Y/N)''
     +              ''  >> '',$)',
     +              CELL_NAMES(IDs(cell)),cell
               OK = ' '
               read (*,'(A)',err=270)OK
               call upper_case(OK,LEN(OK))
            else
               OK = 'Y'
            end if
            if(OK.eq.'N')then   !NO - make sure it is excluded
               excluded(cell) = 1
               goto 2001        !get next cell
            else if(OK.eq.'Y')then !YES - make sure it is included
               excluded(cell)=0
               if((analyzed_cells(IDs(cell)).eq.3).or. !flag this cell as valid 
     +              (analyzed_cells(IDs(cell)).eq.4))
     +              analyzed_cells(IDs(cell))=
     +              analyzed_cells(IDs(cell))-2
               do i2 = 1, MAX_NUM_CODES !flag this cell's pairs as valid
                  if((IDs(i2).ne.0).and.(excluded(i2).eq.0).and.
     +                 (cell.ne.i2))then
                     if((analyzed_pairs(IDs(cell),IDs(i2)).eq.3).or.
     +                    (analyzed_pairs(IDs(cell),IDs(i2)).eq.4))
     +                    analyzed_pairs(IDs(cell),IDs(i2))=
     +                    analyzed_pairs(IDs(cell),IDs(i2))-2
                     if((analyzed_pairs(IDs(i2),IDs(cell)).eq.3).or.
     +                    (analyzed_pairs(IDs(i2),IDs(cell)).eq.4))
     +                    analyzed_pairs(IDs(i2),IDs(cell))=
     +                    analyzed_pairs(IDs(i2),IDs(cell))-2
                  end if
               end do
               call write_DBSAV(pgm_version,DBSAV,ios,included, !write the db gamesave file
     +              excluded,analyzed_cells,analyzed_pairs,
     +              CELL_NAMES,perturb_applied,
     +              perturb,AA_applied,AA,STA_applied,STA,
     +              resp_type,info_prev_imported,
     +              per_prev_imported,total_num_cells_db,
     +              qdtfilenames)
            else
               goto 270         !force an appropriate answer
            end if
         end if

         CTHsaved = 'y'
         analyzed_cells(IDs(cell)) = 2 !tag this cell as ANALYZED

         call write_DBSAV(pgm_version,DBSAV,ios,included, !write the db gamesave file
     +        excluded,analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)

         call write_preDBG(9,included(cell),ios,date_exp,BDT_FILE, !write the cell's info to the preDBG file
     +        recording,histogram_number,CELL_NAMES(IDs(cell)),cell,
     +        resp_type(IDs(cell)),aAP,aRL,adep,sort,
     +        comm,
     +        card(IDs(cell)),AA_cord,                          
     +        AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +        AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,
     +        AA_new3,STA_phrenic,STA_RLN,
     +        STA_cent_vagus,STA_lumbar,STA_cerv_symp,
     +        STA_ELN,STA_splanch,STA_new1,STA_new2,
     +        STA_new3,carotidCO2_x,vertCO2_x,
     +        hypercap_5_O2_x,hypercap_5_air_x,
     +        hypercap_tbd_x,hypoxia_5_x,hypoxia_12_x,
     +        gasp_x,lobel_x,aorta_cath_x,pe_x,
     +        carotid_occ_x,nitrop_x,vc_cath_x,
     +        dopamine_x,vagus_x,capsaicin_x,pinch_x,
     +        deep_pain_x,codeine_x,nalox_x,
     +        methyserg_x,mucimol_x,dexameth_x,
     +        noinf_x,hyperinf_x,hypervent_x,
     +        pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,
     +        SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +        carotidCO2,vertCO2,
     +        hypercap_5_O2,hypercap_5_air,
     +        hypercap_tbd,hypoxia_12,
     +        hypoxia_5,gasp,lobel,
     +        aorta_cath,pe,carotid_occ,
     +        nitrop,vc_cath,dopamine,
     +        vagus,capsaicin,pinch,
     +        deep_pain,codeine,nalox,
     +        methyserg,mucimol,dexameth,
     +        noinf,hyperinf,hypervent,
     +        pertnew1,carotidCO2_2,cgh,
     +        lcgh,SLNcgh,exp_reflex,
     +        sw,SLNsw)

         if((menu_choice.eq.'Z').and.(old_choice.ne.' '))then
            menu_choice = old_choice
            old_choice = ' '
            if(menu_choice.eq.'X')goto 269
         end if

      else if(menu_choice.eq.' ')then !move to next cell if CTH has been saved or if just looking

      else if(menu_choice.eq.'R')then
         goto 250               !redraw the graphics window
         
      else if(menu_choice.eq.'M')then !return to previous menu
         
      else if(menu_choice.eq.'D')then !enter direct-access mode
         if(CTHsaved.eq.'n')then
            print '(/,T5,''Do you want to save this cell''''s ''
     +           ''data before going to another cell? ''
     +           ''(y/n)  >> '',$)'
            read (*,'(A1)') OK
            call upper_case(OK,LEN(OK))
            if(OK.eq.'Y')then
               old_choice = 'D'
               menu_choice = 'Z'
               goto 269
            else if(OK.eq.'N')then
            else
            end if
         end if
         
      else if(menu_choice.eq.'B')then !backup to previous cell's CTHs
         if(CTHsaved.eq.'n')then
            print '(/,T5,''Do you want to save this cell''''s ''
     +           ''data before backing up? (y/n)  >> '',$)'
            read (*,'(A1)') OK
            call upper_case(OK,LEN(OK))
            if(OK.eq.'Y')then
               old_choice = 'B'
               menu_choice = 'Z'
               goto 269
            else if(OK.eq.'N')then
            else
            end if
         end if

      else if(menu_choice.eq.'T')then !toggle for display of I and E pulses on CTHs
         if(show_pulses.ne.'n')then
            show_pulses='n'
         else
            show_pulses='y'
         end if
         goto 250

      else if(menu_choice.eq.'P'.or.menu_choice.eq.'W' !print or write (or view in ps format) the CTH display
     +        .or.menu_choice.eq.'V')then
         if(menu_choice.eq.'P')task='p'
         if(menu_choice.eq.'W')task='w'
         if(menu_choice.eq.'V')task='v' !'view ps file' is a debug option
         call print_or_write_CTHs(task,c_format,CTH,
     +        phrenic_hist,
     +        norm_CTH,norm_phrenic_hist,cardCCH,
     +        cardiac_hist,BINW,NORM_BINW,BINW2,
     +        bwtext1,bwtext2,bwtext3,bwtext4,bwtext5,
     +        Rcode,Revents_p,Revents_n,Revents_c,
     +        Tcode,Tevents,
     +        BDT_FILE,TAR_ACH,TAR_ACH_1,
     +        CELL_NAMES(IDs(cell)),QDT_FILENAME,
     +        sort,
     +        mode,resp_type(IDs(cell)),aAP,aRL,adep,
     +        AA,AA_results,AA_text1,AA_applied,
     +        STA,STA_results,STA_text1,STA_applied,
     +        STA_res_text,per_results,per_text_abbrev,
     +        zmodsig_6(IDs(cell)),zmodsig2_6(IDs(cell)),
     +        ETA2_6(IDs(cell)),coef(IDs(cell)),
     +        card_type(IDs(cell)),card(IDs(cell)),
     +        DELTA2(IDs(cell)),tedfactor(IDs(cell)),
     +        fiveHT(IDs(cell)),mean_rISI(IDs(cell)),
     +        sd_rISI(IDs(cell)),c_MAX_INT,
     +        num_rej_ISI(IDs(cell)),BINW_1,BINW_3,
     +        Revents_ach1,Revents_ach2,IDs,ITAL,cardiac_pls,
     +        STIM_OFFSET,NORM_OFFSET,mean_E,comm,show_pulses)
         
         goto 251

*
      else if(menu_choice.eq.'E')then !do not accept displayed data - edit it!
c 260     if((import_info.eq.'y').or.(info_prev_imported.eq.'y').or.
c     +        (import_per.eq.'y').or.(per_prev_imported.eq.'y'))then
 260     print '(
     +         //,5x,''You may change the following items.  ''
     +         ''If you wish to change others, you must ''
     +         /,5x,''edit the INFO and/or PER.DB files and ''
     +         ''re-import them.''
     +         //,T10,''R --> Respiratory modulation (firing pattern)''
     +         /,T10,''C --> Cardiac modulation''
     +         /,T10,''S --> Sort evaluation''
     +         /,T10,''A --> change All''
     +         //,T10,''<cr> --> return to the previous menu without ''
     +         ''making any changes''
     +         ///,T10,''>> '',$)'
c         end if
         enter_info=' '         !re-set the variable
         read (*,'(A1)',err=260)enter_info
         call upper_case(enter_info,LEN(enter_info))
         if(enter_info.eq.' ')then !do not make any changes - return to previous menu
            goto 251
         else if((enter_info.eq.'R').or.(enter_info.eq.'S').or.
     +           (enter_info.eq.'C'))then !change cell type
            CTHsaved = 'n'
            analyzed_cells(IDs(cell))=1 !re-tag cell as "not yet analyzed"
            goto 218
         else if(enter_info.eq.'A')then !change all information
            CTHsaved = 'n'
            resp_type(IDs(cell))=' '
            analyzed_cells(IDs(cell))=1 !re-tag cell as "not yet analyzed"
            card(IDs(cell))=' '
            sort=' '
            goto 218            !re-enter all data for this cell
         else
            goto 260            !force an appropriate response
         end if
*
      else if(menu_choice.eq.'Q')then
         call clear_top_slot(fildes)
         call plot_spikes_per_cycle(sp_per_cycle,num_acc_cycles,
     +        Tcode,Tcode,fildes,IDs,400.,445.,50.)
         call make_picture_current(fildes)
         goto 251

      else if(menu_choice.eq.'X')then !set the flag in excluded() -- ignore this cell henceforward
 261     print '(//,5x,''Please confirm -- ''
     +        ''do you wish to exclude this cell ''
     +        ''AND ALL ITS PAIRS ''
     +        ''from the database files? (Y/N)''
     +        ''  >> '',$)'
         OK = ' '
         read (*,'(A)',err=261)OK
         call upper_case(OK,LEN(OK))
         if(OK.eq.'Y')then      !YES - make sure it is excluded
            excluded(cell)=1
            if((analyzed_cells(IDs(cell)).eq.1).or.
     +           (analyzed_cells(IDs(cell)).eq.2))
     +           analyzed_cells(IDs(cell))=
     +           analyzed_cells(IDs(cell))+2
            do m = 1,MAX_NUM_CODES
               if((included(m).ne.0).and.(excluded(m).eq.0))then
                  if((analyzed_pairs(IDs(cell),IDs(m)).eq.1).or.
     +                 (analyzed_pairs(IDs(cell),IDs(m)).eq.2))
     +                 analyzed_pairs(IDs(cell),IDs(m))=
     +                 analyzed_pairs(IDs(cell),IDs(m))+2
                  if((analyzed_pairs(IDs(m),IDs(cell)).eq.1).or.
     +                 (analyzed_pairs(IDs(m),IDs(cell)).eq.2))
     +                 analyzed_pairs(IDs(m),IDs(cell))=
     +                 analyzed_pairs(IDs(m),IDs(cell))+2
               end if
            end do
            call write_DBSAV(pgm_version,DBSAV,ios,included, !write the db gamesave file
     +           excluded,analyzed_cells,analyzed_pairs,
     +           CELL_NAMES,perturb_applied,
     +           perturb,AA_applied,AA,STA_applied,STA,
     +           resp_type,info_prev_imported,
     +           per_prev_imported,total_num_cells_db,
     +           qdtfilenames)
         else if(OK.eq.'N')then
         else
            goto 261            !force an appropriate answer
         end if
      elseif(menu_choice.eq.'C')then !edit the comment
         print '(/,T5,''Enter comment (40 char max)  >> '',$)'
         read '(A40)',comm
         changed_comment='y'
         goto 250
      else
         goto 251               !force a valid response
      end if

*     ***** now look at menu_choice again to see if need to go anywhere: *****

 2001 NAVIGATE: if((menu_choice.eq.'M').or.
     +     (menu_choice.eq.'L'))then
         call clear(fildes)
         close = gclose(fildes)
         goto 1150              !C_P_M_menu

      else if(menu_choice.eq.'D')then !direct access mode
 2020    print '(/,10x,''Enter ID code of cell ..or.. ''
     +        ''<cr> to return to menu  >> '',$)'
         read (*,'(I12)',err=2020)icell
         if(icell.eq.0)goto 251 !back to menu
         if(included(icell).eq.0)then
            print '(//,15x,''You have entered an invalid ID code.'',
     +           ''  Try again.'')'
            goto 2020
         end if
         cell = icell
         goto 210               !proceed with direct access of CTHs
      else if(menu_choice.eq.'B')then !back up to previous valid IDcode
         do i = cell-1, 0, -1
            if(i.eq.0)then
               call clear(fildes)
               text='Now you''ve gone too far - '//
     +              'press <return> to continue'
               call text2d(fildes,.3,.5,text//char(0),VDC_TEXT,0)
               call make_picture_current(fildes)
               read (*,'(A)')OK
               call clear(fildes)
               goto 250
            end if
            if((IDs(i).eq.0).or.(excluded(i).eq.1))cycle !invalid IDcode or an excluded IDcode
            if((mode.ne.'vt').and.
     +           ((analyzed_cells(IDs(i)).eq.0).or.
     +           (analyzed_cells(IDs(i)).eq.3).or.
     +           (analyzed_cells(IDs(i)).eq.4)))cycle
            cell = i            !found the next one going backwards, so display it
            call clear(fildes)
            if(.true.)goto 210  !conditional to suppress warning from fort77
         end do

         call clear(fildes)
         goto 210
      else if((menu_choice.eq.' ').or.(menu_choice.eq.'X').or.
     +        (menu_choice.eq.'S'))then
         do k = cell+1,MAX_NUM_CODES !find the next valid, non-excluded cell
            if((IDs(k).ne.0).and.(excluded(k).ne.1))then
               cell = k
               if((mode.eq.'cr').and. !only interested in seeing unanalyzed cells if creating
     +              (analyzed_cells(IDs(cell)).eq.2))cycle
               goto 210         !included cell - display it
            else
               cycle            !this cell no good - look at the next possible cell
            end if
         end do

         call clear(fildes)
         call text2d(fildes,.3,.5,
     +        "End of data"//char(0),VDC_TEXT,0)
         call make_picture_current(fildes)
         print '(10(/),''End of data.  ''
     +        ''Press <cr> to continue  >> ''$)'
         read (*,'(A)')
         call clear(fildes)
         close = gclose(fildes)
         goto 1150              !C_P_M_menu
      else
      end if NAVIGATE

*
      if(mode.eq.'cr')goto 1150 !C_P_M_menu         
*
*       **********************************************************************
*       **********************************************************************
*       **********************************************************************
*       **********************************************************************
*
*       *****  Get ready to display CCHs for pair analysis and entry of results:  *****
*
 298  if((mode.eq.'ed').and.(pairs_completed.eq.'n'))then
         print '(//,''In case you''''re wondering, ''
     +        ''PAIR ANALYSIS IS INCOMPLETE.  '',/,T5,
     +        ''Would you like to see a list of unanalyzed pairs?''
     +        '' (y/n)  >> '',$)'
         read (*,'(A)') OK
         call upper_case(OK,LEN(OK))
         if(OK.eq.'Y')then
            j2 = 0
            print '(//,''UNANALYZED PAIRS by IDcode:'',/)'
            do i = 1, MAX_NUM_CODES          
               if((IDs(i).eq.0).or.(excluded(i).eq.1))cycle !do not consider excluded cells
c                  if(analyzed_cells(IDs(i)).gt.0)then
               if((analyzed_cells(IDs(i)).eq.1).or.
     +              (analyzed_cells(IDs(i)).eq.2))then
                  do j = i+1, MAX_NUM_CODES
                     if((IDs(j).eq.0).or.(excluded(j).eq.1))cycle !do not consider excluded cells
                     if(analyzed_pairs(IDs(i),IDs(j)).eq.1)then
                        write (c_i,'(I3)') i
                        write (c_j,'(I3)') j
                        text='('//c_i//','//c_j//')'
                        call remove_all_blanks(text,LEN(text))
                        call strlength(text,LEN(text),l)
                        print '(A,''  '',$)',text(1:l)
                        j2 = j2 + 1
                        if(j2.eq.10)then !print 10 pairs, then skip to next line
                           print '(/$)'
                           j2 = 0
                        end if
                     end if   
                  end do
               end if
            end do
         else if(OK.eq.'N')then
         else
            goto 298            !force an appropriate response
         end if
 2198    print '(//,T5,''Would you like to enter data for ''
     +        ''the unanalyzed pairs only? (y/n) >> '',$)'
         read (*,'(A1)') OK
         call upper_case(OK,LEN(OK))
         if(OK.eq.'Y')then
            print '(/,T5,''OK.  Just so you know -- ''
     +           ''You will be shown ''
     +           ''ONLY the unanalyzed pairs ''
     +           /,T30,''until you ''
     +           ''return to the CELLS/PAIRS menu.''
     +           //,T10,''Press <cr> to continue  >> '',$)'
            read (*,'(A)')
            mode_original = mode
            mode = 'cr'         !temporarily switch to "create" mode
         else if(OK.eq.'N')then !leave well enough alone
         else
            goto 2198           !force an appropriate response
         end if
      end if

*     ***** PUT CODE HERE TO EXPLAIN HOT KEYS FOR MAIN CCH WINDOW *****
      print '(//,''You may use the following hot keys for ''
     +     ''the main CCH screen:''
     +     //,T3,''<sp> or n --> Next pair''
     +     T40,''b --> Backup to previous pair''
     +     T71,''p --> Print the window''
     +     /,T11,''s --> next Significant pair''
     +     T40,''d --> Direct access''
     +     T71,''w --> Write the window''
     +     /,T11,''f --> next Flat pair''
     +     T40,''a --> go to first pair''
     +     T71,''v --> preView the window''
     +     /,T5,'', and . --> toggle between qdt files''
     +     T40,''2 --> show +/- 2 stds''
     +     T71,''3 --> show +/- 3 stds'',//)'


c$$$         height=height1*0.8
c$$$     write (c_height,'(I6)') height
c$$$     geometry=' '
c$$$     geometry=c_width//'x'//c_height//'-700-5'
c$$$         call remove_all_blanks(geometry,LEN(geometry))
c$$$         call strlength(geometry,LEN(geometry),l)
c$$$         call strlength(USER,LEN(USER),m)
c$$$             isys=SYSTEM ('xwcreate -wmdir /dev/screen/'//USER(1:m)//
c$$$     +          ' -title '//WINDOW1//
c$$$     +          ' -geometry ='//geometry(1:l)//char(0))
c$$$         fildes=gopen ('/dev/screen/'//USER(1:m)//'/'//
c$$$     +              WINDOW1//char(0),
c$$$     +      OUTINDEV,'sox11'// char(0),INIT)
      fildes=gopen(width,height1*4/5,-700,-5,'XAnalysis'//char(0))
*
*       *****   set screen and window parameters  *****
*
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call shade_mode(fildes,IOR(INIT,CMAP_NORMAL),0)
      call mapping_mode(fildes,1)
      call view_port(fildes,0.,0.,1.,1.)
      call view_window(fildes,1.,1.,1700.,800.)
      call text_font_index(fildes,6)
      call background_color(fildes,1.,1.,1.) !white background
      call clear(fildes)        !display the new background color immediately
      call text_color(fildes,0.,0.,0.) !text and line color is black (usually)
      call line_color(fildes,0.,0.,0.)
      call perimeter_color(fildes,0.,0.,0.)
*
      call echo_type(fildes,0,1,.625,.1,0.0)
      call set_locator(fildes,1,.625,.1,0.0)
      call track(fildes,fildes,1)
*
      call character_width(fildes,0.010)
      call character_height(fildes,0.050)
*

*
*
*               ***** read the database "gamesave" file: *****
*
*               
      if((mode.eq.'ed').or.(mode.eq.'cr').or.(mode.eq.'jl'))then 
         call read_DBSAV(file_version,DBSAV,ios,included,excluded,
     +        analyzed_cells,analyzed_pairs,CELL_NAMES,
     +        perturb_applied,
     +        perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +        info_prev_imported,
     +        per_prev_imported,total_num_cells_db,qdtfilenames)
      end if
*
*               *****  Take a look at all possible pairs and figure out which ones  *****
*               *****  are valid for analysis:  *****
*
      pair_counter=0
*
      GET_REF:  do j=1,MAX_NUM_CODES
         if(included(j).eq.0)cycle GET_REF !this cell not valid - get another one
         GET_TAR: do k=j+1,MAX_NUM_CODES
            if(included(k).eq.0)cycle GET_TAR !this cell not valid - get another one
*
            if((mode.eq.'cr').and.(analyzed_pairs(IDs(j),IDs(k)).eq.2))
     +           cycle GET_TAR
            if(((mode.eq.'cr').or.(mode.eq.'ed').or.(mode.eq.'jl')).and.
     +           (analyzed_pairs(IDs(j),IDs(k)).eq.0).or.
     +           (analyzed_pairs(IDs(j),IDs(k)).ge.3))cycle GET_TAR !analyzed_pairs(j,k):
                                !=0 --> not a valid pair
                                !=1 --> valid, but not analyzed yet
                                !=2 --> valid and analyzed
                                !=3 --> excluded, not analyzed yet
                                !=4 --> excluded and analyzed

*
            REF = j             !"good" code --> define as REFerence cell
            TAR = k             !"good" code --> define as TARget cell
*
*               ***** Load the valid_pair_REF/TAR arrays with ID codes according to *****
*               *****   the location of the pair within the "pair queue". Cells that*****
*               *****   have been designated as EXCLUDED are not included in this   *****
*               *****   process and will not be viewed in the "slide-show" display. *****
*               *****   However, if CCH data exists in the .qdt file for any cell   *****
*               *****   it can be viewed using the "direct access" option.          *****
*
            pair_counter = pair_counter+1 !update the counter
            valid_pair_REF(pair_counter) = REF
            valid_pair_TAR(pair_counter) = TAR
            total_num_pairs = pair_counter !total_num_pairs to be analyzed for this file
*
         end do GET_TAR
      end do GET_REF
*
      total_num_sig = 0
      call filesize (pre_DBP, ifilesize)
      iPairsInFile = ifilesize/iDBP_RECL
      do m = 1, iPairsInFile
         call read_preDBP(10,m,ios,char2,char11,
     +        char30,char4,integer,char4,integer,
     +        prim1,char15,char20,char15,char40,
     +        char8,char8,char8,char8,char8,char8,char40,char8,
     +        integer)
         if((prim1.ne.'Flat').and.(prim1.ne.' '))then
            total_num_sig = total_num_sig + 1
         end if
      end do



      pair_counter = 0          !reset the counter
      call clear(fildes)        !clear display
      auto_mode = 'n'           !set the flag:  y --> show CCHs automatically (slide show)
                                !               n --> pause the slide show for analysis (default)
*
*
*       **************************************************************************
*       **************************************************************************
*       **                                                                      **
*       ** BEGIN HUGE DO LOOP that will read in the histograms computed for     **
*       **  each pair of cells. The user will be prompted to ENTER the RESULTS  **
*       **  of the analysis. This information will be stored for later write-out**
*       **  into the database files.                                            **
*       **                                                                      **
*       **                                                                      **
*       **                                                                      **
*       **************************************************************************
*       **************************************************************************
*
*

      call sc_file (included, DB_FILES)

      show_conf_lim = 'y'       !default display values: show the averaged shift-control +/- 3 standard deviations
      show_2_sd = 'n'
      show_3_sd = 'y'
      show_control = 'y'
      show_avg_shift = 'y'
      show_single_shift = 'n'
      scaledMAIN = 'n'
      istart = 1
      flats_only = 'n'
      speed = 3.
      GET_NEXT_PAIRx: do ii = istart,total_num_pairs !for every possible pair of cells:
         pair_counter = ii
*
 350     if((pair_counter.le.0).or.
     +        (pair_counter.gt.total_num_pairs))then
            call clear1(fildes,2) !clear bottom of window
            textstring='End of CCH data'
            call text2d(fildes,50.,50.,textstring//char(0),
     +           ANNOTATION_TEXT,0)
            textstring='Click to continue'
            call text2d(fildes,60.,30.,textstring//char(0),
     +           ANNOTATION_TEXT,0)
            call make_picture_current(fildes)
            call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
            call clear(fildes)
            close = gclose(fildes)
            sig_only = 'n'
            goto 1150           !C_P_M_menu                              !return to C/P/X menu
         end if
*
         REF = valid_pair_REF(pair_counter) !look in the array to find REF and TAR ID codes
         TAR = valid_pair_TAR(pair_counter)
         if((mode.eq.'cr').and. !allow user to back up to any pair -
     +        (analyzed_pairs(IDs(REF),IDs(TAR)).ne.1).and. ! but when move forward, look only at unanalyzed pairs
     +        (button_choice.ne.'prev_pair'))then ! if the user has specifically chosen that option
            pair_counter = pair_counter + 1
            goto 350
         end if
*
*               ***** re-initialize variables for statistical analysis *****
*
 380     ZK=0.0
         Probk=0.0
         det=0.0
         vis=0.0
         ZLAT=0.0
         HALFWD=0.0
         stats_bw = ' '
         czk=' '
         cprobk=' '
         cdet=' '
         cvis=' '
         czlat=' '
         chalfwd=' '
         prim=' '
         sec=' '
         loc1=' '
         loc2=' '
         rescom=' '
         statcomm=' '
         do i = 1, MAX_PERTURB
            per_results_REF (i)= ' '
            per_results_TAR(i) = ' '
         end do
         
*
*               ***** calculate record number in .qdt file for 1st CCH for this pair *****
*               *****   (see comment lines at beginning of this subroutine)     *****
*
         call sc_pair (REF, TAR, total_num_qdts)
c380     Q_pos = 0
         Q_pos = 0
         if(included(REF).eq.1)goto 390
         do i2 = 1,included(REF)-1                     
            Q_pos = Q_pos + (total_num_cells-i2)         
         end do
 390     Q_pos = Q_pos + (included(TAR) - included(REF))
         rel_loc = (Q_pos*4)-3
         REC_num_x = 3 + (7*total_num_cells) + rel_loc !hold the # of the CCH at first binwidth
         REC_num_control_1 = 3 + (5*total_num_cells)+ !record # of "shift by 1" control CCH @ 1st binwidth
     +        (2*(total_num_cells**2)) + rel_loc !hold the # of the control CCH @ 1st bw
         REC_num_control_2 = 3 + (3*total_num_cells)+ !record # of "shift by 2" control CCH @ 1st binwidth
     +        (4*(total_num_cells**2)) + rel_loc !hold the # of the control CCH @ 1st bw
         REC_num_control_avg = 3 + total_num_cells + !record # of "averaged shift" control CCH @ 1st binwidth
     +        (6*(total_num_cells**2)) + rel_loc !hold the # of the control CCH @ 1st bw
*
*
         if((mode.eq.'cr').or.(mode.eq.'ed').or.(mode.eq.'jl'))then !do not read a DATASAVE file if looking at a .TQDT file
            call read_preDBP(10,Q_pos,ios,recording_x,date_exp,
     +           BDT_FILE_30,cell_name(1),IDs_1(1),cell_name(2),IDs_1(2)
     +           ,prim,loc1,sec,loc2,rescom,czk,cprobk,cdet,cvis ,czlat
     +           ,chalfwd,statcomm,stats_bw,REC_num_xx)
            prim_old = prim
            if(REC_num_x.ne.REC_num_xx)then
               call strlength(pre_DBP,LEN(pre_DBP),l)
               print '(//,T5,
     +              ''***** WARNING ***** DATA MIS-MATCH *****''
     +           //,T5,''Datasave file '',A,'' may be corrupted.''
     +           ''  You need to address this issue before continuing.''
     +           /,T5,''Try re-importing the *.info file.''
     +           //,T10,''Press <cr> to return to main menu  >> '',$)',
     +           pre_DBP(1:l)
               read (*,'(A1)') OK
               return
            end if
         end if

         if((sig_only.eq.'y')
     +        .and.((prim.eq.'Flat').or.(prim.eq.' ')))then !if user has so chosen, show only CCHs with features
            pair_counter = pair_counter + 1                      
            goto 350
         end if

         if((flats_only.eq.'y')
     +        .and.(prim.ne.'Flat'))then !if user has so chosen, show only CCHs with no features
            pair_counter = pair_counter + 1                      
            goto 350
         end if

         write (c_num_sig,'(I8)') total_num_sig
         call remove_all_blanks(c_num_sig,LEN(c_num_sig))
         call strlength(c_num_sig,LEN(c_num_sig),l_s)
         ipercent = (FLOAT(total_num_sig)/FLOAT(total_num_pairs))*100.
         write (c_percent,'(I3)') ipercent
         call strlength(c_percent,LEN(c_percent),l_p)
         print '(/,T5,''This is pair '',I8,'' of '',I8,
     +        '' total pairs.  (approx '',A,'' ('',A,''%)'',
     +        '' significant correlations)'')',
     +        pair_counter,total_num_pairs,c_num_sig(1:l_s),
     +        c_percent(1:l_p)

         if((mode.eq.'cr').or.(mode.eq.'ed').or.(mode.eq.'jl'))then
            histogram_number = (included(REF)*3)+1 !calculate the correct record number - check for inconsistencies
            call read_preDBG(9,included(REF),ios,date_exp,
     +           BDT_FILE_30,recording_x,histogram_number_x,
     +           CELL_NAMES(IDs(REF)),cell,resp_type(IDs(REF)),
     +           aAP,aRL,adep,sort,comm,
     +           card(IDs(REF)),
     +           AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +           AA_rVRG,AA_rtPRG,AA_ltPRG,
     +           AA_new1,AA_new2,AA_new3,
     +           STA_phrenic,STA_RLN,STA_cent_vagus,
     +           STA_lumbar,STA_cerv_symp,STA_ELN,
     +           STA_splanch,STA_new1,STA_new2,STA_new3,
     +           carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +           hypercap_5_air_x, hypercap_tbd_x,hypoxia_12_x,
     +           hypoxia_5_x,gasp_x,lobel_x,
     +           aorta_cath_x,pe_x,carotid_occ_x,
     +           nitrop_x,vc_cath_x,dopamine_x,
     +           vagus_x,capsaicin_x,pinch_x,
     +           deep_pain_x,codeine_x,nalox_x,
     +           methyserg_x,mucimol_x,dexameth_x,
     +           noinf_x,hyperinf_x,hypervent_x,
     +           pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +           exp_reflex_x,sw_x,SLNsw_x,carotidCO2,vertCO2,
     +           hypercap_5_O2,hypercap_5_air,hypercap_tbd,
     +           hypoxia_12,hypoxia_5,gasp,lobel,aorta_cath,pe,
     +           carotid_occ,nitrop,vc_cath,dopamine,vagus,capsaicin,
     +           pinch,deep_pain,codeine,nalox,methyserg,mucimol,
     +           dexameth,noinf,hyperinf,hypervent,pertnew1,
     +           carotidCO2_2,
     +           cgh,lcgh,SLNcgh,exp_reflex,sw,SLNsw,per_results_REF,
     +           AA_results,STA_results)
            REF_coords(1) = aAP
            REF_coords(2) = aRL
            REF_coords(3) = adep
            if(histogram_number.ne.histogram_number_x)then
               call strlength(pre_DBG,LEN(pre_DBG),l)
               print '(//,T5,
     +              ''***** WARNING ***** DATA MIS-MATCH *****''
     +              //,T5,
     +              ''Datasave file '',A,'' may be corrupted.''
     +              ''  You need to address this issue ''
     +              ''before continuing.''
     +              /,T5,''Try re-importing the *.info file.''
     +              //,T10,
     +              ''Press <cr> to return to main menu''
     +              ''  >> '',$)',pre_DBG(1:l)
               read (*,'(A1)') OK
               return
            end if
            histogram_number = (included(TAR)*3)+1
            call read_preDBG(9,included(TAR),ios,date_exp,
     +           BDT_FILE_30,recording_x,histogram_number_x,
     +           CELL_NAMES(IDs(TAR)),cell,resp_type(IDs(TAR)),
     +           aAP,aRL,adep,sort,comm,
     +           card(IDs(TAR)),
     +           AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +           AA_rVRG,AA_rtPRG,AA_ltPRG,
     +           AA_new1,AA_new2,AA_new3,
     +           STA_phrenic,STA_RLN,STA_cent_vagus,
     +           STA_lumbar,STA_cerv_symp,STA_ELN,
     +           STA_splanch,STA_new1,STA_new2,STA_new3,
     +           carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +           hypercap_5_air_x, hypercap_tbd_x,hypoxia_12_x,
     +           hypoxia_5_x,gasp_x,lobel_x,
     +           aorta_cath_x,pe_x,carotid_occ_x,
     +           nitrop_x,vc_cath_x,dopamine_x,
     +           vagus_x,capsaicin_x,pinch_x,
     +           deep_pain_x,codeine_x,nalox_x,
     +           methyserg_x,mucimol_x,dexameth_x,
     +           noinf_x,hyperinf_x,hypervent_x,
     +           pertnew1_x,carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +           exp_reflex_x,sw_x,SLNsw_x,carotidCO2,vertCO2,
     +           hypercap_5_O2,hypercap_5_air,hypercap_tbd,
     +           hypoxia_12,hypoxia_5,gasp,lobel,aorta_cath,pe,
     +           carotid_occ,nitrop,vc_cath,dopamine,vagus,capsaicin,
     +           pinch,deep_pain,codeine,nalox,methyserg,mucimol,
     +           dexameth,noinf,hyperinf,hypervent,pertnew1,
     +           carotidCO2_2,
     +           cgh,lcgh,SLNcgh,exp_reflex,sw,SLNsw,per_results_TAR,
     +           AA_results,STA_results)
            TAR_coords(1) = aAP
            TAR_coords(2) = aRL
            TAR_coords(3) = adep
            if(histogram_number.ne.histogram_number_x)then
               call strlength(pre_DBG,LEN(pre_DBG),l)
               print '(//,T5,
     +              ''***** WARNING ***** DATA MIS-MATCH *****''
     +              //,T5,
     +              ''Datasave file '',A,'' may be corrupted.''
     +              ''  You need to address this issue ''
     +              ''before continuing.''
     +              /,T5,''Try re-importing the *.info file.''
     +              //,T10,
     +              ''Press <cr> to return to main menu''
     +              ''  >> '',$)',pre_DBG(1:l)
               read (*,'(A1)') OK
               return
            end if
         end if
*
*
*
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*               *                                                               *
*               *       DISPLAY CROSS-CORRELOGRAMS AND ENTER ANALYSIS RESULTS   *
*               *                                                               *
*               * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*       
         flag = 0               !flag to tell program that the user has
                                ! chosen to investigate an individual CCH more
                                ! closely, perhaps for statistical analysis of a
                                ! primary feature; 
                                ! if =1 --> display ENTER ANALYSIS RESULTS
                                ! button in red to prompt user to enter results
                                ! if the CCH is not flat

 400     DISPLAY_CCHs:   do izz=1,4
            call sc_izz (izz)
            do i = 1, 6
               ICN(i)=0
            end do

            do i = 1, NUM_BINS
               CCH (i)= 0       !clear the CCH histogram array
            end do
            read(2,c_format,REC=(REC_num_x+(izz-1)),err=401) 
     +           itype,NHW,ICN,CCH !read in the CCH data

 401        if(itype.ne.0)print '(''NOT A CCH'')'
            BINW_x=(FLOAT(NHW)/50.)
            if(izz.eq.1)then
               Rcode=ICN(1)
               Tcode=ICN(3)
               if((Rcode.ne.REF).or.(Tcode.ne.TAR))then
                  call strlength(pre_DBP,LEN(pre_DBP),l)
                  print '(//,T5,
     +              ''***** WARNING ***** DATA MIS-MATCH *****''
     +              //,T5,''Datasave file '',A,'' may be corrupted.''
     +              ''  You must address this issue before continuing.''
     +              //,T10,''Press <cr> to return to main menu  >> '',
     +              $)', pre_DBP(1:l)
                  read (*,'(A1)') OK
                  return
               end if
            end if

            ibin_count=0
            do i5=1,NUM_BINS
               ibin_count=CCH(i5)+ibin_count
            end do
            do i = 1, NUM_BINS
               CONTROL(i)=0
            end do
            read(2,c_format,REC=(REC_num_control_1+(izz-1)),err=402)  !read in the single shift-control CCH
     +           itype,NHW,ICN,CONTROL
            single_shift = 'y'
            ibin_count=0
            do i5=1,NUM_BINS
               ibin_count=CONTROL(i5)+ibin_count
            end do
            if(total_histograms.eq.(8*(total_num_cells**2) -
     +           total_num_cells+3))then !2-cycle and averaged shift-controls are in the qdt file
               do i = 1, NUM_BINS
                  CONTROL_2 (i)= 0
                  CONTROL_AVG (i)= 0
               end do
               avg_shift = 'y'
               read(2,c_format,REC=(REC_num_control_2+(izz-1)),
     +              IOSTAT=ios,err=1401) !read in the 2-cycle shift-control CCH
     +              itype,NHW,ICN,CONTROL_2
 1401          read(2,c_format,REC=(REC_num_control_avg+(izz-1)),
     +              IOSTAT=ios,err=402) !read in the averaged shift-control CCH
     +              itype,NHW,ICN,CONTROL_AVG
            else
               avg_shift = 'n'
               show_avg_shift = 'n' !can't show what's not there
            end if

*       **** calculate record numbers of ACHs and read in ACH data: *****

 402        REC_num=(included(REF)*4)+(3*total_num_cells)
            do i = 1, NUM_BINS
               REF_ACH(i)=0
            end do
            read(2,c_format,REC=(REC_num+(izz-1)),err=403) 
     +           itype,NHW,ICN,REF_ACH
*
 403        REC_num=(included(TAR)*4)+(3*total_num_cells)
            do i = 1, NUM_BINS
               TAR_ACH(i)=0
            end do
            read(2,c_format,REC=(REC_num+(izz-1)),err=404) 
     +           itype,NHW,ICN,TAR_ACH
*
 404        if(izz.eq.1)then    !go and get the normalized CTH plots
               REC_num=(included(REF)*3)+2
               do i = 1, NUM_BINS
                  CTH_REF(i)=0
               end do
               read(2,c_format,REC=REC_num,err=405) itype,NHW,ICN
     +              ,CTH_REF
 405           REC_num=(included(TAR)*3)+2
               do i = 1, NUM_BINS
                  CTH_TAR(i)=0
               end do
               read(2,c_format,REC=REC_num,err=4051) itype,NHW,ICN,
     +              CTH_TAR
 4051          REC_num=(included(REF)*3)+3
               do i = 1, NUM_BINS
                  cardiac_REF(i)=0
               end do
               read(2,c_format,REC=REC_num,err=4052) itype,NHW,ICN,
     +              cardiac_REF
 4052          REC_num=(included(TAR)*3)+3
               do i = 1, NUM_BINS
                  cardiac_TAR(i)=0
               end do
               read(2,c_format,REC=REC_num,err=406) itype,NHW,ICN,
     +              cardiac_TAR
            end if
*
 406        button_choice = ' '
            call showCCHs(fildes,IDs,mouse,mode, CCH,CONTROL,REF_ACH
     +           ,TAR_ACH,ICN,BINW_x,izz,CTH_REF,CTH_TAR
     +           ,norm_phrenic_hist,cardiac_REF,cardiac_TAR
     +           ,cardiac_hist,ETA2_6(IDs(Rcode)),zmodsig_6(IDs(Rcode))
     +           ,zmodsig2_6(IDs(Rcode)),card_type(IDs(Rcode)),
     +           ETA2_6(IDs(Tcode)),zmodsig_6(IDs(Tcode)),
     +           zmodsig2_6(IDs(Tcode)),card_type(IDs(Tcode)),Rcode
     +           ,Tcode,CELL_NAMES(IDs(Rcode)),CELL_NAMES(IDs(Tcode)),
     +           resp_type(IDs(Rcode)),resp_type(IDs(Tcode)),
     +           button_choice,REC_num_x,REC_num_control_1,auto_mode,
     +           speed,coefnum(IDs(TAR)),click,prim,sec,loc1,loc2
     +           ,rescom,Q_pos,czk,cprobk,cdet,cvis,chalfwd,czlat,ZK
     +           ,stats_bw,Probk,det,vis,ZLAT,HALFWD,BDT_FILE
     +           ,QDT_FILENAME,date_exp,recording,NHW,id,flag,WINDOW1
     +           ,WINDOW2,show_control,c_format,statcomm,qdt_files
     +           ,current_qdt,included,total_num_cells,DELTA2(IDs(REF))
     +           ,DELTA2(IDs(TAR)),tedfactor(IDs(REF))
     +           ,tedfactor(IDs(TAR)),per_text_abbrev,per_results_REF
     +           ,per_results_TAR,CONTROL_AVG,REC_num_control_2
     +           ,REC_num_control_avg,single_shift,avg_shift
     +           ,show_single_shift,show_avg_shift,CONTROL_2
     +           ,sp_per_cycle,num_acc_cycles,show_2_sd,show_3_sd
     +           ,show_conf_lim,sig_only,ITAL,cardiac_pls,flats_only
     +           ,total_num_qdts,scaledMAIN,REF_coords,TAR_coords
     +           ,exp_name,NORM_OFFSET,mean_E,NORM_BINW,show_pulses)


            if(button_choice.eq.'redraw')goto 400
            if(button_choice.eq.'leave_stats')goto 400 !re-display same pair

*               ***** transformation of analysis results to character format: *****

            if(izz.eq.4)then    !when all CCHs, ACHs, and CTHs have been displayed:
               if((prim.ne.' ').and.
     +              ((mode.eq.'cr').or.(mode.eq.'ed')))then !transform data only if user has entered
                  if(((prim.eq.'Peak').or.(prim.eq.'Trough')) !  analysis results
     +                 .and.(ZK.ne.0.0))then
                     write (czk,'(F6.2)') ZK
                     write (cprobk,'(F6.2)') Probk
                     write (cvis,('(F6.2)')) vis
                     write (czlat,('(F6.2)')) ZLAT
                     write (chalfwd,('(F6.2)')) HALFWD
                  end if
                  if((prim.eq.'M P & T').or.(prim.eq.'Flat'))then
                     cdet = ' '
                     czk = ' '
                     cprobk = ' '
                     cvis = ' '
                     czlat = ' '
                     chalfwd = ' '
                     statcomm = ' '
                  end if
               end if

               if(prim_old.ne.prim)then
                  if((prim_old.eq.' '.or.prim_old.eq.'Flat')
     +                 .and.(prim.ne.' '.and.prim.ne.'Flat'))
     +                 total_num_sig = total_num_sig + 1
                  if((prim_old.ne.' '.and.prim_old.ne.'Flat')
     +                 .and.(prim.eq.' '.or.prim.eq.'Flat'))
     +                 total_num_sig = total_num_sig - 1
               end if
               

*      
*
*       **************** WRITE TO .pre_dbp FILE ***********************************
*
*
               if((mode.eq.'ed').or.(mode.eq.'cr'))then !update the database gamesave file
                  call write_preDBP(10,Q_pos,ios,recording,date_exp,
     +                 BDT_FILE,CELL_NAMES(IDs(REF)),REF,
     +                 CELL_NAMES(IDs(TAR)),TAR,
     +                 prim,loc1,sec,loc2,rescom,czk,cprobk,
     +                 cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +                 REC_num_x)

                  analyzed_pairs(IDs(REF),IDs(TAR)) = 2 !2 --> this pair now tagged as "analyzed"

                  call write_DBSAV(pgm_version,DBSAV,ios,included,
     +                 excluded,analyzed_cells,analyzed_pairs,
     +                 CELL_NAMES,perturb_applied,
     +                 perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +                 info_prev_imported,
     +                 per_prev_imported,total_num_cells_db,
     +                 qdtfilenames)

                  if(excluded(REF).eq.1)then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' is included.'',
     +                    /,T5,''I will store these data in the ''
     +                    ''meantime. '')',REF
                     goto 420   !skip ahead to checking on the TARget cell
                  end if
                  if(analyzed_cells(IDs(REF)).ne.2)then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' has been analyzed.''
     +                    /,T5,''(ie., respiratory pattern / sort ''
     +                    ''evaluation / visual cardiac ''
     +                    ''modulation have not been entered)''
     +                    /,T5,''I will store these data in the '',
     +                    ''meantime.'')',REF
                  end if
                  if(CELL_NAMES(IDs(REF)).eq.' ')then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' has been assigned a name. ''
     +                    /,T5,''(ie., either the *.info file has '',
     +                    ''not yet been imported or IDcode '',I3,
     +                    '' was not in the *.info file)''
     +                    /,T5,''I will store these data in the '',
     +                    ''meantime.'')',REF,REF
                  end if

 420              if(excluded(TAR).eq.1)then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' is included.'',
     +                    /,T5,''I will store these data in the ''
     +                    ''meantime. '')',TAR
                     goto 430
                  end if
                  if(analyzed_cells(IDs(TAR)).ne.2)then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' has been analyzed.''
     +                    /,T5,''(ie., respiratory pattern / sort ''
     +                    ''evaluation / visual cardiac ''
     +                    ''modulation have not been entered)''
     +                    /,T5,''I will store these data in the '',
     +                    ''meantime.'')',TAR
                  end if
                  if(CELL_NAMES(IDs(TAR)).eq.' ')then
                     print '(//,T5,''*** WARNING ***'',
     +                    //,T5,''Analysis data for this pair will ''
     +                    ''not be written to the ''
     +                    ''database-ready tables until cell '', 
     +                    I3,'' has been assigned a name. ''
     +                    /,T5,''(ie., either the *.info file has '',
     +                    ''not yet been imported or IDcode '',I3,
     +                    '' was not in the *.info file)''
     +                    /,T5,''I will store these data in the '',
     +                    ''meantime.'')',TAR,TAR
                  end if


               end if
*
*
 430           if(button_choice.eq.'prev_pair')then !auto_mode = 'n' for all of these:
                  pair_counter = pair_counter - 1
                  goto 350

               elseif(button_choice.eq.'first_pair')then 
                  pair_counter = 1
                  goto 350
*     
               else if(button_choice.eq.'next_pair')then
                  pair_counter = pair_counter + 1
                  goto 350
*
               else if(button_choice.eq.'main_menu')then
                  call clear(fildes)
                  close = gclose(fildes)
                  goto 1150     !C_P_M_menu
*
               else if(button_choice.eq.'dir_acc')then
                  REF=0
                  TAR=0
 450              print '(//,10x,''Enter ID code of reference cell ''
     +                 ''..or.. <cr> to exit direct access mode >> ''
     +                 ,$)'
                  read (*,fmt='(I3)',err=450)REF
                  if(REF.eq.0) goto 350
                  if(included(REF).eq.0)then
                     print '(/,10x,
     +                    ''You have entered an invalid code.'')'
                     goto 450
                  end if
 460              print '(//,10x,
     +                 ''Enter ID code of target cell  >> '',$)'
                  read (*,fmt='(I3)',err=460)TAR
                  if(included(TAR).eq.0)then
                     print '(/,10x,
     +                    ''You have entered an invalid code.'')'
                     goto 460
                  end if
                  if(REF.eq.TAR)then
                     print '(/,10x,
     +                    ''Target must be different from reference.'')'
                     goto 460
                  end if
                  if(TAR.lt.REF)then !switch 'em!
                     i=REF
                     j=TAR
                     REF=j
                     TAR=i
                  end if
*
*                                                       !if using direct access mode,
*                                                       !check to see if the user has asked to 
*                                                       !  view a "valid" pair
                  REFtemp = REF
                  TARtemp = TAR
                  click = 'y'   !user forced to use "clickable" mode
                  speed=21600.

                  print *,'analyzed_pairs(IDs(',REF,',IDs(',TAR,') = ',
     +                 analyzed_pairs(IDs(REF),IDs(TAR))

                  if((analyzed_pairs(IDs(REF),IDs(TAR)).eq.1).or.
     +                 (analyzed_pairs(IDs(REF),IDs(TAR)).eq.2))
     +                 goto 470 !if "valid", reset pair_counter so can view
                                !  previous and next pairs
                  if((analyzed_pairs(IDs(REF),IDs(TAR)).eq.0).or.
     +                 (analyzed_pairs(IDs(REF),IDs(TAR)).eq.3).or.
     +                 (analyzed_pairs(IDs(REF),IDs(TAR)).eq.4))then !if "invalid", find the valid pair immediately

                     print *,'No can do; invalid pair (a_p = ',
     +                    analyzed_pairs(IDs(REF),IDs(TAR)),')'
                     PREV: do m = REF,1,-1 !  preceding and use it to reset pair_counter
                        do n = TAR,1,-1 !  to enable viewing of previous and next pairs
                           if((IDs(m).ne.0).and.(IDs(n).ne.0))then
                              print *,'a_p(IDs(',m,'), IDs(',n,') = ',
     +                             analyzed_pairs(IDs(m),IDs(n))
                              if((analyzed_pairs(IDs(m),IDs(n)).eq.1).or.
     +                             (analyzed_pairs(IDs(m),IDs(n)).eq.2))
     +                             then
                                 REFtemp = m !NOTE: choices of "previous" and "next" pairs 
                                 TARtemp = n !  will show "valid" pairs ONLY
                                 exit PREV
                              end if
                           end if
                        end do
                     end do PREV
                  end if
*
 470              do m = 1,total_num_pairs
                     if((REFtemp.eq.valid_pair_REF(m)).and.
     +                    (TARtemp.eq.valid_pair_TAR(m)))then
                        pair_counter = m
                        exit
                     end if
                  end do
                  goto 380
*
               else if(button_choice.eq.'auto_view')then
                  pair_counter = pair_counter + 1
                  goto 350

               end if           !end of "if(button_choice.eq....)then"
*
            end if              !end of "if(izz.eq.4) ..."
*
         end do DISPLAY_CCHs
*
*
*       ***** write the database "gamesave" file: *****
*
         if((mode.eq.'ed').or.(mode.eq.'cr'))then
            call write_DBSAV(pgm_version,DBSAV,ios,included,
     +           excluded,analyzed_cells,analyzed_pairs,CELL_NAMES,
     +           perturb_applied,
     +           perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +           info_prev_imported,
     +           per_prev_imported,total_num_cells_db,qdtfilenames)
         end if
*

      end do GET_NEXT_PAIRx
*
      call clear(fildes)
      close = gclose (fildes)
      print '(''end of program'')'

      goto 1150                 !C_P_M_menu

      END

*
*
*       *********************************************************************
*       *********************************************************************
*       *********************************************************************
*
*
*

      end module mod_analyze_data
