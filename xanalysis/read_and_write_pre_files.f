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

      module mod_read_and_write_pre_files
      contains
*     filename: read_and_write_pre_files.f          (created 16-Dec-2003    lss)

*     date of last revision: 20-Jan-2004        lss

*     subroutines for reading and writing pre_DBG and pre_DBP files for VERSION 5 of x2002

*     pre_DBG files provide information for the UNITS and RESPONSES text files that will be read into the Access database
*     pre_DBP files provide information for the CROSSES text files that will be read into the Access database


      subroutine read_preDBG(iunit,rec_num,ios,date_exp,BDT_FILE_30,
     +     recording,hist_number,
     +     CELL_NAMES,cell,resp_type,aAP,aRL,adep,sort,comm,
     +     card,
     +     AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +     AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +     STA_phrenic,STA_RLN,STA_cent_vagus,
     +     STA_lumbar,STA_cerv_symp,STA_ELN,
     +     STA_splanch,STA_new1,STA_new2,STA_new3,
     +     carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +     hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +     carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,sw_x,
     +     SLNsw_x,
     +     carotidCO2,vertCO2,hypercap_5_O2,hypercap_5_air,
     +     hypercap_tbd,hypoxia_12,hypoxia_5,gasp,lobel,
     +     aorta_cath,pe,carotid_occ,nitrop,vc_cath,dopamine,
     +     vagus,capsaicin,pinch,deep_pain,codeine,nalox,
     +     methyserg,mucimol,dexameth,noinf,hyperinf,hypervent,
     +     pertnew1,carotidCO2_2,cgh,lcgh,SLNcgh,exp_reflex,sw,SLNsw,
     +     per_results,AA_results,STA_results)

      INCLUDE 'x2000parameter.defs'


      integer cell,ios,iunit,rec_num,
     +     hist_number,cell_temp

      character*40 comm
      character*30 BDT_FILE_30  !BDT_FILE in the pre files may be truncated from 200 to 30 char

      character*15 resp_type
      character*12 STA_phrenic,STA_RLN,STA_lumbar,STA_ELN,
     +     STA_cent_vagus,STA_splanch,STA_cerv_symp,
     +     STA_new1,STA_new2,STA_new3,STA_results(MAX_STA)

      character*11 date_exp
      character*9 carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +     sw_x,SLNsw_x,per_results(MAX_PERTURB)
      CHARACTER*5 aAP,aRL,adep
      CHARACTER*4 CELL_NAMES
      CHARACTER*3 AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +     AA_rVRG,AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +     AA_results(MAX_AA)
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
     +     recording
      CHARACTER*1 card,sort
      character*80 iomsg

      ios = 0
      cell_temp = cell
      
      do i = 1, MAX_PERTURB
         per_results(i) = ' '
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
      CELL_NAMES = ' '
      resp_type = ' '
      aAP = ' '
      aRL = ' '
      adep = ' '
      sort =  ' '
      comm = ' '
      card = ' '


      goto 100
 10   if(ios.ne.0)then
c            print '(''I''''ve gone a record too far -- 
c     +           don''''t worry, it''''s OK'')'
      else
         print '(''error reading pre_dbg file for cell '',I12)',
     +        cell_temp
         print '(''error code = '',I10)',ios
            print *,'error msg = ',iomsg
         print *,'rec_num = ',rec_num
      end if
      goto 999

c 100     read (iunit,'(I10,4(A),I5,237(A))',REC=rec_num,
 100  read (iunit,'(I10,4(A),I5,237(A))',REC=rec_num,
     +     IOSTAT=ios,ERR=10,IOMSG=iomsg),
     +     hist_number,         !histogram number of cell resp CTH
     +     date_exp,BDT_FILE_30,recording,CELL_NAMES,
     +     cell,resp_type,aAP,aRL,adep,sort,comm,
     +     card,
     +     AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +     AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +     STA_phrenic,STA_RLN,STA_cent_vagus,
     +     STA_lumbar,STA_cerv_symp,STA_ELN,
     +     STA_splanch,STA_new1,STA_new2,STA_new3,
     +     carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +     hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +     carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,codeine_x,nalox_x,
     +     methyserg_x,mucimol_x,dexameth_x,noinf_x,hyperinf_x,
     +     hypervent_x,pertnew1_x,carotidCO2_2_x,
     +     cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +     carotidCO2,vertCO2,
     +     hypercap_5_O2,
     +     hypercap_5_air,hypercap_tbd,
     +     hypoxia_12,
     +     hypoxia_5,gasp,lobel,
     +     aorta_cath,pe,
     +     carotid_occ,nitrop,
     +     vc_cath,dopamine,vagus,
     +     capsaicin,pinch,deep_pain,
     +     codeine,nalox,
     +     methyserg,mucimol,
     +     dexameth,noinf,
     +     hyperinf,hypervent,
     +     pertnew1,carotidCO2_2,
     +     cgh,lcgh,SLNcgh,
     +     exp_reflex,sw,SLNsw

      per_results(1) = carotidCO2_x
      per_results(2) = vertCO2_x
      per_results(3) = hypercap_5_O2_x
      per_results(4) = hypercap_5_air_x
      per_results(5) = hypercap_tbd_x
      per_results(6) = hypoxia_12_x
      per_results(7) = hypoxia_5_x
      per_results(8) = gasp_x
      per_results(9) = pertnew1_x
      per_results(10) = carotidCO2_2_x
      per_results(11) = lobel_x
      per_results(12) = aorta_cath_x
      per_results(13) = pe_x
      per_results(14) = carotid_occ_x
      per_results(15) = nitrop_x
      per_results(16) = vc_cath_x
      per_results(17) = dopamine_x
      per_results(18) = vagus_x
      per_results(19) = capsaicin_x
      per_results(20) = pinch_x
      per_results(21) = deep_pain_x
      per_results(22) = codeine_x
      per_results(23) = nalox_x
      per_results(24) = methyserg_x
      per_results(25) = mucimol_x
      per_results(26) = dexameth_x
      per_results(27) = noinf_x
      per_results(28) = hyperinf_x
      per_results(29) = hypervent_x
      per_results(30) = cgh_x
      per_results(31) = lcgh_x
      per_results(32) = SLNcgh_x
      per_results(33) = exp_reflex_x
      per_results(34) = sw_x
      per_results(35) = SLNsw_x

      AA_results(1) = AA_cord
      AA_results(2) = AA_RLN
      AA_results(3) = AA_vagus
      AA_results(4) = AA_cVRG
      AA_results(5) = AA_rVRG
      AA_results(6) = AA_rtPRG
      AA_results(7) = AA_ltPRG
      AA_results(8) = AA_new1
      AA_results(9) = AA_new2
      AA_results(10) = AA_new3
      STA_results(1) = STA_phrenic
      STA_results(2) = STA_RLN
      STA_results(3) = STA_cent_vagus
      STA_results(4) = STA_lumbar
      STA_results(5) = STA_cerv_symp
      STA_results(6) = STA_ELN
      STA_results(7) = STA_splanch
      STA_results(8) = STA_new1
      STA_results(9) = STA_new2
      STA_results(10) = STA_new3

 999  return
      end

*     *****************************************************************
*     *****************************************************************
*     *****************************************************************
*     *****************************************************************


      subroutine write_preDBG(iunit,rec_num,ios,date_exp,BDT_FILE_30,
     +     recording,hist_number,
     +     CELL_NAMES,cell,resp_type,aAP,aRL,adep,sort,comm,
     +     card,
     +     AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +     AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +     STA_phrenic,STA_RLN,STA_cent_vagus,
     +     STA_lumbar,STA_cerv_symp,STA_ELN,
     +     STA_splanch,STA_new1,STA_new2,STA_new3,
     +     carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +     carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,sw_x,
     +     SLNsw_x,
     +     carotidCO2,vertCO2,hypercap_5_O2,hypercap_5_air,
     +     hypercap_tbd,hypoxia_12,hypoxia_5,gasp,lobel,
     +     aorta_cath,pe,carotid_occ,nitrop,vc_cath,dopamine,
     +     vagus,capsaicin,pinch,deep_pain,codeine,nalox,
     +     methyserg,mucimol,dexameth,noinf,hyperinf,hypervent,
     +     pertnew1,carotidCO2_2,cgh,lcgh,SLNcgh,exp_reflex,sw,SLNsw)

      INCLUDE 'x2000parameter.defs'

      integer cell,ios,hist_number,iunit,
     +     rec_num

      character*40 comm
      character*30 BDT_FILE_30  !BDT_FILE may be truncated from 200 to 30 char in the pre files

      character*15 resp_type
      character*12 STA_phrenic,STA_RLN,STA_lumbar,STA_ELN,
     +     STA_cent_vagus,STA_splanch,STA_cerv_symp,
     +     STA_new1,STA_new2,STA_new3

      character*11 date_exp
      character*9 carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +     hypoxia_12_x,gasp_x,
     +     lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +     nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,
     +     codeine_x,nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +     noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +     carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +     sw_x,SLNsw_x
      CHARACTER*5 aAP,aRL,adep
      CHARACTER*4 CELL_NAMES
      CHARACTER*3 AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +     AA_rVRG,AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3
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
     +     recording

      CHARACTER*1 card,sort

      if(global_mode.eq.'jl'.and.global_compare_and_convert_flag.eq.' ')
     +     then
         print '
     +(''bug: global_compare_and_convert_flag is blank in write_preDB
     +G'')'
         stop
      end if

      if (global_mode.eq.'jl'.and.global_compare_and_convert_flag.eq
     +     .'n') return

      ios = 0
      
      goto 100
 10   if(ios.eq.922)then
c            print '(''I''''ve gone a record too far -- 
c     +           don''''t worry, it''''s OK'')'
      else
         print '(''error writing to pre_dbg file for cell '',I3)',
     +        cell
         print '(''error code = '',I10)',ios
      end if
      goto 999

 100  write (iunit,'(I10,4(A),I5,237(A))',REC=rec_num,
     +     IOSTAT=ios,ERR=10)     
     +     hist_number,         !histogram number of cell resp CTH
     +     date_exp,BDT_FILE_30,recording,CELL_NAMES,
     +     cell,resp_type,aAP,aRL,adep,sort,comm,card,
     +     AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +     AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +     STA_phrenic,STA_RLN,STA_cent_vagus,
     +     STA_lumbar,STA_cerv_symp,STA_ELN,
     +     STA_splanch,STA_new1,STA_new2,STA_new3,
     +     carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +     hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +     hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +     carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +     capsaicin_x,pinch_x,deep_pain_x,codeine_x,nalox_x,
     +     methyserg_x,mucimol_x,dexameth_x,noinf_x,hyperinf_x,
     +     hypervent_x,pertnew1_x,carotidCO2_2_x,
     +     cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,sw_x,SLNsw_x,
     +     carotidCO2,vertCO2,hypercap_5_O2,
     +     hypercap_5_air,hypercap_tbd,hypoxia_12,
     +     hypoxia_5,gasp,lobel,aorta_cath,pe,
     +     carotid_occ,nitrop,vc_cath,dopamine,vagus,
     +     capsaicin,pinch,deep_pain,codeine,nalox,
     +     methyserg,mucimol,dexameth,noinf,
     +     hyperinf,hypervent,pertnew1,carotidCO2_2,
     +     cgh,lcgh,SLNcgh,exp_reflex,sw,SLNsw


 999  return
      end


*     *****************************************************************
*     *****************************************************************
*     *****************************************************************
*     *****************************************************************


      subroutine read_preDBP(iunit,rec_num,ios,recording,date_exp,
     +     BDT_FILE_30,CELL_NAMES_REF,REF,CELL_NAMES_TAR,TAR,
     +     prim,loc1,sec,loc2,rescom,czk,cprobk,
     +     cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +     hist_number)

      integer REF,TAR,hist_number,iunit,rec_num,ios
      character*40 statcomm,rescom
      character*30 BDT_FILE_30
      character*20 prim,sec
      character*15 loc1,loc2
      character*11 exp_date,date_exp
      CHARACTER*8 czk,cprobk,cdet,cvis,chalfwd,czlat,stats_bw
      character*4 CELL_NAMES_REF,CELL_NAMES_TAR
      character*2 recording
      character*80 iomsg
      LOGICAL stop_flag

      if(.false.)print *,date_exp !suppress unused variable warning

      stop_flag = .FALSE.
      ios = 0

c         print '(''start reading preDBP for pair '',A4,'','',A4)',
c     +        CELL_NAMES_REF,CELL_NAMES_TAR
      goto 100
 10   if(ios.eq.922)then
c            print '(''I''''ve gone a record too far -- 
c     +           don''''t worry, it''''s OK'')'
      else
         print '(''error reading pre_dbp file for pair '',A4,'','',A4)',
     +        CELL_NAMES_REF,CELL_NAMES_TAR
         print '(''error code = '',I10)',ios
         print *,'iomsg = ',iomsg
         print *,'rec = ', rec_num
         stop_flag = .TRUE.
      end if
      if (stop_flag) stop
      goto 999

 100  read(iunit,'(I10,4(A),I5,A,I5,13(A))',REC=rec_num,
     +        IOSTAT=ios,ERR=10,IOMSG=iomsg),
C     +     IOSTAT=ios,ERR=10),
     +     hist_number,exp_date,recording,BDT_FILE_30,
     +     CELL_NAMES_REF,REF,CELL_NAMES_TAR,TAR,
     +     prim,loc1,sec,loc2,rescom,czk,cprobk,
     +     cdet,cvis,czlat,chalfwd,statcomm,stats_bw

c         print '(''finished reading preDBP for pair '',A4,'','',A4)',
c     +        CELL_NAMES_REF,CELL_NAMES_TAR

 999  return
      end


*     *****************************************************************
*     *****************************************************************
*     *****************************************************************
*     *****************************************************************


      subroutine write_preDBP(iunit,rec_num,ios,recording,date_exp,
     +     BDT_FILE_30,CELL_NAMES_REF,REF,CELL_NAMES_TAR,TAR,
     +     prim,loc1,sec,loc2,rescom,czk,cprobk,
     +     cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +     hist_number)


      integer REF,TAR,hist_number,iunit,rec_num,ios
      character*40 statcomm,rescom
      character*30 BDT_FILE_30
      character*20 prim,sec
      character*15 loc1,loc2
      character*11 date_exp
      CHARACTER*8 czk,cprobk,cdet,cvis,chalfwd,czlat,stats_bw
      character*4 CELL_NAMES_REF,CELL_NAMES_TAR
      character*2 recording

      include 'x2000parameter.defs'


      if(global_mode.eq.'jl'.and.global_compare_and_convert_flag.eq.' ')
     +     then
         print '
     +(''bug: global_compare_and_convert_flag is blank in write_preDB
     +P'')'
         stop
      end if

      if (global_mode.eq.'jl'.and.global_compare_and_convert_flag.eq
     +     .'n') return

      ios = 0

c         print '(''start writing preDBP for pair '',A4,'','',A4)',
c     +        CELL_NAMES_REF,CELL_NAMES_TAR
      goto 100
 10   if(ios.eq.922)then
c            print '(T5,''I''''ve gone a record too far -- '',
c     +                  /,T5,''don''''t worry it''''s OK'')'
      else
         print '(''error writing pre_dbp file for pair '',A4,'','',A4)',
     +        CELL_NAMES_REF,CELL_NAMES_TAR
         print '(''error code = '',I10)',ios
      end if
      goto 999

 100  write(iunit,'(I10,4(A),I5,A,I5,13(A))',REC=rec_num,
     +     IOSTAT=ios,ERR=10)
     +     hist_number,date_exp,recording,BDT_FILE_30,
     +     CELL_NAMES_REF,REF,CELL_NAMES_TAR,TAR,
     +     prim,loc1,sec,loc2,rescom,czk,cprobk,
     +     cdet,cvis,czlat,chalfwd,statcomm,stats_bw

c         pause

c         print '(''finished writing preDBP for pair '',A4,'','',A4)',
c     +        CELL_NAMES_REF,CELL_NAMES_TAR
c         pause

 999  return
      end
      subroutine clear_pert_fields(carotidCO2_x,vertCO2_x,
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
     +     sw(5),SLNsw(5)

      carotidCO2_x= ' ' 
      vertCO2_x= ' '
      hypercap_5_O2_x=' '
      hypercap_5_air_x=' '
      hypercap_tbd_x=' '
      hypoxia_12_x= ' '
      hypoxia_5_x= ' '
      gasp_x=' '
      lobel_x= ' '
      aorta_cath_x= ' '
      pe_x= ' '
      carotid_occ_x= ' '
      nitrop_x= ' '
      vc_cath_x= ' '
      dopamine_x=' '
      vagus_x=' '
      capsaicin_x= ' '
      pinch_x= ' '
      deep_pain_x= ' '
      codeine_x= ' '
      nalox_x= ' '
      methyserg_x= ' '
      mucimol_x= ' '
      dexameth_x=' '
      noinf_x= ' '
      hyperinf_x= ' '
      hypervent_x= ' '
      pertnew1_x=' '
      carotidCO2_2_x=' '
      cgh_x= ' '
      lcgh_x= ' '
      SLNcgh_x= ' '
      exp_reflex_x= ' '
      sw_x= ' '
      SLNsw_x= ' '
      do m = 1,5
         carotidCO2(m)= ' '     !perturbation data was read-in before - wipe out obsolete data
         vertCO2(m)= ' '
         hypercap_5_O2(m)=' '
         hypercap_5_air(m)=' '
         hypercap_tbd(m)=' '
         hypoxia_12(m)= ' '
         hypoxia_5(m)= ' '
         gasp(m)=' '
         lobel(m)= ' '
         aorta_cath(m)= ' '
         pe(m)= ' '
         carotid_occ(m)= ' '
         nitrop(m)= ' '
         vc_cath(m)= ' '
         dopamine(m)=' '
         vagus(m)=' '
         capsaicin(m)= ' '
         pinch(m)= ' '
         deep_pain(m)= ' '
         codeine(m)= ' '
         nalox(m)= ' '
         methyserg(m)= ' '
         mucimol(m)= ' '
         dexameth(m)=' '
         noinf(m)= ' '
         hyperinf(m)= ' '
         hypervent(m)= ' '
         pertnew1(m)=' '
         carotidCO2_2(m)=' '
         cgh(m)= ' '
         lcgh(m)= ' '
         SLNcgh(m)= ' '
         exp_reflex(m)= ' '
         sw(m)= ' '
         SLNsw(m)= ' '
      end do
      
      return
      end
      end module mod_read_and_write_pre_files
