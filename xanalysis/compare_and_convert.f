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

      module mod_compare_and_convert
      contains
*     filename: compare_and_convert.f

*     date of last revision: 16-Feb-2018       lss
*       check for changes in marker codes (<100) that should be reflected within excluded()!!

*     This subroutine of Xanalysis (aka x2004) will allow the user to associate existing database
*     files with new qdt files by allowing the addition and/or deletion of cells to a new qdt file.


         subroutine compare_and_convert_if_add_or_delete_units
     +     (included_QDT,included_DB,exclude_QDT,exclude_DB,
     +     IDs,pre_DBG,pre_DBP,total_num_cells_QDT,
     +     analyzed_cells,analyzed_pairs,BDT_FILE,CELL_NAMES,
     +     resp_type,date_exp,recording,qdt_files,qdt_sav,DBSAV,
     +     esc,mode,total_num_qdts,dirname)
      use mod_miscellaneous_subroutines
      use mod_read_and_write_pre_files

         include 'x2000parameter.defs'

         integer included_DB(MAX_NUM_CODES),delete(MAX_NUM_CODES),
     +        exclude_QDT(MAX_NUM_CODES),total_num_qdts,
     +        exclude_DB(MAX_NUM_CODES),add(MAX_NUM_CODES),
     +        included_QDT(MAX_NUM_CODES),
     +        analyzed_cells(MAX_NUM_CHAN),
     +        analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +        IDs_REF,IDs_TAR,IDs(MAX_NUM_CODES)
         integer ios,hist_number,total_num_cells_QDT,
     +        REF,TAR,Q_pos_DB,Q_pos_QDT,total_num_cells_DB,
     +        rel_loc,hist_number_QDT

         character*(*) pre_DBG,pre_DBP,BDT_FILE,date_exp,recording,
     +        qdt_files(MAX_NUM_QDTS),qdt_sav(MAX_NUM_QDTS),DBSAV,
     +        dirname
         character*40 statcomm,comm,rescom,blank40
         character*30 BDT_FILE_30
         character*20 prim,sec,blank20
         character*15 resp_type(MAX_NUM_CHAN),loc1,loc2,blank15
         character*12 STA_phrenic,STA_RLN,STA_lumbar,STA_ELN,
     +        STA_cent_vagus,STA_splanch,STA_cerv_symp,
     +        STA_new1,STA_new2,STA_new3,STA_results(MAX_STA),
     +        blank12
         character*9 carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +        hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +        hypoxia_12_x,gasp_x,
     +        lobel_x,aorta_cath_x,pe_x,carotid_occ_x,
     +        nitrop_x,vc_cath_x,dopamine_x,vagus_x,
     +        capsaicin_x,pinch_x,deep_pain_x,
     +        codeine_x,nalox_x,methyserg_x,mucimol_x,dexameth_x,
     +        noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +        carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,exp_reflex_x,
     +        sw_x,SLNsw_x,per_results(MAX_PERTURB),blank9
         CHARACTER*8 czk,cprobk,cdet,cvis,chalfwd,czlat,stats_bw,
     +        blank8
         CHARACTER*5 aAP,aRL,adep,blank5
         CHARACTER*4 CELL_NAME_REF,CELL_NAME_TAR,
     +        CELL_NAMES(MAX_NUM_CHAN),blank4
         CHARACTER*3 AA_cord,AA_RLN,AA_vagus,AA_cVRG,
     +        AA_rVRG,AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +        AA_results(MAX_AA),blank3
         character*2 carotidCO2(5),vertCO2(5),hypercap_5_O2(5),
     +        hypercap_5_air(5),hypercap_tbd(5),hypoxia_5(5),
     +        hypoxia_12(5),gasp(5),
     +        lobel(5),aorta_cath(5),pe(5),carotid_occ(5),
     +        nitrop(5),vc_cath(5),dopamine(5),vagus(5),
     +        capsaicin(5),pinch(5),deep_pain(5),
     +        codeine(5),nalox(5),methyserg(5),mucimol(5),
     +        dexameth(5),
     +        noinf(5),hyperinf(5),hypervent(5),pertnew1(5),
     +        carotidCO2_2(5),cgh(5),lcgh(5),SLNcgh(5),
     +        exp_reflex(5),
     +        sw(5),SLNsw(5),
     +        blank2x5(5),blank2,mode
         CHARACTER*1 card,sort,different,
     +        blank1,OK,OK1,OK2,esc

         INTEGER*4 SYSTEM
         EXTERNAL SYSTEM

*       ***** suppress unused variable warnings *****
        if(.false.)print *,exclude_DB
        if(.false.)print *,total_num_qdts
        if(.false.)print *,exclude_QDT

         do i = 1, MAX_NUM_CODES
            if(IDs(i).ne.0)print '(''IDs('',I3,'') = '',I3''; a_c = '',
     +           I1)',i,IDs(i),analyzed_cells(IDs(i))
         end do

 5       call strlength(qdt_files(1),LEN(qdt_files(1)),l_QDT)
         call strlength(DBSAV,LEN(DBSAV),l_DB)
         call strlength(pre_DBG,LEN(pre_DBG),l_dbg)
         call strlength(pre_DBP,LEN(pre_DBP),l_dbp)
         esc = 'n'

         different = 'n'
         do i = 1, MAX_NUM_CODES
            delete(i) = 0
            add(i) = 0
         end do

         do i = 1, MAX_NUM_CODES                        !check to see if data from QDTSAV and DBSAV files is different
            if(included_QDT(i).ne.included_DB(i))then
               different = 'y'
               exit
            end if
            if(exclude_QDT(i).ne.exclude_DB(i))then     !check for changes in marker codes
               different = 'y'
               exit
            end if
         end do
         global_compare_and_convert_flag = 'n'
         CHECK_FOR_CONVERT: if(different.eq.'n')then                   !don't need to convert anything
            print '(''no cells or marker codes added/deleted; do not have to convert'')'
            return
         else                                           !units have been added or deleted - CONVERT!
            print '(''CONVERT - cells added/deleted!'')'
            if(mode.eq.'jl')then                        !backup the datasave files and work with temporaries
               call strlength(pre_DBG,LEN(pre_DBG),l_dbg)
               dirname = 'xanalysis_tempXXXXXX'//char(0)
               call mkdtemp(dirname)
               isys = SYSTEM('cp '//pre_DBG(1:l_dbg)//' '//dirname
     +              //char(0))
               call strlength(pre_DBP,LEN(pre_DBP),l_dbp)
               isys = SYSTEM('cp '//pre_DBP(1:l_dbp)//' '//dirname
     +              //char(0))
               call strlength(DBSAV,LEN(DBSAV),l_db)
               isys = SYSTEM('cp '//DBSAV(1:l_db)//' '//dirname
     +              //char(0))
            else
               num_del_cells = 0
               num_added_cells = 0
               do i = 1, MAX_NUM_CODES
                  if((included_QDT(i).eq.0).and.
     +                 (included_DB(i).ne.0))then
                     delete(i) = 1
                     num_del_cells = num_del_cells + 1
                  end if
                  if((included_QDT(i).ne.0).and.
     +                 (included_DB(i).eq.0))then
                     add(i) = 1
                     num_added_cells = num_added_cells + 1
                  end if
               end do
 15            if(num_added_cells.gt.0)then
                 print '(10(/),''*** WARNING ** WARNING ** WARNING ''
     +                 ''** WARNING ** WARNING ** WARNING ***''
     +                 //,T5,A,'' contains some cells ''
     +                 ''that are not in '',A,''.'')',
     +                 qdt_files(1)(1:l_QDT),DBSAV(1:l_DB)
                  if(num_added_cells.eq.1)then
                     print '(/,T5,''You are about to include the ''
     +                    ''following cell in these database files. ''
     +                    /,T5,''You will have to enter information ''
     +                    ''for it and all its pairs.''
     +                    /,T5,''Please be sure you are adding this ''
     +                    ''cell to the correct database files.'',/)'
                     do j = 1, MAX_NUM_CODES
                        if(add(j).eq.1)print '(T10,I3)',j
                     end do
                     print '(/,T5,''Do you want to add this cell?  ''
     +                    ''(y/n..or.. V to enter view-only mode)''
     +                    ''  >> '',$)'
                  else
                     print '(/,T5,''You are about to include the ''
     +                    ''following '',I3,'' cells in these ''
     +                    ''database files. ''
     +                    /,T5,''You will have to enter information ''
     +                    ''for them and all their pairs.''
     +                    /,T5,''Please be sure you are adding these ''
     +                    ''cells to the correct database files.'',/)',
     +                    num_added_cells
                     do j = 1, MAX_NUM_CODES
                        if(add(j).eq.1)print '(T10,I3)',j
                     end do
                     print '(/,T5,''Do you want to add these cells?  ''
     +                    ''(y/n ..or.. V to enter view-only mode)''
     +                    ''  >> '',$)'
                  end if
                  read (*,'(A)'),OK
                  call upper_case(OK,LEN(OK))
                  if(OK.eq.'Y')then
                     print '(''Okey-dokey.  Here we go.'')'
                  elseif(OK.eq.'N')then
                     esc = 'y'
                     return
                  elseif(OK.eq.'V')then
                     mode = 'jl'
                     goto 5
                  else
                     goto 15
                  end if
               end if

 20            if(num_del_cells.gt.0)then
c                  call strlength(qdt_files(1),LEN(qdt_files(1)),l_QDT)
c                  call strlength(DBSAV,LEN(DBSAV),l_DB)
c                  call strlength(pre_DBG,LEN(pre_DBG),l_dbg)
c                  call strlength(pre_DBP,LEN(pre_DBP),l_dbp)
                 print '(10(/),''*** WARNING ** WARNING ** WARNING ''
     +                 ''** WARNING ** WARNING ** WARNING ***''
     +                 //,T5,A,'' is missing some cells ''
     +                 ''that are in '',A,''.'')',
     +                 qdt_files(1)(1:l_QDT),DBSAV(1:l_DB)
                  if(num_del_cells.eq.1)then
                     print '(/,T5,''You are about to ''
     +                    ''DELETE information for the following cell ''
     +                    ''and all its pairs:'',/)'
                  else
                     print '(/,T5,''You are about to ''
     +                    ''DELETE information for the following '',I3,
     +                    '' cells and all their pairs:'',/)',
     +                    num_del_cells
                  end if
                  do j = 1, MAX_NUM_CODES
                     if(delete(j).eq.1)print '(T10,I3)',j
                  end do
                  print '(/,T5,''[If you want to ''''play around'''' ''
     +                 ''with a data subset, I suggest that you copy ''
     +                 /,T5,''the following ''
     +                 ''files into a temporary directory: '',/)'
                  call strlength(qdt_files(1),LEN(qdt_files(1)),l_qdt)
                  print '(T8,A,$)',qdt_files(1)(1:l_qdt)
                  do m = 2, MAX_NUM_QDTS
                     if(qdt_files(m).ne.' ')then
                        call strlength(qdt_files(m),LEN(qdt_files(m)),
     +                       l_qdt)
                        print '(3x,A,$)',qdt_files(m)(1:l_qdt)
                     end if
                  end do
                  call strlength(qdt_sav(1),LEN(qdt_sav(1)),l_sav)
                  print '(T8,A,$)',qdt_sav(1)(1:l_sav)
                  do m = 2, MAX_NUM_QDTS
                     if(qdt_sav(m).ne.' ')then
                        call strlength(qdt_sav(m),LEN(qdt_sav(m)),
     +                       l_sav)
                        print '(3x,A,$)',qdt_sav(m)(1:l_sav)
                     end if
                  end do
                  print '(/,T8,A,5x,A,5x,A,
     +                 //,''and play with them there.]'')',
     +                 DBSAV(1:l_DB),pre_DBG(1:l_dbg),pre_DBP(1:l_dbp)
                  print '(/,T5,''Are you absolutely sure you want to ''
     +                 ''delete data?  ''
     +                 ''(y/n ..or.. V to enter view-only mode)''
     +                 ''  >> '',$)'
                  read (*,'(A)'),OK
                  call upper_case(OK,LEN(OK))
 30               if(OK.eq.'Y')then
                     print '(10(/),T5,''Because if you mess up these ''
     +                    ''files, it could go badly for you -- ''
     +                 /,T5,''reprimand ... dismissal ... that sort ''
     +                    ''of thing.'')'
                     do m = 1, 650000000
                     end do
                     print '(/,T5,''I am not making this up.  ''
     +                    ''You can ask Dr. Lindsey.'')'
                     do m = 1, 650000000
                     end do
                     print '(/,T5,''Are you still sure that deleting ''
     +                    ''these data is the right thing to do?  '',
     +                    //,T8,''(y/n/V for view mode)  >> '',$)'
                     read (*,'(A)'),OK1
                     call upper_case(OK1,LEN(OK1))
 40                  if(OK1.eq.'Y')then
                        print '(10(/),T5,''OK, then.  Here we go.''
     +                       /,T5,''If you ''
     +                       ''later regret these deletions ''
     +                       ''(and want to save your job), ''
     +                       /,T5,''you might ''
     +                       ''be able to restore data files ''
     +                       ''from the ''
     +                       ''./x2000/backups subdirectory.'',
     +                       //,T10,''<cr> to continue ..or.. ''
     +                       ''X to change your mind (last chance)''
     +                       ''  >> '',$)'
                        read (*,'(A)') OK2
                        call upper_case(OK2,LEN(OK2))
                        if(OK2.eq.'X')then
                           esc = 'y'
                           return
                        elseif(OK2.eq.' ')then
                        else
                           goto 40
                        end if
                     elseif(OK1.eq.'N')then
                        esc = 'y'
                        return
                     elseif(OK1.eq.'V')then
                        mode = 'jl'
                        goto 5
                     else
                        goto 30
                     end if
                  elseif(OK.eq.'N')then
                     esc = 'y'
                     return
                  elseif(OK.eq.'V')then
                     mode = 'jl'
                     goto 5
                  else
                     goto 20
                  end if
               end if
            end if

*     **** OK - we're gonna add or delete cells from the DB files to match the QDT files: *****

            global_compare_and_convert_flag = 'y'
            close (unit = 9)                            !close the already-open preDBG file
            close (unit = 10)                            !close the already-open preDBP file
            blank40 = ' '
            blank20 = ' '
            blank15 = ' '
            blank12 = ' '
            blank9 = ' '
            blank8 = ' '
            blank5 = ' '
            blank4 = ' '
            blank3 = ' '
            blank2 = ' '
            blank1 = ' '
            do i = 1, 5
               blank2x5 = ' '
            end do
            do i = 1, MAX_NUM_CHAN
               analyzed_cells(i) = 0                        !now clear the array to get ready for updated info
               CELL_NAMES(i) = ' '
               resp_type(i) = ' '
               do j = 1, MAX_NUM_CHAN
                  analyzed_pairs(j,i) = 0
               end do
            end do

*     ***** convert the pre_DBG file: *****

            call strlength(pre_DBG,LEN(pre_DBG),l)
            isys=SYSTEM('mv '//pre_DBG//' '//pre_DBG(1:l)//
     +           '.trash'//char(0))                                    !move old files out of the way
            OPEN (30,FILE=pre_DBG(1:l)//'.trash',
     +           FORM='FORMATTED',                                     !open the existing pre_dbg file
     +           STATUS='OLD',ACCESS='DIRECT',RECL=DBG_RECL)
            OPEN (31,FILE=pre_DBG,FORM='FORMATTED',               !create the new pre_dbg file
     +           STATUS='NEW',ACCESS='DIRECT',RECL=DBG_RECL)

            do i = 1, MAX_NUM_CODES
               if(included_QDT(i).ne.0)then                            !if this cell is in the new QDT file 
                  ios = 0
                  hist_number_QDT=(included_QDT(i)*3)+1        !calculate the histogram number of resp CTH

                  if(included_DB(i).ne.0)then                  !if this cell is in the DBSAV file, too,  then ...
                     i_tmp = i
                     call read_preDBG(30,included_DB(i),ios,date_exp,  !read its information from the old pre_dbg file and...
     +                    BDT_FILE_30,
     +                    recording,hist_number,
     +                    CELL_NAMES(IDs(i)),i_tmp,resp_type(IDs(i)),
     +                    aAP,aRL,adep,sort,comm,card,
     +                    AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +                    AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +                    STA_phrenic,STA_RLN,STA_cent_vagus,
     +                    STA_lumbar,STA_cerv_symp,STA_ELN,
     +                    STA_splanch,STA_new1,STA_new2,STA_new3,
     +                    carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +                    hypercap_5_air_x,hypercap_tbd_x,hypoxia_12_x,
     +                    hypoxia_5_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +                    carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,
     +                    vagus_x,
     +                    capsaicin_x,pinch_x,deep_pain_x,
     +                    codeine_x,nalox_x,methyserg_x,mucimol_x,
     +                    dexameth_x,
     +                    noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +                    carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +                    exp_reflex_x,sw_x,SLNsw_x,
     +                    carotidCO2,vertCO2,hypercap_5_O2,
     +                    hypercap_5_air,
     +                    hypercap_tbd,hypoxia_12,hypoxia_5,gasp,lobel,
     +                    aorta_cath,pe,carotid_occ,nitrop,vc_cath,
     +                    dopamine,
     +                    vagus,capsaicin,pinch,deep_pain,codeine,nalox,
     +                    methyserg,mucimol,dexameth,noinf,hyperinf,
     +                    hypervent,
     +                    pertnew1,carotidCO2_2,cgh,lcgh,SLNcgh,
     +                    exp_reflex,sw,SLNsw,per_results,AA_results,
     +                    STA_results)

                     if((ios.eq.0).and.(i_tmp.eq.i))then          !read the record from the old preDBG with no errors
                        BDT_FILE_30 =  BDT_FILE
                        call write_preDBG(31,included_QDT(i),ios,    ! write that same info to the new pre_dbg file
     +                   date_exp,BDT_FILE_30,                       !    using its new record number and histogram #
     +                   recording,hist_number_QDT,
     +                   CELL_NAMES(IDs(i)),i,resp_type(IDs(i)),
     +                   aAP,aRL,adep,sort,comm,card,
     +                   AA_cord,AA_RLN,AA_vagus,AA_cVRG,AA_rVRG,
     +                   AA_rtPRG,AA_ltPRG,AA_new1,AA_new2,AA_new3,
     +                   STA_phrenic,STA_RLN,STA_cent_vagus,
     +                   STA_lumbar,STA_cerv_symp,STA_ELN,
     +                   STA_splanch,STA_new1,STA_new2,STA_new3,
     +                   carotidCO2_x,vertCO2_x,hypercap_5_O2_x,
     +                   hypercap_5_air_x,hypercap_tbd_x,hypoxia_5_x,
     +                   hypoxia_12_x,gasp_x,lobel_x,aorta_cath_x,pe_x,
     +                   carotid_occ_x,nitrop_x,vc_cath_x,dopamine_x,
     +                   vagus_x,
     +                   capsaicin_x,pinch_x,deep_pain_x,
     +                   codeine_x,nalox_x,methyserg_x,mucimol_x,
     +                   dexameth_x,
     +                   noinf_x,hyperinf_x,hypervent_x,pertnew1_x,
     +                   carotidCO2_2_x,cgh_x,lcgh_x,SLNcgh_x,
     +                   exp_reflex_x,sw_x,SLNsw_x,
     +                   carotidCO2,vertCO2,hypercap_5_O2,
     +                   hypercap_5_air,
     +                   hypercap_tbd,hypoxia_12,hypoxia_5,gasp,lobel,
     +                   aorta_cath,pe,carotid_occ,nitrop,vc_cath,
     +                   dopamine,
     +                   vagus,capsaicin,pinch,deep_pain,codeine,nalox,
     +                   methyserg,mucimol,dexameth,noinf,hyperinf,
     +                   hypervent,
     +                   pertnew1,carotidCO2_2,cgh,lcgh,SLNcgh,
     +                   exp_reflex,sw,SLNsw)
                        if((resp_type(IDs(i)).ne.' ').and.(sort.ne.' ')
     +                       .and.(card.ne.' '))then                    
                           analyzed_cells(IDs(i)) = 2                  !this cell has been analyzed
                        else
                           analyzed_cells(IDs(i)) = 1                  !this cell has not been analyzed yet
                        end if
                        cycle
                     else
                        print '(''error reading record for '',I3,
     +                       '': ios = '',I5,
     +                       '', but it''''s OK - I''''ll just '',
     +                       ''write a place-holder record'')',i,ios
                     end if
                  end if

*     ***** If get to this point, either the cell is an "added" one, or it does exist in the old preDBG file but its record  *****
*     *****  wasn't read without an error, so write a "place-holder" record to the new preDBG file:   *****

                  analyzed_cells(IDs(i)) = 1                           !tag this cell as valid, but not yet analyzed
                  BDT_FILE_30 = BDT_FILE
                  call write_preDBG(31,included_QDT(i),ios,date_exp,   ! write a "place-holder" record to the new pre_dbg file
     +                   BDT_FILE_30,                                               
     +                   recording,hist_number_QDT,
     +                   blank4,i,blank15,
     +                   blank5,blank5,blank5,blank1,blank40,blank1,
     +                   blank3,blank3,blank3,blank3,blank3,blank3,
     +                   blank3,blank3,blank3,blank3,
     +                   blank12,blank12,blank12,blank12,blank12,
     +                   blank12,blank12,blank12,blank12,blank12,
     +                   blank9,blank9,blank9,blank9,blank9,blank9,
     +                   blank9,blank9,blank9,blank9,blank9,blank9,
     +                   blank9,blank9,blank9,blank9,blank9,blank9,
     +                   blank9,blank9,blank9,blank9,blank9,blank9,
     +                   blank9,blank9,blank9,blank9,blank9,blank9,
     +                   blank9,blank9,blank9,blank9,blank9,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5,
     +                   blank2x5,blank2x5,blank2x5,blank2x5,blank2x5)
                    
               end if
            end do
           close (unit=30)
           close (unit=31)


*     ***** convert the pre_DBP file: *****

           call strlength(pre_DBP,LEN(pre_DBP),l)
           isys=SYSTEM('mv '//pre_DBP(1:l)//' '//pre_DBP(1:l)//
     +          '.trash'//char(0))                                     !move old files out of the way
           OPEN (35,FILE=pre_DBP(1:l)//'.trash',FORM='FORMATTED',
     +          STATUS='OLD',ACCESS='DIRECT',RECL=iDBP_RECL)
           OPEN (36,FILE=pre_DBP,FORM='FORMATTED',
     +          STATUS='NEW',ACCESS='DIRECT',RECL=iDBP_RECL)

           total_num_cells_DB = 0
           do i = 1, MAX_NUM_CODES                              !find out how many cells were included in the original QDT file - THIS MIGHT NOT BE RIGHT!!!  FIND THE GREATEST # INSTEAD OF INCREMENTING?!?!?!?
              if(included_DB(i).ne.0)then                       ! (ie., the cells now in the db files)
                 total_num_cells_DB = total_num_cells_DB + 1
              end if
           end do
           do REF = 1, MAX_NUM_CODES
              do TAR = REF+1, MAX_NUM_CODES
                 if((included_QDT(REF).ne.0).and.               !if REF and TAR are included in the QDTSAV file...
     +                (included_QDT(TAR).ne.0))then

                    prim = ' '                                  !initialize variables
                    sec = ' '
                    loc1 = ' '
                    loc2 = ' '
                    rescom = ' '
                    czk = ' '
                    cprobk = ' '
                    cdet = ' '
                    cvis = ' '
                    czlat = ' '
                    chalfwd = ' '
                    statcomm = ' '
                    stats_bw = ' '
                    CELL_NAME_REF = CELL_NAMES(IDs(REF))
                    CELL_NAME_TAR = CELL_NAMES(IDs(TAR))
                    IDs_REF =  REF
                    IDs_TAR =  TAR
                    ios = 0

                    Q_pos_QDT = 0
                    if(included_QDT(REF).gt.1)then
                       do i = 1, included_QDT(REF)-1 
                          Q_pos_QDT = Q_pos_QDT + 
     +                         (total_num_cells_QDT-i) 
                       end do
                    end if
                    Q_pos_QDT = Q_pos_QDT +                            !calculate the position of this pair in the pair queue
     +                   (included_QDT(TAR) - included_QDT(REF))       !  in the QDT file
                    rel_loc = (Q_pos_QDT*4)-3                          !re-calculate the record number in the preDBP file
                    hist_number_QDT = 3+(7*total_num_cells_QDT)+rel_loc
c                    print '(''hist_number_QDT(pair) = '',I5)',
c     +                   hist_number_QDT
                    
                    if((included_DB(REF).ne.0).and.                    !if REF and TAR are also included in the DBSAV file,
     +                   (included_DB(TAR).ne.0))then                  !  retrieve the existing data
                       Q_pos_DB = 0
                       if(included_DB(REF).gt.1)then
                          do i = 1, included_DB(REF)-1 
                             Q_pos_DB = Q_pos_DB+(total_num_cells_DB-i) 
                          end do
                       end if
                       Q_pos_DB = Q_pos_DB +                           !calculate the position of this pair in the pair queue
     +                      (included_DB(TAR) - included_DB(REF))      !  in the DB file
c                       print '(''REF: '',I3,''; TAR: '',I3,
c     +                      ''; Q_pos_QDT = '',I,
c     +                      ''; Q_pos_DB = '',I)',
c     +                      REF,TAR,Q_pos_QDT,Q_pos_DB
                       call read_preDBP(35,Q_pos_DB,ios,recording,     !read existing DB info from old preDBP file
     +                      date_exp, 
     +                      BDT_FILE_30,CELL_NAME_REF,IDs_REF,
     +                      CELL_NAME_TAR,IDs_TAR,
     +                      prim,loc1,sec,loc2,rescom,czk,cprobk,
     +                      cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +                      ihist_num)
c                       print '(''ios = '',I,''; REF = '',I3,
c     +                      ''; IDs_REF = '',I3,''// TAR = '',I3,
c     +                      ''; IDs_TAR = '',I3)',ios,REF,IDs_REF,
c     +                      TAR,IDs_TAR
c                       if(ios.eq.0)then                               !record read from old preDBP without error, so re-write
                       if((ios.eq.0).and.(REF.eq.IDs_REF).and.
     +                      (TAR.eq.IDs_TAR))then                      !record read from old preDBP without error, so re-write
                          BDT_FILE_30 = BDT_FILE
c                          print '(''call write: hist_num = '',I5)',
c     +                         hist_number_QDT
                          call write_preDBP(36,Q_pos_QDT,ios,recording,
     +                      date_exp,                                  !write the same info to a new place in the new preDBP file
     +                      BDT_FILE_30,CELL_NAME_REF,IDs_REF,
     +                      CELL_NAME_TAR,IDs_TAR,
     +                      prim,loc1,sec,loc2,rescom,czk,cprobk,
     +                      cdet,cvis,czlat,chalfwd,statcomm,stats_bw,
     +                      hist_number_QDT)

                          if(prim.ne.' ')then
                             analyzed_pairs(IDs(REF),IDs(TAR)) = 2     !this pair has been analyzed
                          else
                             analyzed_pairs(IDs(REF),IDs(TAR)) = 1     !this pair has NOT been analyzed
                          end if
                          cycle
                       else
c                          print '(''error reading record for pair ('',
c     +                         I3,'','',I3,
c     +                         ''): ios = '',I5,
c     +                         '', but it''''s OK - I''''ll just '',
c     +                         ''write a place-holder record'')',
c     +                         REF,TAR,ios
c                          read (*,'(A)')
                       end if
                    end if

*     ***** If get to this point, either the pair is an "added" one, or it does exist in the old preDBP file but its record  *****
*     *****  wasn't read without an error, so write a "place-holder" record to the new preDBP file:   *****

                    analyzed_pairs(IDs(REF),IDs(TAR)) = 1              !tag this pair as valid, but not-yet-analyzed
                    BDT_FILE_30 = BDT_FILE

                    call write_preDBP(36,Q_pos_QDT,ios,recording,
     +                   date_exp,                                     !write a place-holder record to the new preDBP file
     +                   BDT_FILE_30,blank4,IDs_REF,
     +                   blank4,IDs_TAR,
     +                   blank20,blank15,blank20,blank15,blank40,
     +                   blank8,blank8,
     +                   blank8,blank8,blank8,blank8,blank40,blank8,
     +                   hist_number_QDT)

                 end if
              end do
           end do
           
           close (unit=35)
           close (unit=36)

           call strlength(pre_DBG,LEN(pre_DBG),l)
           isys=SYSTEM('rm '//pre_DBG(1:l)//'.trash'//char(0))     !remove old files 
           call strlength(pre_DBP,LEN(pre_DBP),l)
           isys=SYSTEM('rm '//pre_DBP(1:l)//'.trash'//char(0))     !remove old files 

           OPEN (9,FILE=pre_DBG,FORM='FORMATTED',
     +          ACCESS='DIRECT',RECL=DBG_RECL)        !re-open the preDBG file
           OPEN (10,FILE=pre_DBP,FORM='FORMATTED',
     +          ACCESS='DIRECT',RECL=iDBP_RECL)       !re-open the preDBP file

        end if CHECK_FOR_CONVERT

        return
        end
      end module mod_compare_and_convert
