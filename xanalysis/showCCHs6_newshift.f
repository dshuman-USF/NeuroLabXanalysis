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


      module mod_showCCHs6_newshift
      contains
*       SUBROUTINE TO CREATE HISTOGRAMS in X11 window
*       date of last revision = 21-Sep-2006     lss
*
*       filename = showCCHs6_newshift.f
*
*       sep-2000
*        option added to allow user to view CCH window without having
*               to make database entries
*
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*       *** INDIRECT POINTERS are now used to access unit dat
*       ***  within the following arrays:
*       ***     DATA_ARRAY, ITAL, card_type, coef, coefnum, TALLY_NORM,
*       ***     ETA2_*, zmodsig_*, zmodsig2_*, analyzed_cells, analyzed_pairs,
*       ***     resp_type, card, CELL_NAMES
*       ***     [these arrays are now dimensioned to MAX_NUM_CHAN;
*       ***      the location of a cell's information (user-assigned ID code = i)
*       ***      within these arrays is array(IDs(i))]
*       ***
*       *** DIRECT POINTERS are used with these arrays: exclude, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       with v3 - contains these subroutines: ENTER_ANALYSIS_RESULTS, 
*
*       new_plot incorporated 01-apr-99
*
*       LINK with x2000_v* code 
*
*
      SUBROUTINE showCCHs (fildes,IDs,mouse,mode,ihs,iohist_1,ilhist,
     +     irhist,ICN,BINW,izz,CTH_REF,CTH_TAR,phrenic_hist,cardiac_REF
     +     ,cardiac_TAR,cardiac_hist,ETA2_REF,zmodsig_REF,zmodsig2_REF
     +     ,card_type_REF,ETA2_TAR,zmodsig_TAR,zmodsig2_TAR
     +     ,card_type_TAR,REFcode,TARcode,REFname,TARname,REFtype
     +     ,TARtype,button_choice,REC_NUM_X,rec_num_control_1,auto_mode
     +     ,speed,coefval,click,prim,sec,loc1,loc2,rescom,Q_pos,czk
     +     ,cprobk,cdet,cvis,chalfwd,czlat,ZK,stats_bw,Probk,det,vis
     +     ,ZLAT,HALFWD,BDT_FILE,QDT_FILENAME,date,recording,NHW,id
     +     ,flag,WINDOW1,WINDOW2,show_controlMAIN,c_format,statcomm
     +     ,qdt_files,current_qdt,included,total_num_cells,DELTA2_REF
     +     ,DELTA2_TAR,tedfactor_REF,tedfactor_TAR,per_text_abbrev
     +     ,per_results_REF,per_results_TAR,iohist_avg
     +     ,rec_num_control_2, rec_num_control_avg, single_shift
     +     ,avg_shift,show_single_shift,show_avg_shift, iohist_2
     +     ,sp_per_cycle,num_acc_cycles,show_2_sd,show_3_sd,
     +     show_conf_lim,sig_only,ITAL,cardiac_pls,flats_only,
     +     total_num_qdts,scaledMAIN,REF_coords,TAR_coords,exp_name,
     +     NORM_OFFSET,mean_E,bw_n,show_pulses)
      use mod_acumsum_nomad_2
      use mod_clear_routines
      use mod_locate_region
      use mod_locate_region2
      use mod_miscellaneous_subroutines
      use mod_new_draw_button
      use mod_new_plot
      use mod_plot_spikes_per_cycle
      use mod_print_and_write_routines
      use mod_stats6
      use mod_test_draw_button
*     
      INCLUDE 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
*     
      DIMENSION ICN(6),ihs(101),ilhist(101),
     +     irhist(101),ICN_control(6),
     +     iohist_1(101),iohist_2(101),
     +     iohist_avg(101),IDs(MAX_NUM_CODES)
*     
      integer*4 mouse
      integer REFcode,TARcode,REC_NUM_X,CCH(101),flag,
     +     CTH_REF(101),CTH_TAR(101),phrenic_hist(101),
     +     cardiac_REF(101),cardiac_TAR(101),cardiac_hist(101),
     +     regionAUTO,regionSPEED,regionSINGLE,regionDIF,
     +     regionMAIN,regionSTATS,Q_pos,
     +     rec_num,ICN_original(6),rec_num_ACH,
     +     subtract,subtract_1,CCH_scaled(101),
     +     CCH_mirror(101),CONTROL_mirror(101),
     +     CONTROL_cl_mirror(101),CONTROL_cl_scaled(101),
     +     CCH_original(101),CONTROL(101),rec_num_control,
     +     rec_num_control_1,rec_num_control_2,rec_num_control_avg,
     +     DIFF(101),scaled_hist(101),
     +     CONTROL_original(101),CONTROL_scaled(101),current_qdt,
     +     included(MAX_NUM_CODES),total_num_cells,
     +     CCHs(4,101),CONTROLs_1(4,101),CONTROLs_AVG(4,101),
     +     CONTROL_CCHs(4,101),CONTROLs_2(4,101),
     +     CONTROLs_cl(4,101),
     +     REF_ACHs(4,101),TAR_ACHs(4,101),CONTROL_cl(101),
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     CONTROL_1(101),CONTROL_2(101),CONTROL_AVG(101),
     +     rec_num_control_single,rec_num_control_two,
     +     rec_num_control_average,ITAL(MAX_NUM_CHAN),
     +     cardiac_pls,big_ACH(101),total_num_qdts
     +     ,temp1(101),temp2(101),subtract_cl
*     
      character*120 title,textstring,label1,title_original,prim_txt,
     +     sec_txt,comm_txt,prim_txt_stats
      character*9 WINDOW1
      character*5 WINDOW2
      character*5 DELTA2_REF,DELTA2_TAR
      character*(*) QDT_FILENAME
      character*(*) exp_name
      character*(*) REF_coords(3),TAR_coords(3)
      character*(*) BDT_FILE
      character*(*) c_format,qdt_files(MAX_NUM_QDTS)
*     
      character*40 rescom,statcomm
      character*38 label5
*     
      character*30 prev_pair,dir_acc,main_menu,stat_anal,
     +     next_pair,auto_stop,auto_view,
     +     slow,fast,medium,return_box,perimeter,qsum_text,
     +     stats,
     +     label2,label3,label4,clickable,user_defined1,
     +     user_defined2,go,view_all,default,
     +     view_1,view_2,view_3,blank30
*     
      character*15 REFtype,TARtype,REFtype1,TARtype1,button_choice
      character*20 prim,sec
      character*12 bwtext2,REFname1,TARname1,BINWtexts(4),
     +     per_text_abbrev(MAX_PERTURB)
      save BINWtexts
      character*11 date
      character*9 per_results_REF(MAX_PERTURB),
     +     per_results_TAR(MAX_PERTURB)
      character*8 czk,cdet,cprobk,cvis,
     +     chalfwd, czlat,stats_bw
      character*15 loc1,loc2,screen
      character*6 bwtext
      character*5 ETA2_REF,ETA2_TAR
      character*4 REFcodetext,TARcodetext,REFcodetext1,TARcodetext1
      character*4 REFname,TARname
      character*3 CCHtype,ACH,CTH,DIF,SPC,ACHtype,
     +     card_type_REF,card_type_TAR,
     +     zmodsig_REF,zmodsig2_REF,zmodsig_TAR,zmodsig2_TAR
      character*2 recording,mode,
     +     tedfactor_REF,tedfactor_TAR
c     +             tsfs_REF(MAX_PERTURB),
c     +             tsfs_TAR(MAX_PERTURB)
      character*1 auto_mode,click,info,show_control,task,
     +     show_controlMAIN,show_conf_lim,
     +     show_control_STATS,show_control_DIFF,mirror,
     +     scaledup,
     +     single_shift,avg_shift,show_single_shift,
     +     show_avg_shift,show_2_sd,show_3_sd,sig_only,
     +     flats_only,scaledMAIN,show_pulses,ign1
*     
      real x,y,x_dc,y_dc,z_dc,coefval,NORM_OFFSET,mean_E,bw_n
      integer*4 valid,keystroke,ICN0(6)

      logical cth_cch

      save CONTROLs_cl
      save CONTROLs_1
      save CONTROLs_2
      save REF_ACHs
      save TAR_ACHs
      save CCHs

      if(.false.)print *,single_shift !suppress unused variable warning
      if(.false.)print *,fildes2 !suppress unused variable warning

      cth_cch = (qdt_file_version.EQ.'8')

*     
*     
*     *****  set Starbase screen parameters  ******************************
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call mapping_mode(fildes,1)
      call view_port(fildes,0.,0.,1.,1.)
      call view_window(fildes,1.,1.,1500.,800.)
c     *********************************************************************
*     *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
      call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
      call ps_mapping_mode(1)
      call ps_view_port(0.,0.,1.,1.)
      call ps_view_window(1.,1.,1500.,800.)
      call ps_geometry(18.,18.,792.-18.,612.-18.)
      call fill_color(fildes,0.,0.,0.)               
*     *********************************************************************

      ICN0 = 0
      timeout = 2e9         !timeout for waiting for a mouse click = 63 years
      max_bin = 0
      min_bin = 0
      info = 'n'
      mirror = 'n'
      call sc_mirror (0)
      scaledup = 'n'
      button_choice=' '
      CCHtype = 'CCH'
      ACHtype = 'ACH'
      ACH = 'ACH' 
      CTH = 'CTH'
      DIF='DIF'
      SPC='SPC'
      REFname1=' '
      REFtype1=' '
      TARname1=' '
      TARtype1=' '
      blank30=' '
      write (REFcodetext,'(I4)') REFcode
      write (TARcodetext,'(I4)') TARcode
      title=' '
      title = REFname//' > '//TARname//
     +     '    ('//REFcodetext//' > '//TARcodetext//')'
      label1 = 'CCHs:'
      label2 = 'ACHs:'
      label3='('//REFcodetext//','//TARcodetext//')'
      call remove_all_blanks(label3,LEN(label3))
      label4 = 'CTHs:'
      write (bwtext,'(f4.1)') BINW
      bwtext2 = ' '
      bwtext2 = bwtext//' ms.'
      if(.false.)then           !define (0,0) of each CCH
      else if(izz.eq.1)then
         call clear_all(fildes)
         call character_height(fildes,.021)
         call character_width(fildes,.007)
         call strlength(qdt_files(current_qdt),
     +        LEN(qdt_files(current_qdt)),l)
         textstring=qdt_files(current_qdt)(1:l)
         call strlength(textstring,LEN(textstring),l)
         call text2d(fildes,1.,760.,
     +        textstring(1:l)//char(0),ANNOTATION_TEXT,0)
         textstring=BDT_FILE
         call strlength(textstring,LEN(textstring),l)
         call text2d(fildes,1.,775.,textstring(1:l)
     +        //char(0),ANNOTATION_TEXT,0)
         call character_height(fildes,.075)
         call character_width(fildes,.025)
         call strlength(title,LEN(title),l)
         call text2d (fildes,150.,760.,title(1:l)
     +        //char(0),ANNOTATION_TEXT,0) !'REF>TAR'
         call character_height(fildes,.030)
         call character_width(fildes,.010)
         call strlength(label3,LEN(label3),l)
         call text2d (fildes,5.,460.,label3(1:l)
     +        //char(0),ANNOTATION_TEXT,0) !'(ref,tar)'
         call text2d(fildes,5.,405.,'(1st 50 bins shown)'//char(0),
     +        ANNOTATION_TEXT,0)
         call text2d(fildes,5.,435.,'(click to'//char(0),
     +        ANNOTATION_TEXT,0)
         call text2d(fildes,10.,420.,'enlarge)'//char(0),
     +        ANNOTATION_TEXT,0)
         call character_height(fildes,.045)
         call character_width(fildes,.015)
         call strlength(label1,LEN(label1),l)
         call text2d (fildes,1.,710.,label1(1:l)
     +        //char(0),ANNOTATION_TEXT,0) !'CCHs:'
         call strlength(label2,LEN(label2),l)
         call text2d (fildes,1.,475.,label2(1:l)
     +        //char(0),ANNOTATION_TEXT,0) !'ACHs:'
         call strlength(label3,LEN(label3),l)
         call text2d (fildes,1.,300.,label4(1:l)
     +        //char(0),ANNOTATION_TEXT,0) !'CTHs:'
         label5 = '[click a'
         call character_height(fildes,.030)
         call character_width(fildes,.010)
         call text2d(fildes,3.,640.,label5
     +        //char(0),ANNOTATION_TEXT,0)
         label5 = ' CCH for '
         call text2d(fildes,3.,625.,label5
     +        //char(0),ANNOTATION_TEXT,0)
         label5 = ' more'
         call text2d(fildes,3.,610.,label5
     +        //char(0),ANNOTATION_TEXT,0)
         label5 = ' options]'
         call text2d(fildes,3.,595.,label5
     +        //char(0),ANNOTATION_TEXT,0)
         if(scaledMAIN.eq.'n')then
            label1 = 'scale up'//char(0)
            label2 = '(original)'//char(0)
         elseif(scaledMAIN.eq.'y')then
            label1 = 'original'//char(0)
            label2 = '(scaled up)'//char(0)
         end if
         call text2d(fildes,5.,675.,label2,ANNOTATION_TEXT,0)
         call character_height(fildes,.033)
         call character_width(fildes,.011)
         call rectangle(fildes,3.,560.,90.,590.)
         call strlength(label1,LEN(label1),l)
         call text2d(fildes,7.,565.,label1(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
         call character_height(fildes,.030)
         call character_width(fildes,.010)
         textstring = ' '
         if(show_single_shift.eq.'y')then
            if(show_2_sd.eq.'y')then
               textstring='(single shift-control +/- 2sd displayed)'
            else if(show_3_sd.eq.'y')then
               textstring='(single shift-control +/- 3sd displayed)'
            else
               textstring='(single shift-control displayed)'
            end if
         else if(show_avg_shift.eq.'y')then
            if (cth_cch)then
               if(show_2_sd.eq.'y')then
                  textstring='(cth control +/- 2sd displayed)'
               else if(show_3_sd.eq.'y')then
                  textstring='(cth control +/- 3sd displayed)'
               else
                  textstring='(cth control displayed)'
               end if
            else
               if(show_2_sd.eq.'y')then
                  textstring=
     +                 '(averaged shift-control +/- 2sd displayed)'
               else if(show_3_sd.eq.'y')then
                  textstring=
     +                 '(averaged shift-control +/- 3sd displayed)'
               else
                  textstring='(averaged shift-control displayed)'
               end if
            end if
         end if
         call text2d(fildes,10.,530.,textstring//char(0),
     +        ANNOTATION_TEXT,0)
         call sc_label (fildes, 400., 530.,show_single_shift
     +        ,show_avg_shift,show_2_sd,show_3_sd)
         call character_height(fildes,.036)
         call character_width(fildes,.012)
         call character_height(fildes,.030)
         call character_width(fildes,.010)
         label4 = 'resp:'
         call text2d (fildes,15.,375.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         if (CTH_REF(101).gt.0) then
            call text2d(fildes,5.,360.,'(ph. norm.)'//char(0),
     +           ANNOTATION_TEXT,0)
         else
            call text2d(fildes,10.,360.,'(norm.)'//char(0),
     +           ANNOTATION_TEXT,0)
         end if
         label4 = 'cardiac:'
         call text2d (fildes,860.,375.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         call character_height(fildes,.021)
         call character_width(fildes,.007)
         label4 = 'IDcode (type):'
         call text2d (fildes,1.,210.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         label4 = 'ANOVA/BINARY:'
         call text2d (fildes,1.,190.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         label4 = 'ETA2:'
         call text2d (fildes,1.,175.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         x=100
         y=550
         do i = 1, NUM_BINS
            do j = 1, 4
               CCHs(j,i) = 0    !initialize arrays that will hold histogram data
               CONTROLs_1(j,i) = 0
               CONTROLs_2(j,i) = 0
               CONTROLs_AVG(j,i) = 0
               CONTROLs_cl(j,i) = 0
               REF_ACHs(j,i) = 0
               TAR_ACHs(j,i) = 0
            end do
         end do
         BINWtexts = ' '
      else if(izz.eq.2)then
         x=450
         y=550
      else if(izz.eq.3)then
         x=800
         y=550
      else if(izz.eq.4)then
         x=1150
         y=550
      end if

      do i = 1, NUM_BINS
         CCHs (izz,i)= ihs(i)   !save histogram data for printing/writing
         CONTROLs_1(izz,i) = iohist_1(i)
         CONTROLs_AVG(izz,i) = iohist_avg(i)
         CONTROLs_2(izz,i) = iohist_2(i)
         REF_ACHs(izz,i) = ilhist(i)
         TAR_ACHs(izz,i) = irhist(i)
      end do
      BINWtexts(izz) = bwtext2

      if(mode.ne.'vt')then      !if viewing TQDTs --> "view-only" mode - no db info displayed
         if(prim.eq.' ')then
            prim_txt = ' '
            prim_txt_stats=' '
         else if(prim.eq.'Flat')then
            prim_txt = 'FLAT'
            prim_txt_stats=' '
         else if(prim.eq.'M P & T')then
            call strlength(loc1,LEN(loc1),i)
            prim_txt = 'M P & T '//loc1(1:i)//char(0)
         else 
            call strlength(loc1,LEN(loc1),i)
            call strlength(prim,LEN(prim),j)
            call remove_all_blanks(cdet,LEN(cdet))
            call strlength(cdet,LEN(cdet),n)
            call strlength(statcomm,LEN(statcomm),k)
            if(cdet.ne.' ')then
c               prim_txt='Prim: '//loc1(1:i)//' '//
c     +              prim(1:j)//' (d.i. = '//
c     +              cdet(1:n)//'; '//
c     +              statcomm(1:k)//')'//char(0)
               prim_txt='Prim: '//loc1(1:i)//' '//
     +              prim(1:j)
               prim_txt_stats=' (d.i. = '//
     +              cdet(1:n)//'; '//
     +              statcomm(1:k)//')'
            else
               prim_txt='Prim: '//loc1(1:i)//' '//
     +              prim(1:j)//' (no statistics calculated)'
            end if
         end if
         if((sec.ne.'None').and.(sec.ne.' '))then
            call strlength(loc2,LEN(loc2),i)
            call strlength(sec,LEN(sec),j)
            sec_txt = 'Sec: '//loc2(1:i)//' '//
     +           sec(1:j)//char(0)
         else
            sec_txt = ' '
         end if
         if(rescom.ne.' ')then
            comm_txt = 'Comment: '//rescom//char(0)
         else
            comm_txt = ' '
         end if
      else
         prim_txt = ' '
         prim_txt_stats = ' '
      end if
*     
*     
*     **************************************************
*     *****   plot the histograms for this pair:   *****
*     **************************************************
*     
      call character_height(fildes,.045)
      call character_width(fildes,.015)
      call text2d(fildes,(x+100.),(y-155.),bwtext2//char(0),
     +     ANNOTATION_TEXT,0)  
      if((show_single_shift.eq.'y').or.cth_cch)then
         CONTROL = iohist_1     !the "shift by 2" SP is used for calculation of confidence
         CONTROL_cl = iohist_2  !  limits for the "shift by 1" CONTROL
c     rec_num_control = rec_num_control_single
      else if(show_avg_shift.eq.'y')then !in the case of the average shift, use the "shifted by 1" to
         CONTROL = iohist_avg   !  compute the confidence limits
         CONTROL_cl = iohist_1
c     rec_num_control = rec_num_control_average
      end if
      if(scaledMAIN.eq.'y')then
         subtract = 75
         min_bin = ihs(1)
         max_bin = ihs(1)
         do i = 2,101
            min_bin = MIN0(min_bin,ihs(i))   
            max_bin = MAX0(max_bin,ihs(i))
         end do
         subtract_1 = min_bin*(subtract/100.)
         call sc_subtract (subtract_1)
         do i = 1, 101
            CCH_scaled(i) = ihs(i) - subtract_1
            CONTROL_scaled(i) = CONTROL(i) - subtract_1
            subtract_cl = subtract_1
            if (cth_cch) subtract_cl = 0
            CONTROL_cl_scaled(i) = CONTROL_cl(i) - subtract_cl
         end do
         call new_plot(fildes,101,CCH_scaled,CONTROL_scaled,scaled_hist,
     +        CCHtype,x,y,200.,303.,info,'','','','',0,0,          
     +        show_controlMAIN,0,'',ICN0,0,show_conf_lim,'','', !plot the CCH  ...
     +        CONTROL_cl_scaled,show_2_sd,show_3_sd,scaledup,  
     +        min_bin,rmagnify,IDs,ITAL)                       
      else
         call sc_subtract (0)
         call new_plot(fildes,101,ihs,CONTROL,scaled_hist,CCHtype,x,y,  
     +        200.,303.,info,'','','','',0,0,show_controlMAIN,0,'', !plot the CCH  ...
     +        ICN0,0,show_conf_lim,'','',CONTROL_cl,show_2_sd,      
     +        show_3_sd,scaledup,min_bin,rmagnify,IDs,ITAL)    
      end if  

      call new_plot(fildes,50,ilhist,ilhist,scaled_hist,ACH,(x-2.), !the reference ACH and ...
     +     (y-125.),100.,150.,info,'','','','',0,0,'',0,'',ICN0,
     +     0,'','','', [integer::],'','',scaledup,min_bin,rmagnify,IDs,  ! (show the first 50 bins)
     +     ITAL)

      call new_plot(fildes,50,irhist,irhist,scaled_hist,ACH,(x+152.+2.), !the target ACH and ...
     +     (y-125.),100.,150.,info,'','','','',0,0,'',0,'',ICN0,
     +     0,'','','',[integer::],'','',scaledup,min_bin,rmagnify,IDs,   ! (show the first 50 bins)
     +     ITAL)

      if(izz.eq.4)then
         show_control = 'y'
         call new_plot(fildes,100,CTH_REF,phrenic_hist,scaled_hist,CTH,   
     +        100.,230.,150.,225.,info,'','','','',0,0,show_control, !the reference CTH and ...
     +        0,'',ICN0,0,'','','',[integer::],'','',scaledup,
     +        min_bin,rmagnify,IDs, ITAL)                                              
         call new_plot(fildes,100,CTH_TAR,phrenic_hist,scaled_hist,CTH,   
     +        350.,230.,150.,225.,info,'','','','',0,0,show_control, !the target CTH
     +        0,'',ICN0,0,'','','',[integer::],'','',scaledup,
     +        min_bin,rmagnify,IDs, ITAL)                                              
         call new_plot(fildes,100,CTH_REF,CTH_TAR,scaled_hist,CTH,600.,   
     +        230.,150.,225.,info,'','','','',0,0,show_control,0,'',
     +        ICN0, 0,'','','',[integer::],'','',scaledup,
     +        min_bin, rmagnify,IDs,ITAL) !the superimposed CTHs

         call character_height(fildes,.033)
         call character_width(fildes,.011)
         label1 = ' pulses'//char(0)
         call rectangle(fildes,3.,250.,90.,280.)
         call text2d(fildes,7.,255.,label1,ANNOTATION_TEXT,0)
         if(show_pulses.ne.'n')then
            call character_height(fildes,.030)
            call character_width(fildes,.010)
            call line_color(fildes,0.,0.,1.) !blue
            call text_color(fildes,0.,0.,1.)
            call line_type(fildes,DOT)
            call move2d(fildes,100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.)
            call draw2d(fildes,100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+3.)
            call text2d(fildes,100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+5.,'E'//char(0),ANNOTATION_TEXT,0)

            call move2d(fildes,350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.)
            call draw2d(fildes,350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+3.)
            call text2d(fildes,350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+5.,'E'//char(0),ANNOTATION_TEXT,0)

            call move2d(fildes,600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.)
            call draw2d(fildes,600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+3.)
            call text2d(fildes,600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))
     +           ), 230.+150.+5.,'E'//char(0),ANNOTATION_TEXT,0)
            if(mean_E.ne.0.0)then
               call move2d(fildes,100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.)
               call draw2d(fildes,100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+3.)
               call text2d(fildes,100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+5.,'I'//char(0),
     +              ANNOTATION_TEXT, 0)

               call move2d(fildes,350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.)
               call draw2d(fildes,350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+3.)
               call text2d(fildes,350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+5.,'I'//char(0),
     +              ANNOTATION_TEXT,0)

               call move2d(fildes,600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.)
               call draw2d(fildes,600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+3.)
               call text2d(fildes,600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +              (100.*bw_n))),230.+150.+5.,'I'//char(0),
     +              ANNOTATION_TEXT,0)
            end if
            call line_color(fildes,0.,0.,0.)
            call text_color(fildes,0.,0.,0.)
            call line_type(fildes,SOLID)
         end if

         if(cardiac_pls.ne.0)then
            call new_plot(fildes,100,cardiac_REF,cardiac_hist,
     +           scaled_hist,CTH,950.,230.,150.,225.,info,'','','','',0,             !the reference cardiac CCH and ...
     +           0,show_control,0,'',ICN0,0,'','','',[integer::],
     +           '','',scaledup,min_bin, rmagnify,IDs,ITAL)                              

            call new_plot(fildes,100,cardiac_TAR,cardiac_hist,
     +           scaled_hist,CTH,1200.,230.,150.,225.,info,'','','','',           !the target cardiac CCH
     +           0,0,show_control,0,'',ICN0,0,'','','',
     +           [integer::],'','',scaledup,min_bin, rmagnify,IDs,ITAL)                              
         else
            call character_height(fildes,.030)
            call character_width(fildes,.010)
            textstring = 'CARDIAC CCHs NOT CALCULATED'
            call strlength(textstring,LEN(textstring),l)
            call text2d(fildes,1000.,330.,textstring(1:l)//char(0),
     +           ANNOTATION_TEXT,0)
         end if


         call character_height(fildes,.030)
         call character_width(fildes,.010)
         call remove_all_blanks(REFcodetext,LEN(REFcodetext))
         call remove_all_blanks(REFtype,LEN(REFtype))
         call strlength(REFcodetext,LEN(REFcodetext),l)
         call strlength(REFtype,LEN(REFtype),m)
         label4 = REFcodetext(1:l)//' ('//REFtype(1:m)//')'
         call text2d (fildes,150.,210.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         call upper_case(zmodsig_REF,LEN(zmodsig_REF))
         call upper_case(zmodsig2_REF,LEN(zmodsig2_REF))
         label4=zmodsig_REF//'/ '//zmodsig2_REF
         call strlength(label4,LEN(label4),l)
         call text2d(fildes,200.,190.,label4(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
         label4=ETA2_REF
         call text2d(fildes,150.,175.,label4
     +        //char(0),ANNOTATION_TEXT,0)

         call remove_all_blanks(TARcodetext,LEN(TARcodetext))
         call remove_all_blanks(TARtype,LEN(TARtype))
         call strlength(TARcodetext,LEN(TARcodetext),l)
         call strlength(TARtype,LEN(TARtype),m)
         label4 = TARcodetext(1:l)//' ('//TARtype(1:m)//')'
         call strlength(label4,LEN(label4),l)
         call text2d (fildes,400.,210.,label4(1:l)
     +        //char(0),ANNOTATION_TEXT,0)
         call upper_case(zmodsig_TAR,LEN(zmodsig_TAR))
         call upper_case(zmodsig2_TAR,LEN(zmodsig2_TAR))
         label4=zmodsig_TAR//'/ '//zmodsig2_TAR
         call text2d(fildes,450.,190.,label4//char(0),
     +        ANNOTATION_TEXT,0)
         label4=ETA2_TAR
         call text2d(fildes,400.,175.,label4
     +        //char(0),ANNOTATION_TEXT,0)

         call strlength(REFcodetext,LEN(REFcodetext),l)
         call strlength(TARcodetext,LEN(TARcodetext),m)
         label4 = REFcodetext(1:l)//' & '//TARcodetext(1:m)
         call text2d(fildes,710.,210.,label4
     +        //char(0),ANNOTATION_TEXT,0)

         if(cardiac_pls.ne.0)then
            label4=REFcodetext
            call text2d (fildes,1060.,210.,label4
     +           //char(0),ANNOTATION_TEXT,0)
            label4=TARcodetext
            call text2d (fildes,1310.,210.,label4
     +           //char(0),ANNOTATION_TEXT,0)

            label4 = 'DELTA2:'
            call text2d(fildes,860.,190.,label4
     +           //char(0),ANNOTATION_TEXT,0)
            label4 = '  (value; t.f.)'
            call text2d(fildes,860.,175.,label4
     +           //char(0),ANNOTATION_TEXT,0)
            if(card_type_REF.ne.' ')then
               call upper_case(card_type_REF,LEN(card_type_REF))
               call text2d(fildes,1050.,190.,card_type_REF//char(0),
     +              ANNOTATION_TEXT,0) !cardiac-modulated acc'ding to Ted's DELTA**2 test?
               label4='('//DELTA2_REF//', '//tedfactor_REF//')'
               call remove_all_blanks(label4,LEN(label4))
               call strlength(label4,LEN(label4),l)
               call text2d(fildes,1040.,175.,label4(1:l)//char(0),
     +              ANNOTATION_TEXT,0)
            end if
            if(card_type_TAR.ne.' ')then
               call upper_case(card_type_TAR,LEN(card_type_TAR))
               call text2d(fildes,1300.,190.,card_type_TAR//char(0),
     +              ANNOTATION_TEXT,0)
               label4='('//DELTA2_TAR//', '//tedfactor_TAR//')'
               call remove_all_blanks(label4,LEN(label4))
               call strlength(label4,LEN(label4),l)
               call text2d(fildes,1290.,175.,label4(1:l)//char(0),
     +              ANNOTATION_TEXT,0)
            end if
         end if


         if(.true.)goto 300     !skip the next section of code for now
*
         label4='phrenic'
         call character_height(fildes,.036)
         call character_width(fildes,.012)
         call text2d(fildes,100.,300.,REFname
     +        //char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,100.,275.,REFtype
     +        //char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,100.,225.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         call move2d(fildes,50.,310.)
         call draw2d(fildes,90.,310.)
         call move2d(fildes,50.,230.)
         call line_type(fildes,2)
         call draw2d(fildes,90.,230.)
         call line_type(fildes,0)
         call text2d(fildes,600.,300.,TARname
     +        //char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,600.,275.,TARtype
     +        //char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,600.,225.,label4
     +        //char(0),ANNOTATION_TEXT,0)
         call move2d(fildes,550.,310.)
         call draw2d(fildes,590.,310.)
         call move2d(fildes,550.,230.)
         call line_type(fildes,2)
         call draw2d(fildes,590.,230.)
         call line_type(fildes,0)
         call text2d(fildes,1100.,300.,REFname1
     +        //char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,1100.,225.,TARname1
     +        //char(0),ANNOTATION_TEXT,0)
         call move2d(fildes,1050.,310.)
         call draw2d(fildes,1090.,310.)
         call move2d(fildes,1050.,230.)
         call line_type(fildes,2)
         call draw2d(fildes,1090.,230.)
         call line_type(fildes,0)

*       *****  Get ready to draw buttons on bottom of screen *****
*
 300     perimeter = 'perimeter'
         prev_pair = 'PREVIOUS PAIR' !define text strings to label buttons
         dir_acc = 'DIRECT ACCESS'
         main_menu = 'RETURN TO MAIN MENU'
         stat_anal = 'STATISTICAL ANALYSIS'
         next_pair = 'NEXT PAIR'
         auto_stop = 'STOP'
         auto_view = 'AUTOVIEW'
         user_defined1 = 'User-'
         user_defined2 = 'defined'
         slow = '10'
         medium = '5'
         fast = '3        seconds'
         go = 'GO'
         default = '(default)'
         view_all = '   ALL'
         return_box = 'RETURN'
         qsum_text = 'Q_SUM'
         stats = 'STATISTICS'
         clickable = 'CLICKABLE'

         do i = 1, MAX_NUM_QDTS
            x0 = (float(i)*125.)+700. !calculate the x-coordinates for the "switch boxes"
            x1 = x0 + 100.
            if((current_qdt.ne.i).and.(qdt_files(i).ne.' '))then
               call character_height(fildes,.030)
               call character_width(fildes,.010)
               label1='click on one of these to see '
               call text2d(fildes,525.,785.,label1//char(0),
     +              ANNOTATION_TEXT,0)
               label1='  another set of CCHs -->'
               call text2d(fildes,525.,770.,label1//char(0),
     +              ANNOTATION_TEXT,0)
               call strlength(qdt_files(i),LEN(qdt_files(i)),l_q)
               call strlength(exp_name,LEN(exp_name),l_e)
               if(qdt_files(i)(1:l_e).eq.exp_name(1:l_e))then
                  label1 = qdt_files(i)(l_e+1:l_q-4)
                  if(label1.eq.' ')label1 = qdt_files(i)(1:l_q-4)
               else
                  label1 = qdt_files(i)(1:l_q-4)
               end if
               call rectangle(fildes,x0,775.,x1,800.)
               call strlength(label1,LEN(label1),l)
               call text2d(fildes,x0+1.,780.,label1(1:l)//char(0),
     +              ANNOTATION_TEXT,0)
            end if
         end do
         call make_picture_current(fildes)

*
*       ***** draw STOP button on bottom of screen if in automatic mode:   *****
*
         if(auto_mode.eq.'y')then !in automatic mode, so draw STOP button
            call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.,0.
     +           ,blank30,0.,0.)
            call draw_button(fildes,5.,90.,205.,165.,auto_stop,60.,120.
     +           ,blank30,0.,0.) 
            label1='  REDRAW'
            call draw_button(fildes,1295.,90.,1495.,165.,label1,1300.
     +           ,120.,blank30,0.,0.) 
c        label1=' +/- SHIFT'
c         call draw_button(fildes,525.,90.,725.,165.,
c     +                          label1,530.,120.,,,)
            if (.not.cth_cch)then
               label1='     +/-'
               call text2d(fildes,526.,145.,label1//char(0),
     +           ANNOTATION_TEXT,0)
               label1=' SINGLE'
               label2='  SHIFT'
               call draw_button(fildes,525.,90.,620.,165.,
     +              label1,526.,120.,label2,526.,95.)
            end if
            if(avg_shift.eq.'y')then !do not offer the option if no average shift to show
               label1='     +/-'
               call text2d(fildes,631.,145.,label1//char(0),
     +              ANNOTATION_TEXT,0)
               if (cth_cch)then
                  label1=' CTRL'
                  label2=''
               else
                  label1='  AVG.'
                  label2=' SHIFT'
               end if
               call draw_button(fildes,630.,90.,725.,165.,
     +              label1,631.,120.,label2,631.,95.)
            end if
            call sc_buttons (fildes)

            if(click.eq.'y')then
               label1='CLICK HERE'
               label2='TO CONTINUE'
               call draw_button(fildes,1295.,5.,1495.,80.,
     +              label1,1300.,50.,label2,1300.,20.)
            end if
*
            call strlength(prim_txt,LEN(prim_txt),l)
            call text2d(fildes,750.,145.,prim_txt(1:l)
     +           //char(0),ANNOTATION_TEXT,0)
            call strlength(prim_txt_stats,LEN(prim_txt_stats),l)
            call text2d(fildes,800.,125.,prim_txt_stats(1:l)
     +           //char(0),ANNOTATION_TEXT,0)
            call strlength(sec_txt,LEN(sec_txt),l)
            call text2d(fildes,750.,105.,sec_txt(1:l)//char(0),
     +           ANNOTATION_TEXT,0)
            call strlength(comm_txt,LEN(comm_txt),l)
            call text2d(fildes,750.,85.,comm_txt(1:l)//char(0),
     +           ANNOTATION_TEXT,0)

*         
*
*       ***** wait here for a mouse click or for time to run out: ***** 
*
 315        call make_picture_current(fildes)
 316        call request_locator(fildes,1,speed,valid,x_dc,y_dc,z_dc)
            if(valid.eq.0)then  !no button has been clicked - timeout invoked
               button_choice = 'auto_view'
               goto 800
            end if
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
*
*       ***** Check to see if the user clicked the mouse and, if so, determine the      *****
*       *****   chosen part of the window.  If the user did not click or if any part of *****
*       *****   the window other than the STOP, PRINT, WRITE or REDRAW buttons was      *****    
*       *****   chosen, proceed with the automatic showing of cell pairs.               *****
*       *****   If the user did indeed click on STOP, clear the bottom of the screen    *****
*       *****   and draw more buttons to proceed with pair analyis.                     *****
*
            screen = 'auto'
            call locate_region(xloc,yloc,regionAUTO,screen)
*
            if(.false.)then     !STOPorCONTINUE_WINDOW
*
*                       !force an appropriate response
            else if(regionAUTO.eq.5)then !choice = STOP --> proceed with analysis
               flats_only = 'n'
               sig_only = 'n'
               goto 399
*
            else if(regionAUTO.eq.20)then !redraw
               button_choice='redraw'
               goto 315
               
            else if(regionAUTO.eq.21.and.qdt_file_version.lt.'8')then !add or remove the SINGLE SHIFT CONTROL PLOT
               show_avg_shift = 'n'
               if(show_controlMAIN.eq.'y')then !currently displaying one of the shifts
                  if(show_single_shift.eq.'y')then !currently displaying single shift
                     show_controlMAIN = 'n'
                     show_single_shift = 'n'
                     show_conf_lim = 'n'
                  else          !currently displaying averaged shift
                     show_single_shift = 'y'
                  end if
               else             !currently not displaying any shift at all
                  show_controlMAIN = 'y'
                  show_single_shift = 'y'
                  show_conf_lim = 'n'
               end if
               button_choice = 'redraw'
               return

            else if(regionAUTO.eq.121)then !add or remove the AVERAGED SHIFT CONTROL PLOT
               if(avg_shift.ne.'y')goto 316 !can't add or remove something that's not available!
               show_single_shift = 'n'
               if(show_controlMAIN.eq.'y')then
                  if(show_avg_shift.eq.'y')then !currently displaying averaged shift
                     show_controlMAIN = 'n'
                     show_avg_shift = 'n'
                     show_conf_lim = 'n'
                  else          !currently displaying single shift
                     show_avg_shift = 'y'
                  end if
               else             !currently not displaying any shift
                  show_controlMAIN = 'y'
                  show_avg_shift = 'y'
                  show_conf_lim = 'n'
               end if
               button_choice = 'redraw'
               return

            else
               goto 316         !inappropriate button click - get another one
               
            end if              ! STOPorCONTINUE_WINDOW
         end if                 !end code for automatic viewing of CCHs
*
*
*
*               *****   The STOP button WAS chosen or not in automatic mode, so...      *****
*
 399     auto_mode = 'n'        !drop out of automatic viewing mode
*
*               *****  Draw analysis buttons at bottom of data screen:  *****
*
         call clear_bottom(fildes) !clear the bottom part of the window 
         label1 = ' '
         label2 = ' '
         call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.,0.,
     +        blank30,0.,0.)
         label1 = 'STOP'
         label2 = 'VIEWING CCHs'
         call draw_button(fildes,5.,90.,205.,165.,
     +        label1,10.,135.,label2,10.,105.)
         if(id.ne.2) then       !id=2 --> look only at unanalyzed cells
            label1='PREV'
            label2='PAIR'
            call draw_button (fildes,110.,5.,205.,80.,
     +           label1,115.,50.,label2,115.,20.)
            label1='FIRST'
            label2='PAIR'
            call draw_button (fildes,5.,5.,100.,80.,
     +           label1,10.,50.,label2,10.,20.)
         endif
         if(mode.ne.'cr')then
            label1='DIRECT'
            label2='ACCESS'
            call draw_button (fildes,275.,5.,450.,80.,
     +           label1,280.,50.,label2,280.,20.)
         end if
         call text2d(fildes,1350.,62.,'N E X T'//char(0),ANNOTATION_TEXT,0)
         call draw_button (fildes,1295.,5.,1375.,60.,'PAIR'//char(0)
     +        ,1300.,20.,blank30,0.,0.)
         call draw_button (fildes,1377.,5.,1435.,60.,'SIG.'//char(0)
     +        ,1380.,20.,blank30,0.,0.)
         call draw_button (fildes,1437.,5.,1495.,60.,'FLAT'//char(0)
     +        ,1440.,20.,blank30,0.,0.)

         call draw_button (fildes,275.,90.,450.,165.,auto_view,280.,120.
     +        ,blank30,0.,0.)
         label1 = '   PRINT'
         call draw_button (fildes,525.,5.,725.,80.,label1,530.,35.
     +        ,blank30,0.,0.)
         call character_height(fildes,.021)
         call character_width(fildes,.007)
         call text2d(fildes,745.,71.,'P'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,745.,60.,'R'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,745.,49.,'E'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,745.,38.,'V'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,748.,27.,'I'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,745.,16.,'E'//char(0),ANNOTATION_TEXT,0)
         call text2d(fildes,740.,5.,'W'//char(0),ANNOTATION_TEXT,0)
         label1 = '   WRITE'
         call draw_button(fildes,775.,5.,975.,80.,label1,780.,35.,
     +        blank30,0.,0.)
         call character_height(fildes,.036)
         call character_width(fildes,.012)
         label1='  REDRAW'
         call draw_button(fildes,1295.,90.,1495.,165.,label1,1300.,105.,
     +        blank30,0.,0.)                                            
         if (.not.cth_cch)then
            label1='     +/-'
            call text2d(fildes,526.,145.,label1//char(0),
     +           ANNOTATION_TEXT,0)
            label1=' SINGLE'
            label2='  SHIFT'
            call draw_button(fildes,525.,90.,620.,165.,
     +           label1,526.,120.,label2,526.,95.)
         end if
         if(avg_shift.eq.'y')then !do not offer the option if no average shift to show
            label1='     +/-'
            call text2d(fildes,631.,145.,label1//char(0),
     +           ANNOTATION_TEXT,0)
            if (cth_cch)then
               label1=' CTRL'
               label2=''
            else
               label1='  AVG.'
               label2=' SHIFT'
            end if
            call draw_button(fildes,630.,90.,725.,165.,
     +           label1,631.,120.,label2,631.,95.)
         end if
         call sc_buttons (fildes)


         if((mode.eq.'cr').or.(mode.eq.'ed'))then
            label1 = 'ENTER'
            label2 = 'ANALYSIS'
            if(((mode.eq.'cr').or.(mode.eq.'ed')).and.
     +           (prim.eq.' '))then
               call text_color(fildes,1.0,0.0,0.0) !red button
               call perimeter_color(fildes,1.0,0.0,0.0)
            end if
            call draw_button (fildes,1025.,5.,1225.,80.,label1,1030.,55.
     +           ,label2,1030.,35.)
            label1 = 'RESULTS'
            call text2d(fildes,1030.,15.,label1
     +           //char(0),ANNOTATION_TEXT,0)
            if(((mode.eq.'cr').or.(mode.eq.'ed')).and.
     +           (prim.eq.' '))then
               call text_color(fildes,0.,0.,0.) !back to black box and text
               call perimeter_color(fildes,0.,0.,0.)
            end if
         end if
*
         call strlength(prim_txt,LEN(prim_txt),l)
         call text2d(fildes,750.,145.,prim_txt(1:l)
     +        //char(0),ANNOTATION_TEXT,0)
         call strlength(prim_txt_stats,LEN(prim_txt_stats),l)
         call text2d(fildes,800.,125.,prim_txt_stats(1:l)
     +        //char(0),ANNOTATION_TEXT,0)
         call strlength(sec_txt,LEN(sec_txt),l)
         call text2d(fildes,750.,105.,sec_txt(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
         call strlength(comm_txt,LEN(comm_txt),l)
         call text2d(fildes,760.,85.,comm_txt(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
*
*               *** wait for user to click on another button or on a histogram ***
*
 499     call make_picture_current(fildes)
         keystroke = 0
         call request_locator(fildes,keystroke,timeout,valid,
     +        x_dc,y_dc,z_dc)
         if(valid.eq.2)then     !keystroke instead of mouse click
            if(keystroke.eq.44)then !',' = '<'
               if(current_qdt.gt.1.and.current_qdt.le.5)then
                  current_qdt = current_qdt-1
               else
                  goto 499
               end if
            elseif(keystroke.eq.46)then !'.' = '>'
               if(current_qdt.ge.1.and.current_qdt.lt.5.and.
     +              current_qdt.lt.total_num_qdts)then
                  current_qdt = current_qdt+1
               else
                  goto 499
               end if
            elseif(keystroke.eq.32.or.keystroke.eq.110)then !spacebar or 'n'
               regionMAIN = 15
            elseif(keystroke.eq.115)then !'s' --> look at next significant pair
               regionMAIN = 151
            elseif(keystroke.eq.102)then !'f' --> look at next flat pair
               regionMAIN = 152
            elseif(keystroke.eq.98)then !'b' --> look at previous pair (backup)
               regionMAIN = 10
            elseif(keystroke.eq.97)then !'a' --> look at first pair
               regionMAIN = 1010
            elseif(keystroke.eq.100)then !'d' --> direct access to a pair
               regionMAIN = 11
            elseif(keystroke.eq.112)then !'p' --> print the main CCH screen
               regionMAIN = 12
            elseif(keystroke.eq.119)then !'w' --> write the main CCH screen
               regionMAIN = 13
            elseif(keystroke.eq.118)then !'v' --> preview main CCH screen using ghostscript
               regionMAIN = 212
            elseif(keystroke.eq.120)then !'x' --> stop viewing CCHs
               regionMAIN = 5
            elseif(keystroke.eq.50)then !'2' --> show +/- 2 stds
               regionMAIN = 18
            elseif(keystroke.eq.51)then !'3' --> show +/- 3 stds
               regionMAIN = 118
            else
               goto 499         !invalid keystroke = ignore it and look for another one
            end if
            if((keystroke.eq.44).or.(keystroke.eq.46))then
               regionMAIN = current_qdt + 100
            end if
         elseif(valid.eq.1)then !mouseclick
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
            screen = 'main'
            call locate_region(xloc,yloc,regionMAIN,screen)
         else
            goto 499
         end if
         screen = 'main'
*
         if(.false.)then        !MAIN_CCH_WINDOW
*                                                       !no valid button chosen - go back and wait for another click
         else if(regionMAIN.eq.1 !a CCH has been chosen
     +           .or.regionMAIN.eq.2
     +           .or.regionMAIN.eq.3
     +           .or.regionMAIN.eq.4)then
*
            call sc_izz (regionMAIN)
            
            do i = 1, 101
               CCH(i)=0
            end do
            rec_num=REC_NUM_X+(regionMAIN-1)
            read(2,c_format,rec=rec_num,err=501)itype,NHW,ICN,CCH !need to re-read to get ICN data
 501        do i = 1, 101
               CONTROL(i)=0
            end do
            rec_num_control_single=rec_num_control_1+(regionMAIN-1) !calc the record #s for each SP CONTROL histogram
            rec_num_control_two=rec_num_control_2+(regionMAIN-1) !calc the location of the SP needed to calc the conf lim
            rec_num_control_average=rec_num_control_avg+(regionMAIN-1) !  SP; for the averaged shift, use the "shifted by 1"

            read(2,c_format,rec=rec_num_control_single,err=1501)
     +           itype,NHW_control,ICN_control,CONTROL_1
 1501       read(2,c_format,rec=rec_num_control_two,err=1502)
     +           itype,NHW_control,ICN_control,CONTROL_2
 1502       read(2,c_format,rec=rec_num_control_average,err=502)
     +           itype,NHW_control,ICN_control,CONTROL_AVG

 502        CCH_original = CCH  !store original data in CCH_original()
            CONTROL_original = CONTROL
            ICN_original = ICN
            title_original = title
*
            call clear_all(fildes)
            call character_height(fildes,.036)
            call character_width(fildes,.012)
 1520       info = 'y'
            if((show_single_shift.eq.'y').or.cth_cch) then
               if(mirror.eq.'n')then
                  do i = 1, 101
                     CONTROL(i) = CONTROL_1(i) !the "shift by 2" SP is used for calculation of confidence
                     CONTROL_cl(i) = CONTROL_2(i) !  limits for the "shift by 1" CONTROL
                  end do
               elseif(mirror.eq.'y')then
                  do i = 1, 101
                     CONTROL(i) = CONTROL_1(102-i)
                     CONTROL_cl(i) = CONTROL_2(102-i)
                  end do
               end if
               rec_num_control = rec_num_control_single
            else if(show_avg_shift.eq.'y')then
               if(mirror.eq.'n')then !in the case of the average shift, use the "shifted by 1" to
                  do i = 1, 101
                     CONTROL(i) = CONTROL_AVG(i) !  compute the confidence limits
                     CONTROL_cl(i) = CONTROL_1(i)
                  end do
               elseif(mirror.eq.'y')then
                  do i = 1, 101
                     CONTROL(i) = CONTROL_AVG(102-i) !  compute the confidence limits
                     CONTROL_cl(i) = CONTROL_1(102-i)
                  end do
               end if
               rec_num_control = rec_num_control_average
            end if

            CONTROL_original = CONTROL

            if(mode.ne.'vt')then !viewing TQDTs --> "view-only" mode
               if(prim.eq.' ')then
                  prim_txt = ' '
                  prim_txt_stats = ' '
               else if(prim.eq.'Flat')then
                  prim_txt = 'FLAT'
                  prim_txt_stats = ' '
               else if(prim.eq.'M P & T')then
                  call strlength(loc1,LEN(loc1),i)
                  prim_txt = 'M P & T '//loc1(1:i)
               else 
                  call strlength(loc1,LEN(loc1),i)
                  call strlength(prim,LEN(prim),j)
                  call remove_all_blanks(cdet,LEN(cdet))
                  call strlength(cdet,LEN(cdet),l)
                  call strlength(statcomm,LEN(statcomm),k)
                  if(cdet.ne.' ')then
c                     prim_txt='Prim: '//loc1(1:i)//' '//
c     +                    prim(1:j)//' (d.i. = '//
c     +                    cdet(1:l)//'; '//
c     +                    statcomm(1:k)//')'//char(0)
                     prim_txt='Prim: '//loc1(1:i)//' '//
     +                    prim(1:j)
                     prim_txt_stats=' (d.i. = '//
     +                    cdet(1:l)//'; '//
     +                    statcomm(1:k)//')'
                  else
                     prim_txt='Prim: '//loc1(1:i)//' '//
     +                    prim(1:j)//' (no statistics calculated)'
                  end if
               end if
            end if
 1521       call clear_all(fildes)
            x0 = 400.
            y0 = 300.
            height = 300.
            width = 525.
            min_bin = CCH(1)
            max_bin = CCH(1)
            do i = 2,101
               min_bin = MIN0(min_bin,CCH(i))   
               max_bin = MAX0(max_bin,CCH(i))
            end do
            if(scaledup.eq.'n')then
               call sc_subtract (0)
               call new_plot(fildes,101,CCH,CONTROL,scaled_hist,
     +              CCHtype,x0,y0,height,width,
     +              info,BDT_FILE,qdt_files(current_qdt),
     +              date,
     +              recording,NHW,rec_num,show_controlMAIN,             
     +              rec_num_control,title,ICN,max_bin,
     +              show_conf_lim,
     +              show_single_shift,show_avg_shift,
     +              CONTROL_cl,show_2_sd,show_3_sd,scaledup,
     +              min_bin,rmagnify,IDs,ITAL)
            elseif(scaledup.eq.'y')then
               subtract_1 = min_bin*(subtract/100.)
               call sc_subtract (subtract_1)
               subtract_cl = subtract_1
               if (cth_cch) subtract_cl = 0
               do i = 1, 101
                  CCH_scaled(i) = CCH(i) - subtract_1
                  CONTROL_scaled(i) = CONTROL(i) - subtract_1
                  CONTROL_cl_scaled(i) = CONTROL_cl(i) - subtract_cl
               end do
               rmagnify = float(max_bin)/(max_bin-subtract_1)
               call new_plot(fildes,101,CCH_scaled,CONTROL_scaled,
     +              scaled_hist,
     +              CCHtype,x0,y0,height,width,
     +              info,BDT_FILE,qdt_files(current_qdt),
     +              date,
     +              recording,NHW,rec_num,show_controlMAIN,             
     +              rec_num_control,title,ICN,max_bin,
     +              show_conf_lim,
     +              show_single_shift,show_avg_shift,
     +              CONTROL_cl_scaled,show_2_sd,show_3_sd,
     +              scaledup,min_bin,rmagnify,IDs,ITAL)
            end if

            call character_height(fildes,.036)
            call character_width(fildes,.012)
            call strlength(prim_txt,LEN(prim_txt),l)
            call strlength(prim_txt_stats,LEN(prim_txt_stats),l1)
            textstring=prim_txt(1:l)//prim_txt_stats(1:l1)
            print '(''prim_txt = '',A)',prim_txt
            print '(''l = '',I5)', l
            print '(''prim_txt_stats = '',A)',prim_txt_stats
            print '(''l1 = '',I5)',l1
            print '(''textstring = '',A)',textstring
            call strlength(textstring,LEN(textstring),l2)
            print '(''l2 = '',I5)', l2
            call text2d(fildes,x0,y0-120.,textstring(1:l2)//char(0),
     +           ANNOTATION_TEXT,0)

c            read (*,'(A)')OK
        
            call character_height(fildes,.030)
            call character_width(fildes,.010)
            if((mode.eq.'vt').or.(mode.eq.'jl'))then
               textstring = '(view-only mode)'
            else
               textstring = ' '
            end if
            call text2d(fildes,5.,y0+height,textstring//char(0),
     +           ANNOTATION_TEXT,0)

            write (REFcodetext1,'(I4)')ICN_original(1)
            write (TARcodetext1,'(I4)')ICN_original(3)
            call text2d(fildes,1000.,700.,REFcodetext1//char(0),
     +           ANNOTATION_TEXT,0)
            call text2d(fildes,1175.,700.,TARcodetext1//char(0),
     +           ANNOTATION_TEXT,0)
            call text2d(fildes,1325.,700.,'ID code'//char(0),
     +           ANNOTATION_TEXT,0)
            if((REFtype.ne.' ').or.(TARtype.ne.' '))then
               call text2d(fildes,1000.,680.,REFtype//char(0),
     +              ANNOTATION_TEXT,0)
               call text2d(fildes,1175.,680.,TARtype//char(0),
     +              ANNOTATION_TEXT,0)
               call text2d(fildes,1325.,680.,'resp. type'//char(0),
     +              ANNOTATION_TEXT,0)
            end if
            call text2d(fildes,1325.,300.,'AP, RL, depth'//char(0),
     +           ANNOTATION_TEXT,0)
            call text2d(fildes,1350.,285.,'(mm)'//char(0),
     +           ANNOTATION_TEXT,0)
            textstring = ' '
            call strlength(REF_coords(1),LEN(REF_coords(1)),l_1)
            call strlength(REF_coords(2),LEN(REF_coords(2)),l_2)
            call strlength(REF_coords(3),LEN(REF_coords(3)),l_3)
            textstring = REF_coords(1)(1:l_1)//', '//
     +           REF_coords(2)(1:l_2)//', '//
     +           REF_coords(3)(1:l_3)
            call strlength(textstring,LEN(textstring),l)
            call text2d(fildes,1000.,300.,textstring(1:l)//char(0),
     +           ANNOTATION_TEXT,0)
            textstring = ' '
            call strlength(TAR_coords(1),LEN(TAR_coords(1)),l_1)
            call strlength(TAR_coords(2),LEN(TAR_coords(2)),l_2)
            call strlength(TAR_coords(3),LEN(TAR_coords(3)),l_3)
            textstring = TAR_coords(1)(1:l_1)//', '//
     +           TAR_coords(2)(1:l_2)//', '//
     +           TAR_coords(3)(1:l_3)
            call strlength(textstring,LEN(textstring),l)
            call text2d(fildes,1175.,300.,textstring(1:l)//char(0),
     +           ANNOTATION_TEXT,0)
            textstring = ' '
            if(show_single_shift.eq.'y')then
               textstring='(single shift-control displayed)'
            else if(show_avg_shift.eq.'y')then
               textstring='(averaged shift-control displayed)'
               if(cth_cch) textstring='(cth control displayed)'
            end if
            call text2d(fildes,x0+220.,y0-70.,textstring//char(0),
     +           ANNOTATION_TEXT,0)
            if(show_conf_lim.eq.'y')then
               if(show_2_sd.eq.'y')then
                  textstring='(confidence limits:  +/- 2 s.d.)'
               else if(show_3_sd.eq.'y')then
                  textstring='(confidence limits:  +/- 3 s.d.)'
               end if
               call text2d(fildes,x0+220.,y0-90.,textstring//char(0),
     +              ANNOTATION_TEXT,0)
            end if
            call sc_label2 (fildes, x0+220, y0-70,show_single_shift
     +           ,show_avg_shift,show_2_sd,show_3_sd)
            call character_height(fildes,.030)
            call character_width(fildes,.010)
            y=670.
            do i = 1, MAX_PERTURB
               if((per_results_REF(i).ne.' ').or.
     +              (per_results_TAR(i).ne.' '))then
                  y = y - 20.
                  call text2d(fildes,1000.,y,per_results_REF(i)//
     +                 char(0),ANNOTATION_TEXT,0)
                  call text2d(fildes,1175.,y,per_results_TAR(i)//
     +                 char(0),ANNOTATION_TEXT,0)
                  call text2d(fildes,1325.,y,per_text_abbrev(i)//
     +                 char(0),ANNOTATION_TEXT,0)
               end if
            end do

            call plot_spikes_per_cycle(sp_per_cycle,num_acc_cycles,
     +           REFcode,TARcode,fildes,IDs,
     +           400.,745.,50.) !x0, y0, height

            info='n'            !default values
            call clear_bottom(fildes)
            call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.
     +           ,0.,blank30,0.,0.)
            call draw_button(fildes,5.,5.,205.,80.,stats,10.,35.
     +           ,blank30,0.,0.)
            call draw_button(fildes,275.,5.,450.,80.,qsum_text,280.,35
     +           .,blank30,0.,0.)                                     
            label1='DIFFERENCE'
            label2='HISTOGRAM'
            call draw_button(fildes,5.,90.,205.,165.,
     +           label1,10.,135.,label2,10.,105.)
            if (.not.cth_cch)then
               label1='     +/-'
               call text2d(fildes,526.,145.,label1//char(0),
     +              ANNOTATION_TEXT,0)
               label1=' SINGLE'
               label2='  SHIFT'
               call draw_button(fildes,525.,90.,620.,165.,
     +              label1,526.,120.,label2,526.,95.)
            end if
            if(avg_shift.eq.'y')then !do not offer the option if no average shift to show
               label1='     +/-'
               call text2d(fildes,631.,145.,label1//char(0),
     +              ANNOTATION_TEXT,0)
               if (cth_cch)then
                  label1=' CTRL'
                  label2=''
               else
                  label1='  AVG.'
                  label2=' SHIFT'
               end if
               call draw_button(fildes,630.,90.,725.,165.,
     +              label1,631.,120.,label2,631.,95.)
            end if
            call sc_buttons (fildes)

            label1='MIRROR IMAGE'
            call draw_button(fildes,525.,5.,725.,80.,label1,530.,35.
     +           ,blank30,0.,0.)                                       
            label1='SCALE UP'
            call strlength(label1,LEN(label1),l)
            call text2d(fildes,810.,60.,label1(1:l)//char(0),
     +           ANNOTATION_TEXT,0)
            label1='75%'
            call draw_button(fildes,775.,5.,873.,58.,label1,795.,25.
     +           ,blank30,0.,0.)                                       
            label1='custom'
            call draw_button(fildes,877.,5.,975.,58.,label1,887.,25.
     +           ,blank30,0.,0.)                                       

            if((show_avg_shift.eq.'y').or.(show_single_shift.eq.'y'))
     +           then
               label1=' +/- 2'
               label2='std dev'
               call draw_button_color(fildes,775.,90.,870.,165.,
     +              label1,785.,135.,label2,785.,105.,'u')
               
               label1=' +/- 3'
               label2='std dev'
               call draw_button_color(fildes,880.,90.,975.,165.,
     +              label1,885.,135.,label2,885.,105.,'m')
            end if

            label1='PRINT'
            call draw_button(fildes,1025.,90.,1120.,165.,label1,1040.,  
     +           120.,blank30,0.,0.)
            label1='WRITE'
            label2='to disk'
            call draw_button(fildes,1130.,90.,1225.,165.,
     +           label1,1135.,135.,label2,1135.,105.)
            label1='ORIGINAL'
            label2='PLOT'
            call draw_button(fildes,275.,90.,450.,165.,label1,280.,135.,
     +           label2,280.,105.)
            label1='CHANGE TITLE'
            call draw_button(fildes,1025.,5.,1225.,80.,label1,1030.,35.,
     +           blank30,0.,0.)
            call draw_button(fildes,1295.,5.,1495.,80.,return_box,1300.,
     +           35.,blank30,0.,0.)
            label1='  REDRAW'
            call draw_button(fildes,1295.,90.,1495.,165.,label1,1300.,  
     +           120.,blank30,0.,0.)
*
*               *** wait for user to select a control button: ***
*
 525        call make_picture_current(fildes)
            call request_locator(fildes,1,timeout,valid,x_dc,y_dc,z_dc)
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
            screen = 'single'
            call locate_region(xloc,yloc,regionSINGLE,screen)
*
 5526       if(.false.)then     !ENLARGED_CCH_WINDOW
*
            else if(regionSINGLE.eq.20)then !REDRAW
               goto 1520
               
            else if(regionSINGLE.eq.5)then !plot the DIFFERENCE HISTOGRAM
               do i = 1, 101
                  DIFF(i) = CCH(i) - CONTROL(i)
               end do
               info='y'
               show_control_DIFF='y'
               call clear_all(fildes)
               x0 = 200.
               y0 = 300.
               if(scaledup.eq.'y')then !plot CCH with shift-control
                  call new_plot(fildes,101,CCH_scaled,CONTROL_scaled
     +                 ,scaled_hist,CCHtype,x0,y0,300.,525.,info
     +                 ,BDT_FILE,qdt_files(current_qdt),date,recording
     +                 ,NHW,rec_num,show_control_DIFF,rec_num_control
     +                 ,title,ICN,max_bin,'','','',[integer::],'','',
     +                 scaledup,min_bin,rmagnify,IDs,ITAL) 
               else if(scaledup.eq.'n')then !plot CCH with shift-control
                  call new_plot(fildes,101,CCH,CONTROL,scaled_hist
     +                 ,CCHtype,x0,y0,300.,525.,info,BDT_FILE
     +                 ,qdt_files(current_qdt),date,recording,NHW
     +                 ,rec_num,show_control_DIFF,rec_num_control,title
     +                 ,ICN,max_bin,'','','',[integer::],'','',scaledup,
     +                 min_bin,rmagnify,IDs,ITAL)               
               end if
               if(show_single_shift.eq.'y')then
                  textstring='(single shift-control displayed)'
               else if(show_avg_shift.eq.'y')then
                  textstring='(averaged shift-control displayed)'
                  if(cth_cch) textstring='(cth control displayed)'
               else
                  textstring=' '
               end if
               call text2d(fildes,x0+220.,y0-70.,textstring//char(0),
     +              ANNOTATION_TEXT,0)
               info = 'n'
               call new_plot(fildes,101,DIFF,DIFF,scaled_hist,DIF,900., !draw the difference histogram
     +              300.,300.,525.,info,'','','','',0,0,'',0,'',
     +              [integer::],max_bin,'','','',[integer::],'','',
     +              scaledup,min_bin,rmagnify,IDs,ITAL)
 5251          call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.,0.
     +              ,blank30,0.,0.)
               label1='   PRINT'
               call draw_button(fildes,775.,5.,975.,80.,label1,780.,35.
     +              ,blank30,0.,0.)                                       
               label1='WRITE to disk'
               label2=' in ps format'
               call draw_button(fildes,1025.,5.,1225.,80.,label1,
     +              1030.,50.,label2,1030.,20.)
               label1='  REDRAW'
               call draw_button(fildes,1295.,90.,1495.,165.,label1,1300.
     +              ,120.,blank30,0.,0.)
               call draw_button(fildes,1295.,5.,1495.,80.,return_box
     +              ,1300.,35.,blank30,0.,0.)
*
*               *** wait for user to select the return button: ***
*
 526           call make_picture_current(fildes)
               call request_locator(fildes,1,timeout,valid,x_dc,y_dc,
     +              z_dc)
               call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
               screen = 'dif'
               call locate_region(xloc,yloc,regionDIF,screen)
*
               if(.false.)then  !DIFF_HISTOGRAM_WINDOW
*
               else if(regionDIF.eq.13.or.regionDIF.eq.14.or.
     +                 regionDIF.eq.213)then !PRINT or WRITE or VIEW the difference histogram
                  if(regionDIF.eq.13)task = 'p'
                  if(regionDIF.eq.14)task = 'w'
                  if(regionDIF.eq.213)task = 'v'
                  if(scaledup.eq.'y')then
                     call print_or_write_DIFF_CCH(task,CCH_scaled,
     +                    CONTROL_scaled,BDT_FILE, !print the DIFFERENCE HISTOGRAM
     +                    qdt_files(current_qdt),REFcode,TARcode,
     +                    date,recording,
     +                    NHW,rec_num,rec_num_control,title,ICN,
     +                    min_bin,max_bin,scaledup,rmagnify,
     +                    IDs,ITAL)
                  else if(scaledup.eq.'n')then
                     call print_or_write_DIFF_CCH(task,CCH,CONTROL,
     +                    BDT_FILE, !print the DIFFERENCE HISTOGRAM
     +                    qdt_files(current_qdt),REFcode,TARcode,
     +                    date,recording,
     +                    NHW,rec_num,rec_num_control,title,ICN,
     +                    min_bin,max_bin,scaledup,rmagnify,
     +                    IDs,ITAL)
                  end if
                  goto 5251
*
               else if(regionDIF.eq.20)then !REDRAW
                  goto 5526

               else if(regionDIF.eq.15)then !RETURN to previous window
                  goto 1520
               else
                  goto 526
*
               end if           !DIFF_HISTOGRAM_WINDOW
*
            else if(regionSINGLE.eq.6)then !ORIGINAL PLOT
               CCH = CCH_original !reload CCH() with original data
               CONTROL = CONTROL_original !reload CONTROL() with original data
               ICN = ICN_original
               write (REFcodetext,'(I4)')ICN(1)
               write (TARcodetext,'(I4)')ICN(3)
               title = title_original !restore the default title
               mirror = 'n'
               call sc_mirror (0)
               scaledup = 'n'
               goto 1520
*
            else if(regionSINGLE.eq.21.and.qdt_file_version.lt.'8')then !add or remove the SINGLE SHIFT CONTROL PLOT
               show_avg_shift = 'n'
               if(show_controlMAIN.eq.'y')then !currently displaying one of the shifts
                  if(show_single_shift.eq.'y')then !currently displaying single shift
                     show_controlMAIN = 'n'
                     show_single_shift = 'n'
                     show_conf_lim = 'n'
                  else          !currently displaying averaged shift
                     show_single_shift = 'y'
                  end if
               else             !currently not displaying any shift at all
                  show_controlMAIN = 'y'
                  show_single_shift = 'y'
                  show_conf_lim = 'n'
               end if
               goto 1520

            else if((regionSINGLE.eq.121).and.(avg_shift.eq.'y'))then !add or remove the AVERAGED SHIFT CONTROL PLOT
               show_single_shift = 'n'
               if(show_controlMAIN.eq.'y')then
                  if(show_avg_shift.eq.'y')then !currently displaying averaged shift
                     show_controlMAIN = 'n'
                     show_avg_shift = 'n'
                     show_conf_lim = 'n'
                  else          !currently displaying single shift
                     show_avg_shift = 'y'
                  end if
               else             !currently not displaying any shift
                  show_controlMAIN = 'y'
                  show_avg_shift = 'y'
                  show_conf_lim = 'n'
               end if
               goto 1520
            else if(regionSINGLE.eq.206)then !add or remove the SURROGATE CONTROL PLOT
               call sc_control (-1)
               goto 1520
            else if(regionSINGLE.eq.211)then !add or remove the SURROGATE CONTROL CL PLOT
               call sc_cl (-1)
               goto 1520

            else if(regionSINGLE.eq.18)then !redraw with CONFIDENCE LIMITS (mean +/- 2 std devs)
               show_3_sd = 'n'
               if(show_conf_lim.eq.'y')then
                  if(show_2_sd.eq.'y')then !currently displaying 2 s.d.
                     show_conf_lim = 'n'
                     show_2_sd = 'n'
                  else          !currently displaying 3 s.d.
                     show_2_sd = 'y'
                  end if
               else             !currently not displaying any confidence limits
                  show_conf_lim = 'y'
                  show_2_sd = 'y'
               end if
               goto 1520


            else if(regionSINGLE.eq.118)then !redraw with CONFIDENCE LIMITS (mean +/- 3 std devs)
               show_2_sd = 'n'
               if(show_conf_lim.eq.'y')then
                  if(show_3_sd.eq.'y')then !currently displaying 3 s.d.
                     show_conf_lim = 'n'
                     show_3_sd = 'n'
                  else          !currently displaying 2 s.d.
                     show_3_sd = 'y'
                  end if
               else             !currently not displaying any confidence limits
                  show_conf_lim = 'y'
                  show_3_sd = 'y'
               end if
               goto 1520
               
            else if(regionSINGLE.eq.10)then !STATISTICAL ANALYSIS
 527           call clear_bottom(fildes)
               call clear_all(fildes)
               call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.
     +              ,0.,blank30,0.,0.)
               if(show_single_shift.eq.'n')show_avg_shift = 'y'
               if(show_2_sd.eq.'n')show_3_sd = 'y'
               show_conf_lim = 'y'
               binwidth=float(NHW)/50.
               info='y'
               show_control_STATS='y'
               if((show_single_shift.eq.'y').or.cth_cch)then
                  CONTROL = CONTROL_1
                  CONTROL_cl = CONTROL_2
                  rec_num_control = rec_num_control_single
               end if
               if(show_avg_shift.eq.'y'.and..not.cth_cch)then
                  CONTROL = CONTROL_AVG
                  CONTROL_cl = CONTROL_1
                  rec_num_control = rec_num_control_average
               end if
               if(mirror.eq.'y')then
                  do i = 1,101
                     temp1(i) = CONTROL(102-i)
                     temp2(i) = CONTROL_cl(102-i)
                  end do
                  do i = 1,101
                     CONTROL(i) = temp1(i)
                     CONTROL_cl(i) = temp2(i)
                  end do
               end if
               x0 = 500.
               y0 = 300.
               if(scaledup.eq.'y')then
                  call new_plot(fildes,101,CCH_scaled,CONTROL_scaled, !plot the shift-control
     +                 scaled_hist,CCHtype,
     +                 500.,300.,300.,525.,info,BDT_FILE,
     +                 qdt_files(current_qdt),date,recording,NHW,
     +                 rec_num,show_control_STATS,
     +                 rec_num_control,title,ICN,max_bin,
     +                 show_conf_lim,
     +                 show_single_shift,show_avg_shift,
     +                 CONTROL_cl_scaled,show_2_sd,show_3_sd,
     +                 scaledup,min_bin,rmagnify,IDs,ITAL)
               elseif(scaledup.eq.'n')then
                  call new_plot(fildes,101,CCH,CONTROL,scaled_hist, !plot the shift-control
     +                 CCHtype,
     +                 500.,300.,300.,525.,info,BDT_FILE,
     +                 qdt_files(current_qdt),date,recording,NHW,
     +                 rec_num,show_control_STATS,
     +                 rec_num_control,title,ICN,max_bin,
     +                 show_conf_lim,
     +                 show_single_shift,show_avg_shift,
     +                 CONTROL_cl,show_2_sd,show_3_sd,scaledup,
     +                 min_bin,rmagnify,IDs,ITAL)
               end if
               if(show_single_shift.eq.'y')then
                  textstring='(single shift-control displayed)'
               else if(show_avg_shift.eq.'y')then
                  textstring='(averaged shift-control displayed)'
                  if(cth_cch) textstring='(cth control displayed)'
               else
                  textstring = ' '
               end if
               call text2d(fildes,x0+220.,y0-70.,textstring//char(0),
     +              ANNOTATION_TEXT,0)
               if(show_conf_lim.eq.'y')then
                  if(show_2_sd.eq.'y')then
                     textstring='(confidence limits:  +/- 2 s.d.)'
                  else if(show_3_sd.eq.'y')then
                     textstring='(confidence limits:  +/- 3 s.d.)'
                  end if
                  call text2d(fildes,x0+220.,y0-90.,textstring//
     +                 char(0),ANNOTATION_TEXT,0)
               end if
               info = 'n'

               call plot_spikes_per_cycle(sp_per_cycle,num_acc_cycles,
     +              REFcode,TARcode,fildes,IDs,
     +              400.,745.,50.) !x0, y0, height)
               call character_height(fildes,.036)
               call character_width(fildes,.012)
               if (.not.cth_cch)then
                  label1='   USE'
                  call text2d(fildes,526.,145.,label1//char(0),
     +                 ANNOTATION_TEXT,0)
                  label1=' SINGLE'
                  label2='  SHIFT'
                  call draw_button(fildes,525.,90.,620.,165.,
     +                 label1,526.,120.,label2,526.,95.)
               end if
               if((avg_shift.eq.'y').and.(.not.cth_cch))then !do not offer the option if no average shift to show
                  label1='  USE'
                  call text2d(fildes,631.,145.,label1//char(0),
     +                 ANNOTATION_TEXT,0)
                  label1='  AVG.'
                  label2=' SHIFT'
                  call draw_button(fildes,630.,90.,725.,165.,
     +                 label1,631.,120.,label2,631.,95.)
               end if

               call character_height(fildes,0.036)
               call character_width(fildes,.012)
               label1='CALCULATE STATISTICS USING:'
               call text2d(fildes,20.,50.,label1
     +              //char(0),ANNOTATION_TEXT,0)
               if (cth_cch)then
                  label1='CONTROL'
               else
                  label1='SHIFT-CONTROL'
               end if
               label2='   CCH'
               call draw_button (fildes,525.,5.,725.,80.,
     +              label1,530.,50.,label2,530.,20.)
               label1='  RANGE OF'
               label2='    BINS'
               call draw_button (fildes,775.,5.,975.,80.,
     +              label1,780.,50.,label2,780.,20.)
               label1='  CANCEL'
               call draw_button (fildes,1025.,5.,1225.,80.,label1,1030.,   
     +              35.,blank30,0.,0.)
 528           call make_picture_current(fildes)
               call request_locator(fildes,1,timeout,valid,x_dc,y_dc,
     +              z_dc)
               call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
               screen = 'stats'
               call locate_region(xloc,yloc,regionSTATS,screen)
*
*
               STATISTICS_WINDOW: if(.false.)then
*                       !force a valid response
               else if(regionSTATS.eq.14)then !CANCEL
                  goto 1520
               else if(regionSTATS.eq.21.and.qdt_file_version.lt.'8')
     +                 then     !use single shift control 
                  show_single_shift = 'y'
                  show_avg_shift = 'n'
                  goto 527
               else if(regionSTATS.eq.121)then !use averaged shift control
                  show_avg_shift = 'y'
                  show_single_shift = 'n'
                  goto 527
               else if(regionSTATS.eq.12.or.regionSTATS.eq.13)then !proceed with statistical analyses
                  if(regionSTATS.eq.12)then
                     if((show_single_shift.eq.'n').and.
     +                    (show_avg_shift.eq.'n'))goto 527 !must have picked a shift control for this calculation
                  end if
c                if(show_single_shift.eq.'y')then
c                   CONTROL = CONTROL_1  
c                   CONTROL_cl = CONTROL_2
c                else if(show_avg_shift.eq.'y')then      
c                   CONTROL = CONTROL_AVG                
c                   CONTROL_cl = CONTROL_1
c                end if
                  info='y'
                  call clear_bottom(fildes)
                  call clear_all(fildes)
                  call draw_button(fildes,1.,1.,1500.,170.,perimeter
     +                 ,0.,0.,blank30,0.,0.)
                  if(scaledup.eq.'y')then
                     call new_plot(fildes,101,CCH_scaled,CONTROL_scaled,
     +                    scaled_hist,CCHtype,
     +                    500.,300.,300.,525.,info,BDT_FILE,
     +                    qdt_files(current_qdt),date,recording,NHW,
     +                    rec_num,show_control_STATS,
     +                    rec_num_control,title,ICN,max_bin,
     +                    show_conf_lim,
     +                    show_single_shift,show_avg_shift,
     +                    CONTROL_cl_scaled,show_2_sd,show_3_sd,
     +                    scaledup,min_bin,rmagnify,IDs,ITAL)
                  else if(scaledup.eq.'n')then
                     call new_plot(fildes,101,CCH,CONTROL,scaled_hist,
     +                    CCHtype,
     +                    500.,300.,300.,525.,info,BDT_FILE,
     +                    qdt_files(current_qdt),date,recording,NHW,
     +                    rec_num,show_control_STATS,
     +                    rec_num_control,title,ICN,max_bin,
     +                    show_conf_lim,
     +                    show_single_shift,show_avg_shift,
     +                    CONTROL_cl,show_2_sd,show_3_sd,scaledup,
     +                    min_bin,rmagnify,IDs,ITAL)
                  end if
                  pixels_per_bin=525./101. !=width / # bins in CCH
                  call statistics(fildes,mouse,CCH,CONTROL,scaled_hist,
     +                 binwidth,ZK,
     +                 Probk,det,vis,ZLAT,HALFWD,info,
     +                 pixels_per_bin,regionSTATS,NHW,REFcode,
     +                 TARcode,REFtype,TARtype,
     +                 qdt_files(current_qdt),
     +                 BDT_FILE,WINDOW1,
     +                 mode,Q_pos,statcomm,
     +                 czk,cprobk,cdet,cvis,chalfwd,czlat,
     +                 rec_num,rec_num_control,
     +                 show_control,show_conf_lim,show_single_shift, !STATS???
     +                 show_avg_shift,show_2_sd,show_3_sd,
     +                 per_text_abbrev,
     +                 per_results_REF,per_results_TAR,date,recording,
     +                 ICN,title,CONTROL_cl,
     +                 min_bin,max_bin,scaledup,rmagnify,IDs,ITAL,
     +                 CCH_scaled,CONTROL_scaled,CONTROL_cl_scaled,
     +                 stats_bw)

                  goto 1520
               else
                  goto 528
               end if STATISTICS_WINDOW
            else if(regionSINGLE.eq.11)then !Q_SUM
               call clear_all(fildes)
               call new_plot(fildes,101,CCH,CONTROL,scaled_hist,CCHtype,
     +              100.,300.,300.,525.,info,BDT_FILE,
     +              qdt_files(current_qdt),date,recording,NHW,
     +              rec_num,show_control,rec_num_control,title,
     +              ICN,max_bin,
     +              show_conf_lim,
     +              show_single_shift,show_avg_shift,
     +              CONTROL_cl,show_2_sd,show_3_sd,scaledup,
     +              min_bin,rmagnify,IDs,ITAL)
               binwidth=float(NHW)/50.
*
               call qsum(fildes,mouse,CCH,
     +              binwidth,ICN,coefval,
     +              BDT_FILE,qdt_files(current_qdt),date,
     +              recording,NHW,rec_num,title,
     +              REFcode,
     +              TARcode,WINDOW1,WINDOW2,mode,IDs,ITAL)
               goto 1520
*
            else if(regionSINGLE.eq.12)then !MIRROR IMAGE
               if(mirror.eq.'y')then
                  print '(/,''If I mirror a plot that is already ''
     +''mirrored, there''''s no telling what ''
     +''might happen!'')'
                  goto 1521
               end if
               mirror = 'y'
               call sc_mirror (1)
               call clear_top(fildes)
               do i = 1,101     !load CCH into CCH_mirror backwards
                  CCH_mirror(i)=CCH(102-i)
                  CONTROL_mirror(i)=CONTROL(102-i)
                  CONTROL_cl_mirror(i)=CONTROL_cl(102-i)
               end do
               CCH = CCH_mirror
               CONTROL = CONTROL_mirror
               CONTROL_cl = CONTROL_cl_mirror
               ICN(1)=ICN_original(3)
               ICN(2)=ICN_original(4)
               ICN(3)=ICN_original(1)
               ICN(4)=ICN_original(2)
               write (REFcodetext,'(I4)')ICN(1)
               write (TARcodetext,'(I4)')ICN(3)
               title = 'MIRROR IMAGE: '//REFcodetext//' > '//TARcodetext
               info='y'
               goto 1521
*
            else if(regionSINGLE.eq.131)then !auto SCALE UP - subtract 75% of minimum bin
               subtract = 75
               scaledup = 'y'
               goto 1520

            else if(regionSINGLE.eq.13)then !SCALE UP
 540           print '(//,T5,''ENTER percentage of minimum bin ''
     +''to be subtracted from all bins  >> '',$)'
               read (*,fmt='(I3)',err=540)subtract 
               if(subtract.ne.0)scaledup = 'y'
               goto 1520
*
            else if(regionSINGLE.eq.14)then !CHANGE TITLE
               call clear_bottom(fildes)
 541           print '(//,T5,''ENTER new title (30 char. max.) ..or.. ''
     +''x to restore default title''
     +''  >> '',$)'
               title = ' '
               read (*,fmt='(A)',err=541)title
               if(title.eq.'x')title=title_original
               call clear_all(fildes)
               info='y'
               goto 1520
*
            else if(regionSINGLE.eq.19 !PRINT or WRITE the window
     +              .or.regionSINGLE.eq.119.or.regionSINGLE.eq.219)then
               if(regionSINGLE.eq.19)task='p'
               if(regionSINGLE.eq.119)task='w'
               if(regionSINGLE.eq.219)task='v'
               if(scaledup.eq.'n')then
                  call print_or_write_ENLARGED_CCH(task,CCH,CONTROL,
     +                 BDT_FILE,qdt_files(current_qdt),
     +                 REFcode,TARcode,date,recording,
     +                 NHW,rec_num,rec_num_control,title,
     +                 ICN,show_control,
     +                 show_conf_lim,REFtype,TARtype,per_results_REF,
     +                 per_results_TAR,per_text_abbrev,
     +                 show_single_shift,show_avg_shift,CONTROL_cl,
     +                 show_2_sd,show_3_sd,prim_txt,min_bin,max_bin,
     +                 scaledup,rmagnify,IDs,ITAL,
     +                 REF_coords,TAR_coords,prim_txt_stats)
               else if(scaledup.eq.'y')then
                  call print_or_write_ENLARGED_CCH(task,CCH_scaled,
     +                 CONTROL_scaled,
     +                 BDT_FILE,qdt_files(current_qdt),
     +                 REFcode,TARcode,date,recording,
     +                 NHW,rec_num,rec_num_control,title,
     +                 ICN,show_control,
     +                 show_conf_lim,REFtype,TARtype,per_results_REF,
     +                 per_results_TAR,per_text_abbrev,
     +                 show_single_shift,show_avg_shift,
     +                 CONTROL_cl_scaled,
     +                 show_2_sd,show_3_sd,prim_txt,min_bin,max_bin,
     +                 scaledup,rmagnify,IDs,ITAL,
     +                 REF_coords,TAR_coords,prim_txt_stats)
               end if
               goto 1520
*
            else if(regionSINGLE.eq.15)then !RETURN 
               button_choice='leave_stats'
               flag = 1         !set a flag to tell program that the user has
! chosen to investigate an individual CCH more
! closely, perhaps for statistical analysis of a
! primary feature --> display ENTER ANALYSIS RESULTS
! button in red to prompt user to enter results
! if the CCH is not flat
               return           !go back to analyze_data and redisplay the histograms
! for this pair - be careful NOT TO RESUME AUTO MODE
            else
               goto 525
*
            end if              !ENLARGED_CCH_WINDOW
*
*
*
         elseif(regionMAIN.eq.106)then !toggle to scale up the CCHs in the MAIN screen by 75% of min bin
            if(scaledMAIN.eq.'n')then
               scaledMAIN = 'y'
            elseif(scaledMAIN.eq.'y')then
               scaledMAIN = 'n'
            end if
            button_choice = 'redraw'
            return
         elseif(regionMAIN.eq.107)then !toggle to turn on and off the I and E pulse lines in the CTHs
            if(show_pulses.eq.'n')then
               show_pulses = 'y'
            elseif(show_pulses.eq.'y')then
               show_pulses = 'n'
            end if
            button_choice = 'redraw'
            return
         elseif (regionMAIN.ge.41.and.regionMAIN.le.48)then !show enlarged ACH
            if(regionMAIN.eq.41)rec_num_ach=(included(REFcode)*4)+
     +           (3*total_num_cells)
            if(regionMAIN.eq.42)rec_num_ach=(included(TARcode)*4)+
     +           (3*total_num_cells)
            if(regionMAIN.eq.43)rec_num_ach=(included(REFcode)*4)+
     +           (3*total_num_cells)+1
            if(regionMAIN.eq.44)rec_num_ach=(included(TARcode)*4)+
     +           (3*total_num_cells)+1
            if(regionMAIN.eq.45)rec_num_ach=(included(REFcode)*4)+
     +           (3*total_num_cells)+2
            if(regionMAIN.eq.46)rec_num_ach=(included(TARcode)*4)+
     +           (3*total_num_cells)+2
            if(regionMAIN.eq.47)rec_num_ach=(included(REFcode)*4)+
     +           (3*total_num_cells)+3
            if(regionMAIN.eq.48)rec_num_ach=(included(TARcode)*4)+
     +           (3*total_num_cells)+3
            do i = 1, NUM_BINS
               big_ACH(i)=0
            end do
            read(2,c_format,REC=(rec_num_ach),err=403) 
     +           itype,NHW,ICN,big_ACH
 403        call clear(fildes)
            x0 = 400.
            y0 = 300.
            height = 300.
            width = 525.
            min_bin = big_ACH(1)
            max_bin = big_ACH(1)
            do i = 2,100
               min_bin = MIN0(min_bin,big_ACH(i))   
               max_bin = MAX0(max_bin,big_ACH(i))
            end do
            if(MOD(regionMAIN,2).eq.1)then
               title = 'ACH for '//REFname//' (IDcode = '//
     +              REFcodetext//')'
            else
               title = 'ACH for '//TARname//' (IDcode = '//
     +              TARcodetext//')'
            end if
            call new_plot(fildes,100,big_ACH,big_ACH,scaled_hist,
     +           ACH,x0,y0,height,width,
     +           'y',BDT_FILE,qdt_files(current_qdt),
     +           date,
     +           recording,NHW,rec_num_ACH,'n',              
     +           rec_num_control,title,ICN,max_bin,
     +           'n','n','n',CONTROL_cl,'n','n','n',
     +           min_bin,rmagnify,IDs,ITAL)
            label1='    RETURN'
            call draw_button(fildes,5.,90.,205.,165.,label1,10.,120.
     +           ,blank30,0.,0.) 
            call make_picture_current(fildes)
 404        call request_locator(fildes,1,timeout,valid,x_dc,y_dc,
     +           z_dc)
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
            screen = 'main'
            call locate_region(xloc,yloc,regionMAIN,screen)
            if(regionMAIN.eq.5)then
            else
               goto 404
            end if
            button_choice='redraw'
            return


         else if(regionMAIN.eq.5)then !RETURN TO MAIN MENU button - don't forget to set the flag
            button_choice = 'main_menu'
            goto 800            !goto C-P-M menu
*
         else if(regionMAIN.eq.10)then ! PREVIOUS PAIR button
            button_choice='prev_pair' !return to analyze_data (calling code) with a flag set to
! tell the pgm to read in the immediately preceding pair
            sig_only = 'n'
            flats_only = 'n'
            goto 800
         else if(regionMAIN.eq.1010)then ! FIRST PAIR button
            button_choice='first_pair' !return to analyze_data (calling code) with a flag set to
! tell the pgm to read in the first pair
            sig_only = 'n'
            flats_only = 'n'
            goto 800
*
         else if(regionMAIN.eq.11)then !DIRECT ACCESS button - ask user for REF and TAR (console),
! calculate the location of their CCHs and display them
! after clearing the entire window
            if(mode.eq.'cr')goto 499 !cannot directly access a pair if in initial data entry mode
            button_choice='dir_acc'
            goto 800
*
         else if(regionMAIN.eq.12 !PRINT or WRITE the graphics window
     +           .or.regionMAIN.eq.13.or.regionMAIN.eq.212)then
            if(regionMAIN.eq.12)task='p'
            if(regionMAIN.eq.13)task='w'
            if(regionMAIN.eq.212)task='v'
            if((show_single_shift.eq.'y').or.cth_cch)then
               do i = 1, NUM_BINS
                  do j = 1, 4
                     CONTROL_CCHs(j,i) = CONTROLs_1(j,i)
                     CONTROLs_cl(j,i) = CONTROLs_2(j,i)
                  end do
               end do
            else if(show_avg_shift.eq.'y')then
               do i = 1, NUM_BINS
                  do j = 1, 4
                     CONTROL_CCHs(j,i) = CONTROLs_AVG(j,i)
                     CONTROLs_cl(j,i) = CONTROLs_1(j,i)
                  end do
               end do
            end if
            call print_or_write_ALL_CCHs(task,date,recording, BDT_FILE
     +           ,qdt_files(current_qdt), REFcode,TARcode,REFname
     +           ,TARname,REFtype,TARtype, CCHs,CONTROL_CCHs,REF_ACHs
     +           ,TAR_ACHs,BINWtexts, CTH_REF,CTH_TAR,phrenic_hist,
     +           ETA2_REF,ETA2_TAR,zmodsig_REF, zmodsig_TAR,zmodsig2_REF
     +           ,zmodsig2_TAR,cardiac_REF, cardiac_TAR,cardiac_hist
     +           ,card_type_REF,card_type_TAR, prim,sec,loc1,loc2,rescom
     +           ,statcomm, cdet,show_controlMAIN,c_format,included
     +           ,total_num_cells, DELTA2_REF,DELTA2_TAR,tedfactor_REF
     +           ,tedfactor_TAR, show_single_shift,show_avg_shift
     +           ,show_conf_lim, show_2_sd,show_3_sd,CONTROLs_cl,min_bin
     +           ,max_bin,ign1,rmagnify,IDs,ITAL,cardiac_pls,
     +           REF_coords,TAR_coords,scaledMAIN, show_pulses,mean_E
     +           ,NORM_OFFSET,bw_n)

            goto 399
*
         else if(regionMAIN.eq.15)then !NEXT PAIR button
            button_choice='next_pair'
            sig_only = 'n'
            flats_only = 'n'
            goto 800
         else if(regionMAIN.eq.151)then !NEXT SIG FEATURE ONLY  button
            button_choice='next_pair'
            sig_only = 'y'
            flats_only = 'n'
            goto 800
         else if(regionMAIN.eq.152)then !NEXT FLAT ONLY  button
            button_choice='next_pair'
            sig_only = 'n'
            flats_only = 'y'
            goto 800
*
         else if(regionMAIN.eq.101)then !switch to QDT file #1
            if(qdt_files(1).ne.' ')then
               if(c_format.eq.'(109I6)')then !version 2 or lower
                  OPEN (2,FILE=qdt_files(1),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=654)
               else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
                  OPEN (2,FILE=qdt_files(1),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=1082)
               else
               end if
               QDT_FILENAME = qdt_files(1)
               current_qdt = 1
               call sc_qdt (current_qdt)
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.102)then !switch to QDT file #2
            if(qdt_files(2).ne.' ')then
               if(c_format.eq.'(109I6)')then !version 2 or lower
                  OPEN (2,FILE=qdt_files(2),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=654)
               else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
                  print *,qdt_files(2)
                  OPEN (2,FILE=qdt_files(2),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=1082)
               else
               end if
               QDT_FILENAME = qdt_files(2)
               current_qdt = 2
               call sc_qdt (current_qdt)
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.103)then !switch to QDT file #3
            if(qdt_files(3).ne.' ')then
               if(c_format.eq.'(109I6)')then !version 2 or lower
                  OPEN (2,FILE=qdt_files(3),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=654)
               else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
                  OPEN (2,FILE=qdt_files(3),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=1082)
               else
               end if
               QDT_FILENAME = qdt_files(3)
               current_qdt = 3
               call sc_qdt (current_qdt)
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.104)then !switch to QDT file #4
            if(qdt_files(4).ne.' ')then
               if(c_format.eq.'(109I6)')then !version 2 or lower
                  OPEN (2,FILE=qdt_files(4),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=654)
               else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
                  OPEN (2,FILE=qdt_files(4),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=1082)
               else
               end if
               QDT_FILENAME = qdt_files(4)
               current_qdt = 4
               call sc_qdt (current_qdt)
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.105)then !switch to QDT file #5
            if(qdt_files(5).ne.' ')then
               if(c_format.eq.'(109I6)')then !version 2 or lower
                  OPEN (2,FILE=qdt_files(5),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=654)
               else if(c_format.eq.'(2I6,107I10)')then !version 3 or higher
                  OPEN (2,FILE=qdt_files(5),FORM='FORMATTED',
     +                 STATUS='OLD',ACCESS='DIRECT',RECL=1082)
               else
               end if
               QDT_FILENAME = qdt_files(5)
               current_qdt = 5
               call sc_qdt (current_qdt)
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.6)then !AUTOVIEW button - proceed with slide show presentation
            sig_only = 'y'      !default is to look at pairs w/ significant features in CCH
            flats_only = 'n'
 550        call clear_bottom(fildes) !clear the bottom of the window
            if((mode.eq.'jl').or.(mode.eq.'ed'))then
               textstring = 'VIEW:'
               call text2d(fildes,100.,60.,textstring
     +              //char(0),ANNOTATION_TEXT,0)
               call draw_button (fildes,200.,90.,250.,140.,view_all,185. !view all CCHs
     +              ,70.,blank30,0.,0.)
               view_1 = 'significant'
               view_2 = ' features'
               view_3 = '   only'
               call draw_button(fildes,325.,90.,375.,140., !default - view significant CCHs
     +              view_1,295.,70.,
     +              view_2,295.,52.)
               call strlength(view_3,LEN(view_3),l)
               call text2d(fildes,295.,34.,view_3(1:l)//char(0),
     +              ANNOTATION_TEXT,0)
               view_1 = 'FLATs'
               view_2 = ' only'
               call draw_button(fildes,450.,90.,500.,140.,
     +              view_1,450.,70.,
     +              view_2,450.,52.)
               call interior_style(fildes,INT_SOLID,1) !mark the chosen button
               if((sig_only.eq.'n').and.(flats_only.eq.'n'))then !mark the ALL button
                  call rectangle(fildes,210.,100.,240.,130.)
               else if(sig_only.eq.'y')then !mark the SIG FEATURES ONLY button
                  call rectangle(fildes,335.,100.,365.,130.)
               elseif(flats_only.eq.'y')then !mark the FLATS ONLY button
                  call rectangle(fildes,460.,100.,490.,130.)
               end if
               call interior_style(fildes,INT_HOLLOW,1)
            end if
            textstring='CHOOSE A SPEED:'
            call text2d(fildes,550.,60.,textstring
     +           //char(0),ANNOTATION_TEXT,0)
            call draw_button (fildes,1100.,90.,1150.,140.,fast,1120.,70.
     +           ,blank30,0.,0.)
            call draw_button (fildes,1000.,90.,1050.,140.,medium,1020.,   
     +           70.,blank30,0.,0.)
            call draw_button (fildes,900.,90.,950.,140.,slow,915.,70.
     +           ,blank30,0.,0.)                                          
            call draw_button (fildes,800.,90.,850.,140.,
     +           user_defined1,790.,70.,user_defined2,775.,52.)
            call text_color(fildes,0.,.5,0.) !green
            call perimeter_color(fildes,0.,.5,0.)
            call draw_button (fildes,1400.,110.,1450.,160.,go,1410.,90.
     +           ,blank30,0.,0.)                                        
            call text_color(fildes,1.,0.,0.) !red
            call perimeter_color(fildes,1.,0.,0.)
            textstring='CANCEL'
            call draw_button(fildes,1400.,20.,1450.,70.,                  
     +           textstring//char(0),1390.,1.,blank30,0.,0.)
            call text_color(fildes,0.,0.,0.)                    
            call perimeter_color(fildes,0.,0.,0.)
            call interior_style(fildes,INT_SOLID,1) !mark the chosen button
            if(speed.eq.10.)then
               call rectangle(fildes,910.,100.,940.,130.)
            else if(speed.eq.5.)then
               call rectangle(fildes,1010.,100.,1040.,130.)
            else if(speed.eq.3.)then
               call rectangle(fildes,1110.,100.,1140.,130.)
            else 
               call rectangle(fildes,810.,100.,840.,130.)
            end if
            call interior_style(fildes,INT_HOLLOW,1)
            call make_picture_current(fildes)
            call request_locator
     +           (fildes,1,timeout,valid,x_dc,y_dc,z_dc)
            call vdc_to_wc(fildes,x_dc,y_dc,z_dc,
     +           xloc,yloc,zloc) !translate into fp coords.
            screen = 'speed'
            call locate_region2(xloc,yloc,regionSPEED)
            CHOOSE_THE_SPEED: if(.false.)then !wait for a valid click
            else if(regionSPEED.eq.65)then !USER-DEFINED DISPLAY SPEED
 551           print '(/,T10,''Enter display time in seconds ''
     +''(minimum = 0.2 sec)  >> '',$)'
               speed = 0.0
               read (*,fmt='(f5.2)',err=551)speed
               if(speed.eq.0.0)return
               if(speed.lt.0.2)then
                  print '(/,T10,''Sorry -- minimum display time is''
     +'' 0.2 sec.  Please re-enter.'')'
                  goto 551
               end if
            else if(regionSPEED.eq.7)then !SLOW SPEED=10sec.
               speed = 10.
               goto 550
            else if(regionSPEED.eq.8)then !MEDIUM SPEED=5sec.
               speed = 5.
               goto 550
            else if(regionSPEED.eq.9)then !FAST SPEED=3sec.
               speed = 3.
               goto 550
            else if((regionSPEED.eq.92).and.
     +              ((mode.eq.'jl').or.(mode.eq.'ed')))then !view ALL CCHs
               sig_only='n'
               flats_only = 'n'
               goto 550
            else if((regionSPEED.eq.93).and.
     +              ((mode.eq.'jl').or.(mode.eq.'ed')))then !view CCHs with SIGNIFICANT FEATURES
               sig_only='y'
               flats_only = 'n'
               goto 550
            else if((regionSPEED.eq.931).and.
     +              ((mode.eq.'jl').or.(mode.eq.'ed')))then !view CCHs with FLAT FEATURES
               flats_only='y'
               sig_only = 'n'
               goto 550
            else if(regionSPEED.eq.91)then !GO
               click = 'n'      !clear the flag
               button_choice='auto_view'
               auto_mode = 'y'
               goto 800         !proceed with the slide show
            else if(regionSPEED.eq.94)then !cancel out of autoview mode
               button_choice = 'redraw'
               return
            else
               goto 550         !loop until an appropriate click is made
            end if CHOOSE_THE_SPEED
            goto 550            !loop until GO button is clicked
*
         else if(regionMAIN.eq.14)then !ENTER ANALYSIS RESULTS
            if((mode.eq.'cr').or.(mode.eq.'ed'))then
               call ENTER_ANALYSIS_RESULTS(fildes,mouse,prim,sec,loc1,
     +              loc2,rescom,Q_pos,statcomm,cdet)
            else
               goto 499
            end if
            button_choice='leave_stats'
            return                     
*
         else if(regionMAIN.eq.21.and.qdt_file_version.lt.'8')then !add or remove the SINGLE SHIFT CONTROL PLOT
            show_avg_shift = 'n'
            if(show_controlMAIN.eq.'y')then !currently displaying one of the shifts
               if(show_single_shift.eq.'y')then !currently displaying single shift
                  show_controlMAIN = 'n'
                  show_single_shift = 'n'
                  show_conf_lim = 'n'
               else             !currently displaying averaged shift
                  show_single_shift = 'y'
               end if
            else                !currently not displaying any shift at all
               show_controlMAIN = 'y'
               show_single_shift = 'y'
               show_conf_lim = 'n'
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.121)then !add or remove the AVERAGED SHIFT CONTROL PLOT
            if(avg_shift.ne.'y')goto 499 !can't add or remove something that's not available!
            show_single_shift = 'n'
            if(show_controlMAIN.eq.'y')then
               if(show_avg_shift.eq.'y')then !currently displaying averaged shift
                  show_controlMAIN = 'n'
                  show_avg_shift = 'n'
                  show_conf_lim = 'n'
               else             !currently displaying single shift
                  show_avg_shift = 'y'
               end if
            else                !currently not displaying any shift
               show_controlMAIN = 'y'
               show_avg_shift = 'y'
               show_conf_lim = 'n'
            end if
            button_choice = 'redraw'
            return
         else if(regionMAIN.eq.206)then !add or remove the SURROGATE CONTROL PLOT
            call sc_control (-1)
            button_choice = 'redraw'
            return
         else if(regionMAIN.eq.211)then !add or remove the SURROGATE CONTROL CL PLOT
            call sc_cl (-1)
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.20)then !REDRAW the main CCH window
            button_choice='redraw'
            return

         else if(regionMAIN.eq.18)then !redraw with CONFIDENCE LIMITS (mean +/- 2 std devs)
            show_3_sd = 'n'
            if(show_conf_lim.eq.'y')then
               if(show_2_sd.eq.'y')then !currently displaying 2 s.d.
                  show_conf_lim = 'n'
                  show_2_sd = 'n'
               else             !currently displaying 3 s.d.
                  show_2_sd = 'y'
               end if
            else                !currently not displaying any confidence limits
               show_conf_lim = 'y'
               show_2_sd = 'y'
            end if
            button_choice = 'redraw'
            return

         else if(regionMAIN.eq.118)then !redraw with CONFIDENCE LIMITS (mean +/- 3 std devs)
            show_2_sd = 'n'
            if(show_conf_lim.eq.'y')then
               if(show_3_sd.eq.'y')then !currently displaying 3 s.d.
                  show_conf_lim = 'n'
                  show_3_sd = 'n'
               else             !currently displaying 2 s.d.
                  show_3_sd = 'y'
               end if
            else                !currently not displaying any confidence limits
               show_conf_lim = 'y'
               show_3_sd = 'y'
            end if
            button_choice = 'redraw'
            return


         else
            goto 499
*
         end if !MAIN_CCH_WINDOW
*
*       ***** If the user has viewed a CCH screen and not made a point of entering *****
*       *****   analysis results, this program assumes that the CCH was FLAT.  *****
*       ***** But first, BE SURE THAT DATA FOR THIS CCH HAVE NOT BEEN ENTERED PREVIOUSLY: *****
*
 800     if(prim.eq.' ')then
            prim = 'Flat'       !load analysis database fields with the defaults
            sec = 'None'
            loc1 = ' '
            loc2 = ' '
            statcomm = ' '
            cdet=' '
            det = 0.0           !reset all statistical variables
            ZK=0.0
            Probk=0.0
            vis=0.0
            ZLAT=0.0
            HALFWD=0.0
         end if
*
*
      end if
*
*
*
*
      RETURN
      END
*
*
*
*
*       ********************************************************************************
*       ********************************************************************************
*
*               SUBROUTINE ENTER_ANALYSIS_RESULTS
*
*       ********************************************************************************
*       ********************************************************************************
*
*
      subroutine ENTER_ANALYSIS_RESULTS (fildes,mouse,prim,sec,
     +     loc1,loc2,rescom,Q_pos,statcomm,cdet)
*
*
      use mod_clear_routines
      use mod_locate_region1
      use mod_new_draw_button
      INCLUDE 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
*
      integer*4 mouse
      integer regionRESULTS,feature(21),Q_pos
*
      real x_dc,y_dc,z_dc,xloc,yloc,zloc,
     +     x0(21),x1(21),y0(21),y1(21),
     +     x_label(21),y_label(21)
*
      character*120 textstring
*
      character*40 rescom,statcomm
*
      character*30 label(21),perimeter,
     +     label10,blank30
*
      character*20 prim,sec
      character*15 loc1,loc2,screen
      character*8 cdet
*
*
*
*       *****  set screen parameters  ***************************************
*        call mapping_mode(fildes,1)
*        call view_window(fildes,1.,1.,1500.,800.)
*       *********************************************************************
*
*
      DATA x0/4*75.,4*400.,4*800.,4*1150.,5.,275.,500.,650.,800./
      DATA y0/90.,70.,50.,30.,90.,70.,50.,30.,90.,70.,50.,30.,
     +     90.,70.,50.,30.,2*130.,3*5./
      DATA x1/4*90.,4*415.,4*815.,4*1165.,205.,475.,600.,750.,900./
      DATA y1/105.,85.,65.,45.,105.,85.,65.,45.,105.,85.,65.,45.,
     +     105.,85.,65.,45.,2*165.,3*28./
      DATA x_label/4*95.,4*420.,4*820.,4*1170.,75.,300.,515.,655.,815/
      DATA y_label/90.,70.,50.,30.,90.,70.,50.,30.,90.,70.,50.,30.,
     +     90.,70.,50.,30.,2*140.,3*10./
*
*
*       ***** suppress unused variable warnings *****
      if(.false.)print *,mouse
      if(.false.)print *,fildes2 !suppress unused variable warning
      if(.false.)print *,Q_pos

      blank30 = ' '
      timeout = 604800.0        !timeout for waiting for a mouse click = 7 days
      do i = 1, 21
         feature(i)=0           !initialize arrays
      end do
      call fill_color(fildes,0.,0.,0.)                
*
*
      perimeter = 'perimeter'
      label(1) = 'PEAK'         !define button labels:
      label(2) = 'TROUGH'
      label(3) = 'MULT. PK & TR'
      label(4) = 'SKIP'
      label(5) = 'CENTRAL'
      label(6) = 'RIGHT'
      label(7) = 'LEFT'
      label(8) = 'SKIP'
      label(9) = 'NONE'
      label(10) = 'PEAK'  
      label(11) = 'TROUGH'
      label(12) = 'MULT. PK & TR'
      label(13) = 'CENTRAL'
      label(14) = 'RIGHT'
      label(15) = 'LEFT'
      label(16) = 'SKIP'
      label(17) = 'FLAT'
      label(18) = 'SKIP ALL'
      label(19) = 'OK'
      label(20) = 'CLEAR'
      label(21) = 'COMM:'
*
*                       ***** Draw analysis entry boxes: *****
*
*
*
 100  if(.false.)then           !see if any of the boxes should be already checked
      else if(prim.eq.'Peak')then
         feature(1)=1
      else if(prim.eq.'Trough')then
         feature(2)=1
      else if(prim.eq.'M P & T')then
         feature(3)=1
      else if(prim.eq.'Skip')then
         feature(4)=1
         feature(18)=1
      else if(prim.eq.'Flat')then
         feature(17)=1
      end if
*
      if(.false.)then
      else if(loc1.eq.'Central')then
         feature(5)=1
      else if(loc1.eq.'Offset right')then
         feature(6)=1
      else if(loc1.eq.'Offset left')then
         feature(7)=1
      else if(loc1.eq.'Skip')then
         feature(8)=1
      end if
*
      if(.false.)then
      else if(sec.eq.'None')then
         feature(9)=1
      else if(sec.eq.'Peak')then
         feature(10)=1
      else if(sec.eq.'Trough')then
         feature(11)=1
      else if(sec.eq.'M P & T')then
         feature(12)=1
      end if
*
      if(.false.)then
      else if(loc2.eq.'Central')then
         feature(13)=1
      else if(loc2.eq.'Offset right')then
         feature(14)=1
      else if(loc2.eq.'Offset left')then
         feature(15)=1
      else if(loc2.eq.'Skip')then
         feature(16)=1
      end if
*
*
 101  call clear_bottom(fildes) !clear bottom of window
      textstring='Choose all that apply to this pair'
      call text2d(fildes,5.,5.,textstring
     +     //char(0),ANNOTATION_TEXT,0)
*
      do i = 1,21               !draw results choices - fill in those that
         label10=label(i)
         call draw_button(fildes,x0(i),y0(i),x1(i),y1(i),label10, !have been selected
     +        x_label(i),y_label(i),blank30,0.,0.)
         if(feature(i).eq.1)then
            call interior_style(fildes,INT_SOLID,1)
            call rectangle(fildes,(x0(i)+5.),(y0(i)+5.),
     +           (x1(i)-5.),(y1(i)-5.))
            call interior_style(fildes,INT_HOLLOW,1)
            if(.false.)then
            else if(i.eq.17.or.i.eq.18)then
               call text_color(fildes,1.0,1.0,0.0) !black text
               call text2d(fildes,x_label(i),y_label(i),
     +              label(i)
     +              //char(0),ANNOTATION_TEXT,0)
               call text_color(fildes,0.,0.,0.) !back to black text
            end if
         end if
      end do
*
      textstring='PRIMARY FEATURE'
      call text2d(fildes,60.,110.,textstring
     +     //char(0),ANNOTATION_TEXT,0)
*
      textstring='LOCATION'
      call text2d(fildes,385.,110.,textstring
     +     //char(0),ANNOTATION_TEXT,0)
*
      textstring='SECONDARY FEATURE'
      call text2d(fildes,785.,110.,textstring
     +     //char(0),ANNOTATION_TEXT,0)
*
      textstring='LOCATION'
      call text2d(fildes,1135.,110.,textstring
     +     //char(0),ANNOTATION_TEXT,0)
*
      if(rescom.ne.' ')then
         call text2d(fildes,925.,10.,rescom
     +        //char(0),ANNOTATION_TEXT,0)
      end if
*
      call draw_button(fildes,1.,1.,1500.,170.,perimeter,0.,0.
     +     ,blank30,0.,0.)
*
*               *** wait for user to select a control button: ***
*
 110  call make_picture_current(fildes)
      call request_locator(fildes,1,timeout,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
      screen = 'results'
      call locate_region1(xloc,yloc,regionRESULTS)
*
*
      if(.false.)then           !ENTER_RESULTS_OF_ANALYSIS
*                       !no valid button chosen yet -- sample mouse again
      else if(regionRESULTS.eq.16 !primary feature button chosen
     +        .or.regionRESULTS.eq.17
     +        .or.regionRESULTS.eq.18
     +        .or.regionRESULTS.eq.19)then
         prim = ' '             !reset the value
         feature(1:4) = 0
         feature(17:18) = 0
*
         if(regionRESULTS.eq.16)then
            prim = 'Peak'
            feature(1)=1
         end if
*
         if(regionRESULTS.eq.17)then
            prim = 'Trough'
            feature(2)=1
         end if
*       
         if(regionRESULTS.eq.18)then
            prim = 'M P & T'
            feature(3)=1
         end if
*       
         if(regionRESULTS.eq.19)then
            prim = 'Skip'
            feature(4)=1
         end if
*
      else if(regionRESULTS.eq.20 !location of primary feature button chosen
     +        .or.regionRESULTS.eq.21
     +        .or.regionRESULTS.eq.22
     +        .or.regionRESULTS.eq.23)then
         loc1 = ' '             !reset the value
         feature(5:8)=0
*       
         if(regionRESULTS.eq.20)then
            loc1 = 'Central'
            feature(5)=1
         end if
*       
         if(regionRESULTS.eq.21)then
            loc1 = 'Offset right'
            feature(6)=1
         end if
*       
         if(regionRESULTS.eq.22)then
            loc1 = 'Offset left'
            feature(7)=1
         end if
*       
         if(regionRESULTS.eq.23)then
            loc1 = ' '
            feature(8)=1
         end if
*       
      else if(regionRESULTS.eq.24 !secondary feature button chosen
     +        .or.regionRESULTS.eq.25
     +        .or.regionRESULTS.eq.26
     +        .or.regionRESULTS.eq.27)then
         sec = ' '              !reset the value
         feature(9:12)=0
*       
         if(regionRESULTS.eq.24)then
            sec = 'None'
            feature(9)=1
         end if
*       
         if(regionRESULTS.eq.25)then
            sec = 'Peak'
            feature(10)=1
         end if
*       
         if(regionRESULTS.eq.26)then
            sec = 'Trough'
            feature(11)=1
         end if
*       
         if(regionRESULTS.eq.27)then
            sec = 'M P & T'
            feature(12)=1
         end if
*       
      else if(regionRESULTS.eq.28 !location of secondary feature button chosen
     +        .or.regionRESULTS.eq.29
     +        .or.regionRESULTS.eq.30
     +        .or.regionRESULTS.eq.31)then
         loc2 = ' '             !reset the value
         feature(13:16)=0
*       
         if(regionRESULTS.eq.28)then
            loc2 = 'Central'
            feature(13)=1
         end if
*
         if(regionRESULTS.eq.29)then
            loc2 = 'Offset right'
            feature(14)=1
         end if
*       
         if(regionRESULTS.eq.30)then
            loc2 = 'Offset left'
            feature(15)=1
         end if
*       
         if(regionRESULTS.eq.31)then
            loc2 = ' '
            feature(16)=1
         end if
*       
      else if(regionRESULTS.eq.32)then !FLAT
         feature(1:20)=0        !clear all results flags
         prim = 'Flat'
         loc1 = ' '
         sec = 'None'
         loc2 = ' '
         feature(17)=1
*
      else if(regionRESULTS.eq.33)then !SKIP ALL
         feature(1:20)=0        !clear all results flags
         prim = 'Skip'
         loc1 = ' '
         sec = ' '
         loc2 = ' '
         feature(18)=1
         feature(4)=1
*       
      else if(regionRESULTS.eq.34)then !OK --> write to *.dbp file
*
         if((prim.eq.' ').or.(prim.eq.'Skip').or.
     +        (loc1.eq.'Skip').or.(loc2.eq.'Skip'))then !empty data field for primary feature is not allowed
            call clear_bottom(fildes) !clear the bottom of the screen
            textstring='You must select an option for primary '//
     +           'feature.  No skipping allowed!'
            call text2d(fildes,50.,50.,textstring
     +           //char(0),ANNOTATION_TEXT,0)
            textstring='Click to continue.'
            call text2d(fildes,80.,30.,textstring
     +           //char(0),ANNOTATION_TEXT,0)
*
            call make_picture_current(fildes)
            call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
            goto 100            !redraw all analysis choices and start over
         end if
*
         return                 !go back to showCCHs and write to *.pre_dbp file
*
      else if(regionRESULTS.eq.35)then !CLEAR all entries and start over
         feature(1:20)=0        !clear all results flags
         prim = ' '
         loc1 = ' '
         sec = ' '
         loc2 = ' '
         rescom = ' '
         statcomm = ' '
         cdet = ' '
         goto 100
*
      else if(regionRESULTS.eq.36)then !enter COMMENT
         rescom = ' '           !clear a previous entry
 120     print '(5x,''CLICK HERE'',
     +/,5x,''Comment (no commas) >> '',$)'
         read (*,'(A)',err=120)rescom
         do i = 1,LEN(rescom)
            if(rescom(i:i).eq.',')rescom(i:i)=';'
         end do
*
      else
         goto 110
*
      end if     !ENTER_RESULTS_OF_ANALYSIS !end of analysis entry section
*
      goto 101                  !allow user to choose analysis results until OK is chosen
*
*
c        return
      end


      end module mod_showCCHs6_newshift
