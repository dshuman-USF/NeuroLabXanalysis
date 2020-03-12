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

      module mod_print_and_write_routines
      contains
*     filename: print_and_write_routines.f  (created 02-Feb-2004         lss)

*     This file contains subroutines for printing or saving graphics widows in x2004.
*     To do this, window displays are re-written to a postscript file, and printed or saved from there.

*     Subroutines included in this file:        print_or_write_ENLARGED_CCH
*                                               print_or_write_DIFF_CCH
*                                               print_or_write_ALL_CCHs
*                                               print_or_write_CTHs
*                                               print_or_write_QSUM
*                                               print_or_write_STATS

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or save the ENLARGED CCH

*     date of last revision: 16-Mar-2004        lss


        subroutine print_or_write_ENLARGED_CCH(task,CCH,CONTROL,
     +                       BDT_FILE,                               !print the ENLARGED CCH
     +                       QDT_FILE,REF,TAR,date,recording,NHW,
     +                       rec_num,rec_num_control,title,ICN,
     +                       show_control,show_conf_lim,
     +                       REFtype,TARtype,per_results_REF,
     +                       per_results_TAR,per_text_abbrev,
     +                       show_single_shift,show_avg_shift,
     +                       CONTROL_cl,show_2_sd,show_3_sd,prim_txt,
     +                       min_bin,max_bin,scaledup,rmagnify,IDs,ITAL,
     +                       REF_coords,TAR_coords,prim_txt_stats)
      use mod_miscellaneous_subroutines
      use mod_ps_plot_1

        include 'x2000parameter.defs'

        integer*4 isys,SYSTEM
        integer CCH(101),CONTROL(101),scaled_hist(101),
     +       rec_num,rec_num_control,ICN(6),REF,
     +       TAR,CONTROL_cl(101),IDs(MAX_NUM_CODES),
     +       ITAL(MAX_NUM_CHAN)
         
        character*120 title,text,prim_txt,prim_txt_stats
        character*100 ps_file
        character*(*) BDT_FILE,QDT_FILE,REF_coords(3),TAR_coords(3)
        character*15 REFtype,TARtype
        character*12 per_text_abbrev(MAX_PERTURB)
        character*11 date
        character*10 bwtext
        character*9 per_results_REF(MAX_PERTURB),
     +              per_results_TAR(MAX_PERTURB)
        character*2 recording
        character*3 CCHtype,REFtext,TARtext
        character*1 info,show_control,task,show_conf_lim,
     +              show_single_shift,show_avg_shift,
     +              show_2_sd,show_3_sd,scaledup

*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
                 call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
                 call ps_mapping_mode(1)
                 call ps_view_port(0.,0.,1.,1.)
                 call ps_view_window(1.,1.,1500.,800.)
                 call ps_geometry(18.,18.,792.-18.,612.-18.)
*       *********************************************************************

        write (bwtext,'(f7.1)') (float(NHW)/50.)
        write (REFtext,'(I3)') REF
        write (TARtext,'(I3)') TAR
        CCHtype = 'CCH'

        if(task.eq.'w')then
           call strlength(QDT_FILE,LEN(QDT_FILE),l)
           ps_file=QDT_FILE(1:l)//'_'//REFtext//'x'//TARtext//'_'//
     +          bwtext//'.ps'
           if(scaledup.eq.'y')
     +          ps_file=QDT_FILE(1:l)//'_'//REFtext//'x'//TARtext//'_'//
     +          bwtext//'_scale.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if
        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))

                 info='y'
              call ps_character_height(.036)
              call ps_character_width(.012)
              x0 = 400.
              y0 = 300.
              call ps_plot(101,CCH,CONTROL,scaled_hist,CCHtype,        !plot CCH 
     +                      x0,y0,300.,525.,info,BDT_FILE,
     +                      QDT_FILE,date,recording,NHW,
     +                      rec_num,show_control,
     +                      rec_num_control,title,ICN,max_bin,
     +                      show_conf_lim,
     +                      show_single_shift,show_avg_shift,
     +                      CONTROL_cl,show_2_sd,show_3_sd,
     +                      scaledup,min_bin,rmagnify,IDs,ITAL)

              call ps_text2d(1000.,700.,REFtext//char(0),
     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(1175.,700.,TARtext//char(0),
     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(1325.,700.,'ID code'//char(0),
     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(1000.,680.,REFtype//char(0),
     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(1175.,680.,TARtype//char(0),
     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(1325.,680.,'resp. type'//char(0),
     +                  PS_ANNOTATION_TEXT)

              text = ' '
              call strlength(REF_coords(1),LEN(REF_coords(1)),l_1)
              call strlength(REF_coords(2),LEN(REF_coords(2)),l_2)
              call strlength(REF_coords(3),LEN(REF_coords(3)),l_3)
              text = REF_coords(1)(1:l_1)//', '//
     +                  REF_coords(2)(1:l_2)//', '//
     +                  REF_coords(3)(1:l_3)
              call strlength(text,LEN(text),l)
              call ps_text2d(950.,300.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              text = ' '
              call strlength(TAR_coords(1),LEN(TAR_coords(1)),l_1)
              call strlength(TAR_coords(2),LEN(TAR_coords(2)),l_2)
              call strlength(TAR_coords(3),LEN(TAR_coords(3)),l_3)
              text = TAR_coords(1)(1:l_1)//', '//
     +                  TAR_coords(2)(1:l_2)//', '//
     +                  TAR_coords(3)(1:l_3)
              call strlength(text,LEN(text),l)
              call ps_text2d(1175.,300.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)

              if(show_single_shift.eq.'y')then
                 text='(single shift-control displayed)'
              else if(show_avg_shift.eq.'y')then
                 text='(averaged shift-control displayed)'
                 if(qdt_file_version.eq.'8')
     +                text='(cth control displayed)'
              else
                 text = ' '
              end if
              call strlength(text,LEN(text),l)
              call ps_text2d(x0+220.,y0-70.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              if(show_conf_lim.eq.'y')then
                 if(show_2_sd.eq.'y')then
                    text='(confidence limits:  +/- 2 s.d.)'
                 else if(show_3_sd.eq.'y')then
                    text='(confidence limits:  +/- 3 s.d.)'
                 else
                    text = ' '
                 end if
                 call strlength(text,LEN(text),l)
                 call ps_text2d(x0+220.,y0-90.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              end if

              call sc_label2_ps (fildes, x0+220, y0-70,show_single_shift
     +             ,show_avg_shift,show_2_sd,show_3_sd)

              call ps_character_height(.036)
              call ps_character_width(.012)
              
              call strlength(prim_txt,LEN(prim_txt),l)
c              call ps_text2d(x0,y0-150.,prim_txt(1:l)//char(0),
c     +                  PS_ANNOTATION_TEXT)
              call ps_text2d(x0,y0-150.,prim_txt(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              call strlength(prim_txt_stats,LEN(prim_txt_stats),l)
              call ps_text2d(x0+20,y0-170.,prim_txt_stats(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)

              y=670.
              do i = 1, MAX_PERTURB
                 if((per_results_REF(i).ne.' ').or.
     +              (per_results_TAR(i).ne.' '))then
                    y = y - 20.
                    call strlength(per_results_REF(i),
     +                   LEN(per_results_REF(i)),l)
                    call ps_text2d(1000.,y,per_results_REF(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
                    call strlength(per_results_TAR(i),
     +                   LEN(per_results_TAR(i)),l)
                    call ps_text2d(1175.,y,per_results_TAR(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
                    call strlength(per_text_abbrev(i),
     +                   LEN(per_text_abbrev(i)),l)
                    call ps_text2d(1325.,y,per_text_abbrev(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
                 end if
              end do
              call strlength(ps_file,LEN(ps_file),l)
              call ps_fclose(ps_file(1:l))
              if((task.eq.'p').or.(task.eq.'v'))then
                 if(task.eq.'p')then
                    isys=SYSTEM('lpr '//                !print a hard copy of the ps file
     +                  ps_file(1:l)//char(0)) 
                 else
                    isys=SYSTEM('gv temp.ps'//char(0))          !pop-up display of the ps file
                 end if
                 isys=SYSTEM('rm temp.ps'//char(0))             !remove the temporary file
              end if

              return
              end





*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or save the DIFFERENCE CCH

*     date of last revision: 02-Feb-2004        lss


        subroutine print_or_write_DIFF_CCH(task,CCH,CONTROL,BDT_FILE,  !re-plot the DIFFERENCE HISTOGRAM
     +                       QDT_FILE,REF,TAR,date,recording,NHW,
     +                       rec_num,rec_num_control,title,ICN,                            
     +                       min_bin,max_bin,scaledup,rmagnify,
     +                       IDs,ITAL)
      use mod_miscellaneous_subroutines
      use mod_ps_plot_1

        include 'x2000parameter.defs'

        integer*4 isys,SYSTEM
        integer CCH(101),CONTROL(101),DIFF(101),scaled_hist(101),
     +       rec_num,rec_num_control,ICN(6),REF,
     +       TAR,IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
        

        character*120 title
        character*100 ps_file
        character*(*) BDT_FILE,QDT_FILE
        character*11 date
        character*10 bwtext
        character*2 recording
        character*3 CCHtype,DIF,REFtext,TARtext
        character*1 info,show_control_DIFF,task,scaledup

*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
                 call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
                 call ps_mapping_mode(1)
                 call ps_view_port(0.,0.,1.,1.)
                 call ps_view_window(1.,1.,1500.,800.)
                 call ps_geometry(18.,18.,792.-18.,612.-18.)
*       *********************************************************************

        write (bwtext,'(f7.1)') (float(NHW)/50.)
        write (REFtext,'(I3)') REF
        write (TARtext,'(I3)') TAR

        if(task.eq.'w')then
           call strlength(QDT_FILE,LEN(QDT_FILE),l)
           ps_file=QDT_FILE(1:l)//'_'//REFtext//'x'//TARtext//'_DIF_'//
     +          bwtext//'.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if
        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))

        DIFF = CCH - CONTROL
        info='y'
        show_control_DIFF='y'
        CCHtype = 'CCH'
        DIF = 'DIF'
        call ps_plot(101,CCH,CONTROL,scaled_hist,CCHtype,200.,300.,300.,!plot CCH with shift-control
     +       525.,info,BDT_FILE,QDT_FILE,date,recording,NHW,    
     +       rec_num,show_control_DIFF,rec_num_control,title,   
     +       ICN,max_bin,'','','',[integer::],'','',scaledup,min_bin,
     +       rmagnify,IDs,ITAL)                                          
        info = 'n'
        call ps_plot(101,DIFF,DIFF,scaled_hist,DIF,900.,300.,300.,525.,!draw the difference histogram
     +       info,BDT_FILE,QDT_FILE,date,recording,NHW,rec_num, 
     +       show_control_DIFF,rec_num_control,title,ICN, max_bin,'','',
     +       '',[integer::],'','',scaledup,min_bin,rmagnify,IDs, ITAL)                                              

        call strlength(ps_file,LEN(ps_file),l)
        call ps_fclose(ps_file(1:l)//char(0))

        if((task.eq.'p').or.(task.eq.'v'))then
           if(task.eq.'p')then
              isys=SYSTEM('lpr temp.ps'//char(0))
           else
              isys=SYSTEM('gv temp.ps'//char(0))
           end if
           isys=SYSTEM('rm temp.ps'//char(0))
        end if
        
        return
        end

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or save the MAIN CCH window containing CCHs and ACHs at 4 binwidths,
*       as well as respiratory and cardiac CTHs for both cells of the pair

*     date of last revision: 17-April-2006        lss


        subroutine print_or_write_ALL_CCHs(task,date,recording,
     +     BDT_FILE,QDT_FILE,REFcode,TARcode,
     +     REFname,TARname,REFtype,TARtype,
     +     CCHs,CONTROLs,REF_ACHs,TAR_ACHs,BINWtexts,
     +     CTH_REF,CTH_TAR,CTH_phrenic,ETA2_REF,ETA2_TAR,zmodsig_REF,
     +     zmodsig_TAR,zmodsig2_REF,zmodsig2_TAR,cardiac_REF,
     +     cardiac_TAR,cardiac_hist,card_type_REF,card_type_TAR,
     +     prim,sec,loc1,loc2,rescom,statcomm,
     +     cdet,show_controlMAIN,c_format,included,total_num_cells,
     +     DELTA2_REF,DELTA2_TAR,tedfactor_REF,tedfactor_TAR,
     +     show_single_shift,show_avg_shift,show_conf_lim,show_2_sd,
     +     show_3_sd,CONTROLs_cl,min_bin,max_bin,scaled_up,rmagnify,
     +     IDs,ITAL,cardiac_pls,REF_coords,TAR_coords,scaledMAIN,
     +     show_pulses,mean_E,NORM_OFFSET,bw_n)
      use mod_miscellaneous_subroutines
      use mod_ps_plot_1

        include 'x2000parameter.defs'

        integer*4 isys,SYSTEM,subtract_1
        integer CTH_REF(NUM_BINS),CTH_TAR(NUM_BINS),
     +          CTH_phrenic(NUM_BINS),
     +          cardiac_REF(NUM_BINS),cardiac_TAR(NUM_BINS),
     +          cardiac_hist(NUM_BINS),
     +          CCHs(4,NUM_BINS),CONTROLs(4,NUM_BINS),
     +          CONTROLs_cl(4,NUM_BINS),
     +          CCHtemp(NUM_BINS),CONTROLtemp(NUM_BINS),
     +          CONTROL_cltemp(NUM_BINS),
     +          REF_ACHs(4,NUM_BINS),TAR_ACHs(4,NUM_BINS),
     +          scaled_hist(NUM_BINS),included(MAX_NUM_CODES),
     +          IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),cardiac_pls

        integer REFcode,TARcode,total_num_cells,subtract

        real mean_E,NORM_OFFSET,bw_n

        character*200 REFline,TARline
        character*120 title,text,label1,label1a
        character*100 ps_file
        character*40 rescom,statcomm
        character*30 label4
        character*20 prim,sec,c_format
        character*15 loc1,loc2,REFtype,TARtype
        character*5  REF_coords(3),TAR_coords(3)
        character*(*) BDT_FILE,QDT_FILE
        character*30 label2,label3
        character*10 c_ITAL
        character*11 date
        character*12 BINwtexts(4)
        character*8 cdet
        character*5 ETA2_REF,ETA2_TAR,DELTA2_REF,DELTA2_TAR
        character*4 REFcodetext,TARcodetext,REFname,TARname
        character*3 CCHtype,CTH,card_type_REF,
     +              card_type_TAR,zmodsig_REF,zmodsig_TAR,zmodsig2_REF,
     +              zmodsig2_TAR,ACH
        character*2 recording,tedfactor_REF,tedfactor_TAR
        character*1 info,show_controlMAIN,task,show_control,
     +                  show_single_shift,show_avg_shift,show_conf_lim,
     +                  show_2_sd,show_3_sd,scaled_up,scaledMAIN,
     +                  show_pulses,ign1
        integer*4 ICN0(6)

*       ***** suppress unused variable warnings *****
        if(.false.)print *,c_format
        if(.false.)print *,included
        if(.false.)print *,total_num_cells

*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
        call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
        call ps_mapping_mode(1)
        call ps_view_port(0.,0.,1.,1.)
        call ps_view_window(1.,1.,1500.,800.)
        call ps_geometry(18.,18.,792.-18.,612.-18.)
*       *********************************************************************

        ICN0 = 0
        CCHtype = 'CCH'
        CTH = 'CTH'
        ACH = 'ACH'

        write (REFcodetext,'(I4)') REFcode
        write (TARcodetext,'(I4)') TARcode
        if(task.eq.'w')then
           call strlength(QDT_FILE,LEN(QDT_FILE),l)
           ps_file=QDT_FILE(1:l)//'_'//REFcodetext//'x'//TARcodetext//
     +             '_ALL_CCHs.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if
        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))


*     ***** draw to the ps_file: *****


        title=' '
        title = REFname//' > '//TARname//
     +    '    ('//REFcodetext//' > '//TARcodetext//')'
        call ps_character_height(.075)
        call ps_character_width(.025)
        call strlength(title,LEN(title),l)
        call ps_text2d(300.,760.,title(1:l)
     +       //char(0),PS_ANNOTATION_TEXT)                             !'REF>TAR'
        label1 = 'CCHs:'
        label2 = 'ACHs:'
        call strlength(REFcodetext,LEN(REFcodetext),ir)
        call strlength(TARcodetext,LEN(TARcodetext),it)
        label3='('//REFcodetext(1:ir)//','//TARcodetext(1:it)//')'
        call remove_all_blanks(label3,LEN(label3))
        label4 = 'CTHs:'
        info = 'n'
        show_control = 'n'

        DISPLAY_CCHs:   do izz=1,4
           call sc_izz (izz)

           if(.false.)then                                            !define (0,0) of each CCH
              else if(izz.eq.1)then
                 x=100
                 y=550
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

*
*       **************************************************
*       *****   plot the histograms for this pair:   *****
*       **************************************************
*
           call ps_character_height(.045)
           call ps_character_width(.015)
           call ps_text2d((x+100.),(y-155.),BINWtexts(izz)//char(0),
     +          PS_ANNOTATION_TEXT) 
        if(scaledMAIN.eq.'y')then
           subtract = 75
           min_bin = CCHs(izz,1)
           max_bin = CCHs(izz,1)
           do i = 2,101
              min_bin = MIN0(min_bin,CCHs(izz,i))   
              max_bin = MAX0(max_bin,CCHs(izz,i))
           end do
           subtract_1 = min_bin*(subtract/100.)
           call sc_subtract (subtract_1)
           do i = 1, 101
              CCHtemp(i) = CCHs(izz,i) - subtract_1
              CONTROLtemp(i) = CONTROLs(izz,i) - subtract_1
              CONTROL_cltemp(i) = CONTROLs_cl(izz,i) - subtract_1
           end do
        else
           call sc_subtract (0)
           do i = 1, 101
              CCHtemp(i) = CCHs(izz,i)
              CONTROLtemp(i) = CONTROLs(izz,i)
              CONTROL_cltemp(i) = CONTROLs_cl(izz,i)
           end do
        end if
           call ps_plot(101,CCHtemp,CONTROLtemp,scaled_hist,CCHtype,x,y,
     +       200.,303.,info,'','','','',0,0,show_controlMAIN,0,'',
     +       ICN0,0,show_conf_lim,'','',CONTROL_cltemp,show_2_sd, !plot the CCH  ...
     +       show_3_sd,ign1,0,0.,[integer::],[integer::])                                
c           call ps_plot(101,CCHs(izz,1:NUM_BINS),
c     +          CONTROLs(izz,1:NUM_BINS),
c     +          scaled_hist,CCHtype,x,y,200.,303.,                    !plot the CCH  ...
c     +          info,,,,,,,show_controlMAIN,,,,,show_conf_lim,,,
c     +          CONTROLs_cl(izz,1:NUM_BINS),show_2_sd,show_3_sd,,,)    
c           call ps_plot(100,REF_ACHs(izz,1:NUM_BINS),
c     +          REF_ACHs(izz,1:NUM_BINS),
c     +          scaled_hist,ACH,                                      !the reference ACH and ...
c     +          (x+25.),(y-125.),100.,101.,info,,,,,,,show_control,,,,
c     +          max_bin,,,,,,,scaled_up,min_bin,rmagnify)
           call ps_plot(50,REF_ACHs(izz,1:NUM_BINS),
     +          REF_ACHs(izz,1:NUM_BINS),scaled_hist,ACH,(x-2.),
     +          (y-125.),100.,150.,info,'','','','',0,0, !the reference ACH and ...
     +          show_control,0,'',ICN0,max_bin,'','','',
     +          [integer::],'','', ! (show the first 50 bins)
     +          scaled_up,min_bin,rmagnify,[integer::],[integer::])
c           call ps_plot(100,TAR_ACHs(izz,1:NUM_BINS),
c     +          TAR_ACHs(izz,1:NUM_BINS),
c     +          scaled_hist,ACH,                                      !the target ACH and ...
c     +          (x+152.+2.),(y-125.),100.,101.,info,,,,,,,show_control,,,,
c     +          max_bin,,,,,,,scaled_up,min_bin,rmagnify)
           call ps_plot(50,TAR_ACHs(izz,1:NUM_BINS),                    
     +          TAR_ACHs(izz,1:NUM_BINS),scaled_hist,ACH,       
     +          (x+152.+2.),(y-125.),100.,150.,info,'','','','',0,0, !the target ACH and ...
     +          show_control,0,'',ICN0,max_bin,'','','',
     +          [integer::],'','',scaled_up, ! (show the first 50 bins)
     +          min_bin,rmagnify,[integer::],[integer::])                               
        end do DISPLAY_CCHs

              info = 'n'
              show_control = 'y'
              call ps_plot(100,CTH_REF,CTH_phrenic,scaled_hist,CTH,100.,
     +             225.,150.,225.,info,'','','','',0,0,show_control, !the reference CTH and ...
     +             0,'',ICN0,max_bin,'','','',[integer::],'','',
     +             scaled_up,min_bin, rmagnify,[integer::],[integer::])
              call ps_plot(100,CTH_TAR,CTH_phrenic,scaled_hist,CTH,350.,
     +             225.,150.,225.,info,'','','','',0,0,show_control, !the target CTH
     +             0,'',ICN0,max_bin,'','','',[integer::],'','',
     +             scaled_up,min_bin, rmagnify,[integer::],[integer::])
              call ps_plot(100,CTH_REF,CTH_TAR,scaled_hist,CTH,600.,    
     +             225.,150.,225.,info,'','','','',0,0,show_control, !the superimposed CTHs
     +             0,'',ICN0,max_bin,'','','',[integer::],'','',
     +             scaled_up,min_bin, rmagnify,[integer::],[integer::])

        if(show_pulses.ne.'n')then
c          call ps_character_height(.030)
c          call ps_character_width(.010)
           call ps_character_height(.021)
           call ps_character_width(.007)
           call ps_move2d(100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.-5.)
           call ps_draw2d(100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.+5.)
           call ps_text2d(100.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          210.,'E'//char(0),PS_ANNOTATION_TEXT)

           call ps_move2d(350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.-5.)
           call ps_draw2d(350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.+5.)
           call ps_text2d(350.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          210.,'E'//char(0),PS_ANNOTATION_TEXT)

           call ps_move2d(600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.-5.)
           call ps_draw2d(600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          225.+5.)
           call ps_text2d(600.+(225.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          210.,'E'//char(0),PS_ANNOTATION_TEXT)
           if(mean_E.ne.0.0)then
              call ps_move2d(100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.-5.)
              call ps_draw2d(100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.+5.)
              call ps_text2d(100.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),210.,'I'//char(0),PS_ANNOTATION_TEXT)

              call ps_move2d(350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.-5.)
              call ps_draw2d(350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.+5.)
              call ps_text2d(350.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),210.,'I'//char(0),PS_ANNOTATION_TEXT)

              call ps_move2d(600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.-5.)
              call ps_draw2d(600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),225.+5.)
              call ps_text2d(600.+(225.*((ABS(NORM_OFFSET)+mean_E)/
     +          (100.*bw_n))),210.,'I'//char(0),PS_ANNOTATION_TEXT)
           end if
        end if



         if(cardiac_pls.ne.0)then
              call ps_plot(100,cardiac_REF,cardiac_hist,scaled_hist,CTH,
     +           950.,225.,150.,225.,info,'','','','',0,0, !the reference cardiac CCH and ...
     +           show_control,0,'',ICN0,max_bin,'','','',
     +           [integer::],'','', scaled_up,min_bin,rmagnify,IDs,ITAL)         
              call ps_plot(100,cardiac_TAR,cardiac_hist,scaled_hist,CTH,
     +             1200.,225.,150.,225.,info,'','','','',0,0, !the target cardiac CCH
     +             show_control,0,'',ICN0,max_bin,'','','',
     +             [integer::],'','',scaled_up,min_bin,rmagnify,IDs,
     +             ITAL)
           else
              call ps_character_height(.030)
              call ps_character_width(.010)
              text = 'CARDIAC CCHs NOT CALCULATED'
              call strlength(text,LEN(text),l)
              call ps_text2d(1000.,330.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
           end if
              
              
              call ps_character_height(.075)
              call ps_character_width(.025)
c              call strlength(text,LEN(text),l)
c              call ps_text2d (300.,760.,title(1:l)
c     +                //char(0),PS_ANNOTATION_TEXT)                   !'REF>TAR'
              call ps_character_height(.030)
              call ps_character_width(.010)
              call strlength(label3,LEN(label3),l)
              call ps_text2d(5.,430.,label3(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)                    !'(ref,tar)'
              call ps_text2d(10.,400.,'(1st 50 bins shown)'//char(0),
     +                  PS_ANNOTATION_TEXT)
              if(scaledMAIN.eq.'y')then
                 call ps_text2d(5.,675.,'scaled'//char(0),
     +                  PS_ANNOTATION_TEXT)
                 call ps_text2d(10.,660.,'up 75%'//char(0),
     +                  PS_ANNOTATION_TEXT)
              else
                 call ps_text2d(5.,675.,'original'//char(0),
     +                  PS_ANNOTATION_TEXT)
              end if

              call ps_character_height(.045)
              call ps_character_width(.015)
              call strlength(label1,LEN(label1),l)
              call ps_text2d (1.,710.,label1(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)                    !'CCHs:'
              call strlength(label2,LEN(label2),l)
              call ps_text2d (1.,475.,label2(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)                    !'ACHs:'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (1.,300.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)                    !'CTHs:'
              call ps_character_height(.036)
              call ps_character_width(.012)
              call ps_character_height(.030)
              call ps_character_width(.010)

              if((qdt_file_version.eq.'8')
     +             .and.(show_avg_shift.eq.'y')) then
               if(show_2_sd.eq.'y')then
                  text='(cth control +/- 2sd displayed)'
               else if(show_3_sd.eq.'y')then
                  text='(cth control +/- 3sd displayed)'
               else
                  text='(cth control displayed)'
               end if
              else if(show_single_shift.eq.'y')then
               if(show_2_sd.eq.'y')then
                  text='(single shift-control +/- 2sd displayed)'
               else if(show_3_sd.eq.'y')then
                  text='(single shift-control +/- 3sd displayed)'
               else
                  text='(single shift-control displayed)'
               end if
              else if(show_avg_shift.eq.'y')then
               if(show_2_sd.eq.'y')then
                  text='(averaged shift-control +/- 2sd displayed)'
               else if(show_3_sd.eq.'y')then
                  text='(averaged shift-control +/- 3sd displayed)'
               else
                  text='(averaged shift-control displayed)'
               end if
              else
                 text = ' '
              end if
              call sc_label_ps (fildes, 400., 530.,show_single_shift
     +             ,show_avg_shift,show_2_sd,show_3_sd)
              call ps_character_height(.030)
              call ps_character_width(.010)
              call strlength(text,LEN(text),l)
              call ps_text2d(10.,530.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)

              label4 = 'resp:'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (25.,375.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
              if (CTH_REF(101).gt.0) then
                 call ps_text2d(0.,350.,'(ph. norm.)'//char(0),
     +                PS_ANNOTATION_TEXT)
              else
                 call ps_text2d(15.,350.,'(norm.)'//char(0),
     +                PS_ANNOTATION_TEXT)
              end if
              label4 = 'cardiac:'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (860.,375.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
              call ps_character_height(.021)
              call ps_character_width(.007)
              label4 = 'ANOVA/BINARY:'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (1.,175.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
              label4 = 'ETA2:'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (1.,160.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)

              call ps_character_height(.030)
              call ps_character_width(.010)
              
              label4 = REFcodetext                                     !label the respiratory CTHs:
              call remove_all_blanks(REFcodetext,LEN(REFcodetext))
              call remove_all_blanks(REFtype,LEN(REFtype))
              call strlength(REFcodetext,LEN(REFcodetext),l)
              call strlength(REFtype,LEN(REFtype),m)
              label4 = REFcodetext(1:l)//' ('//REFtype(1:m)//')'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (150.,195.,label4(1:l)
     +             //char(0),PS_ANNOTATION_TEXT)
              call upper_case(zmodsig_REF,LEN(zmodsig_REF))
              call upper_case(zmodsig2_REF,LEN(zmodsig2_REF))
              label4=zmodsig_REF//'/ '//zmodsig2_REF
              call strlength(label4,LEN(label4),l)
              call ps_text2d(200.,175.,label4(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              label4 = ETA2_REF
              call strlength(label4,LEN(label4),l)
              call ps_text2d(150.,160.,label4(1:l)
     +             //char(0),PS_ANNOTATION_TEXT)
              
              label4=TARcodetext
              call remove_all_blanks(TARcodetext,LEN(TARcodetext))
              call remove_all_blanks(TARtype,LEN(TARtype))
              call strlength(TARcodetext,LEN(TARcodetext),l)
              call strlength(TARtype,LEN(TARtype),m)
              label4 = TARcodetext(1:l)//' ('//TARtype(1:m)//')'
              call strlength(label4,LEN(label4),l)
              call ps_text2d (400.,195.,label4(1:l)
     +             //char(0),PS_ANNOTATION_TEXT)
              call upper_case(zmodsig_TAR,LEN(zmodsig_TAR))
              call upper_case(zmodsig2_TAR,LEN(zmodsig2_TAR))
              label4=zmodsig_TAR//'/ '//zmodsig2_TAR
              call strlength(label4,LEN(label4),l)
              call ps_text2d(450.,175.,label4(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
              label4=ETA2_TAR
              call strlength(label4,LEN(label4),l)
              call ps_text2d(400.,160.,label4(1:l)
     +             //char(0),PS_ANNOTATION_TEXT)
              
              label4 = REFcodetext//' & '//TARcodetext
              call strlength(label4,LEN(label4),l)
              call ps_text2d(710.,195.,label4(1:l)
     +             //char(0),PS_ANNOTATION_TEXT)
              
              if(cardiac_pls.ne.0)then
                 label4=REFcodetext                             !label the cardicac CCCs
                 call strlength(label4,LEN(label4),l)
                 call ps_text2d (1060.,195.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
                 label4=TARcodetext
                 call strlength(label4,LEN(label4),l)
                 call ps_text2d (1310.,195.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
                 label4 = 'DELTA2:'
                 call strlength(label4,LEN(label4),l)
                 call ps_text2d (860.,175.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
                 label4 = '  (value; t.f.)'
                 call strlength(label4,LEN(label4),l)
                 call ps_text2d(860.,160.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
                 if(card_type_REF.ne.' ')then
                    call upper_case(card_type_REF,LEN(card_type_REF))
                    call ps_text2d(1050.,175.,card_type_REF
     +                //char(0),PS_ANNOTATION_TEXT)                    !cardiac-modulated acc'ding to Ted's DELTA**2 test?
                    label4='('//DELTA2_REF//', '//tedfactor_REF//')'
                    call remove_all_blanks(label4,LEN(label4))
                    call strlength(label4,LEN(label4),l)
                    call ps_text2d(1040.,160.,label4(1:l)//char(0),
     +                PS_ANNOTATION_TEXT)
                 end if
                 if(card_type_TAR.ne.' ')then
                    call upper_case(card_type_TAR,LEN(card_type_TAR))
                    call ps_text2d(1300.,175.,card_type_TAR//
     +                   char(0),PS_ANNOTATION_TEXT)
                    label4='('//DELTA2_TAR//', '//tedfactor_TAR//')'
                    call remove_all_blanks(label4,LEN(label4))
                    call strlength(label4,LEN(label4),l)
                    call ps_text2d(1290.,160.,label4(1:l)//char(0),
     +                   PS_ANNOTATION_TEXT)
                 end if
              end if

              call ps_character_height(.036)
              call ps_character_width(.012)
              text = 'date of experiment: '//date
              call strlength(text,LEN(text),l)
              call ps_text2d (20.,130.,text(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
              label4 = 'recording: '//recording
              call strlength(label4,LEN(label4),l)
              call ps_text2d (20.,110.,label4(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)
              call strlength(QDT_FILE,LEN(QDT_FILE),l)
              call ps_text2d(20.,70.,
     +                QDT_FILE(1:l)//char(0),PS_ANNOTATION_TEXT)
              text=BDT_FILE
              call strlength(text,LEN(text),l)
              call ps_text2d(20.,90.,text(1:l)
     +                //char(0),PS_ANNOTATION_TEXT)

              call remove_all_blanks(REFcodetext,LEN(REFcodetext))
              call strlength(REFcodetext,LEN(REFcodetext),l_code)
              write (c_ITAL,'(I10)') ITAL(IDs(REFcode))
              call remove_all_blanks(c_ITAL,LEN(c_ITAL))
              call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)
              call strlength(REF_coords(1),LEN(REF_coords(1)),l_1)
              call strlength(REF_coords(2),LEN(REF_coords(2)),l_2)
              call strlength(REF_coords(3),LEN(REF_coords(3)),l_3)

              REFline='R: ID='//REFcodetext(1:l_code) //
     +             ' (# events in spike train = '//
     +             c_ITAL(1:l_ITAL)//'; coords = '//
     +             REF_coords(1)(1:l_1)//' AP, '//
     +             REF_coords(2)(1:l_2)//' RL, '//
     +             REF_coords(3)(1:l_3)//' D)'

              call strlength(REFline,LEN(REFline),l)
              call ps_text2d (20.,45.,REFline(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

              call remove_all_blanks(TARcodetext,LEN(TARcodetext))
              call strlength(TARcodetext,LEN(TARcodetext),l_code)
              write (c_ITAL,'(I10)') ITAL(IDs(TARcode))
              call remove_all_blanks(c_ITAL,LEN(c_ITAL))
              call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)
              call strlength(TAR_coords(1),LEN(TAR_coords(1)),l_1)
              call strlength(TAR_coords(2),LEN(TAR_coords(2)),l_2)
              call strlength(TAR_coords(3),LEN(TAR_coords(3)),l_3)
              TARline='T: ID='//TARcodetext(1:l_code)//
     +          ' (# events in spike train = '//
     +          c_ITAL(1:l_ITAL)//'; coords = '//
     +             TAR_coords(1)(1:l_1)//' AP, '//
     +             TAR_coords(2)(1:l_2)//' RL, '//
     +             TAR_coords(3)(1:l_3)//' D)'
              call strlength(TARline,LEN(TARline),l)
              call ps_text2d(20.,25.,TARline(1:l)//char(0),
     +             PS_ANNOTATION_TEXT)

            if(prim.eq.'Flat')then
               label1 = 'FLAT'
               call strlength(label1,LEN(label1),l)
               call ps_text2d(760.,120.,label1(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)
               if(rescom.ne.' ')then
                  label1 = 'Comment: '//rescom
               call strlength(label1,LEN(label1),l)
                  call ps_text2d(600.,70.,label1(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)
               end if
            else if(prim.ne.' ')then
               call strlength(loc1,LEN(loc1),i)
               call strlength(prim,LEN(prim),j)
               if(prim.ne.'Flat')then
                  if((prim.eq.'Peak').or.(prim.eq.'Trough'))then
                     call strlength(statcomm,LEN(statcomm),k)
                     call remove_all_blanks(cdet,LEN(cdet))
                     call strlength(cdet,LEN(cdet),l)
                     if(cdet.ne.' ')then
c                        label1='Primary feature: '//loc1(1:i)//' '//
c     +                       prim(1:j)//' (d.i. = '//
c     +                       cdet(1:l)//'; '//
c     +                       statcomm(1:k)//')'
                        label1='Primary feature: '//loc1(1:i)//' '//
     +                       prim(1:j)
                        label1a=' (d.i. = '//
     +                       cdet(1:l)//'; '//
     +                       statcomm(1:k)//')'
                     else
                        label1='Primary feature: '//loc1(1:i)//' '//
     +                       prim(1:j)//' (no statistics calculated)'
                        end if
                  else
                     label1='Primary feature: '//loc1(1:i)//' '//
     +                  prim(1:j)
                  end if
                  call strlength(label1,LEN(label1),l)
                  call ps_text2d(500.,130.,label1(1:l)//char(0),
     +                 PS_ANNOTATION_TEXT)
                  call strlength(label1a,LEN(label1a),l)
                  call ps_text2d(520.,110.,label1a(1:l)//char(0),
     +                 PS_ANNOTATION_TEXT)
               end if
               if((sec.ne.'None').and.(sec.ne.' '))then
                  do i = LEN(loc2),1,-1
                     if(loc2(i:i).ne.' ')exit
                  end do
                  do j = LEN(sec),1,-1
                     if(sec(j:j).ne.' ')exit
                  end do
                  label1 = 'Secondary feature: '//loc2(1:i)//' '//
     +                      sec(1:j)
                  call strlength(label1,LEN(label1),l)
                  call ps_text2d(500.,90.,label1(1:l)//char(0),
     +                 PS_ANNOTATION_TEXT)
               else
               end if
               if(rescom.ne.' ')then
                  label1 = 'Comment: '//rescom
                  call strlength(label1,LEN(label1),l)
                  call ps_text2d(500.,70.,label1(1:l)//char(0),
     =                 PS_ANNOTATION_TEXT)
               end if
            end if
           
      call strlength(ps_file,LEN(ps_file),l)
      call ps_fclose(ps_file(1:l)//char(0))

        if((task.eq.'p').or.(task.eq.'v'))then
           if(task.eq.'p')then
              isys=SYSTEM('lpr temp.ps'//char(0))       !hard copy of the ps file
           else
              isys=SYSTEM('gv temp.ps'//char(0))                !pop-up display of the ps file
           end if
           isys=SYSTEM('rm temp.ps'//char(0))
        end if

        return
        end

*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or write the CTH window

*     date of last revision: 05-Feb-2004        lss


        subroutine print_or_write_CTHs(task,c_format,cellCTH,
     +                  phrenicCTH,
     +                  norm_cellCTH,norm_phrenicCTH,cardCCH,
     +                  card_overlay,bw_p,bw_n,bw_c,
     +                  bwtext1,bwtext2,bwtext3,
     +                  bwtext4,bwtext5,
     +                  REFcode,REFevents_p,REFevents_n,REFevents_c,
     +                  TARcode,TARevents,
     +                  BDT_FILE,TAR_ACH,TAR_ACH_1,TARname,
     +                  QDT_FILENAME,sort,
     +                  mode,resp_type,AP,RL,dep,
     +                  AA,AA_results,AA_text1,AA_applied,
     +                  STA,STA_results,STA_text1,STA_applied,
     +                  STA_res_text,per_results,per_text_abbrev,
     +                  zmodsig,zmodsig2,ETA2,coef,card_type,card,
     +                  DELTA2,tedfactor,fiveHT,mean_rISI,sd_rISI,
     +                  c_MAX_INT,num_rej_ISI,bw_ach1,bw_ach2,
     +                  REFevents_ach1,REFevents_ach2,IDs,ITAL,
     +                  cardiac_pls,STIM_OFFSET,NORM_OFFSET,mean_E,comm,
     +                  show_pulses)
*
      use mod_miscellaneous_subroutines
      use mod_ps_plot_1
        integer*4 isys,SYSTEM
        integer cellCTH(101),phrenicCTH(101),norm_cellCTH(101),
     +            norm_phrenicCTH(101),cardCCH(101),card_overlay(101),
     +            TAR_ACH(101),TAR_ACH_1(101),scaled_hist(101),
     +            NTOP1,NTOP2,NTOP3,NTOP4,NTOP5,cardiac_pls
*
        include 'x2000parameter.defs'
*
        integer REFcode,TARcode,REFevents_p,REFevents_n,
     +          REFevents_c,TARevents,REFevents_ach1,REFevents_ach2,
     +          IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
        real STIM_OFFSET,NORM_OFFSET,mean_E
        character*10 fiveHT,mean_rISI,sd_rISI,c_MAX_INT
*
        character*(*) BDT_FILE,c_format,QDT_FILENAME,comm
        character*12 bwtext1,bwtext2,bwtext3,bwtext4,bwtext5,
     +                  AA_text1(MAX_AA),STA_text1(MAX_STA),
     +                  STA_results(MAX_STA),
     +                  per_text_abbrev(MAX_PERTURB)
        character*10 num_rej_ISI,c_ITAL
        character*10 c_rate_p,c_rate_n,c_rate_c,c_rate_ach1,c_rate_ach2
        character*60 title,message
        character*100 ps_file
        character*50 text,text2
        character*80 text80
        character*100 REFline,TARline
        character*15 resp_type,STA_res_text(122)
        character*9 per_results(MAX_PERTURB)
        character*7 REFeventstext,TAReventstext
        character*5 AP,RL,dep,ETA2,coef,DELTA2
        character*4 TARname
        character*4 TARcodetext,REFcodetext
        character*3 CTH,ACH,CCC,AA_results(MAX_AA),
     +              zmodsig,zmodsig2,card_type
        character*2 mode,tedfactor
        character*1 show_control,AA(MAX_AA),STA(MAX_STA),
     +                  AA_applied,STA_applied,
     +                  sort,card,task,scaled_up,info,show_pulses

*       ***** suppress unused variable warnings *****
        if(.false.)print *,c_format
        if(.false.)print *,STA_applied

*
*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
        call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
        call ps_mapping_mode(1)
        call ps_view_port(0.,0.,1.,1.)
        call ps_view_window(1.,1.,1700.,500.)
        call ps_geometry(18.,18.,792.-18.,350.-18.)
*       *********************************************************************
*
        write (REFcodetext,'(I4)') REFcode
        write (TARcodetext,'(I4)') TARcode
        if(task.eq.'w')then
           call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
           ps_file=QDT_FILENAME(1:l)//'_CTH_'//TARcodetext//'.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if
c        print '(''ps_file = '',A)',ps_file

        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))
c        print '(''ps_file opened'')'

*     ***** draw to the ps_file: *****

        min_bin=0
        max_bin=0
        scaled_up='n'
        info = 'n'
        rmagnify = 0.0
        show_control = 'y'                              !default value - always show the overlay for a CTH

        write (TARcodetext,'(I4)') TARcode
        write (TAReventstext,'(I7)') TARevents
        write (REFeventstext,'(I7)') REFevents_p
        title=' '
        title = 'CTHs for cell '//TARname
        message = '(ID code = '//TARcodetext//')'
        if(TARname.eq.' ')title = 
     +       'CTHs for cell with ID code = '//TARcodetext
        TARline=' '
        CTH = 'CTH'
        ACH = 'ACH'
        CCC = 'CCC'
*
*
*       *****   label the CTHs according to binwidths, input filename, & codes: *****
*
        call ps_character_width(0.020)
        call ps_character_height(0.1)
        call strlength(title,LEN(title),l)
        call ps_text2d (10.,450.,title(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
        if(TARname.ne.' ')then
           call ps_character_height(.060)
           call ps_character_width(.015)
           call strlength(message,LEN(message),l)
           call ps_text2d(10.,425.,message(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
        end if
*
        call ps_character_width(0.010)
        call ps_character_height(0.050)
        call ps_line_type(DOT)
        text = '<-- binwidths -->'
        call strlength(text,LEN(text),l)
        call ps_text2d (300.,80.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
         call strlength(bwtext1,LEN(bwtext1),l)
         call ps_text2d (600.,80.,bwtext1(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
         call strlength(bwtext2,LEN(bwtext2),l)
         call ps_text2d (1020.,80.,bwtext2(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
         if(cardiac_pls.ne.0)then
            call strlength(bwtext3,LEN(bwtext3),l)
            call ps_text2d (1450.,80.,bwtext3(1:l)//char(0),
     +           PS_ANNOTATION_TEXT)
         end if
         call remove_all_blanks(bwtext4,LEN(bwtext4))
         call strlength(bwtext4,LEN(bwtext4),l)
         call ps_text2d (75.,82.,bwtext4(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
         call strlength(bwtext5,LEN(bwtext5),l)
         call ps_text2d (180.,82.,bwtext5(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)

         call remove_all_blanks(REFcodetext,LEN(REFcodetext))
         call strlength(REFcodetext,LEN(REFcodetext),l_code)

         write (c_ITAL,'(I10)') ITAL(IDs(REFcode))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)
        REFline='R: ID='//REFcodetext(1:l_code) //
     +          ' (# E pulses in spike train = '//
     +          c_ITAL(1:l_ITAL)//')'
        call strlength(REFline,LEN(REFline),l)
        call ps_text2d (1.,45.,REFline(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

         call remove_all_blanks(TARcodetext,LEN(TARcodetext))
         call strlength(TARcodetext,LEN(TARcodetext),l_code)
         write (c_ITAL,'(I10)') ITAL(IDs(TARcode))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)
         TARline='T: ID='//TARcodetext(1:l_code) //
     +          ' (# events in spike train = '//
     +          c_ITAL(1:l_ITAL)//')'
        call strlength(TARline,LEN(TARline),l)
        call ps_text2d(1.,25.,TARline(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

         call remove_all_blanks(REFeventstext,LEN(REFeventstext))
         call strlength(REFeventstext,LEN(REFeventstext),l_Revents)
         call remove_all_blanks(TAReventstext,LEN(TAReventstext))
         call strlength(TAReventstext,LEN(TAReventstext),l_Tevents)
        text = 'sum of all bins = '//TAReventstext(1:l_Tevents)//
     +          '; # cycles shown = '//REFeventstext(1:l_Revents)
        call strlength(text,LEN(text),l)
        call ps_text2d(1.,5.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

        call strlength(BDT_FILE,LEN(BDT_FILE),l_bdt)
        call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
        text = BDT_FILE(1:l_bdt)//' / '//QDT_FILENAME(1:l)
        call strlength(text,LEN(text),l)
        call ps_text2d(1.,65.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

        if(sort.eq.'c')then
           text='ACHs: (all-order; clean)'
        else if(sort.eq.'m')then
           text='ACHs: (all-order; messy)'
        else
           text='ACHs:'
        end if   
        call strlength(text,LEN(text),l)
        call ps_text2d (10.,220.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

        text = 'RESPIRATORY CTH'
        call strlength(text,LEN(text),l)
        call ps_text2d (490.,60.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
        text = '(phrenic overlay)'
        call strlength(text,LEN(text),l)
        call ps_text2d (520.,40.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
        if (norm_cellCTH(101).gt.0) then
           text = 'PHASE NORMALIZED RESP. CTH'
        else
           text = 'NORMALIZED RESP. CTH'
        end if
        call strlength(text,LEN(text),l)
        call ps_text2d (875.,60.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
        if (norm_cellCTH(101).gt.0) then
           text = '(phase norm. phrenic overlay)'
        else
           text = '(norm. phrenic overlay)'
        end if
        call strlength(text,LEN(text),l)
        call ps_text2d (900.,40.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
        if(cardiac_pls.ne.0)then
           text = 'CARDIAC CCH'
           call strlength(text,LEN(text),l)
           call ps_text2d (1400.,60.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
           text = '(cardiac pulse overlay)'
           call strlength(text,LEN(text),l)
           call ps_text2d (1350.,40.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
        else
           text = 'CARDIAC CCHs NOT CALCULATED'
           call strlength(text,LEN(text),l)
           call ps_text2d(1300.,250.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
        end if
        
*       ***** find the maximum bin in each CTH and calculate the peak firing rate: *****
        
        NTOP1=0
        NTOP2=0
        NTOP3=0
        NTOP4=TAR_ACH(1)
        NTOP5=TAR_ACH_1(1)
        do I=1,100
          NTOP1=MAX0(NTOP1,cellCTH(I))
          NTOP2=MAX0(NTOP2,norm_cellCTH(I))
          NTOP3=MAX0(NTOP3,cardCCH(I))
          NTOP4=MAX0(NTOP4,TAR_ACH(I))
          NTOP5=MAX0(NTOP5,TAR_ACH_1(I))
       end do
       rate = (NTOP1*(1000.0/bw_p))/REFevents_p
          write (c_rate_p,'(F9.1)')rate
          call remove_all_blanks(c_rate_p,LEN(c_rate_p))
          call strlength(c_rate_p,LEN(c_rate_p),l)
          text = c_rate_p(1:l)//' spikes/sec'
          call strlength(text,LEN(text),l)
          call ps_text2d (420.,405.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
       rate = (NTOP2*(1000.0/bw_n))/REFevents_n
          write (c_rate_n,'(F9.1)')rate
          call remove_all_blanks(c_rate_n,LEN(c_rate_n))
          call strlength(c_rate_n,LEN(c_rate_n),l)
          text = c_rate_n(1:l)//' spikes/sec'
          call strlength(text,LEN(text),l)
          call ps_text2d (840.,405.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
       if(cardiac_pls.ne.0)then
          rate = (NTOP3*(1000.0/bw_c))/REFevents_c
          write (c_rate_c,'(F9.1)')rate
          call remove_all_blanks(c_rate_c,LEN(c_rate_c))
          call strlength(c_rate_c,LEN(c_rate_c),l)
          text = c_rate_c(1:l)//' spikes/sec'
          call strlength(text,LEN(text),l)
          call ps_text2d (1270.,405.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
       end if
       rate = (NTOP4*(1000.0/bw_ach1))/REFevents_ach1
          write (c_rate_ach1,'(F9.1)')rate
          text = c_rate_ach1
          call remove_all_blanks(text,LEN(text))
          call strlength(text,LEN(text),l)
          call ps_text2d (40.,203.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)
       rate = (NTOP5*(1000.0/bw_ach2))/REFevents_ach2
          write (c_rate_ach2,'(F9.1)')rate
          call remove_all_blanks(c_rate_ach2,LEN(c_rate_ach2))
          call strlength(c_rate_ach2,LEN(c_rate_ach2),l)
          text = c_rate_ach2(1:l)//' spikes/sec'
          call strlength(text,LEN(text),l)
          call ps_text2d (160.,203.,text(1:l)//char(0),
     +         PS_ANNOTATION_TEXT)

*       *****   scale the CTHs: *****

      call ps_plot(100,cellCTH,phrenicCTH,scaled_hist,CTH,420.,100.,   !plot resp CTH and overlay
     +         300.,404.,info,'','','','',0,0,show_control,0,'',
     +         [integer::], max_bin,'','','',[integer::],'','',
     +         scaled_up, min_bin, rmagnify, [integer::],[integer::])
        if(show_pulses.ne.'n')then
           call ps_move2d(420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +          100.-5.)
           call ps_draw2d(420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +          100.+5.)
           call ps_text2d(420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +          80.,'E'//char(0),PS_ANNOTATION_TEXT)
        end if

      call ps_plot(100,norm_cellCTH,norm_phrenicCTH,scaled_hist,CTH,   !plot norm CTH and overlay
     +       840.,100.,300.,404.,info,'','','','',0,0,show_control,0,'',
     +       [integer::], max_bin,'','','',[integer::],'','',
     +       scaled_up, min_bin, rmagnify, [integer::],[integer::])
        if(show_pulses.ne.'n')then
           call ps_move2d(840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          100.-5.)
           call ps_draw2d(840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          100.+5.)
           call ps_text2d(840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +          80.,'E'//char(0),PS_ANNOTATION_TEXT)
           if(mean_E.ne.0 0)then
              call ps_move2d(840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +             (100.*bw_n))),100.-5.)
              call ps_draw2d(840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +             (100.*bw_n))),100.+5.)
              call ps_text2d(840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +             (100.*bw_n))),80.,'I'//char(0),PS_ANNOTATION_TEXT)
           end if
        end if

        if(cardiac_pls.ne.0)then
           call ps_plot(100,cardCCH,card_overlay,scaled_hist,CCC,1270.,!plot cardiac CCH and overlay
     +          100.,300.,404.,info,'','','','',0,0,show_control,0,'',
     +          [integer::], max_bin,'','','',[integer::],'','',
     +          scaled_up, min_bin, rmagnify, [integer::],[integer::])                             
        end if

      call ps_plot(100,TAR_ACH,TAR_ACH,scaled_hist,ACH,40.,100.,100.,  !plot ACHs (smallest and next-to-greatest BW)
     +       101.,info,'','','','',0,0,show_control,0,'',
     +       [integer::], max_bin,'','','',[integer::],'','',
     +       scaled_up, min_bin, rmagnify, [integer::],[integer::])
      call ps_plot(100,TAR_ACH_1,TAR_ACH_1,scaled_hist,ACH,160.,100.,   
     +     100.,101.,info,'','','','',0,0,show_control,0,'',
     +     [integer::], max_bin,'','','',[integer::],'','',
     +     scaled_up, min_bin, rmagnify, [integer::],[integer::])
*
           call ps_character_height(.060)
           call ps_character_width(.015)
           call strlength(resp_type,LEN(resp_type),l)
           call ps_text2d(10.,400.,resp_type(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)
           call ps_character_height(.050)
           call ps_character_width(.010)

           if(mode.ne.'vt')then                                 !do not display the following info if only looking:
              text = 'AP, RL, D: '//AP//', '//
     +                RL//', '//dep//' mm.'
              call strlength(text,LEN(text),l)
              call ps_text2d(10.,380.,text(1:l)//char(0),
     +                    PS_ANNOTATION_TEXT)
              y = 360.
              call ps_text2d(10.,y,'AA:'//char(0),
     +             PS_ANNOTATION_TEXT)
              do ii = 1,MAX_AA
                 if(AA(ii).eq.'y')then
                    text=AA_results(ii)//'  '//
     +                AA_text1(ii)
                    call strlength(text,LEN(text),l)
                    call ps_text2d(50.,y,text(1:l)//char(0),
     +                PS_ANNOTATION_TEXT)
                    y = y-15.
                 end if
              end do

              if(AA_applied.eq.'y')then
                 y = y - 5.
              else
                 y = y-20.
              end if
              call ps_text2d(10.,y,'STA:'//char(0),
     +             PS_ANNOTATION_TEXT)
              do ii = 1,MAX_STA
                 if(STA(ii).eq.'y')then
                    do jj=LEN(STA_results(ii)),1,-1
                       if(STA_results(ii)(jj:jj).ne.' ')exit
                    end do
                    text=STA_results(ii)(1:jj)//' '//
     +               STA_res_text(ichar(STA_results(ii)(1:1)))
                    call strlength(text,LEN(text),l)
                    call ps_text2d(60.,y,text(1:l)//char(0),
     +                          PS_ANNOTATION_TEXT)
                    call strlength(STA_text1(ii),LEN(STA_text1(ii)),l)
                    call ps_text2d(250.,y,STA_text1(ii)(1:l)//char(0),
     +                          PS_ANNOTATION_TEXT)
                    y = y-15.
                 end if
              end do
 
              y = 480.
              call ps_text2d(300.,y,'Responses:'//char(0),
     +                  PS_ANNOTATION_TEXT)
              m = 0
              x = 400.
              do ii = 1, MAX_PERTURB
                 if(per_results(ii).ne.' ')then
                    m=m+1
                    if(m.eq.1)x=450.
                    if(m.ne.1)x=x+250.
                    do jj=LEN(per_text_abbrev(ii)),1,-1
                       if(per_text_abbrev(ii)(jj:jj).ne.' ')exit
                    end do
                    do kk=LEN(per_results(ii)),1,-1
                       if(per_results(ii)(kk:kk).ne.' ')exit
                    end do
                    text=per_text_abbrev(ii)(1:jj)//': '
     +                   //per_results(ii)
                    call strlength(text,LEN(text),l)
                    call ps_text2d(x,y,text(1:l)//char(0),
     +                          PS_ANNOTATION_TEXT)
                 end if
                 if(m.eq.4)then
                    m = 0       !reset
                    y = y-15.
                 end if
              end do

              call ps_text2d(300.,425.,'Comments:'//char(0),
     +                  PS_ANNOTATION_TEXT)
              call strlength(comm,LEN(comm),l)
              call ps_text2d(450.,425.,comm(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)

              if(cardiac_pls.ne.0)then
                 if(card.eq.'c')then
                    text='cardiac mod.:  visual --> YES'
                 else if(card.eq.'n')then
                    text='cardiac mod.:  visual --> NO'
                 else
                    text = 'cardiac mod.:  visual -->'
                 end if
                 call strlength(text,LEN(text),l)
                 call ps_text2d(1250.,5.,text(1:l)//char(0),
     +                PS_ANNOTATION_TEXT)
              end if
           end if
           
           if((zmodsig.eq.'r').or.(zmodsig.eq.'R'))then
              text='resp. mod.:  ANOVA  --> YES'
           else if((zmodsig.eq.'n').or.(zmodsig.eq.'N'))then
              text='resp. mod.:  ANOVA  --> NO'
           else if((zmodsig.eq.'ned').or.(zmodsig.eq.'NED'))then
              text='resp. mod.:  ANOVA  --> NO (but insuff. data)'
           else
              text='resp. mod.:  ANOVA  -->'
           end if
           call strlength(text,LEN(text),l)
           call ps_text2d(860.,20.,text(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)

           if((zmodsig2.eq.'r').or.(zmodsig2.eq.'R'))then
              text='resp. mod.:  BINARY --> YES'
           else if((zmodsig2.eq.'n').or.(zmodsig2.eq.'N'))then
              text='resp. mod.:  BINARY --> NO'
           else if((zmodsig2.eq.'ned').or.(zmodsig2.eq.'NED'))then
              text='resp. mod.:  BINARY --> NO (but insuff. data)'
           else
              text='resp. mod.:  BINARY --> '
           end if
           call strlength(text,LEN(text),l)
            call ps_text2d(860.,5.,text(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)

           text='ETA = '//ETA2//
     +          ';  c.o.v. = '//coef
            call ps_text2d(500.,20.,text//char(0),
     +              PS_ANNOTATION_TEXT)
            if(cardiac_pls.ne.0)then
               text='(using cycles in CTH)'
               call strlength(text,LEN(text),l)
               call ps_text2d(505.,5.,text(1:l)//char(0),
     +              PS_ANNOTATION_TEXT)
*
               if(DELTA2.ne.' ')then                                   !do not display if no value for DELTA2
                  if((card_type.eq.'c').or.(card_type.eq.'C'))then
                     text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +                    ' --> YES'
                  else if((card_type.eq.'n').or.(card_type.eq.'N'))then
                     text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +                    ' --> NO'
                  else
                     text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +                    ' --> insuff. data'
                  end if
                  call strlength(text,LEN(text),l)
                  call ps_text2d(1250.,20.,text(1:l)//char(0),
     +                 PS_ANNOTATION_TEXT)
               end if
            end if
            if(fiveHT.ne.' ')then                               !do not display this stuff if no value for fiveHT
               text2 = '('//fiveHT//')'                         !display 5HT and relativeISI vlaues
               call remove_all_blanks(text2,LEN(text2))
               do i = 1,LEN(fiveHT)
                  if(fiveHT(i:i).ne.' ')exit                    !find the first "real" character in the string
               end do
               if(fiveHT(i:i).eq.'-')then
                  text80='5HT? --> YES '//text2
               else
                  text80='5HT? --> NO '//text2
               end if
               call upper_case(text80,LEN(text80))
               call strlength(text80,LEN(text80),l)
               call ps_text2d(1390.,480.,text80(1:l),
     +                     PS_ANNOTATION_TEXT)
               text2=mean_rISI//' +/- '//sd_rISI
               call remove_all_blanks(text2,LEN(text2))
               text80='mean rel. ISI = '//text2
               call strlength(text80,LEN(text80),l)
               call ps_text2d(1390.,460.,text80(1:l),
     +                     PS_ANNOTATION_TEXT)
               text2 = '(>'//c_MAX_INT//'ms.)'
               call remove_all_blanks(text2,LEN(text2))
               text80=num_rej_ISI//' ISI rejected '//text2
               call remove_leading_blanks(text80,LEN(text80))
               call strlength(text80,LEN(text80),l)
               call ps_text2d(1390.,440.,text80(1:l),
     +                     PS_ANNOTATION_TEXT)
            end if

        call strlength(ps_file,LEN(ps_file),l)
        call ps_fclose(ps_file(1:l)//char(0))

        if((task.eq.'p').or.(task.eq.'v'))then
           if(task.eq.'p')then
              isys=SYSTEM('lpr temp.ps'//char(0))              !hard copy of the ps file
           else
              isys=SYSTEM('gv temp.ps'//char(0))              !pop-up display of the ps file
           end if
           isys=SYSTEM('rm temp.ps'//char(0))
        end if

*
*
        call ps_view_window(1.,1.,1500.,500.)
*
*
        RETURN
        END





*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or write the QSUM window

*     date of last revision: 11-Feb-2004        lss


*       11-Feb-2004 code added (courtesy of ro'c) that will properly scale the view window for the postscript files.
*               We were having problems with line width for the CUSUM plots, and didn't want to have to force the
*               line width to zero to avoid later problems with CorelDraw.

        subroutine print_or_write_QSUM(task,IHIST,BINVAL,ICN,
     +                    COEF,BDT_FILE,QDT_FILENAME,date,
     +                    recording,NHW,rec_num,title,
     +                    REFcode,
     +                    TARcode,first_bin,last_bin,IDs,ITAL)
*
*
      use mod_miscellaneous_subroutines
        include 'x2000parameter.defs'
        PARAMETER(iSOLID=0,iDOT=1)
        
        integer*4 isys,SYSTEM
        integer IHIST(101),KDAT(101),ICN(6),IDs(MAX_NUM_CODES),
     +       ITAL(MAX_NUM_CHAN)
*
        real COEF,BINVAL,CDAT(101)
*
        integer*4 first_bin,last_bin,total_bins
        integer rec_num,REFcode,TARcode
*
        character*160 text160
        character*120 title,text
        character*100 ps_file
        character*(*) BDT_FILE,QDT_FILENAME
        character*11 date
        character*10 c_ITAL
        character*7 c_ICN(6),c_binwd
        character*6 c_ATOP,c_ABOT
        character*5 c_rec_num
        character*4 c_NHW
        character*3 c_first_bin,c_last_bin
        character*2 recording
        character*1 task

*       ***** suppress unused variable warnings *****
        if(.false.)print *,REFcode
        if(.false.)print *,rec_num
        if(.false.)print *,TARcode
*
*        write (c_rec_num,'(I5)')rec_num
         write (c_NHW,'(I4)')NHW
         binwd = (float(NHW)/50.)
         write (c_binwd,'(F7.1)')binwd
         do i = 1,6                             !convert to character format
            write (c_ICN(i),'(6I7)')ICN(i)
            call remove_all_blanks(c_ICN(i),LEN(c_ICN(i)))
         end do

        if(task.eq.'w')then
           call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
           ps_file=QDT_FILENAME//'_'//c_ICN(1)//'x'//c_ICN(3)//
     +             '_QSUM_'//c_binwd//'.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if

c        print '(''ps_file: '',A)', ps_file
        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))


*     ***** draw to the ps_file: *****

        SUMX=0.0
        SUMC=0.0
*
        KDAT=IHIST              !load CCH data into KDAT for later manipulation
*
*
        total_bins=(last_bin-first_bin)+1       !determine total number of bins included
                                                ! within the range
        do j = first_bin,last_bin               !calc total of bin values within the range
         SUMX=SUMX + KDAT(j)
         end do
        average = SUMX/total_bins               !average value of range bin
        DO I=1,101
         SUMC=SUMC+(float(KDAT(I))-average)     !subtract the average from each bin of the CCH
         CDAT(I)=SUMC                           ! and keep a running total of the adjusted bin
*                                               ! values; store the running totals in CDAT()
        end do
*
*
*       scale everything used from here on to # of impulses/trigger event...
*
        do  I=1,101
         CDAT(I) = CDAT(I)/float(ICN(2))                !ICN(2) = total number of reference spikes
         end do
        average = average/float(ICN(2))
*
*
*       Q_SUM CALC COMPLETE..DEFINE Y AXIS AND PLOT
*
        TOP=CDAT(1)
        BOT=CDAT(1)
*
        DO I=2,101
         TOP=AMAX1(TOP,CDAT(I))         !TOP = maximum value in CDAT
         BOT=AMIN1(BOT,CDAT(I))         !BOT = minimum value in CDAT
         end do
*
        BIG=AMAX1(ABS(TOP),ABS(BOT))    !BIG = the larger of the 2 absolute values for
                                        !   TOP and BOT
        BORDER=0.1*(DIM(BIG,-BIG))      !define a border according to the total spread
                                        !   to be used for the Q_SUM plot 
        ATOP=BIG+BORDER                 !maximum and minimum excursions for the Q_SUM 
        ABOT=-BIG-BORDER                !   plot itself (border included)
c        print '(''ATOP: '',f7.2,''; ABOT: '',f7.2)',ATOP,ABOT
*
*

*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
*       *****    (machinations for the view_window call avoid line width scaling problems)
                 call ps_vdc_extent(0.0,0.0,0.0,1.25,1.0,0.0)
                 call ps_mapping_mode(1)
                 call ps_view_port(.15*1.25,.15,.85*1.25,.75)
c                call ps_view_window(85.,ABOT,590.,ATOP)        !old view_window call
                 call ps_geometry(18.,18.,792.-18.,612.-18.)
                 f = (792.-18.-18.)/(590.-82.)                  !f=window width in pixels/(x2-x1)
                 g = (612.-18.-18.)/(ATOP-ABOT)                 !g=window height in pixels/(y2-y1)
                 call ps_view_window(f*85.,g*ABOT,f*590.,g*ATOP)
*       *********************************************************************

         call ps_character_height(.035)
         call ps_character_width(.012)
         write (c_ATOP,'(f5.3)') ATOP
         write (c_ABOT,'(f5.3)') ABOT
         call strlength(c_ATOP,LEN(c_ATOP),l)
         call ps_text2d(f*45.,g*ATOP-5.,c_ATOP(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
         call strlength(c_ABOT,LEN(c_ABOT),l)
         call ps_text2d(f*45.,g*ABOT+5.,c_ABOT(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)


C       *****  draw axes from 85,ATOP to 85,ABOT to 590,ABOT***
          call ps_line_type(iSOLID)
        call ps_move2d(f*85.,g*ATOP)
        call ps_draw2d(f*85.,g*ABOT)
        call ps_draw2d(f*590.,g*ABOT)
c       *****  draw a broken line from 338,200 to 338,700 (middle bin) ****
        call ps_line_type(iDOT)
        call ps_move2d(f*338.,g*ABOT)
        call ps_draw2d(f*338.,g*ATOP)
        call ps_line_type(iSOLID)

       call ps_move2d(f*85.,g*ABOT)
        IZ=101
        IX=85
       DO I=1,IZ
        J=I+1
        IF(J.EQ.102)cycle
        TEMP=CDAT(J)
        IX2=IX+5
        call ps_move2d(f*float(IX),g*CDAT(I))
        call ps_draw2d(f*float(IX2),g*CDAT(I))
        call ps_draw2d(f*float(IX2),g*TEMP)
        IX=IX2
        end do
        
c
c
C       PLOT CONFIDENCE LIMITS ON Q_SUM
C
C
        IF (COEF.GT.0.9)
     +      print '(''Coef. of var. .GT. 0.9.'')' 
        ZIX=85.+((last_bin-1)*5.)       !physical location (within the display window) of the last
                                        ! bin of the range
        YC1=0.0
        YCM1=0.0
        LPCT=101-last_bin       !# OF BINS TO HAVE CONFIDENCE BAND
        U=(total_bins*BINVAL)   !TOTAL TIME OF CONTROL PERIOD (BINS SPANNED IN SEC.)
        ZMEAN= BINVAL/average
        C4=(1.0/6.0)-(COEF**4/6.0)
        C2=COEF**2
        UM=U*ZMEAN
        ISKIP=0  !WHEN =1 SKIP V2 CALC BELOW
        DO 1500 J=1,LPCT
          V1=ICN(2)*(((C2*J*BINVAL)/ZMEAN)
     +       +C4+(((COEF*J*BINVAL)**2)/UM))
          V2=(ICN(2)*(J*BINVAL))/ZMEAN
          IF(V1.LT.V2)ISKIP=1
          VX=V2
          IF(ISKIP.EQ.1)VX=V1
          if (COEF.GT.0.9) VX=V2
          CURSIG=((3*SQRT(VX))/ICN(2)) 
C       
C       DRAW CONFIDENCE BAND ELEMENTS ABOVE AND BELOW THE JTH BIN IN LOOP
C
c
          ZIX2=ZIX+5.
          YC2=CURSIG
          YCM2=-CURSIG
          if (CURSIG .gt. ATOP) goto 1000
          call ps_move2d(f*ZIX,g*YC1)
          call ps_draw2d(f*ZIX2,g*YC2)          !draw line above 0
          call ps_move2d(f*ZIX,g*YCM1)
          call ps_draw2d(f*ZIX2,g*YCM2)         !draw line below 0
          ZIX=ZIX2
          YC1=YC2
          YCM1=YCM2
1500     CONTINUE

1000     call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)    !get ready to print the accompanying text
         call ps_view_port(0.,0.,1.,1.)
         call ps_view_window(85.,10.,590.,500.)         !xmin,ymin,xmax,ymax
         call ps_geometry(18.,18.,792.-18.,612.-18.)
         call ps_character_height(.035)
         call ps_character_width(.012)
*
*
         text='change in'
         call strlength(text,LEN(text),l)
         call ps_text2d(95.,300.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
         text='# spikes'
         call strlength(text,LEN(text),l)
         call ps_text2d(95.,285.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
         text='per'
         call strlength(text,LEN(text),l)
         call ps_text2d(95.,270.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)
         text='REF event'
         call strlength(text,LEN(text),l)
         call ps_text2d(95.,255.,text(1:l)//char(0),
     +          PS_ANNOTATION_TEXT)

         text='.bdt filename:  '//BDT_FILE
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,480.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
         text='.qdt filename:  '//
     +        QDT_FILENAME(1:l)
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,465.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='date of experiment:  '//date
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,450.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='recording #:  '//recording
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,435.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='record # in .qdt file:  '//c_rec_num
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,65.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call strlength(c_ICN(1),LEN(c_ICN(1)),l_1)
         write (c_ITAL,'(I10)') ITAL(IDs(ICN(1)))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

         text='R: ID = '//c_ICN(1)(1:l_1)//
     +          ' (# events in spike train = '//
     +          c_ITAL(1:l_ITAL)//')'
         call strlength(text,LEN(text),l)
c        text='ref:  ID = '//c_ICN(1)//'  ('//c_ICN(2)//' events)'
c          call strlength(text,LEN(text),l)
         call ps_text2d(85.,50.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         call strlength(c_ICN(3),LEN(c_ICN(3)),l_3)
         call strlength(c_ICN(4),LEN(c_ICN(4)),l_4)
         write (c_ITAL,'(I10)') ITAL(IDs(ICN(3)))
         call remove_all_blanks(c_ITAL,LEN(c_ITAL))
         call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

         text='T: ID = '//c_ICN(3)(1:l_3)//
     +        ' (# events in spike train = '//
     +        c_ITAL(1:l_ITAL)//')'
         call strlength(text,LEN(text),l)
c        text='tar:  ID = '//c_ICN(3)//'  ('//c_ICN(4)//' events)'
c          call strlength(text,LEN(text),l)
         call ps_text2d(85.,35.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='binwidth:  '//c_binwd//' msec.'
          call strlength(text,LEN(text),l)
         call ps_text2d(85.,20.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text='-'//c_NHW
          call strlength(text,LEN(text),l)
         call ps_text2d(120.,75.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)
*
         text=c_NHW//' msec'
          call strlength(text,LEN(text),l)
         call ps_text2d(500.,75.,text(1:l)//char(0),
     +        PS_ANNOTATION_TEXT)

         write(c_first_bin,'(I3)') first_bin
         write(c_last_bin,'(I3)') last_bin
         text = '('//c_first_bin//' , '//c_last_bin//')'
         call remove_all_blanks(text,LEN(text))
         text160 = 'bin range used to calculate average: '//text
         call strlength(text160,LEN(text160),l)
         call ps_text2d(320.,35.,text160(1:l),
     +        PS_ANNOTATION_TEXT)
*
         call ps_character_height(.060)
         call ps_character_width(.020)
         text160 = 'QSUM: '//title
         call strlength(text160,LEN(text160),l)
         call ps_text2d (250.,420.,
     +        text160(1:l),PS_ANNOTATION_TEXT)            !'QSUM: REF>TAR'
*
        call strlength(ps_file,LEN(ps_file),l)
        call ps_fclose(ps_file(1:l)//char(0))

        if((task.eq.'p').or.(task.eq.'v'))then
           if(task.eq.'p')then
              isys=SYSTEM('lpr temp.ps'//char(0))       !send the ps file to the default printer
           else
              isys=SYSTEM('gv temp.ps'//char(0))                !pop-up display of the ps file
           end if
           isys=SYSTEM('rm temp.ps'//char(0))
        end if


        return
        end






*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

*     code to print or write the CCH statistics window

*     date of last revision: 14-Apr-2004        lss

        subroutine print_or_write_STATS(task,CCH,CONTROL,
     +                       BDT_FILE,                               !print the STATISTICS CCH
     +                       QDT_FILE,REF,TAR,date,recording,NHW,
     +                       rec_num,rec_num_control,title,ICN,
     +                       show_control,show_conf_lim,
     +                       REFtype,TARtype,
     +                       show_single_shift,show_avg_shift,
     +                       CONTROL_cl,show_2_sd,show_3_sd,
     +                       czk,cdet,cvis,cprobk,czlat,
     +                       chalfwd,
     +                       IPKS,IPKF,IBKS,IBKF,per_text_abbrev,
     +                       per_results_REF,per_results_TAR,
     +                       min_bin,max_bin,scaledup,rmagnify,
     +                       IDs,ITAL)
      use mod_miscellaneous_subroutines
      use mod_ps_plot_1




        INCLUDE 'x2000parameter.defs'

        integer*4 isys,SYSTEM
        integer CCH(101),CONTROL(101),scaled_hist(101),
     +       rec_num,rec_num_control,ICN(6),REF,
     +       TAR,CONTROL_cl(101),total_num_bins,
     +       IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)

        real x0,y0,                                      !lower-left coordinates (0,0) of histogram axes
     +       height, width
         
        character*120 title,text
        character*100 ps_file
        character*(*) BDT_FILE,QDT_FILE
        character*15 REFtype,TARtype
        character*12 per_text_abbrev(MAX_PERTURB)
        character*11 date
        character*10 bwtext
        character*9 per_results_REF(MAX_PERTURB),
     +              per_results_TAR(MAX_PERTURB)
        character*8 czk,cprobk,cdet,cvis,czlat,chalfwd
        character*2 recording
        character*3 CCHtype,REFtext,TARtext,c_IPKF,c_IBKF

        character*1 info,show_control,task,show_conf_lim,
     +              show_single_shift,show_avg_shift,
     +              show_2_sd,show_3_sd,scaledup

*       *****  set histogram printing and saving parameters per Russ's ps_starbase library: *********
                 call ps_vdc_extent(0.0,0.0,0.0,1.0,1.0,0.0)
                 call ps_mapping_mode(1)
                 call ps_view_port(0.,0.,1.,1.)
                 call ps_view_window(1.,1.,1500.,800.)
                 call ps_geometry(18.,18.,792.-18.,612.-18.)
*       *********************************************************************

        write (bwtext,'(f7.1)') (float(NHW)/50.)
        write (REFtext,'(I3)') REF
        write (TARtext,'(I3)') TAR
        CCHtype = 'CCH'

        if(task.eq.'w')then
           call strlength(QDT_FILE,LEN(QDT_FILE),l)
           ps_file=QDT_FILE(1:l)//'_'//REFtext//'x'//TARtext//'_'//
     +          'STATS_'//bwtext//'.ps'
           do i = 1,LEN(ps_file)-3
              if(ps_file(i:i+3).eq.'.qdt')then
                 ps_file(i:i+3)='    '
                 exit
              end if
           end do
           call remove_all_blanks(ps_file,LEN(ps_file))
        else
           ps_file = 'temp.ps'
        end if

        call strlength(ps_file,LEN(ps_file),l)
        call ps_fopen(ps_file(1:l)//char(0))
        info = 'y'
        x0 = 500.
        y0 = 300.
        height = 300.
        width = 525.
        total_num_bins = 101

        call ps_plot(total_num_bins,CCH,CONTROL,scaled_hist,CCHtype,   !plot the CCH 
     +       x0,y0,height,width,info,BDT_FILE,
     +       QDT_FILE,date,recording,NHW,
     +       rec_num,show_control,
     +       rec_num_control,title,ICN,max_bin,
     +       show_conf_lim,
     +       show_single_shift,show_avg_shift,
     +       CONTROL_cl,show_2_sd,show_3_sd,
     +          scaledup,min_bin,rmagnify,IDs,ITAL)

        pixels_per_bin =width/FLOAT(total_num_bins)             !calculate width of bin


*       ***** designate the bins used as feature range for statistics: *****
*
c       call ps_fill_color(0.1,0.1,0.1)                 !fill color = gray
c       call ps_text_color (0.1,0.1,0.1)
        call ps_line_type(DOT)
c       call ps_interior_style(INT_SOLID,1)
        call ps_rectangle((x0+((IPKS-1)*pixels_per_bin)),
     +           y0+50.,(x0+(IPKF*pixels_per_bin)),360.)
        write (text,'(I3)') IPKS
        write (c_IPKF,'(I3)') IPKF
        call ps_character_height(.030)
        call ps_character_width(.010)
        call ps_text2d((450.+((IPKS-1)*pixels_per_bin)),
     +                 y0+50.,text//char(0),PS_ANNOTATION_TEXT)
        call ps_text2d((x0+(float(IPKF)*pixels_per_bin)),
     +                  y0+50.,c_IPKF//' (selected feature)'//
     +                  char(0),PS_ANNOTATION_TEXT)

*       ***** designate the bins used as background range for statistics: *****
*
c       call ps_fill_color(0.5,0.5,0.5)                 !fill color = gray
c       call ps_text_color(0.5,0.5,0.5)
c       call ps_line_color (0.0,0.0,1.0)
c       call ps_interior_style(INT_SOLID,1)
        call ps_rectangle((x0+((IBKS-1)*pixels_per_bin)),
     +           330.,(x0+(IBKF*pixels_per_bin)),340.)
c       call ps_interior_style(INT_HOLLOW,1)
        write (text,'(I3)') IBKS
        write (c_IBKF,'(I3)') IBKF
        call ps_character_height(.030)
        call ps_character_width(.010)
        call ps_text2d((450.+((IBKS-1)*pixels_per_bin)),
     +                 330.,text//char(0),PS_ANNOTATION_TEXT)
        call ps_text2d((x0+(float(IBKF)*pixels_per_bin)),
     +                  330.,c_IBKF//' (background)'//
     +                  char(0),PS_ANNOTATION_TEXT)
        call ps_line_type(SOLID)


        call ps_character_height (.03)
        call ps_character_width (.01)
        goto 111
c        text='(binwidth = '//stats_bw//' ms)'
cc        print '(A)',text
c        call strlength(text,LEN(text),l)
c          call ps_text2d(300.,150.,text(1:l)//char(0),
c     +       PS_ANNOTATION_TEXT)
 111    text='k value  = '//czk
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,150.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)
        text='detectability index (d) = '//cdet
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,130.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)
        text='visibility index = '//cvis
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,110.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)
        text='significance of k value = '//cprobk
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,90.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)
        text='time lag to feature = '//czlat//' msec.'
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,70.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)
        text='half width of feature = '//chalfwd//' msec.'
        call strlength(text,LEN(text),l)
          call ps_text2d(50.,50.,text(1:l)//char(0),
     +       PS_ANNOTATION_TEXT)

        call ps_character_height(.036)
        call ps_character_width(.012)

        call ps_text2d(1000.,700.,REFtext//char(0),
     +                  PS_ANNOTATION_TEXT)
        call ps_text2d(1175.,700.,TARtext//char(0),
     +                  PS_ANNOTATION_TEXT)
        call ps_text2d(1325.,700.,'ID code'//char(0),
     +                  PS_ANNOTATION_TEXT)
        call ps_text2d(1000.,680.,REFtype//char(0),
     +                  PS_ANNOTATION_TEXT)
        call ps_text2d(1175.,680.,TARtype//char(0),
     +                  PS_ANNOTATION_TEXT)
        call ps_text2d(1325.,680.,'resp. type'//char(0),
     +                  PS_ANNOTATION_TEXT)

        call ps_character_height(.030)
        call ps_character_width(.010)
        if(show_single_shift.eq.'y')then
           text='(single shift-control displayed)'
        else if(show_avg_shift.eq.'y')then
           text='(averaged shift-control displayed)'
           if(qdt_file_version.eq.'8') text='(cth control displayed)'
        else
           text = ' '
        end if
        call strlength(text,LEN(text),l)
        call ps_text2d(x0+220.,y0-70.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
        if(show_conf_lim.eq.'y')then
           if(show_2_sd.eq.'y')then
              text='(confidence limits:  +/- 2 s.d.)'
           else if(show_3_sd.eq.'y')then
              text='(confidence limits:  +/- 3 s.d.)'
           else
              text = ' '
           end if
           call strlength(text,LEN(text),l)
           call ps_text2d(x0+220.,y0-90.,text(1:l)//char(0),
     +                  PS_ANNOTATION_TEXT)
        end if
        
        call sc_label2_ps (fildes,x0+220,y0-70,show_single_shift
     +       ,show_avg_shift,show_2_sd,show_3_sd)
        call ps_character_height(.030)
        call ps_character_width(.010)
        y=670.
        do i = 1, MAX_PERTURB
           if((per_results_REF(i).ne.' ').or.
     +              (per_results_TAR(i).ne.' '))then
              y = y - 20.
              call strlength(per_results_REF(i),
     +                   LEN(per_results_REF(i)),l)
              call ps_text2d(1000.,y,per_results_REF(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
              call strlength(per_results_TAR(i),
     +                   LEN(per_results_TAR(i)),l)
              call ps_text2d(1175.,y,per_results_TAR(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
              call strlength(per_text_abbrev(i),
     +                   LEN(per_text_abbrev(i)),l)
              call ps_text2d(1325.,y,per_text_abbrev(i)(1:l)//
     +                  char(0),PS_ANNOTATION_TEXT)
           end if
        end do

        call strlength(ps_file,LEN(ps_file),l)
        call ps_fclose(ps_file(1:l)//char(0))
          if((task.eq.'p').or.(task.eq.'v'))then
             if(task.eq.'p')then
                isys=SYSTEM('lpr '//                    !print a hard copy of the ps file
     +               ps_file(1:l)//char(0))       
             else
                isys=SYSTEM('gv temp.ps'//char(0))              !pop-up display of the ps file
             end if
             isys=SYSTEM('rm temp.ps'//char(0))                 !remove the temporary file
          end if


          return
          end

      end module mod_print_and_write_routines
