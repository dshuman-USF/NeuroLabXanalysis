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


      module mod_showCTHs_v3
      contains
*       SUBROUTINE TO CREATE HISTOGRAMS in starbase window
*       using starbase library calls... derived from autodis.f
*
*       filename = showCTHs_v3.f
*
*       date of latest revision = 21-Sep-2006  lss
*
*       may-2002        lss
*        cardiac pulses now calculated as a CCH with no offset
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
*       *** DIRECT POINTERS are used with these arrays: exclude, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*
*
*
*       new_plot incorporated  01-apr-99
*
*       LINK WITH x2000_v2 code
*
*       This subroutine will display a cell's respiratory CTH with phrenic CTH overlay,
*               the normalized respiratory CTH with normalized phrenic overlay, and
*               the cardiac CCH with cardiac pulse overlay.
*
*       def'ns of iflag:        0 --> display titles and labels for overlays only
*                               1 --> not used
*                               2 --> show only phrenic overlay
*                               3 --> show only cardiac overlay
*
      SUBROUTINE showCTHs (iflag,c_format,fildes,cellCTH,
     +     phrenicCTH,
     +     norm_cellCTH,norm_phrenicCTH,cardCCH,
     +     card_overlay,bw_p,bw_n,bw_c,
     +     bwtext1,bwtext2,bwtext3,
     +     bwtext4,bwtext5,
     +     REFcode,REFevents_p,REFevents_n,REFevents_c,
     +     TARcode,TARevents,
     +     BDT_FILE,TAR_ACH,TAR_ACH_1,TARname,
     +     QDT_FILENAME,sort,
     +     mode,resp_type,AP,RL,dep,
     +     AA,AA_results,AA_text1,AA_applied,
     +     STA,STA_results,STA_text1,STA_applied,
     +     STA_res_text,per_results,per_text_abbrev,
     +     zmodsig,zmodsig2,ETA2,coef,card_type,card,
     +     DELTA2,tedfactor,fiveHT,mean_rISI,sd_rISI,
     +     c_MAX_INT,num_rej_ISI,bw_ach1,bw_ach2,
     +     REFevents_ach1,REFevents_ach2,IDs,ITAL,
     +     cardiac_pls,STIM_OFFSET,NORM_OFFSET,mean_E,
     +     comm,show_pulses)
*
      use mod_clear_routines
      use mod_miscellaneous_subroutines
      use mod_new_plot
      INCLUDE 'x2000parameter.defs'

      integer cellCTH(101),phrenicCTH(101),norm_cellCTH(101),
     +     norm_phrenicCTH(101),cardCCH(101),card_overlay(101),
     +     TAR_ACH(101),TAR_ACH_1(101),scaled_hist(101),
     +     NTOP1,NTOP2,NTOP3,NTOP4,NTOP5,
c     +           sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      integer cardiac_pls
      real STIM_OFFSET,NORM_OFFSET,mean_E
*
*
      INCLUDE 'gopen_type.defs'
      integer REFcode,TARcode,REFevents_p,REFevents_n,
     +     REFevents_c,TARevents,REFevents_ach1,REFevents_ach2
*
      character*(*) fiveHT,mean_rISI,sd_rISI,c_MAX_INT
      character*(*) BDT_FILE
      character*(*) QDT_FILENAME,comm
      character*(*) bwtext1,bwtext2,bwtext3,bwtext4,bwtext5,
     +     AA_text1(MAX_AA),STA_text1(MAX_STA),
     +     STA_results(MAX_STA),
     +     per_text_abbrev(MAX_PERTURB)
      character*(*) num_rej_ISI
      character*(*) resp_type,STA_res_text(122)
      character*(*) per_results(MAX_PERTURB)
      character*(*) AP,RL,dep,ETA2,coef,DELTA2
      character*(*) TARname
      character*(*) AA_results(MAX_AA),zmodsig,zmodsig2,card_type
      character*(*) mode,tedfactor
      character*(*) AA(MAX_AA),STA(MAX_STA), AA_applied,STA_applied,
     +     sort, card,show_pulses

      integer*4 ICN0(6)
      character*1 show_control
      character*3 CTH,ACH,CCC
      character*4 c_TARcode,c_REFcode
      character*7 REFeventstext,TAReventstext
      character*20 c_format
      character*10 c_ITAL
      character*150 REFline,TARline
      character*80 text80
      character*50 text,text2
      character*10 c_rate_p,c_rate_n,c_rate_c,c_rate_ach1,c_rate_ach2
      character*60 title,message
*
*
*       ***** suppress unused variable warnings *****

      if(.false.)print *,c_format
      if(.false.)print *,STA_applied
      if(.false.)print *,fildes2 !suppress unused variable warning

      ICN0 = 0
      min_bin=0
      max_bin=0
      rmagnify=0.0
      show_control = 'y'        !default value - always show the overlay for a CTH

*       *****  set screen parameters  ***************************************
      call text_font_index(fildes,6)
      call vdc_extent(fildes,0.0,0.0,0.0,1.0,1.0,0.0)
      call mapping_mode(fildes,1)
      call view_window(fildes,1.,1.,1700.,500.)
*       *********************************************************************
*
      write (TAReventstext,'(I7)') TARevents
      write (REFeventstext,'(I7)') REFevents_p
      write (c_REFcode,'(I4)') REFcode
      write (c_TARcode,'(I4)') TARcode
      title=' '
      title = 'CTHs for cell '//TARname
      message = '(ID code = '//c_TARcode//')'
      if(TARname.eq.' ')title = 
     +     'CTHs for cell with ID code = '//c_TARcode
      if(iflag.eq.0)then
         if(cardiac_pls.ne.0)then
            title = 'CTHs of phrenic, norm. phrenic, & cardiac pulse'
         else
            title = 'CTHs of phrenic & norm. phrenic'
         end if
      end if
      if(iflag.eq.2)title = 'CTH of phrenic'
      if(iflag.eq.3)title = 'CCH of cardiac pulse'
      TARline=' '
      CTH = 'CTH'
      ACH = 'ACH'
      CCC = 'CCC'
*
*
*       *****   label the CTHs according to binwidths, input filename, & codes: *****
*
      call text_font_index (fildes,6)
      call text_color(fildes,0.0,0.0,0.0) !text color is black
      call background_color(fildes,1.,1.,1.) !background color is white
      call character_width(fildes,0.020)
      call character_height(fildes,0.1)
      call clear(fildes)        !clear the display
      call text2d (fildes,10.,450.,title//char(0),
     +     ANNOTATION_TEXT,0)
      if(TARname.ne.' ')then
         call character_height(fildes,.060)
         call character_width(fildes,.015)
         call text2d(fildes,10.,425.,message//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      if(iflag.eq.2)goto 10     !show phrenic overlay only
      if(iflag.eq.3)goto 20     !show cardiac overlay only
*
      call character_width(fildes,0.010)
      call character_height(fildes,0.030)
      call line_type(fildes,DOT)
      call move2d (fildes,1.,500.)
      call draw2d (fildes,1700.,500.)
      call character_width(fildes,0.005)
      call character_height(fildes,0.025)
      text = '<-- binwidths -->'
      call strlength(text,LEN(text),l)
      call text2d (fildes,300.,80.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call character_width(fildes,0.010)
      call character_height(fildes,0.030)
      call strlength(bwtext1,LEN(bwtext1),l)
      call text2d (fildes,600.,80.,bwtext1(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call strlength(bwtext2,LEN(bwtext2),l)
      call text2d (fildes,1020.,80.,bwtext2(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      if(cardiac_pls.ne.0)then
         call strlength(bwtext3,LEN(bwtext3),l)
         call text2d (fildes,1450.,80.,bwtext3(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      call character_width(fildes,0.005)
      call character_height(fildes,0.025)
      call remove_all_blanks(bwtext4,LEN(bwtext4))
      call strlength(bwtext4,LEN(bwtext4),l)
      call text2d (fildes,75.,82.,bwtext4(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call strlength(bwtext5,LEN(bwtext5),l)
      call text2d (fildes,180.,82.,bwtext5(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call character_width(fildes,0.010)
      call character_height(fildes,0.030)

      call remove_all_blanks(c_REFcode,LEN(c_REFcode))
      call strlength(c_REFcode,LEN(c_REFcode),l_code)
      write (c_ITAL,'(I10)') ITAL(IDs(REFcode))
      call remove_all_blanks(c_ITAL,LEN(c_ITAL))
      call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

      REFline='R: ID='//c_REFcode(1:l_code) //
     +     ' (# E pulses in spike train = '//c_ITAL(1:l_ITAL)//')'
      call strlength(REFline,LEN(REFline),l)
      call text2d (fildes,1.,45.,REFline(1:l)//char(0),
     +     ANNOTATION_TEXT,0)

      call remove_all_blanks(c_TARcode,LEN(c_TARcode))
      call strlength(c_TARcode,LEN(c_TARcode),l_code)
      call remove_all_blanks(TAReventstext,LEN(TAReventstext))
      call strlength(TAReventstext,LEN(TAReventstext),l_events)
      write (c_ITAL,'(I10)') ITAL(IDs(TARcode))
      call remove_all_blanks(c_ITAL,LEN(c_ITAL))
      call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)
      TARline='T: ID='//c_TARcode(1:l_code) //
     +     ' (# events in spike train = '//c_ITAL(1:l_ITAL)//')'
      call strlength(TARline,LEN(TARline),l)
      call text2d (fildes,1.,25.,TARline(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call remove_all_blanks(REFeventstext,LEN(REFeventstext))
      call strlength(REFeventstext,LEN(REFeventstext),la_events)
      text = 'sum of all bins = '//TAReventstext(1:l_events)//
     +     '; # cycles shown = '//REFeventstext(1:la_events)
      call strlength(text,LEN(text),l)
      call text2d (fildes,1.,5.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      

      call strlength(BDT_FILE,LEN(BDT_FILE),l_bdt)
      call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
      text = BDT_FILE(1:l_bdt)//' / '//QDT_FILENAME(1:l)
      call strlength(text,LEN(text),l)
      call text2d(fildes,1.,65.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)

      if(sort.eq.'c')then
         text='ACHs: (all-order; clean)'
      else if(sort.eq.'m')then
         text='ACHs: (all-order; messy)'
      else
         text='ACHs:'
      end if   


      if(iflag.eq.0)text=' '
      call text2d (fildes,10.,220.,text//char(0),
     +     ANNOTATION_TEXT,0)

      text = 'RESPIRATORY CTH'
      call text2d (fildes,490.,60.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text = '(phrenic overlay)'
      call text2d (fildes,520.,40.,text//char(0),
     +     ANNOTATION_TEXT,0)
      if (norm_cellCTH(101).gt.0) then
         text = 'PHASE NORMALIZED RESP. CTH'
      else
         text = 'NORMALIZED RESP. CTH'
      end if
      call text2d (fildes,875.,60.,text//char(0),
     +     ANNOTATION_TEXT,0)
      if (norm_cellCTH(101).gt.0) then
         text = '(phase norm. phrenic overlay)'
      else
         text = '(norm. phrenic overlay)'
      end if
      call text2d (fildes,900.,40.,text//char(0),
     +     ANNOTATION_TEXT,0)
      if(cardiac_pls.ne.0)then
         text = 'CARDIAC CCH'
         call strlength(text,LEN(text),l)
         call text2d (fildes,1400.,60.,text(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
         text = '(cardiac pulse overlay)'
         call strlength(text,LEN(text),l)
         call text2d (fildes,1350.,40.,text(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
      else
         text = 'CARDIAC CCHs NOT CALCULATED'
         call strlength(text,LEN(text),l)
         call text2d (fildes,1300.,250.,text(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
      end if  
*       ***** find the maximum bin in each CTH and ACH and calculate the peak firing rate: *****
      
      NTOP1=cellCTH(1)
      NTOP2=norm_cellCTH(1)
      NTOP3=cardCCH(1)
      NTOP4=TAR_ACH(1)
      NTOP5=TAR_ACH_1(1)
      do I=2,100
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
      call text2d (fildes,420.,405.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      rate = (NTOP2*(1000.0/bw_n))/REFevents_n
      write (c_rate_n,'(F9.1)')rate
      call remove_all_blanks(c_rate_n,LEN(c_rate_n))
      call strlength(c_rate_n,LEN(c_rate_n),l)
      text = c_rate_n(1:l)//' spikes/sec'
      call strlength(text,LEN(text),l)
      call text2d (fildes,840.,405.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      if(cardiac_pls.ne.0)then
         rate = (NTOP3*(1000.0/bw_c))/REFevents_c
         write (c_rate_c,'(F9.1)')rate
         call remove_all_blanks(c_rate_c,LEN(c_rate_c))
         call strlength(c_rate_c,LEN(c_rate_c),l)
         text = c_rate_c(1:l)//' spikes/sec'
         call strlength(text,LEN(text),l)
         call text2d (fildes,1270.,405.,text(1:l)//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      call character_width(fildes,0.005)
      call character_height(fildes,0.025)
      rate = (NTOP4*(1000.0/bw_ach1))/REFevents_ach1
      write (c_rate_ach1,'(F9.1)')rate
      text = c_rate_ach1
      call remove_all_blanks(text,LEN(text))
      call strlength(text,LEN(text),l)
      call text2d (fildes,40.,203.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      rate = (NTOP5*(1000.0/bw_ach2))/REFevents_ach2
      write (c_rate_ach2,'(F9.1)')rate
      call remove_all_blanks(c_rate_ach2,LEN(c_rate_ach2))
      call strlength(c_rate_ach2,LEN(c_rate_ach2),l)
      text = c_rate_ach2(1:l)//' spikes/sec'
      call strlength(text,LEN(text),l)
      call text2d (fildes,160.,203.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call character_width(fildes,0.010)
      call character_height(fildes,0.030)

*       *****   scale the CTHs and display:     *****

 10   call new_plot(fildes,100,cellCTH,phrenicCTH,scaled_hist,CTH,420., !plot resp CTH and overlay
     +     100.,300.,404.,'','','','','',0,0,show_control,0,'',
     +     ICN0,0,'', '','',[integer::],'','','',0,0.,
     +     [integer::], [integer::])
      if(show_pulses.ne.'n')then
         call line_color(fildes,0.,0.,1.) !blue
         call text_color(fildes,0.,0.,1.)
         call line_type(fildes,DOT)
         call move2d(fildes,420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +        100.)
         call draw2d(fildes,420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +        100.+300.)
         call text2d(fildes,420.+(404.*(ABS(STIM_OFFSET)/(100.*bw_p))),
     +        80.,'E'//char(0),ANNOTATION_TEXT,0)
         call line_color(fildes,0.,0.,0.)
         call text_color(fildes,0.,0.,0.)
         call line_type(fildes,SOLID)
      end if
      if(iflag.eq.2)goto 100   

      call new_plot(fildes,100,norm_cellCTH,norm_phrenicCTH,scaled_hist, !plot norm CTH and overlay
     +     CTH,840.,100.,300.,404.,'','','','','',0,0,show_control, 
     +     0,'',ICN0,max_bin,'','','',[integer::],'','','',
     +     min_bin,0.,[integer::],[integer::])          
      if(show_pulses.ne.'n')then
         call line_color(fildes,0.,0.,1.) !blue
         call text_color(fildes,0.,0.,1.)
         call line_type(fildes,DOT)
         call move2d(fildes,840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +        100.)
         call draw2d(fildes,840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +        100.+300.)
         call text2d(fildes,840.+(404.*(ABS(NORM_OFFSET)/(100.*bw_n))),
     +        80.,'E'//char(0),ANNOTATION_TEXT,0)
         if(mean_E.ne.0 0)then
            call move2d(fildes,840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +           (100.*bw_n))),100.)
            call draw2d(fildes,840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +           (100.*bw_n))),100.+300.)
            call text2d(fildes,840.+(404.*((ABS(NORM_OFFSET)+mean_E)/
     +           (100.*bw_n))),80.,'I'//char(0),ANNOTATION_TEXT,0)
         end if
         call line_color(fildes,0.,0.,0.)
         call text_color(fildes,0.,0.,0.)
         call line_type(fildes,SOLID)
      end if

 20   if(cardiac_pls.ne.0)then
         call new_plot(fildes,100,cardCCH,card_overlay,scaled_hist,CCC, !plot cardiac CCH and overlay
     +        1270.,100.,300.,404.,'','','','','',0,0,show_control, 
     +        0,'',ICN0,max_bin,'','','',[integer::],'','','',
     +        min_bin,0.,[integer::],[integer::])       
      end if
      if(iflag.eq.3)goto 100
      if(iflag.eq.0)goto 100    !don't draw an ACH for phrenic
      call new_plot(fildes,100,TAR_ACH,TAR_ACH,scaled_hist,ACH,40.,100., !plot ACHs (smallest and next-to-greatest BW)
     +     100.,101.,'','','','','',0,0,show_control,0,'',ICN0,
     +     max_bin,'','','',[integer::],'','','',min_bin,0.,[integer::],
     +     [integer::])                        
      call new_plot(fildes,100,TAR_ACH_1,TAR_ACH_1,scaled_hist,ACH,160.,
     +     100.,100.,101.,'','','','','',0,0,show_control,0,'',
     +     ICN0,max_bin,'','','',[integer::],'','','',min_bin,0.,
     +     [integer::],[integer::])


      call character_height(fildes,.030)
      call character_width(fildes,.010)
      if(mode.eq.'cr')text='(data entry mode)'
      if(mode.eq.'ed')text='(edit mode)'
      if((mode.eq.'vt').or.(mode.eq.'jl'))
     +     text='(view-only mode)'
      call text2d(fildes,10.,480.,text//char(0),
     +     ANNOTATION_TEXT,0)          
      call character_height(fildes,.060)
      call character_width(fildes,.015)
      call text2d(fildes,10.,400.,resp_type//char(0),
     +     ANNOTATION_TEXT,0)
      call character_height(fildes,.030)
      call character_width(fildes,.010)

      if(mode.ne.'vt')then      !do not display the following info if only looking:
         text = 'AP, RL, D: '//AP//', '//
     +        RL//', '//dep//' mm.'
         call text2d(fildes,10.,380.,text//char(0),
     +        ANNOTATION_TEXT,0)
         y = 360.
         call text2d(fildes,10.,y,'AA:'//char(0),
     +        ANNOTATION_TEXT,0)
         do ii = 1,MAX_AA
            if(AA(ii).eq.'y')then
               text=AA_results(ii)//'  '//
     +              AA_text1(ii)
               call text2d(fildes,50.,y,text//char(0),
     +              ANNOTATION_TEXT,0)
               y = y-15.
            end if
         end do

         if(AA_applied.eq.'y')then
            y = y - 5.
         else
            y = y-20.
         end if
         call text2d(fildes,10.,y,'STA:'//char(0),
     +        ANNOTATION_TEXT,0)
         do ii = 1,MAX_STA
            if(STA(ii).eq.'y')then
               call strlength(STA_results(ii),
     +              LEN(STA_results(ii)),jj)
               text=STA_results(ii)(1:jj)//' '//
     +              STA_res_text(ichar(STA_results(ii)(1:1)))
               call text2d(fildes,60.,y,text//char(0),
     +              ANNOTATION_TEXT,0)
               call text2d(fildes,250.,y,STA_text1(ii)//char(0),
     +              ANNOTATION_TEXT,0)
               y = y-15.
            end if
         end do
         
         y = 480.
         call text2d(fildes,300.,y,'Responses:'//char(0),
     +        ANNOTATION_TEXT,0)
         m = 0
         x = 400.
         do ii = 1, MAX_PERTURB
            if(per_results(ii).ne.' ')then
c                   print '(''per_results('',I2,'') = '',A)',
c     +                 ii,per_results(ii)
               m=m+1
               if(m.eq.1)x=450.
               if(m.ne.1)x=x+250.
               call strlength(per_text_abbrev(ii),
     +              LEN(per_text_abbrev(ii)),jj)
               call strlength(per_results(ii),
     +              LEN(per_results(ii)),kk)
               text=per_text_abbrev(ii)(1:jj)//': '//per_results(ii)
               call text2d(fildes,x,y,text//char(0),
     +              ANNOTATION_TEXT,0)
            end if
            if(m.eq.4)then
               m = 0            !reset - display 4 responses per line
               y = y-15.
            end if
         end do

         call text2d(fildes,300.,425.,'Comments:'//char(0),
     +        ANNOTATION_TEXT,0)
         call strlength(comm,LEN(comm),l)
         call text2d(fildes,450.,425.,comm(1:l)//char(0),
     +        ANNOTATION_TEXT,0)

         if(cardiac_pls.ne.0)then
            if((card.eq.'c').or.(card.eq.'C'))then
               text='cardiac mod.:  visual --> YES'
            else if((card.eq.'n').or.(card.eq.'N'))then
               text='cardiac mod.:  visual --> NO'
            else
               text = 'cardiac mod.:  visual -->'
            end if
            call text2d(fildes,1250.,5.,text//char(0),
     +           ANNOTATION_TEXT,0)
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
      call text2d(fildes,820.,20.,text//char(0),
     +     ANNOTATION_TEXT,0)

      if((zmodsig2.eq.'r').or.(zmodsig2.eq.'R'))then
         text='resp. mod.:  BINARY --> YES'
      else if((zmodsig2.eq.'n').or.(zmodsig2.eq.'N'))then
         text='resp. mod.:  BINARY --> NO'
      else if((zmodsig2.eq.'ned').or.(zmodsig2.eq.'NED'))then
         text='resp. mod.:  BINARY --> NO (but insuff. data)'
      else
         text='resp. mod.:  BINARY --> '
      end if
      call text2d(fildes,820.,5.,text//char(0),
     +     ANNOTATION_TEXT,0)

      text='ETA = '//ETA2//
     +     ';  c.o.v. = '//coef
      call text2d(fildes,460.,20.,text//char(0),
     +     ANNOTATION_TEXT,0)
      if(cardiac_pls.ne.0)then
         text='(using cycles in CTH)'
         call text2d(fildes,465.,5.,text//char(0),
     +        ANNOTATION_TEXT,0)
      end if
*
      if(DELTA2.ne.' ')then     !do not display if no value for DELTA2
         if((card_type.eq.'c').or.(card_type.eq.'C'))then
            text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +           ' --> YES'
         else if((card_type.eq.'n').or.(card_type.eq.'N'))then
            text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +           ' --> NO'
         else
            text='DELTA2 ='//DELTA2//'; ted f. = '//tedfactor//
     +           ' --> insuff. data'
         end if
         call text2d(fildes,1250.,20.,text//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      if(fiveHT.ne.' ')then     !do not display this stuff if no value for fiveHT
         text2 = '('//fiveHT//')' !display 5HT and relativeISI values
         call remove_all_blanks(text2,LEN(text2))
         do i = 1,LEN(fiveHT)
            if(fiveHT(i:i).ne.' ')exit !find the first "real" character in the string
         end do
         if(fiveHT(i:i).eq.'-')then
            text80='5HT? --> YES '//text2
         elseif(fiveHT.eq.'No data')then
            text80='5HT?  (No data)'
         else
            text80='5HT? --> NO '//text2
         end if
         call upper_case(text80,LEN(text80))
         call text2d(fildes,1400.,480.,text80,
     +        ANNOTATION_TEXT,0)
         if(mean_rISI.eq.'No data')then
            text80='mean rel. ISI: No data'
         else
            text2=mean_rISI//' +/- '//sd_rISI
            call remove_all_blanks(text2,LEN(text2))
            text80='mean rel. ISI = '//text2
         end if
         call text2d(fildes,1400.,460.,text80,
     +        ANNOTATION_TEXT,0)
         text2 = '(>'//c_MAX_INT//'ms.)'
         call remove_all_blanks(text2,LEN(text2))
         text80=num_rej_ISI//' ISI rejected '//text2
         call remove_leading_blanks(text80,LEN(text80))
         call text2d(fildes,1400.,440.,text80,
     +        ANNOTATION_TEXT,0)
      end if

      call make_picture_current(fildes)
*
*
 100  call view_window(fildes,1.,1.,1500.,500.)
*
*
      call make_picture_current(fildes)
      RETURN
      END
      end module mod_showCTHs_v3
