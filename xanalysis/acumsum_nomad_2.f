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

      module mod_acumsum_nomad_2
      contains
*       filename = acumsum_nomad.f
*
*       date of last revision = 07-aug-2003     lss
*
*       link with x2000 or x2000_v2 code
*
*       This subroutine of analyze_data (which itself is a subroutine of x2002) asks
*               the user for a range of bins to be used to define the average value.
*               The user enters the extent of this range with two clicks of the mouse --
*               once on the beginning bin and then again on the ending bin of the range.
*
*
*
      SUBROUTINE qsum(fildes,mouse,IHIST,BINVAL,ICN,
     +     COEF,
     +     BDT_FILE,QDT_FILENAME,date,
     +     recording,NHW,rec_num,title,
     +     REFcode,
     +     TARcode,WINDOW1,WINDOW2,mode,
     +     IDs,ITAL)
*
*
      use mod_clear_routines
      use mod_locate_region
      use mod_miscellaneous_subroutines
      use mod_new_draw_button
      use mod_print_and_write_routines
      include 'x2000parameter.defs'

      integer IHIST(101),KDAT(101),ICN(6)
*
      real COEF,BINVAL,x_dc,y_dc,z_dc,xloc,yloc,zloc,
     +     CDAT(101)
*
      integer*4 mouse,first_bin,last_bin,total_bins
      integer*4 close
      INCLUDE 'gopen_type.defs'
      
      integer region,rec_num,REFcode,TARcode,IDs(MAX_NUM_CODES),
     +     ITAL(MAX_NUM_CHAN)
*
      character*120 directions,title,text
      character*130 text130
      character*30 label1,label2
      character*(*) BDT_FILE,QDT_FILENAME
      character*15 screen
      character*11 date
      character*10 c_ITAL
      character*6 c_ATOP,c_ABOT
      character*7 c_ICN(6),c_binwd
      character*5 c_rec_num
      character*4 c_NHW
      character*2 recording,mode
      character*1 task
      character*5 WINDOW1
      character*5 WINDOW2
*

*       ***** suppress unused variable warnings *****

      if(.false.)print *,mode
      if(.false.)print *,mouse
      if(.false.)print *,window1
*
*
      binwd = (float(NHW)/50.)
      SUMX=0.0
      SUMC=0.0
*
c       PRINT '(2X,''COEF:'',F10.4,'' ICN(2)='',I8)',
c     +         COEF,ICN(2)
      KDAT=IHIST                !load CCH data into KDAT for later manipulation
*
      
*      directions='SELECT FIRST BIN FOR DEFINING RANGE FOR AVERAGE'
*      call text2d(fildes,100.,55.,directions//char(0),
*     +     ANNOTATION_TEXT,0)
*      call make_picture_current (fildes)
*      call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
*      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
*      first_bin=INT((xloc-100.)/5)+1
*      if(first_bin.lt.1)first_bin=1
      first_bin=1
      call clear1(fildes,2)     !clear bottom of window
*
 100  directions='SELECT LAST BIN FOR DEFINING RANGE FOR AVERAGE'
      call text2d(fildes,100.,55.,directions//char(0),
     +     ANNOTATION_TEXT,0)
      call make_picture_current(fildes)
      call request_locator(fildes,1,2e9,valid,x_dc,y_dc,z_dc)
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc)
      last_bin=INT((xloc-100.)/5)+1
      if(last_bin.gt.101)last_bin=101
      if(last_bin.lt.1)last_bin=1
      if(last_bin.eq.1) goto 100
      call clear1(fildes,2)     !clear bottom of window
*
*
*       ***** color the bins used as control range for Q_SUM: *****
*
      call fill_color(fildes,0.,0.,0.) !fill color = black
      call interior_style(fildes,INT_SOLID,1)
      call rectangle(fildes,(100.+((first_bin-1)*5.)),
     +     310.,(100.+(last_bin*5.)),350.)
      call make_picture_current(fildes)
      call interior_style(fildes,INT_HOLLOW,1)
*
      total_bins=(last_bin-first_bin)+1 !determine total number of bins included
                                ! within the range
      do j = first_bin,last_bin !calc total of bin values within the range
         SUMX=SUMX + KDAT(j)
      end do
c        print '(''SUMX='',f7.2)',SUMX
*
      average = SUMX/total_bins !average value of range bin
c        print '(''average  '',f9.1)',average
      DO I=1,101
         SUMC=SUMC+(float(KDAT(I))-average) !subtract the average from each bin of the CCH
         CDAT(I)=SUMC           ! and keep a running total of the adjusted bin
*        IDAT(I)=AINT(CDAT(I))                  ! values; store the running totals in CDAT()
      end do
c       print '(''SUMC='',f7.2)',SUMC
*
*
*       scale everything used from here on to # of impulses/trigger event...
*
      do  I=1,101
         CDAT(I) = CDAT(I)/float(ICN(2)) !ICN(2) = total number of reference spikes
c        print '(''CDAT('',I3,'')= '',F9.7)',I,CDAT(I)
      end do
      average = average/float(ICN(2))
c       print '(''scaled average = '',f9.5)',average
*
*
*
*       Q_SUM CALC COMPLETE..DEFINE Y AXIS AND PLOT
*
      TOP=CDAT(1)
      BOT=CDAT(1)
*
      DO I=2,101
         TOP=AMAX1(TOP,CDAT(I)) !TOP = maximum value in CDAT
         BOT=AMIN1(BOT,CDAT(I)) !BOT = minimum value in CDAT
      end do
      print '(''TOP = '',f10.5)',TOP
      print '(''BOT = '',f10.5)',BOT
*
      BIG=AMAX1(ABS(TOP),ABS(BOT)) !BIG = the larger of the 2 absolute values for
                                !   TOP and BOT
      BORDER=0.1*(DIM(BIG,-BIG)) !define a border according to the total spread
                                !   to be used for the Q_SUM plot 
      ATOP=BIG+BORDER           !maximum and minimum excursions for the Q_SUM 
      ABOT=-BIG-BORDER          !   plot itself (border included)
*
c$$$    call strlength(USER,LEN(USER),m)
c$$$    isys=SYSTEM ('xwcreate -wmdir /dev/screen/'//USER(1:m)//
c$$$     +                      ' -title '//WINDOW2//
c$$$     +                      ' -geometry =750x500-5-200 -r'//char(0))
c$$$            fildes2=gopen ('/dev/screen/'//USER(1:m)//'/'//
c$$$     +                      WINDOW2//char(0),
c$$$     +              OUTINDEV,'sox11'//char(0),INIT)
      fildes2=gopen (750,500,-5,-200,WINDOW2//char(0))
*
*       set the view_window to Q_SUM settings:
      call vdc_extent(fildes2,0.0,0.0,0.0,1.25,1.0,0.0)
      call shade_mode(fildes2,IOR(INIT,CMAP_NORMAL),0)
c        call mapping_mode(fildes,1)
      call mapping_mode(fildes2,1)
      call view_port(fildes2,.15*1.25,.17,.85*1.25,.75)
      call view_window(fildes2,85.,ABOT,590.,ATOP) !xmin,ymin,xmax,ymax
      call background_color(fildes2,1.,1.,1.) !background color = white
      call clear(fildes2)
      call line_color(fildes2,0,0.,0.) !black
      call text_color(fildes2,0.,0.,0.) !black
*       call make_picture_current(fildes2)
*        print '(''ATOP='',f9.5,'' ABOT='',f12.5)',ATOP,ABOT
*
*        scale = 250./(ATOP-ABOT)
*        ht=float(AINT((ATOP-ABOT)*scale))
*        call clear1(fildes,1)
*        if(BOT.lt.0)IDAT=IDAT+AINT(ABS(BOT))
*
      call character_height(fildes2,.035)
      call character_width(fildes2,.012)
      write (c_ATOP,'(f5.3)') ATOP
      write (c_ABOT,'(f5.3)') ABOT
      call strlength(c_ATOP,LEN(c_ATOP),l)
      call text2d(fildes2,45.,ATOP,c_ATOP(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
      call strlength(c_ABOT,LEN(c_ABOT),l)
      call text2d(fildes2,45.,ABOT,c_ABOT(1:l)//char(0),
     +     ANNOTATION_TEXT,0)

C       *****  draw axes from 85,ATOP to 85,ABOT to 590,ABOT***
      call line_type(fildes2,SOLID)
c       call line_color(fildes2,0.,0.,0.)               !black
c       call background_color(fildes2,1.,1.,1.)         !background color is white
      call move2d(fildes2,85.,ATOP)
      call draw2d(fildes2,85.,ABOT)
      call draw2d(fildes2,590.,ABOT)
c       *****  draw a broken line from 1253,200 to 1253,700 (middle bin) ****
      call line_type(fildes2,DOT)
      call move2d(fildes2,338.,ABOT)
      call draw2d(fildes2,338.,ATOP)
      call line_type(fildes2,SOLID)
C
      call move2d(fildes2,85.,ABOT)
      IZ=101
      IX=85
      DO I=1,IZ
         J=I+1
         IF(J.EQ.102)cycle
         TEMP=CDAT(J)
         IX2=IX+5
         call move2d(fildes2,float(IX),CDAT(I))
         call draw2d(fildes2,float(IX2),CDAT(I))
         call draw2d(fildes2,float(IX2),TEMP)
         IX=IX2
      end do
      call make_picture_current(fildes2)
c
c
C       PLOT CONFIDENCE LIMITS ON Q_SUM
C
C
      IF (COEF.GT.0.9)
     +     print '(''Coef. of var. .GT. 0.9.'')' 
      ZIX=85.+((last_bin-1)*5.) !physical location (within the display window) of the last
                                ! bin of the range
      YC1=0.0
      YCM1=0.0
      LPCT=101-last_bin         !# OF BINS TO HAVE CONFIDENCE BAND
      U=(total_bins*BINVAL)     !TOTAL TIME OF CONTROL PERIOD (BINS SPANNED IN SEC.)
      ZMEAN= BINVAL/average
      C4=(1.0/6.0)-(COEF**4/6.0)
      C2=COEF**2
      UM=U*ZMEAN
      ISKIP=0                   !WHEN =1 SKIP V2 CALC BELOW
      DO 1500 J=1,LPCT
         V1=ICN(2)*(((C2*J*BINVAL)/ZMEAN)
     +        +C4+(((COEF*J*BINVAL)**2)/UM))
         V2=(ICN(2)*(J*BINVAL))/ZMEAN
         IF(V1.LT.V2)ISKIP=1
         VX=V2
         IF(ISKIP.EQ.1)VX=V1
         if (COEF.GT.0.9) VX=V2
         CURSIG=((3*SQRT(VX))/ICN(2)) 
C       
C       DRAW CONFIDENCE BAND ELEMENTS ABOVE AND BELOW THE JTH BIN IN LOOP
C
c       DEBUG LINES
c
c
c       PRINT 3000,COEF,U,average,ICN(2),BINVAL,V2,CURSIG
c3000   FORMAT (2X,'COEF:',F10.4,'U:',F10.4,'AV:',F10.4,'ICN(2):',I5,/,
c     + 2X,'BINVAL:',F10.4,'V2:',F10.4,'CURSIG:',F10.4)
c
c
c
c
         ZIX2=ZIX+5.
         YC2=CURSIG
         YCM2=-CURSIG
         if (CURSIG .gt. ATOP) goto 1000
         call move2d(fildes2,ZIX,YC1)
         call draw2d(fildes2,ZIX2,YC2) !draw line above 0
         call move2d(fildes2,ZIX,YCM1)
         call draw2d(fildes2,ZIX2,YCM2) !draw line below 0
         ZIX=ZIX2
         YC1=YC2
         YCM1=YCM2
 1500 CONTINUE
*1000    call make_picture_current(fildes2)
 1000 call vdc_extent(fildes2,0.0,0.0,0.0,1.0,1.0,0.0)
      call view_port(fildes2,0.,0.,1.,1.)
      call view_window(fildes2,85.,10.,590.,500.) !xmin,ymin,xmax,ymax
c        call text_color(fildes2,1.,0.,1.)              !black
*

      text='change in'
      call text2d(fildes2,95.,300.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='# spikes'
      call text2d(fildes2,95.,285.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='per'
      call text2d(fildes2,95.,270.,text//char(0),
     +     ANNOTATION_TEXT,0)
      text='REF event'
      call text2d(fildes2,95.,255.,text//char(0),
     +     ANNOTATION_TEXT,0)



      write (c_rec_num,'(I5)')rec_num
*
      write (c_NHW,'(I4)')NHW
*
      write (c_binwd,'(F7.1)')binwd
      do i = 1,6                !convert to character format
         write (c_ICN(i),'(6I7)')ICN(i)
         call remove_all_blanks(c_ICN(i),LEN(c_ICN(i)))
      end do
*
      text='.bdt filename:  '//BDT_FILE
      call text2d(fildes2,85.,480.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l)
      text='.qdt filename:  '//
     +     QDT_FILENAME(1:l)
      call text2d(fildes2,85.,465.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      text='date of experiment:  '//date
      call text2d(fildes2,85.,450.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      text='recording #:  '//recording
      call text2d(fildes2,85.,435.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      text='record # in .qdt file:  '//c_rec_num
      call text2d(fildes2,85.,65.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      call strlength(c_ICN(1),LEN(c_ICN(1)),l_1)
      write (c_ITAL,'(I10)') ITAL(IDs(ICN(1)))
      call remove_all_blanks(c_ITAL,LEN(c_ITAL))
      call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

      text='R: ID = '//c_ICN(1)(1:l_1)//
     +     ' (# events in spike train = '//
     +     c_ITAL(1:l_ITAL)//')'
      call strlength(text,LEN(text),l)
      call text2d(fildes2,85.,50.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
*
      call strlength(c_ICN(3),LEN(c_ICN(3)),l_3)
      call strlength(c_ICN(4),LEN(c_ICN(4)),l_4)
      write (c_ITAL,'(I10)') ITAL(IDs(ICN(3)))
      call remove_all_blanks(c_ITAL,LEN(c_ITAL))
      call strlength(c_ITAL,LEN(c_ITAL),l_ITAL)

      text='T: ID = '//c_ICN(3)(1:l_3)//
     +     ' (# events in spike train = '//
     +     c_ITAL(1:l_ITAL)//')'
      call strlength(text,LEN(text),l)
      call text2d(fildes2,85.,35.,text(1:l)//char(0),
     +     ANNOTATION_TEXT,0)
*
      text='binwidth:  '//c_binwd//' msec.'
      call text2d(fildes2,85.,20.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      text='-'//c_NHW
      call text2d(fildes2,161.,75.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      text=c_NHW//' msec'
      call text2d(fildes2,500.,75.,text//char(0),
     +     ANNOTATION_TEXT,0)
*
      call character_height(fildes2,.060)
      call character_width(fildes2,.020)
      text130 = 'QSUM: '//title
      call strlength(text130,LEN(text130),l)
      call text2d (fildes2,250.,420.,
     +     text130(1:l),ANNOTATION_TEXT,0) !'REF>TAR'
      call character_height(fildes2,.030)
      call character_width(fildes2,.010)
      call make_picture_current(fildes2)
*
 3900 label1='CLICK HERE'
      label2='TO CONTINUE'
      call draw_button(fildes,5.,5.,205.,80.,
     +     label1,10.,50.,label2,10.,20.)
      label1='   PRINT'
      call draw_button(fildes,5.,90.,205.,165.,label1,10.,120.,'',0.,0.) 
      label1='   WRITE'
      call draw_button(fildes,275.,90.,475.,165.,label1,280.,120.,'',0.,
     +     0.)                                             
 4000 call make_picture_current(fildes)
      call request_locator(fildes,1,2e9,ivalid,x_dc,y_dc,z_dc)
      if (ivalid.eq.4) goto 5000 !close button was pressed on fildes2
      call vdc_to_wc(fildes,x_dc,y_dc,z_dc,xloc,yloc,zloc) !translate into fp coords.
      screen = 'qsum'
      call locate_region(xloc,yloc,region,screen)
*         print '(''region= '',I3)',region
      if(.false.)then
*
      else if(region.eq.5       !PRINT or WRITE the Q_SUM window
     +        .or.region.eq.6.or.region.eq.205)then
         if(region.eq.5)task='p'
         if(region.eq.6)task='w'
         if(region.eq.205)task='v'
         call print_or_write_QSUM(task,IHIST,BINVAL,ICN,COEF,
     +        BDT_FILE,QDT_FILENAME,date,
     +        recording,NHW,rec_num,title,
     +        REFcode,
     +        TARcode,first_bin,last_bin,IDs,ITAL)

         call clear_bottom(fildes)
         goto 3900
      else if(region.eq.10)then
         goto 5000
      else
         goto 4000
      end if
*
 5000 call clear_bottom(fildes)
      close = gclose(fildes2)
      RETURN
      END
      
      end module mod_acumsum_nomad_2
