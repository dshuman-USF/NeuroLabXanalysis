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

C       SUBROUTINE TO CREATE HISTOGRAMS in X11 window
c               using starbase library calls... derived from autodisx11.f
c
*       date of lastest revision = 26-Feb-99   lss
*
c       displays cycle duration distribution histograms during the cycle
*               selection process (selected cycles will be used to compute
*               both regular and normalized CTHs)
c
C       LINK WITH MAIN PROGRAM x200*
*
*
      SUBROUTINE DISPLAY_SELECTED_CYCLES (fildes,iohist,ihs,BINW,left)

      include 'x2000parameter.defs'

      integer ihs(101),iohist(101),IH(101)
      INCLUDE 'gopen_type.defs'
      real BINW,x0,x1,y0,y1
      character left
      if(.false.)print *,fildes2 !suppress unused variable warning

*
*       ***** set default screen values: (to draw on right side of screen)*****
*
      x0 = 1000.
      x1 = 1404.
      y0 = 90.
      y1 = 390.
*
      if(left.eq.'y')then       !draw histograms on left side of screen
         x0 = 1.
         x1 = 405.
         y0 = 90.
         y1 = 390.
      end if
*
c       *****  set screen parameters  ***************************************

*       call shade_mode(fildes,IOR(INIT,CMAP_NORMAL),0)
      call vdc_extent(fildes,0.0,0.0,0.0,1.25,1.0,0.0)
      call mapping_mode(fildes,1)
      call view_port(fildes,.05*1.25,.10,.95*1.25,.95)
c        call view_window(fildes,1.,85.,1700.,800.)
      call view_window(fildes,1.,85.,1700.,600.)
c       *********************************************************************
      if(BINW.eq.0.0)then
         call clear2(fildes)
         return
      end if
*
*
c
c       main  histogram
c       SCLFC=1.0
c
      call DISSCZ1 (iohist,IH,0,NTOP,SCLFC)
c
      call VECOVL21(IH,fildes,x0)
      CALL DISSCZ1(ihs,IH,1,NTOP,SCLFC)
      CALL VECHIZ1(IH,fildes,x0,y0,x1,y1) !main plot
      call make_picture_current(fildes)
c
c
c
c
      RETURN
      END
C
c
C       THIS ROUTINE DRAWS SIMPLE FRAME FOR main HISTOGRAM

      SUBROUTINE VECHIZ1(IH,fildes,x0,y0,x1,y1)
c
      include 'x2000parameter.defs'

      integer IH(101)
      INCLUDE 'gopen_type.defs'
      real x0,x1,y0,y1
      if(.false.)print *,fildes2 !suppress unused variable warning
c
c       *****  clear the screen  *****************************************
*        call clear2(fildes)
c
*
*       *****  draw axes from x0,y1 to x0,y0 to x1,y0  *****************
      call line_color(fildes,0.,0.,0.) !black
      call line_type(fildes,SOLID)
      call move2d(fildes,x0,y1)
      call draw2d(fildes,x0,y0)
      call draw2d(fildes,x1,y0)
c
c
      linetype=0
      call dohis21(IH,fildes,linetype,3,x0)
c        call make_picture_current(fildes)
      RETURN
      END
c
c
c       THIS SUBROUTINE DRAWS BARS (POLYGONS) IN HISTOGRAM
      SUBROUTINE dohis21(IH,fildes,linetype,linecolor,x0)
c
c
      include 'x2000parameter.defs'

      INCLUDE 'gopen_type.defs'
      integer IH(101),linecolor
      real x0,R,G,B
      if(.false.)print *,fildes2 !suppress unused variable warning

      if(.false.)then
*               !default color is magenta
      else if(linecolor.eq.0)then !black
         R=0.0
         G=0.0
         B=0.0
      else if(linecolor.eq.1)then !white
         R=1.0
         G=1.0
         B=1.0
      else if(linecolor.eq.2)then !red
         R=1.0
         G=0.0
         B=1.0
      else if(linecolor.eq.3)then !yellow
         R=1.0
         G=1.0
         B=0.0
      else if(linecolor.eq.4)then !green
         R=0.0
         G=1.0
         B=0.0
      else if(linecolor.eq.5)then !cyan
         R=0.0
         G=1.0
         B=1.0
      else if(linecolor.eq.6)then !blue
         R=0.0
         G=0.0
         B=1.0
      else if(linecolor.eq.7)then !magenta
         R=1.0
         G=0.0
         B=1.0
      else
         R=1.0
         G=0.0
         B=1.0
      end if

      call line_color(fildes,R,G,B)
      call line_type(fildes,linetype)
c        call move2d(fildes,85,90)
      IX=x0
      DO 100 I=1,100
         IH(I)=IH(I)+90
         J=I+1
         IF(J.EQ.102)GOTO 100
         ITEMP=IH(J)
         ITEMP=ITEMP+90
         IX2=IX+4
         if(IH(I).eq.90)goto 90 !bin value is 0 - don't draw it
C       **********************************
*        call move2d(fildes,float(IX),float(IH(I)))
         call move2d(fildes,float(IX),90.)
         call draw2d(fildes,float(IX),float(IH(I)))
         call move2d(fildes,float(IX+1),90.)
         call draw2d(fildes,float(IX+1),float(IH(I)))
         call move2d(fildes,float(IX+2),90.)
         call draw2d(fildes,float(IX+2),float(IH(I)))
         call move2d(fildes,float(IX+3),90.)
         call draw2d(fildes,float(IX+3),float(IH(I)))
*        call draw2d(fildes,float(IX2),float(IH(I)))
*        call draw2d(fildes,float(IX2),float(ITEMP))
 90      IX=IX2
 100  CONTINUE
      call line_color(fildes,0.,0.,0.) !black
c        call make_picture_current(fildes)
      RETURN
      END
c
c
C       ROUTINE DRAWS ONLY OVERLAY IN dotted line
      SUBROUTINE VECOVL21(IH,fildes,x0)

      include 'x2000parameter.defs'

      INCLUDE 'gopen_type.defs'
      integer IH(101)
      real x0
      if(.false.)print *,fildes2 !suppress unused variable warning
      call clear2(fildes)
*       call line_color_index(fildes,7)
      linetype=2                !dotted line
      call dohis21(IH,fildes,linetype,7,x0)
c        call make_picture_current(fildes)
      RETURN
      END
c
c
c
C       ROUTINE TO CARRY OUT FINAL SCALING BEFORE PLOT
      SUBROUTINE DISSCZ1(ihs,IH,ISTAT,NTOP,SCLFC)
      integer ihs(101),IH(101)
      if(.false.)print *,ISTAT  !suppress unused variable warning
      NTOP=0
      DO 10 I=1,100
         NTOP=MAX0(NTOP,ihs(I))
 10   CONTINUE
      IF(NTOP.EQ.0)SCLFC=1.0
      IF(NTOP.EQ.0)GOTO 30
      SCLFC=300.0/NTOP
 30   DO 20 J=1,100
         ZVAL=ihs(J)*SCLFC
         IH(J)=IFIX(ZVAL)
         IF(IH(J).GT.300)GOTO 100
         IF(IH(J).LT.0)GOTO 100
         GOTO 20
 100     PRINT 110
 110     FORMAT(2X,'IH VALUE ERROR')
 20   CONTINUE
      IH(101)=0
      RETURN
      END
c
c
c
c

c
      SUBROUTINE clear2 (fildes)
c
c
      include 'x2000parameter.defs'


      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
      call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
      RETURN
      END
c
