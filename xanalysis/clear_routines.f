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

      module mod_clear_routines
      contains
*
*       filename = clear_routines.f
*
*       date of last revision:  08-Aug-2003     lss
*
*       This subroutine of x200* contains several subroutines for clearing the display window.
*
*       * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
*
*
*
*
*
*       ************************************************************
*       ************************************************************
*
*
      SUBROUTINE clear (fildes)
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
      call clear_control(fildes,CLEAR_DISPLAY_SURFACE)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
*
      RETURN
      END
*
*
*
*
*       *****************************************************
*       ***** clear the desired section of the screen:  *****
*       *****************************************************
*       
      subroutine clear1(fildes,flag)
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      integer*4 flag
      if(.false.)print *,fildes2 !suppress unused variable warning
*
*
      if(flag.eq.0)call view_port(fildes,0.,0.,1.,1.) !entire window
      if(flag.eq.1)call view_port(fildes,0.,.21,1.,1.) !top of window
      if(flag.eq.2)call view_port(fildes,0.,0.,1.,.20) !bottom of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      if(flag.ne.0)call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
*
*
*
*       *****************************************************
*       ***** clear the top section of the screen:      *****
*       *****************************************************
*       
      subroutine clear_top(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
*
      call view_port(fildes,0.,.21,1.,1.) !top of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
*
*
*
*       **************************************************************************************
*       ***** clear the bottom section of the screen for selecting resp cycles for CTHs: *****
*       **************************************************************************************
*       
      subroutine clear_bottom_CTH(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
*
      call view_port(fildes,0.,0.,1.,.25) !clear bottom of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end


*       *****************************************************
*       ***** clear the bottom section of the screen:   *****
*       *****************************************************
*       
      subroutine clear_bottom(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
*
      call view_port(fildes,0.,0.,1.,.21) !clear bottom of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
*
*
*
*       *****************************************************
*       ***** clear the whole window:                   *****
*       *****************************************************
*       
      subroutine clear_all(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
*
      call view_port(fildes,0.,0.,1.,1.) !entire window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
*
*
      return
      end
*
*
*
*       ***************************************************     
*       ***** clear the right quarter of the bottom:  *****
*       ***************************************************
*
      subroutine clear_quarter(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
      call view_port(fildes,.45,0.,1.,.25) !clear bottom of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
*
*
*
*       ***************************************************     
*       ***** clear the right half of the screen:  *****
*       ***************************************************
*
      subroutine clear_right_half(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
      call view_port(fildes,.50,0.,1.,1.) !clear right half of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
*
*
*
*       ***************************************************     
*       ***** clear the top slot of the CTH screen:  *****
*       ***************************************************
*
      subroutine clear_top_slot(fildes)
*
*
      include 'x2000parameter.defs'
      INCLUDE 'gopen_type.defs'
      if(.false.)print *,fildes2 !suppress unused variable warning
*
      call view_port(fildes,.176,.9,.818,1.) !clear right half of window
*
      call clear_control(fildes,CLEAR_VIEWPORT)
      call clear_view_surface(fildes)
      call make_picture_current(fildes)
*
      call view_port(fildes,0.,0.,1.,1.) !reset port to entire window
*
*
      return
      end
      end module mod_clear_routines
