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

      module mod_new_draw_button
      contains
        subroutine draw_button(fildes,x0,y0,x1,y1,label11,
     +     x_label1,y_label1,label12,x_label2,y_label2)
*
*       filename = new_draw_button.f
*
*       date of last revision = 18-mar-2005     lss
*
*       link with x200*
*
*
*
        include 'x2000parameter.defs'

        real x0,y0,x1,y1,x_label1,y_label1,x_label2,y_label2
*
        character*(*) label11,label12
*
        INCLUDE 'gopen_type.defs'
        if(.false.)print *,fildes2 !suppress unused variable warning
*
*
        call line_type(fildes,SOLID)
        call interior_style(fildes,INT_HOLLOW,1)
*
        call character_height(fildes,.036)
        call character_width(fildes,.012)
*

        if(label11.eq.'perimeter')then
          call rectangle(fildes,x0,y0,x1,y1)
          call rectangle(fildes,x0+1.,y0+1.,x1-1.,y1-1.)
          call rectangle(fildes,x0+2.,y0+2.,x1-2.,y1-2.)
*         goto 10
          return
          end if
*
        call rectangle(fildes,x0,y0,x1,y1)
*
        if(label11.ne.' ')then
         call text2d(fildes,x_label1,y_label1,label11//char(0),
     +          ANNOTATION_TEXT,0)
         end if
*
        if(label12.ne.' ')then
           call text2d(fildes,x_label2,y_label2,label12//char(0),
     +          ANNOTATION_TEXT,0)
         end if
*
*
*
        return
        end
      end module mod_new_draw_button