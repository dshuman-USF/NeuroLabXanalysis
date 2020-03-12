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


      module mod_plot_spikes_per_cycle
      contains
*     filename: plot_spikes_per_cycle.f (17-Mar-2004    lss)

*     date of last revision:  27-Oct-2004       lss

      subroutine plot_spikes_per_cycle(sp_per_cycle,num_acc_cycles,
     +     REFcode,TARcode,fildes,IDs,
     +     x0,y0,height)
      use mod_miscellaneous_subroutines

      include 'x2000parameter.defs'

      INCLUDE 'gopen_type.defs'

      integer :: sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +     REFcode,TARcode,num_acc_cycles,IDs(MAX_NUM_CODES),
     +     sp_per_cycle_REF(MAX_NUM_ACC_CYCLES),
     +     sp_per_cycle_TAR(MAX_NUM_ACC_CYCLES),
     +     temp_num_cycles = 0

      character*120 textstring
      character*130 text130
      character*4 REFcodetext,TARcodetext
      character*10 c_max_bin, c_num_acc_cycles
      if(.false.)print *,fildes2 !suppress unused variable warning

      call character_height(fildes,.015)
      call character_width(fildes,.005)

      write (REFcodetext,'(I4)') REFcode
      write (TARcodetext,'(I4)') TARcode

      do i = 1, MAX_NUM_ACC_CYCLES !initialize the temporary array
         sp_per_cycle_REF(i)=0
         sp_per_cycle_TAR(i)=0
      end do

      n = num_acc_cycles / 1000
      j = 0
      do i = 1, num_acc_cycles,n+1
         j = j+1
         do k = i,i+n
            sp_per_cycle_REF(j)= sp_per_cycle_REF(j)
     +           +sp_per_cycle(IDs(REFcode),k)
            sp_per_cycle_TAR(j)= sp_per_cycle_TAR(j)
     +           +sp_per_cycle(IDs(TARcode),k)
         end do
      end do
      temp_num_cycles = j

*       ***** draw the spikes-per-cycle histogram - show the REF and TAR cell activities: *****

      call line_color(fildes,0.,0.,0.)
      call line_type (fildes,SOLID)
c               x0 = 400.
c               y0 = 745.
c               height = 50.
      width = float(temp_num_cycles)
      call move2d (fildes,x0,(y0+height)) !Y axis
      call draw2d (fildes,x0,y0)
      call draw2d(fildes,(x0+width),y0) !X axis
      call text_color(fildes,0.,0.,1.)
      textstring = REFcodetext
      call text2d(fildes,x0-40.,y0+20.,textstring//char(0),
     +     ANNOTATION_TEXT,0)
      if(TARcode.ne.REFcode)then
         textstring = TARcodetext
         call text_color(fildes,1.,0.,0.)
         call text2d(fildes,x0-40.,y0,textstring//char(0),
     +        ANNOTATION_TEXT,0)
      end if
      call text_color(fildes,0.,0.,0.)
      textstring = '# spikes'
      call text2d(fildes,x0-80.,y0+height-10.,
     +     textstring//char(0),
     +     ANNOTATION_TEXT,0)

*     ***** need to know which spike train has the highest number of spikes in any cycle: *****

      max_bin_val_REF = 0
      do i = 1,temp_num_cycles  !find the max bin over both histograms
         max_bin_val_REF=MAX0
     +        (max_bin_val_REF,(sp_per_cycle_REF(i)))
      end do
      max_bin_val_TAR = 0
      do i = 1,temp_num_cycles  !find the max bin over both histograms
         max_bin_val_TAR=MAX0
     +        (max_bin_val_TAR,(sp_per_cycle_TAR(i)))
      end do
      max_bin_value = MAX0(max_bin_val_REF,max_bin_val_TAR)

c               max_bin_value = 0
c               do i = 1,temp_num_cycles                                        !find the max bin over both histograms
c                  max_bin_value=MAX0
c     +             (max_bin_value,(sp_per_cycle_REF(i)))
c                  max_bin_value=MAX0
c     +             (max_bin_value,(sp_per_cycle_TAR(i)))
c               end do
      if(max_bin_value.eq.0)then !avoid a "divide by zero" situation
         s_f = 1.0
      else
         s_f = height/FLOAT(max_bin_value)
      end if
      write (c_max_bin,'(I10)') max_bin_value
      write (c_num_acc_cycles,'(I10)') num_acc_cycles
      textstring = '(max_bin = '//c_max_bin/
     +     /', number_of_cycles = '//c_num_acc_cycles//')'
      call remove_all_blanks(textstring,LEN(textstring))
      text130 = 'cycle #  '//textstring
      call text2d(fildes,x0+10.,y0-15.,text130,
     +     ANNOTATION_TEXT,0)

*       *****   plot the data from the sp_per_cycle array for the REFerence cell:   *****

      x1 = x0 
      pixels_per_bin = 1.
      call line_color(fildes,0.,0.,1.) !blue line
      call move2d(fildes,x1,
     +     float(sp_per_cycle_REF(i))*s_f+y0)
      do i = 1, temp_num_cycles
         call draw2d (fildes,x1,
     +        (float(sp_per_cycle_REF(i))*s_f)+y0)
         x1 = x1 + pixels_per_bin
      end do

*       *****   plot the data from the sp_per_cycle array for the TARGET iff REF <> TAR:   *****

      if(TARcode.ne.REFcode)then
         x1 = x0      
         pixels_per_bin = 1.
         call line_color(fildes,1.,0.,0.) !red line
         call move2d(fildes,x1,
     +        float(sp_per_cycle_TAR(i))*s_f+y0)
         do i = 1, temp_num_cycles
            call draw2d (fildes,x1,
     +           (float(sp_per_cycle_TAR(i))*s_f)+y0)
            x1 = x1 + pixels_per_bin
         end do
      end if

*       ***** re-scale and re-plot the data from the less-active spike train: *****

      if(max_bin_val_REF.lt.max_bin_val_TAR)then
         if(max_bin_val_REF.eq.0)then !avoid a "divide by zero" situation
            s_f = 1.0
         else
            s_f = height/FLOAT(max_bin_val_REF)
         end if
         x1 = x0      
         pixels_per_bin = 1.
         call line_color(fildes,0.,0.,1.) !blue line
         do i = 1, temp_num_cycles
            call move2d (fildes,x1,
     +           (float(sp_per_cycle_REF(i))*s_f)+y0)
            call draw2d (fildes,x1+1.,
     +           (float(sp_per_cycle_REF(i))*s_f)+y0)
            x1 = x1 + pixels_per_bin
         end do
      else
         if(max_bin_val_TAR.eq.0)then !avoid a "divide by zero" situation
            s_f = 1.0
         else
            s_f = height/FLOAT(max_bin_val_TAR)
         end if
         x1 = x0      
         pixels_per_bin = 1.
         call line_type(fildes,iDOT)
         call line_color(fildes,1.,0.,0.) !red line
         do i = 1, temp_num_cycles
            call move2d (fildes,x1,
     +           (float(sp_per_cycle_TAR(i))*s_f)+y0)
            call draw2d (fildes,x1+1,
     +           (float(sp_per_cycle_TAR(i))*s_f)+y0)
            x1 = x1 + pixels_per_bin
         end do
      end if

      call line_color(fildes,0.,0.,0.)
      call line_type(fildes,iSOLID)

*     ***** re-plot the x and y-axes just in case they have been drawn-over in a different color: *****

      call move2d (fildes,x0,(y0+height)) !Y axis
      call draw2d (fildes,x0,y0)
      call draw2d(fildes,(x0+width),y0) !X axis

      return
      end
      end module mod_plot_spikes_per_cycle
