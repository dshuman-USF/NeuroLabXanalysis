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

      module mod_serotonergic_discr_func
      contains
*     filename: serotonergic_discr_func.f
*     created 04-Nov-2004       lss

*     date of last revision: 10-Dec-2003        lss

*     derived from Peggy Mason, Physiological identification of pontomedullary serotonergic neurons in the rat.  
*                       J Physiol.,1087-1098, 1997.

*     This subroutine of x2002 will calculate the sertotonergic discriminant score (y) as
*     a measure of the likelihood that a cell is serotonergic-like based upon its firing rate and anatomic location.

*               y                = 146 - meanISI + 0.98sdISI
*                (meanISI-sdISI)
*       The discriminant function defines two conditions that must be met for a cell to be classified as serotonergic:
*               1.  because sdISI must be >0, then meanISI > 146ms (ie., the discharge rate must be <7 Hz)
*               2.  the discharge must be sufficiently regular such that sdISI < meanISI by >= 146 ms.

*                so...  if y < 0 --> serotonergic-like activity
*                       if y > 0 --> not serotonergic-like


      subroutine serotonergic_discr_score(SPIKETIMES,IDs,ITAL,cell,
     +     meanISI_text,sdISI_text,y_text,MAX_INT,num_rej_ISI)
      use mod_mean_and_sd_routines

      INCLUDE 'x2000parameter.defs'

      integer cell,rejected
      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      double precision SPIKE,nextSPIKE,ISI(ITAL(IDs(cell))),
     +     MAX_INT,meanISI,sdISI,y
      character*10 y_text
      character*10 meanISI_text,sdISI_text,num_rej_ISI

      meanISI_text = ' '
      sdISI_text =  ' '
      sum = 0.0
      num_rej_ISI = ' '
      rejected = 0
      do i = 1,ITAL(IDs(cell))-1
         ISI(i) = 0.0
      end do

      j = 0
      do i = 1,ITAL(IDs(cell))-1
         SPIKE = SPIKETIMES(i,IDs(cell))
         nextSPIKE = SPIKETIMES(i+1,IDs(cell))
         if((nextSPIKE - SPIKE).le.MAX_INT)then !acceptable interval --> save it
            j = j+1
            ISI(j) = nextSPIKE - SPIKE
            sum = sum + ISI(j)
         else
            rejected = rejected + 1
         end if
      end do

      if(j.lt.2)then            !need at least 2 acceptable ISIs
         print '(''ID '',I3,'': not enough ISIs to calculate PMS!'')',
     +        cell
         meanISI_text = 'No data'
         sdISI_text = 'No data'
         y_text = 'No data'
         write (num_rej_ISI,'(I10)') rejected
         return
      end if
      
c        print '(I3,''; num_rej_ISI = '',I6)',cell,rejected

      smean= sum/j
c        print '(''smean = '',f7.2)',smean

      call dp_real_mean_and_sd(ISI,j,meanISI,sdISI)

      y = 146 - meanISI + 0.98*sdISI !discriminant score

c        print '(I3,'': meanISI = '',f9.2,''; sdISI = '',f9.2)',
c     +       cell,meanISI,sdISI
c       if(y.gt.0)print '(''serotonergic disc score: '',f7.2,
c     +       '' --> NO'')',y
c       if(y.lt.0)print '(''serotonergic disc score: '',f7.2,
c     +       '' --> YES'')',y
c       if(y.eq.0)print '(''serotonergic disc score: '',f7.2,
c     +       '' --> on the line'')',y


      write (y_text,'(f7.2)') y

      write (y_text,'(ES10.2E2)')y

      write (meanISI_text,'(f7.2)') meanISI
      write (sdISI_text,'(f7.2)') sdISI
      write (num_rej_ISI,'(I10)') rejected

      return
      end

c
c
c
c       filename: relative_ISI.f
*     created 04-Dec-2003       lss

*     date of last revision: 10-Dec-2003        lss
c       
c       from Wang, Pizzonia, and Richerson. Chemosensitivity of rat medullary raphe nerones in primary
c               tissue culture, J. Physiol., 433-450, 1998.
c
c       This code provides a measure of the variability of the interspike interval (ISI) -- the relative ISI is
c       calculated for each spike as the ratio of ISI after the spike to ISI before the spike.  The relative ISI
c       will reflect the beat-to-beat regularity of a neuron whose firing rate is slowly ramping up
c       or down (the CV of the absolute ISI for such a cell may not reveal the cell's highly regular firing pattern).

c       The user may place restrictions on relative ISI values so that 'long' gaps in cell activity will not contribute
c       to the values of the mean rISI and the sd rISI.

c       call rISI(SPIKETIMES,IDs,ITAL,min_rISI,max_rISI,mean_rISI_text(IDs(cell)),sd_rISI_text(IDs(cell)),cell)

      subroutine calc_rISI(SPIKETIMES,IDs,ITAL,
     +     mean_rISI_text,sd_rISI_text,cell,MAX_INT,
     +     num_rej_rISI)
      use mod_mean_and_sd_routines

      INCLUDE 'x2000parameter.defs'

      integer cell,rejected
      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      double precision SPIKE,preSPIKE,postSPIKE,
     +     rISI(ITAL(IDs(cell))),mean_rISI,
     +     sd_rISI, MAX_INT            
      character*10 mean_rISI_text,sd_rISI_text,num_rej_rISI

      num_rej_rISI = ' '
      rejected = 0
      j=0

      do i = 2, ITAL(IDs(cell))-1 !careful -- the rejected count could be inflated because 
         SPIKE = SPIKETIMES(i,IDs(cell)) ! if(postSPIKE-SPIKE)>MAX_INT, then it will be rejected again
         preSPIKE = SPIKETIMES(i-1,IDs(cell)) ! as (SPIKE-preSPIKE)
         if((SPIKE-preSPIKE).le.MAX_INT)then
c              print '(''interval A = '',f)',SPIKE-preSPIKE
            postSPIKE = SPIKETIMES(i+1,IDs(cell))
            if((postSPIKE-SPIKE).le.MAX_INT)then
c                 print '(''interval B = '',f)',postSPIKE-SPIKE
c                 print '(''pre: '',f,''  SPIKE: '',f,''  post: '',f)',
c     +                preSPIKE,SPIKE,postSPIKE
               j=j+1
               rISI(j) = (postSPIKE-SPIKE)/(SPIKE-preSPIKE)
               if(SPIKE.eq.preSPIKE.or.postSPIKE.eq.SPIKE)then
                  print '(//)'
                  print *,'FATAL ERROR'
                  print '(A,I3,A,F11.1,A)',' Cell ',cell,
     +                 ' has more than one spike at ',SPIKE,' ms'
                  print *,'You need to fix this with ',
     +                 'bigfix and then start all over.'
                  print '(//)'
                  stop
c     read (*,'(A)')
               end if
            else
               rejected = rejected + 1
            end if
         else
            rejected = rejected + 1
         end if
      end do

      if(j.lt.2)then            !need at least 2 acceptable ISIs
         mean_rISI_text = 'No data'
         sd_rISI_text = 'No data'
         write (num_rej_rISI,'(I10)') rejected
         return
      end if

c        print '(I3,''; num_rej_rISI = '',I6)',cell,rejected


      call dp_real_mean_and_sd(rISI,j,mean_rISI,sd_rISI)
      write (num_rej_rISI,'(I10)') rejected
      write (mean_rISI_text,'(f7.2)') mean_rISI
      write (sd_rISI_text,'(f7.2)') sd_rISI

      return
      end


      end module mod_serotonergic_discr_func
