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


      module mod_calcCCH_newshift_2
      contains
*       
*       filename: calcCCH_newshift_2.f
*
*       date of last revision = 15-dec-2004     lss
*
*       May-2003        lss
*         this "new" subroutine calculates both the CCHs and the shift-predictor
*           controls; the REF train is shifted within calc_and_write_newshift and
*           the resulting shifted times are stored in SPIKESHIFT; 
*           the value of calc_shift tells the code which histogram to calculate.
*
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*
*
*       link with x2002newshift code
*       The following is a stab at an explanation of why I have monkeyed around with the spike times: 

*               Double precision arithmetic is a treacherous thing.  DP floats may not be represented exactly within
*               the binary depths of the computer.  So ... it behooves one to do as much integer arithmetic as
*               possible.  Since our spike times are initially expressed as milliseconds (xxx.x) in SPIKETIMES(),
*               we can get rid of the decimal point by multiplying by 10 (*10D0) and then rounding to the nearest 
*               integer (IDNINT):  
*                       TIME = IDNINT(TIME * 10D0)
*               Now each time is expressed in tenths of a millisecond, and we can do "integer" math with it.
*       
*               OK, now for the binwidth thing and how to know which bin of the histogram to increment:
*               Our binwidths are also expressed as milliseconds (xx.x).  That's not OK to use with the IDNINT'ed
*               spike times, so we have to do BINW = IDNINT(BINW*10D0).  

*               Now we can calculate the bin by: BIN# = (DELT/BINW_rnd) + 1     [BIN# is defined as LAG in this subroutine]
*
      SUBROUTINE calcCCH_newshift (calc_shift,SPIKETIMES,
     +     SPIKESHIFT,IRPT,ITPT,ITAL,
     +     BINW,IHIST,OFSET,
     +     LAST_REF,LAST_TAR,IDs,count,ishift)
*
      include 'x2000parameter.defs'
*
      double precision SPIKESHIFT(:,:)
      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN)
      integer*4 IHIST(NUM_BINS)
      integer count,low,high,probe

      INTEGER*4 LAST_REF,LAST_TAR,bin_ID,err
*
      DOUBLE PRECISION REF, TAR, first_time,
     +     max_delta,last_time,
     +     high_time,low_time

      character*1 calc_shift


      if(.false.)print *,OFSET  !suppress unused variable warning

*
*       ***** initialize the arrays and variables *****
*
      err = 0
      IHIST=0
      LAST_REF=0
      LAST_TAR=0
      max_delta = IDNINT((((NUM_BINS/2)*BINW) + BINW/2.0)*100D0)
*
*       ************************************************
*
      if(calc_shift.eq.'y')then
         num_of_REFs = count
      else
         num_of_REFs = ITAL(IDs(IRPT))
      end if

c ITB guarantee: either the spike at ITB is in range of REF, or it isn't
c and neither is anything to the left.  If this is true of ITB for one
c REF spike, it is true for all subsequent REF spikes.  This means that
c if ITB is past the range, we don't need to look to the left for a
c starting point in the range.

      ITB = 1
      do IR = 1, num_of_REFs
         if(calc_shift.eq.'y')then
            REF = IDNINT(SPIKESHIFT(IR,ishift)*10D0)*10D0
         else
            REF = IDNINT(SPIKETIMES(IR,IDs(IRPT))*10D0)*10D0
         end if
         
c For each REF there is a target range from and including first_time to
c and excluding last_time

         first_time = REF - max_delta !TAR must not be less than first_time
         last_time = REF + max_delta !TAR must not be greater than last_time

         TAR = IDNINT(SPIKETIMES(ITB,IDs(ITPT))*10D0)*10D0
         
         if (TAR.lt.first_time) then
            goto 500  !try binary search
         else if (TAR.ge.last_time) then
            cycle     !no TAR in range of this REF
         else
            goto 1000 !ITB is in range
         end if

 500     itar = ITAL(IDs(ITPT))
         TAR = IDNINT(SPIKETIMES(itar,IDs(ITPT))*10D0)*10D0
         if (TAR.lt.first_time) then
            exit       !this CCH is done
         else if (TAR.ge.last_time) then
            goto 600   !binary search
         else
            ITB = itar !itar is in range
            goto 1000
         end if

 600     high = itar !perform a BINARY SEARCH to find the location of the first target event
         low = ITB              ! that will be included in the histogram wrt this new reference event
         k=0
         do while (high-low.gt.1) !when low and high are adjacent, the search is finished, and low will
            k=k+1
            low_time=IDNINT(SPIKETIMES(low,IDs(ITPT))*10D0)*10D0
            high_time=IDNINT(SPIKETIMES(high,IDs(ITPT))*10D0)*10D0
            probe=DMAX1((((first_time-low_time)/
     +           (high_time-low_time))*
     +           (high-low) + low),low+1D0)
c first_time > low_time and high_time > low_time and high_time >
c first_time, so the fraction must be less than 1 and DMAX1 truncates,
c so probe must be less than high.
            if(probe.ge.high)goto 2001 !BUG
            TAR = IDNINT(SPIKETIMES(probe,IDs(ITPT))*10D0)*10D0
            if (TAR.lt.first_time) then
               low = probe
            else if (TAR.ge.last_time) then
               high = probe
            else
               ITB = probe !probe is in range
               goto 1000
            end if
         end do
         ITB = high
         cycle !no TAR in range of this REF

c ITB is in range at this point
 1000    do IT = ITB, 1, -1     !look backwards
            TAR = IDNINT(SPIKETIMES(IT,IDs(ITPT))*10D0)*10D0
            if(TAR.lt.first_time) exit
            bin_ID=INT((TAR-first_time)/IDNINT(BINW*100D0))+1
            if((bin_ID.lt.1).or.(bin_ID.gt.NUM_BINS)) goto 2002 !BUG
            IHIST(bin_ID) = IHIST(bin_ID) + 1
            IF(IT.gt.LAST_TAR)LAST_TAR=IT
            IF(IR.gt.LAST_REF)LAST_REF=IR
         end do
         do IT = ITB+1, ITAL(IDs(ITPT)), 1 !OK - have gone backwards far enough -- now go forwards
            TAR = IDNINT(SPIKETIMES(IT,IDs(ITPT))*10D0)*10D0
            if(TAR.ge.last_time) exit
            bin_ID=INT((TAR-first_time)/IDNINT(BINW*100D0))+1
            if((bin_ID.lt.1).or.(bin_ID.gt.NUM_BINS)) goto 2003 !BUG
            IHIST(bin_ID) = IHIST(bin_ID) + 1
            IF(IT.gt.LAST_TAR)LAST_TAR=IT
            IF(IR.gt.LAST_REF)LAST_REF=IR
         end do
         ITB = LAST_TAR
      end do 
      RETURN
 2003 err = err + 1
 2002 err = err + 1
 2001 err = err + 1
      print *,'bin_ID: ', bin_ID
      print *,'TAR: ', TAR
      print *,'first_time: ', first_time
      print *,'BINW: ', BINW
      print *,'ITB: ', ITB
      print *,'IT: ', IT
      print *,'BUG',err,' in calcCCH_newshift, aborting'
      print *,'SPIKETIMES(IT,IDs(ITPT)): ', SPIKETIMES(IT,IDs(ITPT))
      print *,'SPIKETIMES(ITB,IDs(ITPT)): ',SPIKETIMES(ITB,IDs(ITPT))
      print *, 'IDs(ITPT): ',IDs(ITPT) 
      print *,SPIKETIMES(119,IDs(ITPT)),SPIKETIMES(120,IDs(ITPT))
      stop
      END
      end module mod_calcCCH_newshift_2
