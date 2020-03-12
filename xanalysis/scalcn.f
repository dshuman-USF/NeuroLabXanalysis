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

      module mod_scalcn
      contains
      SUBROUTINE SCALCN (begin_E,end_E,middle_I,ebins,ibins,IDs,
     +     ITAL,NORM_BW,NORM_OFFSET,SPIKETIMES,TOTAL_SELECTED_CYCLES,
     +     cell,JHIST,spike_count)

      include 'x2000parameter.defs'

      integer ebins,ibins,IDs(MAX_NUM_CODES),ITAL(MAX_NUM_CHAN),
     +     TOTAL_SELECTED_CYCLES,JHIST(NUM_BINS),spike_count,cell
      real NORM_BW,NORM_OFFSET
      double precision
     +     begin_E(:),
     +     end_E(:),
     +     middle_I(:),
     +     SPIKETIMES(:,:)
      
      integer i,j,bin
      double precision edur,idur,t,ebindur,ibindur,ehist(ebins),
     +     ihist(ibins)

      ebindur = ebins * NORM_BW
      ibindur = ibins * NORM_BW
      ehist = 0
      ihist = 0
      JHIST = 0
      spike_count = 0
      i = 1
      edur = middle_I(i) - begin_E(i)
      idur = end_E(i) - middle_I(i)
      do j = 1, ITAL(IDs(cell))
         t = SPIKETIMES(j,IDs(cell)) 
         do while (i.le.TOTAL_SELECTED_CYCLES.and.t.ge.end_E(i))
            i = i + 1
            if (i.gt.TOTAL_SELECTED_CYCLES) exit
            edur = middle_I(i) - begin_E(i)
            idur = end_E(i) - middle_I(i)
         end do
         if (i.gt.TOTAL_SELECTED_CYCLES) exit
         if (t.lt.begin_E(i)) cycle
         if (t.lt.middle_I(i)) then
            bin = INT ((t - begin_E(i)) / edur * ebins) + 1
            ehist(bin) = ehist(bin) + ebindur / edur
         else
            bin = INT ((t - middle_I(i)) / idur * ibins) + 1
            ihist(bin) = ihist(bin) + ibindur / idur
         end if
         spike_count = spike_count + 1
      end do
      
      j = 1 + NINT (NORM_OFFSET / NORM_BW)
      do while (j.lt.1)
         j = j + ebins + ibins
      end do
      
      do i=1,NUM_BINS-1
         if (j.gt.ebins+ibins) j = 1
         if (j.le.ebins) then
            JHIST(i) = NINT (ehist(j))
         else
            JHIST(i) = NINT (ihist(j - ebins))
         end if
         j = j + 1
      end do
      JHIST(101) = ebins
      end subroutine
      end module mod_scalcn
