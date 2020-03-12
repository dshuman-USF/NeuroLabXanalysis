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


      module mod_cardsig_delta2
      contains
*       filename = cardsig_delta2f.
*
*       date of latest revision = 24-Nov-2003   lss

*       delta-squared value derived as indication of a cell's cardiac modulation

C       derived from TEDS3.F  24-Nov-2003         lss
*
*       mar-2000        lss
*        modified to allow ID codes to range from 1 to 999;     
*        maximum number of units allowed still = 120;           
*        *.defs files inserted                                  
*        (MAX_NUM_CODES=999; MAX_NUM_CHAN=120)                  
*       *** INDIRECT POINTERS are now used to access unit data
*       ***  within the following arrays:
*       ***     SPIKETIMES, ITAL, card_type, coef, coefnum, TALLY_NORM,
*       ***     ETA2_*, zmodsig_*, zmodsig2_*, analyzed_cells, analyzed_pairs,
*       ***     resp_type, card, CELL_NAMES
*       ***     [these arrays are now dimensioned to MAX_NUM_CHAN;
*       ***      the location of a cell's information (user-assigned ID code = i)
*       ***      within these arrays is array(IDs(i))]
*       ***
*       *** DIRECT POINTERS are used with these arrays: exclude, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       link with x2000_v2 code
*
*
*
*       
*               increased max # spike trains and # events
*
*      IMAXCYC is in msec but integer 
*      subroutine FOR DETERMINATION OF A CELL'S cardiac
*      MODULATION, WILL RETURN PROBABILITY  
*      OF THE NULL HYPOTHESIS (THAT THE CELL IS NOT CARDIAC MODULATED)
*
*
      subroutine cardsig_delta2(SPIKETIMES,ICELPT,ICYPT,IMAXCYC,
     +     card_type,IDs,tedfactor,DELTA2)
      use mod_betai


      include 'x2000parameter.defs'

      double precision SPIKETIMES(:,:)
*
      integer IDs(MAX_NUM_CODES)
      PARAMETER (SIGLVL=.05,ZSIGLVL=1.645,itrt=5,
     +     isubmax=51)

      DIMENSION IY(itrt,isubmax), JSUM(itrt)
      REAL*4 :: PROB = 0.0
      DOUBLE PRECISION CYST,CYND,CLLK
      character*5 DELTA2(MAX_NUM_CHAN)
      CHARACTER*3 card_type(MAX_NUM_CHAN)
      character*2 tedfactor(MAX_NUM_CHAN)
      real :: zero = 0
*       real*16 bicokmq
*      CLEAR ALL ARRAYS,INITIALIZE VARIABLES,ETC HERE 
*
*
      max_num_events = size (spiketimes, 1)
c      print '(''max_num_events='',I7)', max_num_events
      DELTA2temp = 0/zero
      itedfactor = 1
      print '(''cell ID for cardsig: '',I3)',ICELPT
 1    iskipcnt=0
      isub=50
      RMSB=0
      RMSW=0
      itedcounter = 0
      DO 184, J=1,itrt
         JSUM(J)=0
         DO 182, I=1,isub
            IY(J,I)=0
 182     CONTINUE
 184  CONTINUE
      IYSUM=0
      I=1
      J=1
      ICELPT2=1
      ICYPT2=1
 188  CYST=SPIKETIMES(ICYPT2,IDs(ICYPT))
      CYND=SPIKETIMES(ICYPT2+1,IDs(ICYPT))
      print '(''CYST = '',F20.5)',CYST
      print '(''CYND = '',F20.5)',CYND
      print '(''CYND-CYST = '',F20.5)',(CYND-CYST)
      print '(''imaxcyc = '',I10)',IMAXCYC
      IF (CYND-CYST.GT.IMAXCYC) THEN
         ICYPT2=ICYPT2+1
         IF ((ICYPT2.GE.MAX_NUM_EVENTS).OR.(CYND.EQ.0)) THEN
            GOTO 300
         ENDIF
         GOTO 188
      ENDIF
C      find the first cell event after the first cycle start
 190  CLLK=SPIKETIMES(ICELPT2,IDs(ICELPT))
      print '(''CLLK = '',f20.5)',CLLK
      IF (CLLK.LT.CYST) THEN
         ICELPT2=ICELPT2+1
         if(ICELPT2.ge.max_num_events) print '(''ICELPT2='',I7)', ICELPT2
         if(ICELPT2.eq.max_num_events) print '(''Warning: There may be ''
     +        ''an issue with cell '',I4,'' - too few spikes, etc.'')',
     +        ICELPT
         if(ICELPT2.eq.max_num_events) goto 300
         GOTO 190
      ENDIF
C      count the cell events in each 1/itrt th of cycle
 200  IF (CLLK.LT.(CYST+J*(CYND-CYST)/itrt)) THEN
         IY(J,I)=IY(J,I)+1      !counts events in array subscripted
C                                                by 1/itrt ths and cycle #
         IYSUM=IYSUM+1          !sum of events in all cycles
         JSUM(J)=JSUM(J)+1      !sum of events is each 1/trt th
         iskipcnt=iskipcnt+1    !sum of events in cycle
         ICELPT2=ICELPT2+1
         CLLK=SPIKETIMES(ICELPT2,IDs(ICELPT))
         IF ((ICELPT2.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
            GOTO 300
         ENDIF
         GOTO 200
      ELSE                      !it's the next 1/itrt th
         J=J+1
         IF (J.LE.itrt) GOTO 200 !is it time for the next cycle?
      ENDIF
 205  ICYPT2=ICYPT2+1
      IF ((ICYPT2.GE.MAX_NUM_EVENTS).OR.(CYND.EQ.0)) THEN
         GOTO 300
      ENDIF
      CYST=SPIKETIMES(ICYPT2,IDs(ICYPT))
      CYND=SPIKETIMES(ICYPT2+1,IDs(ICYPT))
      IF (CYND-CYST.GT.IMAXCYC) GOTO 205 
 207  IF (CLLK.LT.CYST) THEN
         ICELPT2=ICELPT2+1
         CLLK=SPIKETIMES(ICELPT2,IDs(ICELPT))
         IF ((ICELPT2.GE.MAX_NUM_EVENTS).OR.(CLLK.EQ.0)) THEN
            GOTO 300
         ENDIF
         GOTO 207
      ENDIF
      J=1
C       we only count cycles with events now
      IF (ISKIPCNT.NE.0) THEN
         I=I+1
         ISKIPCNT=0
      ENDIF
      IF (I.LE.isub) GOTO 200   !for the next cycle
      itedcounter = itedcounter + 1
      if(itedcounter.lt.itedfactor)then
         I = 1
         J = 1
         goto 200
      end if
      AVGY=IYSUM/float(isub*itrt) !global average
      DO 210, J=1,itrt
         RMSB=RMSB+(float(JSUM(J))/float(isub)-AVGY)**2
 210  CONTINUE
      RMSB=isub*RMSB/(itrt-1)
      DO 230, I=1,isub
         DO 220, J=1,itrt
            RMSW=RMSW+(float(IY(J,I))-
     +           float(JSUM(J))/float(isub))**2
 220     CONTINUE
 230  CONTINUE
      RMSW=RMSW/(itrt*(isub-1))
      IF (RMSW.NE.0) THEN
         F=RMSB/RMSW
      ELSE
         IF(RMSB.eq.0)then
            F=1.0
         ELSE
            F=9999.0
         ENDIF
      ENDIF

      DELTA2temp = (itrt-1)*F/((itrt-1)*F+(itrt*(isub-1)))
      PROB = BETAI(0.5*(itrt*(isub-1)),
     +     0.5*(itrt-1),(itrt*(isub-1))/
     +     ((itrt*(isub-1))+(itrt-1)*F))

      if(itedfactor.eq.1)then
         itedfactor = 10
         goto 1
      end if
      if(itedfactor.eq.10)then
         itedfactor = 20
         goto 1
      end if
      if(itedfactor.eq.20)then
         itedfactor = 50
         goto 1
      end if
*
 300  if(itedfactor.ge.20)then
         if(PROB.lt.SIGLVL)card_type(IDs(ICELPT))='C'
         if(PROB.ge.SIGLVL)card_type(IDs(ICELPT))='N'
      else                      !tedfactor = 1 or 10 --> not enough data for a clear call on modulation
         card_type(IDs(ICELPT))='ned'
      end if

      write (DELTA2(IDs(ICELPT)),'(f5.2)') DELTA2temp
      write (tedfactor(IDs(ICELPT)),'(I2)') itedfactor
c       print '(''ID '',I3,'',: DELTA2 = '',A5,'' ('',f5.2,
c     +            ''); tedfactor = '',I2,
c     +            ''('',A2,
c     +            '');  card_type = '',A3)',
c     +      ICELPT,DELTA2(IDs(ICELPT)),DELTA2temp,itedfactor,
c     +      tedfactor(IDs(ICELPT)),card_type(IDs(ICELPT))

      
      return
      END 

      end module mod_cardsig_delta2
