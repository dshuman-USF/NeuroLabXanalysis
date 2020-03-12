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

      module mod_qdttxt_write
      contains

*       filename = qdttxt_write.f
*
*       date of last revision = 19-APR-2006             lss
*
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
*       *** DIRECT POINTERS are used with these arrays: excluded, included
*       ***     (these arrays are dimensioned to MAX_NUM_CODES)
*
*       This subroutine of x2000 produces text files containing database information
*               that will eventually be imported into the CELLS and PAIRS databases maintained
*               with Microsoft Access.  Specifications for these text files is as follows:
*
*
*               *****************************************************************
*               *                                                               *
*               * FORMAT OF OUTPUT ASCII FILE THAT WILL BE LOADED INTO THE      *
*               *   ACCESS DATABASE 'CELLS':                                    *
*               *                                                               *
*               * -- each field in the output file corresponds to a data        *
*               *    item for a cell (eg., name, type, coordinates) or for      *
*               *    a cell pair (eg., primary feature of histogram)            *
*               *                                                               *
*               * -- each record in the output file contains all the data       *
*               *    for one cell (each record contains ____ fields) or for     *
*               *    one cell pair (each record contains ____ fields)           *
*               *                                                               *
*               * -- records must be separated by <CR><LF>                      *
*               *                                                               *
*               * -- the number of records in an output file equals the         *
*               *    number of simultaneously recorded cells in that group      *
*               *    plus the number of pairs derived from those cells          *
*               *                                                               *
*               * -- text strings containing blanks must be enclosed in         *
*               *    quotation (") marks                                        *
*               *                                                               *
*               * -- all fields may be enclosed in quotation marks              *
*               *                                                               *
*               * -- fields must be separated from one another by either a      *
*               *    comma or a blank space                                     *
*               *                                                               *
*               * -- the output filename must be lower case letters             *
*               *                                                               *
*               *****************************************************************
*
*
        subroutine write_QDTTXT(QDT_FILENAME,shortest,
     +     longest)
      use mod_miscellaneous_subroutines
      use mod_read_and_write_DBSAV


        INCLUDE 'x2000parameter.defs'

        integer BNDRY,total_histograms,
     +     I_pulse,E_pulse,phrenic,cardiac_pls,
     +     included(MAX_NUM_CODES),excluded(MAX_NUM_CODES),
     +     IDs(MAX_NUM_CODES),
     +     sp_per_cycle(MAX_NUM_CHAN,MAX_NUM_ACC_CYCLES),
     +       TOTAL_NUM_SHIFTS,total_num_cells,ITAL(MAX_NUM_CHAN)

        real STIM_OFFSET,BINW,BINW_1,BINW_2,BINW_3,BINW_4,BINW2,
     +       start_time,end_time,coefnum(MAX_NUM_CHAN),
     +       NORM_BW,NORM_OFFSET,shortest,longest
     
        character*(*) QDT_FILENAME
        character*200 protocol,BDT_FILE
        character*175 exp_name
        character*11 date_exp
        character*10 meanISI(MAX_NUM_CHAN),sdISI(MAX_NUM_CHAN),
     +               fiveHT(MAX_NUM_CHAN),mean_rISI(MAX_NUM_CHAN),
     +               sd_rISI(MAX_NUM_CHAN),num_rej_ISI(MAX_NUM_CHAN),
     +               num_rej_rISI(MAX_NUM_CHAN),c_MAX_INT
        character*9 TODAY
        character*2 recording,tedfactor(MAX_NUM_CHAN),version
        character*60 QDTTXT,QDTSAV
        character*1 BOUNDARY
        CHARACTER*5 coef(MAX_NUM_CHAN),ETA2_1(MAX_NUM_CHAN),
     +   ETA2_2(MAX_NUM_CHAN),ETA2_3(MAX_NUM_CHAN),
     +   ETA2_4(MAX_NUM_CHAN),ETA2_5(MAX_NUM_CHAN),
     +   ETA2_6(MAX_NUM_CHAN),DELTA2(MAX_NUM_CHAN)
        character*3 zmodsig_1(MAX_NUM_CHAN),zmodsig_2(MAX_NUM_CHAN),
     +     zmodsig_3(MAX_NUM_CHAN),zmodsig_4(MAX_NUM_CHAN),
     +     zmodsig_5(MAX_NUM_CHAN),zmodsig_6(MAX_NUM_CHAN),
     +     zmodsig2_1(MAX_NUM_CHAN),zmodsig2_2(MAX_NUM_CHAN),
     +     zmodsig2_3(MAX_NUM_CHAN),zmodsig2_4(MAX_NUM_CHAN),
     +     zmodsig2_5(MAX_NUM_CHAN),zmodsig2_6(MAX_NUM_CHAN),
     +     card_type(MAX_NUM_CHAN)
      character*20 cntrl

      zmodsig_1 = ' '
      zmodsig_2 = ' '
      zmodsig_3 = ' '
      zmodsig_4 = ' '
      zmodsig_5 = ' '
      zmodsig_6 = ' '
      zmodsig2_1 = ' '
      zmodsig2_2 = ' '
      zmodsig2_3 = ' '
      zmodsig2_4 = ' '
      zmodsig2_5 = ' '
      zmodsig2_6 = ' '
      card_type = ' '

        call DATE(TODAY)

        call remove_all_blanks(QDT_FILENAME,LEN(QDT_FILENAME))
        call strlength(QDT_FILENAME,LEN(QDT_FILENAME),l_qdt)

        QDTSAV = QDT_FILENAME(1:l_qdt)//'.sav'
        QDTTXT =  QDT_FILENAME(1:l_qdt)//'.txt'

        call remove_all_blanks(QDTSAV,LEN(QDTSAV))
        call remove_all_blanks(QDTTXT,LEN(QDTTXT))
        call strlength(QDTTXT,LEN(QDTTXT),l_txt)
        call strlength(QDTSAV,LEN(QDTSAV),l_sav)

*
*       ***** read the QDT.SAV file: ****
*
         call read_QDTSAV(version,QDTSAV,ios,
     +     date_exp,recording,protocol,BDT_FILE,QDT_FILENAME,
     +     IDs,excluded,included,total_num_cells,
     +     I_pulse,E_pulse,phrenic,BNDRY,cardiac_pls,
     +     BINW,NORM_BW,STIM_OFFSET,NORM_OFFSET,BINW2,
     +     BINW_1,BINW_2,BINW_3,BINW_4,BOUNDARY,start_time,
     +     end_time,icycles,
     +     total_histograms,ETA2_1,ETA2_2,ETA2_3,ETA2_4,
     +     ETA2_5,ETA2_6,zmodsig_1,zmodsig_2,zmodsig_3,
     +     zmodsig_4,zmodsig_5,zmodsig_6,zmodsig2_1,
     +     zmodsig2_2,zmodsig2_3,zmodsig2_4,zmodsig2_5,
     +     zmodsig2_6,coef,coefnum,card_type,exp_name,
     +     DELTA2,tedfactor,meanISI,sdISI,fiveHT,mean_rISI,sd_rISI,
     +     num_rej_ISI,num_rej_rISI,c_MAX_INT,TOTAL_NUM_SHIFTS,
     +     num_acc_cycles,sp_per_cycle,ITAL)
         
        call strlength(BDT_FILE,LEN(BDT_FILE),l_bdt)
        call strlength(protocol,LEN(protocol),l_protocol)
         
        do i = 1, MAX_NUM_CHAN                                         !force all these to UPPER CASE LETTERS
         call upper_case(zmodsig_1(i),LEN(zmodsig_1(i)))
         call upper_case(zmodsig_2(i),LEN(zmodsig_2(i)))
         call upper_case(zmodsig_3(i),LEN(zmodsig_3(i)))
         call upper_case(zmodsig_4(i),LEN(zmodsig_4(i)))
         call upper_case(zmodsig_5(i),LEN(zmodsig_5(i)))
         call upper_case(zmodsig_6(i),LEN(zmodsig_6(i)))
         call upper_case(zmodsig2_1(i),LEN(zmodsig2_1(i)))
         call upper_case(zmodsig2_2(i),LEN(zmodsig2_2(i)))
         call upper_case(zmodsig2_3(i),LEN(zmodsig2_3(i)))
         call upper_case(zmodsig2_4(i),LEN(zmodsig2_4(i)))
         call upper_case(zmodsig2_5(i),LEN(zmodsig2_5(i)))
         call upper_case(zmodsig2_6(i),LEN(zmodsig2_6(i)))
         call upper_case(card_type(i),LEN(card_type(i)))
        end do

        OPEN (UNIT=400,FILE=QDTTXT,FORM='FORMATTED',
     +       ACCESS='SEQUENTIAL')
        rewind (400)

        write (400,'(T2,''This file ('',A,'') was created on '',A,
     +       '' and describes parameters '',/,
     +    T2,''used to create data ''
     +       ''contained within '',A,'' and '',A,''.'',//,
     +    T2,''Date of experiment: '',A,/,T2,''Recording: '',A,/,
     +    T2,''Experimental protocol: '',A,//,
     +    T2,''Input data file = '',A)')
     +       QDTTXT(1:l_txt),TODAY,QDT_FILENAME(1:l_qdt),
     +       QDTSAV(1:l_sav),
     +       date_exp,recording,
     +       protocol(1:l_protocol),BDT_FILE(1:l_bdt)


        write (400,'(/,T2,''The following unit IDcodes were included ''
     +                  ''in histogram calculation:'')')
        j = 0
        do k = 1, MAX_NUM_CODES
            if((IDs(k).ne.0).and.(included(k).ne.0).and.
     +          (excluded(k).eq.0))then
               write (400,'(I3,''  '',$)') k                        !print the codes used to generate histograms
               j = j+1
            end if
            if(j.eq.20)then                             
               write (400,'(/,$)')                       !start a new line
               j = 0
            end if
         end do

         write (400,'(//,T2,''Total number of included cells = '',
     +        I12)')
     +        total_num_cells
         write (400,'(/,T2,''Total number of histograms in '',A,
     +        '' = '',I12)')QDT_FILENAME(1:l_qdt),total_histograms
         write (400,'(/,T2,''The following parameters were used to ''
     +        ''calculate the respiratory CTHs:'',/,
     +     T5,''I pulse = IDcode '',I3,/,
     +     T5,''E pulse = IDcode '',I3,'' (CTH trigger)''/,
     +     T5,''phrenic = IDcode '',I3,/,
     +     T5,''resp CTHs: binwidth = '',f7.1,
     +        '' msec.; offset = '',f9.1,'' msec.'',/,
     +     T5,''normalized resp CTHs: binwidth = '',f7.1,
     +        '' msec.; offset = '',f9.1,'' msec.'')')
     +        I_pulse,E_pulse,phrenic,
     +        BINW,STIM_OFFSET,NORM_BW,NORM_OFFSET
         write (400,'(T5,''cycle set (or a subset thereof):'')')
         if(BOUNDARY.eq.'A')write(400,'(T7,''all cycles'')')
         if(BOUNDARY.eq.'F')write(400,'(T7,''the first 10000 cycles'')')
         if(BOUNDARY.eq.'B')write(400,'(T7,''cycles contained ''
     +          ''within the boundaried control period'',/,
     +      T15,''(boundary ID code = '',I4,'')'')')BNDRY
         if(BOUNDARY.eq.'C')write(400,'(T7,''the first '',I3,
     +          '' cycles'')')icycles
         if(BOUNDARY.eq.'T')write(400,'(T7,''cycles occurring ''
     +          ''between '',f5.1,'' and '',f5.1,'' minutes ''
     +          ''into the recording'')')start_time,end_time
         if(shortest.ne.0.0)
     +        write(400,'(T7,''shortest cycle is '',f9.1,'' ms.'')')
     +        shortest
         if(longest.ne.0.0)
     +        write(400,'(T7,''longest cycle is '',f9.1,'' ms.'')')
     +        longest

         if(cardiac_pls.gt.0)then
            write (400,'(/,T2,''The following parameters were used to ''
     +        ''calculate the cardiac CCHs:'',/,
     +     T5,''cardiac pulse = IDcode '',I3,/,
     +     T5,''binwidth = '',f7.1,'' msec.'')')cardiac_pls,BINW2
            write (400,'(T5, ''Cardiac CTHs were generated using ''
     +                      ''cardiac pulses which occurred '',/,
     +                   T5,''within selected respiratory cycles.'')')
         else
            write (400,'(/,T2,''No cardiac CCHs were computed.'')')
         end if

         write (400,'(/,T2,''Binwidths used to calculate CCHs ''
     +                   ''and ACHs = '',
     +        3(f7.1,'', ''),f7.1,'' msec.'')')BINW_1,BINW_2,BINW_3,
     +        BINW_4

      call getenv ("XANALYSIS_CONTROL",cntrl)
      IF((trim(cntrl).ne.'').and.(cntrl.NE.'cth_cch'))then
            write (400,'(/,T2,''A total of '',I2,'' shift-control ''
     +''correlograms are generated for each CCH.'',/,
     +T5,''A shift-control correlogram is calculated by ''
     +''moving each REF event '',/,
     +T5,''to the corresponding location in the next respiratory ''
     +''phase before matching it with the TAR train.'',/)')
     +           TOTAL_NUM_SHIFTS
         else
            write (400,'(/,T2
     +,''The control correlograms are calculated based on CTHs of the re
     +ference and target cells.'',/)')
         end if

         write (400,'(T2,''Respiratory significance (expressed by the''
     +       '' ETA**2 value and the ANOVA and BINARY tests)'',/,
     +    T3,''is evaluated using several different sets ''
     +       ''of respiratory cycles: '',/,
     +    T5,''1. using the entire data set (case 1)'',/,
     +    T8,''a. the user may, if desired, eliminate cycles ''
     +          ''whose duration is'',/,
     +    T11,  ''greater than a user-defined value (=IMAXCYC) ''
     +          ''(case 2)'',/,
     +    T5,''2. using events which occur only within the ''
     +          ''selected control period'',/,
     +    T8,   ''(NOTE: if a control period was not designated ''
     +          ''during histogram calculation,'',/,
     +    T14,  ''these next three values will not have been ''
     +          ''calculated.)'',/,
     +    T8,''a. use all events (mean_control_cycle = mean ''
     +          ''control cycle duration) (case 3)'',/,
     +    T8,''b. cycles whose duration is greater than the ''
     +          ''mean + 2 SD of '',/,
     +    T11,  ''all control cycles (=control_max_value) ''
     +          ''(case 4)'',/,
     +    T8,''c. the user may eliminate cycles whose duration ''
     +          ''is greater than'',/,
     +    T11,  ''a user-defined value (=user_max_value) ''
     +          ''(case 5)'',/,
     +    T5,''3. using only those control cycles actually used to ''
     +          ''compute CTHs (case 6)'',//,
     +    T2,''The ETA**2 value and ANOVA test are calculated using ''
     +       ''the first 50 acceptable cycles '',/,
     +    T3, ''of the control period (an acceptable cycle is one ''
     +        ''during which the cell is active).'',/,
     +    T3, ''If 50 acceptable cycles are not available, the user ''
     +        ''will be informed that the '',/,
     +    T3, ''test has been performed using less than the ''
     +        ''desired number of respiratory cycles.'',//, 
     +    T2,''The BINARY test splits the respiratory cycle into 2 ''
     +       ''parts according to the activity '',/,  
     +    T3, ''pattern of the cell so that the difference between ''
     +        ''the two parts is maximal.'',/,
     +    T3, ''The BINARY test is most robust if 50 acceptable ''
     +        ''cycles are used to form the hypothesis'',/,
     +    T3, ''of a respiratory pattern, with the remainder of ''
     +        ''the cycles (preferably at least 50) used to'',/, 
     +    T3, ''confirm or deny that hypothesis.  If at least 100 '' 
     +        ''acceptable cycles are not available,'',/,
     +    T3, ''the BINARY test will be performed using half of the ''
     +        ''available cycles to form the hypothesis,'',/,
     +    T3, ''and the remaining half to test the hypothesis.  ''
     +        ''In this case, the user will be warned'',/,
     +    T3, ''that 100 cycles were not available for a robust ''
     +        ''BINARY test.'',//,
     +    T2, ''Results of significance testing for respiratory ''
     +        ''modulation are coded as follows:'',/,
     +    T5, ''r   -- respiratory-modulated'',/,
     +    T5, ''n   -- not respiratory-modulated'',/,
     +    T5, ''ned -- not enough data to confidently say the cell ''
     +               ''is not respiratory-modulated'')')

         if(cardiac_pls.gt.0)then
            write (400,'(/,T2,''The presence or absence of cardiac ''
     +        ''modulation was evaluated using the DELTA**2 statistic ''
     +        ''(Morris & Dick).'',/,
     +    T2, ''Results of significance testing for cardiac ''
     +        ''modulation are coded as follows:'',/,
     +    T5, ''  c -- cardiac-modulated'',/,
     +    T5, ''  n -- not cardiac-modulated'',/,
     +    T5, ''ned -- not enough data to confidently say whether or ''
     +           ''not the cell''''s activity is cardiac-modulated'')')
         end if

         write(400,'(/,T2,''AA responses are imported from a *.info''
     +            '' file and are coded as follows:'',
     +       /,T5,''aa  -- antidromically activated'',
     +       /,T5,''naa -- not antidromically activated'',
     +       /,T5,''nt  -- not tested from this site'')')
    
         write(400,'(/,T2,''STA responses are imported from a *.info''
     +            '' file and are coded as follows:'',
     +       /,T5,''0 -- nothing'',
     +       /,T5,''1 -- peak right - sharp'',
     +       /,T5,''   M -- peak right - sharp --> ''
     +                     ''suggest motoneuron'',
     +       /,T5,''   P -- peak right - sharp --> ''
     +                     ''suggest premotoneuron'',
     +       /,T5,''2 -- peak right - broad'',
     +       /,T5,''3 -- peak center'',
     +       /,T5,''   D -- peak center, followed by dip'',
     +       /,T5,''L -- peak left'',
     +       /,T5,''4 -- trough right'',
     +       /,T5,''5 -- trough center'',
     +       /,T5,''T -- trough left'',
     +       /,T5,''6 -- multiple peaks and troughs'',
     +       /,T5,''7 -- other'',
     +       /,T5,''8 -- not tested using this site'',
     +       /,T5,''9 -- cardiac EKG'')')

         write (400,'(/,T2,''Pertubation responses are imported from''
     +                    '' a *.per.db file and are coded as ''               
     +                    ''follows:''
     +       /,T5,''NC -- no change in activity in response to ''
     +                  ''perturbation''
     +       /,T5,''INC -- increase in activity''
     +       /,T5,''DEC -- decrease in activity''
     +       /,T5,''INC/DEC -- increase then decrease in activity''
     +       /,T5,''DEC/INC -- decrease then increase in activity''
     +       /,T5,''Complex -- complex changes in activity''
     +       /,T5,''No data -- cell is not active during perturbation ''
     +                       ''run or is not active enough to ''
     +                       ''make a call'')')

         close (unit=400)
         return
         end
      end module mod_qdttxt_write
