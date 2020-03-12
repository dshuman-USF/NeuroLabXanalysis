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


      module mod_picktypetest2
      contains
*
*
*       subroutine picktypetest1
*
*       filename: picktypetest2.f
*
*       date of latest revision = 30-Mar-2006    lss
*
*
*       This subroutine of x200_v2 is used to select designation of firing
*               pattern for respiratory cells.  The user may elect to use a "tree"
*               method, or may select unit type from a list.
*
*               4 -- I-AUG-P    17 -- E-AUG-P    28 -- I-I/E-AUG-P    49 -- I-E/I-AUG-P
*               5 -- I-AUG-T    18 -- E-AUG-T    29 -- I-I/E-AUG-T    50 -- I-E/I-AUG-T
*               6 -- I-DEC-P    19 -- E-DEC-P    30 -- I-I/E-DEC-P    51 -- I-E/I-DEC-P,
*               7 -- I-DEC-T    20 -- E-DEC-T    31 -- I-I/E-DEC-T    52 -- I-E/I-DEC-T
*               8 -- I-PLAT-P   21 -- E-PLAT-P   32 -- I-I/E-OTHER-P  53 -- I-E/I-OTHER-P
*               9 -- I-PLAT-T   22 -- E-PLAT-T   33 -- I-I/E-OTHER-T  54 -- I-E/I-OTHER-T
*               10 -- I-LR-P    23 -- E-LR-P     34 -- E-I/E-AUG-P    55 -- E-E/I-AUG-P
*               11 -- I-LR-T    24 -- E-LR-T     35 -- E-I/E-AUG-T    56 -- E-E/I-AUG-T
*               12 -- I-OTHER-P 25 -- E-OTHER-P  36 -- E-I/E-DEC-P    57 -- E-E/I-DEC-P
*               13 -- I-OTHER-T 26 -- E-OTHER-T  37 -- E-I/E-DEC-T    58 -- E-E/I-DEC-T
*               14 -- I-RECRUIT 27 -- E-RECRUIT  38 -- E-I/E-OTHER-P  59 -- E-E/I-OTHER-P
*               15 -- I-PHR-P                    39 -- E-I/E-OTHER-T  60 -- E-E/I-OTHER-T
*               16 -- I-PHR-T                    40 -- I/E-TRANS-P    61 -- E/I-TRANS-P
*                                                41 -- I/E-TRANS-T    62 -- E/I-TRANS-T
*               70 -- I-LATE-P                   42 -- I/E-PLAT-P     63 -- E/I-PLAT-P
*               71 -- I-LATE-t                   43 -- I/E-PLAT-T     64 -- E/I-PLAT-T
*                                                44 -- I/E-LR-P       65 -- E/I-LR-P
*                                                45 -- I/E-LR-T       66 -- E/I-LR-T
*                                                46 -- I/E-OTHER-P    67 -- E/I-OTHER-P
*                                                47 -- I/E-OTHER-T    68 -- E/I-OTHER-T
*                                                48 -- I/E-RECRUIT    69 -- E/I-RECRUIT
*
*
*
*
        subroutine picktypetest1 (pattern,cell_name,cell_ID)
*
*
      use mod_miscellaneous_subroutines
        character*15 pattern
        character*4 cell_name
        character*1 choice1,choice2,choice3,choice4,choice5
        integer cell_ID
        character*2 expert
*
*
        expert=' '                      !initialize variables
        pattern=' '
        choice1=' '
        choice2=' '
        choice3=' '
        choice4=' '
        choice5=' '
*
*
*       *****   Display all possible choices.  The user may use the       *****
*       *****           "tree" to designate a unit type by entering <CR>. *****
*
*
11      print'(20(/),''Choose a firing pattern from the following menu''
     +          '' ..or.. <cr> for sub-menus'',//,
     +    T3,''NON-RESP:'',T13,''1 - NRM-SPONT'',T33,
     +       ''2 - NRM-RECRUIT'',
     +    T60, ''OTHER:'',T70,''3 - OTHER'',/,
     +    T3,''---------'',T60,''------'',/,
*
     +  ''I CELLS    P    T     E CELLS   P     T   '',
     +         ''  I/E CELLS     P    T     E/I CELLS     P    T'',/,
     +  ''-------    -    -     -------   -     -   '',
     +         ''  ---------     -    -     ---------     -    -'',/,
     +  ''I-AUG      4    5     E-AUG     17   18   '',
     +         ''  I-I/E-AUG    28   29     I-E/I-AUG    49   50'',/,
     +  ''I-DEC      6    7     E-DEC     19   20   '',
     +         ''  I-I/E-DEC    30   31     I-E/I-DEC    51   52'',/,
     +  ''I-PLAT     8    9     E-PLAT    21   22   '',
     +         ''  I-I/E-OTHER  32   33     I-E/I-OTHER  53   54'',/,
     +  ''I-LowRate 10   11     E-LowRate 23   24   '',
     +         ''  E-I/E-AUG    34   35     E-E/I-AUG    55   56'',/,
     +  ''I-OTHER   12   13     E-OTHER   25   26   '',
     +         ''  E-I/E-DEC    36   37     E-E/I-DEC    57   58'',/,
     +  ''I-RECRUIT    14       E-RECRUIT    27     '',
     +         ''  E-I/E-OTHER  38   39     E-E/I-OTHER  59   60'',/,
     +  ''I-PHR     15   16                         '',
     +         ''  I/E-TRANS    40   41     E/I-TRANS    61   62'',/,
     +  ''                                          '',
     +         ''  I/E-PLAT     42   43     E/I-PLAT     63   64'',/,
     +  ''I-LATE                                    '',
     +         ''  I/E-LowRate  44   45     E/I-LowRate  65   66'',/,
     +  ''          70   71                         '',
     +         ''  I/E-OTHER    46   47     E/I-OTHER    67   68'',/,
     +  ''                                          '',
     +         ''  I/E-RECRUIT     48       E/I-RECRUIT     69'',//,
     +  ''Enter firing pattern'',
     +  '' ..or.. X to exclude this cell from further analysis'',
     +  '' ..or.. <cr> for sub-menus'',/,T15,''>> '',$)'
*
        expert = ' '                                    !default value --> use explorer mode
        read (*,fmt='(A)',err=11),expert
*
        if(expert.eq.' ')then
           goto 10                                      !go to sub-menu tree
c       if(expert.gt.'69')goto 11                               !force an appropriate answer
*
        elseif(expert.eq.'1')then
           pattern='NRM-SPONT'
        elseif(expert.eq.'2')then
           pattern='NRM-RECRUIT'
        elseif(expert.eq.'3')then
           pattern='OTHER'
        elseif(expert.eq.'4')then
           pattern='I-AUG-P'
        elseif(expert.eq.'5')then
           pattern='I-AUG-T'
        elseif(expert.eq.'6')then
           pattern='I-DEC-P'
        elseif(expert.eq.'7')then
           pattern='I-DEC-T'
        elseif(expert.eq.'8')then
           pattern='I-PLAT-P'
        elseif(expert.eq.'9')then
           pattern='I-PLAT-T'
        elseif(expert.eq.'10')then
           pattern='I-LR-P'
        elseif(expert.eq.'11')then
           pattern='I-LR-T'
        elseif(expert.eq.'12')then
           pattern='I-OTHER-P'
        elseif(expert.eq.'13')then
           pattern='I-OTHER-T'
        elseif(expert.eq.'14')then
           pattern='I-RECRUIT'
        elseif(expert.eq.'15')then
           pattern='I-PHR-P'
        elseif(expert.eq.'16')then
           pattern='I-PHR-T'
        elseif(expert.eq.'70')then
           pattern='I-LATE-P'
        elseif(expert.eq.'71')then
           pattern='I-LATE-T'
*
        elseif(expert.eq.'17')then
           pattern='E-AUG-P'
        elseif(expert.eq.'18')then
           pattern='E-AUG-T'
        elseif(expert.eq.'19')then
           pattern='E-DEC-P'
        elseif(expert.eq.'20')then
           pattern='E-DEC-T'
        elseif(expert.eq.'21')then
           pattern='E-PLAT-P'
        elseif(expert.eq.'22')then
           pattern='E-PLAT-T'
        elseif(expert.eq.'23')then
           pattern='E-LR-P'
        elseif(expert.eq.'24')then
           pattern='E-LR-T'
        elseif(expert.eq.'25')then
           pattern='E-OTHER-P'
        elseif(expert.eq.'26')then
           pattern='E-OTHER-T'
        elseif(expert.eq.'27')then
           pattern='E-RECRUIT'
*
        elseif(expert.eq.'28')then
           pattern='I-I/E-AUG-P'
        elseif(expert.eq.'29')then
           pattern='I-I/E-AUG-T'
        elseif(expert.eq.'30')then
           pattern='I-I/E-DEC-P'
        elseif(expert.eq.'31')then
           pattern='I-I/E-DEC-T'
        elseif(expert.eq.'32')then
           pattern='I-I/E-OTHER-P'
        elseif(expert.eq.'33')then
           pattern='I-I/E-OTHER-T'
        elseif(expert.eq.'34')then
           pattern='E-I/E-AUG-P'
        elseif(expert.eq.'35')then
           pattern='E-I/E-AUG-T'
        elseif(expert.eq.'36')then
           pattern='E-I/E-DEC-P'
        elseif(expert.eq.'37')then
           pattern='E-I/E-DEC-T'
        elseif(expert.eq.'38')then
           pattern='E-I/E-OTHER-P'
        elseif(expert.eq.'39')then
           pattern='E-I/E-OTHER-T'
        elseif(expert.eq.'40')then
           pattern='I/E-TRANS-P'
        elseif(expert.eq.'41')then
           pattern='I/E-TRANS-T'
        elseif(expert.eq.'42')then
           pattern='I/E-PLAT-P'
        elseif(expert.eq.'43')then
           pattern='I/E-PLAT-T'
        elseif(expert.eq.'44')then
           pattern='I/E-LR-P'
        elseif(expert.eq.'45')then
           pattern='I/E-LR-T'
        elseif(expert.eq.'46')then
           pattern='I/E-OTHER-P'
        elseif(expert.eq.'47')then
           pattern='I/E-OTHER-T'
        elseif(expert.eq.'48')then
           pattern='I/E-RECRUIT'
*       
        elseif(expert.eq.'49')then
           pattern='I-E/I-AUG-P'
        elseif(expert.eq.'50')then
           pattern='I-E/I-AUG-T'
        elseif(expert.eq.'51')then
           pattern='I-E/I-DEC-P'
        elseif(expert.eq.'52')then
           pattern='I-E/I-DEC-T'
        elseif(expert.eq.'53')then
           pattern='I-E/I-OTHER-P'
        elseif(expert.eq.'54')then
           pattern='I-E/I-OTHER-T'
        elseif(expert.eq.'55')then
           pattern='E-E/I-AUG-P'
        elseif(expert.eq.'56')then
           pattern='E-E/I-AUG-T'
        elseif(expert.eq.'57')then
           pattern='E-E/I-DEC-P'
        elseif(expert.eq.'58')then
           pattern='E-E/I-DEC-T'
        elseif(expert.eq.'59')then
           pattern='E-E/I-OTHER-P'
        elseif(expert.eq.'60')then
           pattern='E-E/I-OTHER-T'
        elseif(expert.eq.'61')then
           pattern='E/I-TRANS-P'
        elseif(expert.eq.'62')then
           pattern='E/I-TRANS-T'
        elseif(expert.eq.'63')then
           pattern='E/I-PLAT-P'
        elseif(expert.eq.'64')then
           pattern='E/I-PLAT-T'
        elseif(expert.eq.'65')then
           pattern='E/I-LR-P'
        elseif(expert.eq.'66')then
           pattern='E/I-LR-T'
        elseif(expert.eq.'67')then
           pattern='E/I-OTHER-P'
        elseif(expert.eq.'68')then
           pattern='E/I-OTHER-T'
        elseif(expert.eq.'69')then
           pattern='E/I-RECRUIT'
        elseif(expert.eq.'X'.or.expert.eq.'x')then
           pattern='exclude'
        else
           goto 11
        end if
*
        return
*
*
*       ***** EXPLORER MODE OF TYPE ENTRY *****
*
*
10      choice1=' '
        PRINT '(20(/),T3,''Cell = '',A,'' (ID = '',I3,''):'',//,
     +          T5,''Basic firing pattern:'',
     +          T30,''I -- Inspiratory'',/,
     +          T30,''E -- Expiratory'',/,
     +          T30,''P -- I/E phase spanning'',/,
     +          T30,''S -- E/I phase spanning'',/,
     +          T30,''O -- Other respiratory pattern'',/,
     +          T30,''N -- Not respiratory modulated'',//,
     +          T20,''Enter selection  >> '',$)',cell_name,cell_ID
        read (*,fmt='(A)',err=10),choice1
        call upper_case(choice1,LEN(choice1))
*
        if((choice1.ne.'I').and.                !force a valid response
     +     (choice1.ne.'E').and.
     +     (choice1.ne.'P').and.
     +     (choice1.ne.'S').and.
     +     (choice1.ne.'O').and.
     +     (choice1.ne.'N'))goto 10

*
*       *****   for an I cell:  *****
*
*
        if(choice1.eq.'I')then
           call level3(choice1,choice2,choice3)
           if(choice3.eq.'A')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-AUG-P'
              if(choice4.eq.'T')pattern='I-AUG-T'
           endif
           if(choice3.eq.'D')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-DEC-P'
              if(choice4.eq.'T')pattern='I-DEC-T'
           endif
           if(choice3.eq.'P')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-PLAT-P'
              if(choice4.eq.'T')pattern='I-PLAT-T'
           endif
           if(choice3.eq.'L')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-LATE-P'
              if(choice4.eq.'T')pattern='I-LATE-T'
           endif
           if(choice3.eq.'B')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-LR-P'
              if(choice4.eq.'T')pattern='I-LR-T'
           endif
           if(choice3.eq.'F')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-PHR-P'
              if(choice4.eq.'T')pattern='I-PHR-T'
           endif
           if(choice3.eq.'O')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='I-OTHER-P'
              if(choice4.eq.'T')pattern='I-OTHER-T'
           endif
           if(choice3.eq.'R')pattern='I-RECRUIT'
           return
        endif
*
*
*       *****   for an E cell:  *****
*
*
        if(choice1.eq.'E')then
           call level3(choice1,choice2,choice3)
           if(choice3.eq.'A')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='E-AUG-P'
              if(choice4.eq.'T')pattern='E-AUG-T'
           endif
           if(choice3.eq.'D')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='E-DEC-P'
              if(choice4.eq.'T')pattern='E-DEC-T'
           endif
           if(choice3.eq.'P')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='E-PLAT-P'
              if(choice4.eq.'T')pattern='E-PLAT-T'
           endif
           if(choice3.eq.'L')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='E-LR-P'
              if(choice4.eq.'T')pattern='E-LR-T'
           endif
           if(choice3.eq.'O')then
              call level4(choice4)
              if(choice4.eq.'P')pattern='E-OTHER-P'
              if(choice4.eq.'T')pattern='E-OTHER-T'
           endif
           if(choice3.eq.'R')pattern='E-RECRUIT'
           return
        endif
*
*
*       *****   for an I/E spanning cell:       *****
*
*
        if(choice1.eq.'P')then          
           call level2(choice1,choice2)
           if(choice2.eq.'I')then                               !peak during I
              call level3(choice1,choice2,choice3)
              if(choice3.eq.'A')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-I/E-AUG-P'
                 if(choice4.eq.'T')pattern='I-I/E-AUG-T'
              endif
              if(choice3.eq.'D')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-I/E-DEC-P'
                 if(choice4.eq.'T')pattern='I-I/E-DEC-T'
              endif
              if(choice3.eq.'O')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-I/E-OTHER-P'
                 if(choice4.eq.'T')pattern='I-I/E-OTHER-T'
              endif
           endif
*       
           if(choice2.eq.'E')then                               !peak during E
              call level3(choice1,choice2,choice3)
              if(choice3.eq.'A')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-I/E-AUG-P'
                 if(choice4.eq.'T')pattern='E-I/E-AUG-T'
              endif
              if(choice3.eq.'D')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-I/E-DEC-P'
                 if(choice4.eq.'T')pattern='E-I/E-DEC-T'
              endif
              if(choice3.eq.'O')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-I/E-OTHER-P'
                 if(choice4.eq.'T')pattern='E-I/E-OTHER-T'
              endif
           endif
*
           if(choice2.eq.'T')then       !peak during transition
              call level4(choice4)
              if(choice4.eq.'P')pattern='I/E-TRANS-P'
              if(choice4.eq.'T')pattern='I/E-TRANS-T'
           endif
*       
           if(choice2.eq.'P')then       !plateau
              call level4(choice4)
              if(choice4.eq.'P')pattern='I/E-PLAT-P'
              if(choice4.eq.'T')pattern='I/E-PLAT-T'
           endif
*       
           if(choice2.eq.'L')then       !low rate cell
              call level4(choice4)
              if(choice4.eq.'P')pattern='I/E-LR-P'
              if(choice4.eq.'T')pattern='I/E-LR-T'
           endif
*       
           if(choice2.eq.'O')then       !other type
              call level4(choice4)
              if(choice4.eq.'P')pattern='I/E-OTHER-P'
              if(choice4.eq.'T')pattern='I/E-OTHER-T'
           endif
*       
           if(choice2.eq.'R')pattern='I/E-RECRUIT' !recruited
*       
           return
        endif
*
*
*       *****   for an E/I spanning cell:       *****
*
*
        if(choice1.eq.'S')then
           call level2(choice1,choice2)
*       
           if(choice2.eq.'I')then !peak during I
              call level3(choice1,choice2,choice3)
              if(choice3.eq.'A')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-E/I-AUG-P'
                 if(choice4.eq.'T')pattern='I-E/I-AUG-T'
              endif
              if(choice3.eq.'D')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-E/I-DEC-P'
                 if(choice4.eq.'T')pattern='I-E/I-DEC-T'
              endif
              if(choice3.eq.'O')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='I-E/I-OTHER-P'
                 if(choice4.eq.'T')pattern='I-E/I-OTHER-T'
              endif
           endif
*       
           if(choice2.eq.'E')then !peak during E
              call level3(choice1,choice2,choice3)
              if(choice3.eq.'A')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-E/I-AUG-P'
                 if(choice4.eq.'T')pattern='E-E/I-AUG-T'
              endif
              if(choice3.eq.'D')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-E/I-DEC-P'
                 if(choice4.eq.'T')pattern='E-E/I-DEC-T'
              endif
              if(choice3.eq.'O')then
                 call level4(choice4)
                 if(choice4.eq.'P')pattern='E-E/I-OTHER-P'
                 if(choice4.eq.'T')pattern='E-E/I-OTHER-T'
              endif
           endif
*       
           if(choice2.eq.'T')then !peak during transition
              call level4(choice4)
              if(choice4.eq.'P')pattern='E/I-TRANS-P'
              if(choice4.eq.'T')pattern='E/I-TRANS-T'
           endif
*       
           if(choice2.eq.'P')then !plateau
              call level4(choice4)
              if(choice4.eq.'P')pattern='E/I-PLAT-P'
              if(choice4.eq.'T')pattern='E/I-PLAT-T'
           endif
*       
           if(choice2.eq.'L')then !low rate
              call level4(choice4)
              if(choice4.eq.'P')pattern='E/I-LR-P'
              if(choice4.eq.'T')pattern='E/I-LR-T'
           endif
*       
           if(choice2.eq.'O')then !other type
              call level4(choice4)
              if(choice4.eq.'P')pattern='E/I-OTHER-P'
              if(choice4.eq.'T')pattern='E/I-OTHER-T'
           endif
*       
           if(choice2.eq.'R')pattern='E/I-RECRUIT'
*       
           return
        endif
*
*
*       *****   for a cell that is respiratory, but does not fit any    *****
*       *****           of these categories:                            *****
*
*
        if(choice1.eq.'O')pattern='OTHER'
*
*
*       *****   for a cell that is not respiratory modulated:   *****
*
*
        if(choice1.eq.'N')then
 50        print '(///,T5,''This NRM cell is:'',
     +     T30,''S -- spontaneously active (not stimulus-related)'',/,
     +     T30,''R -- recruited by stimulus'',//,
     +     T20,''Enter selection  >> '',$)'
*
           read (*,fmt='(A)',err=50),choice5
           call upper_case(choice5,LEN(choice5))
           if(choice5.eq.'S')then
              pattern='NRM-SPONT'
           elseif(choice5.eq.'R')then
              pattern='NRM-RECRUIT'
           else
              goto 50
           endif
*       
           return
        endif   
*       
*       
        return
        end
*
*
*
*
*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>
*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>
*
*       ******************************************************************
*       *                                                                *
*       *       subroutine level2 is called only for phase-spanning cells*
*       *                                                                *
*       ******************************************************************
*
*
        subroutine level2(choice1,choice2)
      use mod_miscellaneous_subroutines

        character*1 choice2,choice1

10      choice2 = ' '
*
        if(choice1.eq.'P')then
          print '(/,T30,''I -- peak during inspiration'',/,
     +          T30,''E -- peak during expiration'',/,
     +          T30,''T -- peak during transition from I to E'',/,
     +          T30,''P -- Plateau'',/,
     +          T30,''R -- Recruited'',/,
     +          T30,''L -- Low rate'',/,
     +          T30,''O -- Other'',//,
     +          T20,''Enter selection  >> '',$)'
          endif
*
        if(choice1.eq.'S')then
          print '(/,T30,''I -- peak during inspiration'',/,
     +          T30,''E -- peak during expiration'',/,
     +          T30,''T -- peak during transition from E to I'',/,
     +          T30,''P -- Plateau'',/,
     +          T30,''R -- Recruited'',/,
     +          T30,''L -- Low rate'',/,
     +          T30,''O -- Other'',//,
     +          T20,''Enter selection  >> '',$)'
          endif
*
        read (*,fmt='(A)',err=10),choice2
        call upper_case(choice2,LEN(choice2))
*
        if((choice2.ne.'I').and.        !force a valid response
     +     (choice2.ne.'E').and.
     +     (choice2.ne.'P').and.
     +     (choice2.ne.'T').and.
     +     (choice2.ne.'R').and.
     +     (choice2.ne.'L').and.
     +     (choice2.ne.'O'))goto 10
*
        return
        end


*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>
*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>

        subroutine level3(choice1,choice2,choice3)
      use mod_miscellaneous_subroutines
        character*1 choice1,choice2,choice3
*
*
10      choice3 = ' '
*
*
*       *****   for an I cell:  *****
*
*
        if(choice1.eq.'I')then  
           print '(/,T5,''Rate:'',
     +                  T30,''A -- Augmenting'',/,
     +                  T30,''D -- Decrementing'',/,
     +                  T30,''P -- Plateau'',/,
     +                  T30,''L -- Late'',/,
     +                  T30,''R -- Recruited'',/,
     +                  T30,''B -- Low rate'',/,
     +                  T30,''F -- Phrenic'',/,
     +                  T30,''O -- Other'',//,
     +                  T20,''Enter selection  >> '',$)'
*
           read (*,fmt='(A)',err=10),choice3
           call upper_case(choice3,LEN(choice3))
*
           if((choice3.ne.'A').and. !force a valid response
     +       (choice3.ne.'D').and.
     +       (choice3.ne.'P').and.
     +       (choice3.ne.'R').and.
     +       (choice3.ne.'L').and.
     +       (choice3.ne.'F').and.
     +       (choice3.ne.'O'))goto 10
*
*
           return
        endif
*
*
*       *****   for an E cell:  *****
*
*
        if(choice1.eq.'E')then  
           print '(/,T5,''Rate:'',
     +                  T30,''A -- Augmenting'',/,
     +                  T30,''D -- Decrementing'',/,
     +                  T30,''P -- Plateau'',/,
     +                  T30,''R -- Recruited'',/,
     +                  T30,''L -- Low rate'',/,
     +                  T30,''O -- Other'',//,
     +                  T20,''Enter selection  >> '',$)'
*
           read (*,fmt='(A)',err=10),choice3
           call upper_case(choice3,LEN(choice3))
*
           if((choice3.ne.'A').and. !force a valid response
     +       (choice3.ne.'D').and.
     +       (choice3.ne.'P').and.
     +       (choice3.ne.'R').and.
     +       (choice3.ne.'L').and.
     +       (choice3.ne.'O'))goto 10
*
*
           return
        endif
*
*
*       *****   for a phase-spanning cell:      *****
*
*
c       if(choice1.eq.'P')then
c        if(choice2.eq.'I')then         !for I and E cells
        if((choice1.eq.'P').or.(choice1.eq.'S'))then
           if((choice2.eq.'I').or.(choice2.eq.'E'))then !for I and E cells
              print '(/,T5,''Rate:'',                                   
     +                  T30,''A -- Augmenting'',/,
     +                  T30,''D -- Decrementing'',/,
     +                  T30,''O -- Other'',//,
     +                  T20,''Enter selection  >> '',$)'
*
              read (*,fmt='(A)',err=10),choice3
              call upper_case(choice3,LEN(choice3))
*       
              if((choice3.ne.'A').and. !force a valid response
     +          (choice3.ne.'D').and.
     +          (choice3.ne.'O'))goto 10
              return
           endif
*
           return
        endif
*
*
        return
        end

*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>
*       <<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>

        subroutine level4(choice4)
      use mod_miscellaneous_subroutines
        character*1 choice4

10      print '(/,T30,''P -- Phasic'',/,
     +          T30,''T -- Tonic'',//,
     +          T20,''Enter selection  >> '',$)'
*
        read (*,fmt='(A)',err=10),choice4
        call upper_case(choice4,LEN(choice4))
*
        if((choice4.ne.'P').and.(choice4.ne.'T'))goto 10
*
*
        return
        end
      end module mod_picktypetest2
