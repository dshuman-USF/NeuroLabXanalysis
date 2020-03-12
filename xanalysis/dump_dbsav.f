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

      program dump_dbsav
      use mod_read_and_write_DBSAV

      include 'x2000parameter.defs'
      character*2 file_version
      integer included(MAX_NUM_CODES),
     +     excluded(MAX_NUM_CODES),analyzed_cells(MAX_NUM_CHAN),
     +     analyzed_pairs(MAX_NUM_CHAN,MAX_NUM_CHAN),
     +     perturb(50),
     +     total_num_cells_db

      CHARACTER*4 CELL_NAMES(MAX_NUM_CHAN)

      character*1 perturb_applied,
     +     AA_applied,STA_applied,AA(10),STA(10),
     +     info_prev_imported,per_prev_imported

      CHARACTER*15 resp_type(MAX_NUM_CHAN)

      character*220 DBSAV
      character*250 qdtfilenames

      read *,DBSAV

      call read_DBSAV(file_version,DBSAV,ios,included,excluded,
     +     analyzed_cells,analyzed_pairs,CELL_NAMES,
     +     perturb_applied,
     +     perturb,AA_applied,AA,STA_applied,STA,resp_type,
     +     info_prev_imported,
     +     per_prev_imported,total_num_cells_db,qdtfilenames)

      print *,'file_version: ',file_version
      print *,'included: ',included
      print *,'excluded: ',excluded
      print *,'analyzed_cells: ',analyzed_cells
      print *,'analyzed_pairs: ',analyzed_pairs
      print *,'CELL_NAMES: ',CELL_NAMES
      print *,'perturb_applied: ',perturb_applied
      print *,'perturb: ',perturb
      print *,'AA_applied: ',AA_applied
      print *,'AA: ',AA
      print *,'STA_applied: ',STA_applied
      print *,'STA: ',STA
      print *,'resp_type: ',resp_type
      print *,'info_prev_imported: ',info_prev_imported
      print *,'per_prev_imported: ',per_prev_imported
      print *,'total_num_cells_db: ',total_num_cells_db
      print *,'qdtfilenames: ',qdtfilenames

      end
