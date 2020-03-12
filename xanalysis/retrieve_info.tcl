#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

option add *Font {Helvetica -18 bold} startup

set ofile [lindex $argv 0]

set ifile [tk_getOpenFile                            \
	   -filetypes {{"excel files" {.xls .XLS} }} \
	   -initialdir /oberon/experiments                  \
	   -title "SELECT INFO FILE"             \
	  ]

if {$ifile==""} exit
label .lstat -text "Creating $ofile"
pack .lstat
update idletasks
exec echo $ifile > infofilename.txt
exec XLS2CSV.pl $ifile > $ofile

exit
