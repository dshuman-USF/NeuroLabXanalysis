#!/usr/bin/perl -i

while (<>) {
    print;
    last if /LAST U-TURN BEFORE HISTOGRAM FILE IS CREATED.  CONTINUE/;
}
exit 1 if eof;

print <<'EOT';
send -- "y\r"
expect {
# With the following code, all hosts will use the maximum number of
# slots specified in /etc/openmpi/openmpi-default-hostfile.  To use a
# smaller number for some host, uncomment the last line of this
# comment and replace HOSTNAME and SLOTS with appropriate values.  To
# use a smaller number for other hosts, add more lines immediately
# after.
#    -re {	How many processors on HOSTNAME .* >> } {send -- "SLOTS\r" ; exp_continue}

    -re {	How many processors on .* >> } {send -- "\r" ; exp_continue}
    -ex "Press <cr> to continue  >> "          {send -- "\r"}
}
expect -exact "                              >> "
send -- "x\r"
expect eof
EOT

$_ = <>;
exit 0 if /\"y/i;
exit 1
