#!/usr/bin/perl

if (@ARGV == 0) {
    print <<EOT;

    usage: $0 PS_FILE.. > OUTFILE

    extracts information from xanalysis-generated Postscript files.

    PS_FILE... is a space-separated list of one or more Postscript
    files generated by xanalysis from the single-binwidth CCH display.

    OUTFILE is a text file reporting the experiment name, cell IDs,
    binwidth, min rate, and max rate for each PS_FILE, one per line,
    with tab-separated values.

EOT

    exit (0);
  }

print "\n\n\n";

printf "%s\t%s\t%s\t%s\t%s\t%s\n",
        'experiment', 'ref', 'tar', 'bw', 'maxrate', 'minrate';
print "\n";
for $f (@ARGV) {
    open (F, $f);
    $exp = $f;
    $exp =~ s/\/.*//;
    while (<F>) {
        if (/R: ID = (\d+).* events in spike train = (\d+)/) {
            $ref = $1;
            $refcnt = $2;
        }
        elsif (/T: ID = (\d+).* events in spike train = (\d+)/) {
            $tar = $1;
#            $tarcnt = $1;
        }
        elsif (/max bin = (\d+); min bin = (\d+)/) {
            $max = $1;
            $min = $2;
        }
        elsif (/binwidth: *([0-9.]+)/) {
            $bw = $1;
#            print "$f $max $min $bw $refcnt $tarcnt\n";
            printf "%-10s\t%d\t%d\t%4.1f\t%4.1f\t%4.1f\n",
                        $exp, $ref, $tar, $bw,
                    1000. * $max / $bw / $refcnt,
                    1000. * $min / $bw / $refcnt;
        }
    }
    close F;
}
