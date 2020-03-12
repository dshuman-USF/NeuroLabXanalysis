#!/usr/bin/perl

while (<>) {
    if (($xy) = /^([-\d.]+ [-\d.]+) moveto$/) {
        if ($xy ne $curr) {
            if ($inline) {
                print "stroke\n";
                $inline = 0;
            }
            print;
            $curr = $xy;
        }
    }
    elsif (($gray,$xy) = /^([\d.]+ setgray )?([-\d.]+ [-\d.]+) 2 copy lineto stroke moveto$/) {
        print "$gray$xy lineto\n";
        $inline = 1;
        $curr = $xy;
    }
    else {
        if ($inline) {
            print "stroke\n";
            $inline = 0;
        }
        print;
        $curr = "";
    }
}
