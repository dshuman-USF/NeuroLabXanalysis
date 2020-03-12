#!/usr/bin/perl

if (@ARGV == 0) {
    print "\n";
    print "usage: $0 spikefile on_code off_code > output_file\n\n";
    print "Keeps spikes between each on_code and the first subsequent off code.\n";
    print "Drops analog.\n";
    print "\n";
    exit;
}

$off = pop @ARGV;
$on = pop @ARGV;

print scalar <>;
print scalar <>;

while (<>) {
    ($code) = /^(.....)/;
    next if $code >= 1000;
    if ($code == $on) {
        $do_print = 1;
        print;
    }
    elsif ($code == $off) {
        $do_print = 0;
        print;
    }
    elsif ($do_print) {
        print;
    }
}
