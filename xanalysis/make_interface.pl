#!/usr/bin/env perl

$var = "MAX_NUM_EVENTS|MAX_NUM_CHAN|NUM_BINS";

for (`cat x2000parameter.defs`) {
    
    next unless /^ /;
    while (s/($var)\s*=\s*(\d+)//) {
        $pmap{$1} = $2;
    }
}

open F, ">interface.defs";
open C, ">interface.h";

%fmap = (int => "INTEGER", double => "DOUBLE PRECISION", float => "REAL");

for (`cat interface.txt`) {
    next if /^\s*$/;
    if (/^subroutine (\S+)/) {
        $name = $1;
        @{$arg{$name}} = ();
    }
    else {
        ($type, $var, @rest) = split;
        push @{$arg{$name}}, [$type, $var, @rest];
    }
}

for $name (sort keys %arg) {
    print C "void $name (";
    $argcnt = @{$arg{$name}};
    $argn = 0;
    for $a (@{$arg{$name}}) {
        my ($type, $var, @dims) = @$a;
        if (@dims == 0) {
            print C "$type *$var";
        }
        else {
            print C "$type $var";
            for $d (reverse @dims) {
                print C "[$d]";
            }
        }
        print C ", " unless ++$argn == $argcnt;
    }
    print C ");\n";
}

print F "      INTERFACE\n";
for $name (sort keys %arg) {
    print F "         SUBROUTINE $name (";
    $argcnt = @{$arg{$name}};
    $argn = 0;
    for $a (@{$arg{$name}}) {
        my ($type, $var, @dims) = @$a;
        print F "$var";
        print F "," unless ++$argn == $argcnt;
    }
    print F ")\n";
    for $a (@{$arg{$name}}) {
        my ($type, $var, @dims) = @$a;
        for (@dims) {$_ = $pmap{$_} if $pmap{$_};}
        $ftyp = $fmap{$type};
        print F "            $fmap{$type} $var";
        if (@dims > 0) {
            print F " (";
            $, = ",";
            print F @dims;
            print F ")";
        }
        print F "\n";
    }
    print F "         END SUBROUTINE $name\n";
}
print F "      END INTERFACE";

