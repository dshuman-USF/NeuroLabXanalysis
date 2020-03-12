#!/usr/bin/env perl

for (`cat x2000parameter.defs`) {
    
    next unless /^ /;
    while (s/([A-Za-z_]+)\s*=\s*(\d+)//) {
        push @txt, "#define $1 $2\n";
    }
}
open F, ">params.h";
print F "#ifndef PARAMS_H\n";
print F "#define PARAMS_H\n";
print F @txt;
print F "#endif\n";
close F;
