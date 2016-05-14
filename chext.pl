#!/usr/bin/env perl
# Change the extension of an input filename
use strict;
use warnings;
# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Error: one string argument required.\n";
    die;
}

my $OFN=$ARGV[0];
$OFN=~ s/(.+)\..+/$1\.out/; # greedy matching will catch right most dot
print $OFN ."\n";

# now lets say we want several output files from he one arugment name
my ($I, @OFNS);
for $I (0..2) { # (..) idiom actually perl's range operator
    push @OFNS, $ARGV[0];
    $OFNS[$I]=~ s/(.+)\..+/$1\.out$I/;
}
 
for $I (@OFNS) {
    print $I ." ";
}
print "\n";
