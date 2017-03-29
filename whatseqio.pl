#!/usr/bin/env perl
# starting with bioperl: what does Bio::SeqIO give
# # The answer is a Seq object, an instanceof a dats structure with a series of members/
use warnings;
use strict;
 
use Bio::SeqIO;
 
# Argument accounting
if($#ARGV != 1) {
    print "Usage error: this script expects two fasta file names as its argument.\n";
    die;
}

my $sf1 = shift @ARGV;
my $sf2 = shift @ARGV;
my $sqio_obj = Bio::SeqIO->new(-file => $sf1);
my $sq_obj   = $sqio_obj->next_seq;
print $sq_obj ."\n"; # this will actually give a hash of course. Use the ->seq member if you want to see the string peraining to the object
