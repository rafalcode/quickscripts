#!/usr/bin/env perl
# wat2fa.pl Due to not being able to convert to fasat, AND as I'm already parsing it, why not?
use strict;
use warnings;

# subroutines with prototypes up first
# Special function to print Exon separation string in terms of qury string, not alignment strinh, which means gaps in target not accounted for.
# stpt Mne: STring PrinT maybe?
# # could it be the query separation string
sub stpt($;$@) {
    my($TSZ, $taref)= @_;
    my ($I, $J);
    my @SEA; # start end array
    for($I=0; $I<$TSZ; $I++) {
        push @SEA, ($taref->[$I][8], $taref->[$I][9]);
    }
    print "Exon separation string:\n";
    for($I=0; $I<$TSZ; $I++) {
        $J=$I*2;
        if($I==$TSZ-1) {
            printf "<< e%02d:%d-%d >>\n", $I+1, $SEA[$J], $SEA[$J+1];
        } else {
            printf "<< e%02d:%d-%d >> %d ", $I+1, $SEA[$J], $SEA[$J+1], $SEA[$J+2] - $SEA[$J+1];
        }
    }
}

# Argument accounting ... say how many the script should have.
if(scalar @ARGV != 2) {
    print "Sorry, this script requires a target sequence filename first and then second a fasta file with multiple query sequences\n";
    die;
}

use Bio::SeqIO;
my $tin = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
my $qin = Bio::SeqIO->new(-file => $ARGV[1], '-format' => 'Fasta');

# prepare output files here:
my ($I, $J);
my @OFNS; # Mnemonic Output File NameS
my $ANT=0; #preordained number of target sequences.
for $I (0..$ANT) { # (..) is actually perl's range operator. Variabel works yes, despite syntax highlighting not working so well here.
    push @OFNS, $ARGV[0];
    $OFNS[$I]=~ s/(.+)\..+/$1\.out$I/; # in bash this would be ${i%.*}.out{0,1}
}
#
my $tcou=0; # target count ... we should be secretly prepared for having more than one, though in theory eeach context should only have one target, for the sake of semantics
my $qcou=0; # query count.
my(@T, @Q); # target and query arrays. Query can be multisequence, while we expect target to be a single sequence, often much bigger than the individual query sequences.

# Look at target sequence file
while ( my $sq = $tin->next_seq() ) {
    push @T, $sq;
}
if(scalar @T != 1) {
    print "Sorry, only want 1 target sequence at the moment\n";
    die;
}

# and now get the query sequences.
my $TQL=0;
my @QSNS; # query sequence names
while ( my $sq = $qin->next_seq() ) {
    push @Q, $sq;
    push @QSNS, $sq->id;
    $qcou++;
    $TQL += $sq->length();
}
my $NQS=scalar @Q; # number of query sequences.
print "$NQS\n";

use Bio::Factory::EMBOSS;
my $f = Bio::Factory::EMBOSS->new();
my $water = $f->program('water');

my $K;

# For each target:
my $ENDRANGE=0;
for $K (0..$ENDRANGE) {
	# this is some sort of alt options
	$water->run({-asequence => $T[$K], -bsequence => \@Q, -gapopen => '10.0', -gapextend => '2.5', -outfile => $OFNS[$K]});
}
