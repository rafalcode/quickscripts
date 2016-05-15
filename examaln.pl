#!/usr/bin/env perl
# Examine alignments using bioperl and emboss' water aligner. Small scale.
# used for finding gene in Anguilla scaffolds.
use strict;
use warnings;
# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 2) {
    print "Sorry, this script requires a target sequence first and then a multi sequence file of query sequences\n";
    die;
}

use Bio::SeqIO;
my $tin = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
my $qin = Bio::SeqIO->new(-file => $ARGV[1], '-format' => 'Fasta');

# prepare output files here:
my ($I, $J);
my @OFNS;
for $I (0..1) { # (..) idiom actually perl's range operator
    push @OFNS, $ARGV[0];
    $OFNS[$I]=~ s/(.+)\..+/$1\.out$I/;
}
#
my $tcou=0; # target coutn ... we should be secretly prepared for avng more than one
my $qcou=0;
my(@T, @Q); # target and query. Query can be multisequence, while we expect target to be a single sequence, often much bigger than the individual query sequences.

# Look at target sequence file
while ( my $sq = $tin->next_seq() ) {
    push @T, $sq;
}
if(scalar @T != 1) {
    print "Sorry, only want 1 target sequence at the moment\n";
    die;
}
# we also want the reverse complement of the target
push @T, $T[0]->revcom();
$tcou++;
# Just make sure that revcom worked.
# print substr($T[0]->seq, 1, 10) ."\n";
# print substr($T[1]->seq, 1, 10) ."\n";

my $TSL = $T[0]->length(); # will be the same for both

# and now the query sequences.
my $TQL=0;
while ( my $sq = $qin->next_seq() ) {
    push @Q, $sq;
    $qcou++;
    $TQL += $sq->length();
}
my $NQS=scalar @Q; # number of query sequences.

use Bio::Factory::EMBOSS;
my $f = Bio::Factory::EMBOSS->new();
my $water = $f->program('water');

my $K;
use Bio::AlignIO;

my $AIN; # marker ina an alignment
my $ACOU; # number of alignments
my (@AT, @A, @B); # temporary array. per sequence accumulation array, overall start/stop array.
my @ISG; # Idnetity Similarity Gap array.
my $TXTENC = ":encoding(UTF-8)"; # swanky
my $FH; # file handle
my @TAA; # table array of arrays.
my @TAAS; # table array of arrays sorted
my($PET, $PEQ);
my $TSC; # total aln score
my $TSZ;
my @TAH = ("SFI", "ALEN", "SCORE", "IDEN", "IPT", "SIM", "GAPS", "GPT", "TSC", "TEC", "PET", "QSC", "QEC", "QLN", "PEQ"); # header for table
my $R1Z=scalar @TAH; # number of columns
my @FRSTR=("forward-sense", "reverse-sense");
my $alignio_fmt = "emboss";
my $inaln;
my ($TS, $TE); 


my $ENDRANGE=1;
for $K (0..$ENDRANGE) {
    $water->run({-asequence => $T[$K], -bsequence => \@Q, -gapopen => '10.0', -gapextend => '2.5', -outfile => $OFNS[$K]});

    # OK, let's analyze that alignment. It's been written out ot file.
    my $inaln = new Bio::AlignIO(-format => $alignio_fmt, -file => $OFNS[$K]);

    # I want the start and end coords for each of the sequences.
    # by default bioperl doesn't give this, though it is in the output file.
    undef $FH; # this is the file handle we will use.
    open($FH, "< $TXTENC", $OFNS[$K]);
    $AIN=0; # marker for fileptr being in an alignment
    $ACOU=0; # number of alignments
    undef @AT;
    undef @A;
    undef @B;

    while(<$FH>) {
        if( @AT = ($_ =~ /^# Identity:\s+(\d+)/) ) { 
            $ISG[3*$ACOU]=$AT[0];
        } elsif( @AT = ($_ =~ /^# Similarity:\s+(\d+)/) ) { 
            $ISG[3*$ACOU+1]=$AT[0];
        } elsif( @AT = ($_ =~ /^# Gaps:\s+(\d+)/) ) { 
            $ISG[3*$ACOU+2]=$AT[0];
        } elsif( @AT = ($_ =~ /^[^# ]\S+\s+(\d+)\D+(\d+)/)) { 
            if($AIN==0) {
                $AIN=1;
                $ACOU++;
            }
            push @A, @AT;
        } elsif( $AIN && /^#/) {
            $J=4*($ACOU-1);
            $B[$J]=$A[0]; # start index of 1st sequence
            $B[$J+1]=$A[$#A-2]; # end index of 1st sequence
            $B[$J+2]=$A[2]; # start index of 2nd sequence
            $B[$J+3]=$A[$#A]; # end index of 2nd sequence
            $AIN=0;
            undef @A;
        }
    }
    close($FH);

    # Now arrange to store a table in an array of array references
    $ACOU=0; # reusing
    for($J=0; $J<$R1Z; $J++) {
        undef @{$TAA[$J]}
    }
    undef @TAA;
    # TEC target end Coordinate, IPT iden percentage, SFI source file index (whihc is also the alignment index
    $TSC=0;
    # capture table in array of array refs first.
    while ( my $aln = $inaln->next_aln ) {
        $I=3*$ACOU;
        $J=4*$ACOU;
        $PET=100*($B[$J+1] - $B[$J] +1)/$T[0]->length();
        $PEQ=100*($B[$J+3] - $B[$J+2] +1)/$Q[$ACOU]->length();
        # for reverse strand we want coords converted to forward strand.
        $TS=($K%2)? $T[$K]->length() - $B[$J]+1 : $B[$J];
        $TE=($K%2)? $T[$K]->length() - $B[$J+1]+1 : $B[$J+1];
        push @TAA, [$ACOU+1, $aln->length(), $aln->score(), $ISG[$I], 100*$ISG[$I]/$aln->length(), $ISG[$I+1], $ISG[$I+2], 100*$ISG[$I+2]/$aln->length(), ($K%2)? $TE : $TS, ($K%2)? $TS : $TE, $PET, $B[$J+2], $B[$J+3], $Q[$ACOU]->length(), $PEQ];
        $ACOU++;
        $TSC += $aln->score();
    }

    # sort
    undef @TAAS;
    @TAAS = sort {$a->[8] <=> $b->[8]} @TAA;

    # Now print them out
    $TSZ=scalar @TAAS;
    # print header
    for($J=0; $J<$R1Z; $J++) {
        printf "%s", $TAH[$J];
        ($J==$R1Z-1)? print "\n" : print "\t";
    }
    # now print rest of rows.
    for($I=0; $I<$TSZ; $I++) {
        printf "%d\t%d\t%4.1f\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%d\t%3.1f\n", $TAAS[$I][0], $TAAS[$I][1], $TAAS[$I][2], $TAAS[$I][3], $TAAS[$I][4], $TAAS[$I][5], $TAAS[$I][6], $TAAS[$I][7], $TAAS[$I][8], $TAAS[$I][9], $TAAS[$I][10], $TAAS[$I][11], $TAAS[$I][12], $TAAS[$I][13], $TAAS[$I][14];
    }
    printf "Score for %d query sequences (total %d bp) against %s target (%d bp) = %4.2f\n", $NQS, $TQL, $FRSTR[$K%2], $TSL, $TSC;
    if($K==$ENDRANGE) {
        print "Key: SFI src file idx, ALEN aln length, SCORE aln score, IDEN identical bases, IPT percent iden, SIM similar bases, GAPS num gaps, GPT gap percent\n";
        print "\tTSC target start query, TEC target end coord, PET percent of target, QSC Query start coord, QEC query end coord, QLN query aln length, PEQ percent of query\n";
    }
}
