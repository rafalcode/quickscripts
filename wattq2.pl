#!/usr/bin/env perl
# wattq2.pl Due to not being able to convert to fasat, AND as I'm already parsing it, why not?
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
# we may also want the reverse complement of the target, which would add another target
# # in which case we add another target. But we won't do this in this case.
# push @T, $T[0]->revcom();
# $tcou++;
# Just make sure that revcom worked.
# print substr($T[0]->seq, 1, 10) ."\n";
# print substr($T[1]->seq, 1, 10) ."\n";

my $TSL = $T[0]->length(); # will be the same for both, if all targets are the same. A big assumption that may work.

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
my $FH;
my $FH2; # file handle
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
my $EXN; # exon name

# For each target:
my $ENDRANGE=0;
for $K (0..$ENDRANGE) {
	# this is some sort of alt options
	$water->run({-asequence => $T[$K], -bsequence => \@Q, -gapopen => '10.0', -gapextend => '2.5', -outfile => $OFNS[$K]});

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
        } elsif( @AT = ($_ =~ /^[^# ]\S+\s+(\d+)\s+(\S+)\s+(\d+)/)) { 
			# note above how we captured the sequence too
            if($AIN==0) {
                $AIN=1;
                $ACOU++;
            }
            push @A, @AT;
        } elsif( $AIN && /^#/) {
            $J=4*($ACOU-1);
			# for my $y (@A) {
			# 	print "$y\n";
			# }

            $B[$J]=$A[0]; # start index of 1st sequence
            $B[$J+1]=$A[$#A-3]; # end index of 1st sequence
            $B[$J+2]=$A[3]; # start index of 2nd sequence
            $B[$J+3]=$A[$#A]; # end index of 2nd sequence
            $AIN=0;

			# 
			my $SZA=scalar @A;
			open($FH2, ">> $TXTENC", $OFNS[$K]."fa");
			# for sequence A:
			for(my $i=1; $i<$SZA; $i+=6) {
				print $FH2 "$A[$i]\n";
			}
			# for sequence B:
			for(my $i=4; $i<$SZA; $i+=6) {
				print $FH2 "$A[$i]\n";
			}
			close($FH2);
			undef $FH2;
            undef @A;
        }
    }
    close($FH);
}
