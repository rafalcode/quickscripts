#!/usr/bin/env perl
# wattq0.pl first version of template type code to use Bioperl's EMBOSS bindings to align a single target and a multifasta of queries.
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
my $EXN; # exon name

# For each target:
my $ENDRANGE=0;
for $K (0..$ENDRANGE) {
	# this is some sort of alt options
	$water->run({-asequence => $T[$K], -bsequence => \@Q, -gapopen => '10.0', -gapextend => '2.5', -outfile => $OFNS[$K]});
	# the following is in the docs and it doesn't work!
	# $water->run({-sequencea => $T[$K], -seqall => \@Q, -gapopen => '10.0', -gapextend => '2.5', -outfile => $OFNS[$K]});

    # OK, let's analyze that alignment. It's been written out ot file.
    # OK, that alignment has been writtenout to a file, let's analyze it alignment. It's been written out ot file.
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
	my @QNA; # query name array.
	my $DSPN;
	my $C2=0; # for counting numseqs in each alignment. WIll be 2 for water and needle. Unfort correspnding names (id) will be the EMBOSS shortened one.
    while ( my $aln = $inaln->next_aln ) {
		# exploring SimpleAlign .. this each seq is weird. OK, workied in out .. this is because emboss shortens the sequence names
		foreach my $sq ($aln->each_seq) {
			if($C2==1) {
				push @QNA, $sq->id();
			}
			$C2=$C2+1;
		}
		$C2=0;
		# also did not work
		# foreach my $sq (@QNA) {
		# 	print "$sq\n";
		# }
		# $DSPN=$aln->displayname();
		# print "$DSPN\n";
        $I=3*$ACOU;
        $J=4*$ACOU;
        $PET=100*($B[$J+1] - $B[$J] +1)/$T[0]->length();
        $PEQ=100*($B[$J+3] - $B[$J+2] +1)/$Q[$ACOU]->length();
        # for reverse strand we want coords converted to forward strand.
        $TS=($K%2)? $T[$K]->length() - $B[$J]+1 : $B[$J];
        $TE=($K%2)? $T[$K]->length() - $B[$J+1]+1 : $B[$J+1];
        $EXN=sprintf("e%02d", $ACOU+1);
        push @TAA, [$QSNS[$ACOU], $aln->length(), $aln->score(), $ISG[$I], 100*$ISG[$I]/$aln->length(), $ISG[$I+1], $ISG[$I+2], 100*$ISG[$I+2]/$aln->length(), ($K%2)? $TE : $TS, ($K%2)? $TS : $TE, $PET, $B[$J+2], $B[$J+3], $Q[$ACOU]->length(), $PEQ];
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
        # printf "%s", $TAH[$J];
        # ($J==$R1Z-1)? print "\n" : print "\t";
    }
    # now print rest of rows, but first, header
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "QN", "ALEN", "SCORE", "IDEN", "IPT", "SIM", "GAPS", "GPT", "TSC", "TEC", "PET", "QSC", "QEC", "QLN", "PEQ";
    for($I=0; $I<$TSZ; $I++) {
        printf "%s\t%d\t%4.1f\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%d\t%3.1f\n", $TAAS[$I][0], $TAAS[$I][1], $TAAS[$I][2], $TAAS[$I][3], $TAAS[$I][4], $TAAS[$I][5], $TAAS[$I][6], $TAAS[$I][7], $TAAS[$I][8], $TAAS[$I][9], $TAAS[$I][10], $TAAS[$I][11], $TAAS[$I][12], $TAAS[$I][13], $TAAS[$I][14];
    }
    printf "Score for %d query sequences (total %d bp) against %s target (%d bp) = %4.2f\n", $NQS, $TQL, $FRSTR[$K%2], $TSL, $TSC;
    if($K==$ENDRANGE) {
        print "Legend:\n";
		print "QN Query name/ ALEN alignment length / SCORE alignment score / IDEN # identical bases / IPT %identity over alnlen / SIM similar bases / GAPS #gaps, GPT %gap\n";
        print "TSC target start coord / TEC target end coord / PET %target in aln / QSC Query start coord / QEC query end coord / QLN query aln length / PEQ %query in aln\n";
		print "SIM and IDEN the same in DNA, may be different in proteins.";
    }
	# stpt($TSZ, \@TAAS);
}

