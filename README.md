# quickscripts

scripts to satisfy an urgent need.

* chext.pl
trivial script to change extension of filename given as argument

* examaln.pl
Examine alignment. Chiefly takes a large target sequence in a fasta file (first argument)
and a query multiple sequence fasta file (second argument) and runs Smith Waterman local
alignment for each query sequence on the target file.

It also reverses the target sequence and runs the sequences on it so it represents the reverse strand.

* wattq0.pl
Due to emamaln.pl being a bit complicated, wattq0.pl pares down the extras (no revrsing of target)
and just concentrates on basic functionality.
* needtq0.pl

* whatseqio.pl
A toy program.

* COmparing needle and water.
By far the biggest difference is that NW considers the entire target as game, so that starting and ending gaps are counted
so you get a high gap count. Waterman ignores, these, but you also get a lower PET (% of target). In fact almost by definition,
% of target on needle will byy 100% always because of its nature
