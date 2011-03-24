#!/usr/bin/perl
use strict; use warnings;
use FAlite;
use Getopt::Std;
use vars qw($opt_r $opt_s $opt_m $opt_n);
getopts('s:m:n:r:');

my $READS = 1000;
my $SEED  = 1;
my $MIN   = 100;
my $MAX   = 6400;

die "
usage: blastoff.pl [options] <reference.gz> <assembly.gz>
opitons:
  -m <int> mimimum read size [$MIN]
  -n <int> maximum read size [$MAX]
  -r <int> reads [$READS]
  -s <int> seed [$SEED]
" unless @ARGV == 2;

my ($REFERENCE, $ASSEMBLY) = @ARGV;

$READS = $opt_r if $opt_r;
$SEED  = $opt_s if $opt_s;
$MIN   = $opt_m if $opt_m;
$MAX   = $opt_n if $opt_n;

die "bad seed" unless $SEED == int $SEED and $SEED > 0 and $SEED < 10;
srand($SEED);

# format BLAST databases if not already done
unless (-s "$REFERENCE.xni") {system("xdformat -n -I $REFERENCE") == 0 or die}
unless (-s "$ASSEMBLY.xni")  {system("xdformat -n -I $ASSEMBLY")  == 0 or die}

# find sequence lengths
open(my $fh, "gunzip -c $REFERENCE |") or die;
my $fasta = new FAlite($fh);
my $total_length = 0;
my %length;
while (my $entry = $fasta->nextEntry) {
	my ($name) = $entry->def =~ /^>(\S+)/;
	my $len = length($entry->seq);
	$length{$name} = $len;
	$total_length += $len;
}
print STDERR scalar keys %length, " contigs in reference of $total_length bp\n";

# generate fragment files if necessary
my %generated;
for (my $r = $MIN; $r <= $MAX; $r*=2) {
	my %fragname;
	my $frags = "fragments.$SEED.$r.$READS";
	if (-s $frags) {
		$generated{$r} = `grep -c ">" $frags`;
		next;
	}
	print STDERR "generating $READS $r bp fragments with seed $SEED\n";
	open(my $out, ">$frags") or die;
	foreach my $name (keys %length) {
		my $frac = $length{$name} / $total_length;
		my $reads = int 0.5 + $READS * $frac;
		for (my $i = 0; $i < $reads; $i++) {
			my $pos = 1 + int rand($length{$name} - $r);
			my $end = $pos + $r -1;
			my ($def, @seq) = `xdget -n -a $pos -b $end $REFERENCE $name`;
			$def =~ s/\s//g;
			if (defined ($fragname{$def})) {redo}
			else                           {$fragname{$def} = 1}
			chomp @seq;	
			print $out $def, "\n", @seq, "\n";
			$generated{$r}++;
		}
	}
	close $out;
}

#  blasts
for (my $r = $MIN; $r <= $MAX; $r*=2) {
	my $frags = "fragments.$SEED.$r.$READS";
	my $minscore = int $r * 0.9; # 95% identity
	my %hit;
	open(my $blast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags |") or die;
	while (<$blast>) {
		#print;
		my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
			$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;
		next unless $len >= 0.95 * $r; # 95% length
		$hit{$qid}++;
	}
	my $count = keys %hit;
	printf "%d\t%.4f\n", $r, $count / $generated{$r};
}

