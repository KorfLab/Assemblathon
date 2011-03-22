#!/usr/bin/perl
use strict; use warnings;
use FAlite; use DataBrowser;
use Getopt::Std;
use vars qw($opt_r $opt_s $opt_m $opt_n);
getopts('s:m:n:r:');

my $READS = 1000;
my $SEED  = 1;
my $MIN   = 100;
my $MAX   = 25600;

die "
usage: blastoff2.pl [options] <reference.gz> <assembly.gz>
opitons:
  -m <int> mimimum read distance [$MIN]
  -n <int> maximum read distance [$MAX]
  -r <int> read pairs [$READS]
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

# generate 100 bp paired fragment files if necessary
my %generated;
for (my $r = $MIN; $r <= $MAX; $r*=2) {
	my $frags = "fragments2.$SEED.$r.$READS";
	next if -s $frags;
	print STDERR "generating $READS pairs, $r bp apart, with seed $SEED\n";
	open(my $out, ">$frags") or die;
	foreach my $name (keys %length) {
		my $frac = $length{$name} / $total_length;
		my $reads = int 0.5 + $READS * $frac;
		for (my $i = 0; $i < $reads; $i++) {
			my $pos1 = 1 + int rand($length{$name} - $r - 200);
			my $end1 = $pos1 + 99;
			my ($pos2, $end2) = ($pos1 + $r, $end1 + $r);
			my ($def1, @seq1) = `xdget -n -a $pos1 -b $end1 $REFERENCE $name`;
			my ($def2, @seq2) = `xdget -n -a $pos2 -b $end2 $REFERENCE $name`;
			$def1 =~ s/\s//g;
			chomp @seq1;
			chomp @seq2;
			$generated{$r}++;
			print $out ">L-$generated{$r}\n", @seq1, "\n";
			print $out ">R-$generated{$r}\n", @seq2, "\n";
		}
	}
	close $out;
}
unless (%generated) {
	for (my $r = $MIN; $r <= $MAX; $r*=2) {
		my $count = `grep -c ">" fragments2.$SEED.$r.$READS`;
		$generated{$r} = $count / 2;
	}
}

#  blasts
for (my $r = $MIN; $r <= $MAX; $r*=2) {
	my $frags = "fragments2.$SEED.$r.$READS";
	my $minscore = 90; # 95% identity
	my %hit;
	open(my $blast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags |") or die;
	while (<$blast>) {
		#print;
		my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
			$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se) = split;
		next unless $len >= 95;
		my ($side, $num) = split("-", $qid);
		push @{$hit{$num}{$side}}, {
			parent => $sid,
			start  => $ss,
			end    => $se,
			strand => $qs,
		}
	}
	
	my $tolerance = $r * 0.9;
	my $count = 0;
	OUTER: foreach my $num (keys %hit) {
		my $left  = $hit{$num}{L};
		my $right = $hit{$num}{R};
		next unless defined $left and defined $right;
		my $found = 0;
		foreach my $hsp1 (@$left) {
			foreach my $hsp2 (@$right) {
				next if $hsp1->{parent} ne $hsp2->{parent};
				next if $hsp1->{strand} ne $hsp2->{strand};
				my $distance = abs($hsp1->{start} - $hsp2->{start});
				my $diff = abs($distance - $r);
				if ($diff <= $tolerance) {
					$count++;
					next OUTER;
				}
			}
		}
	}
	printf "%d\t%.4f\n", $r, $count / $generated{$r};
}

