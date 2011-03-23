#!/usr/bin/perl
use strict; use warnings;
use FAlite; use DataBrowser;
use Getopt::Std;
use vars qw($opt_r $opt_s $opt_m $opt_n $opt_c $opt_o);
getopts('m:n:r:sco');

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
  -c	   creates CSV file [default OFF]
  -o	   saves blast output file [default OFF]
" unless @ARGV == 2;

my ($REFERENCE, $ASSEMBLY) = @ARGV;

$READS = $opt_r if $opt_r;
$SEED  = $opt_s if $opt_s;
$MIN   = $opt_m if $opt_m;
$MAX   = $opt_n if $opt_n;
my $CSV = $opt_c ? 1 : 0;
my $SAVE = $opt_o ? 1 : 0;

## Variables used for saving csv and blast output
# P1
my ($assembly_tag) = $ASSEMBLY =~ m/(\w\d+)_/;
# contigs or scaffolds
my ($name) = $ASSEMBLY =~ m/_(\w+)/ if ($assembly_tag);
# A2
my ($ref_tag) = $REFERENCE =~ m/(\w\d?)/;

process_csv_output() if ($CSV);

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
	my $frags = $ref_tag ? $ref_tag . ".fragments2.$SEED.$r.$READS" : "fragments2.$SEED.$r.$READS";
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
		my $count = $ref_tag ? `grep -c ">" $ref_tag.fragments2.$SEED.$r.$READS` : `grep -c ">" fragments2.$SEED.$r.$READS` ;
		$generated{$r} = $count / 2;
	}
}

#  blasts
for (my $r = $MIN; $r <= $MAX; $r*=2) {
	my $frags = $ref_tag ? $ref_tag . ".fragments2.$SEED.$r.$READS" : "fragments2.$SEED.$r.$READS";
	my $minscore = 90; # 95% identity
	my %hit;
	my $blast;
	
	if ($SAVE) {
		my $blast_file;
		if ($ref_tag && $assembly_tag && $name) {
			$blast_file =  "$assembly_tag" . "_" . $name . ".$ref_tag.$SEED.$r.$READS.blast.out";
		} else {
			$blast_file = "$ASSEMBLY.$REFERENCE.$SEED.$r.$READS.blast.out";
		}
		unless (-e "$blast_file"){ 
			system("qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags > $blast_file") == 0 or die "can't run qstack.pl";
		}
		open($blast, "<$blast_file") or die "can't open $blast_file";
	} 
	else { open($blast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags |") or die; }
	
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
	print CSV ",$count" if ($CSV); 
}

if ($assembly_tag && $name && $ref_tag eq 'A2' && $CSV) {
	print CSV "\n";
}

sub process_csv_output {
	
	my $csv_file = $assembly_tag ? $assembly_tag . "_"  . $name . "_" . $SEED . "_" . "paired_end_fragments.csv" 
								 : $ASSEMBLY . "_" . $SEED . "_" . "paired_end_fragments.csv"; 
								
	if (-s $csv_file) {
		open(CSV, ">>", $csv_file) or die ("can't open $csv_file");
	} else {
		open(CSV, ">", $csv_file) or die ("can't open $csv_file");
		print CSV "Assembly,Samples";
		for my $genome qw(A A1 A2) {
			for (my $r = $MIN; $r<= $MAX; $r *= 2) {
			print CSV ",$genome"."_"."$r";
			}
		}
		print CSV "\n";
		$assembly_tag ? print CSV "$assembly_tag,$READS" : print CSV "$ASSEMBLY,$READS";
	}
}
			