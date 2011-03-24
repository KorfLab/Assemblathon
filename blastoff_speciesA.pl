#!/usr/bin/perl
# A Slightly modified version of blastoff.pl:
#   Blasts fragments of A A1 and A2. 
#   Finds whether non-matching fragments overlap a repeat region.
#   Added options to save blast output and produce csv file. 
#
# Written by Ian Korf, Ken Yu, and Keith Bradnam
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use FAlite;
use Getopt::Std;
use vars qw($opt_r $opt_s $opt_m $opt_n $opt_c $opt_o);
getopts('m:n:r:s:co');

my $READS = 1000;
my $SEED  = 1;
my $MIN   = 100;
my $MAX   = 6400;

die "
usage: blastoff_2.pl [options] <assembly.gz> <species_dir>
options:
  -m <int>	 mimimum read size [$MIN]
  -n <int> 	 maximum read size [$MAX]
  -r <int>       reads [$READS]
  -s <int>       seed [$SEED]
  -c 	         creates CSV file [default off]
  -o	         saves blast output [default off]
  
(species_dir should include Species A A1 and A2 genome fasta and gff files) 
" unless @ARGV == 2;

my $date = `date`;
print STDERR "Starts: $date\n";

my ($ASSEMBLY, $species_dir) = @ARGV;

$READS = $opt_r if $opt_r;
$SEED  = $opt_s if $opt_s;
$MIN   = $opt_m if $opt_m;
$MAX   = $opt_n if $opt_n;
my $CSV = $opt_c ? 1 : 0;
my $SAVE = $opt_o ? 1 : 0;

die "bad seed" unless $SEED == int $SEED and $SEED > 0 and $SEED < 11;
srand($SEED);

# hash to store repeat info. for all genomes
my %rpt;

# hard-coding the applicable genomes..
my @haplotypes = qw(A A1 A2);

# check if all genome and gff files are present (3 files per genome, 9 files total)
check_files();

# create csv file if requested
process_output_file() if ($CSV);

read_repeatsGFF();

for my $genome (@haplotypes) {
	
	my $REFERENCE = "$genome.seq.masked.fa.gz";

	print "Searching $REFERENCE against $ASSEMBLY\n";
	
	# format BLAST databases if not already done
	unless (-s "$species_dir/$REFERENCE.xni") {system("xdformat -n -I $species_dir/$REFERENCE") == 0 or die}
	unless (-s "$ASSEMBLY.xni")  {system("xdformat -n -I $ASSEMBLY")  == 0 or die}

	# find sequence lengths
	open(my $fh, "gunzip -c $species_dir/$REFERENCE |") or die;
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
	for (my $r = $MIN; $r <= $MAX; $r *= 2 ) {
		my %fragname;
		my $frags = "$species_dir/" . $genome . ".fragments.$SEED.$r.$READS";
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
				my ($def, @seq) = `xdget -n -a $pos -b $end $species_dir/$REFERENCE $name`;
				$def =~ s/\s//g;
				if (defined ($fragname{$def})) {redo}
				else						   {$fragname{$def} = 1}
				chomp @seq;	
				print $out $def, "\n", @seq, "\n";
				$generated{$r}++;
			}
		}
		close $out;
	}
	
	#  blasts
	print "Size\tMatch\t%Match\tNon-match\tOverlap\t%Overlap\n";
	for (my $r = $MIN; $r <= $MAX; $r *= 2) {
		my $frags = "$species_dir/" . $genome . ".fragments.$SEED.$r.$READS";
		my $minscore = int $r * 0.9; # 95% identity
		my %hit;
		my $blast;
	
		if ($SAVE) {
			my $blast_file = "$ASSEMBLY.$genome.$SEED.$r.$READS.blast.out";
			unless (-e "$blast_file"){ 
				system("qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags > $blast_file") == 0 or die "can't run qstack.pl";
			}
			open($blast, "<$blast_file") or die "can't open $blast_file";
		} else {
			open($blast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $frags |") or die;
		}

		while (<$blast>) {
			#print;
			my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
				$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;
			next unless $len >= 0.95 * $r; # 95% length
			$hit{$qid}++;
		}
	
		my $count = keys %hit;
		printf "%d\t%d\t%.2f\t", $r, $count, $count / $generated{$r};
		print OUT ",$count" if ($CSV);
	
		
		# number of fragments that did not match assembly
		my $miss = $READS - $count;
		# hash to store fragments overlapping A/A1/A2 repeats
		my %overlap;
		# to be sure that we're only dealing with missed fragments
		my %check;
		
		open(my $in, "<$frags") or die "can't open $frags";
		$fasta = new FAlite($in);
		
		while (my $entry = $fasta->nextEntry) {
			
			my ($name) = $entry->def =~ /^>(\S+)/;

			# only process headers that aren't hits
			if (exists $hit{$name}){
				next;
			}
			else {
				$check{$name}++;
			}
			
			# chromosome number, start and stop coordinates of missed fragment
			my ($chrom) = $name =~ m/chr(\d+)/;
			my ($frag_start) = $name =~ m/SQ(\d+)/;
			my ($frag_stop) = $name =~ m/-(\d+)/;

			my $part_ID = "$genome" . "_" . "$chrom" . "_"; #A1_0_
			
			# determine repeat coverage in terms of the 'index' (defined when parsing the repeat gff)
			my ($lowIndex, $highIndex);
			{
				use integer;
				($lowIndex, $highIndex) = ($frag_start/100000, $frag_stop/100000);
			}
			
			# loop through coverage
			for (my $i = $lowIndex; $i <= $highIndex; $i++) {
				
				# complete hash key to access coordinates
				my $ID = $part_ID . $i;
				
				my $ref = \%rpt;
				
				# get the number of pairs of coordinates
				my $size = (keys %{$ref->{$ID}});
				
				# loop through coordinates and test for overlap (by 1 or more basepair)
				for (my $j = 0; $j <= $size - 1; $j++) {
					
					# start and stop coordinates of repeat region
					my $rpt_start = $rpt{$ID}{$j}[0];
					my $rpt_stop = $rpt{$ID}{$j}[-1];
					
						# case where repeat is larger than fragment
						# -------- repeat
						#   ----   fragment
					if ($frag_start >= $rpt_start && $frag_start <= $rpt_stop || # start of fragment is within repeat or
						$frag_stop >= $rpt_start && $frag_stop <= $rpt_stop || # stop of fragment is within repeat
						# case where fragment is larger than repeat
						# -------- fragment
						#   ----   repeat
						$rpt_start >= $frag_start && $rpt_start <= $frag_stop || # vice versa of above
						$rpt_stop >= $frag_start && $rpt_stop <= $frag_stop) {
					
						$overlap{$name}++;
						#print "$name overlapped with $rpt_start $rpt_stop\n";
					}
				}
			}
		}
		# number of missed fragments
		my $check_count = keys %check;
		# number of missed fragments that overlapped a repeat
		my $lap_count = keys %overlap;
		printf "%d\t\t%d\t%.2f\n", $check_count, $lap_count, $lap_count/$miss; 			
	}
	print "\n";
}
print OUT "\n" if ($CSV);

$date = `date`;
print STDERR "Ends: $date\n";

#################
## Subroutines ##
#################

sub check_files {
		
	opendir (DIR, $species_dir) or die "Can't open $species_dir\n";
	my @files = readdir(DIR);

	# Check if species A genome and gff files are present
	foreach my $genome_type (@haplotypes){
		die "Can't find ${genome_type}.annots.gff in $species_dir\n"    unless ("${genome_type}.annots.gff"    ~~ @files);
		die "Can't find ${genome_type}.seq.masked.fa.gz in $species_dir\n" unless ("${genome_type}.seq.masked.fa.gz" ~~ @files);
		die "Can't find ${genome_type}.repmask.gff in $species_dir\n" unless ("${genome_type}.repmask.gff" ~~ @files);
	}		
}	

		
sub read_repeatsGFF {	
	
	foreach my $genome (@haplotypes) {
	
		# Process GFF file containing known genes of Species A
		open GFF, "<$species_dir/$genome.repmask.gff" or die "Can't open $genome.repmask.gff\n";
		
		my $previndex = 0;
		
		# creates numerically ascending keys (starting from 0) for each index
		# used for looping through each pair of coordinates
		my $currentID = 0;
		
		while (<GFF>) {

			my @gff_line = split(/\t/, $_);

			# extract chromosome number to form ID
			my ($chrom) = $gff_line[0] =~ m/chr(\d+)/;
			
			# each 100kb region belongs in one bracket 
			my $index; 
			{
				use integer;
				$index = $gff_line[4] / 100000;
			}
			
			$currentID = 0 if $previndex ne $index;
			
			my $id = $genome . "_" . $chrom . "_" . $index; # id = A1_0_13
			
			# add coordinate info to hash for each repeat
			push(@{$rpt{$id}{$currentID}}, $gff_line[3], $gff_line[4]);
			
			#print "$id  $currentID  $gff_line[3] $gff_line[4]\n";
			
			$currentID++;
			$previndex = $index;
			
		}	
		close GFF;
		
	}
}

sub process_output_file {
	

	# creates output file if doesn't exist
	my ($assem) = $ASSEMBLY =~ m/(\w\d+)/;

	#contigs or scaffolds
	my ($name) = $ASSEMBLY =~ m/_(\w+)/;
	
	my $file_name;
	if ($assem) {
		$file_name = $assem . "_" . $name . "_" . "$SEED" . "fragments.csv";
	}
	else {
		$file_name = $ASSEMBLY . "_" . "$SEED" . "fragments.csv";
	}
		
	# Append to results file if exists, else create one
	if (-e $file_name) {
		open(OUT, ">>", $file_name) or die "Can't open $file_name";
	} else {
		open(OUT, ">", $file_name) or die "Can't open $file_name";
		print OUT "Assembly,Samples";
		for my $genome (@haplotypes) {
			for (my $r = $MIN; $r<= $MAX; $r *= 2) {
				print OUT ",$genome"."_"."$r";
			}
		}
		print OUT "\n";
	}
	$assem ? print OUT "$assem,$READS" : print OUT "$ASSEMBLY,$READS";
}

	
