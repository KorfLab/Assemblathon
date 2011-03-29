#!/usr/bin/perl
# An extention of blastoff.pl:
#   Creates two fragment files for each of haplotype A A1 and A2: one with only repeats and the other non-repeats
#   Blasts those files against the assembly to see how well the assembly deals with repeats
#   (the lesser the difference between the two blasts, the higher the 'repeat tolerance' of the assembly)
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
usage: blastoff_repeats_nonrepeats.pl [options] <assembly.gz> <species_dir>
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
my $SAVE = $opt_o ? 1: 0;

die "bad seed" unless $SEED == int $SEED and $SEED > 0 and $SEED < 10;
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
	
	# used for generating output file names
	#my ($ref) = $REFERENCE =~ m/(A\d?)/;

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
	my %reps_generated;
	my %nonreps_generated;
	for (my $r = $MIN; $r <= $MAX; $r *= 2 ) {
		my %fragname;
		my $repfrags = "$species_dir/" . $genome . ".repeat.fragments.$SEED.$r.$READS";
		my $nonrepfrags = "$species_dir/" . $genome . ".nonrepeat.fragments.$SEED.$r.$READS";
		
		$reps_generated{$r} = `grep -c ">" $repfrags` if (-s $repfrags);
		$nonreps_generated{$r} = `grep -c ">" $nonrepfrags` if (-s $nonrepfrags);
		next if ($reps_generated{$r} && $nonreps_generated{$r});

		print STDERR "generating $READS $r bp fragments with seed $SEED\n";
		open(my $outreps, ">$repfrags") or die;
		open(my $outnonreps, ">$nonrepfrags") or die;
		foreach my $name (keys %length) {
			my $frac = $length{$name} / $total_length;
			my $reads = int 0.5 + $READS * $frac;
			$reps_generated{$r} = 0; $nonreps_generated{$r} = 0;
			while ($reps_generated{$r} < $reads || $nonreps_generated{$r} < $reads) {
				my $pos = 1 + int rand($length{$name} - $r);
				my $end = $pos + $r -1;
				my ($def, @seq) = `xdget -n -a $pos -b $end $species_dir/$REFERENCE $name`;
				$def =~ s/\s//g;
				if (defined ($fragname{$def})) {redo}
				else						   {$fragname{$def} = 1}
				chomp @seq;	
				
				my $overlap = check_overlap($def, $genome); #returns true if overlap, false if no overlap
				
				if ($overlap && $reps_generated{$r} < $reads) {
					print $outreps $def, "\n", @seq, "\n";
					$reps_generated{$r}++;
				}
				elsif (!$overlap && $nonreps_generated{$r} < $reads) {
					print $outnonreps $def, "\n", @seq, "\n";
					$nonreps_generated{$r}++;
				}
			}
		}
		close $outreps;
		close $outnonreps;
	}

	
	my $sum_of_diff = 0;
	my ($sumrep, $sumnonrep) = (0,0);
	#  blasts
	print "Size\tRepeatMatch\tNonRepeatMatch\tSumOfDifferences\n";
	for (my $r = $MIN; $r <= $MAX; $r *= 2) {

		my $repfrags = "$species_dir/" . $genome . ".repeat.fragments.$SEED.$r.$READS";
		my $nonrepfrags = "$species_dir/" . $genome . ".nonrepeat.fragments.$SEED.$r.$READS";
		
		my $minscore = int $r * 0.9; # 95% identity
		
		my ($repblast, $nonrepblast);
	
		if ($SAVE) {
			# repeats blast 
			my $blast_file = "$ASSEMBLY.$genome.repeat.$SEED.$r.$READS.blast.out";
			blast($blast_file, "repeats", $minscore, $repfrags);
			open($repblast, "<$blast_file") or die "can't open $blast_file";
			
			# nonrepeats blast 
			$blast_file = "$ASSEMBLY.$genome.nonrepeat.$SEED.$r.$READS.blast.out";
			blast($blast_file, "nonrepeats", $minscore, $nonrepfrags);
			open($nonrepblast, "<$blast_file") or die "can't open $blast_file";	
				
		} else {
			open($repblast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $repfrags |") or die;
			open($nonrepblast, "qstaq.pl -h 0 -s $minscore $ASSEMBLY $nonrepfrags |") or die;
		}
		
		my %rephit = process_blast_output ($repblast, $r);
		my %nonrephit = process_blast_output ($nonrepblast, $r);

		my $repcount = keys %rephit;
		my $nonrepcount = keys %nonrephit;
		$sum_of_diff += $nonrepcount - $repcount;
		$sumrep += $repcount;
		$sumnonrep += $nonrepcount;
		printf "%d\t%d\t\t%d\t\t%d\n", $r, $repcount, $nonrepcount, $sum_of_diff;
		print OUT ",$nonrepcount,$repcount" if ($CSV);
	}
	printf OUT ",%.2lf", ($sumrep/$sumnonrep*100) if ($CSV);
	
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
		
		my $previn = 0;
		
		# creates numerically ascending keys (starting from 0) for each index
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
			
			$currentID = 0 if $previn ne $index;
			
			my $id = $genome . "_" . $chrom . "_" . $index; # id = A1_0_13
			
			# add coordinate info to hash for each repeat
			push(@{$rpt{$id}{$currentID}}, $gff_line[3], $gff_line[4]);
			
			#print "$id  $currentID  $gff_line[3] $gff_line[4]\n";
			
			$currentID++;
			$previn = $index;
			
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
				print OUT ",$genome"."_N_"."$r,$genome"."_R_"."$r";
			}
			print OUT ",$genome"."_%_Repeat_Tolerance";
		}
		print OUT "\n";
	}
	$assem ? print OUT "$assem,$READS" : print OUT "$ASSEMBLY,$READS";
}

	
sub check_overlap {
	
	my $name = shift;
	my $genome = shift;
	my %overlap;
	
	# chromosome number, start and stop coordinates of missed fragment
	my ($chrom) = $name =~ m/chr(\d+)/;
	my ($frag_start) = $name =~ m/SQ(\d+)/;
	my ($frag_stop) = $name =~ m/-(\d+)/;

	my $part_ID = "$genome" . "_" . "$chrom" . "_"; #A1_0_
			
	my ($lowIndex, $highIndex);
	{
		use integer;
		($lowIndex, $highIndex) = ($frag_start/100000, $frag_stop/100000);
	}

	for (my $i = $lowIndex; $i <= $highIndex; $i++) {
	
		my $ID = $part_ID . $i;
				
		my $ref = \%rpt;
		my $size = (keys %{$ref->{$ID}});
		
		for (my $j = 0; $j <= $size - 1; $j++) {	
			## test for any overlap
			# case where fragment is smaller than repeat
			if ($frag_start >= $rpt{$ID}{$j}[0] && $frag_start <= $rpt{$ID}{$j}[-1] ||
				$frag_stop >= $rpt{$ID}{$j}[0] && $frag_stop <= $rpt{$ID}{$j}[-1] ||
				# case where repeat is smaller than fragment
				$rpt{$ID}{$j}[0] >= $frag_start && $rpt{$ID}{$j}[0] <= $frag_stop ||
				$rpt{$ID}{$j}[-1] >= $frag_start && $rpt{$ID}{$j}[-1] <= $frag_stop) {
					
				$overlap{$name}++;
				#print "$name overlapped with $rpt{$ID}{$j}[0] $rpt{$ID}{$j}[-1]\n";
			}
		}
	}
	
	exists $overlap{$name} ? return 1 : return 0;
	
}

sub blast {
	my ($file, $type, $minscore, $frag) = @_;# shift; my $type = shift; my $minscore = shift;
	unless (-e "$file") {
		if ($type eq "repeats") {  
			system("qstaq.pl -h 0 -s $minscore $ASSEMBLY $frag > $file") == 0 or die "can't run qstack.pl";
		}
		else {
			system("qstaq.pl -h 0 -s $minscore $ASSEMBLY $frag > $file") == 0 or die "can't run qstack.pl";
		}
	}
}

sub process_blast_output {
	my $blast_type = shift; my $r = shift;
	my %hits;
	while (<$blast_type>) {
		my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
		$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;
		next unless $len >= 0.95 * $r; # 95% length
		$hits{$qid}++;
	}
	return %hits;
}