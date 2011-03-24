#!/usr/bin/perl
# Script to calculate the amount of E.coli contamination in an assembly 
#
# Written by Keith Bradnam
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use FAlite;
use Getopt::Long;

###############################################
# 
#  C o m m a n d   l i n e   o p t i o n s
#
###############################################

my $csv;   # produce CSV output file of results
GetOptions ("csv"     => \$csv);

# check we have a suitable input file
my $usage = "Usage: assemblathon_stats.pl [options] <E. coli genome FASTA file> <gzipped assembly file>
options:
        -csv         produce a CSV output file of all results
";

die "$usage" unless (@ARGV == 2);

my ($ECOLI, $ASSEMBLY) = @ARGV;
my ($assembly_ID) = $ASSEMBLY =~ m/^(\w\d+)/;

# hard coded length - change if you are using a different E. coli genome sequence
my $ECOLI_LENGTH = 4747819;


# want to keep track of number of sequences + total length of assembly
my ($seq_counter, $total_length) = (0, 0);

# will store sequence lengths of scaffolds in hash (using FASTA def as hash key)
my %length;

open(my $fh, "gunzip -c $ASSEMBLY |") or die "Can't open Assembly\n";
my $fasta = new FAlite($fh);
while (my $entry = $fasta->nextEntry) {
	my ($name) = $entry->def =~ /^>(\S+)/;
	my $len = length($entry->seq);
	$length{$name} = $len;
	$total_length += $len;
	$seq_counter++;
}



# format BLAST databases if not already done
unless (-s "$ASSEMBLY.xni")  {system("xdformat -n -I $ASSEMBLY")  == 0 or die "Can't create blast database for $ASSEMBLY\n"}




#  Run BLAST
# add wink option if supported?
my $params = "-S 50 -M 1 -N -1 -Q 3 -R 3 -W 15 -mformat 2 -kap";

# make suitable output file name for blast & csv results
my $output = $ASSEMBLY;
$output =~ s/fa.gz//;
$output .= "blast.out";

my $csv_out;
if($csv){
	my $csv_file = $ASSEMBLY;	
	$csv_file =~ s/.fa.gz//;
	$csv_file .= "_ecoli.csv";
	open($csv_out, ">", "$csv_file") or die "Can't open $csv_file\n";	
	print $csv_out "Assembly,Number of E. coli nt in assembly,% of assembly that is E. coli,Number of contaminated scaffolds,% of scaffolds with contamination,% of E. coli genome present in assembly\n";
}

# only run blast if no blast output file exists
unless (-s $output) {system("blastn $ASSEMBLY $ECOLI $params -o $output") && die "Can't run blastn\n"}




# process BLAST output

# will store positions within each scaffold that are matches to E. coli by masking a text string that will be set to 
# be the same length as the scaffold.'0' will indicate unmatched base and '1' indicates a match. E.g. 
#  00011100000 - represents a scaffold sequence where positions 4-6 had matches (at 95%) to E. coli
my $scaffold;

# will also store the same information from the point of view of E. coli. I.e. what fraction of the E. coli genome
# matches the assembly
my $ecoli .= '0' x $ECOLI_LENGTH;

# keep track of number of sequences that match (at 95%) and how many bases of each scaffold are actually E. coli
my ($counter, $contaminant_bases) = (0, 0);

# need a way of knowing when we are dealing with HSPs from the 'next' match, i.e. when there is a new subject in
# the BLAST output
my $subject = "";

# want matches at 95% identity
my $identity = 95;

print STDERR "\nFinding scaffolds with E. coli contamination\n\n";

# main loop over BLAST output
open(my $blast, "<", "$output") or die "Can't open $output\n";

while (<$blast>) {

	my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos,
		$qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;

	# skip matches < 95% identity
	next unless $pct >= $identity;
	
	# need to keep track of when we move to statistics for a new subject sequence
	# do this differently if we are looking at the first sequence in a file, vs subsequent ones
	if($subject ne $sid && $counter == 0){
		$counter++;
		$subject = $sid;
		$scaffold = '0' x $length{$sid};
	} elsif ($subject ne $sid){

		# print results for current scaffold sequence and reset other variables
		print_contamination_details($scaffold, 'scaffold');
		$counter++;
		$subject = $sid;
		$scaffold = '0' x $length{$sid};
	}
	
	# need to get start/end coordinates and length of match, factoring in that start coordinate is $qe/$se 
	# if match is on negative strand
	my ($qstart, $qend) = ($qs, $qe);
       ($qstart, $qend) = ($qe, $qs) if ($qe < $qs);
	my ($sstart, $send) = ($ss, $se);
       ($sstart, $send) = ($se, $ss) if ($se < $ss);

	my $qlen = ($qend - $qstart) + 1;
	my $slen = ($send - $sstart) + 1;

	# now mask out regions of $scaffold and $ecoli with '1' corresponding to matching region
	substr($ecoli,   $qstart-1,$qlen) = ("1" x $qlen);
	substr($scaffold,$sstart-1,$slen) = ("1" x $slen);
}
close($blast);

if ($scaffold) {
	# need to deal with last subject sequence in blast output
	print_contamination_details($scaffold, 'scaffold');
}
else {
	print "No E. coli contamination\n";
}
# final stats
print "\n";


my $percent1 = sprintf("%.2f", ($contaminant_bases / $total_length) * 100);
my $percent2 = sprintf("%.2f", ($counter / $seq_counter) * 100);

print "$contaminant_bases nt from $total_length nt (%$percent1) in the assembly were contaminated with E. coli\n";
print "\n$counter of $seq_counter scaffolds (%$percent2) contained E. coli contamination\n";



if($csv){
	print $csv_out "$assembly_ID,$contaminant_bases,$percent1,$counter,$percent2,";
}

# can now calculate how much of E. coli genome was in assembly
print_contamination_details($ecoli, 'ecoli');

close($csv_out) if ($csv);
exit;

##################
# subroutines
##################

sub print_contamination_details{
	my ($seq, $type) = @_;
	my $seq_length = length($seq);

	# can calculate number of E. coli bases by using tr to substitute '1' for '1'
	my $ecoli_bases = $seq =~ tr/1/1/;
	my $percent = sprintf("%.2f", ($ecoli_bases/$seq_length) * 100);
	
	# if we are processing scaffolds, we do things a little different than with E. coli
	if($type eq 'scaffold'){
		print "$counter)\tL=$seq_length\t%E. coli=$percent\n";		
		$contaminant_bases += $ecoli_bases;
	} elsif($type eq 'ecoli'){
		print "E. coli genome\tL=$seq_length\t%of genome in scaffolds=$percent\n";
		if($csv){
			print $csv_out "$percent\n";
		}
	}
}

