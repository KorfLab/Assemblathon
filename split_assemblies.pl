#!/usr/bin/perl
#
# split_assemblies.pl
#
# a script to process 1 Assemblathon file into 2 files (scaffolds + contigs)
#
# Last updated by: $Author$
# Last updated on: $Date$
##############################################

use strict;
use warnings;
use Keith;
use FAlite;

die "Usage: split_assemblies.pl <assembly scaffolds file.gz>\n" unless (@ARGV == 1);
my $seqs = shift @ARGV;

# how many Ns are we using to split scaffolds into contigs?
my $n = 25;

# will keep track of number of input sequences and output contigs (scaffolded + unscaffolded)
my $seq_count = 0;
my $scaffolded_contig_count = 0;
my $unscaffolded_contig_count = 0;

# need to create output file for contigs
my ($contigs) = $seqs =~ m/(\S+).gz/;
$contigs .= "_contigs.fa";

open(my $output, ">", "$contigs") or die "Can't write to $contigs file\n";



##########################################
#
#    M A I N  loop through FASTA file
#
##########################################

open(my $input, "gunzip -c $seqs |") or die "Can't open $seqs\n";
my $fasta = new FAlite($input);

my $contig_counter = 0;

while(my $entry = $fasta->nextEntry){
	$seq_count++;
    my $seq    = uc($entry->seq);
	my $header = $entry->def;
	
	# if there are not at least 25 consecutive Ns in the sequence we need to split it into contigs
	# otherwise the sequence must be a contig itself and it still needs to be put into a separate file
	if ($seq =~ m/N{25}/){
		# can now split into contigs and keep count of how many we make

		foreach my $contig (split(/N{25,}/, $seq)){
			$contig_counter++;
			$scaffolded_contig_count++;
			
			# can now tidy sequence and suitably modify FASTA header (to make it unique by appending contig counter)
			my $tidied_seq = Keith::tidy_seq($contig);
			my $new_header = $header;

			# take anything up to the first whitespace or boundary and just append '.n' where n is the contig count of this sequence
			$new_header =~ s/\b(.*)\b/$1.$contig_counter/;
			print $output "$new_header\n$tidied_seq\n";
		}
	} else {
		# must be here if the scaffold is actually just a contig (or is a scaffold with < 25 Ns)
		$unscaffolded_contig_count++;

		$contig_counter++;

		# can now tidy sequence and suitably modify FASTA header (to make it unique by appending contig counter)
		my $tidied_seq = Keith::tidy_seq($seq);
		my $new_header = $header;

		# take anything up to the first whitespace or boundary and just append '.n' where n is the contig count of this sequence
		$new_header =~ s/\b(.*)\b/$1.$contig_counter/;
		print $output "$new_header\n$tidied_seq\n";
	}	
}

# just need to tidy up and go home now
my $total_contigs = $scaffolded_contig_count + $unscaffolded_contig_count;
print "Produced $total_contigs contigs ($scaffolded_contig_count scaffolded + $unscaffolded_contig_count unscaffolded) from $seq_count input sequences\n";

close($input);
close($output);

# zip file
system("gzip $contigs") && die "Can't gzip $contigs file\n";
exit(0);
