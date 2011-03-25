#!/usr/bin/perl
# A script to calculate coding sequence and individual exon coverage: 
# create gene and exon fasta files from Species A if not already created, and blast those files against an assembly.
#
# Written by Ian Korf, Ken Yu, and Keith Bradnam
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Getopt::Std;
use List::Util qw(min max);
use FAlite;
use vars qw($opt_i $opt_w $opt_W $opt_l $opt_e);
getopts('i:l:w:W:e:');

# set default values
$opt_i = 95 if (!$opt_i);
$opt_w = 13 if (!$opt_w);
$opt_W = 13 if (!$opt_W);
$opt_l = 95 if (!$opt_l);
$opt_e = 30 if (!$opt_e);

die " 
Usage: blast_genes.pl [blast_options] <assembly_file> <species_dir>
blast_options:
 -i <minimum percentage identity required (80-100%)> - default = $opt_i
 -l <minimum percentage length of gene that must match assembly contig (0-100%)> - default = $opt_l
 -w <minimum word size> - default = $opt_w
 -W <wink size> - default = $opt_W
 -e <minimum exon size> - default = $opt_e (we can't BLAST very short exons, so just ignore those below min size)

(species_dir should include Species A A1 and A2 genome fasta and gff files) 
" unless (@ARGV == 2);

my ($assembly_file, $species_dir) = @ARGV;

# check if all genome and gff files are present (2 files per genome, 6 files total)
check_files();

# Convert command-line options to variable names
my ($wink, $word_min, $identity, $min_percent_length, $min_exon_length) = ($opt_W, $opt_w, $opt_i, $opt_l, $opt_e);


# will store genomes as a hash of hashes, first key is genome type (A, A1, A2), second key is chromosome number (0,1,2), value is sequence
my %genomes;

# may want to skip very short exons 
my %exons_to_skip;

# genomes that that we will work with. 2 haplotypes + 1 reference version of the two haplotypes
my @haplotypes = qw(A A1 A2);

# keep track of exon and gene count across genome
my %genomes_to_count;

# store CDS coordinates in a hash of hash structure, hash key is gene ID, second level gives strand and chromosome information
# also stores CDS coordinates as an array of number associated with hash key 'coords'
my %cds_details = ();


# Hash to store the results
my %results;

# total gene and exon lengths for A A1 A2
my %seqlengths;

my %exon_details = ();
my %exon_results;


# keep track of exon and gene lengths
my %seq_to_length;

# keep track of which genome exons and genes are from
my %seq_to_genome;



###########################
#
#   M A I N  L O O P
#
###########################


# Read exisitng genes and exons files if they already exist, otherwise extract sequences from chromosome files
my $genes_file = "all_genes.fa";
my $exons_file = "all_exons.fa";

if (-s $genes_file){
	print STDERR "Will use existing $genes_file file\n";
	read_sequences($genes_file, "gene");
} 

if (-s $exons_file){
	print STDERR "Will use existing $exons_file file\n";
	read_sequences($exons_file, "exon");
}

unless(-s $genes_file && -s $exons_file){
	print STDERR "Need to generate $genes_file and/or $exons_file\n";
	extract_sequences($genes_file, $exons_file);	
}

##################
#				 #
# BLAST ANALYSIS #
#				 #
##################

blast($exons_file, "exon");
blast($genes_file, "gene");

print_output("exon");
print_output("gene");	

exit();
	
	

###############
#			  #
# Subroutines #
#			  #
###############
	
sub check_files {
		
	opendir (DIR, $species_dir) or die "Can't open $species_dir\n";
	my @files = readdir(DIR);

	# Check if species A genome and gff files are present
	foreach my $genome_type (@haplotypes){
		die "Can't find ${genome_type}.annots.gff in $species_dir\n" unless ("${genome_type}.annots.gff" ~~ @files);
		
		unless (("${genome_type}.seq.masked.fa" ~~ @files) or ("${genome_type}.seq.masked.fa.gz" ~~ @files)){
			die "Can't find ${genome_type}.seq.masked.fa or {genome_type}.seq.masked.fa.gz in $species_dir\n" 			
		}
	}		
}	


# Either exons or genes
sub read_sequences {
	
	my ($file, $type) = @_;
		
	open SEQS, "<$file" or die "Can't open $file\n";
	my $fasta_file = new FAlite(\*SEQS);
	while(my $entry = $fasta_file->nextEntry){
		
		my $length = length($entry->seq);
		my ($id, $genome) = $entry->def =~ m/>($type\d+) CDS=(A\d?)_/; 				

		# need to skip of short exons but also keep track of them so we can skip them later
		# when we process the BLAST output. This is necessary if someone changes the value of -e 
		# but uses pre-existing BLAST output
		if($type eq 'exon' && $length < $min_exon_length){
			$exons_to_skip{$id} = 1;
			next;
		}
			
		# track length of sequence and haplotype details for the current sequence
		$seq_to_length{$id} = $length;
		$seq_to_genome{$id} = $genome;
		
		# also keep track of sum length of this type of sequence
		$seqlengths{$type}{$genome} += $length;

		# keep track of how many genes/exons there were for this genome
		$genomes_to_count{$type}{$genome}++; 				
	}	
	close SEQS;
}


# will possibly want to create two files
# 1) sequences of full-length genes (excluding UTRs but including introns)
# 2) sequences of individual CDSs (exons)

sub extract_sequences {

	my ($genes_file, $exons_file) = @_;
	
	# now need to read all of the chromosome sequences for each genome
	read_genome();

	# which files do we need to create?
	my ($need_genes, $need_exons) = (0, 0);
	$need_genes = 1 unless (-s $genes_file);
	$need_exons = 1 unless (-s $exons_file);
	
	# will have to create at least one of the following
	open GENES, ">$genes_file" or die "can't open $genes_file\n" if ($need_genes);
	open EXONS, ">$exons_file" or die "can't open $exons_file\n" if ($need_exons);

	my $exon_count = 0;
	
	foreach my $genome (@haplotypes) {
			
		# Process GFF file containing known genes of Species A
		open GFF, "<$species_dir/$genome.annots.gff" or die "Can't open $genome.annots.gff\n";
	
		# keep track of unique genes and exons in each genome
		my %genes;
		my %exons;
		
		while (<GFF>) {
			# only want to consider CDS entries that belong to each gene
			next unless ($_ =~ m/CDS/);
		
			my @gff_line = split(/\t/, $_);

			# extract gene number to form ID
			my ($gene) = $gff_line[8] =~ m/gene_index (\d+)/;
		
			# keep track of unique genes
			$genes{$gene} = 1;
			
			my $id = $genome . "_" . $gene; # id = A1_35
		
			# add chromosome, strand, length, and coordinate info to hash for each CDS
			my ($chr) = $gff_line[0] =~ m/chr(\d+)/;
			my $strand = $gff_line[6];
			my ($start, $end) = ($gff_line[3], $gff_line[4]);
			my $length = $end - $start + 1;

			$cds_details{$id}{chr}    = $chr;
			$cds_details{$id}{strand} = $strand;
			push(@{$cds_details{$id}{coords}}, $start, $end); 
			
			# can now write to exon file if needed
			if ($need_exons){
				# keep track of unique exons for current haploptype (using coordinates as key) and total number of exons (for all haplotypes)
				$exons{"$start $end"} = 1;
				$exon_count++;

				# store details of current length, and total length as well as which haplotype this sequence belongs to
				my $exonid = "exon$exon_count";

				# skip if exons are shorter than a minimum length which might interfere with BLAST
				# but also keep track of exon ID
				if ($length < $min_exon_length){
					$exons_to_skip{$exonid} = 1;
					next;
				}
						
				$seq_to_length{$exonid} = $length;
				$seqlengths{exon}{$genome} += $length;
				$seq_to_genome{$exonid} = $genome;
				
				# extract sequence
				my $seq = substr ($genomes{$genome}{$chr}, $start - 1, $length);

				# reverse complement sequence?
				$seq = revcomp($seq) if ($strand eq '-');

				# tidy sequence
				$seq = tidy_seq($seq);
				print EXONS ">$exonid CDS=$id location=chr$chr $start-$end $strand\n$seq\n";
			}
			
			# keep track of how many sequences there were for this haplotype
			$genomes_to_count{gene}{$genome} = keys %genes;
			$genomes_to_count{exon}{$genome} = keys %exons;
			
		}	
		close GFF;
	}

	close(EXONS) if ($need_exons);
	
	# if we get here and the genes file already exists, then we don't need to go any further
	return unless ($need_genes);

	my $gene_count = 0;

	for my $cds (keys %cds_details) {
		# cds = "A_15";
		
		# to store results, split ID into separate genome and CDS ID
		my ($genome) = $cds =~ m/(A\d?)_/; 
		
		# get lowest & highest coordinates for each CDS
		my $min = min(@{$cds_details{$cds}{coords}});
		my $max = max(@{$cds_details{$cds}{coords}});	
		my $cds_length = ($max - $min + 1);
		$seqlengths{gene}{$genome} += $cds_length;
					
		# extract sequence from chromosome
		my $sequence = substr ($genomes{$genome}{$cds_details{$cds}{chr}}, $min - 1, $cds_length);
		
		# reverse complement sequence?
		$sequence = revcomp($sequence) if ($cds_details{$cds}{strand} eq '-');
				
		# tidy sequence
		$sequence = tidy_seq($sequence);
		
		$gene_count++;
		my $geneid = "gene$gene_count";				
		$seq_to_length{$geneid} = $cds_length;
		$seq_to_genome{$geneid} = $genome;

		print GENES ">$geneid CDS=$cds $cds_details{$cds}{chr} $min-$max $cds_details{$cds}{strand}\n$sequence\n";
		
	}
	close GENES;
}


sub blast {
	
	my ($file, $type) = @_;
	
	# format BLAST databases if not already done
	unless (-s "$file.xni")          {system("xdformat -n -I $file")           == 0 or die "Can't run xdformat on $file"}
	unless (-s "$assembly_file.xni") {system("xdformat -n -I $assembly_file")  == 0 or die "Can't run xdformat on $assembly_file"}
		
	my ($assembly_id) = $assembly_file =~ m/([A-Z]\d{1,2})_/;
	my $blast_file = "$assembly_id.$type.blast.out";

	my $date = `date`;
	unless (-e "$blast_file"){ 
		print STDERR "Starting BLAST job at $date\n";

		# will want to use a lower BLAST score if searching exons
		my $score = 50;
		$score = $min_exon_length if ($type eq 'exon');
		
		# ideally want to specify -m, -n, -q, and -r as well but will have to change qstack.pl first
#		my $params = "-s $score -M 1 -N -1 -Q 3 -R 3 -d -g -h1 -i 80 -w $word_min -W $wink";
		my $params = "-s $score -m 1 -n -1 -q 3 -r 3 -d -g -h1 -i 80 -w $word_min";
		
		system("qstack.pl $params $assembly_file $file > $blast_file") == 0 or die "Can't run qstack.pl\n";
		$date = `date`;
		print STDERR "Finished BLAST job at $date\n";
	}


	# Process BLAST output file
	open(my $blast, "<$blast_file") or die "can't open $blast_file";
	
	# then open BLAST output
	while (<$blast>) {
		my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos, $qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;

		# is this a short exon that we need to skip?
		next if exists $exons_to_skip{$qid};
		
		my $length = $seq_to_length{$qid};
		my $genome = $seq_to_genome{$qid};

		# do we have a match over 95% of the length of the gene/exon with at least 95% identity?
		# if so just keep track of how many queries matched, and the total length of the match
		if ( ($len >= ($min_percent_length / 100) * $length) && ($pct >= $identity) ){	
			$results{$type}{$genome}{blast}++;
			$results{$type}{$genome}{length} += $len;
		}
	}
	close($blast);
}


sub print_output {
	
	my ($type) = @_;

	my ($assembly_id) = $assembly_file =~ m/(\w\d{1,2})_/;	
	my $file_name = "${assembly_id}.$type.csv";
	

	open(OUT, ">", "$file_name") or die "Can't create $file_name\n";

	# print out header line for CSV file
	print OUT "Assembly name,";
	foreach my $genome (@haplotypes) {
		print OUT "Number of ${type}s found in $genome,% of maximum possible ${type}s in $genome,Total length of ${type}s in $genome,% of maximum possible $type length in $genome,"; 
	}
	print OUT "\n";			


	# Print results
	print OUT "$assembly_id";	

	print "\n\n$type results for assembly $assembly_id\n\n";
	foreach my $genome (@haplotypes) {
		print "Processing genome $genome\n";

		print "Total number of ${type}s present in genome $genome: $genomes_to_count{$type}{$genome}\n";
		my $percent1 = sprintf("%.2f", $results{$type}{$genome}{blast} / $genomes_to_count{$type}{$genome} * 100);
		print "Number of ${type}s from genome $genome found in assembly: $results{$type}{$genome}{blast} (%$percent1)\n";
		
		print "Total length of all possible ${type}s in $genome: $seqlengths{$type}{$genome} bp\n";
		my $percent2 = sprintf("%.2f", $results{$type}{$genome}{length} / $seqlengths{$type}{$genome} * 100);
		print "Total length of ${type}s from $genome found in assembly: $results{$type}{$genome}{length} bp (%$percent2)\n";
		
		# print out results to CSV file (but only if there are results)
		if (exists $results{$type}{$genome}{blast}) {
			print OUT ",$results{$type}{$genome}{blast},$percent1,$results{$type}{$genome}{length},$percent2"; 
		} else {
			print OUT ",0,0,0,0";
		}		
		print "\n";
	}
	print OUT "\n";
	close(OUT);

}


# will always need to process process the known genome of species A
# this consists of two haplotype geneomes (A1 + A2), each of which has 3 chromosomes
# but there is also a reference genome too (an averaged version of A1 and A2)
sub read_genome {
	
	foreach my $genome (@haplotypes){
		my $file = "$species_dir/$genome.seq.masked.fa";

		# If filename doesn't exist but a gzipped version does exist, then use that instead
		if (! -e $file && -e "$file.gz"){
			$file = $file . ".gz";
		}
		
		# if dealing with gzipped FASTA file, treat differently
		my $input;
	    if($file =~ m/\.gz$/){
	            open($input, "gunzip -c $file |") or die "Can't open a pipe to $file\n";
	    } else{
	            open($input, "<", "$file") or die "Can't open $file\n";
	    }

		my $fasta_file = new FAlite(\*$input);

		while(my $entry = $fasta_file->nextEntry){
			my ($chr) = $entry->def =~ m/\.chr(\d+)/;
			$genomes{$genome}{$chr} = uc $entry->seq; #everything after '.' on fasta header
		}
		close $input;
	}
}	


#adds a new line character every 60 bases
sub tidy_seq{
    my ($seq) = @_;
    $seq =~ s/[\s\n]//g;
    $seq =~ tr/a-z/A-Z/;
    
    my ($output_seq) = "";
    my (@seq2) = "";
    my ($start,$end);

    @seq2 = split(//,$seq);
    my $length = @seq2;
    my ($to_add) = int($length/60);

    $end = $start= 0;

    foreach (1..$to_add){
        $output_seq .= substr($seq,$start,60);
        $output_seq .= "\n";
        $start += 60;
    }
    $output_seq .= substr($seq,$start);
    return ($output_seq);
}

# reverse complement sequence
sub revcomp{
	my ($seq) = @_;
	$seq = uc($seq);
	my $revcomp = reverse $seq;
	$revcomp =~ tr/ACGTRSKBD/TGCAYWMVH/;
	return ($revcomp);	
}
