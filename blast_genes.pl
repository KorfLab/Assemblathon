#!/usr/bin/perl
## A script to calculate coding sequence and individual exon coverage: 
#  create gene and exon fasta files from Species A if not already created, and blast those files against an assembly.

# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Getopt::Std;
use List::Util qw(min max);
use FAlite;
use vars qw($opt_i $opt_w $opt_W $opt_l);
getopts('i:l:w:W:');

# set default values
$opt_i = 95 if (!$opt_i);
$opt_w = 13 if (!$opt_w);
$opt_W = 13 if (!$opt_W);
$opt_l = 95 if (!$opt_l);

die " 
Usage: generate_results.pl [blast_options] <assembly_file> <species_dir>
blast_options:
 -i <minimum percentage identity required (80-100%)> - default = $opt_i
 -l <minimum percentage length of gene that must match assembly contig (0-100%)> - default = $opt_l
 -w <minimum word size> - default = $opt_w
 -W <wink size> - default = $opt_W
 
(species_dir should include Species A A1 and A2 genome fasta and gff files) 
" unless (@ARGV == 2);

my ($assembly_file, $species_dir) = @ARGV;

# Species that we will work with
my @species_A = qw(A A1 A2);
my $num_genes_in_A = 176;

# Convert command-line options to variable names
my ($wink, $word_min, $identity, $min_length) = ($opt_W, $opt_w, $opt_i, $opt_l);

# store CDS coordinates in a hash of hash structure, hash key is gene ID, second level gives strand and chromosome information
# also stores CDS coordinates as an array of number associated with hash key 'coords'
my %cds_details = ();

# store genomes as a hash of hashes, first key is genome type (A, A1, A2), second key is chromosome number (0,1,2), value is sequence
my %genomes;

# check if all genome and gff files are present (2 files per genome, 6 files total)
check_files();

# Hash to store the results
my %results;

# total gene length for A A1 A2
my %totlengths;

my $create_exon_file = 0; 
my %exon_details = ();
my %exon_results;

my $exon_output_file = "all_exons.fa";
if (-s $exon_output_file) {
	print STDERR "$exon_output_file already exists\n";
	read_exons($exon_output_file);
}
else {
	$create_exon_file = 1;
}

my $output_file = "all_genes.fa";
# does genes file exist? If so just need to read contents of that and the associated genome file
if (-s $output_file){
	print STDERR "$output_file already exists\n";
	read_files($output_file);	
}
# otherwise, need to create files
 else{
	print STDERR "Need to create $output_file\n";
	print STDERR "Need to create $exon_output_file\n" if ($create_exon_file);
	extract_genes($output_file, $create_exon_file, $exon_output_file);
}
	

my $date;

##################
#				 #
# BLAST ANALYSIS #
#				 #
##################

$date = `date`;
print STDERR "Blasting EXONS against $assembly_file at $date\n";

blast($exon_output_file);

$date = `date`;
print STDERR "Finished blasting EXONS against $assembly_file at $date\n";


$date = `date`;
print STDERR "Blasting GENES against $assembly_file at $date\n";

blast($output_file);

$date = `date`;
print STDERR "Finished blasting GENES against $assembly_file at $date\n";
print_output();
	
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
	foreach my $genome_type (@species_A){
		die "Can't find ${genome_type}.annots.gff in $species_dir\n"    unless ("${genome_type}.annots.gff"    ~~ @files);
		die "Can't find ${genome_type}.seq.masked.fa in $species_dir\n" unless ("${genome_type}.seq.masked.fa" ~~ @files);
	}		
}	

# If all_genes.fa already exists, read file plus genome file
sub read_files {
	
	my ($genes_file) = @_;
	open GENES, "<$genes_file" or die "can't open $genes_file\n";
	my $fasta_file = new FAlite(\*GENES);
	while(my $entry = $fasta_file->nextEntry){
		my ($id) = $entry->def =~ m/>([\d\w]*)/; 				
		# add sequence details to hash
		#$cds_details{$id}{seq}    = $entry->seq;
		$cds_details{$id}{length} = length($entry->seq);
		
		my ($genome) = $id =~ m/(A\d?)_/; 
		$totlengths{$genome} += length($entry->seq);
	}

	close GENES;

	# Now, process all genome FASTA files
	my $genome_seq = {}; # to contain known genome sequence in an anonymous hash, where key is sequence ID, and value is sequence

	foreach my $genome (@species_A){
		read_genome($genome);
	}
}

sub read_exons {
	
	my ($exons_file) = @_;
	open EXONS, "<$exons_file" or die "can't open $exons_file\n";
	my $fasta_file = new FAlite(\*EXONS);
	while (my $entry = $fasta_file->nextEntry) {
		my ($id) = $entry->def =~ m/>([\d\w]*)/;
		$exon_details{$id}{length} = length($entry->seq);
	}
	close EXONS;
}

sub extract_genes {
		
	#Output the sequences
	my $out_file = shift;
	my $create_exon = shift;
	my $exon_out = shift;
	
	open EXONOUT, ">$exon_out" or die "can't open $exon_out\n" if ($create_exon);
	open OUT, ">$out_file" or die "can't open $out_file\n";			
	
	foreach my $genome (@species_A) {
		
		read_genome($genome);

		# Process GFF file containing known genes of Species A
		open GFF, "<$species_dir/$genome.annots.gff" or die "Can't open $genome.annots.gff\n";
		
		while (<GFF>) {
			next unless ($_ =~ m/CDS/);
	
			my @gff_line = split(/\t/, $_);
	
			# extract gene number to form ID
			my ($gene) = $gff_line[8] =~ m/gene_index (\d+)/;
			
			my $id = $genome . "_" . $gene; # id = A1_35
			
			# add chromosome, strand, length, and coordinate info to hash for each CDS
			($cds_details{$id}{chr})  = $gff_line[0] =~ m/chr(\d+)/;
			$cds_details{$id}{strand} = $gff_line[6];
			push(@{$cds_details{$id}{coords}}, $gff_line[3], $gff_line[4]); 
		}	
		close GFF;
	}
	
	for my $cds (keys %cds_details) {
		# cds = "A_15";
		
		# to store results, split ID into separate genome and CDS ID
		my ($genome) = $cds =~ m/(A\d?)_/; 
		
		if ($create_exon) {
			my $size = (@{$cds_details{$cds}{coords}});
			my $exon_tag = 0; # counts the number of exons in each gene
			for (my $j = 1; $j <= $size; $j+=2) {
				my $i = $j - 1;
				my $id = "$cds" . "_" . $exon_tag;
				my $exon_length = ${$cds_details{$cds}{coords}}[$j] - ${$cds_details{$cds}{coords}}[$i] + 1;
				$exon_details{$id}{length} = $exon_length;
				my $seq = substr ($genomes{$genome}{$cds_details{$cds}{chr}}, ${$cds_details{$cds}{coords}}[$i] - 1, $exon_length);
				$seq = revcomp($seq) if ($cds_details{$cds}{strand} eq '-');
				$seq= tidy_seq($seq);
				print EXONOUT ">$id $cds_details{$cds}{chr} ${$cds_details{$cds}{coords}}[$i]-${$cds_details{$cds}{coords}}[$j]" .
											  " $cds_details{$cds}{strand}\n$seq\n";
				$exon_tag++;
			}
		}
			
			
		
		# get lowest & highest coordiantes for each CDS
		my $min = min(@{$cds_details{$cds}{coords}});
		my $max = max(@{$cds_details{$cds}{coords}});	
		my $cds_length = ($max - $min + 1);
		$cds_details{$cds}{length} = $cds_length;
		
		$totlengths{$genome} += $cds_length;
			
		# extract sequence from chromosome
		my $sequence = substr ($genomes{$genome}{$cds_details{$cds}{chr}}, $min - 1, $cds_length);
		
		# reverse complement sequence?
		$sequence = revcomp($sequence) if ($cds_details{$cds}{strand} eq '-');
				
		# add sequence to hash
		$cds_details{$cds}{seq} = $sequence; 
		#print "$cds\n";
		# tidy sequence
		$sequence = tidy_seq($sequence);
		
		#min - max?
		print OUT ">$cds $cds_details{$cds}{chr} $min-$max $cds_details{$cds}{strand}\n$sequence\n";
	}
	close OUT;
	close EXONOUT if ($create_exon);

}


sub blast {
	
	my ($GENES) = @_;
	# format BLAST databases if not already done
	unless (-s "$GENES.xni")         {system("xdformat -n -I $GENES") == 0 or die}
	unless (-s "$assembly_file.xni") {system("xdformat -n -I $assembly_file")  == 0 or die}
	my $blast_file = ($GENES =~ /genes/) ? $assembly_file . ".gene.blast.out" : $assembly_file . ".exon.blast.out";
	
	unless (-e "$blast_file"){ 
		system("qstack.pl -s 50 -d -g -h1 -i 80 -w $word_min -W $wink $assembly_file $GENES > $blast_file") == 0 or die "Can't run qstack.pl\n";
	}

	open(my $blast, "<$blast_file") or die "can't open $blast_file";
	
	# then open BLAST output
	while (<$blast>) {
		my ($qid, $sid, $E, $N, $s1, $s, $len, $idn, $pos, $sim, $pct, $ppos, $qg, $qgl, $sg, $sgl, $qf, $qs, $qe, $sf, $ss, $se, $gr) = split;
		#print;
		my $gene_length = ($GENES =~ /genes/) ? $cds_details{$qid}{length} : $exon_details{$qid}{length};
		
		# to store results, split QID into separate genome and CDS ID
		my ($genome, $cds_id) = $qid =~ m/(A.?)_(\d+)/; 
		# do we have a match over 95% of the length and of the percent identity of the gene?
		if ( ($len >= ($min_length / 100) * $gene_length) && ($pct >= $identity) ){	
			
			if ($GENES =~ /genes/) {
				$results{$genome}{blast}++;
				$results{$genome}{length} += $gene_length;
			}
			else {
				$exon_results{$genome}{blast}++;
				$exon_results{$genome}{length} += $gene_length;
			}
		}
	}
}


sub print_output {
	
	my ($assem) = $assembly_file =~ m/(\w\d+)/;

	#contigs or scaffolds
	my ($name) = $assembly_file =~ m/_(\w+)/;
	
	my $file_name;
	if ($assem) {
		# Append to results file if exists, else create one
		$file_name = $assem . "_" . $name . "_" . "genes.csv";
	}
	else {
		$file_name = $assembly_file . "_" . "genes.csv";
	}

	if (-e $file_name) {
		open(OUT, ">>", "$file_name") or die "Can't append to $file_name\n";
	} else {
		open(OUT, ">", "$file_name") or die "Can't create $file_name\n";
		print OUT "Assembly name,";
		print OUT "Number of CDSs in A,Total length of CDSs in A,Percent of CDSs in A,Percent of Total length of CDSs in A," . 
				  "Number of Exons in A,Total length of Exons in A,Percent of Exons in A,Percent of Total length of Exons in A,";
		print OUT "Number of CDSs in A1,Total length of CDSs in A1,Percent of CDSs in A1,Percent of Total length of CDSs in A1," .
				  "Number of Exons in A1,Total length of Exons in A1,Percent of Exons in A1,Percent of Total length of Exons in A1,";
		print OUT "Number of CDSs in A2,Total length of CDSs in A2,Percent of CDSs in A2,Percent of Total length of CDSs in A2" .
				  "Number of Exons in A2,Total length of Exons in A2,Percent of Exons in A2,Percent of Total length of Exons in A2\n";
	}
	
	if ($assem) {
		print OUT "$assem";
	}
	else {
		print OUT "$assembly_file";
	}

	foreach my $genome (@species_A) {
		# Print results
		print "Number of matching CDSs in $genome: " . $results{$genome}{blast} . "\n";
		print "Number of ALL CDSs in $genome: " . $num_genes_in_A . "\n";
		print "Total length of matching CDSs in $genome: " . $results{$genome}{length} . "\n";
		print "Total length of ALL CDSs in $genome: " . $totlengths{$genome} . "\n\n";
		
		print "Number of matching exons in $genome: " . $exon_results{$genome}{blast} . "\n";
		print "Number of ALL exons in $genome: " . (keys %exon_details) . "\n";
		print "Total length of matching exons in $genome: " . $exon_results{$genome}{length} . "\n\n";
		if (exists $results{$genome}{blast}) {
			print OUT ",$results{$genome}{blast},$results{$genome}{length},"; 
			printf OUT "%.2lf,%.2lf", ($results{$genome}{blast} / $num_genes_in_A * 100), ($results{$genome}{length} / $totlengths{$genome} * 100);
		} else {
			print OUT ",0,0,0,0";
		}
		if (exists $exon_results{$genome}{blast}) {
			print OUT ",$exon_results{$genome}{blast},$exon_results{$genome}{length},";
			printf OUT "%.2lf, %.2lf", ( $exon_results{$genome}{blast} / (keys %exon_details) ), ($exon_results{$genome}{length} / $totlengths{$genome} * 100);
		}
		else {
			print OUT ",0,0,0,0";
		}
	}
	print OUT "\n";
	close(OUT);

}


sub read_genome {
	
	my $genome = shift;
	
	open GENOME, "<$species_dir/$genome.seq.masked.fa" or die "can't open $genome.seq.masked.fa\n";
	my $fasta_file = new FAlite(\*GENOME);

	while(my $entry = $fasta_file->nextEntry){
		my ($chr) = $entry->def =~ m/\.chr(\d+)/;
		$genomes{$genome}{$chr} = uc $entry->seq; #everything after '.' on fasta header
	}
	close GENOME;
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


sub revcomp{
	my ($seq) = @_;
	$seq = uc($seq);
	my $revcomp = reverse $seq;
	$revcomp =~ tr/ACGTRSKBD/TGCAYWMVH/;
	return ($revcomp);	
}
