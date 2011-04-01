#!/usr/bin/perl
# A script to rank assemblies based on overall performance on our metrics.
# A given metric is excluded, desired at a higher value, or lower value. 
# After ranking each column, a sum of ranks for each assembly was produced 
# (notes on sample input file on the bottom).
#
# Written by Keith Bradnam and Ken Yu
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
#
# Last updated by: $Author$
# Last updated on: $Date$

use strict; use warnings;
use Statistics::RankCorrelation;
use List::Util qw(sum);

# arrays to contain info for rows
my @data;
my @headers;
my @xhl;

while(<>){
	chomp;
	if (m/^Assembly/){
		@headers = split(/,/);
		next;
	}
	# eXclude, High, or Low
	if (m/^XHL/) {
		@xhl = split(/,/);
		next;
	}
	push(@data,[split(/,/)]);
}
my @columns;
foreach my $row (@data) {
    for (my $i = 0; $i < @$row; $i++) {
		push(@{$columns[$i]},$$row[$i]);		
  }
}
my %ranks;
for (my $i = 1; $i < @columns; $i++) {
	
	if ($xhl[$i] eq 'x') {
		#print "Skipped $headers[$i]   $xhl[$i]\n";
		next;
	}
	
	#print "$headers[$i]   $xhl[$i]\n";
	
	my @c = @{$columns[$i]};
	my $ranks = Statistics::RankCorrelation::rank(\@c);
	
	if ($xhl[$i] eq 'h') {
		for (my $i = 0; $i < @$ranks; $i++) { $$ranks[$i] = abs( $$ranks[$i] - (@$ranks + 1) )}
	}
	
	for (my $j=0; $j<@{$ranks}; $j++){
		my $assembly = ${$columns[0]}[$j];
		#print "Assembly = $assembly\tRank = ${$ranks}[$j]\n";
		$ranks{$assembly} += ${$ranks}[$j];
	}
}

foreach my $assembly (sort keys %ranks){
	my $average_rank = $ranks{$assembly} / @columns - 1;
	print "$assembly\tsum of ranks = $ranks{$assembly}\taverage rank = ";
	printf "%.2lf\n", $average_rank;
}

# Sample Input file:
=pod
XHL,x,x,x,h,h
Assembly,scaffold %G,scaffold %T,scaffold %N,scaffold %non-ACGTN,Percentage of assembly in scaffolded contigs,
A1,20.28,29.63,0.26,0,65.5,34.5
B1,20.21,29.64,0.26,0,97.4,2.6
B2,20.02,29.43,1.05,0,86.6,13.4
C1,19.46,29.13,2.79,0,98.8,1.2
C2,17.52,26.27,12.46,0,98.4,1.6
=cut
