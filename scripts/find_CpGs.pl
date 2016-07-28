#!/usr/bin/perl

use strict;
use warnings;

my $verbose = 0;

my $genome = "NGS-map/genomes/papAnu2/papAnu2.fa";

open (GEN, "<$genome")
	or die "ERROR: Could not open genome FASTA, $genome. $!\n";

my $chr;
my $chr_bp = 0;

# Flag if line ends with C
my $end_C_flag = 0;

my @CpG_locs;

while (<GEN>) {
	chomp;
	
	if (/^>(.*)/) {
		$chr = $1;
		print STDERR "Searching in $chr...\n";
		$chr_bp = 0;
		next;
	}
	
	my $dna = $_;
	$dna =~ s/\s//g;
	
	# Check if first character is a G if the last line's last char was a C
	if ($end_C_flag && $dna =~ /^G/) {
		my $CpG_start = $chr_bp - 1;
		my $CpG_bp = $chr . "\t" . $CpG_start . "\t" . ($CpG_start + 2);
		push @CpG_locs, $CpG_bp;
		
		if ($verbose) {
			print STDERR "CpG splitting the line at = " . $CpG_bp . "\n";
			print STDERR $dna . "\n";
			print STDERR "^\n";
		}
	}
	
	while ($dna =~ /CG/g) {
		
		my $CpG_start = $chr_bp + $-[0];
		my $CpG_bp = $chr . "\t" . $CpG_start . "\t" . ($CpG_start + 2);
		push @CpG_locs, $CpG_bp;
		
		if ($verbose) {
			print STDERR "CpG at $chr_bp + $-[0] = " . $CpG_bp . "\n";
			print STDERR $dna . "\n";
			print STDERR ' ' x $-[0];
			print STDERR "^^\n";
		}
	}
	
	$chr_bp += length $dna;
	
	# Check if line ends with C, which could mean a CpG spans this and the next line
	if ($dna =~ /C$/) {
		$end_C_flag = 1;
		print STDERR $dna . "\n" if $verbose;
	} else {
		$end_C_flag = 0;
	}
}

close GEN;

print $_ . "\n" foreach (@CpG_locs);

exit;
