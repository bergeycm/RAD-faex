#!/usr/bin/perl

use strict;
use warnings;

# ========================================================================================
# --- Downsample blood reads to create blood BAMs that are the same size as fecal
# ========================================================================================

my @inds = ("Tx01", "Tx02", "Tx05", "Tx09", "Tx10");

my $grep_prefix = "grep '^Mapped' reports/";
my $grep_suffix = ".PE.bwa.baboon.aln_stats.passed.realn.txt | ";
$grep_suffix   .= 'sed -e "s/Mapped reads: *\([0-9]*\)\t.*/\1/"';

my $bam_prefix = "results/";
my $bam_suffix = ".PE.bwa.baboon.passed.realn.bam";

foreach my $ind (@inds) {

	print STDERR "PROCESSING $ind:\n";

	my $awk_cmd = "awk -F \",\" '{ if (\$5 == \"" . $ind . "\") print \$0 }' ";
	$awk_cmd .= "../data/individual_info.csv";
	
	print STDERR "$awk_cmd";
	
	my $sample_info = `$awk_cmd`;
	
	my @sample_info = split /\n/, $sample_info;
	
	my $blood_id;
	my @fecal_ids;
	my $blood_count;
	my @fecal_count;
	
	my @blood_downsample_percents;
	my @fecal_downsample_percents;
	
	foreach (@sample_info) {
			
		my @line_info = split /,/;
				
		if ($line_info[7] eq "blood") {
		
			$blood_id = $line_info[0];
			
			# Get count of mapped blood reads
			my $grep_cmd = $grep_prefix . $blood_id . $grep_suffix;
			
			$blood_count = `$grep_cmd`;
			
			print STDERR "Blood reads: $blood_count\n";

		} else {
			my $this_fecal_id = $line_info[0];
			push @fecal_ids, $this_fecal_id;
			
			# Get count of mapped fecal reads
			my $grep_cmd = $grep_prefix . $this_fecal_id . $grep_suffix;

			my $this_fecal_count = `$grep_cmd`;
			chomp $this_fecal_count;
			push @fecal_count, $this_fecal_count;
			
		}
		
	}
	
	print STDERR "Blood ID is $blood_id\n";
	
	my $blood_bam = $bam_prefix . $blood_id . $bam_suffix;
	
	for (my $i = 0; $i < scalar @fecal_ids; $i++) {

		my $fecal_bam = $bam_prefix . $fecal_ids[$i] . $bam_suffix;
	
		my $fecal_perc = $fecal_count[$i] / $blood_count;
	
		print STDERR $fecal_ids[$i] . " has " . $fecal_count[$i] . " reads.\n";
		print STDERR "\t" . $fecal_perc . " of the blood reads.\n";
		
		if ($fecal_perc < 1) {
			# More blood reads
			# Blood will be downsampled
			$blood_downsample_percents[$i] = $fecal_perc;
			$fecal_downsample_percents[$i] = 1;
			
			my $downsampled_blood_bam = $bam_prefix . $blood_id;
			$downsampled_blood_bam .= "_samp-";
			$downsampled_blood_bam .= $blood_id;
#			$downsampled_blood_bam .= lc (($fecal_ids[$i] =~ /fecalRAD-(.*)/)[0]);
			$downsampled_blood_bam .= $bam_suffix;
			
			my $sam_cmd = "samtools view -s " . $fecal_perc;
			$sam_cmd .= " -b " . $blood_bam;
			$sam_cmd .= " > " . $downsampled_blood_bam;
			
			print STDERR "Command: [$sam_cmd]\n";
			system($sam_cmd);
			
		} else {
			# More fecal reads
			# Feces will be downsampled
			$blood_downsample_percents[$i] = 1;
			$fecal_downsample_percents[$i] = 1 / $fecal_perc;
			
			my $downsampled_feces_bam = $bam_prefix . $fecal_ids[$i];
			$downsampled_feces_bam .= "_samp-";
			$downsampled_feces_bam .= $fecal_ids[$i];
#			$downsampled_feces_bam .= lc (($blood_id =~ /fecalRAD-(.*)/)[0]);
			$downsampled_feces_bam .= $bam_suffix;
			
			my $sam_cmd = "samtools view -s " . (1 / $fecal_perc);
			$sam_cmd .= " -b " . $fecal_bam;
			$sam_cmd .= " > " . $downsampled_feces_bam;
			
			print STDERR "Command: [$sam_cmd]\n";
			system($sam_cmd);			
		}
	
	}

	print STDERR "=================================\n";
	
}
