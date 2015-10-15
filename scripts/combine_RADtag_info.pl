#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

my $region_size = 5000;

# ----------------------------------------------------------------------------------------
# --- Get RADtag info from original file
# ----------------------------------------------------------------------------------------

my $rad = "results/RADtags.bed";

open (RAD, "<$rad")
	or die "ERROR: Could not open RAD BED file. $!\n";

my %radtags;

# Array to hold original locations, for use when using expanded BED
my @rad_locs;

while (<RAD>) {
	chomp;
	my @info = split /\t/;
	# Store info string keyed by location
	my $loc = "$info[0]:$info[1]-$info[2]";
	$radtags{$loc} = $info[3] . ";";	# Forgot to add final semicolon in ddRAD_sim.R
	
	# Add original location to array
	push @rad_locs, $loc;
}

close RAD;

# ----------------------------------------------------------------------------------------
# --- Associate expanded region coordinates with original RADtag location coordinates
# ----------------------------------------------------------------------------------------

my $expanded_rad = "results/RADtags.slop${region_size}.bed ";

open (EXP, "<$expanded_rad")
	or die "ERROR: Could not open expanded RAD BED file. $!\n";

# Hash to hold original RADtag coordinates, keyed by expanded coords
my %orig_RAD_coord;

my $RADtag_index = 0;

while (<EXP>) {
	chomp;
	my @info = split /\t/;
	my $exp_loc = "$info[0]:$info[1]-$info[2]";

	$orig_RAD_coord{$exp_loc} = $rad_locs[$RADtag_index];
	$RADtag_index++
}

close EXP;

# ----------------------------------------------------------------------------------------
# --- Get GC info on RADtags
# ----------------------------------------------------------------------------------------

my $gc = "results/RADtags.gc.bed";

open (GC, "<$gc")
	or die "ERROR: Could not open RAD BED file with GC percents. $!\n";

while (<GC>) {
	chomp;
	my @info = split /\t/;
	
	my $loc = "$info[0]:$info[1]-$info[2]";
	
	$radtags{$loc} .= "gc_perc:$info[5];";		# GC percentage
	$radtags{$loc} .= "N_count:$info[10];";		# Number of N bases
}

close GC;

# ----------------------------------------------------------------------------------------
# --- Get GC info on expanded region surrounding RADtags
# ----------------------------------------------------------------------------------------

my $gc_regional = "results/RADtags.slop${region_size}.gc.bed";

open (GC_REG, "<$gc_regional")
	or die "ERROR: Could not open RAD BED file with regional GC percents. $!\n";

my $header = <GC_REG>;

while (<GC_REG>) {
	chomp;
	my @info = split /\t/;
	
	my $exp_loc = "$info[0]:$info[1]-$info[2]";
	my $loc = $orig_RAD_coord{$exp_loc};
	
	$radtags{$loc} .= "gc_perc_${region_size}:$info[5];";	# Regional GC percentage
	$radtags{$loc} .= "N_count_${region_size}:$info[10];";	# Regional number of N bases
}

close GC_REG;

# ----------------------------------------------------------------------------------------
# --- Get info on nearest CpG and number of overlapping CpGs
# ----------------------------------------------------------------------------------------

my $cpg_dist = "results/RADtags.nearestCpG.bed";

open (CPG_DIST, "<$cpg_dist")
	or die "ERROR: Could not open RAD BED file with distance to CpG. $!\n";

my %min_CpG_dist;
my %CpG_overlap_count;

while (<CPG_DIST>) {
	chomp;
	my @info = split /\t/;

	my $loc = "$info[0]:$info[1]-$info[2]";
	
	if (exists $min_CpG_dist{$loc}) {
		
		if ($info[7] == 0) {
			$CpG_overlap_count{$loc}++;
		}
		
		if ($info[7] < $min_CpG_dist{$loc}) {
			$min_CpG_dist{$loc} = $info[7];
		}
		
	} else {
	
		$min_CpG_dist{$loc} = $info[7];
		
		if ($info[7] == 0) {
			$CpG_overlap_count{$loc} = 1;
		}
	}
}

close CPG_DIST;

for my $loc (@rad_locs) { 
	
	$radtags{$loc} .= "CpG_dist:" . $min_CpG_dist{$loc} . ";";		# Dist to nearest CpG

	my $CpG_ct = 0;
	if (exists $CpG_overlap_count{$loc}) {
		$CpG_ct = $CpG_overlap_count{$loc};
	}
	
	$radtags{$loc} .= "CpG_ct:" . $CpG_ct . ";";					# Overlap CpG count
}

# ----------------------------------------------------------------------------------------
# --- Get distance to nearest CpG island
# ----------------------------------------------------------------------------------------

my $cpg_is_dist = "results/RADtags.nearestCpGisland.bed";

open (CPG_IS_DIST, "<$cpg_is_dist")
	or die "ERROR: Could not open RAD BED file with distance to CpG island. $!\n";

my %min_CpG_island_dist;

while (<CPG_IS_DIST>) {
	chomp;
	my @info = split /\t/;

	my $loc = "$info[0]:$info[1]-$info[2]";
	
	if (exists $min_CpG_island_dist{$loc}) {
		
		if ($info[8] < $min_CpG_island_dist{$loc}) {
			$min_CpG_island_dist{$loc} = $info[8];
		}
		
	} else {
	
		$min_CpG_island_dist{$loc} = $info[8];
	}
}

close CPG_IS_DIST;

for my $loc (@rad_locs) { 
	
	# Dist to nearest CpG island
	$radtags{$loc} .= "CpG_is_dist:" . $min_CpG_island_dist{$loc} . ";";

}

# ----------------------------------------------------------------------------------------
# --- Get count of CpG sites in expanded region surrounding RADtags
# ----------------------------------------------------------------------------------------

my $cpg_regional = "results/RADtags.slop${region_size}.CpGcount.bed";

open (CPG_REG, "<$cpg_regional")
	or die "ERROR: Could not open RAD BED file with regional CpG counts. $!\n";

while (<CPG_REG>) {
	chomp;
	my @info = split /\t/;
	
	my $exp_loc = "$info[0]:$info[1]-$info[2]";
	my $loc = $orig_RAD_coord{$exp_loc};
	
	$radtags{$loc} .= "CpG_${region_size}:$info[4];";	# Regional number of CpGs
}

close CPG_REG;

# ----------------------------------------------------------------------------------------
# --- Print BED with all gathered info
# ----------------------------------------------------------------------------------------

for my $loc (@rad_locs) { 
	my $tab_loc = $loc;
	$tab_loc =~ s/[:-]/\t/g;
	print $tab_loc . "\t" . $radtags{$loc} . "\n";
}

exit;
