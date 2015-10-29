#!/usr/bin/perl

# ----------------------------------------------------------------------------------------
# --- Quantify allelic dropout by comparing blood and fecal DNA from single individual
# ----------------------------------------------------------------------------------------

# Worth it to save STDERR with:
# perl scripts/quantify_ADO.pl 2> reports/quantify_ADO_stderr.txt

use strict;
use warnings;

my @inds = ("Tx01", "Tx02", "Tx05", "Tx07", "Tx09", "Tx10");

my $ind_info_file = "data/fecalRAD_individual_info.csv";

my $prefix = "../NGS-map/results/";
my $suffix = ".PE.bwa.baboon.passed.realn.flt.vcf";

# ----------------------------------------------------------------------------------------

open (INFO, "<$ind_info_file")
	or die "ERROR: Could not open file with individual info, $ind_info_file. $!\n";

my $header = <INFO>;

# All of the following keyed by NGS_ID
my %sample_id;
my %ind_id;
my %sample_type;

# This connects sample ID to NGS ID
# i.e. %NGS_ID_hash{Tx05B} = "fecalRAD-BC3-BC8"
my %NGS_ID_hash;

# This connects individual ID to lists of fecal samples
# i.e. %poop_list{Tx05} = ["Tx05F", "Tx05H", "Tx05H2"]
my %poop_list;

while (<INFO>) {

	my @info = split /,/;

	my $NGS_ID            = $info[0];
	my $this_sample_id    = $info[3];
	my $this_ind_id       = $info[4];
	my $this_sample_type  = $info[5];
		
	$sample_id{$NGS_ID}   = $this_sample_id;
	$ind_id{$NGS_ID}      = $this_ind_id;
	$sample_type{$NGS_ID} = $this_sample_type;
	
	$NGS_ID_hash{$this_sample_id} = $NGS_ID;
	
	print STDERR "NGS ID $NGS_ID, AKA $sample_id{$NGS_ID}, is from ";
	print STDERR "individual ${ind_id{$NGS_ID}}'s $sample_type{$NGS_ID}.\n";

	if ($this_sample_type eq "feces") {		
		push(@{$poop_list{$this_ind_id}}, $this_sample_id);
	}

}

for my $ind (@inds) {

	print STDERR "=" x 80;
	print STDERR "\n";
	print STDERR "Processing individual $ind...\n";
	
	my $blood_id = $ind . "B";
	print STDERR "\tBlood is $blood_id with NGS ID $NGS_ID_hash{$blood_id}\n";
	
	my @poops = @{$poop_list{$ind}};
	print STDERR "\tPoops are:\n";
	foreach my $poop_id (@poops) {
		
		print STDERR "-" x 80;
		print STDERR "\n";
		print STDERR "\t\tComparing to sample $poop_id (NGS ID $NGS_ID_hash{$poop_id})\n";
		
		my $blood_vcf = $prefix . $NGS_ID_hash{$blood_id} . $suffix;
		my $poop_vcf  = $prefix . $NGS_ID_hash{$poop_id}  . $suffix;
		
		my $vcf_cmd_base = "vcftools --vcf $blood_vcf --diff $poop_vcf ";
		$vcf_cmd_base .= "--out results/diff_${blood_id}_${poop_id} ";

		# This fails. Maybe only in later versions of VCFtools
		# Sites that are common / unique to each file. *.diff.sites_in_files
		#my $vcf_cmd_site = $vcf_cmd_base . "--diff-site";
		
		# Discordance on a site by site basis. *.diff.sites
		my $vcf_cmd_site_discord = $vcf_cmd_base . "--diff-site-discordance";

		# Discordance matrix. *.diff.discordance.matrix
		my $vcf_cmd_discord_matrix = $vcf_cmd_base . "--diff-discordance-matrix";

		print STDERR "CMD: [" . $vcf_cmd_site_discord . "]\n";
		system($vcf_cmd_site_discord);
		print STDERR "CMD: [" . $vcf_cmd_discord_matrix . "]\n";
		system($vcf_cmd_discord_matrix);	
	}
}

exit;
