#!/usr/bin/perl

# ========================================================================================
# --- Quantify allelic dropout by comparing blood and fecal DNA from single individual
# ========================================================================================

# Worth it to save STDERR with:
# perl scripts/quantify_ADO.pl \
#   NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz \
#   2> reports/quantify_ADO_stderr.txt

use strict;
use warnings;

my @inds = ("Tx01", "Tx02", "Tx05", "Tx09", "Tx10");

my $ind_info_file = "data/individual_info.csv";

my $prefix = "NGS-map/results/";
my $suffix = ".PE.bwa.baboon.passed.realn.flt.vcf";

# VCF file with SNP calls from GATK
# e.g. "NGS-map/baboon_snps_indiv/baboon.INDIV_DIPLOID.pass.snp.vcf.gz"
if (!defined $ARGV[0]) {
    die "ERROR: Pass VCF file with SNP calls from GATK."
}
my $in_vcf = shift;
chomp $in_vcf;

my $rand = int(rand(10000));

# Optional suffix for output files
my $out_suffix = "";
if ($in_vcf =~ /INDIV/) {
	$out_suffix = "_indiv";
}

my $vcf_prefix = "vcftools --gzvcf $in_vcf --recode --indv ";

# ----------------------------------------------------------------------------------------
# --- Match IDs of paired samples by parsing file with individual info
# ----------------------------------------------------------------------------------------

open (INFO, "<$ind_info_file")
	or die "ERROR: Could not open file with individual info, $ind_info_file. $!\n";

my $header = <INFO>;

# All of the following keyed by NGS_ID
my %sample_id;
my %ind_id;
my %sample_type;

# This connects sample ID to NGS ID
# i.e. %NGS_ID_hash{Tx05B} = "TB05a"
my %NGS_ID_hash;

# This connects individual ID to lists of fecal samples
# i.e. %poop_list{Tx05} = ["Tx05F", "Tx05H", "Tx05H2"]
my %poop_list;

while (<INFO>) {

	my @info = split /,/;

	my $NGS_ID            = $info[0];
	my $this_sample_id    = $info[3];
	my $this_ind_id       = $info[4];
	my $this_sample_type  = $info[7];

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

# ----------------------------------------------------------------------------------------
# --- Compare blood and fecal pairs, original and downsampled
# ----------------------------------------------------------------------------------------

for my $ind (@inds) {

	print STDERR "=" x 80;
	print STDERR "\n";
	print STDERR "Processing individual $ind...\n";

	my $blood_id = $ind . "B";
	print STDERR "\tBlood is $blood_id with NGS ID $NGS_ID_hash{$blood_id}\n";

	my @poops = @{$poop_list{$ind}};

	foreach my $poop_id (@poops) {

		print STDERR "-" x 80;
		print STDERR "\n";
		print STDERR "\t\tComparing to sample $poop_id (NGS ID $NGS_ID_hash{$poop_id})\n";

		# Pull out individuals to make temporary individual VCF files
		my $tmp_vcf_cmd_1 = $vcf_prefix . $NGS_ID_hash{$blood_id} . ".PE ";
		$tmp_vcf_cmd_1 .= "--out tmp${rand}.blood";
		my $tmp_vcf_cmd_2 = $vcf_prefix . $NGS_ID_hash{$poop_id}  . ".PE ";
		$tmp_vcf_cmd_2 .= "--out tmp${rand}.feces";

		print STDERR "CMD 1: [" . $tmp_vcf_cmd_1 . "]\n";
		system($tmp_vcf_cmd_1);
		print STDERR "CMD 2: [" . $tmp_vcf_cmd_2 . "]\n";
		system($tmp_vcf_cmd_2);

		# Compare individual VCF files

		# Write individual mapping file
		open (MAP, ">tmp${rand}.map")
			or die "ERROR: Could not open mapping file. $!\n";
		print MAP $NGS_ID_hash{$poop_id} . ".PE\t" . $NGS_ID_hash{$blood_id} . ".PE\n";
		close MAP;

		my $blood_vcf = "tmp${rand}.blood.recode.vcf";
		my $poop_vcf  = "tmp${rand}.feces.recode.vcf";

		my $vcf_cmd_base = "vcftools --vcf $blood_vcf --diff $poop_vcf ";
		$vcf_cmd_base .= "--diff-indv-map tmp${rand}.map ";
		$vcf_cmd_base .= "--out results/diff_${blood_id}_${poop_id}${out_suffix} ";

		# This fails. Maybe only in later versions of VCFtools
		# Sites that are common / unique to each file. *.diff.sites_in_files
		#my $vcf_cmd_site = $vcf_cmd_base . "--diff-site";

		# Discordance on a site by site basis. *.diff.sites
		my $vcf_cmd_site_discord = $vcf_cmd_base . "--diff-site-discordance";

		# Discordance matrix. *.diff.discordance.matrix
		my $vcf_cmd_discord_matrix = $vcf_cmd_base . "--diff-discordance-matrix";

		print STDERR "CMD 3: [" . $vcf_cmd_site_discord . "]\n";
		system($vcf_cmd_site_discord);
		print STDERR "CMD 4: [" . $vcf_cmd_discord_matrix . "]\n";
		system($vcf_cmd_discord_matrix);

		# Delete temporary VCF files
		system("rm $blood_vcf");
		system("rm $poop_vcf");

		# --------------------------------------------------------------------------------

		# Do comparison again using downsampled samples

		print STDERR "." x 80;
		print STDERR "\n";

		my $poop_to_compare;
		my $blood_to_compare;

		# IDs later suffixed with "-samp" if necessary
		my $poop_id_s  = $poop_id;
		my $blood_id_s = $blood_id;

		# See if downsampled blood exists for this pair
		my $ds_fecal_id = ($NGS_ID_hash{$poop_id} =~ /(.+)/)[0];
		my $ds_blood = $NGS_ID_hash{$blood_id} . "_samp-" . $ds_fecal_id;
		my $ds_blood_file = $prefix . $ds_blood . $suffix;
		print STDERR "Looking for downsampled blood: [$ds_blood_file]\n";
		if (-e $ds_blood_file) {
			$blood_to_compare = $ds_blood;
			if ($blood_to_compare !~ /PE$/) {
				$blood_to_compare .= ".PE";
			}
			$poop_to_compare = $NGS_ID_hash{$poop_id} . ".PE";
			$blood_id_s .= "-samp";
		}

		# See if downsampled feces exists for this pair
		my $ds_blood_id = ($NGS_ID_hash{$blood_id} =~ /(.+)/)[0];
		my $ds_feces = $NGS_ID_hash{$poop_id} . "_samp-" . $ds_blood_id;
		my $ds_feces_file = $prefix . $ds_feces . $suffix;
		print STDERR "Looking for downsampled feces: [$ds_feces_file]\n";
		if (-e $ds_feces_file) {
			print STDERR "\tFound\n";
			$poop_to_compare = $ds_feces;
			if ($poop_to_compare !~ /PE$/) {
				$poop_to_compare .= ".PE";
			}
			$blood_to_compare = $NGS_ID_hash{$blood_id} . ".PE";
			$poop_id_s .= "-samp";
		}

		# Die if neither downsampled blood nor downsample feces exists,
		# because that means something bad happened with downsampling script.
		if ((! -e $ds_feces_file) && (! -e $ds_blood_file)) {
			print STDERR "ERROR: For this pair, neither downsampled blood nor ";
			print STDERR "downsampled feces exists. Aborting.";
			exit;
		}

		print STDERR "\t\tComparing blood $blood_to_compare to poop $poop_to_compare.\n";

		# Pull out individuals to make temporary individual VCF files
		$tmp_vcf_cmd_1 = "";
		$tmp_vcf_cmd_2 = "";
		$tmp_vcf_cmd_1 = $vcf_prefix . $blood_to_compare . " --out tmp${rand}.blood";
		$tmp_vcf_cmd_2 = $vcf_prefix . $poop_to_compare  . " --out tmp${rand}.feces";

		print STDERR "CMD 5: [" . $tmp_vcf_cmd_1 . "]\n";
		system($tmp_vcf_cmd_1);
		print STDERR "CMD 6: [" . $tmp_vcf_cmd_2 . "]\n";
		system($tmp_vcf_cmd_2);

		# Compare individual VCF files

		# Write individual mapping file
		open (MAP, ">tmp${rand}.map")
			or die "ERROR: Could not open mapping file. $!\n";
		print MAP $poop_to_compare . "\t" . $blood_to_compare . "\n";
		close MAP;

		$blood_vcf = "";
		$poop_vcf  = "";
		$blood_vcf = "tmp${rand}.blood.recode.vcf";
		$poop_vcf  = "tmp${rand}.feces.recode.vcf";

		$vcf_cmd_base = "";
		$vcf_cmd_base = "vcftools --vcf $blood_vcf --diff $poop_vcf ";
		$vcf_cmd_base .= "--diff-indv-map tmp${rand}.map ";
		$vcf_cmd_base .= "--out results/diff_${blood_id_s}_${poop_id_s}${out_suffix} ";

		# This fails. Maybe only in later versions of VCFtools
		# Sites that are common / unique to each file. *.diff.sites_in_files
		#my $vcf_cmd_site = $vcf_cmd_base . "--diff-site";

		# Discordance on a site by site basis. *.diff.sites
		$vcf_cmd_site_discord = "";
		$vcf_cmd_site_discord = $vcf_cmd_base . "--diff-site-discordance";

		# Discordance matrix. *.diff.discordance.matrix
		$vcf_cmd_discord_matrix = "";
		$vcf_cmd_discord_matrix = $vcf_cmd_base . "--diff-discordance-matrix";

		print STDERR "CMD 7: [" . $vcf_cmd_site_discord . "]\n";
		system($vcf_cmd_site_discord);
		print STDERR "CMD 8: [" . $vcf_cmd_discord_matrix . "]\n";
		system($vcf_cmd_discord_matrix);

		# Delete temporary VCF files
		system("rm $blood_vcf");
		system("rm $poop_vcf");
		system("rm tmp${rand}.map");
	}
}

exit;
