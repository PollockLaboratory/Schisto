#Perl -Sx "{0}"; Exit {Status}
#!perl

#This script was copied from /data/house_programs/schisto_ddRAD/schisto_analysis/scripts/var_sharing_sampling_rares_v02.pl. It has been modified on 29Apr2019 by JAS with comments to explain the program.

#This script originally only work properly if there was a variant call for every sample at every variant analyzed.
#It has been modidified on 27Oct2016 by JAS to be able to use sites with missing samples as well.
# This version has also been modified to print comparisons across rows and columns instead of just rows

#this version only does comparisons on samples that share a rare variant. 

use warnings;
use strict;
use Data::Dumper;

if ($#ARGV != 2) {die "\n\tUsage: perl vcfSNPs.pl [infile] [# generations] [max_MAF]\n\n"; }

my $infile = $ARGV[0];	# the infile is a vcf
my $gen_count = $ARGV[1];
my $max_maf = $ARGV[2];
my $outfile = $infile."_var_sharing_sampling_".$gen_count."gens_MAF".$max_maf.".txt";

#read in all vars and memorize genotypes for each sample, store var in array
my $samples = []; my $all_vars = []; my $var_count = 0;
open (my $fp, '<', $infile) or die "Couldn't open $infile for reading\n";
while (<$fp>) {chomp;
	if ($_ =~ m/^\#\#/) { next; }	# lines that start with '##' are not used				
	if ($_ =~ m/^\#/) { $samples = getSampleNames ( $samples, $_); next; }	# lines that start with '#' are headers for the rest of the file and contain samples names		
	my @data = split(' ', $_); 
	# these lines separate out fields from the vcf file. $gen_calls is an array ref holding genotpye calls for each sample
	my $chrom = shift(@data); my $pos = shift(@data); my $rsid = shift(@data);
	my $ref = shift(@data); my $alt = shift(@data); my $qual = shift(@data);
	my $filter = shift(@data); my $info = shift(@data); my $format = shift(@data);
	my $gen_calls = \@data; my $freq = get_freq ($info);

	if (isBiallelic($alt) == 1 && $freq <= $max_maf) {	#only use loci that have only two alleles and have a frequency less than or equal to the user-specified maf
		$var_count++;
		# storing locus information and genotypes in a hash
		my $var_info = {};
		$var_info->{'chrom'} = $chrom;
		$var_info->{'pos'} = $pos;
		$var_info->{'samples'} = {};
		my $sample_vars = $var_info->{'samples'};
		for (my $i = 0; $i <= $#$gen_calls; $i++) {
			my $sample_id = $samples->[$i]->{'id'};
			my $var_calls = $gen_calls->[$i];
			my $gen_sum = getGenSum ($var_calls);
			$var_info->{'samples'}->{$sample_id} = $gen_sum;
		}
		push (@$all_vars, $var_info);	#push hash into array containing all the loci used
	}
}
close ($fp);
# $var_count = 3000;
my $shares = {};	#this hash ref stores proportion of alleles shared for each pair
for (my $gen = 1; $gen <= $gen_count; $gen++) {
	print STDERR "\r\e[KGen: $gen";
	my $gen_scores = get_random_vars ($var_count, $all_vars );	#selects a subset of all variants for use in analysis
	for (my $i = 0; $i <= $#$samples-1; $i++) {	#starting a loop to run through all the samples
		my $sample1_info = $samples->[$i];
		my $sample1_id = $sample1_info->{'id'};
		for (my $j= $i; $j <= $#$samples; $j++) {	#now looping through every other sample
			my $sample2_info = $samples->[$j];
			my $sample2_id = $sample2_info->{'id'};		
			compare_genotypes ($gen_scores, $sample1_id, $sample2_id, $shares);
		}
	}
}

my $scores = {};		#stores mean and std deviation of proportion of alleles shared for each pair
for (my $i = 0; $i <= $#$samples; $i++) {
	my $sample1_info = $samples->[$i];
	my $sample1_id = $sample1_info->{'id'};
	for (my $j= $i; $j <= $#$samples; $j++) {
		my $sample2_info = $samples->[$j];
		my $sample2_id = $sample2_info->{'id'};
		my $share1_scores = $shares->{$sample1_id}->{$sample2_id};	# $share1_scores is an array containing $var_count proportions of alleles shared with the other sample in this pair
		my $share2_scores = $shares->{$sample2_id}->{$sample1_id};
		$scores->{$sample1_id}->{$sample2_id} = get_sample_stats ($share1_scores);
		$scores->{$sample2_id}->{$sample1_id} = get_sample_stats ($share2_scores);
	}
}

#print proportion of alleles shared in matrix format
# CURRENTLY ISN'T PRINTING STANDARD DEVIATIONS #
open (my $outfp, '>', $outfile) or die "Couldn't open $outfile for printing\n";
print STDERR "\nNow printing results\n";
for (my $i = 0; $i <= $#$samples; $i++) {
	my $sample1_info = $samples->[$i];
	my $sample1_id = $sample1_info->{'id'};
	print $outfp "\t$sample1_id";
}
print $outfp "\n";
for (my $i = 0; $i <= $#$samples; $i++) {
	my $sample1_info = $samples->[$i];
	my $sample1_id = $sample1_info->{'id'};
	print $outfp "$sample1_id\t";
	for (my $j= 0; $j <= $#$samples; $j++) {
		my $sample2_info = $samples->[$j];
		my $sample2_id = $sample2_info->{'id'};
		print $outfp "$scores->{$sample1_id}->{$sample2_id}->{'mean'}\t";
	}
	print $outfp "\n";
}
close ($outfp);

exit;

sub get_sample_stats { my ($share_scores) = @_;
	my $mean; my $stddv;
	if ($share_scores){
		$mean = getMean( $share_scores );
		$stddv = getStddev ( $share_scores, $mean );
	} else {
		$mean = "NaN";
		$stddv = "NaN";
	}
	return ( {'mean' => $mean, 'stdv' => $stddv} );
}

# counts the proportion of variants shared between two samples as follows: (# variants shared by the two samples)/(# of variants present in one sample)
sub compare_genotypes { my ($gen_scores, $sample1_id, $sample2_id, $shares) = @_;
	my $sample1_count = 0; my $sample2_count = 0; my $sample1_share = 0; my $sample2_share = 0;
	if ( $#{$gen_scores->{$sample1_id}} != $#{$gen_scores->{$sample2_id}} ) { die "Error: Unequal variant list lengths\n"; }
	for (my $i = 0; $i <= $#{$gen_scores->{$sample1_id}}; $i++ ) {		#looping through genotypes of sampled variants
		my $sample1_gen = $gen_scores->{$sample1_id}->[$i];				#gets genotype of first sample
		my $sample2_gen = $gen_scores->{$sample2_id}->[$i]; 			#gets genotype of second sample
		if ( $sample1_gen == -1 || $sample2_gen == -1 ) { next; }		# if either sample is missing a genotype at the locus, skip the locus	
		if ( $sample1_gen >= 1) {										# genotpyes are 0, 1, or 2. If the genotype is greater than 1, it means that it has a copy of the variant. Note that this means that we're only using sites where the sample has a rare variant- two samples might have the same genotype at a locus, but if it doesn't involve a rare variant we're not keeping track of that site
			$sample1_count++;
			if ($sample2_gen >= 1) { $sample1_share++; }				#if second sample also has genotype >= 1, it also has at least one copy of the variant. Count the number of variants shared between the two samples.
		}	
		if ( $sample2_gen >= 1) {
			$sample2_count++;
			if ($sample1_gen >= 1) { $sample2_share++; }
		}
	}
	# print STDERR "\n$sample1_count\t$sample1_share\n";
	store_share_results ($shares, $sample1_count, $sample1_share, $sample1_id, $sample2_id);
	# if ( $sample1_id eq $sample2_id ) { next; }
	store_share_results ($shares, $sample2_count, $sample2_share, $sample2_id, $sample1_id);
}

sub store_share_results { my ($shares, $count, $shared, $id1, $id2) = @_;
	if ( $count >0 ) {
		my $perc_shared = $shared/$count;
		push (@{$shares->{$id1}->{$id2}}, $perc_shared);
	}
}

# randomly samples variants with replacement.
sub get_random_vars { my ($var_count, $all_vars) = @_;
	my $gen_scores = {};
	for (my $i = 0; $i< $var_count; $i++) {
		my $index = int(rand(@{$all_vars}));			# choose an index from array containing all possible variants at random
		my $rand_var_info = $all_vars->[$index];
		my $sample_vars = $rand_var_info->{'samples'};	# get genotypes of randomly chosen locus
		foreach my $sample (keys %$sample_vars) {		# store genotypes of randomly chosen variant for each sample
			push ( @{$gen_scores->{$sample}}, $sample_vars->{$sample});		
		}
	}
	return $gen_scores;
}

sub getMean { my $numbers = shift;
	my $count = 0;
	my $sum = 0;
	foreach my $num (@$numbers) {
		$sum += $num;
		$count++;
	}
	return ($sum/$count);
}

# Calculate the Mean of your data set.
# For each number: subtract the Mean and square the result.
# Calculate the Mean of those squared differences. The result is called the Variance.
# Calculate the square root of the Variance. The result is the Standard Deviation.
# Important note- If the data set is only a sample of the whole population, change the formula to divide by one less when calculating the Variance. This is called Bessel's correction.
sub getStddev { my $numbers = shift; my $mean = shift;
	my $dev = 0; my $count = 0;
	foreach my $num (@$numbers) {
		$dev += ($num-$mean)**2;
		$count++;
	}
	my $variance = sqrt($dev);
	if ($count > 1) { $variance = $dev/($count-1); }
	return ( sqrt($variance) );
}

# checks to see if locus has more than 1 alternate allele
sub isBiallelic { my $alt = shift;
	my @alts = split(',', $alt);
	if ( $#alts > 0 ) { return 0; }
	else {return 1; }
}

# Returns the sum of genotype values at a locus. Ex: vcf genotype= 0/1 (heterozygous), then sum is 1. Missing loci are returned as -1
sub getGenSum { my $var_calls = shift;
	my @alleles = split(':', $var_calls);
	my @genotypes = split('\/', $alleles[0]);
	my $gen_sum;
	if (isGenoCalled($genotypes[0]) && isGenoCalled($genotypes[1]) ) {
		$gen_sum = $genotypes[0] + $genotypes[1];
	} else { $gen_sum = -1; }
	if ($gen_sum > 2) {die "Is this really bi_allelic?\n@alleles\n"; }
	return $gen_sum;
}

# checks to see of a genotype is called or missing 
sub isGenoCalled { my $number = shift;
	if ($number eq ".") {
		# print STDERR "$number is not a number\n"; 
		return 0;
	}
	else { 
		# print STDERR "$number is a number\n"; 
		return 1; 
	}
}

# returns an array containing names of all the samples in the vcf
sub getSampleNames { my $samples = shift; my $line = shift;
	my @fields = split(' ', $line);
	for (my $i = 9; $i <= $#fields; $i++) {
		my $sample_info = {};
		$sample_info->{'id'} = $fields[$i];
		push (@{$samples}, $sample_info);
	}
	return $samples;
}

sub get_freq { my $info = shift;
	# print STDERR "$info\n";
	# $info =~ m/.*;AF\=([0-9]\.[0-9]*).*/;
	$info =~ m/.*;AF\=([0-9]\.[0-9]*e-[0-9]*|0\.[0-9]*|1\.00).*/;
	if (!$1) { die "ERROR: $info\n in get_freq sub\n"; }
	return ($1);
}
