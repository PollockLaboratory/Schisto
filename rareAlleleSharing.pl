#Perl -Sx "{0}"; Exit {Status}
#!perl

# computes genotype similarity for all pairwise sample combinations in a vcf
# requires as input: vcf

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
# use diagnostics;

my $vcf="fake.txt";
my $min_var;

####Defaults set here ###
my $num_vars = 500;
my $gens = 30;
my $max_freq = 0.1;
#########################

GetOptions ("vcf=s" => \$vcf,
			#"max_freq:f" => $max_freq,
			"num_gens:i" => \$gens,
			"num_vars:i" => \$num_vars) 
or die("Error in command line arguments\nUsage:\nperl rareAlleleSharing.pl --vcf [vcf] [ --num_gens [int] --num_vars [int] ]\n");


#### Global Variables ###

#########################
my $outfile = $vcf."_maxFreq".$max_freq."_numVars".$num_vars."_numGens".$gens.".rareAlleleSharing";

print STDERR "Running with $gens generations of $num_vars variants\n";

my ($vars, $samples) = getVcfVars($vcf, $max_freq);
if ($num_vars > scalar (keys %$vars)) {die "Maximum allowed number of vars to sample is".scalar (keys %$vars).".\n"; }
alleleSharingDistance ($vars, $samples, $num_vars, $gens);
printAlleleSharingDistanceMatrix ($samples, $outfile."Matrix");
printAlleleSharingDistancePairs ($samples, $outfile."Pairs");

exit;

sub printAlleleSharingDistancePairs { my ($samples, $file) = @_;
	open (my $fp, '>', $file) or die "couldn't open $file for printing\n";
	foreach my $sample1 (sort keys %$samples) {
		foreach my $sample2 (sort keys %$samples) {
			print $fp "$sample1\t$sample2\t";
			foreach my $score ( @{$samples->{$sample1}->{'sharing'}->{$sample2}} ) {
				print $fp "$score\t";
			} print $fp "\n";
		}
	}
	close($fp);
}

sub printAlleleSharingDistanceMatrix { my ($samples, $file) = @_;
	open (my $fp, '>', $file) or die "couldn't open $file for printing\n";
	print $fp "Sample\t";
	foreach my $sample1 (sort keys %$samples) {
		print $fp "$sample1\t";
	} print $fp "\n";
	foreach my $sample1 (sort keys %$samples) {
		print $fp "$sample1\t";
		foreach my $sample2 (sort keys %$samples) {
			my $mean  = getMean($samples->{$sample1}->{'sharing'}->{$sample2});
			print $fp "$mean\t";
		}
		print $fp "\n";
	}
	close($fp);
}

sub getMean {my $array = shift;
	my $sum = 0;
	foreach my $number (@$array) { $sum += $number; }
	return ($sum/($#$array+1));
}

# get pairwise allele sharing distance by randomly sampling called loci with replacement
sub alleleSharingDistance { my ($vars, $samples, $var_count, $gens) = @_;
	print STDERR "Entering subroutine alleleSharingDistance\n";
	foreach my $sample1 (sort keys %$samples) {
		my $sample1_info = $samples->{$sample1};
		my @var_list = keys %{$sample1_info->{'vars'}};
		if (!$sample1_info->{'sharing'}) { $sample1_info->{'sharing'} = {}; }
		foreach my $sample2 (sort keys %$samples) { if ($sample1_info->{'sharing'}->{$sample2}) { next; }
			print STDERR "\r\e[K$sample1 $sample2";
			$sample1_info->{'sharing'}->{$sample2} = []; my $sharing = $sample1_info->{'sharing'}->{$sample2};
			for (my $gen = 1; $gen <= $gens; $gen++) { 
				my $sample2_info = $samples->{$sample2}; my $scores = [];
				while ($#$scores < $var_count) {
					# print STDERR "\r\e[K$sample1 $sample2\tGen: $gen Var count: $#$scores";
					my $index = int(rand(@var_list)); my $var = $var_list[$index];
					my ($allele1, $allele2) = @{$sample1_info->{'vars'}->{$var}};
					my ($allele3, $allele4) = @{$sample2_info->{'vars'}->{$var}};
					if ( $allele1 eq "." || $allele2 eq "." || $allele3 eq "." || $allele4 eq "." || ($allele1 == 0 && $allele2 == 0) ) {} #do nothing
					else { 
						my $sum1 = $allele1 + $allele2; my $sum2 = $allele3 + $allele4;
						my $similarity;
						if ($sum1 == 0) {
							die "Sum1 is $sum1. The script should never come here.\n";
						} elsif ($sum1 ==1) {
							if ($sum2 == 0) { $similarity = 0; }
							elsif ($sum2 == 1) { $similarity = 1; }
							elsif ($sum2 == 2) { $similarity = 0.5; }
						} elsif ($sum1 == 2) {
							if ($sum2 == 0) { $similarity = 0; }
							elsif ($sum2 == 1) { $similarity = 0.5; }
							elsif ($sum2 == 2) { $similarity = 1; }
						} else {die "How is the sum of $allele1 and $allele2 $sum1?\n";}
						push (@{$scores}, $similarity);
					}
				}
				push(@$sharing, getMean($scores));
			}
		}
	}
	print STDERR "\n";
}

#save info and genotypes from each locus in vcf
sub getVcfVars { my $vcf = shift; $max_freq = shift; my $samples = {}; my $vcf_vars = {}; my $linecount = 0;
	my $sample_order = [];
	open (my $fp, '<', $vcf) or die "Couldn't open $vcf for reading\n";
	print STDERR "Opened $vcf for reading in variants\n";
	burnVcfHeader($fp, $sample_order);
	while (my $line = <$fp>) { chomp $line;
		my $var_info = getVarInfo($line, $sample_order, $samples);
		$var_info->{'lineID'} = $linecount; 
		my $var_id = $var_info->{'id'};
		my $freq = calcFreq ($var_info->{'genos'});
		if (!$vcf_vars->{$var_id} && ($freq <= $max_freq ) ) { # could also add in '|| $freq >= 1-$max_freq' but then 0 becomes the alternate allele and makes things difficult
			$vcf_vars->{$var_id} = $var_info; 
			foreach my $sample (keys %{$var_info->{'genos'}}) {						# set up $samples to point at $var_info
				if (!$samples->{$sample}) {$samples->{$sample} = {}; } my $sample_info = $samples->{$sample};
				if (!$sample_info->{'vars'}) {$sample_info->{'vars'} = {}; }
				$sample_info->{'vars'}->{$var_id} = $var_info->{'genos'}->{$sample};
			}
			$linecount++;
		} #else {print STDERR "Var info for $var_id is already stored\n"; } #or the frequency just wasn't in required range...
		
	}
	close($fp);
	print STDERR "Done reading variants\n";
	return ($vcf_vars, $samples);
}

sub calcFreq { my $genos = shift;
	my $var_count = 0; my $total_count = 0;
	foreach my $sample (keys %$genos) {
		my @calls = @{$genos->{$sample}};
		foreach my $call (@calls) {
			if ($call ne ".") { $total_count++; 
				if ($call == 1 ) { $var_count++; }
			}
		}
	}
	return ($var_count/$total_count);
}

# record the varinfo hash by reading a vcf line
sub getVarInfo { my ($var_line, $sample_names, $samples ) = @_;  my @var_data = split(' ', $var_line);
    my $var_info = {};
    $var_info->{'chrom'} = shift(@var_data); $var_info->{'pos'} = shift(@var_data);
    $var_info->{'id'} = shift(@var_data); $var_info->{'ref'} = shift(@var_data);
    $var_info->{'alt'} = shift(@var_data);
    # next line is just shifting off stuff we don't care about at the moment
    my $qual = shift(@var_data); my $filter = shift(@var_data); my $info = shift(@var_data); my $format = shift(@var_data);
    $var_info->{'genos'} = getVariantCalls (\@var_data, $sample_names );
	$var_info->{'id'} = $var_info->{'chrom'}.".".$var_info->{'pos'};
	return $var_info;
}

# returns an hash ref with variant calls for all samples at a locus
sub getVariantCalls { my $sample_genotypes = shift; my $sample_names = shift;
    my $variants = {};
	if($#$sample_genotypes ne $#$sample_names) { die "Different genotypes and names length\n"; }
    for ( my $i = 0; $i <= $#$sample_genotypes; $i++ ) {
        my $uncut = $sample_genotypes->[$i]; my $sample = $sample_names->[$i];
        my @geno_info = split(':', $uncut);   # split by genotype format fields
        my ($allele1, $allele2) = split('\/', $geno_info[0]); # split on unphased allele
		$variants->{$sample} = [];
        push (@{$variants->{$sample}}, $allele1, $allele2);
    }
    return $variants;
}

# burns through vcf header until final header line and passes final header line to a sub routine designated in $sub
sub burnVcfHeader { my $fp = shift; my $sample_order = shift; 
    while (my $line = <$fp>) {
        if ( $line =~ m/\#\#/) {  }
            elsif ( $line =~ m/\#/) { getSampleNames($line, $sample_order); return; } 
            else { die "Program should never reach this line in sub burnVcfHeader \n$line\n"; return 0; }
    }
}

# returns an array containing the names of all the samples in a vcf
sub getSampleNames { my $line = shift; my $names = shift; chomp $line; my $samplestart = 9;
    my @header = split(' ', $line);
    @{$names} = @header[$samplestart..$#header]; # double check that this works
}