#Perl -Sx "{0}"; Exit {Status}
#!perl

# using a user-defined score threshold, finds clusters of samples exceeding scores
# specifically created to find likely sibling clusters
# requires as input: pairwise scores in list form

#changes in v2:
# added function printClustersExceptOne2File
# added ability to pass in thresholds from command line
# added print statements to state run parameters at beginning of program


use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
# use diagnostics;

####Defaults set here ###
my $score_file = "fake.txt";
my $threshold = "fake.txt";					#this is never used
my $prim_threshold = 0.45;
my $sec_threshold = 0.4;
my $ter_threshold = 0.25;
#########################

GetOptions ("in=s" => \$score_file,
			"threshold1:s" => \$prim_threshold,
			"threshold2:s" => \$sec_threshold,
			"threshold3:s" => \$ter_threshold
			) 
or die("Error in command line arguments\nUsage:\nperl findSibClusters.pl --in [score_file:string] [ --threshold1 [threshold:double]] [ --threshold2 [threshold:double]] [ --threshold3 [threshold:double]]\n");

my $clusters = {};

print STDERR "Using $score_file as input\n";
print STDERR "$prim_threshold is primary threshold\n";
print STDERR "$sec_threshold is secondary threshold\n";
print STDERR "$ter_threshold is tertiary threshold\n";
my $samples = readScoreFile ($score_file);
findClusters($samples, $clusters, $prim_threshold);
getClusterSampleScores($clusters, $samples);
my $subclusters1 = matchCluster($clusters, $samples, $sec_threshold, 1);
my $subclusters2 = matchCluster($clusters, $samples, $ter_threshold, $sec_threshold);
printClustersSubclusters ($clusters, $subclusters1, $subclusters2);
printClustersExceptOne2File ($clusters, $score_file.".sibclst");

#printClustersExceptOne ($clusters);
exit;

# prints all members of a cluster
sub printClustersSubclusters { my ($clusters, @subclusters) = @_;
	foreach my $id (keys %$clusters) {
		if ($id eq $clusters->{$id}->{'id'}) {
			foreach my $sample (keys %{$clusters->{$id}->{'members'}}) {
				print STDOUT "1\t$sample\t$clusters->{$sample}->{'id'}\t$clusters->{$id}->{'members'}->{$sample}\n";
			}
			for (my $i = 0; $i <= $#subclusters; $i++) { my $subcluster = $subclusters[$i];
				my $level = $i + 2;
				if ($subcluster->{$id}) {	#subcluster has some matches to sibling cluster
					foreach my $sample (keys %{$subcluster->{$id}}) {
						print STDOUT "$level\t$sample\t$id\t$subcluster->{$id}->{$sample}\n";
					}
				}
			}
		}
	}
}

# finds partially matching records to each cluster
sub matchCluster { my ($clusters, $samples, $thresh, $prev_thresh) = @_; my $subclusters = {};
	foreach my $id (keys %$clusters) {
		if ($id ne $clusters->{$id}->{'id'}) { next; }
		foreach my $sample (keys %$samples) {
			if ( $clusters->{$id}->{'members'}->{$sample}) { next; }	#if sample is part of the cluster don't bother
			my $average = getAverageMatch($samples, $sample, keys %{$clusters->{$id}->{'members'}});
			if ($average >= $thresh && $average <= $prev_thresh) {
				if (!$subclusters->{$id}) { $subclusters->{$id} = {}; }
				my $subclustinfo = $subclusters->{$id};
				$subclustinfo->{$sample} = $average;
			}
		}
	}
	return $subclusters;
}

# finds how well each sample in a cluster matches to it
sub getClusterSampleScores { my ($clusters, $samples) = @_; 
	foreach my $id (keys %$clusters) {
		if ($id ne $clusters->{$id}->{'id'}) { next; }
		foreach my $sample (keys %{$clusters->{$id}->{'members'}}) {
			my $average = getAverageMatch($samples, $sample, keys %{$clusters->{$id}->{'members'}});
			$clusters->{$id}->{'members'}->{$sample} = $average;
		}
	}
}

sub getAverageMatch { my ($samples, $sample1, @sample2s) = @_; my $sum = 0; my $count = 0; 
	if ($#sample2s == 0 && $sample1 eq $sample2s[0]) { return 1; }	#this is kind of a stupid workaround so that 1 member clusters can calculate their score with no problem
	foreach my $sample2 (@sample2s) {
		if ($sample1 eq $sample2) { next; }
		else { $sum += ($samples->{$sample1}->{$sample2}+$samples->{$sample2}->{$sample1})/2; $count++ }
	}
	return ($sum/$count);
}

# prints all members of a cluster
sub printClusters { my $clusters = shift;
	foreach my $id (keys %$clusters) {
		if ($id eq $clusters->{$id}->{'id'}) {
			if ($#{$clusters->{$id}->{'members'}} > 0) {
				foreach my $sample (@{$clusters->{$id}->{'members'}}) {
					print STDOUT "$sample\t$clusters->{$sample}->{'id'}\n"; 
				}
			}
		}
	}
}

# print all members of a cluster except one to file named after scorefile
sub printClustersExceptOne2File { my $clusters = shift; my $file = shift;
	open (my $fp, '>', $file) or die "Couldn't open $file for writing\n";
	print STDERR "Opened $file to print all members of each primary cluster but one\n";
	foreach my $id (keys %$clusters) {
		if ($id eq $clusters->{$id}->{'id'}) {
			foreach my $sample (keys %{$clusters->{$id}->{'members'}}) {
				if ($sample eq $id) {}
				else { print $fp "$sample\n"; }
			}
		}
	}
	close ($fp);
}

# print all members of a cluster except one
sub printClustersExceptOne { my $clusters = shift;
	foreach my $id (keys %$clusters) {
		if ($id eq $clusters->{$id}->{'id'}) {
			foreach my $sample (@{$clusters->{$id}->{'members'}}) {
				if ($sample eq $id) {}
				else { print STDOUT "$sample\n"; }
			}
		}
	}
}

sub findClusters { my $samples = shift; my $clusters = shift; my $threshold = shift;
	my @sample_names = sort keys %$samples;
	for (my $i = 0; $i < $#sample_names; $i++) {
		my $sample1 = $sample_names[$i];
		if (!$clusters->{$sample1}) { $clusters->{$sample1} = {}; 
			$clusters->{$sample1}->{'id'} = $sample1;
			$clusters->{$sample1}->{'members'} = {};
		} else { next; }
		my $cluster_info = $clusters->{$sample1};
		$cluster_info->{'members'}->{$sample1} = 1;
		for(my $j = $i + 1; $j <= $#sample_names; $j++) {
			my $sample2 = $sample_names[$j];
			my $score = ($samples->{$sample1}->{$sample2} + $samples->{$sample2}->{$sample1})/2;
			if ($score >= $threshold) {
				add2Cluster($clusters, $cluster_info, $sample2, $samples, $threshold)
			} 
		}
	}
}

sub add2Cluster { my ($clusters, $cluster_info, $sample, $samples, $threshold) = @_;
	if (!$clusters->{$sample}) {
		$cluster_info->{'members'}->{$sample} = 1;
		$clusters->{$sample} = $cluster_info;
		foreach my $next_sample (keys %$samples) {
			if ($next_sample eq $sample) { next; }
			my $score = ($samples->{$sample}->{$next_sample} + $samples->{$next_sample}->{$sample})/2;
			if ($score >= $threshold) {
				add2Cluster($clusters, $cluster_info, $next_sample, $samples, $threshold);
			}
		}
	}
	
}

sub readScoreFile { my $file = shift; my $samples = {};
	open (my $fp, '<', $file) or die "Couldn't open $file for reading\n";
	while (<$fp>) { chomp;
		my ($sample1, $sample2, @scores) = split(' ', $_);
		if ($sample1 eq $sample2) { next; }
		if (!$samples->{$sample1}) {$samples->{$sample1} = {}; }
		$samples->{$sample1}->{$sample2} = getMean(@scores);
	}
	close ($fp);
	return ($samples);
}

sub getMean { my @vals = @_;
	my $sum = 0; foreach my $val (@vals) {$sum += $val;}
	return ( $sum/($#vals+1) );
}
