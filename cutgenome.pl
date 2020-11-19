#! /usr/bin/perl

# cutgenome - script to in-silico restriction digest an entire genome
# uses 2 restriction enzymes to cut the DNA and returns information about those pieces cut by both
# Written by Corey Cox 
# Usage:
#	cutgenome.pl cutsitefile genometocut seq_enzyme seq_len outputfile

# set up environment 
use strict;
use warnings;

#define global variables
my $cutsites = {};
my $cutlocs = {};
my $contig_order = [];

# get command line arguments and parse them to script variables
my $cutfile = $ARGV[0];
my $fasta = $ARGV[1];
my $seq_enzyme = $ARGV[2];
my $seq_len = $ARGV[3];
my $output = $ARGV[4];

if ($#ARGV != 4) {die "Usage: cutgenome.pl cutsitefile genometocut seq_enzyme seq_len outputfile\n"; }

readx($cutfile, \&getCutSites, "", $cutsites);
if (!$cutsites->{$seq_enzyme}) { die "\nCheck that the seq_enzyme is in cutfile\n\n"; }
readFasta($fasta, \&getCuts, "", $cutsites, $cutlocs, $contig_order);
printx($output, \&cutsToBed, "", $cutlocs, $contig_order);

sub getCutSites { my $fp = shift; my $cutsites = shift;
    while (<$fp>) { chomp;
		my ($name, $site ) = split (' ', $_); 
		$cutsites->{$name} = $site;
	}
}

sub getCuts { my $name = shift; my $seq = shift; my $cutsites = shift; my $cutlocs = shift; my $contig_order = shift;
    push(@{$contig_order}, $name); my $cuts = getKeyHash($cutlocs, $name);
    foreach my $sitename (keys %{$cutsites}) {
        my $site = $cutsites->{$sitename}; 
		my $cutpos = index($site, '^'); $site =~ s/\^//g;
        while ($seq =~ /$site/ig) { 
			my $loc = pos ($seq) - length ($site) + $cutpos; 
			if ($cuts->{$loc}) { $cuts->{$loc} = "Both"; }
			else { $cuts->{$loc} = $sitename; }
		}
    }
}

sub cutsToBed { my $fp = shift, my $cutsites = shift, my $contig_order = shift;
    foreach my $contig (@{$contig_order}) { my $prev_site = 0, my $prev_loc = 0;
        foreach my $loc (sort {$a <=> $b} keys %{$cutsites->{$contig}}){
            my $cur_site = $cutsites->{$contig}->{$loc};
			my $loc_len = $loc - $prev_loc;
			if (!$prev_site || !$prev_loc) {$prev_site = $cur_site; $prev_loc = $loc; next;}
			if (($prev_site eq $seq_enzyme || $prev_site eq "Both") && ($cur_site ne $seq_enzyme)) {
				if ($loc_len >= $seq_len) { 
					my $loc_stop = $prev_loc + $seq_len;
					print $fp "$contig\t$prev_loc\t$loc_stop\t$prev_site\t$cur_site\t$contig\t$prev_loc\t$loc\n";
				} else { print $fp "$contig\t$prev_loc\t$loc\t$prev_site\t$cur_site\t$contig\t$prev_loc\t$loc\n"; }
			} 
			if ( ($cur_site eq $seq_enzyme || $cur_site eq "Both") && ($prev_site ne $seq_enzyme) ) {
				if ($loc_len >= $seq_len) { 
					my $loc_start = $loc - $seq_len; 
					print $fp "$contig\t$loc_start\t$loc\t$prev_site\t$cur_site\t$contig\t$prev_loc\t$loc\n";
				} else { print $fp "$contig\t$prev_loc\t$loc\t$prev_site\t$cur_site\t$contig\t$prev_loc\t$loc\n"; }
			}
            $prev_site = $cur_site; $prev_loc = $loc;
        }
    }
}

sub readFasta { my $file = shift; my $sub = shift; my $blurb = shift; my @args = @_; my $name = ""; my $seq = "";
    open(my $fp, "<", $file) or die "readFasta cannot open $file for reading: $!";
    while (<$fp>) { if ($_ =~ s/>//) { chomp; $name = $_; }; last; } # get first sequence name
    while (<$fp>) { chomp;
        if ($_ =~ s/>//) { $sub->($name, $seq, @args); $name = $_; $seq = ""; next; }
        $seq .= uc($_);
    }
    $sub->($name, $seq, @args); $name = $_; $seq = "";
    close $fp;
}

sub readx { my ($file, $sub, $blurb, @args) = @_;
    if (open(my $fp, "<", $file)) { my $return = $sub->($fp, @args); close $fp; return $return; }
    else { return -1; }
}

sub printx { my $target = shift; my $sub = shift; my $blurb = shift; my @args = @_;
	my $file_handle; my $handle; 
    if (!$target) { $handle = \*STDERR; }
    if ($target eq \*STDERR || $target eq \*STDOUT) { $handle = $target; }
    else { open($file_handle, ">", $target); $handle = $file_handle; }
    if ($handle) { $sub->($handle, @args); }
    if ($file_handle) { close($file_handle); }
}

sub getKeyHash { my ($hash, $key) = @_; if (!$hash->{$key}) { $hash->{$key} = {}; } return $hash->{$key}; }
