#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 6-22-21
# Last modified: 6-22-21
# Title: splitFastqcData.pl
# Purpose: Split fastqc_data.txt into sub-files based on modules
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $input;
my $outputDir = "";

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "input=s"           => \$input,
            "outputDir=s"       => \$outputDir
      ) or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $outputFH;
my $baseColumn;


##############################
# Code
##############################

if ($outputDir eq "") {
    $input =~ /^(.+\/)/;
    $outputDir = $1;
    mkdir $outputDir . "/split_data"
}

##############################
### Stuff
### More stuff

open my $inputFH, "$input" or die "Could not open first input\n$!";

while (my $input = <$inputFH>){
    chomp $input;
    if ($input =~ /^>>END_MODULE/) {
        # Close file
        close $outputFH;
        undef $baseColumn;
    } elsif ($input =~ /^>>/) {
        # Open new file
        $input =~ s/^>>//;
        $input =~ s/ /_/g;
        $input =~ s/\t.+//;

        open $outputFH, ">", $outputDir . "/split_data/" . $input . ".txt";
    } elsif ($input !~ /^##/) {
        # Deal with data
        ## Find column that has base numbering
        if ($input =~ /Base|Position/ && $input =~ /^#/) {
            # Change the word Position to Base to make downstream plotting easier
            $input =~ s/Position/Base/;
            my @dataArray = split "\t", $input;
            for (my $i = 0; $i < scalar(@dataArray); $i++) {
                if ($dataArray[$i] =~ /^#*Base$/) {
                    $baseColumn = $i;
                }
            }
        }

        # write to file
        if (defined($baseColumn)) {
            my @dataArray = split "\t", $input;
            if ($dataArray[$baseColumn] =~ /-/) {
                my ($lowerNum, $upperNum) = split "-", $dataArray[$baseColumn];
                for (my $i = $lowerNum; $i <= $upperNum; $i++) {
                    my @tempArray = @dataArray;
                    $tempArray[$baseColumn] = $i;
                    print $outputFH join("\t", @tempArray), "\n";
                }
            } else {
                print $outputFH $input, "\n";
            }
        } else {
            print $outputFH $input, "\n";
        }
    }
}


##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    xxxxxx.pl - generates a consensus for a specified gene in a specified taxa
    
Usage:

    perl xxxxxx.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
