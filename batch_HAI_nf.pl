#!/usr/bin/perl

use File::Basename;

use strict;
use warnings;

# --------------------------------------------------------------------------- #

# File:                 batch_HAI_nf.pl
# Date created:         06 June 2022
# Date last modified:   06 June 2022
# Author:               Eliot Stanton (eliot.stanton@state.mn.us)
# Description:          Create list of data analysis deliverables from
#                       processed sequencing files and submit jobs to Slurm.

# --------------------------------------------------------------------------- #

# Define help string to be printed to command-line:
my $var_help    = "batch_HAI.pl [DIRECTORY_IN] [DIRECTORY_OUT] [EMAIL]

Batch submits the pairs of fastq files formatted YYYYEL-#####_R1.fastq.gz
[DIRECTORY IN] Directory containing pairs of FASTQ files for processing
[DIRECTORY OUT] Directory to contain output subdirectories
[EMAIL] Email address for notification (optional)
 ";

# --------------------------------------------------------------------------- #

# Variables passed to this script:
my $dir_in              = $ARGV[0];
my $dir_out             = $ARGV[1];
my $var_email		= $ARGV[2] || "";

# Print $var_help and end script if incorrect number of arguments:
die "$var_help\n" unless scalar @ARGV == 2 || scalar @ARGV == 3;

# Print $var_help and end script if $dir_in is missing:
die "$var_help \nNo input directory found!\n" unless -d $dir_in;

# Ensure $dir_in is accessible for other users:
system ( "chmod 770 $dir_in" );

# Make initial output directory if it doesn't already exist:
unless ( -d $dir_out ) {

	system ( "mkdir $dir_out" );
	system ( "chmod 770 $dir_out" );

}

# Data structures used in this script:
my @array_in;
my @array_out;

# Variables used by this script:
my $var_submit          = dirname(__FILE__);

# Add complete path for slurm submission script:
$var_submit             = `realpath $var_submit`;
chomp $var_submit;

# Add filepaths for submission slurm scripts:
my $var_submit_sh       = "$var_submit\/submit_HAI_nf.sh";

# --------------------------------------------------------------------------- #

print "\tsbatch \\
	-p small \\
	$var_submit_sh \\
	$dir_in \\
	$dir_out\n";

	system ("sbatch --mail-type=END --mail-type=FAIL --mail-user=$var_email -p small $var_submit_sh $dir_in $dir_out");

