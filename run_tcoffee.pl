#!/usr/bin/perl

use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $inputdir;

$inputdir = &overrideDefault("inputfile.dir",'inputdir');

######################################################################
# CODE
######################################################################

opendir(DIR, $inputdir) || die "Can't open directory\n";
my @array = ();
@array = readdir(DIR);

foreach my $file (@array){
	next unless ($file =~ /^\S+.fa$/);
	system ("t_coffee $inputdir/$file");

}
closedir (DIR);

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inputdir|i:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}


=head1 DESCRIPTION

	Run Multiple sequence alignment using T_Coffee in batch mode.

=head1 SYNOPSIS

script.pl  -i [-h]

 [-help -h]                Displays this basic usage information
 [-inputdir -i]            Input directory containing sequence files in fasta format 

