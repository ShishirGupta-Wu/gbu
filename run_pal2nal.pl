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

my $inputproteindir;
my $inputdnadir;
my $outputdir;

$inputproteindir = &overrideDefault("inputdir.aln",'inputproteindir');
$inputdnadir = &overrideDefault("inputdir.dna",'inputdnadir');
$outputdir = &overrideDefault("outputdir.codon",'outputdir');

######################################################################
# CODE
######################################################################

use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);

 system("mkdir $outputdir"); 

my %count = ();
opendir(DIR, $inputproteindir) || die "Can't open directory\n";
my @store_array = ();
@store_array = readdir(DIR);
my $name = '';

foreach my $file (@store_array) {
	next unless ($file =~ /^\S+\.aln$/);
	if ($file =~ /^(\S+)\.aln$/){
		$name = $1;
	}

  system("perl pal2nal.pl $inputproteindir/$file $inputdnadir/$name.fa -output fasta > $outputdir/$name.codon.aln")
}
 
######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inputproteindir|p:s", "inputdnadir|d:s", "outputdir|c:s");
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

	Convert multiple sequence alignments of proteins and the corresponding DNA (or mRNA) sequences into a codon-based DNA alignment.

=head1 SYNOPSIS

script.pl  -p -d -c [-h]

 [-help -h]                 Displays this basic usage information
 [-inputproteindir -p]      Input directory containing multiple sequence alignments of proteins 
 [-inputdnadir -d]          Input directory containing the corresponding DNA sequences
 [-outputdir -c]            A directory of output codon-based alignments
 
