#! /usr/bin/perl
###############################################################################
#
#    runPhyML.pl
#
#	 Reconstruct a maximum likelihood tree for each gene alignment using PhyML.
#
###############################################################################

use Bio::AlignIO;
use Bio::TreeIO;
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
my $outputdir;

$inputdir = &overrideDefault("inputfile.dir",'inputdir');
$outputdir = &overrideDefault("outputfile.dir",'outputdir');

system("mkdir $outputdir"); 

######################################################################
# CODE
######################################################################

opendir(DIR, "$inputdir") or die;
my @array = ();
@array = readdir(DIR);

foreach my $file (@array){
	next unless ($file =~ /^\S+.phylip$/);
	#my $out_file = $file . ".align";		
	
# system ("mpirun -np 80 /artemis/software/bin/phyml-mpi -i $inputdir/$file -d nt -q sequential -m GTR -f m -t e -v e -c 4");	
 system ("/artemis/software/bin/phyml -i $inputdir/$file -d nt -q sequential -m GTR -f m -t e -v e -c 4");       #-d data_type;  
	                                                                                      #-m HKY85| K80 | F81 | GTR | custom (nucleotide-based model)  substitution model
	                                                                                      #-f m : nucleotide sequences: the equilibrium base frequencies are estimated using maximum likelihood
	                                                                                      #-t e: get the maximum likelihood estimate. -t ts/tv_ratio (transition/transversion ratio)
	                                                                                      #-v e: get the maximum likelihood estimate for the proportion of invariable sites.
	                                                                                      #-a gamma: gamma is the value of the gamma shape parameter that is estimated in the maximum likelihood framwork.
	                                                                                      #-c nb_subst_cat: number of relative substituion rate categories.
	                                                                                      #-b int: default is SH-like brance supports alone	
	                                                                                      #-b 0: don't compute branch support                                                                                     
}
closedir (DIR);

my $sym = '';
opendir(DIR, "$inputdir") or die;

@array = readdir(DIR);

foreach my $ele (@array){
	next unless ($ele =~ /_phyml_tree.txt$/);
	($sym) = $ele =~ /^(\S+)_phyml_tree.txt$/;
	#warn "$sym\n";
	open (FH, "$inputdir/$ele")|| die "can't open file:$!\n";
	open (OUT, ">$outputdir/$sym.tree") || die "can't open file:outfile\n";
	while (<FH>){
		chomp;
	$_ =~ s/\)(\d+\.\d+):/\):/g;                                                     #
	print OUT "$_\n";
}
}
closedir (DIR);

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inputdir|i:s", "outputdir|o:s");
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

__DATA__

=head1 NAME

    runPhyML.pl
=head1 DESCRIPTION

	Reconstruct a maximum likelihood tree for each gene alignment using PhyML.

=head1 SYNOPSIS

script.pl -i -o [-h]

 [-help -h]                Displays this basic usage information
 [-inputdir -i]            Input directory containing codon alignments in phylip format
 [-outputdir -o]           Output directory 
 
=cut
