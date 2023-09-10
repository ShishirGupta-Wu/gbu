#!/usr/bin/perl

open my $file_a, '<', 'HPIDB_interaction_file_nrd.txt_notab';
open my $file_b, '<', 'all_vaues_to_map';
open my $file_c, '>', 'primary_interactome_HP';
my %map;
foreach (<$file_b>) {
    chomp;
    my ($id,$dmn) = split;
    push(@{ $map{$id} },$dmn);
}
foreach (<$file_a>) {
    chomp;
    my ($srk,$mrs) = split;
    for my $dmn2 ( @{ $map{$mrs} } ) {
        for my $dmn1 ( @{ $map{$srk} } ) {
            print $file_c "$srk $dmn1 $dmn2 $mrs\n";
        }
    }
}
close $file_a;
close $file_b;
close $file_c;
