#!/usr/bin/perl -w
# extract and count gene*mutations from  the mutation dataset published by
# Davoli et al. (2012) and available from:
# http://elledgelab.med.harvard.edu/wp-content/uploads/2013/11/Mutation_Dataset.txt.zip
#
# Data set described at:
# http://www.cell.com/cell/supplemental/S0092-8674%2813%2901287-7
#
# This script updated 3rd Oct 2014 to include number of times a gene was sequenced
# and to also express position_mutation counts as a percentage of times sequenced. 
# different point mutations at a given position are counted irrespctive of the alt
# amino acid.


use strict;

my $gene_name_dictionary_file = "../resources/gene_name_dictionary_v03.txt";
#my $mutations_file = "../resources/Davoli_mutation_data/Mutation_Dataset.txt";
my $mutations_file = "/Users/jamesc/Documents/113_data_sources/Davoli_TCGA_COSMIC_Alexandrov_set/Mutation_Dataset.txt";
my $output_file = "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/resources/Davoli_mutation_data/Mutation_Dataset.txt.mut_counts.txt";


# format of the data is:
#Gene    Genome.position.hg19    Reference       Mutation        Protein_Change  Mutation_Type   Tumor_Sample    Tumor_Type
#A1BG    19:58861774-58861774    C       -       p.G385fs        Indel Frameshift        TCGA-18-3406-01A        Lung Squamous Cell Carcinoma
#A1BG    19:58864380-58864381    CC      -       p.G85fs Indel Frameshift        TCGA-G2-A3IE-01A        Bladder Carcinoma
#A1BG    19:58862784-58862784    C       T       p.A295T Missense        pfg019T Stomach Adenocarcinoma

# Note that the Gene column contains official HGNC symbols (eg KMT2B, not MLL2 or MLL4...)


# read in the gene name dictionary
my %gene_names;
open GN, "< $gene_name_dictionary_file" or die "unable to open gene name dictionary $gene_name_dictionary_file: $!\n";
while(<GN>){
  next if /^#/;
  my ($alias, $gene_name) = split(/\t/);
  chomp($gene_name);
  $gene_names{$alias} = $gene_name;
}
close GN;


open MUTS, "< $mutations_file" or die "Can't read mutations input file $mutations_file: $!\n";
my %gene_pos_mut_counts;
my %gene_counts;
while(<MUTS>){

  my @fields = split(/\t/);
  my $gene = $fields[0];
  my $prot_pos = $fields[4];
  my $effect = $fields[5];
  
  next unless $effect =~ /Missense/;

  my $position = undef;
  if($prot_pos =~ /^p\.[A-Z]+(\d+)/){
    $position = $1;
  }
  else{
    print "skipping $prot_pos as not a point mutation\n";
    next;
  }

  my $standard_gene = '';
  if(exists($gene_names{$gene})){
  	$standard_gene = $gene_names{$gene};
  }
  else{
  	warn "$gene not found in the gene name dictionary - using the original name.\n";
  	$standard_gene = $gene;
  }
  
#  $gene_pos_mut_counts{$standard_gene . "\t" . $prot_pos . "\t" . $effect} ++;
  $gene_pos_mut_counts{$standard_gene . "\t" . $position . "\t" . $effect} ++;
  
  if(exists($gene_counts{$standard_gene})){
    $gene_counts{$standard_gene} ++;
  }
  else{
    $gene_counts{$standard_gene} = 1;
  }
  
}
close MUTS;


open PROT_POS_COUNTS, "> $output_file" or die "Can't write protein position counts: $!\n";
while(my($prot_pos, $count) = each %gene_pos_mut_counts){
  next if $count < 2;
  
  my ($standard_gene, $prot_pos, $effect) = split(/\t/, $prot_pos);
  
  my $gene_count = $gene_counts{$standard_gene};
  print PROT_POS_COUNTS "$standard_gene\t$prot_pos\t$count\t$gene_count\n";
}
close PROT_POS_COUNTS;



