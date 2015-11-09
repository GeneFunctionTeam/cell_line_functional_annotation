#!/usr/bin/perl -w
# ========================================= #
# process_WTSIcnv_data.pl
# read CNV data from WTSI, count mutations
# by type for each gene and cell line.
# Forked from process_mutation_data.pl on
# 5th Jan 2015 to simplify code management
# jamesc@icr.ac.uk, 25th July 2014
# ========================================= #

use strict;
use Getopt::Long;

my ($help, $wtsi_cnv_data, $genes, $cell_lines, $output_genes_list, $mut_cons_classes, $output);

GetOptions (
  "wtsi_cnv_data=s" => \$wtsi_cnv_data,
  "genes=s" => \$genes,
  "cell_lines=s" => \$cell_lines,
  "output_genes_list=s" => \$output_genes_list,
  "mut_consequence_classes=s" => \$mut_cons_classes,
  "output=s" => \$output,
  "help" => \$help,
 );

# print usage message if requested
if(defined($help)) {
  &usage;
  exit(0);
}

$output = "WTSI_CNV_data_table.txt" unless defined $output;

# ================================================================ #
# this will be a hash of hashes to store the mutation details for
# all gene_celllines in each data set. Primary keys are cell_line
# and gene. Secondary keys are data sources (eg. comsic_mut)
# ================================================================ #
 my %details_table = ();
 
 
# ============================== #
# read cell line name dictionary #
# ============================== #

open CL, "< $cell_lines" or die "Can't read cell line name dictionary $cell_lines: $!\n";
my %cell_lines;
while(<CL>){
  next if /^#/;
  my ($cell_line, $standard_cell_line, $source) = split(/\t/);
  $cell_lines{$cell_line} = $standard_cell_line;
}
close CL;


# ========================= #
# read gene name dictionary #
# ========================= #

open GN, "< $genes" or die "Can't read gene name dictionary $genes: $!\n";
my %genes;
while(<GN>){
  next if /^#/;
  my ($gene_alias, $standard_gene_name) = split(/\t/);
  chomp($standard_gene_name);
  $genes{$gene_alias} = $standard_gene_name;
}
close GN;


# ==================================== #
# read in the list of genes we want to #
# output at the end and standardise    #
# ==================================== #

# we have three identifiers for each gene
# 	the gene symbol (e.g. EIF1AX)
# 	the EntrezGene ID (e.g. 101060318)
# 	the ensembl gene ID(s) (e.g. ENSG00000173674_ENSG00000198692)
# this is followed by the cancer gene classification
# (TSG or OG) at the end.
# These need to be made into hashes so that any of the ID types
# can be used to look up the full set of identifiers...

open OUTGENES, "< $output_genes_list" or die "Can't read list of genes to output from file $output_genes_list: $!\n";
my %output_genes;
my %symbol_to_output_genes;
my %entrez_to_output_genes;
my %ensembl_to_output_genes;
while(<OUTGENES>){
  next if /^#/;
  next if /^[\r\n]/; # skip blank lines
  my ($symbol, $entrez, $ensembl, $type) = split /\t/;
  chomp($type);
  my @ensembls = split(/_/, $ensembl);
  $output_genes{"$symbol\t$entrez\t$ensembl"} = $type;
  $symbol_to_output_genes{$symbol} = "$symbol\t$entrez\t$ensembl";
  $entrez_to_output_genes{$entrez} = "$symbol\t$entrez\t$ensembl";
  foreach my $this_ensembl (@ensembls){
    $ensembl_to_output_genes{$this_ensembl} = "$symbol\t$entrez\t$ensembl";
  }
}
close OUTGENES;


# ============================================== #
# read the mutation consequences classifications #
# ============================================== #

open MCC, "< $mut_cons_classes" or die "Can't read mutation consequence classifications $mut_cons_classes: $!\n";
my %mutation_consequences;
while(<MCC>){
  next if /^#/;
  my ($consequence, $class, $source) = split(/\t/);
  chomp($source);
  $mutation_consequences{$consequence} = $class;
}
close MCC;


# ====================================== #
# Process the WTSI mutation and CNA data #
# ====================================== #


my %master_cell_lines_seen;		# store all cell lines seen
my %master_genes_seen;			# store all genes seen

my %wtsi_cnas;		# store: none/amp/del/gain/loss

if(defined $wtsi_cnv_data){
	# Open and read mutations file
	open WTSI, "< $wtsi_cnv_data" or die "Can't read mutations file $wtsi_cnv_data: $!\n";
	
	my $header = <WTSI>;
	
	while(<WTSI>){
	
	  # get gene, cell line, prot. mut., consequence etc.
	  # [0] ID_VARIANT
	  # [1] SAMPLE_NAME
	  # [2] CHR
	  # [3] GENOME_START
	  # [4] GENOME_STOP
	  # [5] STRAND
	  # [6] ALGORITHM
	  # [7] GENE_NAME
	  # [8] COSMIC_Name
	  # [9] TRANSCRIPT
	  # [10] WT
	  # [11] MT
	  # [12] CDS_SYNTAX
	  # [13] AA_MUT_SYNTAX
	  # [14] DESCRIPTION
	  # [15] QUALITY
	  # [16] ZYGOSITY
	  # [17] VERIF_STATUS_MINT
	  # [18] FATHMM Score
	  # [19] FATHMM Summary
	  # [20] Comment
	
	  my @fields = split(/\t/);
	  my $entrez_gene = $fields[7];  # this is actually a gene symbol
	  my $var_class = $fields[14];
	  my $sample = $fields[1];
	  my $genome_change = $fields[12];
	  my $prot_change = $fields[13];
	  my $mutation_zygosity = $fields[16];

	  # sometimes var_class is blank and we can do nothing with it... skip
	  next if $var_class eq '' || $var_class eq 'Unknown';
	  
	  # lookup gene and cell line name in dictionaries
	
	  $entrez_gene =~ s/_.+//;	# The COSMIC data often has the gene ID concatenated to the 
								# transcript ID. e.g. CEACAM18_ENST00000451626
								# strip out anything after '_' from entrez_gene
	
	  my $standard_gene = $entrez_gene;
	  
	  if(exists($genes{$entrez_gene})){
		$standard_gene = $genes{$entrez_gene};
	  }
	  $entrez_gene = $standard_gene ; # this may look nuts but we are going to munge standard_gene later and want to still keep the gene symbol handy...

	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
	  }
	
	  
	  # skip processing this variant if the gene symbol is not in the set of output genes
	  next unless exists $symbol_to_output_genes{$standard_gene};

	
	  # change $standard_gene from the symbol to the full trio of IDs and keep a record of
	  # the standardised gene symbol for later use in $entrez_gene;
	  $standard_gene = $symbol_to_output_genes{$standard_gene};
	
	  my $wtsi_key = "$standard_cell_line\t$standard_gene";

	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tCNA";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tCNA";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	  my $prot_position = $prot_change;
	  $prot_position =~ s/^p\.[A-Z]+//;
	  $prot_position =~ s/[^\d].*//;
	  my $davoli_key = "$entrez_gene\t$prot_position";

	  if($mutation_consequences{$var_class} eq "deletion"){
		$wtsi_cnas{$wtsi_key} = -2; # hom-del
	  }
	  elsif($mutation_consequences{$var_class} eq "amplification"){
		$wtsi_cnas{$wtsi_key} = 2; # amp
	  }

	  # collect the details for later output:
	  if(exists $details_table{$wtsi_key}{"wtsi_cna"}){
		$details_table{$wtsi_key}{"wtsi_cna"} .=  "; $genome_change";
		$details_table{$wtsi_key}{"wtsi_type"} .=  "; $mutation_consequences{$var_class}";
		$details_table{$wtsi_key}{"wtsi_davoli_freq"} .= "; NA";
	  }
	  else{
		$details_table{$wtsi_key}{"wtsi_cna"} =  $genome_change;
		$details_table{$wtsi_key}{"wtsi_type"} .=  $mutation_consequences{$var_class};
		$details_table{$wtsi_key}{"wtsi_davoli_freq"} = "NA";
	  }
	  
	} # finished reading file, close
	close WTSI;
}



# ==================================== #
# Process the hashes to format outputs
# ==================================== #

# get lists of all genes and cell lines seen in any data set
my @cell_lines_seen = keys %master_cell_lines_seen;

#my @genes_seen = keys %master_genes_seen;
my @genes_seen = keys %output_genes;		# the list of genes from CGC and C5000


# We need to indicate which data sets were available (mut, cna or expression) as a column
# in the output.

my $mutation_matrix = '';			# the functionally relevant changes
my $other_matrix = '';				# all changes regardless of relevance
my $mutation_classification = '';	# del=5,amp=4,trunc_hom=3,trunc_het=2,miss=1,overexpr=6,underexpr=7,wt=0

my $matrix_header = "cell_line";
foreach my $seen_gene (@genes_seen){	# The CGC and C5000s sets
  my $seen_gene_header = $seen_gene;
  $seen_gene_header =~ s/\t/_/g;
  $matrix_header .= "\t$seen_gene_header";
}
$matrix_header .= "\n";


while(my ($seen_cell_line, $cell_line_seen_in_dataset) = each  %master_cell_lines_seen){
  my $datasets = "";
  if($cell_line_seen_in_dataset =~ /mut/){
    $datasets .= 'mut_'
  }
  if($cell_line_seen_in_dataset =~ /CNA/){
    $datasets .= 'CNA_'
  }
  if($cell_line_seen_in_dataset =~ /exprn/){
    $datasets .= 'expression_'
  }
  
  $mutation_matrix .= "$seen_cell_line";
  $other_matrix .= "$seen_cell_line";
  $mutation_classification .= "$seen_cell_line";
  
  foreach my $seen_gene (@genes_seen){
    
    my $hash_key = "$seen_cell_line\t$seen_gene";
    my $hom_del = 0;
    my $loss = 0;
    my $gain = 0;
    my $amp = 0;
    my $trunc = 0;
    my $rec_mis = 0;
    my $other = 0;
    my $exprn = 0;
    
    if(exists($wtsi_cnas{$hash_key})){
      if($wtsi_cnas{$hash_key} == -2){
        $hom_del = 1;
      }
      elsif($wtsi_cnas{$hash_key} == 1){
        $gain = 1;
      }
      elsif($wtsi_cnas{$hash_key} == 2){
        $amp = 1;
        
#
print "found an amp in $seen_gene\n";
        
      }
      elsif($wtsi_cnas{$hash_key} == -1){
        $loss = 1;
      }
    }

    
    # now decide how to classify the cell line / gene...
    my $matrix_value = 0;

    if($output_genes{$seen_gene} eq 'OG'){
      if($amp == 1){
        $matrix_value = 1;
      }
#      elsif($hom_del == 1){
#        $other = 1;
#      }
      #elsif($hom_del == 1 || $loss == 1 || $gain == 1){
      #  $other = 1;
      #}
    }
    elsif($output_genes{$seen_gene} eq 'TSG'){
      if($hom_del == 1){
        $matrix_value = 1;
      }
#      elsif($amp == 1){
#        $other = 1;
#      }
      #elsif($amp == 1 || $loss == 1 || $gain == 1){
      #  $other = 1;
      #}
    }
    
    # add to the mutation/other_matrix depending on the classification
    if($matrix_value == 0){
      $mutation_matrix .= "\t0";
      if($other == 1){
        $other_matrix .= "\t1";
      }
      else{
        $other_matrix .= "\t0";
      }
    }
    else{
      $mutation_matrix .= "\t1";
      $other_matrix .= "\t1";
    }
    
    
    # values for each type of mutation
    # 5: hom del in TSG
    # 4: amp in OG
    # 3: hom trunc in TSG
    # 2: het trunc
    # 1: rec missense
    # 7: underexpressed
    # 6: overexpressed
    # 0: none
    
    # note that this only reports mutations considered to be functionally
    # relevant. We may also want to include non-functionally relevant
    # mutations for a separate plot. This applied especially to CNAs
    
    if($output_genes{$seen_gene} eq 'OG' && $amp == 1){
      $mutation_classification .= "\t4";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $hom_del == 1){
      $mutation_classification .= "\t5";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $loss == 1 && $trunc == 1){
      $mutation_classification .= "\t3";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $trunc == 1){
      $mutation_classification .= "\t2";
    }
    elsif($rec_mis == 1){
      $mutation_classification .= "\t1";
    }
    elsif($output_genes{$seen_gene} eq 'OG' && $exprn == 1){
      $mutation_classification .= "\t6";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $exprn == -1){
      $mutation_classification .= "\t7";
    }
    else{
      $mutation_classification .= "\t0";
    }
    
  }
  
  $mutation_matrix .= "\n";
  $other_matrix .= "\n";
  $mutation_classification .= "\n";
  
}

#
# details table to write
# cell line
# gene
# data source
# consequence (amp/del/missense/trunc/over or under expressed)
# mutation (residue change, expression value etc)
# davoli_freq if missense

open DETAILS, "> $output.details.txt" or die "Can't write details to $output:$! \n";

#print DETAILS "cell_line\tgene_symbol\tEntrez_geneID\tEnsembl_GeneID\tCOSMIC_mutation\tCOSMIC_consequence\tCOSMIC_Davoli\tICR_mutation\tICR_consequence\tICR_Davoli\tBIAKN_mutation\tBIAKN_consequence\tBIANK_Davoli\tWTSI_mutation\tWTSI_consequence\tWTSI_Davoli\tCCLE_CNA\tExprnZ\n";
print DETAILS "cell_line\tgene_symbol\tEntrez_geneID\tEnsembl_GeneID\tdata set\tconsequence\tmutation\tfrequency (if SNV)\n";

my @details_keys = keys %details_table;
foreach my $details_key (@details_keys){
  print DETAILS "$details_key";
  print DETAILS "\tWTSI (CNV)";


# ccle_cna
  if(exists $details_table{$details_key}{"wtsi_cna"}){
    chomp($details_table{$details_key}{"wtsi_cna"});
    
    print "$details_table{$details_key}{wtsi_cna}\n";
    
    if($details_table{$details_key}{"wtsi_type"} eq "deletion"){
      print DETAILS "\thom_del";
    }
    elsif($details_table{$details_key}{"wtsi_type"} eq "amplification"){
      print DETAILS "\tamp";
    }
    
    print DETAILS "\tCN: $details_table{$details_key}{'wtsi_cna'}";
    print DETAILS "\tNA"; # no Davoli data for CNA
    
  }
  else{
    print DETAILS "\tNA\tNA\tNA";
  }

  print DETAILS "\n";
}

close DETAILS;

open MUTMAT, "> $output.functional_mutations.txt" or die "Can't write to functional output file $output: $!\n";
print MUTMAT "$matrix_header$mutation_matrix";
close MUTMAT;

open OTHERMAT, "> $output.nonfunctional_mutations.txt" or die "Can't write to non-functional output file $output: $!\n";
print OTHERMAT "$matrix_header$other_matrix";
close OTHERMAT;

open MUTCLASS, "> $output.mutation_classifications.txt" or die "Can't write to mutation classification output file $output: $!\n";
print MUTCLASS "$matrix_header$mutation_classification";
close MUTCLASS;


sub usage() {
  my $usage =<<END;

# ---------------------------------- #
#  process_WTSIcnv_data.pl
#  James Campbell (jamesc\@icr.ac.uk)
# ---------------------------------- #

Usage:
perl process_WTSIcnv_data.pl [options]

Options:
  --help                    Display this message and quit
  --wtsi_cna_data           Path to the WTSI CNV file [required]
  --genes                   Path to the gene name dictionary [required]
  --cell_lines              Path to the cell line name dictionary [required]
  --output_genes_list       Path to the list of genes requested [required]
  --mut_consequence_classes Path to mutation consequences resource [required]
  --output                  Path to output. [optional]
  
END

  print $usage;
}


