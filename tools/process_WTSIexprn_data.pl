#!/usr/bin/perl -w
# ========================================= #
# process_WTSIexprn_data.pl
# read data from preprocessed Garnet
#  expression for each gene and cell line.
# Forked from process_mutation_data.pl on
# 5th Jan 2015 to simplify code management
# jamesc@icr.ac.uk, 25th July 2014
# ========================================= #

use strict;
use Getopt::Long;

my ($help, $wtsi_expression_z_data, $genes, $cell_lines, $output_genes_list, $mut_cons_classes, $output);

my $exprn_threshold = 3; # approximately p=6.33E-05

GetOptions (
  "wtsi_expression_z_data=s" => \$wtsi_expression_z_data,
  "exprn_threshold=i" => \$exprn_threshold,
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

$output = "wtsi_exprn_data_table.txt" unless defined $output;

if($exprn_threshold <= 0){
  die "Please use a positive integer value for --exprn_threshold\n";
}
my $high_exprn_threshold = $exprn_threshold;
my $low_exprn_threshold = 0 - $exprn_threshold;




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


# =============================== #
# Process the expression z-scores
# =============================== #

my %master_cell_lines_seen;		# store all cell lines seen
my %master_genes_seen;			# store all genes seen

my %wtsi_exprn_z;		# actual z-scores
my %wtsi_exprn_class;	# store: -1/0/1 for under/normal/overexpression

if(defined $wtsi_expression_z_data){
	# Open and read the CCLE CNA file
	open EXPRNZ, "< $wtsi_expression_z_data" or die "Can't read expression z-scores file $wtsi_expression_z_data: $!\n";
	
	
	# gene	expression.z	cell.line
	# CHEK2_11200_ENSG00000183765	NA	SW13_ADRENAL_GLAND
	# CHEK2_11200_ENSG00000183765	-2.00063333901991	NH12_AUTONOMIC_GANGLIA
	# CHEK2_11200_ENSG00000183765	-1.10868384234017	CHP212_AUTONOMIC_GANGLIA
	
	my $header = <EXPRNZ>;
	
	while(<EXPRNZ>){
	
	  my @fields = split(/\t/);
	  my $standard_gene = $fields[0];
	  my $exprnz = $fields[1];
	  my $standard_cell_line = $fields[2];
	  
	  chomp $standard_cell_line;
	  $standard_gene =~ s/^([^_]+)_([^_]+)_/$1\t$2\t/; # replace the underscores following the symbol and EntrezGeneID for tabs
	  next if $exprnz =~ /NA/;
	  
	  # lookup gene and cell line name in dictionaries
	
	  # skip processing this variant if the gene symbol is not in the set of output genes
	  if(!exists $output_genes{$standard_gene}){
        print "Skipping expression for gene: $standard_gene\n";
		next;
	  }
	
	  my $wtsi_key = "$standard_cell_line\t$standard_gene";
	  
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\texprn";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\texprn";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	  
	  # store the actual z-score
	  $wtsi_exprn_z{$wtsi_key} = $exprnz;
	  
      if($exprnz >= $high_exprn_threshold){	# overexpressed
		$wtsi_exprn_class{$wtsi_key} = 1;
	  }
	  elsif($exprnz <= $low_exprn_threshold){ # underexpressed
		$wtsi_exprn_class{$wtsi_key} = -1;
	  }
	  else{
		$wtsi_exprn_class{$wtsi_key} = 0;
	  }
	  
	  # add expression z-scores to the detail table hash
	  if(exists $details_table{$wtsi_key}{'exprnz'}){
		$details_table{$wtsi_key}{'exprnz'} .= "; $exprnz";
	  }
	  else{
		$details_table{$wtsi_key}{'exprnz'} = $exprnz;
	  }
	  
	  
	} # finished reading file, close
	close EXPRNZ;
}


# ==================================== #
# Process the hashes to format outputs
# ==================================== #

# get lists of all genes and cell lines seen in any data set
my @cell_lines_seen = keys %master_cell_lines_seen;

#my @genes_seen = keys %master_genes_seen;
my @genes_seen = keys %master_genes_seen;		# the list of genes from CGC and C5000


# We need to indicate which data sets were available (mut, cna or expression) as a column
# in the output.

my $exprn_zscore_matrix = '';		# expression z-scores
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
  
  $exprn_zscore_matrix .= "$seen_cell_line";
  $mutation_matrix .= "$seen_cell_line";
  $other_matrix .= "$seen_cell_line";
  $mutation_classification .= "$seen_cell_line";
  
  foreach my $seen_gene (@genes_seen){
    
    my $hash_key = "$seen_cell_line\t$seen_gene";
    my $hom_del = 0;
    my $loss = 0;
    my $amp = 0;
    my $trunc = 0;
    my $rec_mis = 0;
    my $other = 0;
    my $exprn = 0;
    
    # write out the actual expression z-score
    
    $exprn_zscore_matrix .= "\t$wtsi_exprn_z{$hash_key}";
    
    # now decide how to classify the cell line / gene...
    my $matrix_value = 0;
    
    if($output_genes{$seen_gene} eq 'OG'){
      if($wtsi_exprn_class{$hash_key} == 1){
        $matrix_value = 1;
      }
      elsif($wtsi_exprn_class{$hash_key} == -1){
        $other = 1;
      }
    }
    elsif($output_genes{$seen_gene} eq 'TSG'){
      if($wtsi_exprn_class{$hash_key} == -1){
        $matrix_value = 1;
      }
      elsif($wtsi_exprn_class{$hash_key} == 1){
        $other = 1;
      }
    }

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
    elsif($output_genes{$seen_gene} eq 'OG' && $wtsi_exprn_class{$hash_key} == 1){
      $mutation_classification .= "\t6";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $wtsi_exprn_class{$hash_key} == -1){
      $mutation_classification .= "\t7";
    }
    else{
      $mutation_classification .= "\t0";
    }
    
  }
  
  $exprn_zscore_matrix .= "\n";
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
#

open DETAILS, "> $output.details.txt" or die "Can't write details to $output:$! \n";

print DETAILS "cell_line\tgene_symbol\tEntrez_geneID\tEnsembl_GeneID\tdata set\tconsequence\tmutation\tfrequency (if SNV)\n";

my @details_keys = keys %wtsi_exprn_class;
foreach my $details_key (@details_keys){

  if($wtsi_exprn_class{$details_key} == 1){
    print DETAILS "$details_key";
    print DETAILS "\tWTSI (exprn)";
    chomp($wtsi_exprn_z{$details_key});
    print DETAILS "\thigh expression";
    print DETAILS "\t$wtsi_exprn_z{$details_key}";
    print DETAILS "\tNA\n";
  }
  elsif($wtsi_exprn_class{$details_key} == -1){
    print DETAILS "$details_key";
    print DETAILS "\tWTSI (exprn)";
    chomp($wtsi_exprn_z{$details_key});
    print DETAILS "\tlow expression";
    print DETAILS "\t$wtsi_exprn_z{$details_key}";
    print DETAILS "\tNA\n";
  }

}

close DETAILS;

open EXPRNZMAT, "> $output.expression_zscores.txt" or die "Can't write to expression z-score output file $output: $!\n";
print EXPRNZMAT "$matrix_header$exprn_zscore_matrix";
close EXPRNZMAT;


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
#  process_WTSIexprn_data.pl
#  James Campbell (jamesc\@icr.ac.uk)
# ---------------------------------- #

Usage:
perl process_WTSIexprn_data.pl [options]

Options:
  --help                    Display this message and quit
  --wtsi_expression_z_data  Path to the WTSI expression z-scores file [required]
  --exprn_threshold         Integer value for extreme expression threshold [optional, default value of 3]
  --genes                   Path to the gene name dictionary [required]
  --cell_lines              Path to the cell line name dictionary [required]
  --output_genes_list       Path to the list of genes requested [required]
  --mut_consequence_classes Path to mutation consequences resource [required]
  --output                  Path to output. [optional]
  
END

  print $usage;
}


