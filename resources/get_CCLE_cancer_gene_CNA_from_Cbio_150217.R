# download CCLE gistic data for kinase genes from cBioPortal
# jamesc@icr.ac.uk, 12th Jan 2015

setwd("/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/resources/")

require(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)


#
# New cancer gene set from CGC + Santarius + Beroukhim
# using HGNC symbols from BioMart, based on EntrezGeneIDs
#

cancer_genes <- c(
	"SEPT6",
	"SEPT9",
	"ABI1",
	"ABL1",
	"ABL2",
	"ACKR3",
	"ACSL3",
	"ACSL6",
	"AFF1",
	"AFF3",
	"AFF4",
	"AKAP9",
	"AKT1",
	"AKT2",
	"AKT3",
	"ALDH2",
	"ALK",
	"AMER1",
	"ANKS1B",
	"APC",
	"AR",
	"ARHGAP26",
	"ARHGEF12",
	"ARID1A",
	"ARID2",
	"ARNT",
	"ARPC1A",
	"ASPSCR1",
	"ASXL1",
	"ATF1",
	"ATIC",
	"ATM",
	"ATP1A1",
	"ATP2B3",
	"ATRX",
	"AURKA",
	"AXIN1",
	"BAG4",
	"BAP1",
	"BCL10",
	"BCL11A",
	"BCL11B",
	"BCL2",
	"BCL2L1",
	"BCL2L2",
	"BCL3",
	"BCL6",
	"BCL7A",
	"BCL9",
	"BCOR",
	"BCR",
	"BIRC2",
	"BIRC3",
	"BLM",
	"BMPR1A",
	"BRAF",
	"BRCA1",
	"BRCA2",
	"BRD3",
	"BRD4",
	"BRIP1",
	"BTG1",
	"BUB1B",
	"C11orf30",
	"C15orf65",
	"C2orf44",
	"C8orf4",
	"CACNA1D",
	"CACNA1E",
	"CALR",
	"CAMTA1",
	"CANT1",
	"CARD11",
	"CARS",
	"CASC5",
	"CASP8",
	"CBFA2T3",
	"CBFB",
	"CBL",
	"CBLB",
	"CBLC",
	"CCDC6",
	"CCNB1IP1",
	"CCND1",
	"CCND2",
	"CCND3",
	"CCNE1",
	"CD274",
	"CD74",
	"CD79A",
	"CD79B",
	"CDC6",
	"CDC73",
	"CDH1",
	"CDH11",
	"CDH13",
	"CDK12",
	"CDK4",
	"CDK6",
	"CDKN2A",
	"CDKN2B",
	"CDKN2C",
	"CDX2",
	"CEBPA",
	"CEP89",
	"CHCHD7",
	"CHD1L",
	"CHEK2",
	"CHIC2",
	"CHN1",
	"CIC",
	"CIITA",
	"CKS1B",
	"CLIP1",
	"CLP1",
	"CLTC",
	"CLTCL1",
	"CNBP",
	"CNOT3",
	"CNTNAP2",
	"CNTRL",
	"COL1A1",
	"COL2A1",
	"COPS3",
	"COX6C",
	"CREB1",
	"CREB3L1",
	"CREB3L2",
	"CREBBP",
	"CRKL",
	"CRLF2",
	"CRTC1",
	"CRTC3",
	"CSF3R",
	"CSMD1",
	"CTNNB1",
	"CUX1",
	"CYLD",
	"DAXX",
	"DCC",
	"DCTN1",
	"DCUN1D1",
	"DDB2",
	"DDIT3",
	"DDX10",
	"DDX5",
	"DDX6",
	"DEK",
	"DICER1",
	"DLG2",
	"DMD",
	"DNM2",
	"DNMT3A",
	"DSCAM",
	"DYRK2",
	"E2F3",
	"EBF1",
	"ECT2L",
	"EEF1A2",
	"EGFR",
	"EIF3E",
	"EIF4A2",
	"EIF5A2",
	"ELF4",
	"ELK4",
	"ELL",
	"ELN",
	"EML4",
	"EP300",
	"EPS15",
	"ERBB2",
	"ERC1",
	"ERCC2",
	"ERCC3",
	"ERCC4",
	"ERCC5",
	"ERG",
	"ETV1",
	"ETV4",
	"ETV5",
	"ETV6",
	"EWSR1",
	"EXT1",
	"EXT2",
	"EZH2",
	"EZR",
	"FADD",
	"FAM46C",
	"FANCA",
	"FANCC",
	"FANCD2",
	"FANCE",
	"FANCF",
	"FANCG",
	"FAS",
	"FBXO11",
	"FBXW7",
	"FCGR2B",
	"FCRL4",
	"FEV",
	"FGFR1",
	"FGFR1OP",
	"FGFR2",
	"FGFR3",
	"FH",
	"FHIT",
	"FIP1L1",
	"FLCN",
	"FLI1",
	"FLT3",
	"FNBP1",
	"FOXA1",
	"FOXL2",
	"FOXO1",
	"FOXO3",
	"FOXO4",
	"FOXP1",
	"FSTL3",
	"FUBP1",
	"FUS",
	"GAS7",
	"GATA1",
	"GATA2",
	"GATA3",
	"GATA6",
	"GMPS",
	"GNA11",
	"GNAQ",
	"GNAS",
	"GOLGA5",
	"GOPC",
	"GPC3",
	"GPC5",
	"GPC6",
	"GPHN",
	"GRB7",
	"H3F3A",
	"H3F3B",
	"HERPUD1",
	"HEY1",
	"HIP1",
	"HIST1H4I",
	"HLA-A",
	"HLF",
	"HMGA1",
	"HMGA2",
	"HMGN2P46",
	"HNF1A",
	"HNRNPA2B1",
	"HOOK3",
	"HOXA11",
	"HOXA13",
	"HOXA9",
	"HOXC11",
	"HOXC13",
	"HOXD11",
	"HOXD13",
	"HRAS",
	"HSP90AA1",
	"HSP90AB1",
	"IDH1",
	"IDH2",
	"IGF1R",
	"IKZF1",
	"IL2",
	"IL21R",
	"IL6ST",
	"IL7R",
	"IRF4",
	"IST1",
	"ITK",
	"JAK1",
	"JAK2",
	"JAK3",
	"JAZF1",
	"JUN",
	"KAT6A",
	"KAT6B",
	"KCNJ5",
	"KDM5A",
	"KDM5C",
	"KDM6A",
	"KDR",
	"KDSR",
	"KIAA1549",
	"KIAA1598",
	"KIF5B",
	"KIT",
	"KLF4",
	"KLF6",
	"KLK2",
	"KMT2A",
	"KMT2C",
	"KMT2D",
	"KRAS",
	"KTN1",
	"LASP1",
	"LCK",
	"LCP1",
	"LHFP",
	"LIFR",
	"LMNA",
	"LMO1",
	"LMO2",
	"LPP",
	"LRIG3",
	"LRP1B",
	"LSM1",
	"LSM14A",
	"LYL1",
	"MACROD2",
	"MAF",
	"MAFB",
	"MAGI2",
	"MALAT1",
	"MALT1",
	"MAML2",
	"MAP2K1",
	"MAP2K2",
	"MAP2K4",
	"MAP3K5",
	"MAX",
	"MCL1",
	"MDM2",
	"MDM4",
	"MECOM",
	"MED12",
	"MED29",
	"MEN1",
	"MET",
	"MITF",
	"MKL1",
	"MLF1",
	"MLH1",
	"MLLT1",
	"MLLT10",
	"MLLT11",
	"MLLT3",
	"MLLT4",
	"MLLT6",
	"MN1",
	"MNX1",
	"MPL",
	"MSH2",
	"MSH6",
	"MSI2",
	"MSN",
	"MTCP1",
	"MTDH",
	"MUC1",
	"MUTYH",
	"MYB",
	"MYC",
	"MYCL",
	"MYCN",
	"MYD88",
	"MYH11",
	"MYH9",
	"MYO5A",
	"NAALADL2",
	"NAB2",
	"NACA",
	"NBN",
	"NCKIPSD",
	"NCOA1",
	"NCOA2",
	"NCOA3",
	"NCOA4",
	"NDRG1",
	"NEGR1",
	"NF1",
	"NF2",
	"NFATC2",
	"NFE2L2",
	"NFIB",
	"NFKB2",
	"NIN",
	"NKX2-1",
	"NKX2-8",
	"NONO",
	"NOTCH1",
	"NOTCH2",
	"NPM1",
	"NR4A3",
	"NRAS",
	"NRG1",
	"NSD1",
	"NT5C2",
	"NTM",
	"NTRK1",
	"NTRK3",
	"NUMA1",
	"NUP214",
	"NUP98",
	"NUTM1",
	"NUTM2A",
	"NUTM2B",
	"OLIG2",
	"OMD",
	"OPCML",
	"P2RY8",
	"PAFAH1B2",
	"PAK1",
	"PALB2",
	"PARD3B",
	"PARK2",
	"PATZ1",
	"PAX3",
	"PAX5",
	"PAX7",
	"PAX8",
	"PAX9",
	"PBRM1",
	"PBX1",
	"PCM1",
	"PCSK7",
	"PDCD1LG2",
	"PDE4D",
	"PDE4DIP",
	"PDGFB",
	"PDGFRA",
	"PDGFRB",
	"PER1",
	"PHF6",
	"PHOX2B",
	"PICALM",
	"PIK3CA",
	"PIK3R1",
	"PIM1",
	"PLA2G10",
	"PLAG1",
	"PLCG1",
	"PML",
	"PMS1",
	"PMS2",
	"POT1",
	"POU2AF1",
	"POU5F1",
	"PPARG",
	"PPFIBP1",
	"PPM1D",
	"PPP2R1A",
	"PRCC",
	"PRDM1",
	"PRDM16",
	"PRF1",
	"PRKAR1A",
	"PRKCI",
	"PRKG1",
	"PRRX1",
	"PSIP1",
	"PTCH1",
	"PTEN",
	"PTK6",
	"PTPN11",
	"PTPRB",
	"PTPRC",
	"PTPRD",
	"PTPRK",
	"PTPRN2",
	"PWWP2A",
	"RAB23",
	"RAB25",
	"RABEP1",
	"RAC1",
	"RAD21",
	"RAD51B",
	"RAF1",
	"RALGDS",
	"RANBP17",
	"RAP1GDS1",
	"RARA",
	"RB1",
	"RBFOX1",
	"RBM15",
	"RECQL4",
	"REL",
	"RET",
	"RHOH",
	"RMI2",
	"RNF213",
	"RNF217-AS1",
	"RNF43",
	"ROS1",
	"RPL10",
	"RPL22",
	"RPL5",
	"RPN1",
	"RPS6KB1",
	"RSPO2",
	"RSPO3",
	"RUNX1",
	"RUNX1T1",
	"RUVBL1",
	"RYR2",
	"SBDS",
	"SDC4",
	"SDHAF2",
	"SDHB",
	"SDHC",
	"SDHD",
	"SDK1",
	"SEPT5",
	"SET",
	"SETBP1",
	"SETD2",
	"SF3B1",
	"SFPQ",
	"SH2B3",
	"SH3GL1",
	"SHC1",
	"SHH",
	"SKP2",
	"SLC34A2",
	"SLC45A3",
	"SMAD4",
	"SMARCA4",
	"SMARCB1",
	"SMARCE1",
	"SMO",
	"SMURF1",
	"SNTG1",
	"SOCS1",
	"SOX2",
	"SPECC1",
	"SRGAP3",
	"SRSF2",
	"SRSF3",
	"SS18",
	"SS18L1",
	"SSX1",
	"SSX2",
	"SSX4",
	"STAG2",
	"STARD3",
	"STAT3",
	"STAT5B",
	"STAT6",
	"STIL",
	"STK11",
	"SUFU",
	"SUZ12",
	"SYK",
	"TAF15",
	"TAL1",
	"TAL2",
	"TBL1XR1",
	"TCEA1",
	"TCF12",
	"TCF3",
	"TCF7L2",
	"TCL1A",
	"TCL6",
	"TERT",
	"TET1",
	"TET2",
	"TFE3",
	"TFEB",
	"TFG",
	"TFPT",
	"TFRC",
	"THRAP3",
	"TLX1",
	"TLX3",
	"TMPRSS2",
	"TNFAIP3",
	"TNFRSF14",
	"TNFRSF17",
	"TOP1",
	"TP53",
	"TPM3",
	"TPM4",
	"TPR",
	"TRAF7",
	"TRIB1",
	"TRIM24",
	"TRIM27",
	"TRIM33",
	"TRIP11",
	"TRRAP",
	"TSC1",
	"TSC2",
	"TSHR",
	"TTL",
	"U2AF1",
	"UBR5",
	"USP6",
	"VHL",
	"VTI1A",
	"WAS",
	"WHSC1",
	"WHSC1L1",
	"WIF1",
	"WRN",
	"WT1",
	"WWOX",
	"WWTR1",
	"XPA",
	"XPC",
	"XPO1",
	"YAP1",
	"YEATS4",
	"YWHAB",
	"YWHAE",
	"YWHAQ",
	"YWHAZ",
	"ZBTB16",
	"ZCCHC8",
	"ZMYM2",
	"ZNF217",
	"ZNF331",
	"ZNF384",
	"ZNF521",
	"ZNF639",
	"ZRSR2",
	"BCL5",
	"IGH",
	"IGL",
	"TRA",
	"TRB",
	"TRD",
	"DUX4",
	"IGK",
	"SNX"
	)



# List the cancer studies at server
getCancerStudies(mycgds)[,1]

# We want the "cellline_ccle_broad" study
ccle_study = "cellline_ccle_broad"

# List the case lists available for data sets
getCaseLists(mycgds,ccle_study)[,1]

# We want the CNA data
ccle_caselist = "cellline_ccle_broad_cna"


# Get available genetic profiles
ccle_geneticprofile = getGeneticProfiles(mycgds,ccle_study)[1,1]

ccle_gistic <- getProfileData(mycgds,cancer_genes,ccle_geneticprofile,ccle_caselist)




## prep ccle_gistic_amps (functionally relevant)
#ccle_gistic_amps <- as.matrix(ccle_gistic)
#ccle_gistic_amps[,] <- 0 
#ccle_gistic_amps[which(ccle_gistic == 2)] <- 1 

## prep ccle_gistic_amps_and_gains (for exclusion)
#ccle_gistic_amps_and_gains <- as.matrix(ccle_gistic)
#ccle_gistic_amps_and_gains[,] <- 0 
#ccle_gistic_amps_and_gains[which(ccle_gistic > 0)] <- 1 


ccle_gistic_t <- t(ccle_gistic)

write.table(
	ccle_gistic_t,
	file="CCLE_cancer_gene_gistic_calls_150217_CGC_Santarius_Beroukhim.txt",
	col.names=TRUE,
	row.names=TRUE,
	sep="\t",
	quote=FALSE
	)
	


#write.table(
#	ccle_gistic_amps,
#	file="CCLE_kinase_gistic_amps_150112.txt",
#	col.names=TRUE,
#	row.names=TRUE,
#	sep="\t",
#	quote=FALSE
#	)
	
	
#write.table(
#	ccle_gistic_amps_and_gains,
#	file="CCLE_kinase_gistic_amps_and_gains_150112.txt",
#	col.names=TRUE,
#	row.names=TRUE,
#	sep="\t",
#	quote=FALSE
#	)


