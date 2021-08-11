#!/usr/bin/perl
use strict; use warnings;

my $output_dir = "HCC_Stem_scRNAseq/data/Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation/Cell_Ranger_filtered_matrices";
mkdir $output_dir unless (-d $output_dir);

##########TO BE FINALIZED
##Download from GEO - 
#my $Day3_Neg = "/Volumes/My_Passport/MaS_10X5RNA_SS-190711-02a/MaS_10X5RNA_SS-190711-02a_AnalysisData/MaS_10X5RNA_SS-190711-02a_AnalysisData/cellranger_v3/tdTneg-day3/count_tdTneg-day3/outs/filtered_feature_bc_matrix";
#my $Day3_Pos = "/Volumes/My_Passport/MaS_10X5RNA_SS-190711-02a/MaS_10X5RNA_SS-190711-02a_AnalysisData/MaS_10X5RNA_SS-190711-02a_AnalysisData/cellranger_v3/tdTpos-day3/count_tdTpos-day3/outs/filtered_feature_bc_matrix";

#my $Day10_Neg = "/Volumes/My_Passport/MaS_10X5RNA_SS-190726-02a/MaS_10X5RNA_SS-190726-02a_AnalysisData/MaS_10X5RNA_SS-190726-02a_AnalysisData/cellranger_v3/tdTneg-day10/count_tdTneg-day10/outs/filtered_feature_bc_matrix";
#my $Day10_Pos = "/Volumes/My_Passport/MaS_10X5RNA_SS-190726-02a/MaS_10X5RNA_SS-190726-02a_AnalysisData/MaS_10X5RNA_SS-190726-02a_AnalysisData/cellranger_v3/tdTpos-day10/count_tdTpos-day10/outs/filtered_feature_bc_matrix";

#my $Day30_Neg = "/Volumes/My_Passport/MaS_10X5RNA_SS-190711-02b/MaS_10X5RNA_SS-190711-02b_AnalysisData/MaS_10X5RNA_SS-190711-02b_AnalysisData/cellranger_v3/Fixed-negative/count_Fixed-negative/outs/filtered_feature_bc_matrix";
#my $Day30_Pos = "/Volumes/My_Passport/MaS_10X5RNA_SS-190711-02b/MaS_10X5RNA_SS-190711-02b_AnalysisData/MaS_10X5RNA_SS-190711-02b_AnalysisData/cellranger_v3/Fixed-positive/count_Fixed-positive/outs/filtered_feature_bc_matrix";

my @dirs = ($Day3_Neg, $Day3_Pos, $Day10_Neg, $Day10_Pos, $Day30_Neg, $Day30_Pos);

for my $dir (@dirs) {
	my ($platform) = $dir =~ m/.*\/(cellranger_v3)\//;
	my ($day) = $dir =~ m/.*day(\d+)\//;
	if (!$day) {
		$day = 30;
	}
	
	if ( (!$day) || (!$platform)) {
		die "Can't find patient or platform\n";
	}
	
	my $type;
	#my $type_long;
	if ($dir =~ m/\/(tdTneg|Fixed-negative).*\//) {
		$type = "Neg";
	} elsif ($dir =~ m/\/(tdTpos|Fixed-positive).*\//) {
		$type = "Pos";
	}
	
	my $out_path = $output_dir . "/Day$day";
	mkdir $out_path unless (-d $out_path);
	#my $metadata = $out_path . "/metadata";
	#mkdir $metadata unless (-d $metadata);
	
	$out_path = $out_path . "/$platform";
	mkdir $out_path unless (-d $out_path);
	$out_path = $out_path . "/$type";
	mkdir $out_path unless (-d $out_path);
	
	my $output_file = $out_path . "/barcodes.tsv";
	open my $OUT, ">$output_file" or die "can't open $output_file\n";
	my $metadata_out = $out_path . "/metadata";
	open my $META, ">$metadata_out" or die "can't open $metadata_out\n";
	print $META "\t" . "type\t" . "day\t" . "day_type\n";
	
	open my $IN, "gunzip -c $dir/barcodes.tsv.gz |" or die "can't open barcodes.tsv.gz\n";
	while (<$IN>) {
		chomp;
		if ($_ =~ /([ACTG]+\-)1$/) {
			#print $1 . "\n";
			#die;
			print $OUT $1 . "$day$type\n";
			
			print $META $1 . "$day$type\t$type\t" . "Day$day\t" . "Day$day" . "_" . "$type\n";
		} else {
			die "Can't match barcode: " . $_ . "\n" . "in $out_path\n";
		}
	}
	close $IN;
	close $OUT;
	close $META;
	
	system("gzip $output_file");
	
	system("cp $dir/features.tsv.gz $out_path/features.tsv.gz");
	system("cp $dir/matrix.mtx.gz $out_path/matrix.mtx.gz");
	#die;
	
}