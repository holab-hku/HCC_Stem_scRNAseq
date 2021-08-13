#!/usr/bin/perl
use strict; use warnings;

###################################
my $bam_dir = "./outputs";
my $subset = "HCC";
###################################

my $dir = "../../../../HCC_Stem_scRNAseq"
my $output_dir = $dir . "/data/Figure_6_and_S6/6E_S6D_scVelo_RNA_velocity/processed_data_for_scVelo";

my $gtf = "$dir/data/Misc/Mus_musculus.GRCm38.93.gtf";

die"
Usage: run_velocyto <day>
" unless @ARGV==1;

my $day = $ARGV[0];
if ($day !~ /^Day\d+$/) {
	die "<day> $day not recognized\n";
}

for my $type (qw(Neg Pos)) {
	
	my $file = "$bam_dir/$subset" . "_" . $day . "_" . $type . "_BAMs/possorted_genome_bam.$subset.$day.$type.bam";
	
	my $cellsorted_file =  "$bam_dir/$subset" . "_" . $day . "_" . $type . "_BAMs/cellsorted_possorted_genome_bam.$subset.$day.$type.bam";
	
	system("samtools sort -t CB -O BAM -@ 12 -o $cellsorted_file $file");
	
	my $sample_id = $day . "_" . $type;
	system("velocyto run -e $sample_id -o $output_dir $file $gtf");
	
}
