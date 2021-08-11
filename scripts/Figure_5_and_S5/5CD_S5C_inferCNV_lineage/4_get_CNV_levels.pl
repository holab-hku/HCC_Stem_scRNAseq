#!/usr/bin/perl
use strict; use warnings;

my $dir = "../../../../HCC_Stem_scRNAseq";
my $figure = "Figure_5_and_S5/5CD_S5C_inferCNV_lineage";

my $prefix = "S5C_inferCNV";

$dir = "$dir/outs/$figure/$prefix";

#modified expression cutoff
my $CUTOFF = 0.05;

#50%
my $PERC_CELLS_IN_GENE_HIGH_CNV_THRESH = 0.5;

my $obs = "$dir/infercnv.observations.txt";
my $refs = "$dir/infercnv.references.txt";

my %num_high_CNV_genes;

my %gene_num_high_CNV_cells;
my $tot_barcodes;
my $tot_genes;

incr_data($obs);
incr_data($refs);

my $output_file = "$dir/$prefix.num_high_CNV_genes.cutoff$CUTOFF.txt";
open my $OUT, ">$output_file" or die "can't open $output_file\n";
print $OUT "\t" . "Num_high_CNV_genes\n";
for my $barcode (keys %num_high_CNV_genes) {
	print $OUT $barcode . "\t" . $num_high_CNV_genes{$barcode} . "\n";
}
close $OUT;

my $num_high_CNV_genes;
my $high_CNV_gene_obs = extract_high_CNV_gene($obs);
my $high_CNV_gene_refs = extract_high_CNV_gene($refs);

print "num_high_CNV_genes: $num_high_CNV_genes\n";

#die;

my %cell_tot_CNV;

incr_CNV($high_CNV_gene_obs);
incr_CNV($high_CNV_gene_refs);

$output_file = "$dir/$prefix.hotspot_region_mean_CNV.cutoff$CUTOFF.txt";
open $OUT, ">$output_file" or die "can't open $output_file\n";
print $OUT "\t" . "Mean_Absolute_Gain_or_Loss_of_Copy_Number\n";
for my $barcode (keys %cell_tot_CNV) {
	my $mean_CNV_level = $cell_tot_CNV{$barcode}/$num_high_CNV_genes;
	print $OUT $barcode . "\t" . $mean_CNV_level . "\n";
}
close $OUT;


%cell_tot_CNV = ();
incr_CNV($obs);
incr_CNV($refs);

$output_file = "$dir/$prefix.mean_CNV.cutoff$CUTOFF.txt";
open $OUT, ">$output_file" or die "can't open $output_file\n";
print $OUT "\t" . "Mean_Absolute_Gain_or_Loss_of_Copy_Number\n";
for my $barcode (keys %cell_tot_CNV) {
	my $mean_CNV_level = $cell_tot_CNV{$barcode}/$tot_genes;
	print $OUT $barcode . "\t" . $mean_CNV_level . "\n";
}
close $OUT;

sub incr_CNV {
	my $file = shift;
	my %num_cols;
	
	open my $IN, "<$file" or die "can't open $file\n";
	my $header = <$IN>;
	$header =~ s/"//g;
	my @header_ar = split(/\s+/, $header);
	chomp(@header_ar);
	print "Num_barcodes: " . scalar(@header_ar) . "\n";

	while (<$IN>) {
		chomp;
		my @line = split(/\s+/, $_);
		
		my $gene = shift @line;
		
		$num_cols{scalar(@line)}++;
		
		for (my $i=0; $i <@line; $i++) {
			my $modified_expr = $line[$i];
			$modified_expr -= 1;
			$modified_expr = abs($modified_expr);
			$cell_tot_CNV{$header_ar[$i]} += $modified_expr;
		}
	}
	close $IN;
	
	for my $num_col (keys %num_cols) {
		print "$num_col" . ":" . " " . $num_cols{$num_col} . " genes\n";
	}
}

sub extract_high_CNV_gene {
	my $file = shift;
	my $output_file = $file . ".high_CNV_gene";
	open my $OUT, ">$output_file" or die "can't open $output_file\n";
	open my $IN, "<$file" or die "can't open $file\n";
	my $header = <$IN>;
	print $OUT $header;
	my $num_genes;
	while (<$IN>) {
		my @line = split(/\s+/, $_);
		
		my $gene = shift @line;
		
		if ( ($gene_num_high_CNV_cells{$gene} / $tot_barcodes) >= $PERC_CELLS_IN_GENE_HIGH_CNV_THRESH) {
			print $OUT $_;
			$num_genes++;
		}
	}
	close $IN;
	close $OUT;
	if (!$num_high_CNV_genes) {
		$num_high_CNV_genes = $num_genes;
	}
	return $output_file;
}

sub incr_data {
	my $file = shift;
	
	my %num_cols;
	
	open my $IN, "<$file" or die "can't open $file\n";
	my $header = <$IN>;
	$header =~ s/"//g;
	my @header_ar = split(/\s+/, $header);
	chomp(@header_ar);
	print "Num_barcodes: " . scalar(@header_ar) . "\n";
	
	$tot_barcodes += scalar(@header_ar);
	
	for my $barcode (@header_ar) {
		$num_high_CNV_genes{$barcode} = 0;
	}
	
	while (<$IN>) {
		chomp;
		my @line = split(/\s+/, $_);
		
		my $gene = shift @line;
		
		
		
		$num_cols{scalar(@line)}++;
		
		for (my $i=0; $i <@line; $i++) {
			my $modified_expr = $line[$i];
			$modified_expr -= 1;
			$modified_expr = abs($modified_expr);
			if ($modified_expr >= $CUTOFF) {
				$num_high_CNV_genes{$header_ar[$i]}++;
				$gene_num_high_CNV_cells{$gene}++;
			}
		}
	}
	close $IN;
	
	for my $num_col (keys %num_cols) {
		print "$num_col" . ":" . " " . $num_cols{$num_col} . " genes\n";
		$tot_genes = $num_cols{$num_col};
	}
	
}