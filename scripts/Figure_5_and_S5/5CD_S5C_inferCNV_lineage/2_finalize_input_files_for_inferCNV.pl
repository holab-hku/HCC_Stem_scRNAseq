#!/usr/bin/perl
use strict; use warnings;

my $dir = "../../../../HCC_Stem_scRNAseq";
my $figure = "Figure_5_and_S5/5CD_S5C_inferCNV_lineage";

my $CNV_dir = "$dir/outs/$figure/lineage_counts_matrix_and_meta";

my $gtf_file = "$dir/data/Misc/Mus_musculus.GRCm38.93.gtf.gz";

open my $IN, "<$CNV_dir/exprMatrix.tsv" or die "can't open $CNV_dir/exprMatrix.tsv";
my $header = <$IN>;
my %genes_seen;
while (<$IN>) {
	my @line = split("\t", $_);
	$genes_seen{$line[0]}++;
}
close $IN;

my $meta_index = 3;
my $input_file = "$CNV_dir/meta.tsv";
my $output_file = "$CNV_dir/meta.noheader.tsv";
open $IN, "<$input_file" or die "can't open $input_file\n";
open my $OUT, ">$output_file" or die "can't open $output_file\n";
$header = <$IN>;
while (<$IN>) {
	chomp;
	my @line = split("\t", $_);
	print $OUT $line[0] . "\t" . $line[$meta_index] . "\n";
}
close $OUT;
close $IN;

my %gene_positions;

my %genes_out;

open $IN, "gunzip -c $gtf_file | " or die "can't open $gtf_file\n";
while (<$IN>) {
	if ($_ =~ /^#/) {
		next;
	}
	my @line = split("\t", $_);
	if ($line[2] eq "gene") {
		my ($gene_name) = $_ =~ /.*gene_name \"(.*?)\"/;
		if (exists $genes_seen{$gene_name}) {
			$gene_positions{$gene_name} = "$line[0]\t$line[3]\t$line[4]";
			$genes_out{$gene_name}++;
		}
	}
}
close $IN;

$output_file = "$CNV_dir/gene_ordering.tsv";
open $OUT, ">$output_file" or die "can't open $output_file\n";
for my $gene (keys %gene_positions) {
	print $OUT $gene . "\t" . $gene_positions{$gene} . "\n";
}
close $OUT;


open $OUT, ">$CNV_dir/exprMatrix.genes_out.tsv";
open $IN, "<$CNV_dir/exprMatrix.tsv";
$header = <$IN>;
print $OUT $header;

while (<$IN>) {
	my @line = split("\t", $_);
	if (exists $genes_out{$line[0]}) {
		print $OUT $_;
	}
}
close $OUT;
close $IN;

#system("gzip $CNV_dir/exprMatrix.genes_out.tsv");
