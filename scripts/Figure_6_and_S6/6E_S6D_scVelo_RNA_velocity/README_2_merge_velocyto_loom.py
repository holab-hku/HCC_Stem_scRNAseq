import loompy

dir = '../../../../HCC_Stem_scRNAseq'
data_dir = dir+'/data/Figure_6_and_S6/6E_S6D_scVelo_RNA_velocity/processed_data_for_scVelo'

files = [data_dir+'/Day3_Pos.loom', data_dir+'/Day10_Pos.loom', data_dir+'/Day30_Pos.loom']
loompy.combine(files, data_dir+'/Pos.loom')

files = [data_dir+'/Day3_Neg.loom', data_dir+'/Day10_Neg.loom', data_dir+'/Day30_Neg.loom']
loompy.combine(files, data_dir+'/Neg.loom')
