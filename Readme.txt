Read me DmojSSRsupplementCode

Pipeline for performing RNA-seq expression analysis: RNAseq_expression_pipeline.txt
	-This pipeline calls on the script "edgeR_expression.R"
		-"edgeR_expression.R" optionally uses "DmojID.csv" and "Or_list.csv" in the folder "sample_files" as input files
	-Target gene sequences used in this study listed in file "targets.fa" included in folder "sample_files"
	-Final annotation of FBgn gene name for each gene in "targets.fa" listed in file "FBgnAnnotations.xlsx" in folder "sample_files"

Pipeline for identifying amino acid substitutions in a given gene of interest: aa_changes_pipeline.txt
	-This pipeline calls on the scripts "translation.py" and "MSA_summary_Dmoj.py" 