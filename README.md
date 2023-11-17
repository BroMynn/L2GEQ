# L2GEQ
L2GEQ tries to quantify a gene expression difference/similarity between two groups in a single value. It was developed to aid a project in Nation Cancer Center (Republic of Korea) of Dr. Youngjoo Lee, M.D., Ph.D., and was a poster presentation at 2023 AACR-KCA Joint Conference on Precision Medicine in Cancer on day 2 (2023-11-17 Fri).
Below is an example of calculating L2GEQ on a RNA-seq experiment data normalized in transcript-per-million (TPM), the first column being gene names the first cell being "Gene_id", and each row representing a gene expression value for each gene of each sample except for the first cell being the name of your sample. Note L2GEQ does not automatically parse the sample names to assign each sample to the appropriate groups. Users need to designate group labels for each sample using .designate_labels() method, accepting a list containing labels for each sample in the order presented in your gene expression file.


## Basic Use
resistant_cancer_cell_lines = L2GEQ(r"YourDirectory\lung_cancer_ensembl_tpm.tsv", sep="\t")
resistant_cancer_cell_lines.designate_label([
    'parent_cancer', 'parent_cancer', 'parent_cancer',
    'resistant_cancer_drug1', 'resistant_cancer_drug1', 'resistant_cancer_drug1',
    'resistant_cancer_drug2', 'resistant_cancer_drug2', 'resistant_cancer_drug2',
    'resistant_cancer_drug3', 'resistant_cancer_drug3', 'resistant_cancer_drug3',
  ]
)
resistant_cancer_cell_lines.calculate_L2GEQ(sig_test='ttest', avg_method='harmonic', weight_method='limit1_circular', weight_denom='multiply')

These three lines of code will give you L2GEQ value between every 2-member combination of experimental groups (parent_cancer vs. resistant_cancer_drug1, parent_cancer vs. resistant_cancer_drug2, ..., resistant_cancer_drug2 vs. resistant_cancer_drug3) with the designated L2GEQ parameters.


## Additional Functionalities
1. If you want to iterate every parameters implemented in the current version of L2GEQ, simply use:
resistant_cancer_cell_lines.parameter_iterator(mode='main', sig_test='ttest', avg_method='harmonic', weight_method='limit1_circular', weight_denom='multiply')

Note that if you wish to perform that ordinary L2GEQ calculation pass mode='main' as an argument, or pass nothing as it is the default value for the function. If you wish to calculate L2GEQ between samples belong to the same label, a mode of L2GEQ calculation which I explain right below, pass mode='within' as an argument.


2. If you wish to calculate L2GEQ values between the samples that belong to the same group (it could be useful if you are checking integrity of your samples beloging to the same group, or trying to determine your threshold of an identical/similar gene expression), try below:
resistant_cancer_cell_lines.calculate_withinL2GEQ(sig_test='ttest', avg_method='harmonic', weight_method='limit1_circular', weight_denom='multiply')

Note that you need at least one label containing at least four samples to perform .calculate_withinL2GEQ() method, and I do not recommend performing this unless you have 6 or more samples in any of the labels, as calculating L2GEQ in any combinations of parameter is not recommended in a group containing less than at least three samples. Note that in performing calculate_withinL2GEQ(), L2GEQ tries to split samples belonging to the same label in two groups allowing no common members or duplicate sample in the group. That's why I recommend to use this method if you have at least 6 samples of the same group.


3. You can perform L2GEQ calculation in the samples of your choice. Use below:
resistant_cancer_cell_lines.L2GEQ_custom_index(label1=[0,1,2], label2=[3,4,5,6,7,8], sig_test='ttest', avg_method='harmonic', weight_method='limit1_circular', weight_denom='multiply')

Use caution as the implementation of assigning comparison groups is somewhat different from other calculate_L2GEQ methods. In .L2GEQ_custom_index() method, users are required to pass the index numbers of labels in a list. It will try to calculate L2GEQ by comparing label1 vs. label2 using the indices you pass.


4. L2GEQ has an additional functionality of extracting differentially expressed genes (DEGs) based on the L2GEQ value previously calculated. Refer to below:
resistant_cancer_cell_lines.calculate_L2GEQ(cutoff=0.0005)
resistant_cancer_cell_lines.set_annotation_db(r"YourDirectory\Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip", sep='\t')
resistant_cancer_cell_lines.extract_DEGs()

You need a previous calculation of L2GEQ to activate .extract_DEGs() method. You also need a conversion table from Ensembl Id to HGNC Gene Symbol, such as one provided by MSigBD, if you are using a gene expression table with gene names in Ensembl ID. It has a cutoff parameter, which filters out the portion of L2GEQ value below that number (If you have a L2GEQ value of 100 and you threshold is 0.01, only genes that contribute to the upper 99 value will be extracted as DEGs)



If you want to contact me, my email is bromynn@ncc.re.kr.
