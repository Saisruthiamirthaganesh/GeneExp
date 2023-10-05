import pandas as pd 

data = pd.read_csv('/home/coderpad/data/mock_gene_counts.csv', index_col = 0)
# There appears to be around 1000 genes, 3 samples under condition 1 and 3 samples under condition 2. Based on the information provided, it can be assumed that the study is interested in understanding the gene expression across two different scenarios - possibly case-control. This experimental design can enable the understanding of genes implicated by comparing the expression levels of a normal individual to potentially identify the incidence of a disease or a perturbed system.

data.describe() # Summary statistics pre-normalization
def normalizeRNA(df): # Function to Normalize for read depth
    temp = df/df.sum() # Sample mean 
    z_score = (temp-temp.mean(axis = 0))/(temp.std(axis = 0))
    return(z_score)

normalized_data = normalizeRNA(data)
normalizeRNA(data).describe() # Post-normalization summary statistics indicating mean zero and variance 1

def Expression(df,col):
    df.rename(columns = {df.columns[0]: 'Baseline'}, inplace = True)
    fc_values = df.iloc[:,1:].div(df.Baseline, axis = 0) # Calculate Fold change for condition against baseline value
    # A fold-Change cutoff of 3 and -3 was arbitrarily chosen for upregulated and downregulated genes respectively
    up_exp = fc_values[col][fc_values[col]> 3]
    down_exp = fc_values[col][fc_values[col]< -3]
    upregulated_genes = []
    downregulated_genes = []
    for row in up_exp.index:
        upregulated_genes.append(row)
    for row in down_exp.index:
        downregulated_genes.append(row)
    return(upregulated_genes,downregulated_genes)

upregulated_genes,downregulates_genes = Expression(normalized_data,'condition_1_2') # Returns expression lists for Condition as specified in the column name

class GeneExpressionAnalyzer:

    def normalizeRNA(df): # Function to Normalize for read depth
        temp = df/df.sum() # Sample mean 
        z_score = (temp-temp.mean(axis = 0))/(temp.std(axis = 0))
        return(z_score)

    def Expression(df,col):
        normalized_data = normalizeRNA(df)
        normalized_data.rename(columns = {normalized_data.columns[0]: 'Baseline'}, inplace = True)
        fc_values = normalized_data.iloc[:,1:].div(normalized_data.Baseline, axis = 0)  
        # Calculate Fold change for condition against baseline value
        # A fold-Change cutoff of 3 and -3 was arbitrarily chosen for upregulated and downregulated genes respectively
        up_exp = fc_values[col][fc_values[col]>3]
        down_exp = fc_values[col][fc_values[col]< -3]
        upregulated_genes = []
        downregulated_genes = []
        for row in up_exp.index:
            upregulated_genes.append(row)
        for row in down_exp.index:
            downregulated_genes.append(row)
        return(upregulated_genes,downregulated_genes)



upregulated_genes,downregulated_genes = GeneExpressionAnalyzer.Expression(data,'condition_1_2') # Returns expression lists for Condition as specified in the column name with normalization and fold change calculation encapsulated within the Class.
print(f'The Upregulated genes for the specified condition are',upregulated_genes)
print(f'The Downregulated genes for the specified condition are',downregulated_genes)