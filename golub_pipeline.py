# data from https://www.openintro.org/data/index.php?data=golub

# Visualize high-dimensional gene expression data using UMAP, with points colored by cancer type.

import pandas as pd
import umap
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('golub.csv')
gene_expression = data.select_dtypes(include=[float])
cancer_type = data.iloc[:, -1]

# Apply UMAP
reducer = umap.UMAP(n_components=2, random_state=42)
embedding = reducer.fit_transform(gene_expression)

# Create a DataFrame with the UMAP results
umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
umap_df['CancerType'] = cancer_type

# Plot the UMAP results
plt.figure(figsize=(10, 8))
sns.scatterplot(x='UMAP1', y='UMAP2', hue='CancerType', data=umap_df, palette='viridis', s=50)
plt.title('UMAP Projection of Gene Expression Data')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.legend(title='Cancer Type', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()


# Perform PCA and K-means clustering on scaled gene expression data, then visualize clusters on a UMAP projection.

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Scale the gene expression data
scaler = StandardScaler()
scaled_gene_expression = scaler.fit_transform(gene_expression)

# Calculate the first 5 principal components
pca = PCA(n_components=5)
principal_components = pca.fit_transform(scaled_gene_expression)

# Apply k-means clustering to the principal components
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(principal_components)

# Add the cluster assignments to the UMAP DataFrame
umap_df['Cluster'] = clusters

# Visualize
plt.figure(figsize=(10, 8))
sns.scatterplot(x='UMAP1', y='UMAP2', hue='Cluster', data=umap_df, palette='viridis')
plt.title('UMAP projection of gene expression data colored by K-means clusters')
plt.show()


# Identify genes differentially expressed between AML and other samples using linear regression, adjusting for sample source and correcting for multiple testing.


import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

# Add a binary variable for AML
data['AML'] = (data['cancer'] == 'aml').astype(int)

# Determine sample source (BM or PB)
data['SampleSource'] = data['tissue.mf'].apply(lambda x: 1 if 'BM' in x else 0)

# Extract gene expression data
gene_expression = data.iloc[:, 6:-2]
gene_names = gene_expression.columns

# Initialize lists to store results
p_values = []

# Perform linear regression for each gene
for gene in gene_names:
    X = sm.add_constant(pd.DataFrame({'AML': data['AML'], 'SampleSource': data['SampleSource']}))
    y = gene_expression[gene]
   
    if len(y) > 0:
        model = sm.OLS(y, X).fit()
        p_values.append(model.pvalues['AML'])
    else:
        p_values.append(1.0)

# Adjust p-values for multiple testing using FDR
reject, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

# Identify differentially expressed genes
diff_expressed_genes = gene_names[reject]
num_diff_expressed_genes = len(diff_expressed_genes)

# Find the gene with the smallest corrected p-value
most_diff_expressed_gene = gene_names[pvals_corrected.argmin()]

# Output the results
print(f"Number of differentially expressed genes: {num_diff_expressed_genes}")
print(f"Most differentially expressed gene: {most_diff_expressed_gene}")