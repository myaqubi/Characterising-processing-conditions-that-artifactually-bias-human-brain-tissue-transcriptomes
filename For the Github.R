### Normalizin the read counst using DESEq2 package
library(DESeq2)
count_data <- read.csv("your_count_data.csv", row.names = 1)  # Replace with your file path
col_data <- read.csv("your_sample_info.csv", row.names = 1)  # Replace with your file path
all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition) 
dds <- estimateSizeFactors(dds)  
normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)
write.csv(normalized_counts, "normalized_counts.csv")


### Using the normlaized values for the heatmap and the PCA

### PCA
normalized_counts <- read.csv("normalized_counts.csv")
row.names(counts) = counts[[1]] ## makes rownames from 1st column
counts <- counts[, -1]
pca <- prcomp(t(counts), scale=TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1], 
                       Y=pca$x[,2])

### Heatmap
heatmap.2(normalized_counts,
          scale = "row",                   # Scale by row (gene-wise Z-scores)
          dendrogram = "both",             # Display dendrogram for both rows and columns
          trace = "none",                  # No trace lines inside the heatmap
          col = bluered(100),              # Color scheme
          margins = c(8, 8),               # Margins for the plot
          cexRow = 0.5,                    # Size of row labels
          cexCol = 0.7,                    # Size of column labels
          key = TRUE,                      # Show color key
          key.title = "Expression",        # Title of the color key
          key.xlab = "Z-score",            # Label of the color key
          main = "Heatmap of Normalized Counts")  # Main title

### Differential expression analysis between groups:
library(DESeq2)
count_data <- read.csv("your_count_data.csv", row.names = 1)  # Replace with your file path
col_data <- read.csv("your_sample_info.csv", row.names = 1)  # Replace with your file path
all(colnames(count_data) %in% rownames(col_data))
all(colnames(count_data) == rownames(col_data))
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)  # Replace 'condition' with your actual column name
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))
res <- res[order(res$pvalue), ]
head(res)
write.csv(as.data.frame(res), "deseq2_results_LMPI_vs_control.csv")
resSig <- subset(res, padj < 0.05)
head(resSig)
write.csv(as.data.frame(resSig), "deseq2_significant_genes_treatment_vs_control.csv")


### Using the DESeq2 outoput file for the GO analysis

library(clusterProfiler)
counts_data <- read.csv("deseq2_significant_genes_treatment_vs_control.csv", check.names = FALSE)
row.names(counts_data) = counts_data[[1]] 
counts_data <- counts_data[, -1]

### GO for the up reg genes
sigs <- counts_data[counts_data$log2FoldChange > 1,]
up_reg_genes <- row.names(sigs)
Go_results_highPMI_vs_control <- enrichGO(gene = up_reg_genes, OrgDb = 'org.Hs.eg.db',
                                          keyType = "SYMBOL", ont = "BP")
tiff('High_PMI_final_samples_vs_surgery_control_up_reg_GO.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
barplot(Go_results_highPMI_vs_control, showCategory = 20)
dev.off()

### GO for the down reg genes
sigs <- counts_data[counts_data$log2FoldChange < 1,]
down_reg_genes <- row.names(sigs)
Go_results_highPMI_vs_control_downreg <- enrichGO(gene = down_reg_genes, OrgDb = 'org.Hs.eg.db',
                                                  keyType = "SYMBOL", ont = "BP")
as.data.frame(Go_results_highPMI_vs_control_downreg)
write.csv(Go_results_highPMI_vs_control_downreg@result, file = "Go_results_highPMI_vs_control_downreg_genes.csv")
tiff('High_PMI_final_samples_vs_surgery_control_down_reg_GO.tiff', units="in", width=10, height=2, res=300, compression = 'lzw')
barplot(Go_results_highPMI_vs_control_downreg, showCategory = 20)
dev.off()

### Single nuclear RNAseq data analysis
The package (scFlow) that was used to analyse the single nuclear data https://github.com/neurogenomics/scFlow


### Using the results of the DESeq2 DEG analysis or the volcano plot in python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the data from the Excel file
file_path = 'DESeq2_output.xlsx'   # Replace with your actual file path
df = pd.read_excel(file_path)

# Assuming the columns are named 'Gene Symbol', 'log2FoldChange', and 'padj'
df.columns = ['Gene Symbol', 'log2FoldChange', 'padj']

# Create a new column for coloring based on criteria
def color_genes(row):
  if row['padj'] < 0.05 and (row['log2FoldChange'] > 1 or row['log2FoldChange'] < -1):
  return 'red'
elif row['padj'] < 0.05:
  return 'blue'
else:
  return 'grey'

df['color'] = df.apply(color_genes, axis=1)

# Create the volcano plot
plt.figure(figsize=(10, 8))
sns.scatterplot(data=df, x='log2FoldChange', y=-np.log10(df['padj']), hue='color', palette={'red':'red', 'blue':'blue', 'grey':'grey'}, legend=False)
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(Adjusted P-value)')
plt.title('Volcano Plot')

# Add a horizontal line at y = -log10(0.05) to denote the p-value threshold
plt.axhline(y=-np.log10(0.05), color='grey', linestyle='--')

# Add vertical lines at x = -1 and x = 1 to denote the fold change thresholds
plt.axvline(x=-1, color='grey', linestyle='--')
plt.axvline(x=1, color='grey', linestyle='--')

# Add a legend in the upper left corner
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', label='p adj < 0.05 & |log2FC| > 1', markerfacecolor='red', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='p adj < 0.05', markerfacecolor='blue', markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='grey', markersize=10)]
plt.legend(handles=legend_elements, loc='upper left')

# Save the plot with high resolution
save_path = 'volcano_plot.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

### Single nucleus RNAseq data analysis
The package (scFlow) that was used to analyse the single nuclear data https://github.com/neurogenomics/scFlow

### Machine learning code
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from tensorflow.keras.layers import Dense, Input, Dropout
from tensorflow.keras.models import Model
import numpy as np
import tensorflow as tf
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt

# Load and preprocess the data
Moein_New = pd.read_csv('C:/input_file.csv', header=0)
Moein_New = Moein_New.transpose()

new_header = Moein_New.iloc[0]  # Get the new header from the first row
Moein_New = Moein_New[1:]  # Remove the first row from the data
Moein_New.columns = new_header  # Set the new header

X = Moein_New.drop('Sample', axis=1)  # Features
y = Moein_New['Sample']  # Target variable

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

y_labels = pd.factorize(y)[0]

rf = RandomForestClassifier(n_estimators=100)
rf.fit(X_scaled, y_labels)
importances = rf.feature_importances_

# Define neural network architecture
def build_model(input_shape):
  input_layer = Input(shape=(input_shape,))
x = Dense(2600, activation='relu')(input_layer)
x = Dropout(0.6)(x)
x = Dense(400, activation='relu')(x)
x = Dropout(0.6)(x)
x = Dense(40, activation='relu')(x)
output_layer = Dense(len(np.unique(y_labels)), activation='softmax')(x)
model = Model(inputs=input_layer, outputs=output_layer)
return model

epochs = 50
batch_size = 64
learning_rate = 0.001

# Function to train, evaluate and save model and features
def train_and_evaluate(X_selected, y_labels, y, feature_count, X_full):
  kfold = StratifiedKFold(n_splits=7, shuffle=True, random_state=64)
best_val_loss = float('inf')
best_model = None
best_model_info = {}

for train_idx, test_idx in kfold.split(X_selected, y_labels):
  X_train, X_test = X_selected[train_idx], X_selected[test_idx]
y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
y_train = pd.get_dummies(y_train)
y_test = pd.get_dummies(y_test)

model = build_model(X_train.shape[1])
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate),
              loss='categorical_crossentropy', metrics=['accuracy'])

history = model.fit(X_train, y_train, epochs=epochs, batch_size=batch_size,
                    validation_split=0.2, verbose=1)

val_loss = history.history['val_loss'][-1]

if val_loss < best_val_loss:
  best_val_loss = val_loss
best_model = model
best_model_info = {
  'Feature_Count': feature_count,
  'Validation_Loss': best_val_loss,
  'Validation_Accuracy': history.history['val_accuracy'][-1],
  'Neural_Architecture': [2600, 400, 40]
}

# Save model
model_save_path = f"C:/best_model_feature_{feature_count}.keras"
best_model.save(model_save_path)

# Save features
selected_features = X_full.columns[np.argsort(importances)[-feature_count:]]
features_save_path = f"C:/selected_features_{feature_count}.csv"
selected_features.to_series().to_csv(features_save_path, index=False, header=False)

return best_model, best_model_info, X_test, y_test

# Feature indices
indices_5000 = np.argsort(importances)[-5000:]
indices_900 = np.argsort(importances)[-900:]
indices_2000 = np.argsort(importances)[-2000:]
indices_3000 = np.argsort(importances)[-3000:]
indices_4000 = np.argsort(importances)[-4000:]

# Feature matrices
X_selected_5000 = X_scaled[:, indices_5000]
X_selected_900 = X_scaled[:, indices_900]
X_selected_2000 = X_scaled[:, indices_2000]
X_selected_3000 = X_scaled[:, indices_3000]
X_selected_4000 = X_scaled[:, indices_4000]

# Train, evaluate, and save models and features
best_model_5000, best_model_info_5000, X_test_5000, y_test_5000 = train_and_evaluate(X_selected_5000, y_labels, y, 5000, X)
best_model_900, best_model_info_900, X_test_900, y_test_900 = train_and_evaluate(X_selected_900, y_labels, y, 900, X)
best_model_2000, best_model_info_2000, X_test_2000, y_test_2000 = train_and_evaluate(X_selected_2000, y_labels, y, 2000, X)
best_model_3000, best_model_info_3000, X_test_3000, y_test_3000 = train_and_evaluate(X_selected_3000, y_labels, y, 3000, X)
best_model_4000, best_model_info_4000, X_test_4000, y_test_4000 = train_and_evaluate(X_selected_4000, y_labels, y, 4000, X)
