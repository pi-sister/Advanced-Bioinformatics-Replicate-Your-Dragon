import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 6)
project_path = "/work/TALC/mdsc519_2026w/students/jamie/Dragon"

# Define columns of interest
columns_of_interest = [
    'Run', 'ReleaseDate', 'LoadDate', 'spots', 'bases', 'spots_with_mates', 
    'avgLength', 'size_MB','InsertDev', 'Platform', 'Model', 
    'SRAStudy', 'BioProject', 'Study_Pubmed_id', 'ProjectID', 'Sample', 
    'BioSample', 'SampleType', 'TaxID', 'ScientificName', 'SampleName', 
    'g1k_pop_code', 'source', 'g1k_analysis_group', 'Subject_ID', 'Sex', 
    'Disease', 'Tumor', 'Affection_Status', 'Analyte_Type', 'Histological_Type', 
    'Body_Site'
]

# Load metadata
metadata_path = Path(f"{project_path}/sra_data/metadata")
metadata_files = list(metadata_path.glob("*.csv")) + list(metadata_path.glob("*.tsv"))

if not metadata_files:
    print("No metadata files found in sra_data/metadata")
    exit(1)

# Load and combine metadata
dfs = []
for file in metadata_files:
    print(f"Loading {file.name}...")
    if file.suffix == '.csv':
        df = pd.read_csv(file)
    else:
        df = pd.read_csv(file, sep='\t')
    dfs.append(df)

metadata = pd.concat(dfs, ignore_index=True)

# Select only columns of interest (if they exist)
available_cols = [col for col in columns_of_interest if col in metadata.columns]
missing_cols = [col for col in columns_of_interest if col not in metadata.columns]

if missing_cols:
    print(f"\nWarning: {len(missing_cols)} requested columns not found: {missing_cols[:10]}")

metadata = metadata[available_cols]

print(f"\nLoaded {len(metadata)} records with {len(available_cols)} columns")
print(metadata.head())
print("\nData types:")
print(metadata.dtypes)
print("\nBasic statistics:")
print(metadata.describe())

# Create output directory
output_dir = Path("outputs")
output_dir.mkdir(exist_ok=True)

# Separate numeric and categorical columns
numeric_cols = metadata.select_dtypes(include=['int64', 'float64']).columns.tolist()
categorical_cols = metadata.select_dtypes(include=['object']).columns.tolist()

# 1. Distribution of Platform and LibraryStrategy
if 'Platform' in categorical_cols:
    plt.figure(figsize=(10, 5))
    metadata['Platform'].value_counts().plot(kind='barh')
    plt.title("Platform Distribution")
    plt.xlabel("Count")
    plt.tight_layout()
    plt.savefig(output_dir / "01_platform_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

if 'LibraryStrategy' in categorical_cols:
    plt.figure(figsize=(10, 5))
    metadata['LibraryStrategy'].value_counts().plot(kind='barh')
    plt.title("Library Strategy Distribution")
    plt.xlabel("Count")
    plt.tight_layout()
    plt.savefig(output_dir / "02_library_strategy.png", dpi=300, bbox_inches='tight')
    plt.close()

# 2. Numeric data distributions
if numeric_cols:
    fig, axes = plt.subplots(len(numeric_cols), 1, figsize=(12, 3*len(numeric_cols)))
    if len(numeric_cols) == 1:
        axes = [axes]
    
    for idx, col in enumerate(numeric_cols):
        axes[idx].hist(metadata[col].dropna(), bins=30, edgecolor='black', alpha=0.7)
        axes[idx].set_title(f"Distribution of {col}")
        axes[idx].set_xlabel(col)
        axes[idx].set_ylabel("Frequency")
    
    plt.tight_layout()
    plt.savefig(output_dir / "03_numeric_distributions.png", dpi=300, bbox_inches='tight')
    plt.close()

# 3. Correlation heatmap for numeric columns
if len(numeric_cols) > 1:
    plt.figure(figsize=(12, 10))
    corr_matrix = metadata[numeric_cols].corr()
    sns.heatmap(corr_matrix, annot=True, fmt='.2f', cmap='coolwarm', center=0)
    plt.title("Correlation Matrix of Numeric Columns")
    plt.tight_layout()
    plt.savefig(output_dir / "04_correlation_matrix.png", dpi=300, bbox_inches='tight')
    plt.close()

# 4. Size distribution (size_MB if available)
if 'size_MB' in numeric_cols:
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    metadata['size_MB'].hist(bins=50, edgecolor='black')
    plt.xlabel("Size (MB)")
    plt.ylabel("Frequency")
    plt.title("File Size Distribution")
    
    plt.subplot(1, 2, 2)
    metadata['size_MB'].plot(kind='box')
    plt.title("File Size Boxplot")
    plt.tight_layout()
    plt.savefig(output_dir / "05_size_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

# 5. If both Platform and Model exist, show relationship
if 'Platform' in categorical_cols and 'Model' in categorical_cols:
    plt.figure(figsize=(12, 6))
    platform_model = metadata.groupby(['Platform', 'Model']).size().unstack(fill_value=0)
    platform_model.plot(kind='bar', stacked=False, figsize=(12, 6))
    plt.title("Platform and Model Distribution")
    plt.xlabel("Platform")
    plt.ylabel("Count")
    plt.legend(title="Model", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / "06_platform_model.png", dpi=300, bbox_inches='tight')
    plt.close()

# 6. Save summary statistics
summary_stats = metadata.describe()
summary_stats.to_csv(output_dir / "summary_statistics.csv")

# 7. Save value counts for categorical columns
for col in ['BioProject', 'Platform', 'LibraryStrategy', 'LibraryLayout']:
    if col in categorical_cols:
        counts = metadata[col].value_counts()
        counts.to_csv(output_dir / f"value_counts_{col}.csv")

# 8. Overall summary report
with open(output_dir / "analysis_report.txt", 'w') as f:
    f.write(f"Metadata Analysis Report\n")
    f.write(f"========================\n")
    f.write(f"Total Records: {len(metadata)}\n")
    f.write(f"Columns Analyzed: {len(available_cols)}\n")
    f.write(f"Numeric Columns: {len(numeric_cols)}\n")
    f.write(f"Categorical Columns: {len(categorical_cols)}\n\n")
    
    if 'size_MB' in numeric_cols:
        f.write(f"Total Size: {metadata['size_MB'].sum():.2f} MB ({metadata['size_MB'].sum()/1024:.2f} GB)\n")
        f.write(f"Average File Size: {metadata['size_MB'].mean():.2f} MB\n")
    
    if 'BioProject' in categorical_cols:
        f.write(f"\nBioProjects: {metadata['BioProject'].nunique()}\n")
        f.write(str(metadata['BioProject'].value_counts()))

print(f"\nAnalysis complete! Graphs and reports saved to {output_dir}/")
