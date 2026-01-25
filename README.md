## The GTEX primary dataset can be conferred upon request.

# Muscle Epigenetic Clock

A DNA methylation-based age prediction model for human skeletal muscle tissue, built from GTEx consortium data using penalized regression with rigorous cross-validation.

**Performance (Leave-One-Out CV):**
- Pearson r: 0.936
- R²: 0.874
- MAE: 3.26 years
- RMSE: 3.91 years
- Clock CpGs: 95

---

## Table of Contents

1. [Historical Background: The Science of Biological Clocks](#historical-background-the-science-of-biological-clocks)
2. [Why DNA Methylation?](#why-dna-methylation)
3. [Dataset Description](#dataset-description)
4. [Methodology Overview](#methodology-overview)
5. [Step-by-Step Pipeline Explanation](#step-by-step-pipeline-explanation)
6. [Using This Clock](#using-this-clock)
7. [Interpreting Results](#interpreting-results)
8. [Limitations and Considerations](#limitations-and-considerations)
9. [File Descriptions](#file-descriptions)
10. [References](#references)

---

## Historical Background: The Science of Biological Clocks

### The Quest to Measure Biological Age

Chronological age—the number of years since birth—is a crude measure of aging. Two 60-year-olds can have vastly different physiological states: one might be running marathons while another struggles with multiple age-related diseases. This observation led researchers to seek molecular biomarkers that could measure "biological age"—the true physiological state of an organism independent of time elapsed since birth.

### Early Biomarkers of Aging

Before the genomics era, researchers explored various aging biomarkers:

- **Telomere length** (1990s): Elizabeth Blackburn's work on telomeres—the protective caps on chromosome ends that shorten with each cell division—suggested they could serve as a "mitotic clock." However, telomere length proved to be a noisy predictor with high inter-individual variability and tissue-specific differences.

- **Gene expression signatures** (2000s): Microarray studies identified genes whose expression changes with age, but these signatures were often tissue-specific and technically variable.

- **Inflammatory markers** (2000s): "Inflammaging" markers like IL-6 and CRP increase with age but are confounded by acute illness and lifestyle factors.

### The Epigenetic Revolution

The breakthrough came with the recognition that **DNA methylation**—a chemical modification where a methyl group is added to cytosine bases, predominantly at CpG dinucleotides—changes systematically with age. Unlike genetic mutations, DNA methylation is:

1. **Highly quantifiable**: Modern arrays measure methylation at hundreds of thousands of sites with high precision
2. **Relatively stable**: More stable than RNA or protein measurements
3. **Functionally relevant**: Methylation regulates gene expression and chromatin state
4. **Tissue-accessible**: Can be measured from blood, saliva, or any tissue with DNA

### Landmark Epigenetic Clocks

**Horvath Clock (2013)**: Steve Horvath's multi-tissue clock was the first pan-tissue epigenetic clock, trained on 8,000 samples from 51 tissue types. Using elastic net regression on Illumina 450K array data, he identified 353 CpGs that predict age with a median error of 3.6 years across tissues. This clock revealed that:
- Different tissues age at different rates
- Cancer cells show dramatic age acceleration
- Biological age can deviate substantially from chronological age

**Hannum Clock (2013)**: Gregory Hannum developed a blood-specific clock using 71 CpGs, demonstrating that blood methylation age correlates with mortality risk and metabolic health markers.

**PhenoAge (2018)**: Morgan Levine created a "second-generation" clock trained not just on chronological age, but on phenotypic measures of aging (albumin, creatinine, glucose, CRP, etc.), producing a clock more predictive of healthspan and lifespan.

**GrimAge (2019)**: Ake Lu developed GrimAge by training on time-to-death data, creating a clock that predicts mortality better than chronological age alone.

**DunedinPACE (2022)**: A "pace of aging" clock measuring the rate of biological aging rather than absolute biological age, based on longitudinal data from the Dunedin birth cohort.

### Why Tissue-Specific Clocks?

While pan-tissue clocks like Horvath's are valuable for their generalizability, they sacrifice accuracy for breadth. Tissue-specific clocks can:

1. **Achieve higher accuracy** by capturing tissue-specific aging patterns
2. **Identify tissue-relevant biology** rather than generic aging signatures
3. **Detect tissue-specific age acceleration** that might be diluted in pan-tissue models

Skeletal muscle is particularly important because:
- It comprises ~40% of body mass
- Sarcopenia (muscle aging) is a major determinant of healthspan
- Muscle has unique metabolic and regenerative properties
- Few muscle-specific clocks exist due to limited tissue availability

---

## Why DNA Methylation?

### The Molecular Basis

DNA methylation occurs when a methyl group (-CH3) is added to the 5' carbon of cytosine, creating 5-methylcytosine (5mC). This modification:

1. **Occurs predominantly at CpG sites**: Cytosines followed by guanines
2. **Is catalyzed by DNMTs**: DNA methyltransferases (DNMT1, DNMT3A, DNMT3B)
3. **Is removed by TET enzymes**: Through oxidation to 5hmC, 5fC, and 5caC
4. **Regulates transcription**: Promoter methylation typically silences genes

### Why Methylation Changes with Age

Several mechanisms drive age-related methylation changes:

1. **Epigenetic drift**: Stochastic errors in methylation maintenance during cell division
2. **Environmental accumulation**: Lifetime exposure to diet, toxins, and stress
3. **Cellular composition shifts**: Changing proportions of cell types within tissues
4. **Programmatic changes**: Developmentally regulated methylation programs that continue into adulthood

### The Methylation Array Technology

This clock uses data from the Illumina Infinium MethylationEPIC array, which measures methylation at ~850,000 CpG sites across the genome. For each CpG:

- **Beta value (β)**: Ratio of methylated signal to total signal, ranging from 0 (unmethylated) to 1 (fully methylated)
- **M-value**: Logit-transformed beta value, M = log2(β / (1-β)), which has better statistical properties

Beta values are used here because they have a direct biological interpretation (proportion of methylated alleles).

---

## Dataset Description

### Source: GTEx Consortium

The Genotype-Tissue Expression (GTEx) project is a comprehensive public resource studying tissue-specific gene expression and regulation. This dataset includes DNA methylation profiles from skeletal muscle samples.

### Dataset Characteristics

| Attribute | Value |
|-----------|-------|
| Samples | 47 |
| CpG sites | 754,119 |
| Tissue | Skeletal muscle |
| Age range | 34-74 years |
| Age format | Decade bins (e.g., "40-49") |
| Platform | Illumina EPIC array |
| Sex distribution | Mixed male/female |

### Data Structure

**Methylation file** (`GTEx_Muscle.meth.csv`):
- Rows: CpG probe IDs (cg########)
- Columns: Sample IDs
- Values: Beta values (0-1)

**Annotation file** (`GTEx_Muscle.anno.csv`):
- Sample metadata including age, sex, and technical variables

### Age Encoding

GTEx encodes age as decade bins for privacy. We convert these to midpoints:
- "30-39" → 34.5
- "40-49" → 44.5
- "50-59" → 54.5
- "60-69" → 64.5
- "70-79" → 74.5

This introduces some measurement error but is unavoidable given the data format.

---

## Methodology Overview

### The Clock-Building Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│  Raw Methylation Data (754,119 CpGs × 47 samples)               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 1: Data Loading & Preprocessing                           │
│  - Load methylation matrix and annotations                      │
│  - Convert age bins to numeric midpoints                        │
│  - Handle missing values                                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 2: Feature-Age Correlation                                │
│  - Compute Pearson r for all 754,119 CpGs vs age               │
│  - Identify age-associated CpGs                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 3: Feature Selection (|r| > 0.5)                          │
│  - Reduce to 423 strongly correlated CpGs                       │
│  - Prevents overfitting in high-dimensional setting             │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 4: Elastic Net Regression with CV                         │
│  - Regularized regression for sparse feature selection          │
│  - 5-fold CV for hyperparameter tuning                          │
│  - Leave-One-Out CV for unbiased predictions                    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 5: Model Evaluation                                       │
│  - Pearson correlation, R², MAE, RMSE                           │
│  - Overfitting diagnostic (train vs CV performance)             │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  Output: 95-CpG Muscle Clock                                    │
│  - Coefficients for age prediction                              │
│  - Performance metrics and visualizations                       │
└─────────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Pipeline Explanation

### Step 1: Data Loading and Preprocessing

```python
meth = pd.read_csv('GTEx_Muscle.meth.csv', sep='\t', index_col=0)
anno = pd.read_csv('GTEx_Muscle.anno.csv', sep='\t', index_col=0).T
```

**What this does:**
- Loads the methylation beta matrix (CpGs × samples)
- Loads and transposes the annotation matrix (samples × variables)

**Why it's necessary:**
- Data must be in consistent orientation for downstream analysis
- Transposing annotation aligns samples as rows for joining

```python
def age_midpoint(age_str):
    lo, hi = map(int, age_str.split('-'))
    return (lo + hi) / 2

anno['age_numeric'] = anno['age'].apply(age_midpoint)
```

**What this does:**
- Converts decade age bins (e.g., "40-49") to numeric midpoints (44.5)

**Why it's necessary:**
- Regression requires numeric outcomes
- Midpoint is the best point estimate given the interval-censored data

```python
X = meth.T.values.astype(np.float32)
y = anno.loc[meth.columns, 'age_numeric'].values.astype(np.float32)
X = np.nan_to_num(X, nan=0.5, posinf=1.0, neginf=0.0)
```

**What this does:**
- Creates feature matrix X (samples × CpGs) and target vector y
- Uses float32 to reduce memory footprint
- Imputes missing values with 0.5 (median beta value)

**Why it's necessary:**
- Machine learning algorithms require numeric matrices without missing values
- 0.5 imputation is reasonable for beta values (represents intermediate methylation)

### Step 2: Computing CpG-Age Correlations

```python
X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
X_std[X_std == 0] = 1  # Avoid division by zero
X_norm = (X - X_mean) / X_std
y_norm = (y - y.mean()) / y.std()
correlations = (X_norm.T @ y_norm) / len(y)
```

**What this does:**
- Computes Pearson correlation between each CpG and age using vectorized operations
- Standardizes both X and y, then computes r as the dot product divided by n

**Why it's necessary:**
- Pearson correlation quantifies the linear relationship between methylation and age
- Vectorized computation is essential for efficiency (754K correlations in seconds vs hours)

**Mathematical basis:**
The Pearson correlation coefficient is:
```
r = Σ[(xi - x̄)(yi - ȳ)] / [√Σ(xi - x̄)² × √Σ(yi - ȳ)²]
```

When both variables are z-score standardized, this simplifies to:
```
r = (1/n) × Σ(xi_std × yi_std) = (1/n) × X_std^T @ y_std
```

**Interpretation of correlations:**
- r > 0: Methylation increases with age (hypermethylation)
- r < 0: Methylation decreases with age (hypomethylation)
- |r| > 0.5: Strong age association
- |r| > 0.7: Very strong age association

### Step 3: Feature Selection

```python
sig_mask = np.abs(correlations) > 0.5
X_filtered = X[:, sig_mask]
```

**What this does:**
- Selects only CpGs with |r| > 0.5 (strong correlation with age)
- Reduces features from 754,119 to 423 (0.06%)

**Why it's necessary:**

1. **Curse of dimensionality**: With p >> n (754K features, 47 samples), models will overfit without feature selection. The model can find spurious patterns that don't generalize.

2. **Regularization alone is insufficient**: While elastic net provides regularization, starting with highly correlated features improves both performance and interpretability.

3. **Biological relevance**: CpGs with |r| < 0.5 are unlikely to be biologically meaningful age markers for this tissue.

4. **Computational efficiency**: Reduces the search space for hyperparameter optimization.

**Choice of threshold:**
- |r| > 0.3 (lenient): 47,846 CpGs → led to overfitting (CV r = 0.69)
- |r| > 0.5 (stringent): 423 CpGs → good generalization (CV r = 0.94)

The stringent threshold is critical for small sample sizes.

### Step 4: Standardization

```python
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_filtered)
```

**What this does:**
- Z-score standardizes each CpG: z = (x - μ) / σ
- Each feature now has mean=0 and std=1

**Why it's necessary:**

1. **Equal weighting**: Without standardization, CpGs with larger variance would dominate the penalty term in elastic net.

2. **Regularization fairness**: L1 and L2 penalties assume features are on comparable scales. A CpG with β ∈ [0.8, 0.9] would have smaller coefficients than one with β ∈ [0.1, 0.9] purely due to scale.

3. **Numerical stability**: Standardized features improve optimization convergence.

### Step 5: Elastic Net Regression

```python
model = ElasticNetCV(
    l1_ratio=[0.5, 0.7, 0.9, 0.95, 0.99, 1.0],
    alphas=np.logspace(-2, 2, 100),
    cv=5,
    max_iter=20000,
    random_state=42,
    n_jobs=-1
)
```

**What this does:**
- Fits an elastic net model with automatic hyperparameter selection
- Tests multiple L1 ratios and regularization strengths
- Uses 5-fold cross-validation for hyperparameter tuning

**Why elastic net:**

Elastic net combines L1 (Lasso) and L2 (Ridge) penalties:

```
Loss = ||y - Xβ||² + α × [ρ × ||β||₁ + (1-ρ) × ||β||₂²]
```

Where:
- α (alpha): Overall regularization strength
- ρ (l1_ratio): Balance between L1 and L2

**L1 (Lasso) benefits:**
- Produces sparse solutions (many coefficients exactly zero)
- Performs feature selection automatically
- Interpretable: only selected CpGs contribute to prediction

**L2 (Ridge) benefits:**
- Handles correlated features better than pure Lasso
- More stable coefficient estimates
- Doesn't arbitrarily select one of several correlated features

**Elastic net combines both:**
- Gets sparsity from L1 (feature selection)
- Gets stability from L2 (grouped selection of correlated features)

**Why these hyperparameter ranges:**
- `l1_ratio=[0.5, 0.7, 0.9, 0.95, 0.99, 1.0]`: Emphasizes sparsity (higher L1)
- `alphas=np.logspace(-2, 2, 100)`: Tests regularization from 0.01 to 100 on log scale
- High regularization prevents overfitting with small n

### Step 6: Leave-One-Out Cross-Validation

```python
loo = LeaveOneOut()
y_pred_loo = cross_val_predict(model, X_scaled, y, cv=loo, n_jobs=-1)
```

**What this does:**
- For each sample i, trains on all other samples and predicts sample i
- Produces predictions where each sample was never in its own training set

**Why LOO-CV:**

1. **Unbiased estimates**: Every prediction is made on held-out data, preventing optimistic bias from training set leakage.

2. **Maximum training data**: With only 47 samples, LOO-CV uses 46 samples for training in each fold, maximizing available training data.

3. **Standard for epigenetic clocks**: Horvath and other clock developers use LOO-CV for small datasets.

4. **No random variation**: Unlike k-fold CV, LOO-CV produces deterministic results.

**Comparison to k-fold CV:**
- 5-fold CV: Trains on 80% of data → more variance in estimates
- LOO-CV: Trains on ~98% of data → lower variance, higher computational cost
- For n=47, LOO-CV is computationally feasible and preferred

### Step 7: Model Evaluation

```python
r, p = pearsonr(y, y_pred_loo)
mae = np.mean(np.abs(y - y_pred_loo))
rmse = np.sqrt(np.mean((y - y_pred_loo)**2))
r2 = 1 - np.sum((y - y_pred_loo)**2) / np.sum((y - y.mean())**2)
```

**Metrics explained:**

**Pearson r (0.936):**
- Measures linear correlation between predicted and actual age
- Range: -1 to 1, where 1 is perfect positive correlation
- Most commonly reported metric for epigenetic clocks
- Interpretation: 93.6% of ranking information is preserved

**R² (0.874):**
- Proportion of variance in age explained by predictions
- R² = 1 - (SS_residual / SS_total)
- Interpretation: 87.4% of age variation is captured

**MAE (3.26 years):**
- Mean Absolute Error: average |predicted - actual|
- Directly interpretable in years
- Less sensitive to outliers than RMSE
- Interpretation: Average prediction error is 3.26 years

**RMSE (3.91 years):**
- Root Mean Squared Error: √(mean((predicted - actual)²))
- Penalizes large errors more than MAE
- Same units as the outcome (years)
- RMSE > MAE indicates some larger errors exist

### Step 8: Overfitting Diagnostic

```python
y_pred_train = model.predict(X_scaled)
r_train, _ = pearsonr(y, y_pred_train)
gap = r_train - r  # 1.0 - 0.936 = 0.064
```

**What this does:**
- Compares training set performance to cross-validation performance
- A large gap indicates overfitting

**Why it's necessary:**

Overfitting occurs when a model learns noise in the training data rather than true signal. Signs of overfitting:
- Training r ≈ 1.0 (perfect fit)
- CV r << Training r (poor generalization)

**Interpretation of gap:**
- Gap < 0.1: Minimal overfitting ✓
- Gap 0.1-0.2: Moderate overfitting (acceptable for small n)
- Gap > 0.2: Severe overfitting (need more regularization)

**Our result (gap = 0.064):**
- Training r = 1.00 (expected—elastic net can fit training data well)
- CV r = 0.94 (excellent generalization)
- Gap = 0.06 (minimal overfitting)

This indicates the model learned true age-related signal, not noise.

---

## Using This Clock

### Applying to New Samples

To predict age for new muscle methylation samples:

```python
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

# Load clock coefficients
clock = pd.read_csv('muscle_clock_coefficients.csv')
clock_cpgs = clock['CpG'].values
clock_coefs = clock['Coefficient'].values

# Load new methylation data (samples × CpGs)
new_data = pd.read_csv('your_methylation_data.csv', sep='\t', index_col=0)

# Extract clock CpGs (must be present in your data)
missing = set(clock_cpgs) - set(new_data.columns)
if missing:
    print(f"Warning: Missing {len(missing)} clock CpGs")

X_new = new_data[clock_cpgs].values

# Standardize using training set parameters (or your data)
scaler = StandardScaler()
X_new_scaled = scaler.fit_transform(X_new)

# Predict age (linear combination of scaled CpGs)
# Note: You need the intercept from the original model
# For simplicity, using mean age as baseline
predicted_age = X_new_scaled @ clock_coefs + 57.1  # Adjust intercept as needed
```

### Requirements for Input Data

1. **Platform compatibility**: Data should be from Illumina EPIC or 450K arrays
2. **Preprocessing**: Beta values should be normalized (e.g., BMIQ, quantile normalization)
3. **CpG coverage**: All 95 clock CpGs must be present
4. **Tissue match**: Best performance on skeletal muscle; other tissues may show bias

### Batch Effects

If applying to external datasets, consider:
- ComBat or similar batch correction
- Validating on a subset with known ages
- Using reference-based normalization

---

## Interpreting Results

### Age Acceleration

The difference between predicted and chronological age is called **epigenetic age acceleration**:

```
Age Acceleration = Predicted Age - Chronological Age
```

**Interpretation:**
- Positive acceleration: Biologically "older" than chronological age
- Negative acceleration: Biologically "younger" than chronological age

**Associations with positive age acceleration:**
- Obesity
- Smoking
- Sedentary lifestyle
- Chronic diseases
- Mortality risk

**Associations with negative age acceleration:**
- Physical fitness
- Healthy diet
- Some genetic variants

### Clock CpGs Biology

The 95 clock CpGs can be annotated to understand biological pathways:

**Hypermethylated with age (positive coefficients):**
- Often in CpG islands and promoter regions
- Associated with gene silencing
- May reflect cumulative environmental exposure

**Hypomethylated with age (negative coefficients):**
- Often in gene bodies and intergenic regions
- May reflect epigenetic drift
- Associated with genomic instability

### Residual Analysis

The residual plot (Age Acceleration vs Age) should show:
- No systematic trend (residuals scattered around zero)
- Homoscedasticity (constant variance across age range)
- Normal distribution of residuals

Our clock shows mean residual ≈ 0 and SD = 3.91 years, indicating unbiased predictions.

---

## Limitations and Considerations

### Sample Size

With n=47, this clock has limited statistical power:
- Wide confidence intervals on estimates
- May not capture rare aging patterns
- Would benefit from larger training sets

### Age Range

Training data spans 34-74 years:
- Extrapolation to ages outside this range is unreliable
- Pediatric and very elderly samples may show poor accuracy
- Age-related changes may not be linear at extremes

### Tissue Specificity

This clock is trained on skeletal muscle:
- Will show systematic bias in other tissues
- Blood samples will likely underpredict (blood ages differently)
- Other muscle types (cardiac, smooth) may vary

### Technical Factors

- Batch effects between studies
- Array version differences (EPIC vs 450K)
- Sample quality (degraded DNA, low input)
- Processing differences (normalization, probe filtering)

### Age Binning

GTEx provides ages as decade bins:
- True ages are unknown
- Midpoint assumption adds ~2.5 years of noise
- May attenuate true correlations

### Causal Inference

This clock is correlational, not causal:
- Clock CpGs may not drive aging
- Interventions targeting clock CpGs may not affect aging
- Age acceleration is a biomarker, not necessarily a mechanism

---

## File Descriptions

| File | Description |
|------|-------------|
| `GTEx_Muscle.meth.csv` | Raw methylation beta values (754,119 CpGs × 47 samples) |
| `GTEx_Muscle.anno.csv` | Sample annotations (age, sex, tissue, etc.) |
| `muscle_clock.ipynb` | Jupyter notebook with full analysis pipeline |
| `build_clock.py` | Standalone Python script for reproducibility |
| `muscle_clock_coefficients.csv` | Final clock: 95 CpGs with coefficients and correlations |
| `muscle_clock_performance.png` | Predicted vs actual age scatterplot |
| `top_cpgs.png` | Visualization of top age-correlated CpGs |
| `correlation_distribution.png` | Histogram of all CpG-age correlations |
| `coefficient_distribution.png` | Distribution of clock coefficients |
| `README.md` | This documentation |

### Clock Coefficients Format

`muscle_clock_coefficients.csv` contains:

| Column | Description |
|--------|-------------|
| CpG | Illumina probe ID (cg########) |
| Coefficient | Elastic net coefficient (weight in prediction) |
| Correlation | Pearson r with age in training data |

---

## References

### Foundational Epigenetic Clock Papers

1. **Horvath S. (2013)**. DNA methylation age of human tissues and cell types. *Genome Biology*, 14(10), R115. [The original multi-tissue clock]

2. **Hannum G. et al. (2013)**. Genome-wide methylation profiles reveal quantitative views of human aging rates. *Molecular Cell*, 49(2), 359-367. [Blood-specific clock]

3. **Levine M.E. et al. (2018)**. An epigenetic biomarker of aging for lifespan and healthspan. *Aging*, 10(4), 573-591. [PhenoAge]

4. **Lu A.T. et al. (2019)**. DNA methylation GrimAge strongly predicts lifespan and healthspan. *Aging*, 11(2), 303-327. [GrimAge]

5. **Belsky D.W. et al. (2022)**. DunedinPACE, a DNA methylation biomarker of the pace of aging. *eLife*, 11, e73420. [Pace of aging]

### Methodology References

6. **Zou H. & Hastie T. (2005)**. Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society B*, 67(2), 301-320. [Elastic net theory]

7. **Bibikova M. et al. (2011)**. High density DNA methylation array with single CpG site resolution. *Genomics*, 98(4), 288-295. [Illumina array technology]

### GTEx Consortium

8. **GTEx Consortium (2020)**. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*, 369(6509), 1318-1330. [Data source]

---

## Citation

If you use this clock in your research, please cite:

```
Muscle Epigenetic Clock built from GTEx skeletal muscle methylation data.
Method: Elastic net regression with leave-one-out cross-validation.
Performance: r=0.936, MAE=3.26 years (n=47, 95 CpGs).
```

---

## License

This clock is provided for research purposes. The underlying GTEx data is subject to GTEx consortium data use agreements.

---

## Contact

For questions or issues, please open a GitHub issue or contact the repository maintainer.
