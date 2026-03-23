# Building a Skeletal Muscle Epigenetic Clock: What I Learned (and What Went Wrong)

A tissue-specific age predictor from 750K CpG methylation sites across 47 vastus lateralis biopsies (GTEx EPIC 850K array). This is a walkthrough of the analysis in `muscle_clock.ipynb`, written for a technical audience familiar with methylation arrays and regularized regression.

---

## The Setup

The goal was straightforward: take beta values (0–1 methylation proportions) from ~754K CpG sites across 47 skeletal muscle samples and build an ElasticNet regression model that predicts chronological age. Standard epigenetic clock methodology — univariate pre-filter to reduce dimensionality, then penalized regression with leave-one-out cross-validation.

The initial results looked excellent: **r = 0.94, MAE = 3.3 years, 95 clock CpGs**.

Then I checked the methodology properly, and the real numbers fell out.

## The Bug Everyone Ships

The original pipeline computed CpG-age correlations on all 47 samples, filtered to |r| > 0.5, scaled the data, and *then* ran LOO-CV. This is textbook feature selection leakage — the held-out sample in each fold had already influenced which features entered the model. It's described in Ambroise & McLachlan (2002) and Simon et al. (2003), and it is *extremely common* in published epigenetic clock papers.

The fix: move everything inside the CV loop. For each of the 47 folds, compute correlations on the 46 training samples only, apply the filter, fit the scaler, tune ElasticNetCV via inner 5-fold CV, then predict the single held-out sample. No information leaks.

**The corrected numbers: r = 0.54, MAE = 7.6 years.**

That's a delta of +0.40 in r and -4.3 years in MAE — entirely attributable to a single methodological error. The side-by-side comparison in the notebook makes this viscerally clear.

## What the 3D PCA Tells Us

PCA on the top 1,000 most variable CpGs shows that **age is not the dominant source of methylation variance in skeletal muscle**. PC1 through PC3 collectively explain ~17% of the variance, and none of them correlate significantly with age (all p > 0.19). The scatter in 3D PC space shows no obvious age gradient — the color mapping is essentially random.

This matters because it sets expectations: the clock is trying to extract a weak, distributed signal from a background dominated by inter-individual variation, sex (28M/19F, no correction applied), and Sentrix batch effects. A corrected r of 0.54 is actually reasonable given this reality.

## What the Clock Did Find

Despite the honest performance being modest, there is real biology here:

- **FHL2** — one of the 108 genes mapped from the 95 clock CpGs — is a canonical Horvath clock gene involved in muscle-specific Wnt signaling. Its presence in a tissue-specific muscle clock built from scratch is reassuring.

- The CpG island distribution of clock sites is enriched for **Open Sea regions** relative to the EPIC array background. This is consistent with the literature: age-associated methylation changes preferentially occur at intergenic and enhancer regions rather than promoter islands.

- SHAP decomposition (exact, not approximate, because ElasticNet is linear) shows clean **directional reversal** between the youngest and oldest samples across the top CpGs. The model learned a genuine age axis, not noise — it just doesn't predict very accurately.

- GO enrichment pointed toward intracellular transport, GTPase signaling, and cell cycle regulation. None survived FDR correction (108 genes is underpowered for enrichment), but the themes are biologically plausible for muscle homeostasis.

## The Hard Limits

**N = 47 is not enough.** The feature count per CV fold ranged from 266 to 1,103 — removing a single sample can quadruple the feature set. This instability means the "clock" is a proof-of-concept, not a biomarker.

**Ages are decade-binned** (5 unique midpoints: 34.5, 44.5, 54.5, 64.5, 74.5). This imposes a ~2.5-year noise floor on MAE and inflates Pearson r by collapsing the target to 5 groups. The corrected MAE of 7.6 years means the model is off by roughly 1.5 age bins on average.

**No covariate adjustment.** Mixed sex, one severe autolysis sample (GTEX-X62O), no Sentrix batch correction. In a 47-sample cohort, any of these can confound the age signal.

**The |r| > 0.5 threshold is arbitrary.** A threshold of 0.4 or 0.6 produces a substantially different clock. Nested CV protects against inflated metrics but doesn't address the researcher degree of freedom in choosing the filter.

## What a Production Clock Needs

To go from this prototype to something publishable:

1. **N > 200** with continuous age measurements (not decade bins)
2. **Batch correction** (ComBat or equivalent) for Sentrix chip and position effects
3. **Sex as a covariate** or sex-stratified models
4. **External validation** on an independent cohort (MESA muscle, or additional GTEx batches)
5. **Head-to-head comparison** with pan-tissue clocks (Horvath, Hannum, PhenoAge, GrimAge) applied to the same samples

The value of this notebook isn't the clock itself — it's the pipeline, and the demonstration that feature selection leakage can inflate r by 0.40 in a real dataset. That's a number worth knowing before your next interview.
