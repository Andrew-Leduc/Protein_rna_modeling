# Modeling Approaches

## Central Question

Can the observed distribution of protein abundance across single cells be explained by the underlying mRNA count distribution under a model of constant translation rate and constant protein degradation? Or does the protein variance require additional sources of noise beyond what mRNA variability provides?

We test this using both **log2 fold change** and **absolute copy number** representations of protein abundance.

---

## 1. mRNA Distribution Modeling

### 1a. Empirical Fit (Negative Binomial)

Fit the observed mRNA count distribution to a negative binomial:

$$P(X = x) = \frac{\Gamma(x + r)}{\Gamma(r)\, x!} \left(\frac{r}{r + \mu}\right)^r \left(\frac{\mu}{r + \mu}\right)^x$$

where $\mu$ is the mean and $r$ is the dispersion (size) parameter. This captures overdispersion relative to Poisson but does not directly inform about the underlying transcriptional mechanism.

### 1b. Mechanistic mRNA Kinetics (Spliced/Unspliced)

Using spliced and unspliced counts, we can learn the kinetic parameters of transcription following the approach of the [Pachter lab](https://www.nature.com/articles/s41592-025-02832-x).

The telegraph model of transcription assumes the gene switches between inactive and active states:

$$\text{OFF} \xrightarrow{k_{on}} \text{ON} \xrightarrow{k_{off}} \text{OFF}$$

In the ON state, mRNA is transcribed at rate $s$. The dynamics of unspliced ($u$) and spliced ($m$) mRNA are:

$$\frac{du}{dt} = s \cdot z(t) - \beta \, u$$

$$\frac{dm}{dt} = \beta \, u - \gamma_m \, m$$

where $z(t) \in \{0, 1\}$ is the gene state, $\beta$ is the splicing rate, and $\gamma_m$ is the mRNA degradation rate.

From the steady-state distributions and the ratio of unspliced to spliced counts, we can estimate:
- **Burst frequency**: $k_{on} / \gamma_m$
- **Burst size**: $s / k_{off}$
- **mRNA half-life**: $t_{1/2}^m = \ln(2) / \gamma_m$

---

## 2. Protein Abundance Models

Given the mRNA distribution (from either 1a or 1b), we model protein levels under constant translation and degradation.

### 2a. Steady-State Copy Number Model

At steady state, the expected protein copy number for a cell with $m$ mRNA molecules is:

$$\langle P \rangle = \frac{k_t}{\gamma_p} \cdot m$$

where $k_t$ is the translation rate (proteins/mRNA/unit time) and $\gamma_p$ is the protein degradation rate.

The full stochastic model gives:

$$\frac{dP}{dt} = k_t \cdot m(t) - \gamma_p \cdot P(t)$$

Because mRNA fluctuates, protein inherits that variability filtered through the protein lifetime $t_{1/2}^p = \ln(2) / \gamma_p$.

### 2b. Simplified Per-Molecule Model

A simpler formulation: each mRNA molecule produces a fixed number of proteins $c$ over its lifetime:

$$c = \frac{k_t}{\gamma_m}$$

giving total protein as:

$$P = c \cdot m + \epsilon$$

This is a snapshot model that ignores temporal averaging but is useful as a null hypothesis.

### 2c. Noise Decomposition

The observed protein variance can be decomposed as:

$$\text{Var}(P) = \underbrace{\left(\frac{k_t}{\gamma_p}\right)^2 \text{Var}(m)}_{\text{mRNA-driven}} + \underbrace{\sigma^2_{\text{extrinsic}}}_{\text{other sources}}$$

If $\sigma^2_{\text{extrinsic}} \approx 0$, then mRNA variability alone explains protein variability under constant translation. If not, additional noise sources (translational bursting, cell-to-cell variation in rates, etc.) are needed.

---

## 3. Fitting Strategy

### Parameters from data
| Parameter | Source |
|-----------|--------|
| $\mu_m, r_m$ (mRNA NB params) | Fit to mRNA count distribution |
| $k_{on}, k_{off}, s, \gamma_m$ | Spliced/unspliced kinetics |
| $\gamma_p$ | Measured protein degradation rate |
| $k_t$ or $\sigma$ | **Fit to match observed protein distribution** |

### Approach

1. Fit the mRNA distribution (NB empirically, or mechanistic via spliced/unspliced)
2. Forward-simulate protein abundance under constant translation + degradation
3. Compare simulated protein distribution to observed (Wasserstein distance or KS test)
4. Optimize $k_t$ (and optionally $\sigma_{\text{extrinsic}}$) to minimize distance
5. Evaluate whether $\sigma_{\text{extrinsic}} \approx 0$ — i.e., is constant translation sufficient?

### Both representations
- **Absolute copies**: Compare simulated $P$ directly to observed copy numbers
- **Log2 fold change**: Compare $\log_2(P) - \overline{\log_2(P)}$ to observed log2FC

---

## References

- Pachter lab mRNA kinetics from spliced/unspliced counts: [Bravo González-Blas et al., Nature Methods (2025)](https://www.nature.com/articles/s41592-025-02832-x)
