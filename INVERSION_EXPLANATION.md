# FLAN: Analytical Inversion Theory & Operational Manual

> **Disclosure**: While the implementation logic leverages standard atmospheric inverse modeling theory, parts of this theoretical explanation and code structure were developed with the assistance of Google Antigravity.

This document details the **4D Analytical Inversion** system implemented in FLAN. It covers the mathematical foundation, the memory-efficient implementation (Implicit Kronecker), and the operational strategy for real-time estimation.

---

## 1. Theoretical Foundation

The model solves for the "true" state of emissions by optimizing the match between observed atmospheric concentrations and modeled transport, subject to prior constraints.

### The Cost Function
We minimize the Bayesian cost function $J(x)$, which represents the trade-off between sticking to our initial guess and fitting the new measurements:

$$J(x) = (x - x_a)^T B^{-1} (x - x_a) + (y - [Hx + y_{bg}])^T R^{-1} (y - [Hx + y_{bg}])$$

Where:
*   $x$: State vector (**Scaling Factors** for emissions).
*   $x_a$: **Prior** state vector (typically 1.0).
*   $B$: **Prior Covariance** (Errors in emission maps, handles spatial/temporal smoothing).
*   $y$: **Observations** (Mixing ratios, usually in ppb).
*   $y_{bg}$: **Background** CO2 (The air concentration before it entered our study domain).
*   $H$: **Sensitivity Matrix** (Jacobian/Footprints conjoined with the Prior Flux).
*   $R$: **Observation Covariance** (Measurement Noise + Model/Transport Error).

### The Analytical Solution
Since transport is linear, we solve for the posterior state $x_{post}$ exactly using the **Kalman Gain ($K$)**:

$$ K = B H^T (H B H^T + R)^{-1} $$
$$ x_{post} = x_a + K (y - [H x_a + y_{bg}]) $$

---

## 2. The 4D State Vector & Convolution

Unlike static inversions that optimize a single map, this system is **Time-Resolved (4D)**.

*   **State Vector Dimensions**: $N_{lon} \times N_{lat} \times N_{time}$ (e.g., $25 \times 20 \times 80$).
*   **Resolution**: Optimized for **3-hourly** fluxes.

### The Lagged Convolution (H Matrix)
An atmospheric observation at time $t$ depends on surface fluxes from time $t$ back to $t - 10$ days. The $H$ matrix represents this Lagrangian convolution:

$$ y_{obs}(t) = \sum_{\tau=0}^{-240h} \text{Footprint}(\tau) \times \text{Flux}(t + \tau) $$

In the code (`main.f90`), this is implemented by element-wise multiplication of the 3D footprint with the 3D prior field, flattened into the row of $H$ corresponding to that observation.

---

This reduces memory usage from **100+ GB** to **~50 MB**, allowing high-resolution 4D inversions on standard hardware.

---

## 4. Error Covariance Treatment (R & B)

### Prior Support ($B$ Matrix)
FLAN supports **Multi-Component Priors**. You can provide separate files for Anthropogenic and Natural (Biogenic) emissions. 
*   **Summation**: The code sums them ($Prior_{tot} = Prior_{anth} + Prior_{nat}$) before building $H$.
*   **Correlation**: The Kronecker smoothing ($B_t \otimes B_s$) ensures that information from an observation at one tower is spread to neighboring pixels and time-steps in a physically consistent way.

### Observation Error ($R$ Matrix)
To handle **Model-Data Mismatch (MDM)**, $R$ is treated as a diagonal matrix where each entry $R_{ii}$ accounts for:
1.  **Measurement Error**: Instrument precision (read from the 9th column of the CSV).
2.  **Model/Transport Error**: Uncertainty in wind fields and representation error (set in `config.nml`).

$$R_{ii} = \sigma_{meas}^2 + \sigma_{transport}^2$$

---

## 4. Operational Real-Time Strategy

The system is designed to run in an operational **Sliding Window** (Fixed-Lag Kalman Smoother) mode.

### The Problem: Aliasing
A single observation at 15:00 provides information about emissions at 15:00 *and* emissions from yesterday. If we strictly inverted for 12:00-15:00 using only current data, we would misattribute yesterday's signals to today's fluxes.

### The Solution: Sequential Cycle
We define a look-back window (e.g., 10 days) matching the footprint length.

#### Cycle A (Current Time $T$)
1.  **Assimilate**: Observations from $T-3h$ to $T$.
2.  **Optimize**: The full 10-day history.
3.  **Result**: 
    *   $T-10d$: Finalized.
    *   $T-24h$: Refined (Seen by 24h of observations).
    *   $T$: Preliminary (Seen only by current observations).

#### Cycle B (Time $T+3h$)
1.  **Shift**: Slide the time window forward by 3 hours.
2.  **Prior Update**: Use the **Posterior from Cycle A** as the **Prior for Cycle B** for the overlapping period.
3.  **Assimilate**: New observations from $T$ to $T+3h$.

This effectively propagates information forward, ensuring that today's preliminary estimate becomes tomorrow's refined estimate.

### Implementation in FLAN
You can trigger this behavior using the `use_prior_sf` toggle in `config.nml`. 
- When enabled, the model reads a 3D scaling factor file from a previous run.
- It performs a **1-step time shift**: the information at Time Index $k$ in Run A becomes the Prior at Time Index $k-1$ in Run B. 
- The very last time step in the new window (the new data) is initialized with the standard prior of **1.0**.

---

## 5. Operational Inputs & Outputs

### Input Requirements
*   **Observations (CSV)**: For a real-time run, provide only the **new** observations (e.g., the latest 3 hours).
*   **Priors & Footprints**: Must cover the full **Look-Back Period** (e.g., 10 days) to allow the convolution to "see" the history of the air parcels.

### Expected Dimensions
*   **Input Prior**: `[Lon, Lat, N_time]` (e.g., 240 distinct 3-hour blocks).
*   **Output Posterior**:
    *   `nc/scaling_factors.nc`: Time-averaged optimized multipliers.
    *   `nc/posterior_flux.nc`: The final physical flux values ($x_{post} \times Prior_{tot}$).
    *   `nc/posterior_uncertainty.nc`: The analytical standard deviation for every pixel and time step.

### Applicability
While designed for resolving the **Diurnal Cycle of Biogenic CO2** (Photosynthesis vs. Respiration), this 3-hourly 4D framework is valid for:
*   **Methane (CH4)**: Resolving intermittent leaks or wetland variability.
*   **Air Quality (NOx/CO)**: Tracking sub-daily urban pollution plumes.
---

## 6. Statistical Validation: The Chi-Square Diagnostic

To ensure the inversion is reliable, FLAN calculates the **Reduced Chi-Square ($\chi_r^2$)**:

$$\chi_r^2 = \frac{J(x_{post})}{N_{obs}}$$

*   **$\chi_r^2 \approx 1.0$**: The "Gold Standard." Your error assumptions ($R$ and $B$) are perfectly balanced with your model residuals.
*   **$\chi_r^2 \gg 1.0$**: Underestimated errors. The model is struggling to fit the data. You may be missing sources (e.g., Natural emissions) or your transport error is too small.
*   **$\chi_r^2 \ll 1.0$**: Overestimated errors. The model is "over-fitting" the noise. Your error bars are too large.

### Unit Consistency
The model expects all concentration units to be in **ppb** (parts per billion). If your CSV is in **ppm**, you must multiply by 1000 before inputting, otherwise the Chi-Square will explode due to the 1000x scale mismatch.
