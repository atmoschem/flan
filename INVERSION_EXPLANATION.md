# FLAN: Analytical Inversion Theory & Operational Manual

> **Disclosure**: While the implementation logic leverages standard atmospheric inverse modeling theory, parts of this theoretical explanation and code structure were developed with the assistance of Google Antigravity.

This document details the **4D Analytical Inversion** system implemented in FLAN. It covers the mathematical foundation, the memory-efficient implementation (Implicit Kronecker), and the operational strategy for real-time estimation.

---

## 1. Theoretical Foundation

The model solves for the "true" state of emissions by optimizing the match between observed atmospheric concentrations and modeled transport, subject to prior constraints.

### The Cost Function
We minimize the Bayesian cost function $J(x)$:

$$J(x) = (x - x_a)^T B^{-1} (x - x_a) + (y - Hx)^T R^{-1} (y - Hx)$$

Where:
*   $x$: State vector (Scaling Factors for emissions).
*   $x_a$: **Prior** state vector (Initial guess, typically 1.0).
*   $B$: **Prior Covariance** matrix (Uncertainty in $x_a$).
*   $y$: **Observations** vector (Mixing ratios in ppm/ppb).
*   $H$: **Sensitivity Matrix** (Jacobian/Footprints).
*   $R$: **Observation Covariance** matrix (Measurement + Model error).

### The Analytical Solution
Since the transport is linear, we solve for the posterior state $x_{post}$ directly using the **Kalman Gain ($K$)**:

$$ K = B H^T (H B H^T + R)^{-1} $$
$$ x_{post} = x_a + K (y - H x_a) $$

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

## 3. Efficient Covariance (Implicit Kronecker)

A major challenge in 4D inversion is the size of the $B$ matrix. For a state vector of size 120,000, $B$ requires ~115 GB of RAM.

**Solution**: We define $B$ as the **Kronecker Product** of temporal and spatial correlations:
$$ B = B_{temporal} \otimes B_{spatial} $$

We avoid constructing $B$ explicitly. Instead, we compute the term $Z = B H^T$ used in the Kalman Gain by exploiting the property:
$$ (B_t \otimes B_s) \cdot \text{vec}(V) = \text{vec}(B_s \cdot V \cdot B_t^T) $$

This reduces memory usage from **100+ GB** to **~50 MB**, allowing high-resolution 4D inversions on standard hardware.

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

---

## 5. Operational Inputs & Outputs

### Input Requirements
*   **Observations (CSV)**: For a real-time run, provide only the **new** observations (e.g., the latest 3 hours).
*   **Priors & Footprints**: Must cover the full **Look-Back Period** (e.g., 10 days) to allow the convolution to "see" the history of the air parcels.

### Expected Dimensions
*   **Input Prior**: `[Lon, Lat, N_time]` (e.g., 240 distinct 3-hour blocks).
*   **Output Posterior**: `[Lon, Lat, N_time]`.
    *   `scaling_factors_3d.nc`: The optimized scaling factors.
    *   `posterior_flux.nc`: The final flux values.

### Applicability
While designed for resolving the **Diurnal Cycle of Biogenic CO2** (Photosynthesis vs. Respiration), this 3-hourly 4D framework is valid for:
*   **Methane (CH4)**: Resolving intermittent leaks or wetland variability.
*   **Air Quality (NOx/CO)**: Tracking sub-daily urban pollution plumes.