# Analytical Inversion Theory & Implementation

This document provides a step-by-step explanation of the analytical inversion logic implemented in the FLAN (Fortran Lagrangian ANalysis) model. The implementation uses a Bayesian framework, specifically the **Kalman Inversion** method, to derive posterior fluxes from prior information and observations.

---

## 1. Theoretical Background

The goal of analytical inversion is to estimate "true" state variables (e.g., emissions/fluxes) by combining **prior knowledge** with **observations** (e.g., atmospheric concentrations).

The cost function minimized in this Bayesian framework is:
$$J(x) = (x - x_a)^T B^{-1} (x - x_a) + (y - Hx)^T R^{-1} (y - Hx)$$

Where:
- $x$: State vector (what we want to find). In this code, these are **scaling factors** for fluxes.
- $x_a$: **Prior state vector**.
- $B$: **Prior error covariance matrix** (uncertainty in $x_a$).
- $y$: **Observations** vector (receptor concentrations).
- $H$: **Observation operator** (Sensitivity/Jacobian matrix).
- $R$: **Observation error covariance matrix** (measurement + model error).

The analytical solution for the posterior state $x_{post}$ is:
$$x_{post} = x_a + G(y - Hx_a)$$
Where $G$ is the **Kalman Gain**:
$$G = B H^T (H B H^T + R)^{-1}$$

---

## 2. Step-by-Step Code Walkthrough

The logic is primarily located in `app/main.f90`.

### Step 1: Building the Sensitivity Matrix ($H$)
The matrix $H$ represents how a unit change in emissions at each grid cell affects the concentration at each receptor.

```fortran
! Lines 126-153 in main.f90
n_grid = size(prior_in, 1) * size(prior_in, 2)
allocate(h_mat(n_obs, n_grid))

do i = 1, size(receptor_path)
   ! Read footprint for receptor i
   call read_3d_netcdf(trim(receptor_path(i)), foot_in, foot_name)
   ! H matrix is the footprint multiplied by the prior, integrated over time
   h_mat(i, :) = reshape(sum(foot_in * prior_in, dim=3), [n_grid])
end do
```
- **Theory**: $H_{ij} = \frac{\partial y_i}{\partial x_j}$. 
- **Code**: Since $x$ represents scaling factors (initial value 1.0), the sensitivity is the total enhancement contributed by grid cell $j$ to receptor $i$ given the prior emissions.

### Step 2: Observation Error Covariance ($R$) (NEEDS TO BE CHANGED)
$R$ accounts for how much we trust our observations and our transport model.

```fortran
! Lines 178-181 in main.f90
allocate(r_diag(size(hsp)))
r_diag = obs_err_sd**2 + model_err_sd**2
r_mat = diag(r_diag)
```
- The code assumes errors are uncorrelated (diagonal matrix).
- The total variance is the sum of measurement error (`obs_err_sd`) and model transport error (`model_err_sd`).

### Step 3: Prior State and Covariance ($x_a$ and $B$)
We start with a "guess" that our prior emissions are correct (scaling factors = 1.0).

```fortran
! Lines 190-198 in main.f90
xa = 1.0_wp ! Prior scaling factors
allocate(b_mat(n_grid, n_grid))
do i = 1, n_grid
  b_mat(i, i) = prior_err_sd**2
end do
```
- $x_a$ is initialized to 1.0.
- $B$ is a diagonal matrix where each entry is the variance of the scaling factors (`prior_err_sd**2`).

### Step 4: Solving the Inversion (Kalman Gain)
This is where the matrix math happens.

```fortran
! Lines 205-210 in main.f90
! S = H B H' + R
s_mat = matmul(h_mat, matmul(b_mat, transpose(h_mat))) + r_mat

! Invert S
call invert_matrix(s_mat)

! gain_k = B H' inv(S)
gain_k = matmul(b_mat, matmul(transpose(h_mat), s_mat))
```
- **$S$ matrix**: This is the innovation covariance. Its dimension is $N_{obs} \times N_{obs}$.
- `invert_matrix`: Uses `stdlib_linalg` to perform the inversion.
- **$K$ (gain_k)**: The key weighting factor.

### Step 5: Updating to Posterior State ($x_{post}$)
We calculate the final scaling factors based on the "mismatch" between observed $(y)$ and modeled $(H x_a)$ concentrations.

```fortran
! Line 216 in main.f90
x_post = xa + matmul(gain_k, (receptor_gas - hsp))
```
- `receptor_gas`: The actual concentration observed at the receptor.
- `hsp`: The modeled concentration ($H \cdot x_a$).
- The difference (innovation) is multiplied by the Kalman Gain to update the scaling factors.

### Step 6: Generating Posterior Fluxes
Finally, the scaling factors are applied to the 3D prior flux maps.

```fortran
! Lines 238-241 in main.f90
do i = 1, n_time
  post_flux(:, :, i) = prior_in(:, :, i) * sf_map
end do
```
- Each time-step of the prior emission field is scaled by the calculated `sf_map` to get the final "optimized" fluxes.

---

## 3. Summary of Variables

| Code Variable | Math Symbol | Description |
| :--- | :---: | :--- |
| `h_mat` | $H$ | Sensitivity Matrix (Jacobian) |
| `xa` | $x_a$ | Prior State Vector (Scaling Factors) |
| `b_mat` | $B$ | Prior Error Covariance |
| `r_mat` | $R$ | Observation Error Covariance |
| `s_mat` | $S$ | Innovation Covariance ($HBH^T + R$) |
| `gain_k` | $K$ | Kalman Gain |
| `x_post` | $x_{post}$ | Posterior State Vector |
| `receptor_gas` | $y$ | Observed Concentrations |

---

## 4. Supporting Utilities
The inversion relies on the `linear_algebra` module in `src/linear_algebra.f90`, which provides:
- **`invert_matrix(A)`**: A wrapper for the high-level `inv()` function from the Fortran Standard Library (`stdlib`).
- **`diag(v)`**: Creates a diagonal matrix from a vector.

##TODO

- Improve dimenisons handling between H and x. H is hourly back in days, x can be annual or monthly, and so on.
- Add techniques such as exponential decays for the spatial autocorelation.