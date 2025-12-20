# Fortran LagrangiAN Inverse Model (FLAN)
My cool new project!


<img src="figs/logo2.png" align="right" alt="" width="220" />

Cool thing, it is an [fpm](https://fpm.fortran-lang.org/index.html) project that already compiles with

- LAPACK
- BLAS
- STDLIB
- NETCDF
- OPENMP

## package dependencies

see [fpm.toml]()

- openmp
- hdf5
- netcdf
- blas
- stdlib =  https://github.com/fortran-lang/stdlib
- datetime = https://github.com/wavebitscientific/datetime-fortran
- fpm = https://github.com/fortran-lang/fpm
- csv-fortran = https://github.com/jacobwilliams/csv-fortran.git

fpm will find libraries in teh system such as openblas, netcdf, df5 and openmp.
Also, it will download and compile 

- stdlib
- datetime to work with time 
- fpm extra functions including strings
- csv-fortran, to read csv files

## Installation

fpm needs access to libraries above mentioned. For instance, for an HPC, you could:

```bash
module load slurm
module load netcdf-c
module load netcdf-fortran
module load hdf5/1.14.3
module load openblas
```
then

## build

fpm build

```
fpm build
fpm_backend_console.f90                done.
filesystem_utilities.c                 done.
fpm_strings.f90                        done.
iscygpty.c                             done.
isatty.c                               done.
csv_kinds.f90                          done.
constants.f90                          done.
version.f90                            done.
token.f90                              done.
M_CLI2.F90                             done.
regex.f90                              done.
version.f90                            done.
shlex_module.f90                       done.
stdlib_ascii.f90                       done.
stdlib_linalg_constants.F90            done.
stdlib_optval.f90                      done.
error.f90                              done.
csv_parameters.f90                     done.
error.f90                              done.
datetime.f90                           done.
io.f90                                 done.
stdlib_error.f90                       done.
stdlib_linalg_blas_aux.f90             done.
stdlib_string_type.f90                 done.
stdlib_blas_constants.f90              done.
fpm_environment.f90                    done.
versioning.f90                         done.
csv_utilities.f90                      done.
utils.f90                              done.
abc.f90                                done.
f08estop.f90                           done.
stdlib_strings.f90                     done.
stdlib_blas.f90                        done.
stdlib_io.f90                          done.
fpm_filesystem.F90                     done.
fpm_release.F90                        done.
csv_module.F90                         done.
terminal.f90                           done.
value.f90                              done.
lexer.f90                              done.
io_manager.f90                         done.
stdlib_codata_type.f90                 done.
stdlib_linalg_state.f90                done.
stdlib_blas_level2_sym.f90             done.
stdlib_strings_to_string.f90           done.
stdlib_blas_level2_tri.f90             done.
stdlib_blas_level2_ban.f90             done.
stdlib_blas_level3_tri.f90             done.
stdlib_blas_level3_gen.f90             done.
stdlib_blas_level3_sym.f90             done.
stdlib_blas_level2_pac.f90             done.
stdlib_string_type_constructor.f90     done.
stdlib_linalg_blas.F90                 done.
stdlib_blas_level2_gen.f90             done.
stdlib_blas_level1.f90                 done.
fpm_os.F90                             done.
fpm_pkg_config.f90                     done.
diagnostic.f90                         done.
sort.f90                               done.
keyval.f90                             done.
node.f90                               done.
map.f90                                done.
list.f90                               done.
stdlib_codata.f90                      done.
stdlib_linalg.f90                      done.
stdlib_linalg_lapack_aux.f90           done.
fpm_command_line.f90                   done.
context.f90                            done.
array_list.f90                         done.
ordered_map.f90                        done.
linear_algebra.f90                     done.
stdlib_constants.f90                   done.
stdlib_lapack_base.f90                 done.
stdlib_linalg_cross_product.f90        done.
stdlib_linalg_diag.f90                 done.
stdlib_linalg_kronecker.f90            done.
stdlib_linalg_outer_product.f90        done.
structure.f90                          done.
lexer.f90                              done.
main.f90                               done.
stdlib_lapack_blas_like_scalar.f90     done.
stdlib_lapack_auxiliary.f90            done.
stdlib_lapack_blas_like_l1.f90         done.
stdlib_lapack_blas_like_base.f90       done.
stdlib_lapack_blas_like_mnorm.f90      done.
stdlib_lapack_solve.f90                done.
stdlib_lapack_givens_jacobi_rot.f90    done.
stdlib_lapack_blas_like_l2.f90         done.
stdlib_lapack_blas_like_l3.f90         done.
stdlib_lapack_householder_reflectors.f9done.
stdlib_lapack_orthogonal_factors.f90   done.
table.f90                              done.
array.f90                              done.
stdlib_lapack_solve_chol_comp.f90      done.
stdlib_lapack_solve_tri_comp.f90       done.
stdlib_lapack_solve_ldl_comp4.f90      done.
stdlib_lapack_solve_ldl.f90            done.
stdlib_lapack_eig_svd_lsq.f90          done.
stdlib_lapack_solve_ldl_comp3.f90      done.
stdlib_lapack_solve_ldl_comp.f90       done.
stdlib_lapack_solve_lu.f90             done.
stdlib_lapack_orthogonal_factors_rz.f90done.
stdlib_lapack_solve_lu_comp.f90        done.
stdlib_lapack_orthogonal_factors_ql.f90done.
stdlib_lapack_orthogonal_factors_qr.f90done.
stdlib_lapack_solve_chol.f90           done.
stdlib_lapack_solve_aux.f90            done.
stdlib_lapack_others.f90               done.
stdlib_lapack_solve_ldl_comp2.f90      done.
type.f90                               done.
stdlib_lapack_eigv_svd_drivers.f90     done.
stdlib_lapack_cosine_sine.f90          done.
stdlib_lapack_eigv_std_driver.f90      done.
stdlib_lapack_eigv_svd_bidiag_dc.f90   done.
stdlib_lapack_eigv_gen2.f90            done.
stdlib_lapack_eigv_tridiag3.f90        done.
stdlib_linalg_lapack.F90               done.
stdlib_lapack_svd_comp.f90             done.
stdlib_lapack_svd_bidiag_qr.f90        done.
stdlib_lapack_lsq_aux.f90              done.
stdlib_lapack_eigv_gen3.f90            done.
stdlib_lapack_eigv_sym.f90             done.
stdlib_lapack_cosine_sine2.f90         done.
stdlib_lapack_lsq.f90                  done.
stdlib_lapack_eigv_comp2.f90           done.
stdlib_lapack_lsq_constrained.f90      done.
stdlib_lapack_others_sm.f90            done.
stdlib_lapack_eigv_svd_drivers2.f90    done.
stdlib_lapack_eigv_svd_drivers3.f90    done.
stdlib_lapack_eigv_gen.f90             done.
stdlib_lapack_eigv_sym_comp.f90        done.
stdlib_lapack_eigv_comp.f90            done.
stdlib_lapack_eigv_tridiag2.f90        done.
stdlib_lapack_svd_comp2.f90            done.
stdlib_lapack_eigv_tridiag.f90         done.
ser.f90                                done.
merge.f90                              done.
keyval.f90                             done.
parser.f90                             done.
ser.f90                                done.
stdlib_linalg_schur.f90                done.
stdlib_linalg_matrix_functions.f90     done.
stdlib_linalg_qr.f90                   done.
stdlib_linalg_inverse.f90              done.
stdlib_linalg_least_squares.f90        done.
stdlib_linalg_determinant.f90          done.
stdlib_linalg_eigenvalues.f90          done.
stdlib_linalg_svd.f90                  done.
stdlib_linalg_solve.f90                done.
stdlib_linalg_norms.f90                done.
stdlib_linalg_pinv.f90                 done.
stdlib_linalg_cholesky.f90             done.
de.f90                                 done.
table.f90                              done.
array.f90                              done.
path.f90                               done.
build.f90                              done.
tomlf.f90                              done.
parser.f90                             done.
jonquil.f90                            done.
downloader.f90                         done.
toml.f90                               done.
fpm_settings.f90                       done.
fpm_compile_commands.F90               done.
git.f90                                done.
fortran.f90                            done.
library.f90                            done.
profiles.f90                           done.
build.f90                              done.
preprocess.f90                         done.
meta.f90                               done.
install.f90                            done.
fpm_compiler.F90                       done.
dependency.f90                         done.
platform.f90                           done.
executable.f90                         done.
test.f90                               done.
example.f90                            done.
feature.f90                            done.
feature_collection.f90                 done.
package.f90                            done.
manifest.f90                           done.
dependency.f90                         done.
fpm_model.f90                          done.
fpm_source_parsing.f90                 done.
fpm_meta_base.f90                      done.
fpm_sources.f90                        done.
fpm_meta_openmp.f90                    done.
fpm_meta_util.f90                      done.
fpm_meta_stdlib.f90                    done.
fpm_meta_minpack.f90                   done.
fpm_targets.f90                        done.
fpm_meta_mpi.f90                       done.
fpm_meta_netcdf.f90                    done.
fpm_meta_blas.f90                      done.
fpm_meta_hdf5.f90                      done.
fpm_meta.f90                           done.
fpm_backend_output.f90                 done.
fpm_backend.F90                        done.
fpm.f90                                done.
publish.f90                            done.
libflan.a                              done.
flan                                   done.
[100%] Project compiled successfully.
```

## run

```
fpm run
 found blas package: openblas
Project is up to date
 ---------------------------------------------------
 Reading configuration from 'namelists/config.nml'...
 input_config: receptor: receptors/receptor_1h.csv
 input_config: prior: /media/sergio/ext6/noaa/co2_nmolm2s_denver_240.nc
 model_config: add_bg: false
 ---------------------------------------------------
 Reading CSV file: receptors/receptor_1h.csv
 --- CSV Summary ---
 Total records:            7
 First Record:
   Time:         2024           9          17           7           0           0
   Gas:    10.100000000000000     
   Bg:    2.0200000000000000     
   Path: /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.7861Nx104.9886Wx00002/hysplit2024x09x17x07x00x39.7861Nx104.9886Wx00002.nc                                                                                                                                    
 Last Record:
   Time:         2024           9          17           7           0           0
   Gas:    3.2000000000000002     
   Bg:   0.64000000000000001     
   Path: /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.4381Nx108.0261Wx00002/hysplit2024x09x17x07x00x39.4381Nx108.0261Wx00002.nc                                                                                                                                    
 First receptor time object: 2024-09-17T07:00:00.000
 Total receptors loaded:            7
 ---------------------------------------------------
 Reading prior
 Prior dimensions (lon, lat, time):           25          20         240
 Prior range (min, max):    10.813813209533691        30470.574218750000     
 Prior sum:    152436660.70426941     
 ---------------------------------------------------
 Reading footprints
 State vector size (N_grid):          500
 Processing footprint            1 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.7861Nx104.9886Wx00002/hysplit2024x09x17x07x00x39.7861Nx104.9886Wx00002.nc
    Footprint            1  range (min, max):    0.0000000000000000        6.5165825188159943E-002
 Processing footprint            2 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x37.3039Nx107.4842Wx00002/hysplit2024x09x17x07x00x37.3039Nx107.4842Wx00002.nc
    Footprint            2  range (min, max):    0.0000000000000000        0.0000000000000000     
 Processing footprint            3 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.9128Nx105.1886Wx00002/hysplit2024x09x17x07x00x39.9128Nx105.1886Wx00002.nc
    Footprint            3  range (min, max):    0.0000000000000000        4.0791746228933334E-002
 Processing footprint            4 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.7795Nx105.0052Wx00002/hysplit2024x09x17x07x00x39.7795Nx105.0052Wx00002.nc
    Footprint            4  range (min, max):    0.0000000000000000        5.4270178079605103E-002
 Processing footprint            5 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.8381Nx104.9498Wx00002/hysplit2024x09x17x07x00x39.8381Nx104.9498Wx00002.nc
    Footprint            5  range (min, max):    0.0000000000000000        4.9182698130607605E-002
 Processing footprint            6 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.7322Nx105.0153Wx00002/hysplit2024x09x17x07x00x39.7322Nx105.0153Wx00002.nc
    Footprint            6  range (min, max):    0.0000000000000000        4.7325719147920609E-002
 Processing footprint            7 : /media/sergio/ext6/2024/09/tmp_2024x09x17x07x00x39.4381Nx108.0261Wx00002/hysplit2024x09x17x07x00x39.4381Nx108.0261Wx00002.nc
    Footprint            7  range (min, max):    0.0000000000000000        0.0000000000000000     
 H matrix dimensions (n_obs, n_grid):            7         500
 ---------------------------------------------------
 Enhancements co2 (ppb)
   1 :      3468.28477
   2 :         0.00000
   3 :        47.05454
   4 :      2402.17249
   5 :      3661.73907
   6 :      1386.99579
   7 :         0.00000
 ---------------------------------------------------
 Skipping background.
 Writing hsp to: hsp/hsp_output.txt
 ---------------------------------------------------
 Creating R matrix (measurement and model error covariance)
 R matrix dimensions (n_obs, n_obs):            7           7
 R diagonal (first 5 elements):    2.0000000000000000        2.0000000000000000        2.0000000000000000        2.0000000000000000        2.0000000000000000     
 ---------------------------------------------------
 Performing Kalman Inversion...
 B matrix dimensions (n_grid, n_grid):          500         500
 B diagonal (variance):   0.25000000000000000     
 S matrix dimensions (n_obs, n_obs):            7           7
 Kalman Gain K dimensions (n_grid, n_obs):          500           7
 Posterior state x_post size:          500
 Posterior scaling factors range:   -1.3543819814210214       to    1.0635400089597165     
 Average scaling factor:   0.94587940100104351     
 ---------------------------------------------------
 Saving posterior outputs...
 Reshaped scaling factors (sf_map) dimensions:           25          20
 Writing scaling factors to: nc/scaling_factors.nc
 Posterior flux dimensions:           25          20         240
 Average posterior flux:    954.40075301462684     
 Writing posterior flux to: nc/posterior_flux.nc
 ---------------------------------------------------
 ("--- Datetime Test ---")
 2025-12-20T01:00:52.803```

