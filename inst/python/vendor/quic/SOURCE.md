# QUIC provenance

QUIC 1.2 by Cho-Jui Hsieh and Matyas A. Sustik, from skggm 0.2.8 source distribution.
Upstream: https://github.com/skggm/skggm/tree/master/inverse_covariance/pyquic
Download: https://pypi.org/project/skggm/0.2.8/
The source header specifies GNU GPL version 3. License text: https://www.gnu.org/licenses/gpl-3.0.txt

Local adaptation: the LAPACK dimension and info declarations use `int` instead of
`ptrdiff_t`, matching the default LP64 macOS Accelerate ABI. No optimization
algorithm, objective, or tuning parameter is altered. The Python ctypes binding
is in ../../quic_solver.py. For redistribution, retain upstream authorship and
the accompanying GPL license. The original source distribution is available
through the download link above; it is not bundled separately here.

The Python binding compiles this source on demand into the caller's output
directory using Accelerate on macOS or LP64 LAPACK/BLAS on Linux. See
`inst/reproduction/EMPIRICAL.md` in the source repository for requirements.
No compiled library is distributed in the R package.
