# fastDiversity 1.3.0

## Changes to existing functions
- `calculate_pi_theta_d()`: reimplemented with a vectorised approach for 
  improved performance. Non-segregating and monomorphic sites are now removed 
  prior to calculation. The function now accepts the same `gt` matrix format 
  used throughout the package. Return value structure has changed: the function 
  now returns `tajima_d`, `num_sites`, `raw_pi`, `watterson_theta`, and 
  `d_stdev` (previously returned `pi`, `S`, `theta`, `D`, and `per_site_stats`).