double ReconfigurableIntelligentSurface::gain(double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad) {
  // fix the corner case in which phi is between -180 and -180 + d_phi/2
  // the gain matrix (if d_phi = 1 degree) would have values from -179 to 180
  // if phi = -179.8, for example, this will fall into the 180 degrees bin, as
  // it is equivalent to 180.2 degrees the 180 degrees bin indeed goes from
  // 179.5 to 180.5 degrees (for d_phi = 1 degree)
  phiRX_rad = fix_azimuth(phiRX_rad);
  phiTX_rad = fix_azimuth(phiTX_rad);

  if (canUseCache(phiTX_rad, thetaTX_rad)) {
    return cachedGain(phiRX_rad, thetaRX_rad);
  }

  CMatrix phase = new_cmatrix(Ts, Ps);

  Matrix PHI = coding;
  double alpha;
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      // compute the phase offset due to the incidence of the signal
      alpha = du_k * (n * sin(thetaTX_rad) * cos(phiTX_rad) +
                      m * sin(thetaTX_rad) * sin(phiTX_rad));
      // compute the phase offsets for all the possible phiRX,thetaRX pairs
      Matrix n_k_du_sin_cos = new_matrix(k_du_sin_cos);
      Matrix m_k_du_sin_sin = new_matrix(k_du_sin_sin);
      // scale by negative m and n, as these matrices need to be subtracted
      gsl_matrix_scale(n_k_du_sin_cos, -n);
      gsl_matrix_scale(m_k_du_sin_sin, -m);
      gsl_matrix_add(n_k_du_sin_cos, m_k_du_sin_sin);
      // add phase offset due to the position of the transmitter
      gsl_matrix_add_constant(n_k_du_sin_cos, alpha);
      // add phase offset due to coding
      gsl_matrix_add_constant(n_k_du_sin_cos, gsl_matrix_get(PHI, m, n));
      CMatrix complex_matrix = new_cmatrix(n_k_du_sin_cos);
      gsl_matrix_complex_scale(complex_matrix, C_0_M1j);
      exp(complex_matrix);
      // compute final phases e^-j(...)
      gsl_matrix_complex_add(phase, complex_matrix);

      gsl_matrix_complex_free(complex_matrix);
      gsl_matrix_free(n_k_du_sin_cos);
      gsl_matrix_free(m_k_du_sin_sin);
    }
  }

  // compute absolute value of all phases, i.e., length of all vectors
  // if length = 0 then all signals were interfering destructively
  // the longer the vector, the more constructively all the signals are summing
  // up together
  Matrix Fa = abs(phase);
  gsl_matrix_complex_free(phase);

  // P = Fa ^ 2
  // compute the power by squaring the absolute values
  matrix_copy(P, Fa);
  gsl_matrix_mul_elements(P, P);

  // (P^T * 1)
  Vector powers_by_theta = matrix_sum(P, true);
  // (P^T * 1) * A(THETA, d_theta, d_phi)
  p_tot = dot_product(powers_by_theta, spherical_elements);

  gsl_vector_free(powers_by_theta);
  gsl_matrix_free(Fa);

  recomputeGainMap = false;
  cached_phiTX = phiTX_rad;
  cached_thetaTX = thetaTX_rad;

  return cachedGain(phiRX_rad, thetaRX_rad);
}
