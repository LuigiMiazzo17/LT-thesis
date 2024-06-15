double ReconfigurableIntelligentSurface::gain(double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad) {
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
      alpha = du_k * (n * sin(thetaTX_rad) * cos(phiTX_rad) + m * sin(thetaTX_rad) * sin(phiTX_rad)); Matrix n_k_du_sin_cos = new_matrix(k_du_sin_cos);
      Matrix m_k_du_sin_sin = new_matrix(k_du_sin_sin);
      gsl_matrix_scale(n_k_du_sin_cos, -n);
      gsl_matrix_scale(m_k_du_sin_sin, -m);
      gsl_matrix_add(n_k_du_sin_cos, m_k_du_sin_sin);
      gsl_matrix_add_constant(n_k_du_sin_cos, alpha);
      gsl_matrix_add_constant(n_k_du_sin_cos, gsl_matrix_get(PHI, m, n));
      CMatrix complex_matrix = new_cmatrix(n_k_du_sin_cos);
      gsl_matrix_complex_scale(complex_matrix, C_0_M1j);
      exp(complex_matrix);
      gsl_matrix_complex_add(phase, complex_matrix);

      gsl_matrix_complex_free(complex_matrix);
      gsl_matrix_free(n_k_du_sin_cos);
      gsl_matrix_free(m_k_du_sin_sin);
    }
  }

  Matrix Fa = abs(phase);
  gsl_matrix_complex_free(phase);

  matrix_copy(P, Fa);
  gsl_matrix_mul_elements(P, P);

  Vector powers_by_theta = matrix_sum(P, true);
  p_tot = dot_product(powers_by_theta, spherical_elements);

  gsl_vector_free(powers_by_theta);
  gsl_matrix_free(Fa);

  recomputeGainMap = false;
  cached_phiTX = phiTX_rad;
  cached_thetaTX = thetaTX_rad;

  return cachedGain(phiRX_rad, thetaRX_rad);
}
