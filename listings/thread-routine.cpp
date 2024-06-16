void *ReconfigurableIntelligentSurface::gain_compute_phase_CPU_routine(void *thread_args) {
  // extract the arguments
  struct thread_gain_args *t_args = (struct thread_gain_args *)thread_args;

  Matrix n_k_du_sin_cos = new_matrix(t_args->ris->k_du_sin_cos);
  Matrix m_k_du_sin_sin = new_matrix(t_args->ris->k_du_sin_sin);
  CMatrix complex_matrix = new_cmatrix(t_args->ris->Ts, t_args->ris->Ps);
  CMatrix tmp_phase = new_cmatrix(t_args->ris->Ts, t_args->ris->Ps);
  *(t_args->tmp_phase) = tmp_phase;

  for (int i = t_args->start; i < t_args->end; i++) {
    // extract m and n from the linerized index i
    int m = i / t_args->ris->N;
    int n = i % t_args->ris->N;

    // compute the phase offset due to the incidence of the signal
    double alpha = t_args->ris->du_k * (n * sin(t_args->thetaTX_rad) * cos(t_args->phiTX_rad) + m * sin(t_args->thetaTX_rad) * sin(t_args->phiTX_rad));
    // compute the phase offsets for all the possible phiRX,thetaRX pairs
    gsl_matrix_memcpy(n_k_du_sin_cos, t_args->ris->k_du_sin_cos);
    gsl_matrix_memcpy(m_k_du_sin_sin, t_args->ris->k_du_sin_sin);
    // scale by negative m and n, as these matrices need to be subtracted
    gsl_matrix_scale(n_k_du_sin_cos, -n);
    gsl_matrix_scale(m_k_du_sin_sin, -m);
    gsl_matrix_add(n_k_du_sin_cos, m_k_du_sin_sin);
    // add phase offset due to the position of the transmitter and add phase
    // offset due to coding
    gsl_matrix_add_constant(n_k_du_sin_cos, alpha + gsl_matrix_get(t_args->ris->coding, m, n));
    matrix_to_cmatrix(complex_matrix, n_k_du_sin_cos);
    gsl_matrix_complex_scale(complex_matrix, t_args->ris->C_0_M1j);
    // compute final phases e^-j(...)
    exp(complex_matrix);

    gsl_matrix_complex_add(tmp_phase, complex_matrix);
  }
  gsl_matrix_free(n_k_du_sin_cos);
  gsl_matrix_free(m_k_du_sin_sin);
  gsl_matrix_complex_free(complex_matrix);
  pthread_exit(NULL);
}
