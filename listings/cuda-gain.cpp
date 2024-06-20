void ReconfigurableIntelligentSurface::gain_compute_phase(CMatrix phase, double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad) {
  int max_threads_per_block = withcuda::get_cuda_max_threads_per_block();
  withcuda::cuda_matrix cuda_k_du_sin_cos;
  withcuda::cuda_matrix cuda_k_du_sin_sin;
  withcuda::cuda_cmatrix cuda_phase;

  withcuda::cuda_matrix_alloc(cuda_k_du_sin_cos, k_du_sin_cos->size1, k_du_sin_cos->size2);
  withcuda::cuda_matrix_alloc(cuda_k_du_sin_sin, k_du_sin_sin->size1, k_du_sin_sin->size2);
  withcuda::cuda_cmatrix_alloc(cuda_phase, phase->size1, phase->size2);

  withcuda::gsl_matrix_to_cuda_matrix(cuda_k_du_sin_cos, k_du_sin_cos);
  withcuda::gsl_matrix_to_cuda_matrix(cuda_k_du_sin_sin, k_du_sin_sin);

  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      // compute the phase offset due to the incidence of the signal
      double alpha = du_k * (n * sin(thetaTX_rad) * cos(phiTX_rad) + m * sin(thetaTX_rad) * sin(phiTX_rad));
      double PHI = gsl_matrix_get(coding, m, n);
      // call the kernel function to compute the phase
      withcuda::gain_compute_phase(max_threads_per_block, cuda_k_du_sin_cos, cuda_k_du_sin_sin, n, m, alpha, PHI, cuda_phase);
    }
  }
  withcuda::cuda_matrix_free(cuda_k_du_sin_cos);
  withcuda::cuda_matrix_free(cuda_k_du_sin_sin);

  withcuda::cuda_cmatrix_to_gsl_cmatrix(phase, cuda_phase);
  withcuda::cuda_cmatrix_free(cuda_phase);
}
