void gain_compute_phase(int max_threads_per_block, const cuda_matrix &k_du_sin_cos, const cuda_matrix &k_du_sin_sin, double n, double m, double alpha, double PHI, cuda_cmatrix &phase) {
  int num_blocks = (k_du_sin_cos.rows * k_du_sin_cos.cols + max_threads_per_block - 1) / max_threads_per_block;
  gain_compute_phase_kernel<<<num_blocks, max_threads_per_block>>>( k_du_sin_cos.data, k_du_sin_sin.data, n, m, alpha, PHI, phase.data_real, phase.data_img, phase.rows * phase.cols);
}