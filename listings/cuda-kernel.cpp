__global__ void gain_compute_phase_kernel(const double *k_du_sin_cos, const double *k_du_sin_sin, double n, double m, double alpha, double PHI, double *phase_real, double *phase_img, size_t size) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < size) {
    double tmp_img = -((-n * k_du_sin_cos[i]) + (-m * k_du_sin_sin[i]) + alpha + PHI);
    phase_real[i] += cos(tmp_img);
    phase_img[i] += sin(tmp_img);
  }
}
