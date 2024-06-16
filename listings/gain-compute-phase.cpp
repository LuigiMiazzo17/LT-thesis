void ReconfigurableIntelligentSurface::gain_compute_phase(CMatrix phase, double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad) {
  // number of threads to use
  int n_threads = this->n_max_threads;
  if (this->n_max_threads > this->M * this->N)
    n_threads = this->M * this->N;

  // parameters for the threads
  vector<thread_gain_args> t_args_list(n_threads);
  vector<pthread_t> gain_threads_list(n_threads);
  vector<CMatrix> phase_list(n_threads);

  // divide the computation of the gain in n_threads threads
  int jobs_per_thread = (N * M) / n_threads;
  int remainder = (N * M) % n_threads;
  int start = 0;

  // start the threads
  for (int i = 0; i < n_threads; i++) {
    thread_gain_args *t_args = &t_args_list[i];
    t_args->ris = this;
    t_args->thetaTX_rad = thetaTX_rad;
    t_args->phiTX_rad = phiTX_rad;
    t_args->thetaRX_rad = thetaRX_rad;
    t_args->phiRX_rad = phiRX_rad;
    t_args->tmp_phase = &phase_list[i];

    t_args->start = start;
    start = t_args->end = start + jobs_per_thread + (i < remainder ? 1 : 0);

    pthread_t ptid;
    pthread_create(&ptid, NULL, &gain_compute_phase_CPU_routine,
                   (void *)t_args);
    gain_threads_list[i] = ptid;
  }

  // wait for all threads to finish and sum up the results
  for (int i = 0; i < n_threads; i++) {
    pthread_join(gain_threads_list[i], NULL);
    gsl_matrix_complex_add(phase, phase_list[i]);
    gsl_matrix_complex_free(phase_list[i]);
  }
}
