struct ReconfigurableIntelligentSurface::thread_gain_args {
    ReconfigurableIntelligentSurface* ris;
    double thetaTX_rad;
    double phiTX_rad;
    double thetaRX_rad;
    double phiRX_rad;
    int start;
    int end;
    CMatrix* tmp_phase;
};
