//This is a classic test system of two coupled, nonlinear equations called
//the Brusselator (Brussels oscillator).
void brusselator (double t_, double *solin, double *fout) {
    (void)t_; //supress unused variable warning
    fout[0] = 1.0 + solin[0]*solin[0]*solin[1] - 4.0*solin[0];
    fout[1] = 3.0*solin[0] - solin[0]*solin[0]*solin[1];
}
