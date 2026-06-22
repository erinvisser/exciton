#ifndef STATE_DENSITY_HH
#define STATE_DENSITY_HH

double factorial(int n);

double theta(double x);

double particleHoleStateDensity(int p_pi, int h_pi, int p_nu, int h_nu, double E_x, int Z, int N, int chi, double V, double n_crit, bool applyPairing = true, double parityFactor = 1.0);

#endif
