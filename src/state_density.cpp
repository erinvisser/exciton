#include "state_density.hh"
#include "preeq_common.hh"

#include <algorithm>
#include <cmath>

double factorial(int n)
{
    return std::tgamma(n + 1.0);
}

double theta(double x)
{
    return (x >= 0.0) ? 1.0 : 0.0;
}

double particleHoleStateDensity(int p_pi, int h_pi, int p_nu, int h_nu, double E_x, int Z, int N, int chi, double V, double n_crit, bool applyPairing, double parityFactor)
{
    // Equation (13.4)
    int n_pi = p_pi + h_pi;
    int n_nu = p_nu + h_nu;
    int n = n_pi + n_nu;
    int p = p_pi + p_nu;
    int h = h_pi + h_nu;

    if (n < 1)
        return 0.0;
    if (E_x <= 0.0)
        return 0.0;

    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);

    double pair_energy = 0.0;
    if (applyPairing)
    {
        int A = Z + N;
        double delta = pairingEnergyDelta(A, chi);
        double P_ph = fuPairingCorrection(delta, E_x, n, n_crit);
        pair_energy = availableExcitationEnergy(E_x, P_ph);
        if (pair_energy <= 0.0)
            return 0.0;
    }
    else
    {
        pair_energy = E_x;
    }

    double A_pauli = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n);

    double E_avail = pair_energy - A_pauli;
    if (E_avail <= 0.0)
        return 0.0;

    double g_product = std::pow(g_p, n_pi) * std::pow(g_n, n_nu);
    double factorial_product = factorial(p_pi) * factorial(h_pi) * factorial(p_nu) * factorial(h_nu) * factorial(n - 1);

    double term_1 = g_product / factorial_product;
    double term_2 = std::pow(E_avail, n - 1);
    double term_3 = finiteWell(p, h, applyPairing ? pair_energy : E_x, V);

    return parityFactor * term_1 * term_2 * term_3;
}
