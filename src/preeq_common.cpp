#include "preeq_common.hh"

#include <algorithm>
#include <cmath>

double spDensityProton(int Z)
{
    // Equation (13.9)
    return Z / 13.0;
}

double spDensityNeutron(int N)
{
    // Equation (13.9)
    return N / 13.0;
}

double binomial(int h, int i)
{
    if (i < 0 || i > h)
        return 0.0;
    if (i == 0 || i == h)
        return 1.0;

    int n = std::min(i, h - i);

    double result = 1;
    for (int j = 0; j < n; ++j)
    {
        result *= (h - j);
        result /= (j + 1);
    }
    return result;
}

double potentialWellDepth(ProjectileType incident, int A, double E_p, int h)
{
    // Equation (13.11)
    if (h > 1)
        return 38.0;

    double E4 = std::pow(E_p, 4);
    double denom = 1.0;

    if (incident == ProjectileType::Proton)
    {
        denom = std::pow(450 / std::cbrt(A), 4);
        return 22.0 + 16.0 * E4 / (E4 + denom);
    }
    else if (incident == ProjectileType::Neutron)
    {
        denom = std::pow(250 / std::cbrt(A), 4);
        return 12.0 + 26.0 * E4 / (E4 + denom);
    }
    else
    {
        return 38.0;
    }
}

double pauliCorrection(int p_pi, int h_pi, int p_nu, int h_nu, double g_p, double g_n)
{
    // Equation (13.8)
    double term_1 = std::pow(std::max(p_pi, h_pi), 2) / g_p;
    double term_2 = std::pow(std::max(p_nu, h_nu), 2) / g_n;
    double term_3 = (std::pow(p_pi, 2) + std::pow(h_pi, 2) + p_pi + h_pi) / (4 * g_p);
    double term_4 = (std::pow(p_nu, 2) + std::pow(h_nu, 2) + p_nu + h_nu) / (4 * g_n);
    return term_1 + term_2 - term_3 - term_4;
}

double finiteWell(int p, int h, double U, double V)
{
    // Equation (13.10)
    if (U <= 0)
        return 1.0;

    int n = p + h;

    if (n < 1)
        return 1.0;

    double f = 1.0;

    for (int i = 1; i <= h; i++)
    {
        double sign = (i % 2 == 0) ? 1.0 : -1.0;
        double diff = U - i * V;

        if (diff >= 0.0)
        {
            f += sign * binomial(h, i) * std::pow(diff / U, n - 1);
        }
    }

    return f;
}

double fuPairingCorrection(double delta, double E_x, int n, double n_crit)
{
    // Equation (13.5)
    if (delta <= 0.0 || E_x <= 0.0)
        return 0.0;

    if (n_crit < 1)
        return delta;

    double energy_ratio = E_x / delta;
    double n_ratio = static_cast<double>(n) / n_crit;

    double threshold = 0.716 + 2.44 * std::pow(n_ratio, 2.17);
    if (energy_ratio >= threshold)
    {
        double base = 0.996 - 1.76 * std::pow(n_ratio, 1.6) / std::pow(energy_ratio, 0.68);
        return delta - delta * std::pow(base, 2);
    }
    else
    {
        return delta;
    }
}

double pairingEnergyDelta(int A, int chi)
{
    // Equation (13.7)
    // chi = 0 for odd-odd, 1 for odd, 2 for even-even
    return chi * 12.0 / std::sqrt(A);
}

double nCrit(int Z, int N, double delta)
{
    // Equation (13.6)
    double g = spDensityProton(Z) + spDensityNeutron(N);

    if (g <= 0.0 || delta <= 0.0)
        return 1;

    double T_crit = 2 * std::sqrt(delta / (0.25 * g)) / 3.5;
    double result = 2 * g * T_crit * std::log(2);

    return std::max(result, 1.0);
}

double availableExcitationEnergy(double E_x, double P_ph)
{
    return E_x - P_ph;
}

double spinDistribution(double J, int A, int n, bool talysCompat)
{
    // Equation (13.59) and (13.61)
    double sigma2 = 0.24 * n * std::pow(A, 2.0 / 3.0);
    double sigma = std::sqrt(sigma2);
    int k = talysCompat ? 2 : n;
    double denom = std::sqrt(PI) * std::pow(static_cast<double>(k), 3.0 / 2.0) * std::pow(sigma, 3.0);
    double expo = -std::pow(J + 0.5, 2) / (k * sigma2);
    return (2.0 * J + 1.0) / denom * std::exp(expo);
}