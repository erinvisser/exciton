#include "transition_rates.hh"
#include "preeq_common.hh"
#include "state_density.hh"

#include "marley/Integrator.hh"
#include "marley/KoningDelarocheOpticalModel.hh"
#include "marley/MassTable.hh"

#include <cmath>
#include <algorithm>

// Volume-averaged imaginary potential: matching TALYS bonetti.f90:104-122
double wvolRadialIntegral(double Wv, double Wd, double Rv, double av,
                          double Rd, double ad)
{
    constexpr int nrbins = 50;
    constexpr double dr = 20.0 / nrbins;

    double sum1 = 0.0, sum2 = 0.0;
    for (int i = 1; i <= nrbins; ++i)
    {
        double rr = (i - 0.5) * dr;

        double fwsvol = 0.0;
        double expo_v = (rr - Rv) / av;
        if (expo_v <= 80.0)
            fwsvol = 1.0 / (1.0 + std::exp(expo_v));

        double fwssurf = 0.0;
        double expo_d = (rr - Rd) / ad;
        if (expo_d <= 80.0)
        {
            double exp_d = std::exp(expo_d);
            fwssurf = -exp_d / (ad * (1.0 + exp_d) * (1.0 + exp_d));
        }

        double term2 = fwsvol * rr * rr * dr;
        double term1 = term2 * (Wv * fwsvol - 4.0 * ad * Wd * fwssurf);
        sum1 += term1;
        sum2 += term2;
    }

    return (sum2 != 0.0) ? sum1 / sum2 : 0.0;
}

double M2(int A_p, int A, double E_tot, int n, double C1, double C2, double C3, PreeqMode mode)
{
    // Equation (13.26)
    double prefactor = C1 * A_p / std::pow(A, 3);
    double denom = 10.7 * C3 + E_tot / (n * A_p);
    double inside = 7.48 * C2 + 4.62e5 / std::pow(denom, 3);
    double result = prefactor * inside;
    if (mode == PreeqMode::Analytical)
    {
        return 1.20 * result;
    }
    return result;
}

double M2_element(MatrixElement type, int A_p, int A, double E_tot, int n,
                  double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                  double C1, double C2, double C3, PreeqMode mode)
{
    // Equation (13.24)
    double matrix_squared = M2(A_p, A, E_tot, n, C1, C2, C3, mode);
    switch (type)
    {
    case MatrixElement::NuNu:
        return R_nu_nu * matrix_squared;
    case MatrixElement::NuPi:
        return R_nu_pi * matrix_squared;
    case MatrixElement::PiPi:
        return R_pi_pi * matrix_squared;
    default:
        return R_pi_nu * matrix_squared;
    }
}

double B(int p_pi, int h_pi, int p_nu, int h_nu, double g_p, double g_n, ConversionType type)
{
    // Equation (13.28)
    double A_current = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n);

    double A_next;

    if (type == ConversionType::ProtonToNeutron)
    {
        A_next = pauliCorrection(p_pi - 1, h_pi - 1, p_nu + 1, h_nu + 1, g_p, g_n);
    }
    else
    {
        A_next = pauliCorrection(p_pi + 1, h_pi + 1, p_nu - 1, h_nu - 1, g_p, g_n);
    }

    return std::max(A_current, A_next);
}

double lambdaNewPairAnalytical(ExcitonType particle, int Z, int N, int A_p,
                               int p_pi, int h_pi, int p_nu, int h_nu,
                               double E_tot, double U, double V,
                               double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                               double C1, double C2, double C3)
{
    // first half of Equation (13.27) — TALYS lambdapiplus / lambdanuplus
    int n_proton = p_pi + h_pi;
    int n_neutron = p_nu + h_nu;
    int p = p_pi + p_nu;
    int h = h_pi + h_nu;
    int n = p + h;
    int A_compound = Z + N + A_p;

    double g_p = spDensityProton(Z);
    int N_comp = N + A_p; // compound N (Z_proj=0 for neutron projectile; N_comp = N + 1)
    double g_n = spDensityNeutron(N_comp);

    double prefactor = 2 * PI / HBAR;
    double current_energy = U - pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n);
    if (current_energy <= 0.0)
        return 0.0;

    double well_term = finiteWell(p + 1, h + 1, U, V);

    double denom = 2 * n * (n + 1) * std::pow(current_energy, n - 1);
    if (denom <= 0)
        return 0.0;

    if (particle == ExcitonType::Proton)
    {
        double next_energy = U - pauliCorrection(p_pi + 1, h_pi + 1, p_nu, h_nu, g_p, g_n);
        if (next_energy <= 0.0)
            return 0.0;
        double ratio = next_energy / current_energy;
        if (ratio < 0.01)
            return 0.0;
        double numerator = std::pow(g_p, 2) * std::pow(next_energy, n + 1);
        double matrix_term = n_proton * g_p * M2_element(MatrixElement::PiPi, A_p, A_compound, E_tot, n) + 2 * n_neutron * g_n * M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n);
        return prefactor * numerator * matrix_term * well_term / denom;
    }
    else
    {
        double next_energy = U - pauliCorrection(p_pi, h_pi, p_nu + 1, h_nu + 1, g_p, g_n);
        if (next_energy <= 0.0)
            return 0.0;
        double ratio = next_energy / current_energy;
        if (ratio < 0.01)
            return 0.0;
        double numerator = std::pow(g_n, 2) * std::pow(next_energy, n + 1);
        double matrix_term = n_neutron * g_n * M2_element(MatrixElement::NuNu, A_p, A_compound, E_tot, n) + 2 * n_proton * g_p * M2_element(MatrixElement::PiNu, A_p, A_compound, E_tot, n);
        return prefactor * numerator * matrix_term * well_term / denom;
    }
}

double lambdaNewPairNumerical(ExcitonType particle, int Z, int N, int A_p,
                              int p_pi, int h_pi, int p_nu, int h_nu,
                              double E_tot, double U, double V,
                              double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                              double C1, double C2, double C3, double M2constant,
                              IntegrationMethod method, int midBins,
                              bool guardBounds, CollisionKernel kernel,
                              int Z_proj)
{
    int n = p_pi + h_pi + p_nu + h_nu;
    double prefactor = 2 * PI / HBAR;

    int A_target = Z + N;
    int A_compound = A_target + A_p;

    int Z_comp = Z + Z_proj;
    int N_comp = N + (A_p - Z_proj);
    double g_p = spDensityProton(Z_comp);
    double g_n = spDensityNeutron(N_comp);
    const auto &mt = marley::MassTable::Instance();
    double sep_p = mt.get_fragment_separation_energy(
        Z_comp, A_compound, 2212);
    double sep_n = mt.get_fragment_separation_energy(
        Z_comp, A_compound, 2112);

    double wompfac_same = 0.0, wompfac_cross = 0.0;
    if (kernel == CollisionKernel::OpticalModel)
    {
        constexpr double C_OMP = 0.55; // TODO: potentially a configurable later on.
        double M2c = M2constant;
        double Rpinu = R_pi_nu / R_pi_pi;
        double denom = 1.0 + 2.0 * Rpinu;
        wompfac_same = M2c * C_OMP / denom;
        wompfac_cross = M2c * C_OMP * 2.0 * Rpinu / denom;
    }

    // turning off pairing bc passing in U, not E_x
    double omega_initial = particleHoleStateDensity(
        p_pi, h_pi, p_nu, h_nu, U, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
    if (omega_initial <= 0.0)
        return 0.0;

    marley::KoningDelarocheOpticalModel kd_omp(Z, A_target);
    marley::Integrator integrator(50);
    double sum = 0.0;

    auto integrate_sub = [&](int t_pi, int t_hi, int t_pn, int t_hn,
                             int r_pi, int r_hi, int r_pn, int r_hn,
                             double M2_val,
                             double g_factor,
                             double L1, double L2,
                             double wompfac = 0.0,
                             int ref_pi = 0, int ref_hi = 0,
                             int ref_pn = 0, int ref_hn = 0,
                             int collider_pdg = 0) -> double
    {
        auto integrand = [&](double e) -> double
        {
            double omega_r = particleHoleStateDensity(
                r_pi, r_hi, r_pn, r_hn, U - e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
            if (kernel == CollisionKernel::MatrixElement)
            {
                double omega_t = particleHoleStateDensity(
                    t_pi, t_hi, t_pn, t_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                return prefactor * M2_val * omega_t * g_factor * omega_r;
            }
            else
            {
                double e_kin = (collider_pdg == 2212) ? e - sep_p : e - sep_n;
                if (e_kin < -20.0)
                    e_kin = -20.0;
                kd_omp.setIncidentEnergyAndFragment(e_kin, collider_pdg);
                double Wv = kd_omp.getWv();
                double Wd = kd_omp.getWd();
                double Rv = kd_omp.getRv();
                double av = kd_omp.getav();
                double Rd = kd_omp.getRd();
                double ad = kd_omp.getad();
                double wvol = wvolRadialIntegral(Wv, Wd, Rv, av, Rd, ad);
                static bool diagnostic_printed = false;
                if (!diagnostic_printed)
                {
                    std::cerr << "state=" << p_pi << "," << h_pi << "," << p_nu << "," << h_nu
                              << " term=" << t_pi << "," << t_hi << "," << t_pn << "," << t_hn
                              << " e=" << e << " e_kin=" << e_kin
                              << " sep_p=" << sep_p << " sep_n=" << sep_n
                              << " Wv=" << Wv << " Wd=" << Wd
                              << " Rv=" << Rv << " av=" << av
                              << " Rd=" << Rd << " ad=" << ad
                              << " wvol=" << wvol
                              << " wompfac=" << wompfac << "\n";
                    diagnostic_printed = true;
                }
                double collision = (2.0 / HBAR) * wompfac * wvol;
                if (ref_pi != 0 || ref_hi != 0 || ref_pn != 0 || ref_hn != 0)
                {
                    double omega_hole = particleHoleStateDensity(
                        t_pi, t_hi, t_pn, t_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                    double omega_particle = particleHoleStateDensity(
                        ref_pi, ref_hi, ref_pn, ref_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                    if (omega_particle > 1.0)
                        collision *= omega_hole / omega_particle;
                }
                return collision * g_factor * omega_r;
            }
        };

        if (method == IntegrationMethod::ClenshawCurtis)
        {
            if (guardBounds)
            {
                double guard_t = std::pow(std::max(t_pi, t_hi), 2) / g_p + std::pow(std::max(t_pn, t_hn), 2) / g_n;
                double guard_r = std::pow(std::max(r_pi, r_hi), 2) / g_p + std::pow(std::max(r_pn, r_hn), 2) / g_n;
                if (L1 < guard_t)
                    L1 = guard_t;
                if (L2 > U - guard_r)
                    L2 = U - guard_r;
            }
            if (L2 <= L1)
                return 0.0;
            return integrator.num_integrate(integrand, L1, L2);
        }
        else
        {
            if (L2 <= L1)
                return 0.0;
            double dx = (L2 - L1) / midBins;
            double sum = 0.0;
            for (int i = 0; i < midBins; ++i)
                sum += integrand(L1 + (i + 0.5) * dx);
            return sum * dx;
        }
    };

    if (particle == ExcitonType::Proton)
    {
        // Equation (13.18)
        double new_pi_pair = pauliCorrection(p_pi + 1, h_pi + 1, p_nu, h_nu, g_p, g_n);

        double less_pi_part = pauliCorrection(p_pi - 1, h_pi, p_nu, h_nu, g_p, g_n);
        double less_pi_hole = pauliCorrection(p_pi, h_pi - 1, p_nu, h_nu, g_p, g_n);

        double less_nu_part = pauliCorrection(p_pi, h_pi, p_nu - 1, h_nu, g_p, g_n);
        double less_nu_hole = pauliCorrection(p_pi, h_pi, p_nu, h_nu - 1, g_p, g_n);

        double L1pip = new_pi_pair - less_pi_part;
        double L2pip = U - less_pi_part;

        double L1pih = new_pi_pair - less_pi_hole;
        double L2pih = U - less_pi_hole;

        double L1nup = new_pi_pair - less_nu_part;
        double L2nup = U - less_nu_part;

        double L1nuh = new_pi_pair - less_nu_hole;
        double L2nuh = U - less_nu_hole;

        double m2_pipi = M2_element(MatrixElement::PiPi, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);
        double m2_nupi = M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        // Equation (13.17)
        sum += integrate_sub(2, 1, 0, 0, p_pi - 1, h_pi, p_nu, h_nu, m2_pipi, g_p, L1pip, L2pip,
                             wompfac_same, 0, 0, 0, 0, 2212);
        sum += integrate_sub(1, 2, 0, 0, p_pi, h_pi - 1, p_nu, h_nu, m2_pipi, g_p, L1pih, L2pih,
                             wompfac_same, 2, 1, 0, 0, 2212);
        sum += integrate_sub(1, 1, 1, 0, p_pi, h_pi, p_nu - 1, h_nu, m2_nupi, g_n, L1nup, L2nup,
                             wompfac_cross, 0, 0, 0, 0, 2112);
        sum += integrate_sub(1, 1, 0, 1, p_pi, h_pi, p_nu, h_nu - 1, m2_nupi, g_n, L1nuh, L2nuh,
                             wompfac_cross, 1, 1, 1, 0, 2112);
    }
    else
    {
        // Equation (13.18)
        double new_nu_pair = pauliCorrection(p_pi, h_pi, p_nu + 1, h_nu + 1, g_p, g_n);

        double less_pi_part = pauliCorrection(p_pi - 1, h_pi, p_nu, h_nu, g_p, g_n);
        double less_pi_hole = pauliCorrection(p_pi, h_pi - 1, p_nu, h_nu, g_p, g_n);

        double less_nu_part = pauliCorrection(p_pi, h_pi, p_nu - 1, h_nu, g_p, g_n);
        double less_nu_hole = pauliCorrection(p_pi, h_pi, p_nu, h_nu - 1, g_p, g_n);

        double L1nup = new_nu_pair - less_nu_part;
        double L2nup = U - less_nu_part;

        double L1nuh = new_nu_pair - less_nu_hole;
        double L2nuh = U - less_nu_hole;

        double L1pip = new_nu_pair - less_pi_part;
        double L2pip = U - less_pi_part;

        double L1pih = new_nu_pair - less_pi_hole;
        double L2pih = U - less_pi_hole;

        double m2_nunu = M2_element(MatrixElement::NuNu, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);
        double m2_pinu = M2_element(MatrixElement::PiNu, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        sum += integrate_sub(0, 0, 2, 1, p_pi, h_pi, p_nu - 1, h_nu, m2_nunu, g_n, L1nup, L2nup,
                             wompfac_same, 0, 0, 0, 0, 2112);
        sum += integrate_sub(0, 0, 1, 2, p_pi, h_pi, p_nu, h_nu - 1, m2_nunu, g_n, L1nuh, L2nuh,
                             wompfac_same, 0, 0, 2, 1, 2112);
        sum += integrate_sub(1, 0, 1, 1, p_pi - 1, h_pi, p_nu, h_nu, m2_pinu, g_p, L1pip, L2pip,
                             wompfac_cross, 0, 0, 0, 0, 2212);
        sum += integrate_sub(0, 1, 1, 1, p_pi, h_pi - 1, p_nu, h_nu, m2_pinu, g_p, L1pih, L2pih,
                             wompfac_cross, 1, 0, 1, 1, 2212);
    }

    // Equation (13.18)
    return sum / omega_initial;
}

double lambdaPairConversionAnalytical(ConversionType conversion, int Z, int N, int A_p,
                                      int p_pi, int h_pi, int p_nu, int h_nu,
                                      double E_tot, double U, double V,
                                      double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                                      double C1, double C2, double C3)
{
    // Second half of Equation (13.27) — TALYS lambdapinu / lambdanupi
    int n_proton = p_pi + h_pi;
    int n_neutron = p_nu + h_nu;
    int p = p_pi + p_nu;
    int h = h_pi + h_nu;
    int n = p + h;
    int A_compound = Z + N + A_p;

    double g_p = spDensityProton(Z);
    int N_comp = N + A_p;
    double g_n = spDensityNeutron(N_comp);

    double pauli = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n);

    double prefactor = 2 * PI / HBAR;

    double denom = U - pauli;
    if (denom <= 0)
        return 0.0;

    if (conversion == ConversionType::ProtonToNeutron)
    {
        double front = M2_element(MatrixElement::PiNu, A_p, A_compound, E_tot, n) * p_pi * h_pi * std::pow(g_n, 2) * finiteWell(p, h, U, V) / n;
        double B_correct = U - B(p_pi, h_pi, p_nu, h_nu, g_p, g_n, ConversionType::ProtonToNeutron);
        if (B_correct <= 0.0)
            return 0.0;
        double frac = std::pow(B_correct / denom, n - 1);
        if (frac < 0.01)
            return 0.0;
        double factor4 = 2.0 * B_correct;
        if (p_pi > 0 && h_pi > 0)
            factor4 += n * std::abs(pauli - pauliCorrection(p_pi - 1, h_pi - 1, p_nu + 1, h_nu + 1, g_p, g_n));
        return prefactor * front * frac * factor4;
    }
    else
    {
        double front = M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n) * p_nu * h_nu * std::pow(g_p, 2) * finiteWell(p, h, U, V) / n;
        double B_correct = U - B(p_pi, h_pi, p_nu, h_nu, g_p, g_n, ConversionType::NeutronToProton);
        if (B_correct <= 0.0)
            return 0.0;
        double frac = std::pow(B_correct / denom, n - 1);
        if (frac < 0.01)
            return 0.0;
        double factor4 = 2.0 * B_correct;
        if (p_nu > 0 && h_nu > 0)
            factor4 += n * std::abs(pauli - pauliCorrection(p_pi + 1, h_pi + 1, p_nu - 1, h_nu - 1, g_p, g_n));
        return prefactor * front * frac * factor4;
    }
}

double lambdaPairConversionNumerical(ConversionType conversion, int Z, int N, int A_p,
                                     int p_pi, int h_pi, int p_nu, int h_nu,
                                     double E_tot, double U, double V,
                                     double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                                     double C1, double C2, double C3, double M2constant,
                                     IntegrationMethod method, int midBins,
                                     bool guardBounds, CollisionKernel kernel,
                                     int Z_proj)
{
    int n = p_pi + h_pi + p_nu + h_nu;

    double prefactor = 2 * PI / HBAR;

    int A_target = Z + N;
    int A_compound = A_target + A_p;

    int Z_comp = Z + Z_proj;
    int N_comp = N + (A_p - Z_proj);
    double g_p = spDensityProton(Z_comp);
    double g_n = spDensityNeutron(N_comp);
    const auto &mt = marley::MassTable::Instance();
    double sep_p_res = mt.get_fragment_separation_energy(
        Z_comp - 1, A_compound - 1, 2212);
    double sep_n_res = mt.get_fragment_separation_energy(
        Z_comp, A_compound - 1, 2112);

    double wompfac_conv = 0.0;
    if (kernel == CollisionKernel::OpticalModel)
    {
        constexpr double C_OMP = 0.55;
        double M2c = M2constant;
        double Rpinu = R_pi_nu / R_pi_pi;
        double denom = 1.0 + 2.0 * Rpinu;
        double wompfac_same = M2c * C_OMP / denom;
        double wompfac_cross = M2c * C_OMP * 2.0 * Rpinu / denom;
        wompfac_conv = 0.5 * (wompfac_same + wompfac_cross);
    }

    double omega_initial = particleHoleStateDensity(
        p_pi, h_pi, p_nu, h_nu, U, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
    if (omega_initial <= 0.0)
        return 0.0;

    marley::KoningDelarocheOpticalModel kd_omp(Z_comp, A_compound);
    marley::Integrator integrator(50);

    auto integrate_conversion = [&](int c_pi, int c_hi, int c_pn, int c_hn,
                                    int a_pi, int a_hi, int a_pn, int a_hn,
                                    int r_pi, int r_hi, int r_pn, int r_hn,
                                    double M2_val,
                                    double L1, double L2,
                                    double wompfac = 0.0,
                                    int ref_pi = 0, int ref_hi = 0,
                                    int ref_pn = 0, int ref_hn = 0,
                                    int collider_pdg = 0) -> double
    {
        auto integrand = [&](double e) -> double
        {
            double omega_r = particleHoleStateDensity(
                r_pi, r_hi, r_pn, r_hn, U - e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
            if (kernel == CollisionKernel::MatrixElement)
            {
                double omega_c = particleHoleStateDensity(
                    c_pi, c_hi, c_pn, c_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                double omega_a = particleHoleStateDensity(
                    a_pi, a_hi, a_pn, a_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                return prefactor * M2_val * omega_c * omega_a * omega_r;
            }
            else
            {
                double omega_a = particleHoleStateDensity(
                    a_pi, a_hi, a_pn, a_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                double e_kin = (collider_pdg == 2212) ? e - sep_p_res : e - sep_n_res;
                if (e_kin < -20.0)
                    e_kin = -20.0;
                kd_omp.setIncidentEnergyAndFragment(e_kin, collider_pdg);
                double Wv = kd_omp.getWv();
                double Wd = kd_omp.getWd();
                double Rv = kd_omp.getRv();
                double av = kd_omp.getav();
                double Rd = kd_omp.getRd();
                double ad = kd_omp.getad();
                double wvol = wvolRadialIntegral(Wv, Wd, Rv, av, Rd, ad);
                static bool conv_diag_printed = false;
                if (!conv_diag_printed)
                {
                    std::cerr << "CONV state=" << p_pi << "," << h_pi << "," << p_nu << "," << h_nu
                              << " term=" << c_pi << "," << c_hi << "," << c_pn << "," << c_hn
                              << " e=" << e << " e_kin=" << e_kin
                              << " sep_p=" << sep_p_res << " sep_n=" << sep_n_res
                              << " Wv=" << Wv << " Wd=" << Wd
                              << " Rv=" << Rv << " av=" << av
                              << " Rd=" << Rd << " ad=" << ad
                              << " wvol=" << wvol
                              << " wompfac=" << wompfac << "\n";
                    conv_diag_printed = true;
                }
                double collision = (2.0 / HBAR) * wompfac * wvol;
                if (ref_pi != 0 || ref_hi != 0 || ref_pn != 0 || ref_hn != 0)
                {
                    double omega_created = particleHoleStateDensity(
                        c_pi, c_hi, c_pn, c_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                    double omega_ref = particleHoleStateDensity(
                        ref_pi, ref_hi, ref_pn, ref_hn, e, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
                    if (omega_ref > 1.0)
                        collision *= omega_created / omega_ref;
                }
                return collision * omega_a * omega_r;
            }
        };

        if (method == IntegrationMethod::ClenshawCurtis)
        {
            if (guardBounds)
            {
                double guard_c = std::pow(std::max(c_pi, c_hi), 2) / g_p + std::pow(std::max(c_pn, c_hn), 2) / g_n;
                double guard_a = std::pow(std::max(a_pi, a_hi), 2) / g_p + std::pow(std::max(a_pn, a_hn), 2) / g_n;
                double guard_r = std::pow(std::max(r_pi, r_hi), 2) / g_p + std::pow(std::max(r_pn, r_hn), 2) / g_n;
                if (L1 < guard_c)
                    L1 = guard_c;
                if (L1 < guard_a)
                    L1 = guard_a;
                if (L2 > U - guard_r)
                    L2 = U - guard_r;
            }
            if (L2 <= L1)
                return 0.0;
            return integrator.num_integrate(integrand, L1, L2);
        }
        else
        {
            if (L2 <= L1)
                return 0.0;
            double dx = (L2 - L1) / midBins;
            double sum = 0.0;
            for (int i = 0; i < midBins; ++i)
                sum += integrand(L1 + (i + 0.5) * dx);
            return sum * dx;
        }
    };

    if (conversion == ConversionType::ProtonToNeutron)
    {
        double less_proton_pair = pauliCorrection(p_pi - 1, h_pi - 1, p_nu, h_nu, g_p, g_n);
        double L1 = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n) - less_proton_pair;
        double L2 = U - less_proton_pair;

        double m2_pinu = M2_element(MatrixElement::PiNu, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        double integral = integrate_conversion(
            0, 0, 1, 1, 1, 1, 0, 0, p_pi - 1, h_pi - 1, p_nu, h_nu, m2_pinu, L1, L2,
            wompfac_conv, 1, 0, 1, 1, 2212);
        return integral / omega_initial;
    }
    else
    {
        double less_neutron_pair = pauliCorrection(p_pi, h_pi, p_nu - 1, h_nu - 1, g_p, g_n);
        double L1 = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n) - less_neutron_pair;
        double L2 = U - less_neutron_pair;

        double m2_nupi = M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n,
                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        double integral = integrate_conversion(
            1, 1, 0, 0, 0, 0, 1, 1, p_pi, h_pi, p_nu - 1, h_nu - 1, m2_nupi, L1, L2,
            wompfac_conv, 1, 1, 1, 0, 2112);
        return integral / omega_initial;
    }
}

double lambdaRate(LambdaType type, PreeqMode mode, int Z, int N, int A_p,
                  int p_pi, int h_pi, int p_nu, int h_nu,
                  double E_tot, double U, double V,
                  double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                  double C1, double C2, double C3, double M2constant,
                  IntegrationMethod method, int midBins,
                  bool guardBounds, CollisionKernel kernel,
                  int Z_proj)
{
    int n = p_pi + h_pi + p_nu + h_nu;
    bool useAnalytical = (mode == PreeqMode::Analytical) || (n == 1);

    double result;

    switch (type)
    {
    case LambdaType::ProtonPairCreation:
        if (useAnalytical)
            result = lambdaNewPairAnalytical(ExcitonType::Proton, Z, N, A_p,
                                             p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                             R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaNewPairNumerical(ExcitonType::Proton, Z, N, A_p,
                                            p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                            R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2constant,
                                            method, midBins, guardBounds, kernel, Z_proj);
        break;

    case LambdaType::NeutronPairCreation:
        if (useAnalytical)
            result = lambdaNewPairAnalytical(ExcitonType::Neutron, Z, N, A_p,
                                             p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                             R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaNewPairNumerical(ExcitonType::Neutron, Z, N, A_p,
                                            p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                            R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2constant,
                                            method, midBins, guardBounds, kernel, Z_proj);
        break;

    case LambdaType::ProtonToNeutronConversion:
        if (useAnalytical)
            result = lambdaPairConversionAnalytical(ConversionType::ProtonToNeutron, Z, N, A_p,
                                                    p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaPairConversionNumerical(ConversionType::ProtonToNeutron, Z, N, A_p,
                                                   p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                   R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2constant,
                                                   method, midBins, guardBounds, kernel, Z_proj);
        break;

    default:
        if (useAnalytical)
            result = lambdaPairConversionAnalytical(ConversionType::NeutronToProton, Z, N, A_p,
                                                    p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaPairConversionNumerical(ConversionType::NeutronToProton, Z, N, A_p,
                                                   p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                   R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2constant,
                                                   method, midBins, guardBounds, kernel, Z_proj);
        break;
    }

    // TALYS matrix.f90 applies 1.20 M² factor only for preeqmode==1.
    // When n==1 forces analytical formula in numerical mode, cancel the factor.
    if (n == 1 && mode == PreeqMode::Numerical)
        result /= 1.20;

    return result;
}