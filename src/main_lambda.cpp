#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstring>

#include "preeq_common.hh"
#include "transition_rates.hh"
#include "marley/KoningDelarocheOpticalModel.hh"

int main(int argc, char *argv[])
{
    PreeqMode mode = PreeqMode::Numerical;
    IntegrationMethod method = IntegrationMethod::ClenshawCurtis;
    int midBins = 50;
    bool guardBounds = true;
    CollisionKernel kernel = CollisionKernel::MatrixElement;

    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--analytical") == 0)
            mode = PreeqMode::Analytical;
        else if (std::strcmp(argv[i], "--no-guards") == 0)
            guardBounds = false;
        else if (std::strcmp(argv[i], "--method") == 0 && i + 1 < argc)
        {
            ++i;
            if (std::strcmp(argv[i], "midpoint") == 0)
                method = IntegrationMethod::Midpoint;
        }
        else if (std::strcmp(argv[i], "--midbins") == 0 && i + 1 < argc)
        {
            midBins = std::atoi(argv[++i]);
        }
        else if (std::strcmp(argv[i], "--kernel") == 0 && i + 1 < argc)
        {
            ++i;
            if (std::strcmp(argv[i], "optical-model") == 0)
                kernel = CollisionKernel::OpticalModel;
            else if (std::strcmp(argv[i], "matrix-element") == 0)
                kernel = CollisionKernel::MatrixElement;
        }
    }

    if (kernel == CollisionKernel::OpticalModel && mode == PreeqMode::Analytical)
    {
        std::cerr << "# warning: OMP kernel requires numerical mode; overriding --analytical\n";
        mode = PreeqMode::Numerical;
    }

    // Diagnostic: OMP parameters matching TALYS bonetti.f90 energy grid
    {
        marley::KoningDelarocheOpticalModel kd_target(18, 40);
        auto print_omp = [&](double e, int k)
        {
            int pdg = (k == 2) ? 2212 : 2112;
            double e_safe = std::max(e, 0.0);
            kd_target.setIncidentEnergyAndFragment(e_safe, pdg);
            double wv = kd_target.getWv();
            double wd = kd_target.getWd();
            double rv = kd_target.getRv();
            double av = kd_target.getav();
            double rd = kd_target.getRd();
            double ad = kd_target.getad();
            double wvol = wvolRadialIntegral(wv, wd, rv, av, rd, ad);
            std::cout << "# OUR_OMP k=" << k << " e=" << e
                      << " w=" << wv << " wd=" << wd
                      << " rw=" << rv << " aw=" << av
                      << " rwd=" << rd << " awd=" << ad
                      << " wvol=" << wvol << "\n";
        };
        for (int nen = -80; nen <= 200; ++nen)
        {
            double e = 0.1 * nen;
            bool print = (nen == -80 || nen == -70 || nen == -60 || nen == 0 || nen == 10 || nen == 20 || (nen > 0 && nen % 20 == 0));
            if (print)
            {
                print_omp(e, 1);
                print_omp(e, 2);
            }
        }
    }

    int Z = 18;
    int N = 22;
    int A_p = 1;
    int A_target = Z + N;
    int A_compound = A_target + A_p;
    int chi = 1;
    double E_comp = 25.60156;
    double V = 38.0;

    double R_nu_nu = 1.5;
    double R_nu_pi = 1.0;
    double R_pi_pi = 1.0;
    double R_pi_nu = 1.0;
    double C1 = 1.0;
    double C2 = 1.0;
    double C3 = 1.0;
    double M2c = 1.0;

    double delta = pairingEnergyDelta(A_target, chi);
    double ncrit = nCrit(Z, N, delta);

    std::cout << "# header: \n";
    std::cout << "#   title: Ar40(n,x) two-component exciton model\n";
    std::cout << "#   source: exciton (C++), "
              << (mode == PreeqMode::Analytical ? "analytical" : "numerical")
              << " mode, "
              << (method == IntegrationMethod::Midpoint ? "midpoint" : "clenshaw-curtis")
              << (method == IntegrationMethod::ClenshawCurtis && guardBounds ? ", guarded" : "")
              << " integration, "
              << (kernel == CollisionKernel::OpticalModel ? "optical-model" : "matrix-element")
              << " kernel\n";
    std::cout << "#   user: Erin Visser\n";
    std::cout << "#   date: 2026-06-16\n";
    std::cout << "#   format: YANDF-0.4\n";
    std::cout << "# target: \n";
    std::cout << "#   Z: " << Z << "\n";
    std::cout << "#   A: 40\n";
    std::cout << "#   nuclide: Ar40\n";
    std::cout << "# parameters: \n";
    std::cout << "#   E-incident [MeV]:  2.000000E+01\n";
    std::cout << "#   E-compound [MeV]:  " << E_comp << "\n";
    std::cout << "#   Fermi well depth [MeV]:  " << V << "\n";
    if (kernel == CollisionKernel::MatrixElement)
    {
        std::cout << "#   Constant for matrix element:  1.000000E+00\n";
        std::cout << "#   p-p ratio for matrix element:  1.000000E+00\n";
        std::cout << "#   n-n ratio for matrix element:  1.500000E+00\n";
        std::cout << "#   p-n ratio for matrix element:  1.000000E+00\n";
        std::cout << "#   n-p ratio for matrix element:  1.000000E+00\n";
    }
    else
    {
        double Rpinu = R_pi_nu / R_pi_pi;
        double denom = 1.0 + 2.0 * Rpinu;
        double ws = M2c * 0.55 / denom;
        double wc = M2c * 0.55 * 2.0 * Rpinu / denom;
        double w0 = 0.5 * (ws + wc);
        std::cout << "#   OMP constant C_omp:  0.55\n";
        std::cout << "#   M2 constant (C1):    " << C1 << "\n";
        std::cout << "#   Rpinu:               " << Rpinu << "\n";
        std::cout << "#   Wompfac(same-type):  " << ws << "\n";
        std::cout << "#   Wompfac(cross-type): " << wc << "\n";
        std::cout << "#   Wompfac(conversion): " << w0 << "\n";
    }
    std::cout << "#   quantity: \n";
    std::cout << "#     type: internal transition rate\n";
    std::cout << "#   datablock: \n";
    std::cout << "#     columns: 8\n";

    std::vector<ExcitonState> states;

    for (int k = 0; k <= 5; ++k)
    {
        for (int i = 0; i <= k; ++i)
        {
            states.push_back({i, i, k + 1 - i, k - i});
        }
    }

    std::cout << "#     entries: " << states.size() << "\n";
    std::cout << "##     p(p)           h(p)           p(n)           h(n)";
    std::cout << "       lambdapiplus   lambdanuplus    lambdapinu     lambdanupi\n";
    std::cout << "##      []             []             []             []";
    std::cout << "          [sec^-1]       [sec^-1]       [sec^-1]       [sec^-1]\n";

    std::cout << std::scientific << std::uppercase;
    std::cout.precision(6);

    for (auto &s : states)
    {
        int n = s.p_pi + s.h_pi + s.p_nu + s.h_nu;
        double P_ph = fuPairingCorrection(delta, E_comp, n, ncrit);
        double U = availableExcitationEnergy(E_comp, P_ph);

        double lam_pi_plus = lambdaRate(LambdaType::ProtonPairCreation,
                                        mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                        R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2c,
                                        method, midBins, guardBounds, kernel, /* Z_proj = */ 0);
        double lam_nu_plus = lambdaRate(LambdaType::NeutronPairCreation,
                                        mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                        R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2c,
                                        method, midBins, guardBounds, kernel, /* Z_proj = */ 0);
        double lam_pi_nu = lambdaRate(LambdaType::ProtonToNeutronConversion,
                                      mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                      R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2c,
                                      method, midBins, guardBounds, kernel, /* Z_proj = */ 0);
        double lam_nu_pi = lambdaRate(LambdaType::NeutronToProtonConversion,
                                      mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                      R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, M2c,
                                      method, midBins, guardBounds, kernel, /* Z_proj = */ 0);

        std::cout << std::setw(14) << s.p_pi
                  << std::setw(14) << s.h_pi
                  << std::setw(14) << s.p_nu
                  << std::setw(14) << s.h_nu
                  << std::setw(15) << lam_pi_plus
                  << std::setw(15) << lam_nu_plus
                  << std::setw(15) << lam_pi_nu
                  << std::setw(15) << lam_nu_pi
                  << "\n";
    }

    return 0;
}
