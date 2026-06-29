#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <ctime>
#include <sys/stat.h>

#include "preeq_common.hh"
#include "transition_rates.hh"

int main(int argc, char *argv[])
{
    PreeqMode mode = PreeqMode::Numerical;
    IntegrationMethod method = IntegrationMethod::ClenshawCurtis;
    CollisionKernel kernel = CollisionKernel::MatrixElement;
    int midBins = 50;
    bool guardBounds = true;
    int Z_proj = 0;

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
        }
        else if (std::strcmp(argv[i], "--diagnostic") == 0)
        {
            set_diagnostic_output(true);
        }
    }

    // Build output filename with timestamp
    std::time_t now = std::time(nullptr);
    char time_buf[32];
    std::strftime(time_buf, sizeof(time_buf), "%Y%m%d_%H%M%S", std::localtime(&now));

    std::string kernel_str = (kernel == CollisionKernel::OpticalModel)
      ? "optical-model" : "matrix-element";
    std::string method_str = (method == IntegrationMethod::Midpoint)
      ? "midpoint" : "clenshaw-curtis";

    std::string outdir = "/Users/erinvisser/Downloads/opencode/exciton/lambda_printouts/OMP/TALYS_Mimic";
    mkdir(outdir.c_str(), 0755); // will fail silently if exists

    std::string outpath = outdir + "/lambda_" + kernel_str + "_" + method_str + "_" + time_buf + ".txt";
    std::ofstream outfile(outpath);
    std::streambuf* old_cout_buf = nullptr;
    if (outfile.is_open()) {
      old_cout_buf = std::cout.rdbuf(outfile.rdbuf());
    }
    std::cerr << "# Output: " << outpath << "\n";

    int Z = 18;
    int N = 23; // compound N (target is 22, projectile adds 1)
    int A_p = 1;
    int A_target = 40; // target mass number for pairing energy (Ar40)
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

    double delta = pairingEnergyDelta(Z + N, chi); // compound mass, matches TALYS pair() call
    double ncrit = nCrit(Z, N, delta);

    std::cout << "# header: \n";
    std::cout << "#   title: Ar40(n,x) two-component exciton model\n";
    std::cout << "#   source: exciton (C++), "
              << (mode == PreeqMode::Analytical ? "analytical" : "numerical")
              << " mode, "
              << (method == IntegrationMethod::Midpoint ? "midpoint" : "clenshaw-curtis")
              << (method == IntegrationMethod::ClenshawCurtis && guardBounds ? ", guarded" : "")
              << (kernel == CollisionKernel::OpticalModel ? ", optical-model kernel" : ", matrix-element kernel")
              << "\n";
    std::cout << "#   user: Erin Visser\n";
    std::cout << "#   date: 2026-06-29\n";
    std::cout << "#   format: YANDF-0.4\n";
    std::cout << "# target: \n";
    std::cout << "#   Z: " << Z << "\n";
    std::cout << "#   A: 40\n";
    std::cout << "#   nuclide: Ar40\n";
    std::cout << "# parameters: \n";
    std::cout << "#   E-incident [MeV]:  2.000000E+01\n";
    std::cout << "#   E-compound [MeV]:  " << E_comp << "\n";
    std::cout << "#   Fermi well depth [MeV]:  " << V << "\n";
    std::cout << "#   Constant for matrix element:  1.000000E+00\n";
    std::cout << "#   p-p ratio for matrix element:  1.000000E+00\n";
    std::cout << "#   n-n ratio for matrix element:  1.500000E+00\n";
    std::cout << "#   p-n ratio for matrix element:  1.000000E+00\n";
    std::cout << "#   n-p ratio for matrix element:  1.000000E+00\n";
    std::cout << "#   quantity: \n";
    std::cout << "#     type: internal transition rate\n";
    std::cout << "#   datablock: \n";
    std::cout << "#     columns: 8\n";

    if (kernel == CollisionKernel::OpticalModel) {
      std::cout << "#   optical-model parameters: \n";
      std::cout << "#     M2c (C1): " << C1 << "\n";
      std::cout << "#     Rpinu (R_pi_nu): " << R_pi_nu << "\n";
      double denom = 1. + 2. * R_pi_nu;
      double wf_same = C1 * 0.55 / denom;
      double wf_cross = C1 * 0.55 * 2. * R_pi_nu / denom;
      double wf_exch = 0.5 * (wf_same + wf_cross);
      std::cout << "#     Wompfac_same (Wompfac1): " << wf_same << "\n";
      std::cout << "#     Wompfac_cross (Wompfac2): " << wf_cross << "\n";
      std::cout << "#     Wompfac_exchange (Wompfac0): " << wf_exch << "\n";
      std::cout << "#     Z_proj: " << Z_proj << "\n";
    }

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
                                        R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                        method, midBins, guardBounds,
                                        kernel, Z_proj);
        double lam_nu_plus = lambdaRate(LambdaType::NeutronPairCreation,
                                        mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                        R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                        method, midBins, guardBounds,
                                        kernel, Z_proj);
        double lam_pi_nu = lambdaRate(LambdaType::ProtonToNeutronConversion,
                                      mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                      R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                      method, midBins, guardBounds,
                                      kernel, Z_proj);
        double lam_nu_pi = lambdaRate(LambdaType::NeutronToProtonConversion,
                                      mode, Z, N, A_p, s.p_pi, s.h_pi, s.p_nu, s.h_nu, E_comp, U, V,
                                      R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                      method, midBins, guardBounds,
                                      kernel, Z_proj);

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

    if (old_cout_buf) std::cout.rdbuf(old_cout_buf);
    return 0;
}
