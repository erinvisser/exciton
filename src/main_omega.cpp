#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "preeq_common.hh"
#include "state_density.hh"

int main()
{
    int Z = 18;
    int N = 23;
    int A = Z + N;
    int chi = 1;
    double E_p = 20.0;

    double delta = pairingEnergyDelta(A, chi);
    double ncrit = nCrit(Z, N, delta);

    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);

    std::vector<ExcitonState> states = {
        {1, 1, 0, 0},
        {0, 0, 1, 1},
        {1, 1, 1, 0},
        {1, 0, 1, 1},
        {2, 1, 0, 0},
        {0, 0, 2, 1},
        {2, 2, 0, 0},
        {0, 0, 2, 2},
        {1, 1, 1, 1},
    };

    std::vector<std::string> labels = {
        "1p1h0p0h",
        "0p0h1p1h",
        "1p1h1p0h",
        "1p0h1p1h",
        "2p1h0p0h",
        "0p0h2p1h",
        "2p2h0p0h",
        "0p0h2p2h",
        "1p1h1p1h",
    };

    std::cout << "# header: \n";
    std::cout << "#   title: Ar40 particle-hole state density\n";
    std::cout << "#   source: main_omega.cpp\n";
    std::cout << "#   user: Erin Visser\n";
    std::cout << "#   date: 2026-06-15\n";
    std::cout << "#   format: YANDF-0.4\n";
    std::cout << "# residual: \n";
    std::cout << "#   Z: " << Z << "\n";
    std::cout << "#   A: " << A << "\n";
    std::cout << "#   nuclide: Ar41\n";
    std::cout << "#   quantity: \n";
    std::cout << "#     type: particle-hole state density\n";
    std::cout << "#   datablock: \n";
    std::cout << "#     columns: 13\n";
    std::cout << "#     entries: 25\n";

    std::cout << "##      Ex            P(3)           g(p)           g(n)          ";
    for (auto &l : labels)
        std::cout << std::setw(14) << l;
    std::cout << "\n";

    std::cout << "##     [MeV]          [MeV]        [MeV^-1]       [MeV^-1]       ";
    for (size_t i = 0; i < labels.size(); ++i)
        std::cout << std::setw(14) << "[MeV^-1]";
    std::cout << "\n";

    std::cout << std::scientific << std::uppercase;
    std::cout.precision(6);

    for (int i = 1; i <= 25; ++i)
    {
        double E_x = static_cast<double>(i);

        double P3 = fuPairingCorrection(delta, E_x, 3, ncrit);
        std::cout << std::setw(14) << E_x;
        std::cout << std::setw(14) << P3;
        std::cout << std::setw(14) << g_p;
        std::cout << std::setw(14) << g_n;

        for (auto &s : states)
        {
            int h = s.h_pi + s.h_nu;
            // double V = potentialWellDepth(ProjectileType::Neutron, A, E_p, h);
            double V = 38.0; // TODO: FIXING (ignoring surface effects) for validation with TALYS.

            double omega = particleHoleStateDensity(
                s.p_pi, s.h_pi, s.p_nu, s.h_nu,
                E_x, Z, N, chi, V, ncrit, false);
            std::cout << std::setw(14) << omega;
        }
        std::cout << "\n";
    }

    // Spin distribution table
    // Parity: each J state splits evenly into π=+1 and π=-1.
    // To obtain Jπ-dependent ω: ω(p,h,E,J,π) = ω(p,h,E) × R_n(J) × 1/2
    int max_exc = 12;
    int max_J = 8;
    int max_J_sum = 30;
    std::cout << "#   quantity: \n";
    std::cout << "#     type: spin distribution\n";
    std::cout << "#     note: parity is 50/50 (factor 1/2 per π state)\n";
    std::cout << "#   datablock: \n";
    std::cout << "#     columns: 11\n";
    std::cout << "#     entries: " << max_exc << "\n";
    std::cout << "##       n";
    for (int J = 0; J <= max_J; ++J)
        std::cout << "         J= " << std::left << J << std::right;
    std::cout << "            Sum\n";
    std::cout << "##      []";
    for (int J = 0; J <= max_J; ++J)
        std::cout << "             []";
    std::cout << "             []\n";

    int A_target = A - 1; // TALYS uses target mass in spin cutoff factor
    for (int n = 1; n <= max_exc; ++n)
    {
        std::cout << "        " << std::setw(4) << n << "  ";
        for (int J = 0; J <= max_J; ++J)
        {
            double value = (2.0 * J + 1.0) * spinDistribution(static_cast<double>(J), A_target, n);
            std::cout << std::setw(15) << value;
        }
        double sum = 0.0;
        for (int J = 0; J <= max_J_sum; ++J)
        {
            double value = (2.0 * J + 1.0) * spinDistribution(static_cast<double>(J), A_target, n);
            sum += value;
        }
        std::cout << std::setw(15) << sum << "\n";
    }

    return 0;
}
