#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <ctime>

#include "preeq_common.hh"
#include "transition_rates.hh"
#include "marley/KoningDelarocheOpticalModel.hh"
#include "marley/MassTable.hh"

enum class ValidationMode { KD03, WellDepths };

int main(int argc, char *argv[])
{
    PreeqMode mode = PreeqMode::Numerical;
    IntegrationMethod method = IntegrationMethod::ClenshawCurtis;
    int midBins = 50;
    bool guardBounds = true;
    CollisionKernel kernel = CollisionKernel::MatrixElement;
    bool quiet = false;
    ValidationMode vmode = ValidationMode::KD03;

    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--analytical") == 0)
            mode = PreeqMode::Analytical;
        else if (std::strcmp(argv[i], "--no-guards") == 0)
            guardBounds = false;
        else if (std::strcmp(argv[i], "--quiet") == 0)
            quiet = true;
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
        else if (std::strcmp(argv[i], "--mode") == 0 && i + 1 < argc)
        {
            ++i;
            if (std::strcmp(argv[i], "well_depths") == 0)
                vmode = ValidationMode::WellDepths;
            else if (std::strcmp(argv[i], "kd03") == 0)
                vmode = ValidationMode::KD03;
        }
    }

    if (kernel == CollisionKernel::OpticalModel && mode == PreeqMode::Analytical)
    {
        std::cerr << "# warning: OMP kernel requires numerical mode; overriding --analytical\n";
        mode = PreeqMode::Numerical;
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

    // Output file setup
    std::string outdir = "lambda_printouts/OMP";
    std::string kd03_dir = outdir + "/KD03_Output";
    std::string rates_dir = outdir + "/Rates_Output";
    std::string wd_dir = outdir + "/Well_Depths_Output";

    std::time_t now = std::time(nullptr);
    char ts[32];
    std::strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", std::localtime(&now));

    std::string mode_str = (mode == PreeqMode::Analytical) ? "analytical" : "numerical";
    std::string method_str = (method == IntegrationMethod::Midpoint) ? "midpoint" : "clenshaw-curtis";
    std::string guard_str = guardBounds ? "_guarded" : "";
    std::string kernel_str = (kernel == CollisionKernel::OpticalModel) ? "optical-model" : "matrix-element";

    std::string kd03_name = kd03_dir + "/kd03_Ar40_n_" + ts + ".log";
    std::string wd_name = wd_dir + "/well_depths_Ar40_n_" + ts + ".log";

    std::ostringstream rates_name;
    rates_name << rates_dir << "/lambda_Ar40_n_" << mode_str << "_" << method_str
               << guard_str << "_" << kernel_str << "_" << ts << ".log";

    std::ofstream kd03_file;
    std::ofstream rates_file;
    std::ofstream wd_file;

    if (vmode == ValidationMode::KD03)
    {
        kd03_file.open(kd03_name);
        rates_file.open(rates_name.str());
        if (!kd03_file.is_open())
        {
            std::cerr << "Error: could not open " << kd03_name << "\n";
            return 1;
        }
        if (!rates_file.is_open())
        {
            std::cerr << "Error: could not open " << rates_name.str() << "\n";
            return 1;
        }
        if (!quiet)
        {
            std::cout << "# wrote KD03 diagnostics to " << kd03_name << "\n";
            std::cout << "# wrote rates to " << rates_name.str() << "\n";
        }
    }
    else if (vmode == ValidationMode::WellDepths)
    {
        wd_file.open(wd_name);
        if (!wd_file.is_open())
        {
            std::cerr << "Error: could not open " << wd_name << "\n";
            return 1;
        }
        if (!quiet)
            std::cout << "# wrote well depths to " << wd_name << "\n";
    }

    // KD03 diagnostic block
    if (vmode == ValidationMode::KD03)
    {
        marley::KoningDelarocheOpticalModel kd(18, 40);
        double A_t = 40.0;
        double A13 = std::pow(A_t, 1.0 / 3.0);

        kd03_file << "# KD03 optical model parameters\n";
        kd03_file << "# Target: Z=18  A=40  nuclide=Ar40\n";
        kd03_file << "# Projectile: n\n";
        kd03_file << "# MARLEY KoningDelarocheOpticalModel\n";
        kd03_file << "#\n";

        // Neutron coefficients
        kd03_file << "# --- Neutron (k=1) coefficients ---\n";
        kd03_file << "# Real volume:  v1n=" << kd.getV1n() << "  v2n=" << kd.getV2n()
                  << "  v3n=" << kd.getV3n() << "  v4n=" << kd.getV4n() << "\n";
        kd03_file << "# Imag volume:  w1n=" << kd.getW1n() << "  w2n=" << kd.getW2n() << "\n";
        kd03_file << "# Imag surf:    d1n=" << kd.getD1n() << "  d2n=" << kd.getD2n()
                  << "  d3n=" << kd.getD3n() << "\n";
        kd03_file << "# Spin-orbit:   vso1n=" << kd.getVso1n() << "  vso2n=" << kd.getVso2n()
                  << "  wso1n=" << kd.getWso1n() << "  wso2n=" << kd.getWso2n() << "\n";
        kd03_file << "# Fermi energy: Efn=" << kd.getEfn() << "\n";
        kd03_file << "# Geometry:\n";
        kd03_file << "#   Rvn=" << kd.getRvn() << "  avn=" << kd.getavn() << "\n";
        kd03_file << "#   Rdn=" << kd.getRdn() << "  adn=" << kd.getadn() << "\n";
        kd03_file << "#   Rso_n=" << kd.getRso_n() << "  aso_n=" << kd.getaso_n() << "\n";
        kd03_file << "# Reduced radii:  rv=" << (kd.getRvn() / A13)
                  << "  rd=" << (kd.getRdn() / A13)
                  << "  rso=" << (kd.getRso_n() / A13) << "\n";
        kd03_file << "#\n";

        // Proton coefficients
        kd03_file << "# --- Proton (k=2) coefficients ---\n";
        kd03_file << "# Real volume:  v1p=" << kd.getV1p() << "  v2p=" << kd.getV2p()
                  << "  v3p=" << kd.getV3p() << "  v4p=" << kd.getV4p() << "\n";
        kd03_file << "# Imag volume:  w1p=" << kd.getW1p() << "  w2p=" << kd.getW2p() << "\n";
        kd03_file << "# Imag surf:    d1p=" << kd.getD1p() << "  d2p=" << kd.getD2p()
                  << "  d3p=" << kd.getD3p() << "\n";
        kd03_file << "# Spin-orbit:   vso1p=" << kd.getVso1p() << "  vso2p=" << kd.getVso2p()
                  << "  wso1p=" << kd.getWso1p() << "  wso2p=" << kd.getWso2p() << "\n";
        kd03_file << "# Fermi energy: Efp=" << kd.getEfp() << "\n";
        kd03_file << "# Coulomb:      Vcbar_p=" << kd.getVcbar_p() << "  Rc=" << kd.getRc() << "\n";
        kd03_file << "# Geometry:\n";
        kd03_file << "#   Rvp=" << kd.getRvp() << "  avp=" << kd.getavp() << "\n";
        kd03_file << "#   Rdp=" << kd.getRdp() << "  adp=" << kd.getadp() << "\n";
        kd03_file << "#   Rso_p=" << kd.getRso_p() << "  aso_p=" << kd.getaso_p() << "\n";
        kd03_file << "#\n";

        // Q-values from MARLEY MassTable
        {
            auto& mt = marley::MassTable::Instance();
            double u = 931.494061; // atomic mass unit (MeV)
            double ME_n = mt.get_particle_mass(marley_utils::NEUTRON) - u;
            double ME_Ar40 = mt.get_mass_excess(18, 40);
            double ME_Ar41 = mt.get_mass_excess(18, 41);
            double ME_K40  = mt.get_mass_excess(19, 40);
            double ME_Cl39 = mt.get_mass_excess(17, 39);
            double ME_Cl38 = mt.get_mass_excess(17, 38);
            double ME_S38  = mt.get_mass_excess(16, 38);
            double ME_S37  = mt.get_mass_excess(16, 37);
            double ME_H1   = mt.get_mass_excess(1, 1);
            double ME_d    = mt.get_mass_excess(1, 2);
            double ME_t    = mt.get_mass_excess(1, 3);
            double ME_h    = mt.get_mass_excess(2, 3);
            double ME_a    = mt.get_mass_excess(2, 4);

            kd03_file << "# --- Q-values [MeV] from MARLEY MassTable ---\n";
            kd03_file << "# (n,gamma):   " << (ME_Ar40 + ME_n - ME_Ar41) << "\n";
            kd03_file << "# (n,n):       0.0\n";
            kd03_file << "# (n,p):       " << (ME_Ar40 + ME_n - ME_H1 - ME_K40) << "\n";
            kd03_file << "# (n,d):       " << (ME_Ar40 + ME_n - ME_d - ME_Cl39) << "\n";
            kd03_file << "# (n,t):       " << (ME_Ar40 + ME_n - ME_t - ME_Cl38) << "\n";
            kd03_file << "# (n,3He):     " << (ME_Ar40 + ME_n - ME_h - ME_S38) << "\n";
            kd03_file << "# (n,alpha):   " << (ME_Ar40 + ME_n - ME_a - ME_S37) << "\n";
            kd03_file << "#\n";
        }

        // Energy-dependent OMP values (TALYS flagomponly format)
        kd03_file << "# --- Energy-dependent OMP values (TALYS flagomponly format) ---\n";
        kd03_file << "# Grid: -8.0 to 20.0 MeV, 0.1 MeV step\n";
        for (int nen = -80; nen <= 200; ++nen)
        {
            double e = 0.1 * nen;
            bool print = (nen == -80 || nen == -70 || nen == -60 || nen == 0
                        || nen == 10 || nen == 20 || (nen > 0 && nen % 20 == 0));
            if (!print) continue;
            for (int k = 1; k <= 2; ++k)
            {
                int pdg = (k == 2) ? 2212 : 2112;
                double e_safe = std::max(e, 0.0);
                kd.setIncidentEnergyAndFragment(e_safe, pdg);

                double Vv = kd.getVv();
                double Wv = kd.getWv();
                double Wd = kd.getWd();
                double Rv = kd.getRv();
                double av = kd.getav();
                double Rd = kd.getRd();
                double ad = kd.getad();
                double rv = Rv / A13;
                double rd = Rd / A13;

                std::string frag_type = (k == 2) ? "proton" : "neutron";
                double Ef = (k == 2) ? kd.getEfp() : kd.getEfn();

                double v1 = (k == 2) ? kd.getV1p() : kd.getV1n();
                double v2 = (k == 2) ? kd.getV2p() : kd.getV2n();
                double v3 = (k == 2) ? kd.getV3p() : kd.getV3n();
                double v4 = (k == 2) ? kd.getV4p() : kd.getV4n();
                double Vcoul = (k == 2) ? kd.getVcbar_p() : 0.0;
                double w1 = (k == 2) ? kd.getW1p() : kd.getW1n();
                double w2 = (k == 2) ? kd.getW2p() : kd.getW2n();
                double d1 = (k == 2) ? kd.getD1p() : kd.getD1n();
                double d2 = (k == 2) ? kd.getD2p() : kd.getD2n();
                double d3 = (k == 2) ? kd.getD3p() : kd.getD3n();
                double vso1 = (k == 2) ? kd.getVso1p() : kd.getVso1n();
                double vso2 = (k == 2) ? kd.getVso2p() : kd.getVso2n();
                double wso1 = (k == 2) ? kd.getWso1p() : kd.getWso1n();
                double wso2 = (k == 2) ? kd.getWso2p() : kd.getWso2n();

                double rvso = (k == 2) ? kd.getRso_p() / A13 : kd.getRso_n() / A13;
                double avso = (k == 2) ? kd.getaso_p() : kd.getaso_n();

                kd03_file << "KD03 OMP parameters for " << frag_type
                          << "  E:" << std::setw(11) << std::fixed << std::setprecision(5) << e
                          << " Ef:   " << std::setw(9) << std::fixed << std::setprecision(5) << Ef << "\n";
                kd03_file << "   rv:" << std::setw(10) << std::fixed << std::setprecision(5) << rv
                          << "   av:" << std::setw(10) << std::fixed << std::setprecision(5) << av
                          << "   v1:" << std::setw(10) << std::fixed << std::setprecision(5) << v1
                          << "   v2:" << std::setw(10) << std::fixed << std::setprecision(5) << v2
                          << "   v3:" << std::scientific << std::setprecision(5) << std::setw(12) << v3
                          << "   v4:" << std::scientific << std::setprecision(5) << std::setw(12) << v4
                          << " Vcoul:" << std::fixed << std::setprecision(5) << std::setw(10) << Vcoul << "\n";
                kd03_file << "   rw:" << std::setw(10) << std::fixed << std::setprecision(5) << rv
                          << "   aw:" << std::setw(10) << std::fixed << std::setprecision(5) << av
                          << "   w1:" << std::setw(10) << std::fixed << std::setprecision(5) << w1
                          << "   w2:" << std::setw(10) << std::fixed << std::setprecision(5) << w2 << "\n";
                kd03_file << "  rwd:" << std::setw(9) << std::fixed << std::setprecision(5) << rd
                          << "  awd:" << std::setw(8) << std::fixed << std::setprecision(5) << ad
                          << "   d1:" << std::setw(10) << std::fixed << std::setprecision(5) << d1
                          << "   d2:" << std::setw(10) << std::fixed << std::setprecision(5) << d2
                          << "   d3:" << std::setw(10) << std::fixed << std::setprecision(5) << d3 << "\n";
                kd03_file << " rvso:" << std::setw(9) << std::fixed << std::setprecision(5) << rvso
                          << " avso:" << std::setw(8) << std::fixed << std::setprecision(5) << avso
                          << " vso1:" << std::setw(10) << std::fixed << std::setprecision(5) << vso1
                          << " vso2:" << std::setw(10) << std::fixed << std::setprecision(5) << vso2 << "\n";
                kd03_file << " rwso:" << std::setw(9) << std::fixed << std::setprecision(5) << rvso
                          << " awso:" << std::setw(8) << std::fixed << std::setprecision(5) << avso
                          << " wso1:" << std::setw(10) << std::fixed << std::setprecision(5) << wso1
                          << " wso2:" << std::setw(10) << std::fixed << std::setprecision(5) << wso2 << "\n";
            }
        }
    }

    // Well-depths diagnostic block
    if (vmode == ValidationMode::WellDepths)
    {
        marley::KoningDelarocheOpticalModel kd(18, 40);
        double A_t = 40.0;
        double A13 = std::pow(A_t, 1.0 / 3.0);

        wd_file << "# Well depths from MARLEY KoningDelarocheOpticalModel\n";
        wd_file << "# Target: Z=18  A=40  nuclide=Ar40\n";
        wd_file << "# Grid: -8.0 to 20.0 MeV, 0.1 MeV step\n";
        wd_file << "#\n";
        wd_file << "# --- Table 1: Energy-dependent potential depths ---\n";
        wd_file << "# Columns: E(MeV)  k  Vv(MeV)  Wv(MeV)  Wd(MeV)  Vso(MeV)  Wso(MeV)\n";

        for (int nen = -80; nen <= 200; ++nen)
        {
            double e = 0.1 * nen;
            bool print = (nen == -80 || nen == -70 || nen == -60 || nen == 0
                        || nen == 10 || nen == 20 || (nen > 0 && nen % 20 == 0));
            if (!print) continue;
            for (int k = 1; k <= 2; ++k)
            {
                int pdg = (k == 2) ? 2212 : 2112;
                double e_safe = std::max(e, 0.0);
                kd.setIncidentEnergyAndFragment(e_safe, pdg);

                double Vv = kd.getVv();
                double Wv = kd.getWv();
                double Wd = kd.getWd();
                double Vso = kd.getVso();
                double Wso = kd.getWso();

                wd_file << std::setw(10) << std::fixed << std::setprecision(5) << e
                        << std::setw(4) << k
                        << std::setw(14) << std::scientific << std::setprecision(6) << Vv
                        << std::setw(14) << Wv
                        << std::setw(14) << Wd
                        << std::setw(14) << Vso
                        << std::setw(14) << Wso << "\n";
            }
        }

        wd_file << "#\n";
        wd_file << "# --- Table 2: Folded geometries ---\n";
        wd_file << "# Columns: E(MeV)  k  Rv(fm)  av(fm)  Rd(fm)  ad(fm)  Rso(fm)  aso(fm)\n";

        for (int nen = -80; nen <= 200; ++nen)
        {
            double e = 0.1 * nen;
            bool print = (nen == -80 || nen == -70 || nen == -60 || nen == 0
                        || nen == 10 || nen == 20 || (nen > 0 && nen % 20 == 0));
            if (!print) continue;
            for (int k = 1; k <= 2; ++k)
            {
                int pdg = (k == 2) ? 2212 : 2112;
                double e_safe = std::max(e, 0.0);
                kd.setIncidentEnergyAndFragment(e_safe, pdg);

                double Rv = kd.getRv();
                double av = kd.getav();
                double Rd = kd.getRd();
                double ad = kd.getad();
                double Rso = (k == 2) ? kd.getRso_p() : kd.getRso_n();
                double aso = (k == 2) ? kd.getaso_p() : kd.getaso_n();

                wd_file << std::setw(10) << std::fixed << std::setprecision(5) << e
                        << std::setw(4) << k
                        << std::setw(14) << std::fixed << std::setprecision(6) << Rv
                        << std::setw(14) << av
                        << std::setw(14) << Rd
                        << std::setw(14) << ad
                        << std::setw(14) << Rso
                        << std::setw(14) << aso << "\n";
            }
        }
    }

    // Rates computation
    if (vmode == ValidationMode::KD03)
    {
        rates_file << "# header: \n";
        rates_file << "#   title: Ar40(n,x) two-component exciton model\n";
        rates_file << "#   source: exciton (C++), "
                  << mode_str << " mode, "
                  << method_str
                  << (method == IntegrationMethod::ClenshawCurtis && guardBounds ? ", guarded" : "")
                  << " integration, "
                  << kernel_str << " kernel\n";
        rates_file << "#   user: Erin Visser\n";
        rates_file << "#   date: 2026-06-23\n";
        rates_file << "#   format: YANDF-0.4\n";
        rates_file << "# target: \n";
        rates_file << "#   Z: " << Z << "\n";
        rates_file << "#   A: 40\n";
        rates_file << "#   nuclide: Ar40\n";
        rates_file << "# parameters: \n";
        rates_file << "#   E-incident [MeV]:  2.000000E+01\n";
        rates_file << "#   E-compound [MeV]:  " << E_comp << "\n";
        rates_file << "#   Fermi well depth [MeV]:  " << V << "\n";
        if (kernel == CollisionKernel::MatrixElement)
        {
            rates_file << "#   Constant for matrix element:  1.000000E+00\n";
            rates_file << "#   p-p ratio for matrix element:  1.000000E+00\n";
            rates_file << "#   n-n ratio for matrix element:  1.500000E+00\n";
            rates_file << "#   p-n ratio for matrix element:  1.000000E+00\n";
            rates_file << "#   n-p ratio for matrix element:  1.000000E+00\n";
        }
        else
        {
            double Rpinu = R_pi_nu / R_pi_pi;
            double denom = 1.0 + 2.0 * Rpinu;
            double ws = M2c * 0.55 / denom;
            double wc = M2c * 0.55 * 2.0 * Rpinu / denom;
            double w0 = 0.5 * (ws + wc);
            rates_file << "#   OMP constant C_omp:  0.55\n";
            rates_file << "#   M2 constant (C1):    " << C1 << "\n";
            rates_file << "#   Rpinu:               " << Rpinu << "\n";
            rates_file << "#   Wompfac(same-type):  " << ws << "\n";
            rates_file << "#   Wompfac(cross-type): " << wc << "\n";
            rates_file << "#   Wompfac(conversion): " << w0 << "\n";
        }
        rates_file << "#   quantity: \n";
        rates_file << "#     type: internal transition rate\n";
        rates_file << "#   datablock: \n";
        rates_file << "#     columns: 8\n";

        std::vector<ExcitonState> states;

        for (int k = 0; k <= 5; ++k)
        {
            for (int i = 0; i <= k; ++i)
            {
                states.push_back({i, i, k + 1 - i, k - i});
            }
        }

        rates_file << "#     entries: " << states.size() << "\n";
        rates_file << "##     p(p)           h(p)           p(n)           h(n)";
        rates_file << "       lambdapiplus   lambdanuplus    lambdapinu     lambdanupi\n";
        rates_file << "##      []             []             []             []";
        rates_file << "          [sec^-1]       [sec^-1]       [sec^-1]       [sec^-1]\n";

        rates_file << std::scientific << std::uppercase;
        rates_file.precision(6);

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

            rates_file << std::setw(14) << s.p_pi
                      << std::setw(14) << s.h_pi
                      << std::setw(14) << s.p_nu
                      << std::setw(14) << s.h_nu
                      << std::setw(15) << lam_pi_plus
                      << std::setw(15) << lam_nu_plus
                      << std::setw(15) << lam_pi_nu
                      << std::setw(15) << lam_nu_pi
                      << "\n";
        }
    }

    return 0;
}
