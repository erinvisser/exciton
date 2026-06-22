#ifndef PREEQ_COMMON_HH
#define PREEQ_COMMON_HH

constexpr double PI = 3.14159265358979323846;
constexpr double HBAR = 6.582119569e-22; // MeV*s

struct ExcitonState
{
    int p_pi, h_pi, p_nu, h_nu;
};

enum class ProjectileType
{
    Proton,
    Neutron,
    Neutrino
};

double binomial(int h, int i);

double spDensityProton(int Z);

double spDensityNeutron(int N);

double potentialWellDepth(ProjectileType incident, int A, double E_p, int h);

double pauliCorrection(
    int p_pi, int h_pi,
    int p_nu, int h_nu,
    double g_p, double g_n);

double finiteWell(
    int p, int h,
    double U, double V);

double pairingEnergyDelta(int A, int chi);

double nCrit(int Z, int N, double delta);

double fuPairingCorrection(double delta, double E_x, int n, double n_crit);

double availableExcitationEnergy(double E_x, double P_ph);

double spinDistribution(double J, int A, int n, bool talysCompat = true);

#endif