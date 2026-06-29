#ifndef TRANSITION_RATES_HH
#define TRANSITION_RATES_HH

#include "preeq_common.hh"

enum class CollisionKernel
{
    MatrixElement,
    OpticalModel
};

enum class ExcitonType
{
    Proton,
    Neutron
};

enum class ConversionType
{
    ProtonToNeutron,
    NeutronToProton
};

enum class MatrixElement
{
    PiPi,
    NuNu,
    PiNu,
    NuPi
};

enum class PreeqMode
{
    Analytical,
    Numerical
};

enum class IntegrationMethod
{
    ClenshawCurtis,
    Midpoint
};

enum class LambdaType
{
    ProtonPairCreation,
    NeutronPairCreation,
    ProtonToNeutronConversion,
    NeutronToProtonConversion
};

double M2(int A_p, int A, double E_tot, int n, double C1, double C2, double C3, PreeqMode mode);

double M2_element(MatrixElement type, int A_p, int A, double E_tot, int n,
                  double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                  double C1 = 1.0, double C2 = 1.0, double C3 = 1.0, PreeqMode mode = PreeqMode::Analytical);

double B(int p_pi, int h_pi, int p_nu, int h_nu, double g_p, double g_n, ConversionType type);

// Analytical-mode helpers
double lambdaNewPairAnalytical(ExcitonType particle, int Z, int N, int A_p,
                               int p_pi, int h_pi, int p_nu, int h_nu,
                               double E_tot, double U, double V,
                               double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                               double C1 = 1.0, double C2 = 1.0, double C3 = 1.0);

double lambdaPairConversionAnalytical(ConversionType conversion, int Z, int N, int A_p,
                                      int p_pi, int h_pi, int p_nu, int h_nu,
                                      double E_tot, double U, double V,
                                      double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                                      double C1 = 1.0, double C2 = 1.0, double C3 = 1.0);

// Numerical-mode helpers
double lambdaNewPairNumerical(ExcitonType particle, int Z, int N, int A_p,
                              int p_pi, int h_pi, int p_nu, int h_nu,
                              double E_tot, double U, double V,
                              double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                              double C1 = 1.0, double C2 = 1.0, double C3 = 1.0,
                              IntegrationMethod method = IntegrationMethod::ClenshawCurtis,
                              int midBins = 20,
                              bool guardBounds = true);

double lambdaPairConversionNumerical(ConversionType conversion, int Z, int N, int A_p,
                                     int p_pi, int h_pi, int p_nu, int h_nu,
                                     double E_tot, double U, double V,
                                     double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                                     double C1 = 1.0, double C2 = 1.0, double C3 = 1.0,
                                     IntegrationMethod method = IntegrationMethod::ClenshawCurtis,
                                     int midBins = 20,
                                     bool guardBounds = true);

double lambdaRate(LambdaType type, PreeqMode mode, int Z, int N, int A_p,
                  int p_pi, int h_pi, int p_nu, int h_nu,
                  double E_tot, double U, double V,
                  double R_nu_nu = 1.5, double R_nu_pi = 1.0, double R_pi_pi = 1.0, double R_pi_nu = 1.0,
                  double C1 = 1.0, double C2 = 1.0, double C3 = 1.0,
                  IntegrationMethod method = IntegrationMethod::ClenshawCurtis,
                  int midBins = 20,
                  bool guardBounds = true,
                  CollisionKernel kernel = CollisionKernel::MatrixElement,
                  int Z_proj = 0);

void set_diagnostic_output(bool on);

#endif
