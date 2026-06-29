#include "transition_rates.hh"
#include "preeq_common.hh"
#include "state_density.hh"

#include "marley/Integrator.hh"
#include "marley/TransitionRates.hh"
#include "marley/MassTable.hh"
#include "marley/marley_utils.hh"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

// ---- TALYS-clone diagnostic mode ------------------------------------------

static bool diagnostic_ = false;
void set_diagnostic_output(bool on) { diagnostic_ = on; }

// ---- TALYS OMP table (wvol) infrastructure --------------------------------
// Mimics bonetti.f90: precomputes wvol(k, nen) for k=1 (neutron), k=2 (proton)
// and nen from NEN_MIN to nenend, where e = 0.1 * nen.

static const int NEN_OFFSET = 300;  // enough to shift nen >= -200 to index >= 0
static const int NUMEN_TALYS = 200;

struct WvolTable {
  double data[2][2000]; // [k-1][nen + NEN_OFFSET], k=1->idx0, k=2->idx1
  int nenend = 0;
  int Z_target = 0;
  int A_target = 0;
  double M2c = 0;
  double Rpinu = 0;
  double Wompfac_same = 0;
  double Wompfac_cross = 0;
  double Wompfac_exch = 0;
  bool built = false;

  double& wvol(int k, int nen) { return data[k-1][nen + NEN_OFFSET]; }
  double wvol(int k, int nen) const { return data[k-1][nen + NEN_OFFSET]; }
};

static WvolTable& wvol_table() {
  static WvolTable tbl;
  return tbl;
}

static void ensure_wvol_table(int Z_comp, int N_comp, int A_p, int Z_proj,
                              double M2c, double Rpinu_val)
{
  WvolTable& tbl = wvol_table();
  int Z_target = Z_comp - Z_proj;
  int N_target = N_comp - (A_p - Z_proj);
  int A_target = Z_target + N_target;

  // Rebuild if parameters changed
  if (tbl.built && tbl.Z_target == Z_target && tbl.A_target == A_target
      && tbl.M2c == M2c && tbl.Rpinu == Rpinu_val)
    return;

  tbl.Z_target = Z_target;
  tbl.A_target = A_target;
  tbl.M2c = M2c;
  tbl.Rpinu = Rpinu_val;

  // Wompfac (TALYS matrix.f90 preeqmode==3)
  double denom = 1. + 2. * Rpinu_val;
  tbl.Wompfac_same = M2c * 0.55 / denom;
  tbl.Wompfac_cross = M2c * 0.55 * 2. * Rpinu_val / denom;
  tbl.Wompfac_exch = 0.5 * (tbl.Wompfac_same + tbl.Wompfac_cross);

  // Build the OMP table using target-like TransitionRates
  // TALYS bonetti.f90: Zix = Zindex(0,0,k0), Nix = Nindex(0,0,k0) -> target
  marley::TransitionRates rates(Z_target, A_target, 0.1);
  rates.set_M2c(M2c);
  rates.set_Rpinu(Rpinu_val);

  // emax = enincmax + S(0,0,k0).  For our case ~20 + S_n(compound).
  // Use a generous bound.
  double emax = 100.;
  int nenend = 10 * std::min(NUMEN_TALYS, static_cast<int>(emax + 1.));
  tbl.nenend = nenend;

  // Build for k=1 (neutron, PDG 2112) and k=2 (proton, PDG 2212)
  for (int k = 0; k < 2; ++k) {
    int pdg = (k == 0) ? marley_utils::NEUTRON : marley_utils::PROTON;
    for (int nen = -80; nen <= nenend; ++nen) {
      double e = 0.1 * static_cast<double>(nen);
      tbl.data[k][nen + NEN_OFFSET] = rates.compute_wvol(pdg, e);
    }
  }

  tbl.built = true;
}

// ---- TALYS channel mapping for OMP (literal indexing) ---------------------
//
// TALYS conventions:
//   k = 1 -> neutron-like (PDG 2112), wvol[1][*]
//   k = 2 -> proton-like (PDG 2212), wvol[2][*]
//   S(Zix, Nix, k) -> separation energy from residual (Z_comp-Zix, N_comp-Nix)

struct TalysOmpChannel {
  int k;      // 1=neutron, 2=proton (TALYS wvol index)
  int Z_res;  // Z of residual for separation energy
  int A_res;  // A of residual for separation energy
  int pdg;    // PDG code for separation energy and wvol
};

// Maps TALYS index j (1..4) to the corresponding OMP channel
// for lambdapiplus (Proton pair creation)
static TalysOmpChannel channel_piplus(int j, int Z_comp, int N_comp) {
  if (j <= 2) {
    // proton-proton terms: same-type, k=2 (proton)
    // S(1,0,2): residual = (Z_comp-1, N_comp)
    int Z = Z_comp - 1;
    int A = Z + N_comp;
    return {2, Z, A, marley_utils::PROTON};
  } else {
    // neutron-proton terms: cross-type, k=1 (neutron)
    // S(0,1,1): residual = (Z_comp, N_comp-1)
    int Z = Z_comp;
    int A = Z + (N_comp - 1);
    return {1, Z, A, marley_utils::NEUTRON};
  }
}

// For lambdanuplus (Neutron pair creation)
static TalysOmpChannel channel_nuplus(int j, int Z_comp, int N_comp) {
  if (j <= 2) {
    // neutron-neutron terms: same-type, k=1 (neutron)
    // S(0,1,1): residual = (Z_comp, N_comp-1)
    int Z = Z_comp;
    int A = Z + (N_comp - 1);
    return {1, Z, A, marley_utils::NEUTRON};
  } else {
    // proton-neutron terms: cross-type, k=2 (proton)
    // S(1,0,2): residual = (Z_comp-1, N_comp)
    int Z = Z_comp - 1;
    int A = Z + N_comp;
    return {2, Z, A, marley_utils::PROTON};
  }
}

// For lambdapinu (ProtonToNeutron conversion)
static TalysOmpChannel channel_pinu(int Z_comp, int N_comp) {
  // Zix=1, Nix=0, k=2 -> proton from (Z_comp-1, N_comp)
  int Z = Z_comp - 1;
  int A = Z + N_comp;
  return {2, Z, A, marley_utils::PROTON};
}

// For lambdanupi (NeutronToProton conversion)
static TalysOmpChannel channel_nupi(int Z_comp, int N_comp) {
  // Zix=0, Nix=1, k=1 -> neutron from (Z_comp, N_comp-1)
  int Z = Z_comp;
  int A = Z + (N_comp - 1);
  return {1, Z, A, marley_utils::NEUTRON};
}

// Separation energy via MARLEY MassTable (TALYS convention)
static double talys_sep(const TalysOmpChannel& ch) {
  return marley::MassTable::Instance().get_fragment_separation_energy(
    ch.Z_res, ch.A_res, ch.pdg, true);
}

// ============================================================================
// TALYS-clone OMP lambda functions (literal translations)
// ============================================================================

static double lambdaNewPairOMP(ExcitonType particle, int Z, int N, int A_p,
    int p_pi, int h_pi, int p_nu, int h_nu,
    double E_tot, double U, double V,
    double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
    double C1, double C2, double C3,
    IntegrationMethod method, int midBins, bool guardBounds,
    int Z_proj)
{
  int n = p_pi + h_pi + p_nu + h_nu;
  double g_p = spDensityProton(Z);
  double g_n = spDensityNeutron(N);
  int Z_comp = Z;
  int N_comp = N;

  // Ensure wvol table is built
  ensure_wvol_table(Z_comp, N_comp, A_p, Z_proj, C1, R_pi_nu);
  const WvolTable& tbl = wvol_table();
  double same = tbl.Wompfac_same;
  double cross = tbl.Wompfac_cross;

  double omega_initial = particleHoleStateDensity(
    p_pi, h_pi, p_nu, h_nu, U, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
  if (omega_initial <= 0.0) return 0.0;

  int nexc = std::max(midBins / 2, 2);  // TALYS: primary? nbins/2 : nbins/4

  struct SubTerm {
    double L1, L2, dEx;
    int t_pi, t_hi, t_pn, t_hn;
    int r_pi, r_hi, r_pn, r_hn;
    int j;            // TALYS j index (1..4)
    double g_factor;  // gsp or gsn
  };

  std::vector<SubTerm> terms;

  if (particle == ExcitonType::Proton) {
    // ---- lambdapiplus-style: 4 subterms ----
    double new_pi = pauliCorrection(p_pi+1, h_pi+1, p_nu, h_nu, g_p, g_n);
    double less_pi_p = pauliCorrection(p_pi-1, h_pi, p_nu, h_nu, g_p, g_n);
    double less_pi_h = pauliCorrection(p_pi, h_pi-1, p_nu, h_nu, g_p, g_n);
    double less_nu_p = pauliCorrection(p_pi, h_pi, p_nu-1, h_nu, g_p, g_n);
    double less_nu_h = pauliCorrection(p_pi, h_pi, p_nu, h_nu-1, g_p, g_n);

    terms.push_back({new_pi - less_pi_p, U - less_pi_p, 0,
      2,1,0,0, p_pi-1,h_pi,p_nu,hnu, 1, g_p});
    terms.push_back({new_pi - less_pi_h, U - less_pi_h, 0,
      1,2,0,0, p_pi,h_pi-1,p_nu,hnu, 2, g_p});
    terms.push_back({new_pi - less_nu_p, U - less_nu_p, 0,
      1,1,1,0, p_pi,h_pi,p_nu-1,hnu, 3, g_n});
    terms.push_back({new_pi - less_nu_h, U - less_nu_h, 0,
      1,1,0,1, p_pi,h_pi,p_nu,hnu-1, 4, g_n});
  } else {
    // ---- lambdanuplus-style: 4 subterms ----
    double new_nu = pauliCorrection(p_pi, h_pi, p_nu+1, h_nu+1, g_p, g_n);
    double less_nu_p = pauliCorrection(p_pi, h_pi, p_nu-1, h_nu, g_p, g_n);
    double less_nu_h = pauliCorrection(p_pi, h_pi, p_nu, h_nu-1, g_p, g_n);
    double less_pi_p = pauliCorrection(p_pi-1, h_pi, p_nu, h_nu, g_p, g_n);
    double less_pi_h = pauliCorrection(p_pi, h_pi-1, p_nu, h_nu, g_p, g_n);

    terms.push_back({new_nu - less_nu_p, U - less_nu_p, 0,
      0,0,2,1, p_pi,h_pi,p_nu-1,hnu, 1, g_n});
    terms.push_back({new_nu - less_nu_h, U - less_nu_h, 0,
      0,0,1,2, p_pi,h_pi,p_nu,hnu-1, 2, g_n});
    terms.push_back({new_nu - less_pi_p, U - less_pi_p, 0,
      1,0,1,1, p_pi-1,h_pi,p_nu,hnu, 3, g_p});
    terms.push_back({new_nu - less_pi_h, U - less_pi_h, 0,
      0,1,1,1, p_pi,h_pi-1,p_nu,hnu, 4, g_p});
  }

  for (auto& t : terms) t.dEx = (t.L2 - t.L1) / nexc;

  double sum = 0.;

  for (int i = 0; i < nexc; ++i) {
    for (size_t ti = 0; ti < terms.size(); ++ti) {
      auto& t = terms[ti];
      double uu = t.L1 + (static_cast<double>(i) + 0.5) * t.dEx;
      if (uu <= 0) continue;

      // TALYS channel for this j
      TalysOmpChannel ch = (particle == ExcitonType::Proton)
        ? channel_piplus(t.j, Z_comp, N_comp)
        : channel_nuplus(t.j, Z_comp, N_comp);

      double sep = talys_sep(ch);
      double eopt = std::max(uu - sep, -20.);
      int nen = std::min(10 * NUMEN_TALYS, static_cast<int>(eopt * 10.));
      if (nen < -80) continue;

      double wv = tbl.wvol(ch.k, nen);
      double Weff = (t.j <= 2) ? same * wv : cross * wv;
      double lambda_col = 2. * Weff / HBAR;

      // Hole terms (j==2 or j==4) get density ratio
      if (t.j == 2) {
        double densh = particleHoleStateDensity(
          1,2,0,0, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
        double densp = particleHoleStateDensity(
          2,1,0,0, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
        double ratio = (densp > 1.) ? densh / densp : 1.;
        lambda_col *= ratio;
      } else if (t.j == 4) {
        double densh = particleHoleStateDensity(
          t.t_pi,t.t_hi,t.t_pn,t.t_hn, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
        double densp = particleHoleStateDensity(
          1,1,1,0, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
        double ratio = (densp > 1.) ? densh / densp : 1.;
        lambda_col *= ratio;
      }

      double omega_r = particleHoleStateDensity(
        t.r_pi, t.r_hi, t.r_pn, t.r_hn, U - uu,
        Z_comp, N_comp, 0, V, 0.0, false, 1.0);

      sum += lambda_col * t.g_factor * t.dEx * omega_r;

      // Diagnostic: print first bin of first term
      if (diagnostic_ && i == 0 && ti == 0) {
        std::cerr << "# DIAG pair_creation particle="
                  << (particle == ExcitonType::Proton ? "proton" : "neutron")
                  << " p_pi=" << p_pi << " h_pi=" << h_pi
                  << " p_nu=" << p_nu << " h_nu=" << h_nu
                  << " j=" << t.j
                  << " k=" << ch.k
                  << " Z_res=" << ch.Z_res << " A_res=" << ch.A_res
                  << " pdg=" << ch.pdg
                  << " sep=" << sep
                  << " L1=" << t.L1 << " L2=" << t.L2 << " dEx=" << t.dEx
                  << " uu=" << uu
                  << " eopt=" << eopt << " nen=" << nen
                  << " wvol=" << wv << " Weff=" << Weff
                  << " lambda_col=" << lambda_col
                  << " g_factor=" << t.g_factor
                  << " omega_r=" << omega_r
                  << " term=" << lambda_col * t.g_factor * t.dEx * omega_r
                  << "\n";
      }
    }
  }

  return sum / omega_initial;
}

static double lambdaPairConversionOMP(ConversionType conversion,
    int Z, int N, int A_p,
    int p_pi, int h_pi, int p_nu, int h_nu,
    double E_tot, double U, double V,
    double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
    double C1, double C2, double C3,
    IntegrationMethod method, int midBins, bool guardBounds,
    int Z_proj)
{
  int n = p_pi + h_pi + p_nu + h_nu;
  double g_p = spDensityProton(Z);
  double g_n = spDensityNeutron(N);
  int Z_comp = Z;
  int N_comp = N;

  ensure_wvol_table(Z_comp, N_comp, A_p, Z_proj, C1, R_pi_nu);
  const WvolTable& tbl = wvol_table();

  double omega_initial = particleHoleStateDensity(
    p_pi, h_pi, p_nu, h_nu, U, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
  if (omega_initial <= 0.0) return 0.0;

  int nexc = std::max(midBins / 2, 2);

  // ---- lambdapinu (ProtonToNeutron) or lambdanupi (NeutronToProton) ----
  TalysOmpChannel ch;
  double L1, L2;
  // Annihilated pair (c in TALYS convention)
  int a_pi, a_hi, a_pn, a_hn;

  if (conversion == ConversionType::ProtonToNeutron) {
    // lambdapinu
    ch = channel_pinu(Z_comp, N_comp);
    double less = pauliCorrection(p_pi-1, h_pi-1, p_nu, h_nu, g_p, g_n);
    L1 = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n) - less;
    L2 = U - less;
    // densh: created pair = (0,0,1,1)
    // densp: (1,0,1,1)
    // annihilated pair for term multiplier = (1,1,0,0)
    a_pi = 1; a_hi = 1; a_pn = 0; a_hn = 0;
  } else {
    // lambdanupi
    ch = channel_nupi(Z_comp, N_comp);
    double less = pauliCorrection(p_pi, h_pi, p_nu-1, h_nu-1, g_p, g_n);
    L1 = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n) - less;
    L2 = U - less;
    // densh: created pair = (1,1,0,0)
    // densp: (1,1,1,0)
    // annihilated pair for term multiplier = (0,0,1,1)
    a_pi = 0; a_hi = 0; a_pn = 1; a_hn = 1;
  }

  double dEx = (L2 - L1) / nexc;
  double sep = talys_sep(ch);

  // Residual state
  int r_pi, r_hi, r_pn, r_hn;
  if (conversion == ConversionType::ProtonToNeutron) {
    r_pi = p_pi - 1; r_hi = h_pi - 1; r_pn = p_nu; r_hn = h_nu;
  } else {
    r_pi = p_pi; r_hi = h_pi; r_pn = p_nu - 1; r_hn = h_nu - 1;
  }

  double sum = 0.;

  for (int i = 0; i < nexc; ++i) {
    double uu = L1 + (static_cast<double>(i) + 0.5) * dEx;
    if (uu <= 0) continue;

    double eopt = std::max(uu - sep, -20.);
    int nen = std::min(10 * NUMEN_TALYS, static_cast<int>(eopt * 10.));
    if (nen < -80) continue;

    double wv = tbl.wvol(ch.k, nen);
    double Weff = tbl.Wompfac_exch * wv;
    double lambda_col = 2. * Weff / HBAR;

    // TALYS density ratio
    // densh = phdens2 of the created pair at uu
    // densp = phdens2 of (1,0,1,1) for pinu or (1,1,1,0) for nupi at uu
    if (conversion == ConversionType::ProtonToNeutron) {
      double densh = particleHoleStateDensity(
        0,0,1,1, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
      double densp = particleHoleStateDensity(
        1,0,1,1, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
      double ratio = (densp > 1.) ? densh / densp : 1.;
      lambda_col *= ratio;
    } else {
      double densh = particleHoleStateDensity(
        1,1,0,0, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
      double densp = particleHoleStateDensity(
        1,1,1,0, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
      double ratio = (densp > 1.) ? densh / densp : 1.;
      lambda_col *= ratio;
    }

    double omega_a = particleHoleStateDensity(
      a_pi, a_hi, a_pn, a_hn, uu, Z_comp, N_comp, 0, V, 0.0, false, 1.0);
    double omega_r = particleHoleStateDensity(
      r_pi, r_hi, r_pn, r_hn, U - uu,
      Z_comp, N_comp, 0, V, 0.0, false, 1.0);

    sum += lambda_col * omega_a * dEx * omega_r;

    // Diagnostic for first bin
    if (diagnostic_ && i == 0) {
      std::cerr << "# DIAG conversion type="
                << (conversion == ConversionType::ProtonToNeutron
                    ? "ProtonToNeutron" : "NeutronToProton")
                << " p_pi=" << p_pi << " h_pi=" << h_pi
                << " p_nu=" << p_nu << " h_nu=" << h_nu
                << " k=" << ch.k
                << " Z_res=" << ch.Z_res << " A_res=" << ch.A_res
                << " pdg=" << ch.pdg
                << " sep=" << sep
                << " L1=" << L1 << " L2=" << L2 << " dEx=" << dEx
                << " uu=" << uu
                << " eopt=" << eopt << " nen=" << nen
                << " wvol=" << wv << " Weff=" << Weff
                << " lambda_col=" << lambda_col
                << " omega_a=" << omega_a
                << " omega_r=" << omega_r
                << " term=" << lambda_col * omega_a * dEx * omega_r
                << "\n";
    }
  }

  return sum / omega_initial;
}

// ============================================================================
// Existing functions (unchanged below this line)
// ============================================================================

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
    // A_compound = Z + N is the compound mass (projectile absorbed)
    int A_compound = Z + N;

    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);

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
        double matrix_term = n_proton * g_p * M2_element(MatrixElement::PiPi, A_p, A_compound, E_tot, n) + 2 * n_neutron * g_n * M2_element(MatrixElement::PiNu, A_p, A_compound, E_tot, n);
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
        double matrix_term = n_neutron * g_n * M2_element(MatrixElement::NuNu, A_p, A_compound, E_tot, n) + 2 * n_proton * g_p * M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n);
        return prefactor * numerator * matrix_term * well_term / denom;
    }
}

double lambdaNewPairNumerical(ExcitonType particle, int Z, int N, int A_p,
                              int p_pi, int h_pi, int p_nu, int h_nu,
                              double E_tot, double U, double V,
                              double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                              double C1, double C2, double C3,
                              IntegrationMethod method, int midBins,
                              bool guardBounds)
{
    int n = p_pi + h_pi + p_nu + h_nu;
    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);
    double prefactor = 2 * PI / HBAR;

    int A_compound = Z + N;

    // turning off pairing bc passing in U, not E_x
    double omega_initial = particleHoleStateDensity(
        p_pi, h_pi, p_nu, h_nu, U, Z, N, 0, V, 0.0, false, 1.0);
    if (omega_initial <= 0.0)
        return 0.0;

    marley::Integrator integrator(50);
    double sum = 0.0;

    auto integrate_sub = [&](int t_pi, int t_hi, int t_pn, int t_hn,
                             int r_pi, int r_hi, int r_pn, int r_hn,
                             double M2_val,
                             double g_factor,
                             double L1, double L2) -> double
    {
        auto integrand = [&](double e) -> double
        {
            double omega_t = particleHoleStateDensity(
                t_pi, t_hi, t_pn, t_hn, e, Z, N, 0, V, 0.0, false, 1.0);
            double omega_r = particleHoleStateDensity(
                r_pi, r_hi, r_pn, r_hn, U - e, Z, N, 0, V, 0.0, false, 1.0);
            return prefactor * M2_val * omega_t * g_factor * omega_r;
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
                                    R_pi_pi, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);
        double m2_nupi = M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n,
                                    R_pi_pi, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        // Equation (13.17)
        sum += integrate_sub(2, 1, 0, 0, p_pi - 1, h_pi, p_nu, h_nu, m2_pipi, g_p, L1pip, L2pip);
        sum += integrate_sub(1, 2, 0, 0, p_pi, h_pi - 1, p_nu, h_nu, m2_pipi, g_p, L1pih, L2pih);
        sum += integrate_sub(1, 1, 1, 0, p_pi, h_pi, p_nu - 1, h_nu, m2_nupi, g_n, L1nup, L2nup);
        sum += integrate_sub(1, 1, 0, 1, p_pi, h_pi, p_nu, h_nu - 1, m2_nupi, g_n, L1nuh, L2nuh);
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
                                    R_pi_pi, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        sum += integrate_sub(0, 0, 2, 1, p_pi, h_pi, p_nu - 1, h_nu, m2_nunu, g_n, L1nup, L2nup);
        sum += integrate_sub(0, 0, 1, 2, p_pi, h_pi, p_nu, h_nu - 1, m2_nunu, g_n, L1nuh, L2nuh);
        sum += integrate_sub(1, 0, 1, 1, p_pi - 1, h_pi, p_nu, h_nu, m2_pinu, g_p, L1pip, L2pip);
        sum += integrate_sub(0, 1, 1, 1, p_pi, h_pi - 1, p_nu, h_nu, m2_pinu, g_p, L1pih, L2pih);
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
    int A_compound = Z + N;

    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);

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
                                     double C1, double C2, double C3,
                                     IntegrationMethod method, int midBins,
                                     bool guardBounds)
{
    int n = p_pi + h_pi + p_nu + h_nu;

    double g_p = spDensityProton(Z);
    double g_n = spDensityNeutron(N);

    double prefactor = 2 * PI / HBAR;

    int A_compound = Z + N;

    double omega_initial = particleHoleStateDensity(
        p_pi, h_pi, p_nu, h_nu, U, Z, N, 0, V, 0.0, false, 1.0);
    if (omega_initial <= 0.0)
        return 0.0;

    marley::Integrator integrator(50);

    auto integrate_conversion = [&](int c_pi, int c_hi, int c_pn, int c_hn,
                                    int a_pi, int a_hi, int a_pn, int a_hn,
                                    int r_pi, int r_hi, int r_pn, int r_hn,
                                    double M2_val,
                                    double L1, double L2) -> double
    {
        auto integrand = [&](double e) -> double
        {
            double omega_c = particleHoleStateDensity(
                c_pi, c_hi, c_pn, c_hn, e, Z, N, 0, V, 0.0, false, 1.0);
            double omega_a = particleHoleStateDensity(
                a_pi, a_hi, a_pn, a_hn, e, Z, N, 0, V, 0.0, false, 1.0);
            double omega_r = particleHoleStateDensity(
                r_pi, r_hi, r_pn, r_hn, U - e, Z, N, 0, V, 0.0, false, 1.0);
            return prefactor * M2_val * omega_c * omega_a * omega_r;
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
                                    R_pi_pi, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        double integral = integrate_conversion(
            0, 0, 1, 1, 1, 1, 0, 0, p_pi - 1, h_pi - 1, p_nu, h_nu, m2_pinu, L1, L2);
        return integral / omega_initial;
    }
    else
    {
        double less_neutron_pair = pauliCorrection(p_pi, h_pi, p_nu - 1, h_nu - 1, g_p, g_n);
        double L1 = pauliCorrection(p_pi, h_pi, p_nu, h_nu, g_p, g_n) - less_neutron_pair;
        double L2 = U - less_neutron_pair;

        double m2_nupi = M2_element(MatrixElement::NuPi, A_p, A_compound, E_tot, n,
                                    R_pi_pi, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3, PreeqMode::Numerical);

        double integral = integrate_conversion(
            1, 1, 0, 0, 0, 0, 1, 1, p_pi, h_pi, p_nu - 1, h_nu - 1, m2_nupi, L1, L2);
        return integral / omega_initial;
    }
}

double lambdaRate(LambdaType type, PreeqMode mode, int Z, int N, int A_p,
                  int p_pi, int h_pi, int p_nu, int h_nu,
                  double E_tot, double U, double V,
                  double R_nu_nu, double R_nu_pi, double R_pi_pi, double R_pi_nu,
                  double C1, double C2, double C3,
                  IntegrationMethod method, int midBins,
                  bool guardBounds,
                  CollisionKernel kernel,
                  int Z_proj)
{
    int n = p_pi + h_pi + p_nu + h_nu;

    if (kernel == CollisionKernel::OpticalModel) {
      // OMP kernel forces numerical integration
      switch (type)
      {
      case LambdaType::ProtonPairCreation:
        return lambdaNewPairOMP(ExcitonType::Proton, Z, N, A_p,
          p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
          R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
          method, midBins, guardBounds, Z_proj);
      case LambdaType::NeutronPairCreation:
        return lambdaNewPairOMP(ExcitonType::Neutron, Z, N, A_p,
          p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
          R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
          method, midBins, guardBounds, Z_proj);
      case LambdaType::ProtonToNeutronConversion:
        return lambdaPairConversionOMP(ConversionType::ProtonToNeutron, Z, N, A_p,
          p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
          R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
          method, midBins, guardBounds, Z_proj);
      default:
        return lambdaPairConversionOMP(ConversionType::NeutronToProton, Z, N, A_p,
          p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
          R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
          method, midBins, guardBounds, Z_proj);
      }
    }

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
                                            R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                            method, midBins, guardBounds);
        break;

    case LambdaType::NeutronPairCreation:
        if (useAnalytical)
            result = lambdaNewPairAnalytical(ExcitonType::Neutron, Z, N, A_p,
                                             p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                             R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaNewPairNumerical(ExcitonType::Neutron, Z, N, A_p,
                                            p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                            R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                            method, midBins, guardBounds);
        break;

    case LambdaType::ProtonToNeutronConversion:
        if (useAnalytical)
            result = lambdaPairConversionAnalytical(ConversionType::ProtonToNeutron, Z, N, A_p,
                                                    p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaPairConversionNumerical(ConversionType::ProtonToNeutron, Z, N, A_p,
                                                   p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                   R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                                   method, midBins, guardBounds);
        break;

    default:
        if (useAnalytical)
            result = lambdaPairConversionAnalytical(ConversionType::NeutronToProton, Z, N, A_p,
                                                    p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                    R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3);
        else
            result = lambdaPairConversionNumerical(ConversionType::NeutronToProton, Z, N, A_p,
                                                   p_pi, h_pi, p_nu, h_nu, E_tot, U, V,
                                                   R_nu_nu, R_nu_pi, R_pi_pi, R_pi_nu, C1, C2, C3,
                                                   method, midBins, guardBounds);
        break;
    }

    // TALYS matrix.f90 applies 1.20 M^2 factor only for preeqmode==1.
    // When n==1 forces analytical formula in numerical mode, cancel the factor.
    if (n == 1 && mode == PreeqMode::Numerical)
        result /= 1.20;

    return result;
}
