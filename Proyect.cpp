#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <utility>

struct Particle {
    double r, t, phi;
    int sign_pr;
    
    void flip_sign() { sign_pr *= -1; }
};

std::pair<double, double> rk4_step(const Particle& p, double dt, 
                                  double M, double a, double E, double L);
double horizon_radius(double M, double a) { return (M + std::sqrt(M*M - a*a)); }

constexpr double EPSILON = 1e-8;
constexpr double MIN_RADIUS = 1e-2;

double safe_sqrt(double x) {
    return x >= 0 ? std::sqrt(x) : std::numeric_limits<double>::quiet_NaN();
}

// Effective potential R(r) derived from Chandrasekar
double compute_R(double r, double M, double a, double E, double L) {
    const double Delta = r*r - 2*M*r + a*a;
    const double term1 = (2*M/r) * (a*E - L) * (a*E - L);
    const double term2 = E*E*r*r + (a*a*E*E - L*L) - Delta;
    return term1 + term2;
}
// Evolution equation for r from Chandrasekar 
double drdtao(double r, double M, double a, double E, double L) {
    if (r < MIN_RADIUS) return std::numeric_limits<double>::quiet_NaN();
    const double R = compute_R(r, M, a, E, L);
    return R >= 0 ? std::sqrt(R) / std::abs(r) : std::numeric_limits<double>::quiet_NaN();
}
// Evolution equation for \phi from Chandrasekar 
double dphidtao(double r, double M, double a, double E, double L) {
    const double Delta = r*r - 2*M*r + a*a;
    if (std::abs(Delta) < EPSILON || r < MIN_RADIUS) 
        return std::numeric_limits<double>::quiet_NaN();
    return (E*a*2*M + L*(r - 2*M)) / (r * Delta);
}

// Compute allowed range for L at a given initial r such that R(r)==0.
// The effective potential R(r)=0 can be rearranged as a quadratic in L:
// A*L^2 + B*L + C = 0.
std::pair<double, double> allowed_L_range(double r, double M, double a, double E) {
    // Expand:
    // R(r) = (2*M/r)*(aE - L)^2 + E^2*r^2 + a^2E^2 - L^2 - (r^2 - 2*M*r + a*a) = 0
    // This gives:
    //   A*L^2 + B*L + C = 0,
    double A = (2*M/r - 1.0);
    double B = - (4*M*a*E / r);
    double C = E*E*r*r + a*a*E*E + (2*M*a*a*E*E)/r - (r*r - 2*M*r + a*a);
    
    double disc = B*B - 4*A*C;
    if(disc < 0)
        throw std::runtime_error("No allowed L values for the given parameters.");
    double sqrt_disc = std::sqrt(disc);
    double L1 = (-B - sqrt_disc) / (2*A);
    double L2 = (-B + sqrt_disc) / (2*A);
    return std::make_pair(std::min(L1, L2), std::max(L1, L2));
}

Particle runge_kutta(Particle p, double M, double a, double E, double L,
                     double dt_initial, double t_end, const std::string& filename) {
    std::ofstream fout(filename);
    fout << "t\tr\tphi\n";
    fout << p.t << "\t" << p.r << "\t" << p.phi << "\n";

    const double r_horizon = horizon_radius(M, a);
    double dt = dt_initial;
    bool valid = true;

    while (valid && p.t < t_end) {
        if (p.r <= r_horizon + EPSILON) {
            fout << "# Horizon crossed at t = " << p.t << "\n";
            break;
        }

        bool step_ok = false;
        double current_dt = std::min(dt, t_end - p.t);
        Particle trial = p;

        for (int attempts = 0; attempts < 50 && !step_ok; ++attempts) {
            try {
                const auto [dr, dphi] = rk4_step(trial, current_dt, M, a, E, L);
                
                // Check new state validity
                if (std::isnan(dr) || std::isnan(dphi) || trial.r + dr < MIN_RADIUS)
                    throw std::runtime_error("Invalid step");
                
                trial.r += dr;
                trial.phi += dphi;
                trial.t += current_dt;

                // Check for turning points
                const double R_new = compute_R(trial.r, M, a, E, L);
                if (R_new < -EPSILON) throw std::runtime_error("Negative R");
                if (R_new < EPSILON) trial.flip_sign();

                step_ok = true;
            }
            catch (...) {
                current_dt *= 0.5;
                trial = p;
                if (current_dt < 1e-6) {
                    valid = false;
                    fout << "# Minimum step reached at t = " << p.t << "\n";
                    break;
                }
            }
        }

        if (step_ok) {
            p = trial;
            dt = std::clamp(current_dt*2.0, 1e-6, 0.1);
            fout << p.t << "\t" << p.r << "\t" << p.phi << "\n";
        }
    }
    return p;
}

std::pair<double, double> rk4_step(const Particle& p, double dt, 
                                  double M, double a, double E, double L) {
    auto deriv = [&](double r, double phi, int sign) {
        const double dr = sign * drdtao(r, M, a, E, L);
        const double dp = dphidtao(r, M, a, E, L);
        if (std::isnan(dr) || std::isnan(dp)) 
            throw std::runtime_error("NaN in derivatives");
        return std::make_pair(dr, dp);
    };

    const auto [k1r, k1p] = deriv(p.r, p.phi, p.sign_pr);
    const auto [k2r, k2p] = deriv(p.r + 0.5*dt*k1r, p.phi + 0.5*dt*k1p, p.sign_pr);
    const auto [k3r, k3p] = deriv(p.r + 0.5*dt*k2r, p.phi + 0.5*dt*k2p, p.sign_pr);
    const auto [k4r, k4p] = deriv(p.r + dt*k3r, p.phi + dt*k3p, p.sign_pr);

    return {
        dt/6.0 * (k1r + 2*k2r + 2*k3r + k4r),
        dt/6.0 * (k1p + 2*k2p + 2*k3p + k4p)
    };
}

int main() {
    try {
        const double M = 1.0, a = 0.2;

        const double E = 0.97;
        if (a >= M) throw std::invalid_argument("a must be < M");

        // Initial radius and particle state
        Particle p{10.0, 0.0, 0.0, -1};
        
        // Precompute allowed L range at the initial r
        auto L_range = allowed_L_range(p.r, M, a, E);
        std::cout << "Allowed L range at r = " << p.r << " is: [" 
                  << L_range.first << ", " << L_range.second << "]\n";
        
       
        double L_direct = a*E;
       
        double L_retro = -L_direct;

        // Run the direct orbit simulation
        runge_kutta(p, M, a, E, L_direct, 0.001, 5000.0, "direct_orbit.txt");
        
        // Reset particle initial state for retrograde simulation
        p = Particle{10.0, 0.0, 0.0, -1};
        runge_kutta(p, M, a, E, L_retro, 0.001, 5000.0, "retrograde_orbit.txt");
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
