#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

namespace {
    // Calculation Constants
    constexpr double MU = 398600.4418;    // Earth's gravitational parameter (km^3/s^2)
    constexpr double R_EARTH = 6378.137;  // Earth's equatorial radius (km)
    constexpr double INV_3600 = 1.0 / 3600.0;  // Precomputed constant for hours conversion
}

class HohmannTransfer {
private:
    const double r1;            // Initial orbital radius (km)
    const double r2;            // Final orbit radius (km)
    const double a_transfer;    // Transfer orbit semi-major axis (km)
    const double v1;            // Initial orbital velocity (km/s)
    const double delta_v_departure;  // Delta-V for departure burn (km/s)
    const double delta_v_arrival;    // Delta-V for arrival burn (km/s)
    const double transfer_time; // Precomputed transfer time (seconds)

    // Calculate all velocities and delta-Vs at construction
    static std::tuple<double, double, double> calculateValues(double r1, double r2, double a_transfer, double v1) {
        // Calculate transfer orbit velocities using precomputed values
        const double mu_2_r1 = MU * 2.0 / r1;
        const double mu_inv_a = MU / a_transfer;
        const double mu_2_r2 = MU * 2.0 / r2;
        
        const double v_transfer_1 = std::sqrt(mu_2_r1 - mu_inv_a);
        const double v_transfer_2 = std::sqrt(mu_2_r2 - mu_inv_a);
        const double v2 = std::sqrt(MU / r2);
        
        // Calculate transfer time once
        const double transfer_time = M_PI * std::sqrt(a_transfer * a_transfer * a_transfer / MU);
        
        return {v_transfer_1 - v1, v2 - v_transfer_2, transfer_time};
    }

public:
    // Constructor with member initializer list
    HohmannTransfer(double initial_altitude, double final_altitude) 
        : r1(initial_altitude + R_EARTH)
        , r2(final_altitude + R_EARTH)
        , a_transfer((r1 + r2) * 0.5)
        , v1(std::sqrt(MU / r1))
        , delta_v_departure([&]() {
            if (initial_altitude <= 0 || final_altitude <= 0) {
                throw std::invalid_argument("Altitudes must be positive");
            }
            auto [dep, arr, time] = calculateValues(r1, r2, a_transfer, v1);
            const_cast<double&>(delta_v_arrival) = arr;
            const_cast<double&>(transfer_time) = time;
            return dep;
        }())
        , delta_v_arrival(0.0)  // Will be set in lambda above
        , transfer_time(0.0)    // Will be set in lambda above
    {}

    // Get precomputed transfer time
    [[nodiscard]] constexpr double getTransferTime() const noexcept { 
        return transfer_time; 
    }

    // Optimized print function with fewer calculations
    void printTransferDetails() const {
        const double total_delta_v = delta_v_departure + delta_v_arrival;
        const double transfer_hours = transfer_time * INV_3600;
        
        std::cout << std::fixed << std::setprecision(2)
                  << "\nHohmann Transfer Orbit Details:\n"
                  << "Initial Orbit Altitude: " << (r1 - R_EARTH) << " km\n"
                  << "Final Orbit Altitude: " << (r2 - R_EARTH) << " km\n"
                  << "Transfer Orbit Semi-Major Axis: " << a_transfer << " km\n"
                  << "Delta-V (Departure Burn): " << delta_v_departure << " km/s\n"
                  << "Delta-V (Arrival Burn): " << delta_v_arrival << " km/s\n"
                  << "Total Delta-V: " << total_delta_v << " km/s\n"
                  << "Transfer Time: " << transfer_hours << " hours\n";
    }

    // Const getters
    [[nodiscard]] constexpr double getDeltaVDeparture() const noexcept { return delta_v_departure; }
    [[nodiscard]] constexpr double getDeltaVArrival() const noexcept { return delta_v_arrival; }
    [[nodiscard]] constexpr double getTransferSemiMajorAxis() const noexcept { return a_transfer; }
    [[nodiscard]] constexpr double getTotalDeltaV() const noexcept { return delta_v_departure + delta_v_arrival; }
};

int main() {
    try {
        double initial_alt, final_alt;
        
        std::cout << "Enter initial orbit altitude (km): ";
        if (!(std::cin >> initial_alt)) {
            throw std::invalid_argument("Invalid initial altitude input");
        }
        
        std::cout << "Enter final orbit altitude (km): ";
        if (!(std::cin >> final_alt)) {
            throw std::invalid_argument("Invalid final altitude input");
        }
        
        const HohmannTransfer transfer(initial_alt, final_alt);
        transfer.printTransferDetails();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    
    return 0;
}
