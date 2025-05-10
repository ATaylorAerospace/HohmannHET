#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

namespace {
    // Calculation Constants
    constexpr double MU = 398600.4418;    // Earth's gravitational parameter (km^3/s^2)
    constexpr double R_EARTH = 6378.137;  // Earth's equatorial radius (km)
}

class HohmannTransfer {
private:
    const double r1;            // Initial orbital radius (km)
    const double r2;            // Final orbit radius (km)
    const double a_transfer;    // Transfer orbit semi-major axis (km)
    const double v1;            // Initial orbital velocity (km/s)
    const double delta_v_departure;  // Delta-V for departure burn (km/s)
    const double delta_v_arrival;    // Delta-V for arrival burn (km/s)

    // Calculate all velocities,  delta-Vs at construction
    static std::tuple<double, double> calculateDeltaVs(double r1, double r2, double a_transfer, double v1) {
        // Calculate transfer orbit velocities
        double v_transfer_1 = std::sqrt(MU * (2.0/r1 - 1.0/a_transfer));
        double v_transfer_2 = std::sqrt(MU * (2.0/r2 - 1.0/a_transfer));
        double v2 = std::sqrt(MU / r2);
        
        return {v_transfer_1 - v1, v2 - v_transfer_2};
    }

public:
    // Constructor with member initializer list
    HohmannTransfer(double initial_altitude, double final_altitude) 
        : r1(initial_altitude + R_EARTH)
        , r2(final_altitude + R_EARTH)
        , a_transfer((r1 + r2) * 0.5)
        , v1(std::sqrt(MU / r1)) {
        if (initial_altitude <= 0 || final_altitude <= 0) {
            throw std::invalid_argument("Altitudes must be positive");
        }
        
        auto [dep, arr] = calculateDeltaVs(r1, r2, a_transfer, v1);
        const_cast<double&>(delta_v_departure) = dep;
        const_cast<double&>(delta_v_arrival) = arr;
    }

    // Calculate transfer time (now inlined)
    [[nodiscard]] inline double calculateTransferTime() const {
        return M_PI * std::sqrt(std::pow(a_transfer, 3) / MU);
    }

    // Print transfer details with StringBuilder pattern for better performance
    void printTransferDetails() const {
        std::cout << std::fixed << std::setprecision(2)
                  << "\nHohmann Transfer Orbit Details:\n"
                  << "Initial Orbit Altitude: " << (r1 - R_EARTH) << " km\n"
                  << "Final Orbit Altitude: " << (r2 - R_EARTH) << " km\n"
                  << "Transfer Orbit Semi-Major Axis: " << a_transfer << " km\n"
                  << "Delta-V (Departure Burn): " << delta_v_departure << " km/s\n"
                  << "Delta-V (Arrival Burn): " << delta_v_arrival << " km/s\n"
                  << "Total Delta-V: " << (delta_v_departure + delta_v_arrival) << " km/s\n"
                  << "Transfer Time: " << (calculateTransferTime() / 3600) << " hours\n";
    }

    // Const getters
    [[nodiscard]] double getDeltaVDeparture() const noexcept { return delta_v_departure; }
    [[nodiscard]] double getDeltaVArrival() const noexcept { return delta_v_arrival; }
    [[nodiscard]] double getTransferSemiMajorAxis() const noexcept { return a_transfer; }
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
        
        HohmannTransfer transfer(initial_alt, final_alt);
        transfer.printTransferDetails();
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
