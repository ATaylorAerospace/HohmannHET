#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>

// Calculation Constants
const double MU = 398600.4418;    // Earth's gravitational parameter (km^3/s^2)
const double R_EARTH = 6378.137;  // Earth's equatorial radius (km)

class HohmannTransfer {
private:
    double r1;            // Initial orbital radius (km)
    double r2;            // Final orbital radius (km)
    double a_transfer;    // Transfer orbit semi-major axis (km)
    double v1;            // Initial orbital velocity (km/s)
    double delta_v_departure;  // Delta-V for departure burn (km/s)
    double delta_v_arrival;    // Delta-V for arrival burn (km/s)

    // Private method to calculate departure delta-V
    double calculateDeltaVDeparture() {
        double v_transfer_1 = std::sqrt(MU * (2.0/r1 - 1.0/a_transfer));
        return v_transfer_1 - v1;
    }

    // Private method to calculate arrival delta-V
    double calculateDeltaVArrival() {
        // Velocity at final orbit
        double v2 = std::sqrt(MU / r2);
        
        // Velocity at transfer orbit apoapsis
        double v_transfer_2 = std::sqrt(MU * (2.0/r2 - 1.0/a_transfer));
        
        return v2 - v_transfer_2;
    }

public:
    // Constructor
    HohmannTransfer(double initial_altitude, double final_altitude) {
        // Validate inputs
        if (initial_altitude <= 0 || final_altitude <= 0) {
            throw std::invalid_argument("Altitudes must be positive");
        }
        
        // Calculate orbital radii
        r1 = R_EARTH + initial_altitude;
        r2 = R_EARTH + final_altitude;
        
        // Calculate transfer orbit semi-major axis
        a_transfer = (r1 + r2) / 2.0;
        
        // Calculate initial orbital velocity
        v1 = std::sqrt(MU / r1);
        
        // Compute delta-V for transfer
        delta_v_departure = calculateDeltaVDeparture();
        delta_v_arrival = calculateDeltaVArrival();
    }

    // Calculate transfer time
    double calculateTransferTime() const {
        return M_PI * std::sqrt(std::pow(a_transfer, 3) / MU);
    }

    // Print transfer details
    void printTransferDetails() const {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "\nHohmann Transfer Orbit Details:\n";
        std::cout << "Initial Orbit Altitude: " << r1 - R_EARTH << " km\n";
        std::cout << "Final Orbit Altitude: " << r2 - R_EARTH << " km\n";
        std::cout << "Transfer Orbit Semi-Major Axis: " << a_transfer << " km\n";
        std::cout << "Delta-V (Departure Burn): " << delta_v_departure << " km/s\n";
        std::cout << "Delta-V (Arrival Burn): " << delta_v_arrival << " km/s\n";
        std::cout << "Total Delta-V: " << delta_v_departure + delta_v_arrival << " km/s\n";
        std::cout << "Transfer Time: " << calculateTransferTime() / 3600 << " hours\n";
    }

    // Getters for key parameters
    double getDeltaVDeparture() const { return delta_v_departure; }
    double getDeltaVArrival() const { return delta_v_arrival; }
    double getTransferSemiMajorAxis() const { return a_transfer; }
};

int main() {
    try {
        double initial_alt, final_alt;

        // User inputs
        std::cout << "Enter initial orbit altitude (km): ";
        std::cin >> initial_alt;

        std::cout << "Enter final orbit altitude (km): ";
        std::cin >> final_alt;
        
        // Create Hohmann transfer object
        HohmannTransfer transfer(initial_alt, final_alt);
        
        // Print transfer details
        transfer.printTransferDetails();
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Please enter valid numerical inputs.\n";
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
