#include <iostream>
#include <cmath>

// Transfer Constants
const double mu = 398600.4418;
const double req = 6378.137;
const double pi = 3.141592653589793;

// Conversion factors
const double rtd = 180 / pi;
const double dtr = pi / 180;

// Function to get input safely
double get_input(const std::string& prompt) {
    double value;
    while (true) {
        std::cout << prompt;
        if (std::cin >> value) {
            return value;
        } else {
            std::cout << "Invalid input. Please try again.\n";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
}

// Current Function for Hohmann transfer calculations
double hohmfunc(double dinc1, double v1, double hn1, double hn2, double hn3, double dinc) {
    double dinc2 = dinc - dinc1;
    return v1 * std::sqrt(1.0 + hn1*hn1 - 2 * hn1 * std::cos(dinc1 * dtr)) +
           v1 * std::sqrt(hn2*hn2 * hn3*hn3 + hn2*hn2 - 2 * hn2*hn2 * hn3 * std::cos(dinc2 * dtr));
}

int main() {
    std::cout << "\nHohmann Orbit Transfer Analysis\n";
    double alt1 = get_input("Please input the initial altitude (kilometers): ");
    double alt2 = get_input("Please input the final altitude (kilometers): ");
    double inc1 = get_input("Please input the initial inclination (degrees) (0 <= inclination <= 180): ") * dtr;
    double inc2 = get_input("Please input the final inclination (degrees) (0 <= inclination <= 180): ") * dtr;

    // Current Calculation section 
    
    return 0;
}
