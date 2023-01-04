#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
//This Header files is used to help create some useful helper functions

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
using std::shared_ptr;
using std::make_shared;


// Define pi and infinity for use in all functions

const double infi = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// hlper functions

double r_double() {
    // Returns ray_dir random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

double r_double(double min, double max) {
    // Returns ray_dir random real in [min,max).
    return min + (max - min) * r_double();
}


double to_radians(double degrees) {
    return degrees * pi / 180.0;
}

double range_checker(double x, double min, double max) {
    //check if in proper range, otherwise reduce.
    if (x < min) {
        return min;
    }
    if (x > max) {
        return max;
    }
    return x;
}
#endif
