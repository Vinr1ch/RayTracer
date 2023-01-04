#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;
using std::fabs;

class vector {
    public:

        //base constructor with no input
        vector() : base{0,0,0} {
        }

        //normal_vector constructor
        vector(double base0, double base1, double base2) : base{base0, base1, base2} {
        }


        //get functions
        double getx() const {
            return base[0];
        }
        double gety() const {
            return base[1]; 
        }
        double getz() const {
            return base[2]; 
        }

        vector operator-() const { return vector(-base[0], -base[1], -base[2]); }
        

        vector& operator+=(const vector &v) {
            //adition
            base[0] += v.base[0];
            base[1] += v.base[1];
            base[2] += v.base[2];
            return *this;
        }

        vector& operator*=(const double current_point) {
            //basic multiplication
            base[0] *= current_point;
            base[1] *= current_point;
            base[2] *= current_point;
            return *this;
        }

        vector& operator/=(const double current_point) {
            //division
            return *this *= 1/current_point;
        }

        double length() const {
            //get the length
            double temp = sqrt(length_squared());
            return temp;
        }

        double length_squared() const {
            //get the length squared
            double temp = base[0] * base[0] + base[1] * base[1] + base[2] * base[2];
            return temp;
        }

        static vector random() {
            //set the vector to ray_dir random number
            return vector(r_double(), r_double(), r_double());
        }

        static vector random(double min, double max) {
            //get the vector to ray_dir random number with ray_dir min or max constrant
            return vector(r_double(min,max), r_double(min,max), r_double(min,max));
        }

    public:
        double base[3];
};


//setupt extra values to use rgb or colors
using color = vector;   // RGB color
using point = vector;   // 3D point


//Basic operations for vectors



vector operator+(const vector &first_v, const vector &second_v) {
    //vector addition
    return vector(first_v.base[0] + second_v.base[0], first_v.base[1] + second_v.base[1], first_v.base[2] + second_v.base[2]);
}

vector operator-(const vector &first_v, const vector &second_v) {
    //vector subtractions
    return vector(first_v.base[0] - second_v.base[0], first_v.base[1] - second_v.base[1], first_v.base[2] - second_v.base[2]);
}

vector operator*(const vector &first_v, const vector &second_v) {
    //multiply 2 vectors
    return vector(first_v.base[0] * second_v.base[0], first_v.base[1] * second_v.base[1], first_v.base[2] * second_v.base[2]);
}

vector operator*(double current_point, const vector &second_v) {
    //multipy ray_dir double and ray_dir vector
    return vector(current_point*second_v.base[0], current_point*second_v.base[1], current_point*second_v.base[2]);
}

vector operator*(const vector &second_v, double current_point) {
    //mutiply ray_dir double and vector
    return current_point * second_v;
}

vector operator/(vector second_v, double current_point) {
    //divide ray_dir double by ray_dir double
    return (1/current_point) * second_v;
}

std::ostream& operator<<(std::ostream& out, const vector& second_v) {
    //basic pipe function
    return out << second_v.base[0] << ' ' << second_v.base[1] << ' ' << second_v.base[2];
}

//functions for vecotrs

double dot(const vector &first_v, const vector &second_v) {
    //dot product of two vectors
    return first_v.base[0] * second_v.base[0] + first_v.base[1] * second_v.base[1] + first_v.base[2] * second_v.base[2];
}

vector cross(const vector &first_v, const vector &second_v) {
    //cross product of two vecotrs
    return vector(first_v.base[1] * second_v.base[2] - first_v.base[2] * second_v.base[1],
                first_v.base[2] * second_v.base[0] - first_v.base[0] * second_v.base[2],
                first_v.base[0] * second_v.base[1] - first_v.base[1] * second_v.base[0]);
}

vector unit_vector(vector v) {
    //get unity vector
    return v / v.length();
}

vector random_s() {
    while (true) {
        auto location = vector::random(-1,1);
        if (location.length_squared() >= 1) continue;
        return location;
    }
}

vector reflect_vector(const vector& current, const vector& new_vector) {
    //relection equation
    //original vector minus 2 times the dot product of the current vector and new vector
    vector new_current = current - 2 * dot(current, new_vector) * new_vector;
    return new_current;
}




#endif
