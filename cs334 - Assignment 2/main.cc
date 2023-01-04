//CS334 - Fundamentals of Computer Graphics
//By Colin Vinarcik

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
#include <cmath>
#include <memory>
#include "helperFunctions.h"
#include "vector.h"

//creation of struct and material for use
struct hit_rec;
class m_type;


class ray {
    //ray object
    //default two vecotrs
    public:
        point o;
        vector d;
    public:
        //default constructor
        ray() {
        }

        //constructor with point and direction
        ray(const point& hold, const vector& there) : o(hold), d(there){
        }

        //basic getter functions
        point origin() const {
            return o; 
        }
        vector direction() const {
            return d; 
        }
        
        point at(double current_point) const {
            return o + current_point * d;
        }       
};

struct hit_rec {
    //vairables to hold the locations, what is hit
    double current_point;
    bool visible;    
    point location;
    vector normal_vector;
    shared_ptr<m_type> type_ptr;
    
    //set function for normal vector.
    inline void set_face_normal(const ray& r, const vector& normal_check) {
        //can we see the face
        visible = dot(r.direction(), normal_check) < 0;
        //if yest then set the normal vector
        normal_vector = visible ? normal_check : -normal_check;
    }
};


class check_hit {
    //helper for ray hit
    public:
        virtual bool ray_hit(const ray& r, double t_min, double t_max, hit_rec& rec) const = 0;

};


class hit_lit : public check_hit {
    //ray_dir class to store the hit values
    public:
        //default constructor
        hit_lit() {
        }
        //take an object as something that is hittable
        hit_lit(shared_ptr<check_hit> object) {
            add(object); 
        }

        //add ray_dir new object to the hit values
        void add(shared_ptr<check_hit> object) {
            objects.push_back(object); 
        }

        //get the rayhit boolean
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;


    public:
        std::vector<shared_ptr<check_hit>> objects;
};

class camera {
    //private varibles to be used 
private:
    point origin;
    point corner_bl;
    vector horizontal;
    vector vertical;
    vector unit_vect, vert_temp, widt_temp;
    double lens_radius;
public:
    //default value for camera object
    camera() : camera(point(0, 0, -1), point(0, 0, 0), vector(0, 1, 0), 40, 1, 1) {}

    //camera construcotr
    camera(point origin, point view, vector default_view, double feild_of_view, double aspect, double focus_dist) {

        //get the unity vecotrs for the viewport
        widt_temp = unit_vector(origin - view);
        unit_vect = unit_vector(cross(default_view, widt_temp));
        vert_temp = cross(widt_temp, unit_vect);


        //get theta value for get the correct feild of view
        auto theta = to_radians(feild_of_view);
        //get true hight
        auto h = tan(theta / 2);
        //widow height
        auto w_height = 2.0 * h;
        //window length
        auto w_width = aspect * w_height;

        //calculate the view distance for the sence
        horizontal = focus_dist * w_width * unit_vect;
        vertical = focus_dist * w_height * vert_temp;
        corner_bl = origin - horizontal / 2 - vertical / 2 - focus_dist * widt_temp;
        //set default lens radius
        lens_radius = 0;

    }

    //get function for ray
    ray get_ray(double s, double current_point) const {
        //get random data from unit vector, and offectdata
        vector rd = lens_radius * vector(0, 0, 0);
        vector offset = unit_vect * rd.getx() + vert_temp * rd.gety();
        //(printf("testing %f %f %f", rd.getx(), rd.getx(), rd.getx())
        //return ray based on origin.
        return ray(origin + offset, corner_bl + s * horizontal + current_point * vertical - origin - offset);
        //exit(0);
    }


};

class m_type {
    public:
        //variables
        double sharp;
        color new_color;
    public:

        //basic constructor.
        m_type(const color& ray_dir, double f) : new_color(ray_dir), sharp(f < 1 ? f : 1) {
        }

        //scatter function, uesed to find if there is ray_dir reflection for the object to have.
        virtual bool scatter(const ray& r_in, const hit_rec& rec, color& current_color, ray& reflect) {
            //get the reflected vector
            vector reflected = reflect_vector(unit_vector(r_in.direction()), rec.normal_vector);
            //make ray_dir ray of the relected object
            //reflect = ray(rec.location, reflected + sharp * 5.0);
            reflect = ray(rec.location, reflected + sharp * random_s());
            //set the objects reflection to the new color found
            current_color = new_color;
            //return if the dot of the new rays directo and hit rec is greater than 0
            return (dot(reflect.direction(), rec.normal_vector) > 0);
        }
};


class sphere : public check_hit {
    //variables
    public:
        point center;
        double radius;
        shared_ptr<m_type> type_ptr;

    public:
        //constructors
        sphere() {
        }

        //sphere(point cen, double curent_ray, shared_ptr<m_type> m): center(cen), radius(curent_ray)
        sphere(point cen, double curent_ray, shared_ptr<m_type> m): center(cen), radius(curent_ray), type_ptr(m) {
        };

        //check if hit
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;

};


class a_a_box : public check_hit {
    //variables
    public:
        point box_min;
        point box_max;
        hit_lit sides;
    public:
        //constructors
        a_a_box() {
        }
        //constructors with given points and how reflective
        a_a_box(const point& min, const point& max, shared_ptr<m_type> ptr);

        //check if hit
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;
   
};


class xy_box : public check_hit {
    //variables
    public:
        shared_ptr<m_type> mp;
        double x0, x1, y0, y1, other_axis;
    public:
        //constructors
        xy_box() {}

        //constructors with given points and how reflective
        xy_box(double minx, double maxx, double miny, double maxy, double _k, shared_ptr<m_type> refl_am)
            : x0(minx), x1(maxx), y0(miny), y1(maxy), other_axis(_k), mp(refl_am) {};


        //check if hit
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;


    
};

class xz_box : public check_hit {
    //variables
    public:
        shared_ptr<m_type> mp;
        double x0, x1, z0, z1, other_axis;
    public:
        //constructors
        xz_box() {}

        //constructors with given points and how reflective
        xz_box(double minx, double maxx, double minz, double maxz, double _k,shared_ptr<m_type> refl_am)
            : x0(minx), x1(maxx), z0(minz), z1(maxz), other_axis(_k), mp(refl_am) {};


        //check if hit
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;
};

class yz_box : public check_hit {
    //variables
    public:
        shared_ptr<m_type> mp;
        double y0, y1, z0, z1, other_axis;
    public:
        //constructors
        yz_box() {}

        //constructors with given points and how reflective
        yz_box(double miny, double maxy, double minz, double maxz, double _k, shared_ptr<m_type> refl_am)
            : y0(miny), y1(maxy), z0(minz), z1(maxz), other_axis(_k), mp(refl_am) {};


        //check if hit
        virtual bool ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const override;



};


bool hit_lit::ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const {
    //temp values
    hit_rec temp_rec;
    //default nothing is hit
    bool anything_hit = false;
    //how far
    double closest_so_far = t_max;

    //check all objects for hit
    for (const auto& object : objects) {
        //if object ray hits 
        if (object->ray_hit(curent_ray, t_min, closest_so_far, temp_rec)) {
            //it hit something, update values accordingly
            anything_hit = true;
            closest_so_far = temp_rec.current_point;
            rec = temp_rec;
        }
    }

    //return bool value for it hits something
    return anything_hit;
}

bool xy_box::ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const {

    //get the point of the object
    auto current_point = (other_axis - curent_ray.origin().getz()) / curent_ray.direction().getz();

    //is current point larger than our max and min 
    if (current_point < t_min || current_point > t_max) {
        return false;
    }

    //get the the x and y points
    auto x = curent_ray.origin().getx() + current_point * curent_ray.direction().getx();
    auto y = curent_ray.origin().gety() + current_point * curent_ray.direction().gety();

    //check x and y values are not larger or smaller than max and min
    if (x < x0 || x > x1 || y < y0 || y > y1) {
        return false;
    }

    //update the current point
    rec.current_point = current_point;

    //get the normal and update it for the ray
    vector normal_check = vector(0, 0, 1);
    rec.set_face_normal(curent_ray, normal_check);

    //update pointers
    rec.type_ptr = mp;
    rec.location = curent_ray.at(current_point);

    //return true for hit
    return true;
}

bool xz_box::ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const {

    //get the point of the object
    auto current_point = (other_axis - curent_ray.origin().gety()) / curent_ray.direction().gety();
    //is current point larger than our max and min 
    if (current_point < t_min || current_point > t_max) {
        return false;
    }

    //get the the x and z points
    auto x = curent_ray.origin().getx() + current_point * curent_ray.direction().getx();
    auto z = curent_ray.origin().getz() + current_point * curent_ray.direction().getz();
    
    //check x and z values are not larger or smaller than max and min
    if (x < x0 || x > x1 || z < z0 || z > z1) {
        return false;
    }

    //update the current point
    rec.current_point = current_point;

    //get the normal and update it for the ray
    vector normal_check = vector(0, 1, 0);
    rec.set_face_normal(curent_ray, normal_check);

    //update pointers
    rec.type_ptr = mp;
    rec.location = curent_ray.at(current_point);

    //return true for hit
    return true;
}

bool yz_box::ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const {

    //get the point of the object
    auto current_point = (other_axis - curent_ray.origin().getx()) / curent_ray.direction().getx();

    //is current point larger than our max and min 
    if (current_point < t_min || current_point > t_max) {
        return false;
    }

    //get the the y and z points
    auto y = curent_ray.origin().gety() + current_point * curent_ray.direction().gety();
    auto z = curent_ray.origin().getz() + current_point * curent_ray.direction().getz();

    //check z and y values are not larger or smaller than max and min
    if (y < y0 || y > y1 || z < z0 || z > z1) {
        return false;
    }

    //update the current point
    rec.current_point = current_point;

    //get the normal and update it for the ray
    vector normal_check = vector(1, 0, 0);
    rec.set_face_normal(curent_ray, normal_check);

    //update pointers
    rec.type_ptr = mp;
    rec.location = curent_ray.at(current_point);

    //return true for hit
    return true;
}



a_a_box::a_a_box(const point& min, const point& max, shared_ptr<m_type> ptr) {

    //set box points min and maxs
    box_min = min;
    box_max = max;

    //add the xy planes
    sides.add(make_shared<xy_box>(min.getx(), max.getx(), min.gety(), max.gety(), max.getz(), ptr));
    sides.add(make_shared<xy_box>(min.getx(), max.getx(), min.gety(), max.gety(), min.getz(), ptr));

    //add the xz planes
    sides.add(make_shared<xz_box>(min.getx(), max.getx(), min.getz(), max.getz(), max.gety(), ptr));
    sides.add(make_shared<xz_box>(min.getx(), max.getx(), min.getz(), max.getz(), min.gety(), ptr));

    //add the yz planes
    sides.add(make_shared<yz_box>(min.gety(), max.gety(), min.getz(), max.getz(), max.getx(), ptr));
    sides.add(make_shared<yz_box>(min.gety(), max.gety(), min.getz(), max.getz(), min.getx(), ptr));
}

bool a_a_box::ray_hit(const ray& curent_ray, double t_min, double t_max, hit_rec& rec) const {
    //return the value from the ray hit sides functions
    return sides.ray_hit(curent_ray, t_min, t_max, rec);
}



void output_ppm(int samples , std::ostream& output_stream, color new_color) {
    //Caluclate sharpness
    double sharpness = 1.0 / samples;
    //get the color and stor in vairables
    double r = new_color.getx();
    double g = new_color.gety();
    double b = new_color.getz();

    
    // find scale value, and sharpen based on samples input
    r = sqrt(sharpness * r);
    g = sqrt(sharpness * g);
    b = sqrt(sharpness * b);
    //Printf("%f  %f  %f" r, g, b);
    //exit(0);
    //outputstream << r  << ' ' <<  g << ' ' << b

    // output to the stream.
    output_stream << static_cast<int>(256 * range_checker(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * range_checker(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * range_checker(b, 0.0, 0.999)) << '\n';
}

bool sphere::ray_hit(const ray& current_ray, double t_min, double t_max, hit_rec& rec) const {

    //get the values for the sphere
    vector sphere_o = current_ray.origin() - center;
    double ray_dir = current_ray.direction().length_squared();
    double dot_s_d = dot(sphere_o, current_ray.direction());
    double constant = sphere_o.length_squared() - radius * radius;

    //get the discrimeant and check is its less than zero
    double discr = dot_s_d * dot_s_d - ray_dir * constant;
    if (discr < 0) {
        //return false on hit
        return false;
    }

    //get the sqrt of the d
    double sqrt_d = sqrt(discr);

    // get the root that is within the proper range
    double root = (-dot_s_d - sqrt_d) / ray_dir;

    //check the new root is larger or smaller than the min and max
    if (root < t_min || t_max < root) {
        //return false on hit
        return false;
       
    }

    //update the current point and its location
    rec.current_point = root;
    rec.location = current_ray.at(rec.current_point);
    
    //update normals
    vector normal_check = (rec.location - center) / radius;
    rec.set_face_normal(current_ray, normal_check);

    //update pointers
    rec.type_ptr = type_ptr;


    //return true for hit
    return true;
}


color get_color_ray(const ray& current,  const check_hit& world, int bounces) {
    hit_rec rec;

    // Check for the number of bounes
    if (bounces <= 0) {
        return color(0, 0, 0);
    }
   
    //
    if (world.ray_hit(current, 0.001, infi, rec)) {
        //default values
        ray reflect;
        color current_color;
        //check the the object being hi and get ray color, otherwise return black
        if (rec.type_ptr->scatter(current, rec, current_color, reflect)) {
            return current_color * get_color_ray(reflect, world, bounces - 1);
        }
        return color(0,0,0);
    }

    //get the direction 
    vector unit_direction = unit_vector(current.direction());
    auto current_point = 0.5*(unit_direction.gety() + 1.0);
    return (current_point)*color(1.0, 1.0, 1.0) ;
}


int main(int argc, char** argv) {

    auto name_file = argv[1];
    
    //temp variables to hold the data for use of the creation of the camera.
    int w, h, m;
    double location = 0;

    
    //list of all check_hit objects within the defined world
    hit_lit world;

    //open the files and prepare for reading.
    //this could be written beter, but works for the need
    std::ifstream inFile(name_file);
    if (inFile) {
        //make ray_dir string holder for each lines
        std::string line{};
        //run through the file line by line.
        while (getline(inFile, line)) {
            if (inFile) {
                //get information about the first character in the line to know what is being added.
                char found = line[0];
                //printf("testing this %c\n", found);
                // printf("%s\n", line.c_str());
                if (found == 'I') {
                    //get the image dimesnions  
                    std::istringstream buf(line);
                    std::istream_iterator<std::string> beg(buf), end;
                    //Add item to vector variable for ease of acess 
                    std::vector<std::string> tokens(beg, end); 
                    int counter = 0;
                    //for loop to set data from stirng value.
                    for (auto& s : tokens) {
                        if (counter == 1) {
                            //get image width and convert to int/double 
                            std::stringstream geek(s);
                            geek >> w;

                        } else if (counter == 2) {
                            //get image height and convert to int/double 
                            std::stringstream geek(s);
                            geek >> h;
                        }
                        counter++;
                    }

                   //printf("image %d, %d\n", w, h);

                }
                else if (found == 'P')
                {
                    //pixel size
                    std::istringstream buf(line);
                    std::istream_iterator<std::string> beg(buf), end;
                    //Add item to vector variable for ease of acess 
                    std::vector<std::string> tokens(beg, end);
                    int counter = 0;
                    //for loop to set data from stirng value.
                    for (auto& s : tokens) {
                        if (counter == 1) {
                            //get image pixel size and convert to int/double 
                            std::stringstream geek(s);
                            geek >> location;
                        }
                        counter++;
                    }

                }
                else if (found == 'M') {
                    //max number of bounces
                    std::istringstream buf(line);
                    std::istream_iterator<std::string> beg(buf), end;
                    //Add item to vector variable for ease of acess 
                    std::vector<std::string> tokens(beg, end);
                    int counter = 0;
                    //for loop to set data from stirng value.
                    for (auto& s : tokens) {
                        if (counter == 1) {
                            //get max number of bounces and convert to int/double
                            std::stringstream geek(s);
                            geek >> m;
                        }
                        counter++;
                    }
                }
                else if (found == 'B')
                {
                    // create ray_dir a_a_box
                    //temp variables to hold values need to create ray_dir a_a_box
                    double minx, miny, minz, maxx, maxy, maxz;
                    double red, green, blue;
                    double refractive = 100;
                    std::istringstream buf(line);
                    std::istream_iterator<std::string> beg(buf), end;
                    //Add item to vector variable for ease of acess 
                    std::vector<std::string> tokens(beg, end);
                    int counter = 0;
                    //for loop to set data from stirng value.
                    for (auto& s : tokens) {
                        if (counter == 1) {
                            //get minx and convert to int/double
                            std::stringstream geek(s);
                            geek >> minx;
                        } else if (counter == 2) {
                            //get miny and convert to int/double
                            std::stringstream geek(s);
                            geek >> miny;
                        } else if (counter == 3) {
                            //get minz and convert to int/double
                            std::stringstream geek(s);
                            geek >> minz;
                        } else if (counter == 4) {
                            //get maxx and convert to int/double
                            std::stringstream geek(s);
                            geek >> maxx;
                        } else if (counter == 5) {
                            //get maxy and convert to int/double
                            std::stringstream geek(s);
                            geek >> maxy;
                        } else if (counter == 6) {
                            //get maxz and convert to int/double
                            std::stringstream geek(s);
                            geek >> maxz;
                        } else if (counter == 7) {
                            //get red and convert to int/double
                            std::stringstream geek(s);
                            geek >> red;
                        } else if (counter == 8) {
                            //get green and convert to int/double
                            std::stringstream geek(s);
                            geek >> green;
                        } else if (counter == 9) {
                            //get blue and convert to int/double
                            std::stringstream geek(s);
                            geek >> blue;
                        } else if (counter == 10) {
                            //get refractive index and convert to int/double
                            std::stringstream geek(s);
                            geek >> refractive;
                        }
                        counter++;
                    }


                    //printf("%f %f  %f  %f   %f   %f   %f   %f   %f ", minx, miny, minz, maxx, maxy, maxz, red, green, blue);
                    //creat the color using red, grean and blue, also set the refractive
                    auto box_color = make_shared<m_type>(color(red, green, blue), 0);
                    //add the object to the world with the given coordinates
                    world.add(make_shared<a_a_box>(point(minx, miny, minz), point(maxx, maxy, maxz), box_color));




                }
                else if (found == 'S')
                {
                    //craete ray_dir sphere
                    //temp varirables to create sphere
                    double centerx, centery, centerz, radius;
                    double red, green, blue;
                    double refractive = 100;

                    std::istringstream buf(line);
                    std::istream_iterator<std::string> beg(buf), end;
                    //Add item to vector variable for ease of acess 
                    std::vector<std::string> tokens(beg, end);
                    int counter = 0;
                    //for loop to set data from stirng value.
                    for (auto& s : tokens) {
                        if (counter == 1) {
                            //get centerx and convert to int/double
                            std::stringstream geek(s);
                            geek >> centerx;
                        } else if (counter == 2) {
                            //get centery and convert to int/double
                            std::stringstream geek(s);
                            geek >> centery;
                        } else if (counter == 3) {
                            //get centerz blue and convert to int/double
                            std::stringstream geek(s);
                            geek >> centerz;
                        } else if (counter == 4) {
                            //get radius and convert to int/double
                            std::stringstream geek(s);
                            geek >> radius;
                        } else if (counter == 5) {
                            //get red and convert to int/double
                            std::stringstream geek(s);
                            geek >> red;
                        } else if (counter == 6) {
                            //get Green and convert to int/double
                            std::stringstream geek(s);
                            geek >> green;
                        } else if (counter == 7) {
                            //get blue and convert to int/double
                            std::stringstream geek(s);
                            geek >> blue;
                        } else if (counter == 8) {
                            //get refractive index and convert to int/double
                            std::stringstream geek(s);
                            geek >> refractive;
                        }
                        counter++;
                    }


                    //printf("%f %f  %f  %f   %f   %f   %f   %f   %f ", minx, miny, minz, maxx, maxy, maxz, red, green, blue);
                    //creat the color using red, grean and blue, also set the refractive
                    auto sphere_color = make_shared<m_type>(color(red, green, blue), 0);
                    //add the object to the world with the given coordinates
                    world.add(make_shared<sphere>(point(centerx, centery, centerz), radius, sphere_color));

                }
                else if (found == 'L')
                {
                    //L <light xpos> <light ypos> <light zpos> <red> <grn> <blue>
                    //light case, not impleneted

                }
                else
                {
                    //std::cout << "No letters on this line \n";
                }
            }
        }
    }
    else
    {
        std::cerr << "Unable to open file \n";
    }

    //exit(0);
 
    // create the camera object
    point location_camera(0,0,0);
    point look_direction(0,0,-1);
    vector default_view(0,1,0);
    auto focus = 0.0005;
   

    //get data for the camera 
    const double aspect = w / h;
    const int image_width = w;
    const int image_height = h;
    const int samples = 10;
    const int max_depth = m;
    const double feild_of_view = atan(aspect) * 2 * 180 / 3.1415;
    color background(0, 0, 0);
    //create camera object
    camera cam(location_camera, look_direction, default_view, feild_of_view, aspect, focus);

    // Render
    //change the output stream ot and out.ppm file
    std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
    //header of ppm files
    ofs << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    //double for loop to render all pixels within the image.
    for (int x = image_height-1; x >= 0; --x) {
        //countdown for user
        std::cerr << "\rRendering the image: " << x << ' ' << std::flush;
        for (int y = 0; y < image_width; ++y) {
            //temp variable for pixel
            color current_pixel(0,0,0);
            for (int s = 0; s < samples; ++s) {
                //used for defraction in the image, getting random data for where ray will go.
                auto tester = (y + r_double()) / (image_width-1);
                auto tester2 = (x + r_double()) / (image_height-1);
                //get the ray
                ray r = cam.get_ray(tester, tester2);
                //update current pixel
                current_pixel += get_color_ray(r,  world, max_depth);
            }
            //write pixel to output file
            output_ppm(samples, ofs, current_pixel);
        }
    }

    std::cerr << "\nDone.\n";
}
