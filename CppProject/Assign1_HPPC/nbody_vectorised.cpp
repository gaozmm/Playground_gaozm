#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <H5Cpp.h>
#include <chrono>
#include <argparse.hpp>

// SI units: Newtons constant of gravity, Solar Mass, distance from Earth to Sun (AU), and JPL unit of velocity, AU/day
const double G = 6.673e-11;
const double SOLAR_MASS = 1.98892e30;
const double AU = 149597870700;
const double AD = AU / (24. * 3600.);
int num_of_asteroids = 1e5;

/** Class that represent bodies, which can be a sun, a planet, or an asteroid */
class Bodies {
public:
    std::vector<double> mass;
    // The position of the bodies in the x, y, and z axis
    std::vector<double> pos_x;
    std::vector<double> pos_y;
    std::vector<double> pos_z;
    // The velocity of the bodies in the x, y, and z axis
    std::vector<double> vel_x;
    std::vector<double> vel_y;
    std::vector<double> vel_z;
    explicit Bodies(int Nx) : mass(Nx), pos_x(Nx), pos_y(Nx), pos_z(Nx), vel_x(Nx), vel_y(Nx), vel_z(Nx) {}
    // The mass of the bodies

#pragma acc routine seq
    void update(int &nv, double &m, double &x, double &y, double &z, double &vx, double &vy, double &vz){
        {
            mass[nv] = m;
            pos_x[nv] = x;
            pos_y[nv] = y;
            pos_z[nv] = z;
            vel_x[nv] = vx;
            vel_y[nv] = vy;
            vel_z[nv] = vz;
        };
    };
};

/** Class that represent a solar system, which consist of a sun, some planets, and many asteroids. */
class SolarSystem {
public:
    Bodies asteroids;
    Bodies sun_and_planets;
    SolarSystem(int noa): asteroids(noa), sun_and_planets(noa) {}
    // The first Body is the sun and the rest are planets
};

/** Function that returns the Kepler velocity -- corresponding to a body in a circular orbit */
double kepler_velocity(const double &pos_x, const double &pos_y, const double &pos_z) {
    double r = std::sqrt(pos_x * pos_x + pos_y * pos_y + pos_z * pos_z);
    return std::sqrt(G * SOLAR_MASS / r);
}

/** Return a new random solar system
 *
 * @param num_of_asteroids The number of asteroids
 * @param seed             The random seed the use
 * @return                 The new solar system
 */
// #pragma acc routine seq
SolarSystem random_system(int num_of_asteroids, int seed) {
    auto begin = std::chrono::steady_clock::now();
    SolarSystem solar_system(num_of_asteroids);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> random_real(0, 1);
    int n = 0;
    // Individually update the celestial body into solar_system
    // Begin with the solar
    double mass = SOLAR_MASS;
    double pos_x = -6.939120887892797E-3*AU; double pos_y =  5.747188451007198E-3*AU;
    double pos_z =  1.148384057799366E-4*AU;
    double vel_x = -6.544433519666372E-06*AD; double vel_y = -6.172720711215085E-06*AD;
    double vel_z =  2.101452452998284E-07*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    //Mercury
    n++;
    mass = SOLAR_MASS / 6023682.;
    pos_x = -3.489435534060263E-01*AU; pos_y =  1.198095410068464E-01*AU; pos_z =  4.080789837136244E-02*AU;
    vel_x = -1.472289494295207E-02*AD; vel_y = -2.549472074554867E-02*AD; vel_z = -7.326786621442656E-04*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Venus
    n++;
    mass = SOLAR_MASS / 408523.72;
    pos_x =  3.581712636091036E-01*AU; pos_y = -6.235356522141564E-01*AU; pos_z = -2.958988935951608E-02*AU;
    vel_x =  1.735247476983314E-02*AD; vel_y =  1.007432834238533E-02*AD; vel_z = -8.631657907430781E-04*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Earth
    mass = SOLAR_MASS / 332946.0487;
    n++;
    pos_x = -8.078442505231266E-01*AU; pos_y =  5.831320274468754E-01*AU; pos_z =  9.273772511744191E-05*AU;
    vel_x = -1.035103163864096E-02*AD; vel_y = -1.403326698993734E-02*AD; vel_z =  1.019838772702252E-06*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Mars
    n++;
    mass = SOLAR_MASS / 3098703.59;
    pos_x =  5.433407676778666E-02*AU; pos_y =  1.568439781218769E+00*AU; pos_z =  3.135948622597734E-02*AU;
    vel_x = -1.346015687592686E-02*AD; vel_y =  1.729635176008936E-03*AD; vel_z =  3.666138547700492E-04*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Jupiter
    n++;
    mass = SOLAR_MASS / 1047.56549688;
    pos_x =  3.284755554859510E+00*AU; pos_y = -3.864876891397587E+00*AU; pos_z =  -5.745522202367392E-02*AU;
    vel_x =  5.656797807092814E-03*AD; vel_y =  5.243733492952428E-03*AD; vel_z = -1.483174252349262E-04*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Saturn
    n++;
    mass = SOLAR_MASS / 3498.76667474;
    pos_x =  5.669212821032890E+00*AU; pos_y = -8.202584462326351E+00*AU; pos_z =  -8.307987098803682E-02*AU;
    vel_x =  4.278650561361220E-03*AD; vel_y =  3.157963341369132E-03*AD; vel_z = -2.250710871565522E-04*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Uranus
    n++;
    mass = SOLAR_MASS / 22905.3426;
    pos_x =  1.523494497890529E+01*AU; pos_y =  1.259053595876394E+01*AU; pos_z =  -1.506088277844468E-01*AU;
    vel_x = -2.534304435254144E-03*AD; vel_y =  2.848473830050199E-03*AD; vel_z =   4.326056025806510E-05*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    // Neptune
    n++;
    mass = SOLAR_MASS / 19416.3129471;
    pos_x =  2.947625865563749E+01*AU; pos_y = -5.092058884278901E+00*AU; pos_z = -5.744488884613864E-01*AU;
    vel_x =  5.139497530346715E-04*AD; vel_y =  3.112318538936126E-03*AD; vel_z = -7.589746140797978E-05*AD;
    solar_system.sun_and_planets.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);

    // Place asteroids in a belt from 2 to 3.5 AU
//     #pragma acc routine seq
//     std::cout << solar_system.asteroids.mass.size() << std::endl;
//     #pragma acc parallel loop
    for (n = 0; n < num_of_asteroids; ++n) {
        double dist = AU * (2. + 1.5 * random_real(gen));
        double theta = M_PI / 2 * random_real(gen);
        pos_x = dist * std::cos(theta);
        pos_y = dist * std::sin(theta);
        pos_z = (random_real(gen) - 0.5) * .1 * dist;

        double v_kep = kepler_velocity(pos_x, pos_y, pos_z);
        vel_x = -1 * std::sin(theta) * v_kep;
        vel_y = std::cos(theta) * v_kep;
        vel_z = 0;
        mass = random_real(gen) * 1e21 + 1e14;
        solar_system.asteroids.update(n, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z);
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "random time: " << (end - begin).count() / 1000000000.0 << " sec"<< std::endl;
    std::cout << "Performance: " << (end - begin).count() / ((9 + num_of_asteroids)) << " nanosecs / particle update"<< std::endl;
    return solar_system;
}

/** Update the velocity of `a` based on `b`
 *
 * @param a  The body to update
 * @param b  The body which act on `a`
 * @param dt The time step size
 */
void update_velocity(Bodies &a, const Bodies &b, const int &i, const double &dt) {
    double vel_x;
    double vel_y;
    double vel_z;
#pragma acc parallel reduction(+:vel_x,vel_y,vel_z)
    {
#pragma acc loop
        for (int k = 0; k < i; k++) {
            vel_x = 0;
            vel_y = 0;
            vel_z = 0;
#pragma acc loop reduction(+:vel_x,vel_y,vel_z)
            for (int h = 0; h < 9; h++){
                double xdiff = b.pos_x[h] - a.pos_x[k];
                double ydiff = b.pos_y[h] - a.pos_y[k];
                double zdiff = b.pos_z[h] - a.pos_z[k];
                double r = std::sqrt(xdiff * xdiff + ydiff * ydiff + zdiff * zdiff);
                double temp = (b.mass[h] / (r * r * r)) * dt;
                vel_x += temp * (xdiff);
                vel_y += temp * (ydiff);
                vel_z += temp * (zdiff);
            }
            a.vel_x[k] = a.vel_x[k] + vel_x;
            a.vel_y[k] = a.vel_y[k] + vel_y;
            a.vel_z[k] = a.vel_z[k] + vel_z;
        }
    };
}
/** Kick a set of bodies forward in time due to their mutual gravitational interaction
 *
 * @param a  The bodies to update
 * @param dt The time step size
 */
void kick_same(Bodies &a, const double &dt) {
    int n = 9;
    update_velocity(a, a, n, dt);
}

/** Kick a set of bodies forward in time due to gravitational interaction with another set of bodies
 *
 * @param a  The bodies to update
 * @param b  The bodies that perturb
 * @param dt The time step size
 */
void kick_other(Bodies &a, const Bodies &b, const double &dt) {
    int n_a = a.pos_x.size();
    update_velocity(a, b, n_a, dt);
}

/** Drift a set of bodies forward in time
 *
 * @param bodies The bodies to update
 * @param dt     The time step size
 */
void drift(Bodies &bodies, const double &dt) {
//     #pragma acc enter data copyout(bodies)
//     #pragma acc data present(bodies.pos_x[0:n]) present(bodies.pos_x[0:n]) present(bodies.pos_z[0:n])
//     #pragma acc data present(bodies.vel_x[0:n]) present(bodies.vel_x[0:n]) present(bodies.vel_z[0:n])
#pragma acc parallel
    {
#pragma acc loop
        for (int i = 0; i < num_of_asteroids; i++) {
            bodies.pos_x[i] += bodies.vel_x[i] * dt;
            bodies.pos_y[i] += bodies.vel_y[i] * dt;
            bodies.pos_z[i] += bodies.vel_z[i] * dt;
        }
    }
}

/** Integrate one time step of the solar system
 *
 * @param solar_system  The solar system to update
 * @param dt            The time step size
 */
void integrate(SolarSystem &solar_system, double dt) {
    // Kick is done twice --> only half dt each time
    double const Ghdt = G*0.5*dt;
    // First kick
    // Update velocity of all bodies
        kick_same(solar_system.sun_and_planets, Ghdt);
        kick_other(solar_system.asteroids,solar_system.sun_and_planets, Ghdt);
    // Drift: Update position of all bodies
        drift(solar_system.sun_and_planets, dt);
        drift(solar_system.asteroids, dt);
    // Second kick
    // Update velocity of all bodies
        kick_same(solar_system.sun_and_planets, Ghdt);
        kick_other(solar_system.asteroids,solar_system.sun_and_planets, Ghdt);
}

/** Write data to a hdf5 file
 *
 * @param group  The hdf5 group to write in
 * @param name   The name of the data
 * @param shape  The shape of the data
 * @param data   The data
 */
void write_hdf5(H5::Group &group, const std::string &name, const std::vector<hsize_t> &shape,
                const std::vector<double> &data) {

    H5::DataSpace dataspace(static_cast<int>(shape.size()), &shape[0]);
    H5::DataSet dataset = group.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
}

/** Write the solar system to a hdf5 file (use `visual.py` to visualize the hdf5 data)
 *
 * @param solar_systems  The solar system to write
 * @param filename       The filename to write to
 */
void write_hdf5(const std::vector<SolarSystem> &solar_systems, const std::string &filename) {

    H5::H5File file(filename, H5F_ACC_TRUNC);

    for (int i = 0; i < solar_systems.size(); ++i) {
        H5::Group group(file.createGroup("/" + std::to_string(i)));
        {
            const std::vector<double> &data = solar_systems[i].sun_and_planets.mass;
            int n = data.size();
            write_hdf5(group, "sun_and_planets_mass", {n}, data);
        }
        {
            std::vector<double> data;
            const Bodies &bodies = solar_systems[i].sun_and_planets;
            int n = bodies.mass.size();
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_x[i]);
            }
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_y[i]);
            }
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_z[i]);
            }
            write_hdf5(group, "sun_and_planets_position", {3, n}, data);
        }
        {
            std::vector<double> data;
            const Bodies &bodies = solar_systems[i].asteroids;
            int n = bodies.mass.size();
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_x[i]);
            }
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_y[i]);
            }
            for (int i = 0; i < n; ++i) {
                data.push_back(bodies.pos_z[i]);
            }
            write_hdf5(group, "asteroids_position", {3, n}, data);
        }

    }
}

/** N-body Solar System simulation
 *
 * @param num_of_iterations  Number of iterations
 * @param num_of_asteroids   Number of asteroids
 * @param seed               Random seed
 * @param filename           Filename to write each time step to. If empty, do data is written.
 */
void simulate(int num_of_iterations, int seed,
              int every, const std::string &filename) {
    double dt = 1e5;
    SolarSystem system = random_system(num_of_asteroids, seed);
    std::vector<SolarSystem> systems;
    systems.push_back(system);
    auto begin = std::chrono::steady_clock::now();
    for (int i = 0; i < num_of_iterations; ++i) {
        integrate(system, dt);
        if (!filename.empty() && (i % every) == 0) {
            systems.push_back(system);
        }
    }
     std::cout << systems.size() << std::endl;
     if (!filename.empty()) {
         write_hdf5(systems, filename);
     }
    auto end = std::chrono::steady_clock::now();
    std::cout << "elapsed time: " << (end - begin).count() / 1000000000.0 << " sec"<< std::endl;
    std::cout << "Performance: " << (end - begin).count() /
                                    (num_of_iterations * (1 + 8 + num_of_asteroids)) << " nanosecs / particle update"<< std::endl;
}

/** Main function that parses the command line and start the simulation */
int main(int argc, char **argv) {
    util::ArgParser args(argc, argv);
    int iterations, every, seed;
    if (args.cmdOptionExists("--iter")) {
        iterations = std::stoi(args.getCmdOption("--iter"));
        if (iterations < 0) {
            throw std::invalid_argument("iter most be a positive integer (e.g. --iter 100)");
        }
    } else {
        throw std::invalid_argument("You must specify the number of iterations (e.g. --iter 100)");
    }
    if (args.cmdOptionExists("--every")) {
        every = std::stoi(args.getCmdOption("--every"));
        if (every < 0) {
            throw std::invalid_argument("every most be a positive integer (e.g. --every 10)");
        }
    } else {
        every = 10;
    }
    if (args.cmdOptionExists("--asteroids")) {
        num_of_asteroids = std::stoi(args.getCmdOption("--asteroids"));
        if (num_of_asteroids < 0) {
            throw std::invalid_argument("asteroids most be a positive integer (e.g. --asteroids 1000)");
        }
    } else {
        throw std::invalid_argument("You must specify the number of asteroids (e.g. --asteroids 1000)");
    }
    if (args.cmdOptionExists("--seed")) {
        seed = std::stoi(args.getCmdOption("--seed"));
        if (seed < 0) {
            throw std::invalid_argument("Seed most be a positive integer (e.g. --seed 42)");
        }
    } else {
        seed = std::random_device{}(); // Default seed is taken from hardware
    }
    const std::string &filename = args.getCmdOption("--out");
    if (!filename.empty()) {
        std::cout << "Writing data to " << filename << std::endl;
        std::cout << "Writing out every " << every << "th iteration to disk (--every "<< every << ")" << std::endl;
    }

    simulate(static_cast<int>(iterations),
             static_cast<int>(seed),
             static_cast<int>(every),
             filename);
    return 0;
}
