/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;
	default_random_engine gen;

	normal_distribution<double> noise_x(x, std[0]);
	normal_distribution<double> noise_y(y, std[1]);
	normal_distribution<double> noise_theta(theta, std[2]);

	for (unsigned int i = 0; i < num_particles; i++) {
		Particle particle;
		// Add Gaussian Noise
		particle.x = noise_x(gen);
		particle.y = noise_y(gen);
		particle.theta = noise_theta(gen);

		particel.weight = 1.0;

		particles.push_back(particle);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// ************ PREDICTION FORMULAS ******************************
	// xf = x_0 + v/theta_dot[sin(theta + theta_dot(dt)) - sin(theta)]
	// yf = y_0 + v/theta_dot[cos(theta) - cos(theta + theta_dot(dt))]
	// thetaf = theta_0 + theta_dot(dt)
	// ***************************************************************

	default_random_engine gen;
	normal_distribution<double> noise_x(0.0, std_pos[0]);
	normal_distribution<double> noise_y(0.0, std_pos[1]);
	normal_distribution<double> noise_theta(0.0, std_pos[2]);

	for (unsigned int i = 0; i < num_particles; i++) {
		particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
		particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
		particles[i].theta += yaw_rate*delta_t;

		//Add Random Gaussian Noise now
		particles[i].x += noise_x(gen);
		particles[i].y += noise_y(gen);
		particles[i].theta += noise_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < observations.size(); i++) {
		double current_closest = numeric_limits<double>::infinity();
		for (unsigned int j = 0; j < predicted.size(); j++) {
			double diff_x = predicted[j].x - observations[i].x;
			double diff_y = predicted[j].y - observations[i].y;
			double distance = sqrt((diff_x*diff_x) + (diff_y*diff_y));

			int assigned_id = -1;
			if (distance < current_closest) {
				assigned_id = predicted[j].id;
				current_closest = distance;
			}

			observations[i].id = assigned_id;
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// transform car observation coordinates to map coordinates
	std::vector<LandmarkObs> observations_map_coordinates;
	for (unsigned int i = 0; i < num_particles; i++) {
		for (unsigned int j = 0; j < observations.size(); j++) {
			LandmarkObs transformed_landmark;
			transformed_landmark.x = particles[i].x + (cos(particles[i].theta * observations[j].x) - (sin(particles[i].theta * observations[j].y)));
			transformed_landmark.y = particles[i].y + (sin(particles[i].theta * observations[j].x) + (cos(particles[i].theta * observations[j].y)));
			transformed_landmark.id = observations[j].id;

			observations_map_coordinates.push_back(transformed_landmark);
		}

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
