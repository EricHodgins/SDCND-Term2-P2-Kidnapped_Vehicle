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

	for (int i = 0; i < num_particles; i++) {
		Particle particle;
		// Add Gaussian Noise
		particle.x = noise_x(gen);
		particle.y = noise_y(gen);
		particle.theta = noise_theta(gen);

		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(particle.weight);
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

	for (int i = 0; i < num_particles; i++) {
		if (fabs(yaw_rate) < 0.0001) {
			particles[i].x += velocity * cos(particles[i].theta) * delta_t;
			particles[i].y += velocity * sin(particles[i].theta) * delta_t;
		} else {
			particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}

		//Add Random Gaussian Noise now
		normal_distribution<double> noise_x(particles[i].x, std_pos[0]);
		normal_distribution<double> noise_y(particles[i].y, std_pos[1]);
		normal_distribution<double> noise_theta(particles[i].theta, std_pos[2]);

		particles[i].x = noise_x(gen);
		particles[i].y = noise_y(gen);
		particles[i].theta = noise_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (size_t i = 0; i < observations.size(); i++) {
		double current_closest = numeric_limits<double>::infinity();
		int assigned_id = -1;
		for (size_t j = 0; j < predicted.size(); j++) {
			double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);//sqrt((diff_x*diff_x) + (diff_y*diff_y));


			if (distance < current_closest) {
				assigned_id = predicted[j].id;
				current_closest = distance;
			}
		}

		observations[i].id = assigned_id;
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

	// Steps:
	// 1. Transform coordinates to map coordinates
	// 2. Associate landmarks
	// 3. Update Weights
	// 		- Determine measurement probability
	// 		- Combine probabilities


	for (int i = 0; i < num_particles; i++) {
		// 1. Transform coordinates
		std::vector<LandmarkObs> transformed_coordinates;

		for (size_t j = 0; j < observations.size(); j++) {
			LandmarkObs transformed_landmark;
			transformed_landmark.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
			transformed_landmark.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);
			transformed_landmark.id = observations[j].id;

			transformed_coordinates.push_back(transformed_landmark);
		}

		// Find landmarks that are only within the sensor range
		vector<LandmarkObs> sensed_landmarks;
		for (size_t s = 0; s < map_landmarks.landmark_list.size(); s++) {
			Map::single_landmark_s current_landmark = map_landmarks.landmark_list[s];
			// Check if within x range
			if (fabs(particles[i].x - current_landmark.x_f)) {
				// Check if within y range
				if (fabs(particles[i].y - current_landmark.y_f)) {
					LandmarkObs sensed_landmark;
					sensed_landmark.id = current_landmark.id_i;
					sensed_landmark.x = current_landmark.x_f;
					sensed_landmark.y = current_landmark.y_f;
					sensed_landmarks.push_back(sensed_landmark);
				}
			}
		}

		// 2. Associate Landmarks
		dataAssociation(sensed_landmarks, transformed_coordinates);

		// 3. Update Weights
		particles[i].weight = 1.0;
		double std_x = std_landmark[0];
		double std_y = std_landmark[1];
		double std_x_2 = 2 * (std_x*std_x);
		double std_y_2 = 2 * (std_y*std_y);
		double noramlizer = 1 / (2*M_PI*std_x*std_y);
		for (size_t k = 0; k < transformed_coordinates.size(); k++) {
			for (size_t l = 0; l < sensed_landmarks.size(); l++) {
				if (transformed_coordinates[k].id == sensed_landmarks[l].id) {
					// Mulitvariate-Gaussian Formula
					double x = transformed_coordinates[k].x;
					double y = transformed_coordinates[k].y;
					double mean_x = sensed_landmarks[l].x;
					double mean_y = sensed_landmarks[l].y;

					double exponent = -(((x - mean_x)*(x - mean_x)/(std_x_2)) + ((y - mean_y)*(y - mean_y)/(std_y_2)));
					double e = exp(exponent);

					double prob = noramlizer * e;

					particles[i].weight *= prob;
					weights[i] = particles[i].weight;
				}
			}
		}

	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::vector<Particle> resampled_particles;
	default_random_engine gen;

	uniform_int_distribution<int> random_int(0, num_particles-1);

	double beta = 0.0;
	double max_w = 2 * (*max_element(weights.begin(), weights.end())); // max_element returns an iterator
	uniform_real_distribution<double> random_double(0, max_w);

	for (int i = 0; i < num_particles; i++) {
		int idx_random = random_int(gen); // -1 ?
		beta = beta + random_double(gen);
		while (weights[idx_random] < beta) {
			idx_random += 1;
			idx_random = idx_random % num_particles;
			beta = beta - weights[idx_random];
		}

		resampled_particles.push_back(particles[idx_random]);
	}

	particles = resampled_particles;
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

    return particle;
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
