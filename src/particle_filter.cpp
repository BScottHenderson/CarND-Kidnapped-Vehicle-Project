/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 *
 *  Modified on: Nov 19, 2018
 *      Author: Scott Henderson
 */

#include <random>
#include <algorithm>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <iterator>
#include <functional>
#include <string>
#include <map>
#include <time.h>

#include "particle_filter.h"

#define M 1000

static std::random_device         rd;
static std::default_random_engine gen(rd());

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  if (is_initialized)
    return;

  // Set the number of particles
  num_particles = M;
  particles = std::vector<Particle>(num_particles);
  weights   = std::vector<double>(num_particles);   // Weights for each particle. Stored redundantly here to make
                                                    // calculations easier.

  // Create normal distributions for x, y, theta.
  std::normal_distribution<double>  dist_x(x, std[0]);
  std::normal_distribution<double>  dist_y(y, std[1]);
  std::normal_distribution<double>  dist_theta(theta, std[2]);

  // For each particle ...
  int id = 0;
  for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p) {
    // Just assign a sequential, 0-based id value.
    p->id = id++;

    // Set the position (x, y) and heading using normal distributions based on the
    // given initial position and heading and standard deviation values.
    p->x = dist_x(gen);
    p->y = dist_y(gen);
    p->theta = dist_theta(gen);

    // Make sure other values are initialized.
    p->weight = 1.0;
    //p->associations.empty();
    //p->sense_x.empty();
    //p->sense_y.empty();
  }

  // Set all weights to 1.0
  std::fill(weights.begin(), weights.end(), 1.0);

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // Create normal distributions for sensor noise.
  std::normal_distribution<double>  dist_x(0, std_pos[0]);
  std::normal_distribution<double>  dist_y(0, std_pos[1]);
  std::normal_distribution<double>  dist_theta(0, std_pos[2]);

  // Cache values not dependent on particle data.
  double  c = fabs(yaw_rate) < 0.0001 ? 0.0 : velocity / yaw_rate;
  double  yaw = yaw_rate * delta_t;

  // For each particle ...
  for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p) {
    // Update the particle position and heading.
    if (fabs(yaw_rate) < 0.0001) {  // Avoid divide-by-zero.
      // Not turning so just add distance traveled at the current heading.
      double  dist = velocity * delta_t;
      p->x += dist * cos(p->theta);
      p->y += dist * sin(p->theta);
      // yaw_rate is 0 so don't update p->theta
    }
    else {
      p->x += c * (sin(p->theta + yaw) - sin(p->theta));
      p->y += c * (cos(p->theta) -  cos(p->theta + yaw));
      p->theta += yaw;
    }

    // Add sensor noise based on measurement standard deviations.
    p->x += dist_x(gen);
    p->y += dist_y(gen);
    p->theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

  // For each sensor observation ...
  for (std::vector<LandmarkObs>::iterator o = observations.begin(); o != observations.end(); ++o) {
    // Examine each landmark to find the closest.
    int     min_distance_id = -1;
    double  min_distance = std::numeric_limits<double>::max();
    for (std::vector<LandmarkObs>::const_iterator l = predicted.begin(); l != predicted.end(); ++l) {
      double  d = dist(o->x, o->y, l->x, l->y);
      if (d < min_distance) {
        min_distance = d;
        min_distance_id = l->id;
      }
    }

    // Update the observation with the closest landmark id.
    o->id = min_distance_id;
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

#if __DEBUG
  time_t  start = clock();
  time_t  local_start;
  time_t  local_elapsed;
  time_t  local_time = 0;
#endif

  //
  // Set the weight for each particle ...
  //

  for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p) {

    //
    // Generate a list of map landmarks that are within sensor range of the vehicle/particle.
    //

#if __DEBUG
    local_start = clock();
#endif
    // Extract a subset of map landmarks to make the nearest neighbor calculation faster.
    // See the dataAssociation() method.
    std::vector<LandmarkObs> predicted;
    int id = 0;
    for (std::vector<Map::single_landmark_s>::const_iterator l = map_landmarks.landmark_list.begin();
          l != map_landmarks.landmark_list.end(); ++l) {
      double  dx = l->x_f - p->x;
      double  dy = l->y_f - p->y;
      if (dist(dx, dy) <= sensor_range) {
        //predicted.push_back({l->id_i, l->x_f, l->y_f});
        // Do not use the actual landmark id value. Since we are building a subset of landmarks
        // anyway we may as well reassign the id values here. The huge advantage of doing this is
        // that we can avoid the "search for landmark matching id" step when updating particle
        // weights (see below). Instead we can just use the id value assigned here as an index
        // into the predicted observation array.
        predicted.push_back({id++, l->x_f, l->y_f});
      }
    }
#if __DEBUG
    local_elapsed = clock() - local_start;
    local_time += local_elapsed;
    times["1-sensor_range"] += local_elapsed;
#endif

    //
    // Transform observations from vehicle/particle coordinates to map coordinates.
    //

#if __DEBUG
    local_start = clock();
#endif
    // Cache expensive calculations.
    double  cos_theta = cos(p->theta);
    double  sin_theta = sin(p->theta);

    std::vector<LandmarkObs> transformed_observations(observations.size());
    for (int o = 0; o < observations.size(); ++o) {
      double map_x = p->x + (cos_theta * observations[o].x) - (sin_theta * observations[o].y);
      double map_y = p->y + (sin_theta * observations[o].x) + (cos_theta * observations[o].y);
      transformed_observations[o] = {0, map_x, map_y};
    }
#if __DEBUG
    local_elapsed = clock() - local_start;
    local_time += local_elapsed;
    times["2-transform"] += local_elapsed;
#endif

    //
    // Nearest Neighbor
    //

#if __DEBUG
    local_start = clock();
#endif
    // Find the closest landmark for each transformed observation.
    // Note: this will update the 'id' member for each observation.
    dataAssociation(predicted, transformed_observations);
#if __DEBUG
    local_elapsed = clock() - local_start;
    local_time += local_elapsed;
    times["3-dataAssociation"] += local_elapsed;

    if (write_data) {
      std::cout << "landmarks: " << predicted.size() << std::endl;
      for (int i = 0; i < predicted.size(); ++i) {
        std::cout << "(" << predicted[i].id << ", " << predicted[i].x << ", " << predicted[i].y << ")" << std::endl;
      }
      std::cout << "observations: " << transformed_observations.size() << std::endl;
      for (int i = 0; i < transformed_observations.size(); ++i) {
        std::cout << "(" << transformed_observations[i].id << ", " << transformed_observations[i].x << ", " << transformed_observations[i].y << ")" << std::endl;
      }
    }
#endif

    //
    // Update particle observation weight based on the closest landmarks.
    //

#if __DEBUG
    local_start = clock();
#endif
    // Reset particle weight.
    p->weight = 1.0;

    // Cache a few values not dependent on observations.
    double  cxy = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
    double  cx = 2.0 * std_landmark[0] * std_landmark[0];
    double  cy = 2.0 * std_landmark[1] * std_landmark[1];

    // For each transformed observation ...
#if __DEBUG
    if (write_data)
      std::cout << "particle: " << p - particles.begin() << std::endl;
#endif
    for (int o = 0; o < transformed_observations.size(); ++o) {

      // Find the predicted observation that we've tagged as nearest neighbor.
      LandmarkObs predicted_observation =
        //*std::find_if(predicted.begin(), predicted.end(),
        //  [=] (const LandmarkObs& l) { return l.id == transformed_observations[o].id; });
        predicted[transformed_observations[o].id];

      double  dx = transformed_observations[o].x - predicted_observation.x;
      double  dy = transformed_observations[o].y - predicted_observation.y;

      double  tx = (dx * dx) / cx;
      double  ty = (dy * dy) / cy;
#if __DEBUG
      if (write_data)
        std::cout << "temp: (" << tx << ", " << ty << ")" << std::endl;
#endif

      // Update the particle weight.
      // ***************************
      // Make sure the observation weight is not zero. This happens sometimes, I'm not sure why.
      // It seems to occur when a predicted landmark is "too far away" and so one of the delta
      // position terms ends up being very large (when squared and divided by the variance)
      // and this causes the exp() function to return 0. In some interations all particles have
      // at least one of these types of observation and the particle weight array ends up being
      // all zero and this is not allowed for the discrete_distribution class (see the resample()
      // method below). The failure does not always happen at the same time step which suggests
      // that it is somehow dependent on random noise added to various values here and there but
      // I still don't really understand how the delta position values end up being so large.
      // ***************************
      double  ow = cxy * exp(-(tx + ty));
      if (ow != 0.0)
        p->weight *= ow;
#if __DEBUG
      if (write_data)
        std::cout << "obs weight: " << ow << std::endl;
      //if (ow == 0) {
      //  std::cout << "temp: (" << tx << ", " << ty << ")" << std::endl;
      //  std::cout << "tobs: (" << transformed_observations[o].x << ", " << transformed_observations[o].y << ")" << std::endl;
      //  std::cout << "pred: (" << predicted_observation.x << ", " << predicted_observation.y << ")" << std::endl;
      //}
#endif
    }
#if __DEBUG
    local_elapsed = clock() - local_start;
    local_time += local_elapsed;
    times["4-weights"] += local_elapsed;
#endif
  }

  //
  // Normalize particle weights
  //

  // Transform particle weights into a probability distribution.
  double  weight_sum =
    std::accumulate(particles.begin(), particles.end(), 0.0,
                    std::bind(std::plus<>(), std::placeholders::_1,
                      std::bind(&Particle::weight, std::placeholders::_2)));
  for (std::vector<Particle>::iterator p = particles.begin(); p != particles.end(); ++p)
    p->weight /= weight_sum;

#if __DEBUG
  time_t  elapsed = (clock() - start) - local_time;
  times["updateWeights"] += elapsed;
#endif
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Construct a discrete_distribution from all particle weights.
  for (int i = 0; i < particles.size(); ++i)
    weights[i] = particles[i].weight;
  std::discrete_distribution<>  dd(weights.begin(), weights.end());

  // Draw a new sample of particles.
  std::vector<Particle> new_particles(particles.size());
  for (int i = 0; i < particles.size(); ++i)
    new_particles[i] = particles[dd(gen)];

  // Replace the old set of particles with the new set.
  particles = new_particles;
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

template<typename T>
std::string getValues(std::vector<T> v) {
  std::stringstream ss;
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(ss, " "));
  std::string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

std::string ParticleFilter::getAssociations(Particle best)
{
  return getValues(best.associations);
}
std::string ParticleFilter::getSenseX(Particle best)
{
  return getValues(best.sense_x);
}
std::string ParticleFilter::getSenseY(Particle best)
{
  return getValues(best.sense_y);
}
