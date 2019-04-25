/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * NOTE: Initialize all the particles
   */

  if (is_initialized) {
    return;
  }
  // Set the number of particles for the particle filter
  num_particles = 100; 

  std::default_random_engine gen;
  // Create a normal (Gaussian) distribution for x, y and thet
  // around the provided GPS position and std deviation
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // Looping through the particles and adding the properties 
  for (int i=0;i< num_particles;++i){
    Particle particle;
    particle.id = i;
    // Randomly picking value from the distribution
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);  
    // Setting weight as 1
    particle.weight = 1.0; 

    particles.push_back(particle);
  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * NOTE: Add measurements to each particle and add random Gaussian noise.
   */

  std::default_random_engine gen;
  // Creates a normal distribution around zero for each of the standard deviations
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  // Looping through the particles
  for (int i=0;i< num_particles;++i){

    // Calculate the change in position if there is an existing yaw rate
    if(fabs(yaw_rate) > 0.00001){
      particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;  
    } else {
      // Calculate the change in position is there is no yaw rate
      particles[i].x += velocity*delta_t*cos(particles[i].theta);
      particles[i].y += velocity*delta_t*sin(particles[i].theta);      
    }

    // Add noise to the prediction
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen); 
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * NOTE: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   */
  
  // looping through all the observations of landmarks
  for(unsigned int i=0; i<observations.size(); ++i){
    double minDistance = 0.0;
    int mapId = -1;
    // looping through all the locations of landmarks in map coordinates
    for(unsigned int j=0; j<predicted.size(); ++j){
      // calculate the distaance between the observation and landmarks
      double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
      // identify the least distance
      if (j==0){
        minDistance = distance;
        mapId = predicted[j].id;
      } else {
        if (distance < minDistance){
          minDistance = distance;
          mapId = predicted[j].id;
        }
      }
    }
    // set the id of the landmark to the associated observation
    observations[i].id = mapId;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * NOTE: Update the weights of each particle using a mult-variate Gaussian 
   *       distribution.
   */
  
  // Loop through the particles
  for (int i=0; i< num_particles;++i){ 
    // Store each of the particles locations
    double x_particle = particles[i].x;
    double y_particle = particles[i].y;
    double theta_particle = particles[i].theta;

    vector<LandmarkObs> rangeLandmarks;
    // Loop through all of the landmarks
    for(unsigned int j=0; j < map_landmarks.landmark_list.size(); ++j){
      double lk_x = map_landmarks.landmark_list[j].x_f;
      double lk_y = map_landmarks.landmark_list[j].y_f;
      int lk_id = map_landmarks.landmark_list[j].id_i;

      // Check if this landmark is within sensor range
      double distance = dist(x_particle,y_particle,lk_x,lk_y);
      if(distance<sensor_range){
        rangeLandmarks.push_back(LandmarkObs{lk_id,lk_x,lk_y});
      }
    }

    vector<LandmarkObs> mapcoord_observations;
    // Convert observations from local co-ordinates to map co-ordinates
    for(unsigned int j=0; j < observations.size(); ++j){
      double mapcoord_x = x_particle + (cos(theta_particle)*observations[j].x) - (sin(theta_particle)*observations[j].y);
      double mapcoord_y = y_particle + (sin(theta_particle)*observations[j].x) + (cos(theta_particle)*observations[j].y);

      mapcoord_observations.push_back(LandmarkObs{observations[j].id,mapcoord_x,mapcoord_y});

    }

    // Associate the observations with the landmarks
    dataAssociation(rangeLandmarks,mapcoord_observations);

    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

    double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
    particles[i].weight = 1.0;
    
    // Loop through the observations to calculate weighted error
    for(unsigned int j=0; j<mapcoord_observations.size(); ++j){

      double x_obs = mapcoord_observations[j].x;
      double y_obs = mapcoord_observations[j].y;
      double mu_x, mu_y;
      int assoc_id = mapcoord_observations[j].id;

      for(unsigned int k=0; k < rangeLandmarks.size(); ++k){
        if(rangeLandmarks[k].id == assoc_id){
          mu_x = rangeLandmarks[k].x;
          mu_y = rangeLandmarks[k].y;
        }

      }
    
      // calculate exponent
      double exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
                  + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
      
    // calculate weight using normalization terms and exponent
      double weight = gauss_norm * exp(-exponent);
      if (weight == 0) {
         particles[i].weight *= 0.00001;
      } else {
         particles[i].weight *= weight;
      }
    }
  } 
}

void ParticleFilter::resample() {
  /**
   * NOTE: Resample particles with replacement with probability proportional 
   *   to their weight. 
   */

  double max_weight = 0.0;
  vector<double> weights;
  vector<Particle> particlesNew;

  // Identify the maximum weight
  for (int i=0; i < num_particles;++i){ 
    weights.push_back(particles[i].weight);

    if (particles[i].weight > max_weight){
      max_weight = particles[i].weight;
    }
  }

  std::default_random_engine gen;
  std::uniform_int_distribution<int> seed(0,num_particles-1);
  std::uniform_real_distribution<double> distribution(0.0,max_weight);

  int index = seed(gen);
  double beta = 0.0;

  // Wheel of weights to resample the particles
  for (int i=0; i < num_particles;++i){ 

    beta += distribution(gen)*2.0;
    while (weights[index] < beta){
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    particlesNew.push_back(particles[index]);
  }

  particles = particlesNew;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}