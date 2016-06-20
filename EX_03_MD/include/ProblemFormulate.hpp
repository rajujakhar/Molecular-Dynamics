#ifndef ProblemFormulate_HPP
#define ProblemFormulate_HPP


#include<vector>
#include<string>

#include "Particle.hpp"
#include "Type.hpp"


class ProblemFormulate
{
	public:
		void readDataInputFile(const std::string&);// to read the particle input data
     		void simulateMD();// to drive the molecular dynamics simulation
		void updateForce();// to update the force vector for each particle
		void updatePosition(); // to update the position vector for each particle
		void saveOldForce();// to copy f_ to fOld_	
		void resetForces();// to reset the f_ vector to 0 after each iteration		
		void updateVelocity(); // to update the velocity vector for each particle
		real calcDistance(const size_t&, const size_t&);// to find distance between two particles
		void applyPeriodicBC(const size_t&);//to apply periodic BCs on particle i		
		void writeOutput(const size_t);// to write the output for i-th timestep


        ProblemFormulate(){;}


	private:
		std::vector<Particle> molecules_;
	/*	
		std::string title = "blocks"; // tmp variable for name
		size_t out_freq = 1600;// tmp variable for vis_space
		real t_o = 0.0; //tmp variable for t_start
		real t_n = 8.0; //tmp variable for t_end
		real del_t = 0.00005; //tmp variable for delta_t
		real x_min = 0.0; //tmp variable for x_min
		real y_min = 0.0; //tmp variable for y_min
		real z_min = 0.0; //tmp variable for z_min
		real x_max = 224.4924; //tmp variable for x_max
		real y_max = 224.4924; //tmp variable for y_max
		real z_max = 224.4924; //tmp variable for z_max
		real r_c = 2.5; //tmp variable for r_cut
		real e = 5.0; //tmp variable for epsilon
		real s = 1.0;// tmp variable for sigma
  	*/
	
		std::string title = "mini"; // tmp variable for name
		size_t out_freq = 1;// tmp variable for vis_space
		real t_o = 0.0; //tmp variable for t_start
		real t_n = 0.001; //tmp variable for t_end
		real del_t = 0.0001; //tmp variable for delta_t
		real x_min = 0.0; //tmp variable for x_min
		real y_min = 0.0; //tmp variable for y_min
		real z_min = 0.0; //tmp variable for z_min
		real x_max = 200; //tmp variable for x_max
		real y_max = 200; //tmp variable for y_max
		real z_max = 200; //tmp variable for z_max
		real r_c = 2.5; //tmp variable for r_cut
		real e = 5.0e-5; //tmp variable for epsilon
		real s = 1.0;// tmp variable for sigma

	
};


#endif