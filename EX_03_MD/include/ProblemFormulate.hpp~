#ifndef ProblemFormulate_HPP
#define ProblemFormulate_HPP


#include<vector>
#include<string>
#include<list>
#include<array>

#include "Particle.hpp"
#include "Type.hpp"


class ProblemFormulate
{
	public:
		void readDataInputFile(const std::string&);// to read the particle input data
     		void simulateMD();// to drive the molecular dynamics simulation
		void updateForce();// to update the force vector for each particle
		void updatePosition(); // to update the position vector for each particle
		void updateCells(const size_t&, const size_t&, const size_t&);//to update cells_ if particle migrates
		void saveOldForce();// to copy f_ to fOld_	
		void resetForces();// to reset the f_ vector to 0 after each iteration		
		void updateVelocity(); // to update the velocity vector for each particle
		real calcDistance(const size_t&, const size_t&);// to find distance between two particles
		void applyPeriodicBC(const size_t&);//to apply periodic BCs on particle i		
		void writeOutput(const size_t);// to write the output for i-th timestep
		void setUpLinkedCells();//to setup the cells 
		
		size_t calcCellPos(const size_t&, const size_t&, const size_t&);//to find the cell position from the coordinates
		void calcCellCoords(const size_t&, std::array<size_t,3>& );//to find cell coords from position and updates in the array passed by reference
		size_t determineParticleLoc(const size_t&);//to determine cell in which particle k resides 
		//void clearCells();//to clear the cells_			
		void determineNeighbourSet(const size_t&, std::vector<size_t>&);//to find the neighbours of a cell		

	private:
		std::vector<Particle> molecules_;
		std::vector<std::list<size_t>> cells_;//contains the ID of particles in a cell
		/*cells are numbered first along -x, then along -y and then along -z direction*/
		
		size_t n_1_c, n_2_c, n_3_c; // number of cells in each of the 3-directions (-x, -y, -z)
		size_t n_tot_c; //total numbe rof cells (= n_1_c* n_2_c* n_3_c)
		real r_1_c, r_2_c, r_3_c; // modified values of r_cut in each of the 3-directions
		


		
		std::string name = "blocks";
		size_t vis_space = 1600;
		real t_start = 0.0;
		real t_end = 8.0; 
		real delta_t = 0.00005;
		real x_min = 0.0;
		real y_min = 0.0;
		real z_min = 0.0;
		real x_max = 224.4924;
		real y_max = 224.4924;
		real z_max = 224.4924;
		real r_cut = 2.5;
		real e = 5.0; //tmp variable for epsilon
		real s = 1.0; // tmp variable for sigma
  	
	/*
		std::string name = "mini";
		size_t vis_space = 1;
		real t_start = 0.0; 
		real t_end = 0.001; 
		real delta_t = 0.0001;
		real x_min = 0.0;
		real y_min = 0.0;
		real z_min = 0.0;
		real x_max = 200;
		real y_max = 200;
		real z_max = 200;
		real r_cut = 2.5;
		real e = 5.0e-5; //tmp variable for epsilon
		real s = 1.0;// tmp variable for sigma
	*/
	
};


#endif
