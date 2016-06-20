
#include "ProblemFormulate.hpp"
#include "ParameterReader.hpp"

#include<iostream>
#include<fstream> // ifstream, ofstream
#include<sstream> //ostringstream 
#include<cmath>	 //pow(x,n)

/*
This function takes the filename as an arguemeent and reads the input from the file and intializes the molecules_ data structure
*/
void ProblemFormulate::readDataInputFile(const std::string& filename)
{
	std::ifstream file_input;
	real m,x,y,z,u,v,w;

	file_input.open(filename);
	while(file_input.good())
	{
		file_input>>m;
		file_input>>x;
		file_input>>y;
		file_input>>z;
		file_input>>u;
		file_input>>v;
		file_input>>w;
		
		Particle tmpParticle(m,x,y,z,u,v,w);
		molecules_.push_back(tmpParticle);
	}
	
	std::cout<<"Finished reading the data input file, read "<<molecules_.size()<<" lines"<<std::endl;
	
//this approach leads to repetition of the last line
//The last entry must be deleted
	molecules_.erase(molecules_.end());

	std::cout<<"There are "<<molecules_.size()<<" particles in the system"<<std::endl;

//displaying the data as stored in molecules_
/*
	const size_t loopLimit = molecules_.size();

	for(size_t i =0; i<loopLimit; ++i)
	{
		molecules_[i].displayParticle(i);
	}
*/
	
}

/*
This function drives the overall simulation
*/
void ProblemFormulate::simulateMD()
{
	const size_t n_steps = (t_n-t_o)/del_t;
	std::cout<<"This simulation has "<<n_steps<<" number of timesteps"<<std::endl;	
	std:: cout<<"Starting the simulation"<<std::endl;	

//writing the intial data to the output file
	writeOutput(0);

//Updating the initial force vector
	updateForce();		

//running the simulationg for n_tsteps
	for(size_t i=1; i<=n_steps; ++i)
	{
		std::cout<<"currently processing timestep: "<<i<<std::endl;		
		updatePosition();		
		saveOldForce();
		resetForces();
		updateForce();
		updateVelocity();
		
		if(i%out_freq ==0)
			writeOutput(i);
		
	}
	std::cout<<"Simulation Complete !!!!!"<<std::endl;
}

/*
This function calculates and updates the force vector for all the particle
*/
void ProblemFormulate::updateForce()
{
	const size_t N = molecules_.size();
	real r_ij, powVal, del_x_1, del_x_2, del_x_3;
	
	for(size_t i = 0; i<N; ++i)
	{
//computing the total force vector  acting on particle i 	
		for(size_t j=0; j<N; ++j)
		{

		if(i!=j)
		{		
			r_ij = calcDistance(i,j);

			if(r_ij>0)
			{	
				del_x_1 = (molecules_[j].x_[0] - molecules_[i].x_[0]);
				del_x_2 = (molecules_[j].x_[1] - molecules_[i].x_[1]);
				del_x_3 = (molecules_[j].x_[2] - molecules_[i].x_[2]);	
				powVal = pow((s/r_ij),6);
	
				molecules_[i].f_[0] += 24*e*(1/(r_ij*r_ij)) * powVal * (1- (2*powVal))*del_x_1;
				molecules_[i].f_[1] += 24*e*(1/(r_ij*r_ij)) * powVal * (1- (2*powVal))*del_x_2;
				molecules_[i].f_[2] += 24*e*(1/(r_ij*r_ij)) * powVal * (1- (2*powVal))*del_x_3;
	
			}
		}
		}
	}
	
}

/*
This function calculates and updates the position vector for all the particle
*/
void ProblemFormulate::updatePosition()
{
	const size_t N = molecules_.size();

	for(size_t i =0; i<N; ++i)
	{
//updating the position of particle i
		for(size_t k =0; k<3; ++k)
		{
			molecules_[i].x_[k] += (del_t*molecules_[i].v_[k]) + (0.5*(molecules_[i].f_[k])*del_t*del_t/(molecules_[i].m_));
		}
//Position update for particle i completed, time to apply periodic BC
		applyPeriodicBC(i);	
	}
}

/*
This function is for copying the forces f_ into the old forces vector fOld_
*/
void ProblemFormulate::saveOldForce()
{
	const size_t N = molecules_.size();

	for(size_t i =0; i<N; ++i)
	{
		for(size_t k =0; k<3; ++k)
			molecules_[i].fOld_[k] = molecules_[i].f_[k];	
	}
}

/*
This function is used to reset the force vector f_ to 0 after each iteration
*/
void ProblemFormulate::resetForces()
{
	const size_t N = molecules_.size();

	for(size_t i =0; i<N; ++i)
	{
		for(size_t k =0; k<3; ++k)
			molecules_[i].f_[k] = 0.0;	
	}
	

}


/*
This function calculates and updates the velocity vector for all the particle
*/
void ProblemFormulate::updateVelocity()
{
	const size_t N = molecules_.size();

	for(size_t i =0; i<N; ++i)
	{
		for(size_t k =0; k<3; ++k)
		{
			molecules_[i].v_[k] += 0.5*del_t*(molecules_[i].f_[k] + molecules_[i].fOld_[k])/molecules_[i].m_;
		}
		
	}
}

/*
This function returns the distance between particles i and j
*/
real ProblemFormulate::calcDistance(const size_t& i, const size_t& j)
{
	real diff = 0.0;	

	for(size_t k=0;k<3; ++k)
	{
		diff += pow((molecules_[i].x_[k] - molecules_[j].x_[k]),2);
	}
 
	return (std::sqrt(diff));
}

/*
This function applies periodic Boundary Condition to particle i, once its position has been updated
*/
void ProblemFormulate::applyPeriodicBC(const size_t& i)
{
	real x = molecules_[i].x_[0];	real y = molecules_[i].x_[1];	real z = molecules_[i].x_[2];

	if(x<x_min)
		x = x_max;

	if(x>x_max)
		x = x_min;
	
	if(y<y_min)
		y = y_max;

	if(y>y_max)
		y = y_min;
		
	if(z<z_min)
		z = z_max;

	if(z>z_max)
		z = z_min;
}




/*
This function writes the output for the n-th timestep to a .vtk file 
*/
void ProblemFormulate::writeOutput(const size_t n)
{
	const size_t N = molecules_.size();
	std::ofstream file_output;
	std::ostringstream nameVar;
	std::string filename;
	
//generating the outputfile name
		nameVar << title << n << ".vtk";
		filename = nameVar.str();
		
//once the filename is generated, we need to reset the ostringstream		
		nameVar.str("");
		//nameVar.clear();

		std::cout<<"Now writing the output to the file "<<filename<<std::endl;
//now opening the file for writing the output
		file_output.open(filename);
		file_output<<"# vtk DataFile Version 3.0"<<std::endl;
		file_output<<"SiWiR-2 Molecular Dynamics Visualization"<<std::endl;
		file_output<<"ASCII"<<std::endl;
		file_output<<"DATASET UNSTRUCTURED_GRID"<<std::endl;

// wrting the positions data		
		file_output<<"POINTS "<<N<<" DOUBLE"<<std::endl;
		for(size_t i=0;i<N;++i)
		{
			file_output<<molecules_[i].x_[0]<<" "<<molecules_[i].x_[1]<<" "<<molecules_[i].x_[2]<<std::endl;
		}

// writing the mass data
		file_output<<"POINT_DATA "<<N<<std::endl;
		file_output<<"SCALARS mass double"<<std::endl;
		file_output<<"LOOKUP_TABLE default"<<std::endl;
		for(size_t i=0;i<N;++i)
		{
			file_output<<molecules_[i].m_<<std::endl;
		}

// writing the force data
		file_output<<"VECTORS force double"<<std::endl;
		for(size_t i=0;i<N;++i)
		{
			file_output<<molecules_[i].f_[0]<<" "<<molecules_[i].f_[1]<<" "<<molecules_[i].f_[2]<<std::endl;
		}

// writing the velocity data
		file_output<<"VECTORS velocity double"<<std::endl;
		for(size_t i=0;i<N;++i)
		{
			file_output<<molecules_[i].v_[0]<<" "<<molecules_[i].v_[1]<<" "<<molecules_[i].v_[2]<<std::endl;
		}

		file_output.close();
	
	std::cout<<"Completed Writing the Output"<<std::endl;
}




