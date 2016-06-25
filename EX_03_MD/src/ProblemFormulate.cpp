
#include "ProblemFormulate.hpp"
#include "ParameterReader.hpp"

#include<iostream>
#include<fstream> // ifstream, ofstream
#include<sstream> //ostringstream 
#include<cmath>	 //pow(x,n)
#include<assert.h>

/*
 The construcor which sets the parameter variables
 */
ProblemFormulate::ProblemFormulate(const std::string& filename)
{
    ParameterReader param;
    param.readParameters(filename);
    param.getParameter(std::string("name"), name);
    param.getParameter(std::string("vis_space"), vis_space);
    param.getParameter(std::string("t_start"), t_start);
    param.getParameter(std::string("t_end"), t_end);
    param.getParameter(std::string("delta_t"), delta_t);
    param.getParameter(std::string("x_min"),x_min );
    param.getParameter(std::string("y_min"),y_min );
    param.getParameter(std::string("z_min"),z_min );
    param.getParameter(std::string("x_max"),x_max );
    param.getParameter(std::string("y_max"),y_max );
    param.getParameter(std::string("z_max"),z_max );
	param.getParameter(std::string("r_cut"),r_cut );
    param.getParameter(std::string("epsilon"),e );
    param.getParameter(std::string("sigma"), s);

    std::cout << "name is: " << name << std::endl;
  

	
    std::cout << "sigma: " << s << std::endl;
}



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
	const size_t n_steps = (t_end-t_start)/delta_t;
	std::cout<<"This simulation has "<<n_steps<<" number of timesteps"<<std::endl;	
	std:: cout<<"Starting the simulation"<<std::endl;	

//writing the intial data to the output file
	writeOutput(0);

//setting up the linked-cell framework
	setUpLinkedCells();


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
		
		if(i%vis_space ==0)
			writeOutput(i);
		
	}
	std::cout<<"Simulation Complete !!!!!"<<std::endl;


}

/*
This function calculates and updates the force vector for all the particle
*/
void ProblemFormulate::updateForce()
{
	const size_t n_c = cells_.size();
	std::vector<size_t> neighbourSet;
	size_t limitL, currCellNum,i,j;
	std::list<size_t>::iterator it_1, it_1_beg, it_1_end, it_2, it_2_beg, it_2_end;
	real r_ij, powVal, del_x_1, del_x_2, del_x_3;	

//looping over all the cells
	for(size_t k=0; k<n_c; ++k)
	{

//first check if this cell is not empty and then only proceed further
		if(!(cells_[k].empty()))
		{
//Now determing the neighbours of this cell and storing them locally
			neighbourSet.clear();
			determineNeighbourSet(k, neighbourSet); 
			limitL = neighbourSet.size();
			
//now looping over all perticles in cell-k
			it_1_beg = cells_[k].begin();
			it_1_end = cells_[k].end();
			for(it_1=it_1_beg; it_1!=it_1_end; ++it_1)
			{
				i = (*it_1);
//now looping over all cells in neighbourset of cell-k 				
				for(size_t l=0; l<limitL; ++l)
				{
					currCellNum = neighbourSet[l];
					it_2_beg = cells_[currCellNum].begin();
					it_2_end = cells_[currCellNum].end();
					for(it_2=it_2_beg; it_2!=it_2_end; ++it_2)
					{
					if(it_1 !=it_2 )
					{
						
						j = (*it_2);
						r_ij = calcDistance(i,j);

						if(r_ij<=r_cut)
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
		}
	}
}



/*
This function determines the r_i_c, n_i_c and the particles lying in each cell
*/
void ProblemFormulate::setUpLinkedCells()
{
	std::cout<<"Now Setting up the Linked Cells"<<std::endl;

// Now determining the values of n_i_c and r_i_c 
	const real l1 = fabs(x_max-x_min);
	const real l2 = fabs(y_max-y_min);
	const real l3 = fabs(z_max-z_min);
	
	std::cout << "x_min: " << x_min << std::endl;
	std::cout << "x_max: " << x_max << std::endl;
	
    std::cout << "y_min: " << y_min << std::endl;
	std::cout << "y_max: " << y_max << std::endl;
	
    std::cout << "z_min: " << z_min << std::endl;
	std::cout << "z_max: " << z_max << std::endl;

	n_1_c = l1/r_cut;	n_2_c = l2/r_cut;	n_3_c = l3/r_cut;		
	r_1_c = l1/n_1_c;	r_2_c = l2/n_2_c;	r_3_c = l3/n_3_c;



	std::cout<<"The cell parameters are: "<<std::endl;
	std::cout<<"l1 : "<< l1<<" l2 : "<<l2<<" l3 is : "<<l3<<std::endl;
	std::cout<<"n_1_c : "<< n_1_c<<" n_2_c : "<<n_2_c<<" n_3_c is : "<<n_3_c<<std::endl;
	std::cout<<"r_1_c : "<< r_1_c<<" r_2_c : "<<r_2_c<<" r_3_c is : "<<r_3_c<<std::endl;

		assert(n_1_c>0);	assert(n_2_c>0);	assert(n_3_c>0);

//Now resizing the cells_ data structure
	n_tot_c  = n_1_c * n_2_c * n_3_c;
	cells_.resize(n_tot_c);
	
	std::cout<<"size of cells is "<<cells_.size()<<std::endl;


// Now determing the particles in each cell and updating the cells_ data structure
	const size_t N = molecules_.size();
	size_t n;	
	for(size_t i = 0; i<N; ++i)
	{
		n = determineParticleLoc(i);
		cells_[n].push_back(i);
	//std::cout<<"x:\t"<<x<<"\ty:\t"<<y<<"\tz:\t"<<z<<"\tn1:\t"<<n1<<"\tn2:\t"<<n2<<"\tn3:\t"<<n3<<"\tn:\t"<<n<<std::endl;
	}

}



/*
This fuctions calculates the position of a cell in the cells_ vector from its coordinates
*/
size_t ProblemFormulate::calcCellPos(const size_t& n1, const size_t& n2, const size_t& n3)
{
	return ((n_1_c * n_2_c * n3) + (n_1_c*n2) + n1);

}


/*
This fuctions calculates the coordinates of a cell from its position in the cells_ vector 
*/
void ProblemFormulate::calcCellCoords(const size_t& n, std::array<size_t,3>& tmpArr)
{	
		
	tmpArr[2] = n/(n_1_c*n_2_c);
	tmpArr[1] = (n%(n_1_c*n_2_c))/n_1_c; 
	tmpArr[0] = (n%(n_1_c*n_2_c))%n_1_c;
	

	if(tmpArr[0]>=n_1_c)
	{
	std::cout<<"Before Update, n1 exceeds limit. The value of n is "<<n<<" The value of n1 is "<<tmpArr[0]<<std::endl;
	exit(EXIT_FAILURE);
	}	
	if(tmpArr[1]>=n_2_c)
	{
	std::cout<<"Before Update, n2 exceeds limit. The value of n is "<<n<<" The value of n2 is "<<tmpArr[1]<<std::endl;
	exit(EXIT_FAILURE);
	}
	if(tmpArr[2]>=n_3_c)
	{
	std::cout<<"Before Update, n3 exceeds limit. The value of n is "<<n<<" The value of n3 is "<<tmpArr[2]<<std::endl;
	exit(EXIT_FAILURE);
	}


	/*	
	assert(tmpArr[0]<n_1_c && tmpArr[0]>=0);	
	assert(tmpArr[1]<n_2_c && tmpArr[1]>=0);
	assert(tmpArr[2]<n_3_c && tmpArr[2]>=0);
	*/
}

/*
This function determines the position of the cell in which particle k lies 
*/
size_t ProblemFormulate::determineParticleLoc(const size_t& k)
{
	real x,y,z;
	size_t n1,n2,n3;	
	x = molecules_[k].x_[0];	y = molecules_[k].x_[1];	z = molecules_[k].x_[2];

//determining the coordinates of the cell in which the particle lies
/*		
		assert(x>=x_min);	assert(y>=y_min);	assert(z>=z_min); 
		assert(x<=x_max);	assert(y<=y_max);	assert(z<=z_max);
*/
		n1 = (fabs(x-x_min))/r_1_c;		n2 = fabs((y-y_min))/r_2_c;		n3 = fabs((z-z_min))/r_3_c;

		if(x==x_max)
			n1 = n_1_c-1; //correction required if particle is on the end boundary
		if(y==y_max)
			n2 = n_2_c-1;
		if(z==z_max)
			n3 = n_3_c-1;

/*
		assert(n1<n_1_c);
		assert(n2<n_2_c);
		assert(n3<n_3_c);
				
*/


/*
		assert(n1<n_1_c && n1>=0);
		assert(n2<n_2_c && n2>=0);
		assert(n3<n_3_c && n3>=0);
*/		


//calculating and returning the position of the cell in the cells_ vector
		return (calcCellPos(n1,n2,n3));
}





/*
This function calculates and updates the position vector for all the particle
*/
void ProblemFormulate::updatePosition()
{
	const size_t N = molecules_.size();
	size_t oldLoc, newLoc;
	real x;		
	real y;
	real z;
	std::array<size_t,3> oldCellCords, newCellCords;
	for(size_t i =0; i<N; ++i)
	{
// determing the old location of this particle	
		
		x = molecules_[i].x_[0];	
		y = molecules_[i].x_[1];	
		z = molecules_[i].x_[2];
		
		if(x<x_min)
		{
			std::cout<<"Before position update, x is less than x_min for particle "<<i<<" The value of x is "<<x<<std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(y<y_min)
		{
			std::cout<<"Before position update, y is less than y_min for particle "<<i<<" The value of y is "<<y<<std::endl;
			exit(EXIT_FAILURE);
		}
		if(z<z_min)
		{
			std::cout<<"Before position update, z is less than z_min for particle "<<i<<" The value of z is "<<z<<std::endl;
			exit(EXIT_FAILURE);
		}

						
		oldLoc = determineParticleLoc(i);

//determining the coordinates of this cell and asserting
		calcCellCoords(oldLoc, oldCellCords);
		if (oldCellCords[0]>=n_1_c)
		{		
			std::cout<<"Before Update, n1 exceeds limit for particle "<<i<<" The value of n is "<<oldLoc<<" The value of n1 is "<<oldCellCords[0]<<std::endl;
			exit(EXIT_FAILURE);		
		}		
		if(oldCellCords[1]>=n_2_c)
		{		
			std::cout<<"Before Update, n2 exceeds limit for particle "<<i<<" The value of n is "<<oldLoc<<" The value of n2 is "<<oldCellCords[1]<<std::endl;
			exit(EXIT_FAILURE);		
		}		
		if(oldCellCords[2]>=n_3_c)
		{		
			std::cout<<"Before Update, n3 exceeds limit for particle "<<i<<" The value of n is "<<oldLoc<<" The value of n3 is "<<oldCellCords[2]<<std::endl;
			exit(EXIT_FAILURE);
		}
//updating the position of particle i
		for(size_t k =0; k<3; ++k)
		{
			molecules_[i].x_[k] += (delta_t*molecules_[i].v_[k]) + (0.5*(molecules_[i].f_[k])*delta_t*delta_t/(molecules_[i].m_));
		}
		
		x = molecules_[i].x_[0];	
		y = molecules_[i].x_[1];	
		z = molecules_[i].x_[2];
		
		if(x<x_min)
		{
			std::cout<<"After position update, x is less than x_min for particle "<<i<<" The value of x is "<<x<<std::endl;
		}
	


//Position update for particle i completed, time to apply periodic BC
		applyPeriodicBC(i);

// determing the new location of this particle	
		newLoc = determineParticleLoc(i);

		x = molecules_[i].x_[0];	
		y = molecules_[i].x_[1];	
		z = molecules_[i].x_[2];
		
		if(x<x_min)
		{
			std::cout<<"After applying Periodic BC, x is less than x_min for particle "<<i<<" The value of x is "<<x<<std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(x>x_max)
		{
			std::cout<<"After applying Periodic BC, x is greater than x_max for particle "<<i<<" The value of x is "<<x<<std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(y<y_min)
		{
			std::cout<<"After applying Periodic BC, y is less than y_min for particle "<<i<<" The value of y is "<<y<<std::endl;
			exit(EXIT_FAILURE);
		}
		

		if(y>y_max)
		{
			std::cout<<"After applying Periodic BC, y is greater than y_max for particle "<<i<<" The value of y is "<<y<<std::endl;
			exit(EXIT_FAILURE);
		}
		
		if(z<z_min)
		{
			std::cout<<"After applying Periodic BC, z is less than z_min for particle "<<i<<" The value of z is "<<z<<std::endl;
			exit(EXIT_FAILURE);
		}

		if(z>z_max)
		{
			std::cout<<"After applying Periodic BC, z is greater than z_max for particle "<<i<<" The value of z is "<<z<<std::endl;
			exit(EXIT_FAILURE);
		}
/*
		assert(x>=x_min);	
		assert(y>=y_min);	
		assert(z>=z_min); 
		assert(x<=x_max);	
		assert(y<=y_max);	
		assert(z<=z_max);
*/


//determining the coordinates of this cell and asserting
		calcCellCoords(newLoc, newCellCords);
		if(newCellCords[0]>=n_1_c)
		{
			std::cout<<"After applying Periodic BC, n1 exceeded limit for particle "<<i<<" The value of n is "<<newLoc<<" The value of n1 is "<<newCellCords[0]<<std::endl;
			exit(EXIT_FAILURE);
		}		
		if(newCellCords[1]>=n_2_c)
		{		
			std::cout<<"After applying Periodic BC, n2 exceeded limit for particle "<<i<<" The value of n is "<<newLoc<<" The value of n2 is "<<newCellCords[1]<<std::endl;
			exit(EXIT_FAILURE);
		}
		if(newCellCords[2]>=n_3_c)
		{
			std::cout<<"After applying Periodic BC, n3 exceeded limit for particle "<<i<<" The value of n is "<<oldLoc<<" The value of n3 is "<<oldCellCords[2]<<std::endl;
			exit(EXIT_FAILURE);	
		}	

//checking if this particle has migrated and updating cells_ accordingly
		updateCells(i, oldLoc, newLoc);
	}
}

/*
This function takes the locations of particle-k before and after position update and updates the cells_ data strusture if the particle migrates from one cell to another
*/
void ProblemFormulate::updateCells(const size_t& k, const size_t& old_, const size_t& new_)
{
	if(new_ !=old_)
		{
			cells_[old_].remove(k);
			cells_[new_].push_back(k);

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
			molecules_[i].v_[k] += 0.5*delta_t*(molecules_[i].f_[k] + molecules_[i].fOld_[k])/molecules_[i].m_;
		}
		
	}
}

/*
This function returns the distance between particles i and j
*/
real ProblemFormulate::calcDistance(const size_t& i, const size_t& j)
{
	real diff = 0.0, r =0.0;	
	
	for(size_t k=0;k<3; ++k)
	{
		diff += pow((molecules_[i].x_[k] - molecules_[j].x_[k]),2);
	}
	assert(diff>0);
 	r = std::sqrt(diff);
	
	return r;
}

/*
This function applies periodic Boundary Condition to particle i, once its position has been updated
*/
void ProblemFormulate::applyPeriodicBC(const size_t& i)
{
	
	if(molecules_[i].x_[0]<x_min)
	{
		//std::cout<<"While applying Periodic BC, x is less than x_min for particle "<<i<<" The value of x is "<<molecules_[i].x_[0]<<std::endl;	
		molecules_[i].x_[0] = x_max;
		//std::cout<<"x has been updated to "<<molecules_[i].x_[0]<<std::endl;
	}

	if(molecules_[i].x_[0]>x_max)
		molecules_[i].x_[0] = x_min;
	

	if(molecules_[i].x_[1]<y_min)
		molecules_[i].x_[1] = y_max;

	if(molecules_[i].x_[1]>y_max)
		molecules_[i].x_[1] = y_min;
		
	if(molecules_[i].x_[2]<z_min)
	{
		//std::cout<<"While applying Periodic BC, z is less than z_min for particle "<<i<<" The value of z is "<<molecules_[i].x_[2]<<std::endl;	
		molecules_[i].x_[2] = z_max;
		//std::cout<<"z has been updated to "<<molecules_[i].x_[2]<<std::endl;
	}

	if(molecules_[i].x_[2]>z_max)
		molecules_[i].x_[2] = z_min;
}



/*
This function determines the 27 elements in neighbourSet of a cell 
*/
void ProblemFormulate::determineNeighbourSet(const size_t& i, std::vector<size_t>& NSet)
{
	std::array<size_t,3> cellCoords;
	size_t c_x, c_y, c_z;
	size_t tmp_n, tmp_c_x, tmp_c_y, tmp_c_z;

//first determining the coordinates of this cell
	calcCellCoords(i, cellCoords);

	/////////////////////////////////////////////////////////////////////////////
	/////Now determing the neighbours of this cell and storing them locally /////
	/////////////////////////////////////////////////////////////////////////////
	c_x = cellCoords[0];	c_y = cellCoords[1];	c_z = cellCoords[2];

// 1- inserting the MIDDLE centre cell to the neighbourSet 	
	NSet.push_back(i);
// 2 - determining and inserting the East Neighbour (here only x-coord changes)
	tmp_c_x = (c_x+1)%n_1_c;
	tmp_n = calcCellPos(tmp_c_x, c_y, c_z);
	NSet.push_back(tmp_n);
// 3 - determining and inserting the North-East Neighbour (here both x- and y-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 4 - determining and inserting the north Neighbour (here only y-coord changes)
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 5 - determining and inserting the North-West Neighbour (here both x- and y-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 6 - determining and inserting the West Neighbour (here only x-coord changes)
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;
	tmp_n = calcCellPos(tmp_c_x, c_y, c_z);
	NSet.push_back(tmp_n);
// 7 - determining and inserting the South-West Neighbour (here both x- and y-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 8 - determining and inserting the South Neighbour (here only y-coord changes)
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 9 - determining and inserting the South-East Neighbour (here both x- and y-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, c_z);
	NSet.push_back(tmp_n);
// 10 - determining and inserting the TOP centre Neighbour (here only z-coord changes)	
	tmp_c_z = (c_z+1)%n_3_c;	
	tmp_n = calcCellPos(c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 11 - determining and inserting the TOP East Neighbour (here both x- and z-coords change)
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_z = (c_z+1)%n_3_c;	
	tmp_n = calcCellPos(tmp_c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 12 - determining and inserting the TOP North-East Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 13 - determining and inserting the TOP north Neighbour (here both y- and z-coords change)
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 14 - determining and inserting the TOP North-West Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 15 - determining and inserting the TOP West Neighbour (here both x- and z-coords change)
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 16 - determining and inserting the TOP South-West Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 17 - determining and inserting the TOP South Neighbour (here both y- and z-coords change)
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 18 - determining and inserting the TOP South-East Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (c_z+1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 19 - determining and inserting the Bottom centre Neighbour (here only z-coord changes)	
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;	
	tmp_n = calcCellPos(c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 20 - determining and inserting the BOTTOM East Neighbour (here both x- and z-coords change)
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;	
	tmp_n = calcCellPos(tmp_c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 21 - determining and inserting the BOTTOM North-East Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 22 - determining and inserting the BOTTOM north Neighbour (here both y- and z-coords change)
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 23 - determining and inserting the BOTTOM North-West Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (c_y+1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 24 - determining and inserting the BOTTOM West Neighbour (here both x- and z-coords change)
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 25 - determining and inserting the BOTTOM South-West Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (n_1_c + c_x - 1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 26 - determining and inserting the BOTTOM South Neighbour (here both y- and z-coords change)
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
// 27 - determining and inserting the BOTTOM South-East Neighbour (here all x-, y- and z-coords change)	
	tmp_c_x = (c_x+1)%n_1_c;	
	tmp_c_y = (n_2_c + c_y - 1)%n_2_c;
	tmp_c_z = (n_3_c + c_z - 1)%n_3_c;
	tmp_n = calcCellPos(tmp_c_x, tmp_c_y, tmp_c_z);
	NSet.push_back(tmp_n);
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
		nameVar << name << n << ".vtk";
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




