#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "Type.hpp"

#include <cstddef> // for size_t

class Particle
{
	public:
		
	Particle(const real&, const real&, const real&, const real&, const real&, const real&, const real&); // constrcutor to intialize the particle parameters		
		~Particle(); //destructor
		void displayParticle(const size_t&);//to display the particle details
	
	//private:
		real m_;
		real x_[3];
		real v_[3];
		real f_[3];
		real fOld_[3];
		
};

#endif
