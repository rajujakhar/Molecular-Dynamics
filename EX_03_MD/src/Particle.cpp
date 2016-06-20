#include"Particle.hpp"

#include<iostream>



/*
This is the constructor of the Prticle class to intialize the particle properties
*/
Particle::Particle(const real& m, const real& x, const real& y, const real& z, const real& u, const real& v, const real&w )
{
	m_ = m;
	x_[0] = x;	x_[1] = y;	x_[2] = z;
	v_[0] = u;	v_[1] = v;	v_[2] = w;
	f_[0] = 0;	f_[1] = 0;	f_[2] = 0;
	fOld_[0] = 0;	fOld_[1] = 0;	fOld_[2] = 0;

}

/*
This id the destructor of the Particle class
*/
Particle::~Particle()
{

}

/*
This function is used to display the details of a Particle
*/
void Particle::displayParticle(const size_t& i)
{
	std::cout<<i<<"\t mass: "<<m_<<" x1: "<<x_[0]<<" x2: "<<x_[1]<<" x3: "<<x_[2]<<" v1: "<<v_[0]<<" v2: "<<v_[1]<<" v3: "<<v_[2]<<std::endl;
}
