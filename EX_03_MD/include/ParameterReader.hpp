#ifndef PARAM_READER_HPP
#define PARAM_READER_HPP

#include"Type.hpp"

#include <cstddef> // for size_t
#include <iostream>
#include <string>
#include <vector>
#include <fstream>  //std::ifstream
#include <sstream>
#include <map>
#include <algorithm> //find
#include <assert.h>
//#include <utility>

class ParameterReader
{
public :
    

	ParameterReader();
	void readParameters(const std::string& filename);
	inline bool isDefined(const std::string& key);
	void displayMap();
	//template<typename Type>
	//inline void getParameter(const std::string& key, double &value) const;

	void setParams();

	//private:
		std::vector<std::string> keys_;
		std::map<std::string, double> globalMap_;
		
		/*
		const static std::string name = "blocks";
		const static size_t vis_space = 1600;
		const static real t_start = 0.0;
		const static real t_end = 8.0;
		const static real delta_t = 0.00005;
		const static real x_min = 0.0;
		const static real y_min = 0.0;
		const static real z_min = 0.0;
		const static real x_max = 224.4924;
		const static real y_max = 224.4924;
		const static real z_max = 224.4924;
		const static real r_cut = 2.5;
		const static real epsilon = 5.0;
		const static real sigma = 1.0;
		*/

};
#endif
