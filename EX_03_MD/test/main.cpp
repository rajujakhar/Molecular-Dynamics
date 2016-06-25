#include<iostream>
#include<string>

#include "ProblemFormulate.hpp"

int main(int argc, char*argv[])
{

	if(argc<3)
		std::cout<<"Unsufficient number of input parameters"<<std::endl;
	
	std::string paramsFileName = argv[1];
	std::string dataFileName = argv[2];

	std::cout<<"The params file is "<<paramsFileName<<std::endl;
	std::cout<<"The data file is "<<dataFileName<<std::endl;

	ProblemFormulate MDsim(paramsFileName);
	MDsim.readDataInputFile(dataFileName);
	MDsim.simulateMD();
	
	return 0;	


}
