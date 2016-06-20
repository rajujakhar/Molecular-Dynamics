#include "ParameterReader.hpp"

int main()
{
    ParameterReader param;
    param.readParameters(std::string("blocks.par"));
    param.displayMap();
}
