#ifndef PARAM_READER_HPP
#define PARAM_READER_HPP

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
private:
    std::vector<std::string> keys_;
    std::map<std::string, double> globalMap_;

public:
    ParameterReader();
    void readParameters(const std::string& filename);
    inline bool isDefined(const std::string& key);
    void displayMap();
    template<typename Type>
    inline void getParameter(const std::string& key, const Type& val);

};
#endif
