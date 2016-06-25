#ifndef PARAM_READER_HPP
#define PARAM_READER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>   //std::ifstream
#include <sstream>
#include <map>
#include <algorithm>    //find
#include <assert.h>
//#include <type_traits>   // std::is_same

class ParameterReader
{
private:
    std::vector<std::string> keys_;
    std::map<std::string, std::string> globalMap_;

public:
    ParameterReader();
    void readParameters(const std::string&);
    inline bool isDefined(const std::string&);
    void displayMap();
    std::string trim(std::string&);

    void getParameter(const std::string&, std::string& );
    void getParameter(const std::string&, double& );
    void getParameter(const std::string&, size_t& );
    //template<typename Type>
    //void getParameter(const std::string& key, Type& value);
    //void getParameter(const std::string&, std::string& );

   /*template<typename Type>
    void getParameter(const std::string& key, Type& value);
    {
        if (std::is_same<Type, double>::value)
        {
            std::cout << "Type is double\n";
            value = std::stod(globalMap_[key]);
        }
        else if(std::is_same<Type, std::string>::value)
        {
            std::cout << "Type is string\n";
           value = globalMap_[key];
            //value = "blocks";
        }

    }*/

};


/*template<typename Type>
void ParameterReader::getParameter(const std::string& key, Type& value)
{
    value = globalMap_[key];
}*/

#endif