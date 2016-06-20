#include "ParameterReader.hpp"

//fill the keys_ of the Parameter
ParameterReader::ParameterReader()
{
	keys_.push_back("name");
	keys_.push_back("vis_space");
	keys_.push_back("t_start");
	keys_.push_back("t_end");
	keys_.push_back("delta_t");
	keys_.push_back("x_min");
	keys_.push_back("y_min");
	keys_.push_back("z_min");
	keys_.push_back("x_max");
	keys_.push_back("y_max");
	keys_.push_back("z_max");
	keys_.push_back("r_cut");
	keys_.push_back("epsilon");
	keys_.push_back("sigma");
}

//Check is the given key is defind in the existing database
bool ParameterReader::isDefined(const std::string &key)
{
    auto iter = std::find(keys_.begin(), keys_.end(), key);

    if(iter== keys_.end())
    {	
        std::cout << "Key " << key << "is not valid\n";
        return 0;    //key not found
    }
    else
        return 1;    // key found
}

// Read input parameters from file
void ParameterReader::readParameters(const std::string &file)
{
    std::string line, key, val;
    char delim = ' ';
    std::ifstream input_file(file);

    assert(input_file.fail() == 0);

    getline(input_file, line);
    //std::istringstream ss(line);
    //getline(ss, key, delim);
    //getline(ss, val);

    // Insert the first line in the map
    //getParameter(key, val);

    while(getline(input_file, line))
    {

        std::istringstream ss(line);              //create a input string stream

        getline(ss, key, delim);       // find key
        getline(ss, val);             // find val

        //std::cout << "key is: " << key << std::endl;
        //std::cout << "Value is: " << val << std::endl;
        getParameter(key, std::stod(val));
        //getParameter(key, val);
    }
     input_file.close();
}

template<typename Type>
void ParameterReader::getParameter(const std::string& key, const Type& value)
{
    if(isDefined(key))
    {
        globalMap_.insert(std::pair<std::string, Type> (key,value));
    }
}

void ParameterReader::displayMap()
{
    std::cout << "Displaying the contents of map\n";
    std::cout << "map size: " << globalMap_.size() <<  std::endl;

    for(auto iter = globalMap_.begin(); iter!=globalMap_.end(); ++iter)
    {
        std::cout << iter->first << "\t" << iter->second << std::endl;
    }

}





