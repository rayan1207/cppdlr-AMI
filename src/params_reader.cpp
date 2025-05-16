#include "dlr_ami.hpp"
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    if (first == std::string::npos || last == std::string::npos)
        return "";
    return str.substr(first, (last - first + 1));
}

// Function to parse the parameter name and value from a line
void parseLine(const std::string& line, std::string& paramName, std::string& paramValue) {
    size_t equalPos = line.find('=');
    if (equalPos != std::string::npos) {
        paramName = trim(line.substr(0, equalPos));
        paramValue = trim(line.substr(equalPos + 1));
    }
}

// Function to fill the struct from the provided text file
void params_loader(const std::string& filename, params_param& params) {


    // Open the file and read its contents
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::string paramName, paramValue;
        parseLine(line, paramName, paramValue);
        if (paramName == "Uval")
            params.Uval = std::stod(paramValue);
        else if (paramName == "Emax")
            params.Emax = std::stod(paramValue);
        else if (paramName == "eps")
            params.eps = std::stod(paramValue);
        else if (paramName == "beta")
            params.beta = std::stod(paramValue);
        else if (paramName == "L")
            params.L = std::stoi(paramValue);
		else if (paramName == "iter")
            params.iter = std::stoi(paramValue);
    
    }

    inputFile.close();
	std::cout << "Loaded params file with values: " << " Beta="  << params.beta << ", Emax="<<  params.Emax 
	<< ", eps = " << params.eps << ", L= " << params.L << " and SCS iteration = " << params.iter;
}