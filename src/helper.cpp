#include "dlr_ami.hpp"
using namespace cppdlr;


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




void mDLR::write_data_momenta(const std::string& filename,
                            Bz_container& data,
                            nda::array<dcomplex,1>& mfreq  )
{
    // Extract the NDA 1D array for the (i,j) momenta
 
    int size = data[0][0].size();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Optional: write header
    //file << "# wn    qx     qy     Re       Im\n";
    for (int i =0; i<kl;++i){
		for (int j =0; j<kl;++j){		
			for (int f = 0; f < size; ++f) {
				double wn_imag = mfreq[f].imag();
				double qx = kvals[i];
				double qy = kvals[j];
				double re = data[i][j](f).real();
				double im = data[i][j](f).imag();

			
				// Write to file
				file << wn_imag << "  "
					 << qx << "  "
					 << qy << "  "
					 << re << "  "
					 << im << "\n";
					}
		}
	}

    file.close();
    std::cout << "Data written to " << filename << std::endl;
}
	
	
	
void mDLR::write_data_ij_momenta(const std::string& filename,
                            Bz_container& data,
                            nda::array<dcomplex,1>& mfreq,
                            std::pair<int, int> ij)
{
    // Extract the NDA 1D array for the (i,j) momenta
    auto array1d = data[ij.first][ij.second];
    int size = array1d.size();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Optional: write header
    file << "# wn_imag    qx     qy     Re[G]       Im[G]\n";

    for (int i = 0; i < size; ++i) {
        double wn_imag = mfreq[i].imag();
        double qx = kvals[ij.first];
        double qy = kvals[ij.second];
        double re = array1d[i].real();
        double im = array1d[i].imag();

        // Print to console
        std::cout << wn_imag << "  "
                  << qx << "  "
                  << qy << "  "
                  << re << "  "
                  << im << std::endl;

        // Write to file
        file << wn_imag << "  "
             << qx << "  "
             << qy << "  "
             << re << "  "
             << im << "\n";
    }

    file.close();
    std::cout << "Data written to " << filename << std::endl;
}


AmiBase::g_prod_t construct_example2(){

	AmiBase::g_prod_t g;

	// Setting up G array
	// defining alpha's /// std::vector<int>


	AmiBase::alpha_t alpha_1={1,0,0};
	AmiBase::alpha_t alpha_2={0,1,0};
	AmiBase::alpha_t alpha_3={-1,1,1};

	//defining epsilon's
	AmiBase::epsilon_t epsilon_1={1,0,0};
	AmiBase::epsilon_t epsilon_2={0,1,0};
	AmiBase::epsilon_t epsilon_3={0,0,1};

	AmiBase::g_struct g1(epsilon_1,alpha_1);
	AmiBase::g_struct g2(epsilon_2,alpha_2);
	AmiBase::g_struct g3(epsilon_3,alpha_3);

	// AmiBase::g_prod_t R0={g1,g2,g3};

	// OR
	AmiBase::g_prod_t R0;
	R0.push_back(g1);
	R0.push_back(g2);
	R0.push_back(g3);
	return R0;

}

std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec) {
    std::vector<std::complex<double>> cplx_vec;
    cplx_vec.reserve(vec.size()); 

    for (auto elem : vec) {
        cplx_vec.emplace_back(elem, 0.0);
    }

    return cplx_vec;
}


nda::array<dcomplex,1> mDLR::evaluate_auxillary_weights(nda::array<double,1> &energy) {
	nda::array<dcomplex ,1> weights;
	std::vector<nda::array<dcomplex,1>> G_dlr_list;
	std::vector<nda::array<dcomplex,1>> G_dlr_w_list;
	for (int i =0; i< N;i++){
		auto dlr_R0 =  multiple_dlr_structs[i];
		auto gdlr_R0 = generate_nda_Gdlr_from_energy(dlr_R0.if_ops, energy[i]);
		auto weights = dlr_R0.if_ops.vals2coefs(beta, gdlr_R0);
		G_dlr_list.push_back(gdlr_R0);
		G_dlr_w_list.push_back(weights);
	}
	
	
	nda::array<dcomplex,1> full_weights (CN);
	for (int i=0; i < CN;i++){
		auto const& combo  = cartesian_combo_list[i];
	    auto result = dcomplex(1,0);
		for (int j =0; j< N ; j++){
			result = -1*result * G_dlr_w_list[j](combo[j]);
		}
		full_weights(i) = result;	
	}
	
	return full_weights;
}


double hubbard_dispersion(double kx, double ky){
  return -2.0*(std::cos(kx) + std::cos(ky));
}


