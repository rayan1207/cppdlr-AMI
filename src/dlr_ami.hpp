#ifndef DLR_AMI_HPP
#define DLR_AMI_HPP
#include "ami_calc.hpp"
#include "amigraph.hpp"
#include <cassert>
#include <experimental/filesystem>
#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <unistd.h>
#include <mpi.h>
#include <chrono>
#include <thread>
#include <sstream>
#include <boost/random/sobol.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <typeinfo> 
#include <string>
#include <utility> 
#include <locale>
#include <stdexcept>
#include <tuple>
#include <Eigen/Dense> 
#include <complex>
#include <iostream>
#include <iomanip>
#include <cppdlr/cppdlr.hpp>




AmiBase::g_prod_t construct_example2();
AmiBase::ami_vars construct_ext_example2();


struct dlr_obj{
	cppdlr::imfreq_ops if_ops;
	//Eigen::VectorXcd im_freqs;
    std::vector<std::complex<double>> im_freqs;
	std::vector<double> pole_locs;
	AmiBase::g_struct ginfo;
	int pole_num;
	std::vector<std::vector<double>> evec;
	std::vector<std::vector<nda::array<dcomplex,1>>> dlrW_in_square;
	
};
class mDLR{
	public:
	std::vector<dlr_obj> multiple_dlr_structs;
	double beta; double eps; double Emax; AmiBase::g_prod_t R0;
	size_t N; ///num of greens function
	size_t CN; ///total number of cartesian, pole_num1* pole_num2* ...pole numN 
	size_t kl;///total number of momentum k grid;
	size_t kN;///total number of cartesian momenta, kl_1^2* kl_2^2.....
	std::vector<double> kvals;
	size_t master_pole_num;
	cppdlr::imfreq_ops master_if_ops;
	std::vector<std::vector<nda::array<dcomplex,1>>> master_dlrW_in_square;
	std::vector<std::vector<int>> cartesian_combo_list;
	std::vector<std::vector<int>> cartesian_k_combo_list;
	
	std::vector<std::vector<double>> auxillary_energy_list;
	mDLR(double _beta, double _eps, double _Emax,size_t kl,AmiBase::g_prod_t _R0);
	//////// methods ////////////
	void create_DLR_master_if_ops();
	
	nda::array<dcomplex,1> generate_nda_Gdlr_from_energy( cppdlr::imfreq_ops &ops, double &energy);
	void create_multiple_gstruct();
	void generate_cartesian_list();
	void generate_auxillary_energy_list();
	nda::array<dcomplex,1> evaluate_auxillary_energies(std::complex<double> &imfreq);
	nda::array<dcomplex,1> evaluate_auxillary_weights( nda::array<double,1> &energy);
	void populate_master_dlrW_from_G0();
	void reshape_dlrW_square_per_kgrid();
	nda::array<dcomplex,1> recover_dlro_G_from_master_weights(nda::array<dcomplex,1> &master_weights, std::vector<std::complex<double>> &dlro_if);
	void transfer_master_DLR_weights_to_dlrR0_elements();
	void generate_momenta_cartesian_combo();
};

dlr_obj create_dlr_obj(double beta, double eps, double Emax,AmiBase::g_struct R0_element);

std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec);




nda::array<dcomplex,1> recover_G_from_poles_n_weights(dlr_obj& dlr_R0, nda::array<dcomplex,1> weights, std::vector<std::complex<double>> imfreqs);

double hubbard_dispersion(double kx, double ky);


template<typename T>
std::vector<T> sumVectors(const std::vector<std::vector<T>>& vecs);


////////////////////////////////////// Print stuff ////////////////////////////


template<typename T>
inline void print2d( std::vector< std::vector<T>> vec)
{
    for ( auto row : vec) {
        for ( auto elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
	std::cout << std::endl;
	
}


template<typename T> 
inline void print1d( std::vector<T>& vec) {
  std::cout << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i != vec.size() - 1) {
      std::cout << ", ";
    }
  }
  std::cout << "]\n";
}

template<typename T>
inline void cprint1d(const std::vector<std::complex<T>>& vec)
{
    for (const auto& elem : vec) {
        std::cout << "(" << elem.real() << "," << elem.imag() << ") ";
    }
    std::cout << std::endl;
}
template<typename T>
inline void cprint2d(const std::vector<std::vector<std::complex<T>>>& vec)
{
    for (const auto& row : vec) {
        cprint1d(row); // Use cprint1d to print each row
    }
    std::cout << std::endl;
}




#endif // DLR_AMI_HPP