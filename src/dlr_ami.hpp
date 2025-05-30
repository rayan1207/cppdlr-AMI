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
#include <omp.h>
#include <format>
#include <string>




AmiBase::g_prod_t construct_example2();
AmiBase::ami_vars construct_ext_example2();




//////params loader functions/////////////
using Bz_container =  std::vector<std::vector<nda::array<dcomplex,1>>>;
struct params_param {
	double Uval;
	double beta;
	double Emax;
	double eps;
	int iter;
	int L;
};

std::string trim(const std::string& str);
void parseLine(const std::string& line, std::string& paramName, std::string& paramValue);
void params_loader(const std::string& filename, params_param& params);





/////////// DLR stuff //////////////
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
	private:
	double two_pi     ;
    double inv_two_pi ;
    double inv_dk     ;
    double* kvals_ptr;
	public:
	std::vector<dlr_obj> multiple_dlr_structs;
	double beta; double eps; double Emax;double Uval; AmiBase::g_prod_t R0;
	size_t N; ///num of greens function
	size_t CN; ///total number of cartesian, pole_num1* pole_num2* ...pole numN 
	size_t kl;///total number of momentum k grid;
	size_t kN;///total number of cartesian momenta, kl_1^2* kl_2^2.....
	double dk;
	int ord;
	
	std::vector<double> kvals;/// kgrid vals
	size_t master_pole_num;/// number of poles in master DLR
	cppdlr::imfreq_ops master_if_ops; /// master dlr obj
	std::vector<std::vector<nda::array<dcomplex,1>>> master_dlrW_in_square; //master dlr weight
	std::vector<std::vector<int>> cartesian_combo_list; /// contains the cartesian product indices for poles
	std::vector<std::vector<int>> cartesian_k_combo_list; ///co
	
	std::vector<std::vector<double>> auxillary_energy_list;
	mDLR(double _beta,double _Uval, double _eps, double _Emax,size_t kl,AmiBase::g_prod_t _R0);
	//////// methods ////////////
	void create_DLR_master_if_ops();
	
	nda::array<dcomplex,1> generate_nda_Gdlr_from_energy( cppdlr::imfreq_ops &ops, double &energy);
	void create_multiple_gstruct();
	void generate_cartesian_list();
	void generate_auxillary_energy_list();
	nda::array<dcomplex,1> evaluate_auxillary_energies(nda::dcomplex &imfreq);
	nda::array<dcomplex,1> evaluate_auxillary_weights( nda::array<double,1> &energy);
	void populate_master_dlrW_from_G0();
	void reshape_dlrW_square_per_kgrid();
	nda::array<dcomplex,1> recover_dlro_G_from_master_weights(nda::array<dcomplex,1> &master_weights, std::vector<std::complex<double>> &dlro_if);
	void transfer_master_DLR_weights_to_dlrR0_elements();
	void generate_momenta_cartesian_combo();
	inline nda::dcomplex compute_momenta_one_kCN_kernel(double kx_ext,double ky_ext,const int* combo_ptr,const int* kcombo_ptr);
	nda::array<nda::dcomplex,1> compute_momenta_kernel_qext(double kx_ext,double ky_ext);
    Bz_container compute_momenta_kernel_bz();
	Bz_container vdot_freq_momenta_kernel_M(const std::vector<std::vector<nda::array<dcomplex,1>>> mk, const std::vector<nda::array<dcomplex,1>> fk);
	Bz_container G_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq);
	void repopulate_master_dlrW_from_G(Bz_container &G );
	void write_data_momenta(const std::string& filename,Bz_container& data, nda::array<dcomplex,1>& mfreq);
	void write_data_ij_momenta(const std::string& filename,
                            Bz_container& data,
                            nda::array<dcomplex,1>& mfreq,
                            std::pair<int, int> ij);
};

dlr_obj create_dlr_obj(double beta, double eps, double Emax,AmiBase::g_struct R0_element);

std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec);
template<typename T>
void triangle_to_square(std::vector<std::vector<T>>& M);
template<typename T>
std::vector<std::vector<T>> data_to_full_bz(const std::vector<std::vector<T>>& M);
template<typename T>
std::vector<T> sumVectors(const std::vector<std::vector<T>>& vecs);
nda::array<dcomplex,1> recover_G_from_poles_n_weights(dlr_obj& dlr_R0, nda::array<dcomplex,1> weights, std::vector<std::complex<double>> imfreqs);

double hubbard_dispersion(double kx, double ky);





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

template<typename T>
void triangle_to_square(std::vector<std::vector<T>>& M) {
    int n = M.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                M[j][i] = M[i][j];
            }
        }
    }
}	

template<typename T>
std::vector<std::vector<T>> data_to_full_bz(const std::vector<std::vector<T>>& M) {
    int n = M.size();
    int N = 2 * n - 1;
    int s_ind = n - 1;    
    int b_ind = N - 1;    
    std::vector<std::vector<T>> F(N, std::vector<T>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i > s_ind && j <= s_ind) F[i][j] = M[b_ind - i][j];
            else if (j > s_ind && i <= s_ind) F[i][j] = M[i][b_ind - j];
            else if (i > s_ind && j > s_ind) F[i][j] = M[b_ind - i][b_ind - j];
            else                               F[i][j] = M[i][j];
        }
    }
    return F;
}


	
template<typename T>
std::vector<T> sumVectors(const std::vector<std::vector<T>>& vecs) {
    if (vecs.empty()) return {};
    
    size_t N = vecs[0].size();
    // Verify all inner vectors have the same size
    for (const auto& v : vecs) {
        if (v.size() != N) {
            throw std::invalid_argument("All vectors must have the same length");
        }
    }
    
    std::vector<T> result(N, T{});
    for (const auto& v : vecs) {
        for (size_t i = 0; i < N; ++i) {
            result[i] += v[i];
        }
    }
    return result;
}





#endif // DLR_AMI_HPP