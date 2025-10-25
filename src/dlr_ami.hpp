#ifndef DLR_AMI_HPP
#define DLR_AMI_HPP
#include "ami_calc.hpp"
#include "amigraph.hpp"
#include <cassert>
#include <filesystem>
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

inline AmiGraph g(AmiBase::Sigma, 0);


//////params loader functions/////////////
using Bz_container =  std::vector<std::vector<nda::array<dcomplex,1>>>;
using Fk_container =  std::vector<nda::array<dcomplex,2>>;
Bz_container sum_containers(const std::vector<Bz_container>& all);
struct params_param {
	double Uval;
	double beta;
	double Emax;
	double eps;
	int iter;
	int L;
	double mu;
	double tp;
	double target_n;
	std::string graph;
	int ord_max;
	int ord_min;
	double mu_L;
	double mu_R;
	int DCA;
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



struct MPI_INFO {
	int rank,size,chunk,remainder,start,end,count;
};




class mDLR{
	private:
	double two_pi     ;
    double inv_two_pi ;
    double inv_dk     ;
    std::vector<int> kcombo_element;
	

	public:
	MPI_INFO MPI_obj;
	std::vector<dlr_obj> multiple_dlr_structs;
	double beta; double eps; double tp; double Emax;double Uval; AmiGraph::graph_t graph; AmiBase::g_prod_t R0; cppdlr::imfreq_ops master_if_ops;
	size_t N; ///num of greens function
	size_t CN; ///total number of cartesian, pole_num1* pole_num2* ...pole numN 
	size_t kl;///total number of momentum k grid;
	size_t kN;///total number of cartesian momenta, kl_1^2* kl_2^2.....
	double dk;
	double prefactor;
	int ord;
	int total_num=1;
	std::vector< int> num_pole_each_dlr;
	std::vector< int> num_k_each_dlr;

	
	std::vector<double> kvals;/// kgrid vals
	size_t master_pole_num;/// number of poles in master DLR
	nda::array<double,1> master_poles;
	nda::array<dcomplex,1> fd_master_poles;
	std::vector<std::vector<nda::array<dcomplex,1>>> master_dlrW_in_square; //master dlr weight
	mDLR(double _beta,double _Uval, double _eps, double _Emax,size_t _kl, double _tp, AmiGraph::graph_t _graph,cppdlr::imfreq_ops _master_if_ops );
	//////// methods ////////////
	
	AmiBase::g_prod_t create_R0_from_graph();
	void create_DLR_master_if_ops();
	double hubbard_dispersion(double kx, double ky,double mu);
	nda::array<dcomplex,1> generate_nda_Gdlr_from_energy( cppdlr::imfreq_ops &ops, double &energy);
	void fill_dlro_pole_info();
	void fill_dlro_momenta_info();
	void create_multiple_gstruct();
	// void generate_cartesian_list();
	std::vector<int> generate_single_CN(int index);
	// void generate_auxillary_energy_list();
	AmiBase::energy_t generate_auxillary_energy(std::vector<int> &combo) ;
	nda::array<dcomplex,1> evaluate_auxillary_energies(nda::dcomplex &imfreq);
	
	// nda::array<dcomplex,1> evaluate_auxillary_weights( nda::array<double,1> &energy);
	
	
	
	
	void populate_master_dlrW_from_G0(double mu);
	void reshape_dlrW_square_per_kgrid();
	nda::array<dcomplex,1> recover_dlro_G_from_master_weights(nda::array<dcomplex,1> &master_weights, std::vector<std::complex<double>> &dlro_if);
	void transfer_master_DLR_weights_to_dlrR0_elements();
	inline void generate_ith_momenta_cartesian_combo(int i, std::vector<int>& result );
	
	
	inline nda::dcomplex compute_momenta_one_kCN_kernel(double kx_ext,double ky_ext,std::vector<int> &combo,const std::vector<int> &kcombo);
	nda::array<nda::dcomplex,1> compute_momenta_kernel_qext(double kx_ext,double ky_ext);
	nda::dcomplex patch_avg_one_GF(double Kx,double Ky, double Nc, double mu, nda::dcomplex iw, nda::dcomplex SE_K  );	
    Bz_container G_from_DLR_SE_M_DCA(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu,double NC);
	Bz_container G_from_DLR_SE_M_DMFT(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu);
	Bz_container G_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu);
	Bz_container compute_momenta_kernel_bz();
	Bz_container vdot_freq_momenta_kernel_M( std::vector<std::vector<nda::array<dcomplex,1>>> &mk,  std::vector<nda::array<dcomplex,1>> &fk);
	Bz_container MPI_vdot_freq_momenta_kernel_M(Bz_container &mk, nda::array<dcomplex,2> &fk);
	nda::array<nda::dcomplex,1>  LocalG_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq,  double mu);
	Bz_container SE_mixer(Bz_container SE_old, Bz_container SE_new, double alpha);
	nda::array<dcomplex,1> fd_on_master_poles();
	double compute_density_from_SE(Bz_container &SE,nda::array<dcomplex,1> &mfreq, double mu);
	double adjust_chemical_potential_bisc(params_param &params, Bz_container &SE,nda::array<dcomplex,1> &mfreq, int max);
	double non_interacting_density(nda::array<dcomplex,1> &mfreq,double mu);
	void repopulate_master_dlrW_from_G(Bz_container &G );
	void write_data_momenta(const std::string& filename,Bz_container& data, nda::array<dcomplex,1>& mfreq);
	void write_data_ij_momenta(const std::string& filename,
                            Bz_container& data,
                            nda::array<dcomplex,1>& mfreq,
                            std::pair<int, int> ij);
							
							
};



class ggm_mDLR{
private:
int rank;

public:
 ggm_mDLR(const params_param& _params,
             const AmiGraph::gg_matrix_t& _ggm,
             const cppdlr::imfreq_ops& _master_if_ops);
params_param params; 
AmiGraph::gg_matrix_t ggm;
std::vector<AmiGraph::graph_t> graphlist;
std::vector<std::string> graphlist_names;
std::vector<mDLR> mDLR_list;
size_t graph_size;
cppdlr::imfreq_ops master_if_ops;
void ggm_to_graphlist();
void graphlist_to_mDLRlist();
nda::array<dcomplex,1> master_mfreq;
Fk_container  Fk_ggm;

void reshape_Fk_ggm();
void ggm_Fk_solver();
void intialize_ggm_DLR_W();
Bz_container generate_SE(mDLR &_mDLR, nda::array<dcomplex,2> &fk);	
	
};





dlr_obj create_dlr_obj(double beta, double eps, double Emax,AmiBase::g_struct R0_element);

MPI_INFO create_MPI_obj(int total_num);

std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec);
template<typename T>
void triangle_to_square(std::vector<std::vector<T>>& M);
template<typename T>
std::vector<std::vector<T>> data_to_full_bz(const std::vector<std::vector<T>>& M);
template<typename T>
std::vector<T> sumVectors(const std::vector<std::vector<T>>& vecs);
nda::array<dcomplex,1> recover_G_from_poles_n_weights(dlr_obj& dlr_R0, nda::array<dcomplex,1> weights, std::vector<std::complex<double>> imfreqs);


double fermi_distribution(double energy, double beta );





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