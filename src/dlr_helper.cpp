#include "dlr_ami.hpp"
using namespace cppdlr;
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


dlr_obj create_dlr_obj(double beta, double eps, double Emax,AmiBase::g_struct R0_element) {
    dlr_obj dlro;
    dlro.ginfo = R0_element;
	std::cout << " Created a DLR object with Epsilon and Alpha: ";
	print1d(dlro.ginfo.eps_);print1d(dlro.ginfo.alpha_);
		
	double lambda = beta*Emax;
    auto dlr_rf = build_dlr_rf(lambda, eps);
	dlro.if_ops = imfreq_ops(lambda, dlr_rf, Fermion); //create an ifops objects inside dlrp
	auto if_nodes = dlro.if_ops.get_ifnodes();
	for (int i=0; i < if_nodes.size();i++){  //// filling in fermionic matsubara frequncy
		auto dlr_n = if_nodes(i);
		double wn = (2*dlr_n +1)*M_PI/beta;
		dlro.im_freqs.push_back(std::complex<double>(0,wn));
	}
	auto all_poles = dlro.if_ops.get_rfnodes()/beta;
	dlro.pole_num = all_poles.size();
	int eps_size = dlro.ginfo.eps_.size();
	for (int i = 0; i< dlro.pole_num; i++){          /// filling in pole locations 
		dlro.pole_locs.emplace_back(all_poles[i]);	
		std::vector<double> tmp;
		tmp.reserve(eps_size);
		for (auto element : dlro.ginfo.eps_ ){ tmp.emplace_back( -1*element* all_poles[i]);}
		dlro.evec.emplace_back(tmp);	
	}
	
	std::cout<< " The pole weights are :";
	print2d(dlro.evec);
	return dlro;
}








///////part of mdlr 
mDLR::mDLR(double _beta, double _eps, double _Emax,size_t _kl,AmiBase::g_prod_t _R0):beta(_beta),eps(_eps),Emax(_Emax),kl(_kl), R0(_R0){
	auto t0 = std::chrono::high_resolution_clock::now();
	N = R0.size();
	std::cout<<"-_-_-_-_-_-_-_-_-  Constructing multiple DLR Object for num(G):" << N << "  -_-_-_-_-_-_-_-_-  \n";
	create_multiple_gstruct();
	std::cout<< "Done\n";
	std::cout << " Generating cartesian list from auxillary poless from G DLR rep \n";
	generate_cartesian_list();
	std::cout<<std::endl;
	CN = cartesian_combo_list.size();
	std::cout<<std::endl;
	std::cout <<" Generating cartesian momenta grid for flat momenta sampling \n";
	generate_momenta_cartesian_combo();
	kN = cartesian_k_combo_list.size();
	
	std::cout <<" Total number of auxillary epsilon_t required to be computed num(epsilon_t): " << CN <<std::endl;
	std::cout <<" Total number of cartesian momemta grid t required to be sampled Npoints: " << kN <<std::endl;
	std::cout << " Populating the auxillary energy lists from cartesian list \n";
	generate_auxillary_energy_list();
	std::cout << "done\n";
	std::cout<<std::endl;
	kvals.resize(kl);
	for(size_t i=0; i<kl; i++){kvals[i] = 2*M_PI * double(i) / double(kl - 1);}
	std::cout << "Operating on the grid\n " ;
	print1d(kvals);
	std::cout<<" creating master DLR object \n";
	create_DLR_master_if_ops();
	std::cout<< " the matsubara frequency nodes of master DLR object is \n ";
	std::cout<< master_if_ops.get_ifnodes() << std::endl;
	master_pole_num = master_if_ops.get_ifnodes().size();
	
	
	master_dlrW_in_square.resize(kl);
	for (size_t i = 0; i < kl; ++i)
	master_dlrW_in_square[i].resize(kl);


	std::cout <<"Populating the 2-dim master_dlrW_in_square grid with weights from G0 with L: " << kl <<std::endl;
	populate_master_dlrW_from_G0();
	std::cout << "Done";
	reshape_dlrW_square_per_kgrid();
	transfer_master_DLR_weights_to_dlrR0_elements();
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	std::cout << " Construction of mDLR took time: " <<duration.count() << " ms \n";
}


void mDLR::create_multiple_gstruct(){
	for (int i =0; i< R0.size(); i++){
		double new_Emax = Emax - static_cast<double>(i)/2.0;
		multiple_dlr_structs.push_back(create_dlr_obj(beta,eps, new_Emax, R0[i]));
	}
}



void mDLR::generate_cartesian_list(){
	int total_num=1;
	int num_dlr = multiple_dlr_structs.size();
	std::vector<int> num_pole_each_dlr;
	
	for (auto dlr_R0 : multiple_dlr_structs){
		total_num =total_num*dlr_R0.pole_num;
		num_pole_each_dlr.push_back(dlr_R0.pole_num);	
	}
	// std::cout<<"pole num for each dlr";
	// print1d(num_pole_each_dlr);
	cartesian_combo_list.reserve(total_num);
	int i =0;
	while ( i < total_num ){
		std::vector<int> tmp;
		int previous = 1;
		for (int j =0;j<  num_pole_each_dlr.size();j++){
			tmp.push_back( i/previous %( num_pole_each_dlr[j] ));
			previous = previous*num_pole_each_dlr[j];    
		}
		cartesian_combo_list.emplace_back(tmp);	
		i++;
	}		
}


void mDLR::generate_auxillary_energy_list(){
	for (auto const& combo : cartesian_combo_list){
		std::vector<std::vector<double>> tmp;
		for (int i = 0; i< N;i++){
			tmp.emplace_back(multiple_dlr_structs[i].evec[combo[i]]);		
		}
		auxillary_energy_list.push_back(sumVectors(tmp));
	}
}


nda::array<dcomplex,1> mDLR::evaluate_auxillary_energies(std::complex<double> &imfreq){
	nda::array<dcomplex,1> frequency_kernel(CN);
	AmiBase ami;
	AmiBase::frequency_t frequency;
	for(int i=0;i<N-1;i++){ frequency.push_back(std::complex<double>(0,0));}
	frequency.push_back(imfreq);
	
	
	for (int i =0; i<CN;i++){
		AmiBase::energy_t energy =  convertToComplex(auxillary_energy_list[i]);
		AmiBase::ami_vars external(energy, frequency,beta);
	
	// Storage objects for S,P,R 
	AmiBase::S_t S_array;
	AmiBase::P_t P_array;
	AmiBase::R_t R_array;

	// Integration/Evaluation parameters
	double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
	int N_INT=2;  // Number of Matsubara sums to perform
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);

	//Construction Stage
	ami.construct(test_amiparms, R0, R_array, P_array, S_array);
	std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  external);
	//std::cout<<"Result was "<< calc_result<<std::endl;
	
	frequency_kernel(i) = calc_result;

	}
	return frequency_kernel;
	

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


void mDLR::create_DLR_master_if_ops(){
	double lambda = beta*Emax;
    auto dlr_rf = build_dlr_rf(lambda, 1e-10 );
    master_if_ops = imfreq_ops(lambda, dlr_rf, Fermion);	
}

void mDLR::populate_master_dlrW_from_G0(){
	auto nodes = master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	for (size_t i =0; i<nodes.size();i++){ mfreq[i]=dcomplex(0,(2*nodes[i]+1)/beta);}
	
	for (int i =0;i< kl;i++){
		for (int j=0;j<kl;j++){
        double e = hubbard_dispersion(kvals[i],kvals[j]);
		
		auto master_gdlr = generate_nda_Gdlr_from_energy(master_if_ops,e);
		auto weights= master_if_ops.vals2coefs(beta,master_gdlr);
		master_dlrW_in_square[i][j]=weights;
		}	
	}	
}

nda::array<dcomplex,1> mDLR::generate_nda_Gdlr_from_energy( cppdlr::imfreq_ops &ops,
    double &energy)
{
    auto nodes = ops.get_ifnodes();                
    size_t N   = nodes.size();


    auto G = nda::zeros<dcomplex>( N);

    for (size_t i = 0; i < N; ++i) {
        dcomplex mf(0.0, (2*nodes(i) + 1) * M_PI / beta);
        G(i) = 1.0 / ( mf - dcomplex(energy, 0.0) );
    }
    return G;
}

void mDLR::reshape_dlrW_square_per_kgrid(){
	for (int i =0; i < N;i++){
		multiple_dlr_structs[i].dlrW_in_square.resize(kl);
		for (size_t j = 0; j < kl; j++){multiple_dlr_structs[i].dlrW_in_square[j].resize(kl);}
				
	}
}


nda::array<dcomplex,1> mDLR::recover_dlro_G_from_master_weights(nda::array<dcomplex,1> &master_weights, std::vector<std::complex<double>> &dlro_if){
	int size =dlro_if.size();
	auto recovered_G = nda::zeros<dcomplex> (size);
	auto master_dlr_poles = master_if_ops.get_rfnodes()/beta;
	
	for (int i=0; i < size;i++){
		auto iw = dlro_if[i];
		for (int j=0; j < master_weights.size(); j++){
			recovered_G(i) += master_weights(j)/(iw - master_dlr_poles[j]);		
		}			
	}
	return recovered_G;
	
}


void mDLR::transfer_master_DLR_weights_to_dlrR0_elements(){
	for (auto &dlr_R0 : multiple_dlr_structs){
		auto wn_list = dlr_R0.im_freqs;
		//dlr_R0.dlrW_in_square.clear();
		for (int i =0; i <kl;i++){
			for (int j =0; j < kl;j++){
				auto master_weights = master_dlrW_in_square[i][j];
				auto dlr_R0_g = recover_dlro_G_from_master_weights(master_weights, wn_list);
				dlr_R0.dlrW_in_square[i][j] = dlr_R0.if_ops.vals2coefs(beta, dlr_R0_g);
				//std::cout << dlr_R0.dlrW_in_square[i][j] <<std::endl;
			}	
		}	
	}
}


void mDLR::generate_momenta_cartesian_combo(){
	
	int kC_num =  std::pow(std::pow(kl,2),(N-1));

	std::vector<int> num_k_each_dlr;
	
	for (int i = 0; i< 2*(N-1);i++){
		num_k_each_dlr.push_back(kl);	
	}
	print1d(num_k_each_dlr);
	
	cartesian_k_combo_list.reserve(kC_num);
	int i =0;
	while ( i < kC_num ){
		std::vector<int> tmp;
		int previous = 1;
		for (int j =0;j<  num_k_each_dlr.size();j++){
			tmp.push_back( i/previous %( num_k_each_dlr[j] ));
			previous = previous*num_k_each_dlr[j];    
		}
		cartesian_k_combo_list.emplace_back(tmp);	
		i++;
	}
	
}



	






	




