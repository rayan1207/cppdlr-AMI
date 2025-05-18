#include "dlr_ami.hpp"
using namespace cppdlr;



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
///////part of mdlr ///////////////////////////////////////////
mDLR::mDLR(double _beta,double _Uval, double _eps, double _Emax,size_t _kl,AmiBase::g_prod_t _R0):beta(_beta),Uval(_Uval),eps(_eps),Emax(_Emax),kl(_kl), R0(_R0){
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
	dk = 2*M_PI/(kl-1);
	std::cout <<" Total number of auxillary epsilon_t required to be computed num(epsilon_t): " << CN <<std::endl;
	std::cout <<" Total number of cartesian momemta grid t required to be sampled Npoints: " << kN <<std::endl;
	std::cout << " Populating the auxillary energy lists from cartesian list \n";
	generate_auxillary_energy_list();
	std::cout << "done\n";
	std::cout<<std::endl;
	kvals.resize(kl);
	for(size_t i=0; i<kl; i++){kvals[i] = dk* double(i);}
	std::cout << "Operating on the grid\n " ;
	print1d(kvals);
	std::cout<<" creating master DLR object \n";
	create_DLR_master_if_ops();
	std::cout<< " the matsubara frequency nodes of master DLR object is \n ";
	std::cout<< master_if_ops.get_ifnodes() << std::endl;
	master_pole_num = master_if_ops.get_ifnodes().size();

    two_pi     = 2.0 * M_PI;
    inv_two_pi = 1.0 / two_pi;
    inv_dk     = 1.0 / dk;
    internal_dof = N-1;
    kvals_ptr = kvals.data();
	master_dlrW_in_square.resize(kl);
	for (size_t i = 0; i < kl; ++i)
	master_dlrW_in_square[i].resize(kl);
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


nda::array<dcomplex,1> mDLR::evaluate_auxillary_energies(nda::dcomplex &imfreq){
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
	int N_INT=N-1;  // Number of Matsubara sums to perform
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);

	//Construction Stage
	ami.construct(test_amiparms, R0, R_array, P_array, S_array);
	std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  external);
	//std::cout<<"Result was "<< calc_result<<std::endl;
	
	frequency_kernel(i) = calc_result;

	}
	return frequency_kernel;
	

}

void mDLR::create_DLR_master_if_ops(){
	double lambda = beta*Emax;
    auto dlr_rf = build_dlr_rf(lambda,eps );
    master_if_ops = imfreq_ops(lambda, dlr_rf, Fermion);	
}

void mDLR::populate_master_dlrW_from_G0(){
	auto nodes = master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
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
		dlr_R0.dlrW_in_square.clear();
		dlr_R0.dlrW_in_square.resize(kl);
		for (size_t i = 0; i < kl; ++i){
		dlr_R0.dlrW_in_square[i].resize(kl);}
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




inline nda::dcomplex mDLR::compute_momenta_one_kCN_kernel(double kx_ext,double ky_ext,const int* combo_ptr,const int* kcombo_ptr)
{
    


    nda::dcomplex val{1.0, 0.0};

    for (int i = 0; i < N; ++i) {

        const auto& info      = multiple_dlr_structs[i].ginfo;
        const int* alpha   = info.alpha_.data();
        int   asz      = info.alpha_.size();
        double qx = alpha[asz - 1] * kx_ext;
        double qy = alpha[asz - 1] * ky_ext;
        for (int j = 0; j < internal_dof; ++j) {
            double a = static_cast<double>(alpha[j]);
            qx += a * kvals_ptr[ kcombo_ptr[2*j    ] ];
            qy += a * kvals_ptr[ kcombo_ptr[2*j + 1] ];
        }


        qx -= std::floor(qx * inv_two_pi) * two_pi;
        qy -= std::floor(qy * inv_two_pi) * two_pi;

        int idx1 = static_cast<int>(qx * inv_dk);
        int idx2 = static_cast<int>(qy * inv_dk);

     
        const auto& wlist =
          multiple_dlr_structs[i].dlrW_in_square[idx1][idx2];
        val = -val * wlist[ combo_ptr[i] ];
    }

    return val;
}
        
 nda::array<nda::dcomplex,1> mDLR::compute_momenta_kernel_qext(double kx_ext,double ky_ext)
{
    const int C = CN;
  
    nda::array<nda::dcomplex,1> kernel = nda::zeros<nda::dcomplex>(C);

    const int K = static_cast<int>(cartesian_k_combo_list.size());
    
    for (int c = 0; c < C; ++c) {
        const int* combo_ptr = cartesian_combo_list[c].data();
        auto sum = nda::dcomplex(0,0);

        for (int k = 0; k < K; ++k) {
            const int* kcombo_ptr = cartesian_k_combo_list[k].data();
            sum += compute_momenta_one_kCN_kernel(
                       kx_ext, ky_ext,
                       combo_ptr,
                       kcombo_ptr);
        }

        kernel(c) = sum;
    }

    return kernel;
}

Bz_container mDLR::compute_momenta_kernel_bz(){
	int n = (kl+1)/2;
	auto reduced_kgrid = nda::array<double,1> (n);
	Bz_container M(n,std::vector<nda::array<dcomplex,1>>(n));
	for(size_t i=0; i<n; i++){reduced_kgrid[i] = dk* double(i);}
	int data = n*(n+1)/2;
	int count =1;
	for (int i =0;i< n; i ++){
		for (int j=0;j<=i;j++){
			double qx = reduced_kgrid(i);
			double qy = reduced_kgrid(j);
			std::cout<< "Computing -> " << count <<"/" << data << " data point \n";
			auto momenta_kernel = mDLR::compute_momenta_kernel_qext(qx,qy);
			M[i][j] = momenta_kernel;
			count++;
		}		
	}
	triangle_to_square(M);
	
	return data_to_full_bz(M);
	 
}	
Bz_container mDLR::vdot_freq_momenta_kernel_M(Bz_container mk, std::vector<nda::array<dcomplex,1>> fk){
	
	
	Bz_container result(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(fk.size())));
	for (int i =0;i<kl;i++){
		for (int j=0;j<kl;j++){
			auto const &momenta_kernel = mk[i][j];
			
			for (int k; k<fk.size();k++){
				result[i][j](k) = std::pow(Uval,internal_dof)*nda::dotc(fk[k],momenta_kernel)/(-1*std::pow((double) kl,internal_dof*internal_dof));
			}			
		}	
	}
	return result;	
}


Bz_container mDLR::G_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq){
	Bz_container G(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(mfreq.size())));
	for (int i =0; i<kl; i++) {
		for (int j =0; j <kl; j++){
			double kx = kvals[i]; double ky = kvals[j];
			double e = hubbard_dispersion(kx,ky);
			for (int f =0; f< mfreq.size();f++) {
				G[i][j](f) = 1/(mfreq(f) - e - SE[i][j](f));			
			}
		}	
	}	
	return G;	
}

void mDLR::repopulate_master_dlrW_from_G(Bz_container &G ){
	master_dlrW_in_square.clear();
	master_dlrW_in_square.resize(kl);
	for (size_t i = 0; i < kl; ++i){
	master_dlrW_in_square[i].resize(kl);}
	
	for (int i =0;i< kl;i++){
		for (int j=0;j<kl;j++){
		auto master_gdlr = G[i][j];
		auto weights= master_if_ops.vals2coefs(beta,master_gdlr);
		master_dlrW_in_square[i][j]=weights;
		}	
	}	
}

	

	 
	 
	 
	 
	 
	 
	 
	 
	 
 


	






	




