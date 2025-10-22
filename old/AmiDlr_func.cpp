#include "dlr_ami.hpp"
using namespace cppdlr;


dlr_obj create_dlr_obj(double beta, double eps, double Emax,AmiBase::g_struct R0_element) {
    dlr_obj dlro;
    dlro.ginfo = R0_element;
	//std::cout << " Created a DLR object with Epsilon and Alpha: ";
	//print1d(dlro.ginfo.eps_);print1d(dlro.ginfo.alpha_);
		
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
		//std::cout << all_poles[i] << std::endl;
		std::vector<double> tmp;
		tmp.reserve(eps_size);
		for (auto element : dlro.ginfo.eps_ ){ tmp.emplace_back(element* all_poles[i]);}///could negative
		dlro.evec.emplace_back(tmp);	
	}
	
	return dlro;
}


MPI_INFO create_MPI_obj(int total_num){
	MPI_INFO  MPI_obj;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &MPI_obj.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPI_obj.size);
	MPI_obj.chunk = total_num/MPI_obj.size;
	MPI_obj.remainder =  total_num% MPI_obj.size;
	MPI_obj.start = MPI_obj.rank*MPI_obj.chunk + std::min(MPI_obj.rank, MPI_obj.remainder);
	MPI_obj.count = MPI_obj.chunk + (MPI_obj.rank < MPI_obj.remainder ? 1 : 0);
	MPI_obj.end   = MPI_obj.start + MPI_obj.count; 
	
	auto s = std::format("created a DLR with rank={}, with start_index={}, and end_index={}, perfoming {}/{} of  Freq kernel  \n",MPI_obj.rank,
	MPI_obj.start, MPI_obj.end, MPI_obj.count, total_num);
	std::cout << s;
	
	return MPI_obj;
}
///////part of mdlr ///////////////////////////////////////////
mDLR::mDLR(double _beta,double _Uval, double _eps, double _Emax,size_t _kl,double _tp, AmiGraph::graph_t _graph,cppdlr::imfreq_ops _master_if_ops ):beta(_beta),Uval(_Uval),eps(_eps),Emax(_Emax),kl(_kl),tp (_tp), graph(_graph),master_if_ops(_master_if_ops){
	
	auto t0 = std::chrono::high_resolution_clock::now();
	g.label_systematic(graph);
	R0 = create_R0_from_graph();
	
	//R0 = construct_example2();
	N = R0.size();
	ord = g.graph_order(graph);
	prefactor =  g.get_prefactor(graph,ord);
	std::cout<<"-_-_-_-_-_-_-_-_-  Constructing multiple DLR Object for num(G):" << N << "  -_-_-_-_-_-_-_-_-  \n";
	create_multiple_gstruct();
	fill_dlro_pole_info();
	fill_dlro_momenta_info();
	MPI_obj =create_MPI_obj(total_num);
	std::cout<< "Done\n";
	std::cout << " Generating cartesian list from auxillary poless from G DLR rep \n";
	generate_cartesian_list();
	std::cout<<std::endl;
	CN = cartesian_combo_list.size();
	std::cout<<std::endl;

	dk = 2*M_PI/(kl);
	kvals.resize(kl);
	for(size_t i=0; i<kl; i++){kvals[i] = dk* double(i);}
	std::cout << "Operating on the grid\n " ;
	two_pi     = 2.0 * M_PI;
    inv_two_pi = 1.0 / two_pi;
    inv_dk     = 1.0 / dk;
	print1d(kvals);

	
	
	
	std::cout <<" Total number of auxillary epsilon_t required to be computed num(epsilon_t): " << CN <<std::endl;
	std::cout <<" Total number of cartesian momemta grid t required to be sampled Npoints: " << kN <<std::endl;
	std::cout << " Populating the auxillary energy lists from cartesian list \n";
	generate_auxillary_energy_list();
	std::cout << "done\n";
	std::cout<<std::endl;
	
	
	std::cout<<" creating master DLR object \n";
	//create_DLR_master_if_ops();
	std::cout<< " the matsubara frequency nodes of master DLR object is \n ";
	std::cout<< master_if_ops.get_ifnodes() << std::endl;
	master_pole_num = master_if_ops.get_ifnodes().size();
	master_poles = master_if_ops.get_rfnodes()/beta;
	fd_master_poles = fd_on_master_poles();
	std::cout << "real poles in master dlr is: \n";
	std::cout << master_poles;
   
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
		double new_Emax = Emax - static_cast<double>(i)/1.8;
		multiple_dlr_structs.push_back(create_dlr_obj(beta,eps, new_Emax, R0[i]));
	}
}

void mDLR::create_DLR_master_if_ops(){
	double E = 10; double eps=1e-10;
	double lambda = beta*E;
    auto dlr_rf = build_dlr_rf(lambda,eps );
    master_if_ops = imfreq_ops(lambda, dlr_rf, Fermion);	
}




void mDLR::generate_cartesian_list(){
	cartesian_combo_list.reserve(MPI_obj.count);
	for (int i=MPI_obj.start; i < MPI_obj.end;i++){
		std::vector<int> tmp;
		int previous = 1;
		for (int j =0;j<  num_pole_each_dlr.size();j++){
			tmp.push_back( i/previous %( num_pole_each_dlr[j] ));
			previous = previous*num_pole_each_dlr[j];    
		}
		cartesian_combo_list.emplace_back(tmp);		
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
	//print2d(auxillary_energy_list);
}


nda::array<dcomplex,1> mDLR::evaluate_auxillary_energies(nda::dcomplex &imfreq){
	nda::array<dcomplex,1> frequency_kernel(CN);
	AmiBase ami;
	AmiBase::frequency_t frequency;
	for(int i=0;i<ord;i++){ frequency.push_back(std::complex<double>(0,0));}
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
	int N_INT=ord;  // Number of Matsubara sums to perform
	AmiBase::ami_parms test_amiparms(N_INT, E_REG);

	//Construction Stage
	ami.construct(test_amiparms, R0, R_array, P_array, S_array);
	std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  external);
	//std::cout<<"Result was "<< calc_result<<std::endl;
	
	frequency_kernel(i) = dcomplex(-1,0)*calc_result;

	}
	return frequency_kernel;
	

}

nda::array<dcomplex,1> mDLR::generate_nda_Gdlr_from_energy( cppdlr::imfreq_ops &ops,
    double &energy)
{
    auto nodes = ops.get_ifnodes();                
    size_t N   = nodes.size();


    auto G = nda::zeros<dcomplex>( N);

    for (size_t i = 0; i < N; ++i) {
        dcomplex mf(0.0, (2*nodes(i) + 1) * M_PI / beta);
        G(i) = 1.0 / ( mf -dcomplex(energy, 0.0) );
    }
    return G;
}

void mDLR::populate_master_dlrW_from_G0(double mu){
	auto nodes = master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	for (int i =0;i< kl;i++){
		for (int j=0;j<kl;j++){
        double e = hubbard_dispersion(kvals[i],kvals[j],mu);
		
		auto master_gdlr = generate_nda_Gdlr_from_energy(master_if_ops,e);
		auto weights= master_if_ops.vals2coefs(beta,master_gdlr);
		master_dlrW_in_square[i][j]=weights;
		}	
	}	
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
			recovered_G(i) +=master_weights(j)/(iw - master_dlr_poles[j]);		
		}			
	}
	return recovered_G;
	
}


void mDLR::transfer_master_DLR_weights_to_dlrR0_elements(){
	for (auto &dlr_R0 : multiple_dlr_structs){
		auto wn_list = dlr_R0.im_freqs;
		dlr_R0.dlrW_in_square.clear();
		dlr_R0.dlrW_in_square.resize(kl);
		for (size_t i = 0; i < kl; ++i){dlr_R0.dlrW_in_square[i].resize(kl);}
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


Bz_container mDLR::vdot_freq_momenta_kernel_M(Bz_container mk, std::vector<nda::array<dcomplex,1>> fk){
	
	
	Bz_container result(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(fk.size())));
	for (int i =0;i<kl;i++){
		for (int j=0;j<kl;j++){
			auto const &momenta_kernel = mk[i][j];
			
			for (int k; k<fk.size();k++){
				result[i][j](k) = prefactor*std::pow(Uval,ord)*nda::dotc(fk[k],momenta_kernel)/(std::pow((double) kl*kl,ord));
			}			
		}	
	}
	return result;	
}




Bz_container mDLR::G_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu){
	Bz_container G(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(mfreq.size())));
	for (int i =0; i<kl; i++) {
		for (int j =0; j <kl; j++){
			double kx = kvals[i]; double ky = kvals[j];
			double e = hubbard_dispersion(kx,ky,mu);
			for (int f =0; f< mfreq.size();f++) {
				G[i][j](f) = 1/(mfreq(f) - e - SE[i][j](f));			
			}
		}	
	}	
	return G;	
}





Bz_container mDLR::G_from_DLR_SE_M_DMFT(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu){
	std::cout << "APPLYING DMFT SCHEMES\n";
	
	auto SE_loc = nda::zeros<dcomplex>(mfreq.size());
	auto G_loc = nda::zeros<dcomplex>(mfreq.size());
	
	////create local SE
	for (int f =0; f< mfreq.size();f++) {
		for (int i =0; i<kl; i++) {
			for (int j =0; j <kl; j++){
				SE_loc(f) +=  SE[i][j](f);			
			}
		}	
	}
	
	SE_loc = SE_loc/nda::dcomplex(kl*kl,0);
	
	
	for (int i =0; i < mfreq.size();i++){
		std::cout <<i<< " :" << mfreq(i) <<", " << SE_loc(i) << std::endl;
	}
	for (int f =0; f< mfreq.size();f++) {
		for (int i =0; i<kl; i++) {
			for (int j =0; j <kl; j++){
				double kx = kvals[i]; double ky = kvals[j];
				double e = hubbard_dispersion(kx,ky,mu);
				G_loc(f) +=  1/( mfreq(f) - e- SE[i][j](f));			
			}
		}
	}
	
	G_loc = G_loc/nda::dcomplex(kl*kl,0);
	
	
	
	
	Bz_container G(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(mfreq.size())));
	
	for (int i =0; i<kl; i++) {
		for (int j =0; j <kl; j++){
			double kx = kvals[i]; double ky = kvals[j];
			for (int f =0; f< mfreq.size();f++) {
					G[i][j](f)=	1.0/(1.0/G_loc(f) - SE_loc(f));
			}
		}	
	}
	
	return G;

}



Bz_container mDLR::G_from_DLR_SE_M_DCA(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu,double NC){
	if (MPI_obj.rank == 0){
	std::cout << " Applying  DCA patch averaging "<<std::endl;}
	Bz_container G(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(mfreq.size())));
	for (int i =0; i<kl; i++) {
		for (int j =0; j <kl; j++){
			double kx = kvals[i]; double ky = kvals[j];
			for (int f =0; f< mfreq.size();f++) {
					G[i][j](f)=	patch_avg_one_GF(kx,ky,NC,mu,mfreq[f],SE[i][j](f));
			}
		}	
	}	
	return G;	
}



nda::dcomplex  mDLR::patch_avg_one_GF(double Kx,double Ky, double Nc, double mu, nda::dcomplex iw, nda::dcomplex SE_K  ){
	double kx_top = Kx + dk/2;
	double kx_bottom = Kx-dk/2;
	double ky_top = Ky + dk/2;
	double ky_bottom = Ky-dk/2;

	double dK = dk/(Nc-1);
	auto GF_loc = nda::dcomplex(0,0);
	
	for (int i = 0; i < Nc; i++){
		for (int j = 0; j < Nc; j++){
			//std::cout << "kx : " << kx_bottom+ i*dK << " ky : "<< ky_bottom+j*dK <<std::endl;
			GF_loc += 1/(iw - hubbard_dispersion(kx_bottom+ i*dK,ky_bottom+j*dK,mu) - SE_K);	
		}
		
	}
	
	GF_loc = GF_loc/nda::dcomplex(Nc*Nc,0);
	
	nda::dcomplex g0_inv = 1/GF_loc + SE_K;
	
	
	
	return 1/g0_inv;
	


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

	


	 
	 
	 
	 
	 
	 
	 
	 
 


	






	




