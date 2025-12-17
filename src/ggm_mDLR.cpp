#include "dlr_ami.hpp"



ggm_mDLR::ggm_mDLR(const params_param& _params, AmiBase::graph_type _baseType,
                   const AmiGraph::gg_matrix_t& _ggm,
                   const cppdlr::imfreq_ops& _master_if_ops)
    : params(_params), baseType (_baseType), ggm(_ggm), master_if_ops(_master_if_ops)
{

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ggm_to_graphlist();
    graphlist_to_mDLRlist();
	graph_size = graphlist.size();
	auto master_nodes = master_if_ops.get_ifnodes();
	master_mfreq.resize(master_nodes.size());
	
	for (int i =0; i < master_nodes.size(); i++){
		master_mfreq(i) = nda::dcomplex( 0,(2*master_nodes(i)+1)*M_PI/params.beta);
	}
	if (rank==0){
	std::cout << " \n \n DLR master frequencies are :" << master_nodes <<std::endl;}
	reshape_Fk_ggm();
	
}
	
	
void ggm_mDLR::ggm_to_graphlist(){
	int gmin;int gmax;
    if (baseType == AmiBase::Sigma){
		gmin  = params.ord_min;
		gmax  = params.ord_max;
	}

	 else{
		gmin  = params.ord_min_2p;
		gmax  = params.ord_max_2p;
	}


	std::cout << " gmin is " << gmin << " , gmax is " << gmax;
	
	int count =0;
	for (int i = gmin; i < gmax+1; ++i){
	 for (int j= 0; j< ggm[i].size(); ++j){
		 for (int k=0; k <ggm[i][j].graph_vec.size();++k){
			 std::cout <<i <<" " << j << " " << k <<std::endl;
			 std::string name = "o"+ std::to_string(i) +"_g"+std::to_string(j) +"_n"+std::to_string(k);
			 AmiGraph::graph_t graph = ggm[i][j].graph_vec[k];
			 graphlist_names.push_back(name);
			 graphlist.push_back(graph);
			 std::cout<< "Sampling graph " << name << std::endl;
			 count++;
			}
		}
	}
	// std::cout<< " Loading completed \n";
}


void ggm_mDLR::graphlist_to_mDLRlist(){
	for (int i=0; i < graphlist.size(); i++){
		if (rank ==0){
		std::cout<<"-_-_-_-_-_-_-_-_-  Constructing multiple_DLR Object for :" <<  graphlist_names[i] << "  -_-_-_-_-_-_-_-_-  \n";}
		mDLR_list.emplace_back(params,baseType, graphlist[i], master_if_ops);			
	}	
}


void ggm_mDLR::reshape_Fk_ggm(){
	Fk_ggm.resize(graph_size);
	for (int i =0; i < graph_size; i++){
		 Fk_ggm[i] = nda::array<dcomplex, 2>(master_mfreq.size(), mDLR_list[i].MPI_obj.count);
	}
}


void ggm_mDLR::ggm_Fk_solver(){
	auto t0 = std::chrono::high_resolution_clock::now();
	for (int i =0; i < graph_size; i++){
	if (rank ==0){
	std::cout << " Computing Frequency kernel for the Graph: " << graphlist_names[i] << std::endl;}
	auto& mDLR = mDLR_list[i];
		for (int j=0; j < master_mfreq.size();j++){
			if (rank==0){std::cout<< "Computing -> " << j+1 <<"/" << master_mfreq.size() << " frequency point\n";}
				
				Fk_ggm[i](j,nda::range::all) = mDLR.evaluate_auxillary_energies(master_mfreq(j));
			}
			auto t1 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
			if (rank==0){
				std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";}
				t0 = t1;
		}
	}



	
	
	
Bz_container ggm_mDLR::generate_SE(mDLR &_mDLR, nda::array<dcomplex,2> &fk  ){                  
	auto momenta_kernel = _mDLR.compute_momenta_kernel_bz();
	if (rank ==0) {std::cout << "Momenta kernel computed \n";}
    return _mDLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, fk);
}


Bz_container ggm_mDLR::generate_summed_SE(int iter,bool write){
	std::vector<Bz_container> SE_list(graph_size);
	for (int i =0; i< graph_size; i++){
		auto SE =  generate_SE( mDLR_list[i], Fk_ggm[i]);
		mDLR_list[0].denoise_FS_points(SE);

		if (write){
		auto name = graphlist_names[i];
		std::string file_SE ="data/Individual_data/SE_" + name + "_i" + std::to_string(i) + ".txt";
		mDLR_list[i].write_data_momenta(file_SE,SE, master_mfreq);	
		SE_list[i] = std::move(SE);}
		
	}

	return sum_containers(SE_list);
}



void ggm_mDLR::intialize_ggm_DLR_W(){
	for (auto &mDLR : mDLR_list){
		mDLR.populate_master_dlrW();
	    mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
	}	
}

void ggm_mDLR::intialize_ggm_DLR_W(Bz_container &GF){
	for (auto &mDLR : mDLR_list){
		mDLR.populate_master_dlrW(GF);
	    mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
	}	
}


void ggm_mDLR::transfer_ggm_DLR_W(Bz_container &GF){
	for (auto &mDLR : mDLR_list){
		mDLR.repopulate_master_dlrW_from_G(GF);
		mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
	}	
}


void ggm_mDLR::ScPT_solver(int max_iters, Bz_container &final_SE, Bz_container &final_GF){
	double TOLERANCE = 1e-4;
	int    BISECT_STEPS    = 1000;
    double alpha=0.3;

	double mu =params.mu;
	Bz_container SE_old;
	Bz_container SE;
	Bz_container GF;
    int iter_idx =0;
	bool converged = false;
	while (iter_idx < max_iters ) {

	auto SE_new = generate_summed_SE(iter_idx,true);
	// if (params.mu ==0){    mDLR_list[0].denoise_FS_points(SE_new);}
	//mDLR_list[0].symmetrize_fermionic_DLR_Bz(SE_new);
	// 2) Mix old SE with new SE 
	SE = ( iter_idx < 1) ? SE_new : mDLR_list[0].SE_mixer(SE_old,SE_new, alpha);
	mDLR_list[0].denoise_FS_points(SE);
	mDLR_list[0].symmetrize_fermionic_DLR_Bz(SE);
	mDLR_list[0].denoise_FS_points(SE);
	
	double density = mDLR_list[0].compute_density_from_SE(SE, master_mfreq, mu);
	double mu_new = (std::abs(params.target_n - density) < 1e-3)
						? mu
						: mDLR_list[0].adjust_chemical_potential_bisc(params, SE, master_mfreq, BISECT_STEPS);

	double density_adj = mDLR_list[0].compute_density_from_SE(SE, master_mfreq, mu_new);
    // mu_new =0;
	if (rank == 0) {
			std::cout << "Iteration: " << iter_idx
					  << ", Initial-density = " << density <<" ,Adjusted-density = " << density_adj<< std::endl;
		}

	
	Bz_container GF;
	GF = mDLR_list[0].G_from_DLR_SE_M(SE, master_mfreq, mu_new);
	// if (params.DCA ==1){
	// 	GF = mDLR_list[0].G_from_DLR_SE_M_DCA(SE, master_mfreq, mu_new,params.patch_N);
	// 		}
	// else {
	// 	GF = mDLR_list[0].G_from_DLR_SE_M(SE, master_mfreq, mu_new);
	// 		}

	

	
	std::string file_SE = "data/Summed_data/" + std::to_string(iter_idx) + "i_shot_SE.txt";
    std::string file_GF = "data/Summed_data/" + std::to_string(iter_idx) + "i_shot_GF.txt";

	mDLR_list[0].write_data_momenta(file_SE, SE, master_mfreq);
	mDLR_list[0].write_data_momenta(file_GF, GF, master_mfreq);
	Bz_container GF_bar;
	Bz_container GF_weiss;
    if (params.DCA==1){
		GF_bar = mDLR_list[0].G_from_DLR_SE_M_DCA(SE, master_mfreq, mu_new,params.patch_N);
		GF_weiss = mDLR_list[0].generate_weiss_G(master_mfreq,SE, GF_bar);
		transfer_ggm_DLR_W(GF_weiss);
	}

	else {
		transfer_ggm_DLR_W(GF);
	}
   
	converged = SE_converged_abs(SE,
                      SE_old,
                      TOLERANCE);
	
	SE_old = SE;
	final_GF = GF;
	final_SE = SE;
	mu = mu_new;
	iter_idx++;
	}
	
	std::cout << "RAN for : " << iter_idx;

}


void ggm_mDLR::ggm_Fk_2p_solver(nda::array<dcomplex,1> &master_bfreq){
	auto t0 = std::chrono::high_resolution_clock::now();
	for (int i =0; i < graph_size; i++){
	if (rank ==0){
	std::cout << " Computing Frequency kernel for the Graph: " << graphlist_names[i] << std::endl;}
	auto& mDLR = mDLR_list[i];
		for (int j=0; j < master_bfreq.size();j++){
			if (rank==0){std::cout<< "Computing -> " << j+1 <<"/" << master_bfreq.size() << " frequency point\n" ;}
				
				Fk_ggm[i](j,nda::range::all) = mDLR.evaluate_auxillary_ph_energies(master_bfreq(j));
			}
			auto t1 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
			if (rank==0){
				std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";}
				t0 = t1;
		}
	}
Bz_container ggm_mDLR::generate_2p_bz(mDLR &_mDLR, nda::array<dcomplex,2> &fk ){                  
	auto momenta_kernel = _mDLR.compute_momenta_kernel_2p_bz();
	if (rank ==0) {std::cout << "Ph Momenta kernel computed \n";}
    return _mDLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, fk);
}


Bz_container ggm_mDLR::generate_summed_2p_bz(nda::array<dcomplex,1> &master_bfreq){
	std::vector<Bz_container> ph_list(graph_size);
	for (int i =0; i< graph_size; i++){
		auto ph =  generate_2p_bz( mDLR_list[i], Fk_ggm[i]);
		auto name = graphlist_names[i];
		std::string file_ph ="data_ph/Individual_data/ph_" + name + "_i" + std::to_string(i) + ".txt";
		mDLR_list[i].write_data_momenta(file_ph,ph,master_bfreq);	
	ph_list[i] = std::move(ph);
		
	}

	return sum_containers(ph_list);
}


void ggm_mDLR::single_shot_2p_solver_bz(nda::array<dcomplex,1> &master_bfreq,Bz_container &chi){
	chi = generate_summed_2p_bz(master_bfreq);
    std::string file_ph = "data_ph/Summed_data/single_shot_ph.txt";
	mDLR_list[0].write_data_momenta(file_ph, chi, master_bfreq);
}


Bz_container ggm_mDLR::generate_2p_qext(double qx, double qy, mDLR &_mDLR, nda::array<dcomplex,2> &fk ){                  
	auto momenta_kernel = _mDLR.compute_momenta_kernel_2p_qext(qx,qy);
	if (rank ==0) {std::cout << "Ph Momenta kernel computed \n";}
    return _mDLR.MPI_vdot_freq_momenta_kernel(momenta_kernel, fk);
}
Bz_container ggm_mDLR::generate_summed_2p_qext(double qx, double qy,nda::array<dcomplex,1> &master_bfreq){
	std::vector<Bz_container> ph_list(graph_size);
	for (int i =0; i< graph_size; i++){
		auto ph =  generate_2p_qext( qx, qy,mDLR_list[i], Fk_ggm[i]);
		auto name = graphlist_names[i];
		std::string file_ph ="data_ph/Individual_data/ph_" + name + "_i" + std::to_string(params.iter) + ".txt";
		mDLR_list[i].write_data_momenta(file_ph,ph,master_bfreq);	
		ph_list[i] = std::move(ph);
		
	}

	return sum_containers(ph_list);
}
void ggm_mDLR::single_shot_2p_solver_qext(double qx, double qy,nda::array<dcomplex,1> &master_bfreq,Bz_container &chi){
	chi = generate_summed_2p_qext( qx, qy,master_bfreq);
    std::string file_ph = "data_ph/Summed_data/single_shot_ph.txt";
	mDLR_list[0].write_data_momenta(file_ph, chi, master_bfreq);
	
}










	
	

	
	