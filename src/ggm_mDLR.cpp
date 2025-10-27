#include "dlr_ami.hpp"
using namespace cppdlr;


ggm_mDLR::ggm_mDLR(const params_param& _params,
                   const AmiGraph::gg_matrix_t& _ggm,
                   const cppdlr::imfreq_ops& _master_if_ops)
    : params(_params), ggm(_ggm), master_if_ops(_master_if_ops)
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
	// std::cout<< "Loading graphs \n";
	// std::cout << "the  ggm size is :"<< ggm.size();
	int count =0;
	for (int i = params.ord_min; i < params.ord_max+1; ++i){
	 for (int j= 0; j< ggm[i].size(); ++j){
		 for (int k=0; k <ggm[i][j].graph_vec.size();++k){
			 // std::cout <<i <<" " << j << " " << k <<std::endl;
			 std::string name = "o"+ std::to_string(i) +"_g"+std::to_string(j) +"_n"+std::to_string(k);
			 AmiGraph::graph_t graph = ggm[i][j].graph_vec[k];
			 graphlist_names.push_back(name);
			 graphlist.push_back(graph);
			 // std::cout<< "Sampling graph " << name << std::endl;
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
		mDLR_list.emplace_back(params.beta, params.Uval, params.eps, params.Emax,params.L,params.tp,graphlist[i], master_if_ops);			
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
	
void ggm_mDLR::intialize_ggm_DLR_W(){
	for (auto &mDLR : mDLR_list){
		mDLR.populate_master_dlrW_from_G0(params.mu);
	    mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
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

		if (write){
		auto name = graphlist_names[i];
		std::string file_SE = std::format("data/Individual_data/SE_{}_i{}.txt", name,iter);
		mDLR_list[i].write_data_momenta(file_SE,SE, master_mfreq);	
		SE_list[i] = std::move(SE);}
		
	}

	return sum_containers(SE_list);
}

void ggm_mDLR::transfer_ggm_DLR_W(Bz_container &GF){
	for (auto &mDLR : mDLR_list){
		mDLR.repopulate_master_dlrW_from_G(GF);
		mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
	}	
}


void ggm_mDLR::ScPT_solver(int max_iters){
	double DENSITY_TOLERANCE = 1e-4;
	int    BISECT_STEPS    = 1000;
    double alpha=0.3;
	double mu =params.mu;
	Bz_container SE_old;
	Bz_container SE;

	for (int iter_idx = 1; iter_idx < max_iters; ++iter_idx) {

	auto SE_new = generate_summed_SE(iter_idx,true);
	// 2) Mix old SE with new SE 
	SE = ( iter_idx < 2) ? SE_new : mDLR_list[0].SE_mixer(SE_old,SE_new, alpha);
	double density = mDLR_list[0].compute_density_from_SE(SE, master_mfreq, mu);
	double mu_new = (std::abs(params.target_n - density) < DENSITY_TOLERANCE)
						? mu
						: mDLR_list[0].adjust_chemical_potential_bisc(params, SE, master_mfreq, BISECT_STEPS);

	double density_adj = mDLR_list[0].compute_density_from_SE(SE, master_mfreq, mu_new);

	if (rank == 0) {
			std::cout << "Iteration: " << iter_idx
					  << ", Initial-density = " << density <<" ,Adjusted-density = " << density_adj<< std::endl;
		}

	
	Bz_container GF;
	if (params.DCA ==1){
		GF = mDLR_list[0].G_from_DLR_SE_M_DCA(SE, master_mfreq, mu_new,params.patch_N);
			}
	else {
		GF = mDLR_list[0].G_from_DLR_SE_M(SE, master_mfreq, mu_new);
			}

	std::string file_SE = std::format("data/Summed_data/{}i_shot_SE.txt", iter_idx);
	std::string file_GF = std::format("data/Summed_data/{}i_shot_GF.txt", iter_idx);

	mDLR_list[0].write_data_momenta(file_SE, SE, master_mfreq);
	mDLR_list[0].write_data_momenta(file_GF, GF, master_mfreq);

    transfer_ggm_DLR_W(GF);
	SE_old = SE;
	mu = mu_new;
	}
}
	
	

	
	