#include "dlr_ami.hpp"
using namespace cppdlr;


ggm_mDLR::ggm_mDLR(const params_param& _params,
                   const AmiGraph::gg_matrix_t& _ggm,
                   const cppdlr::imfreq_ops& _master_if_ops)
    : params(_params), ggm(_ggm), master_if_ops(_master_if_ops)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::cout << "Initializing ggm_mDLR\n";
    // std::cout << "ggm top-level size = " << ggm.size() << std::endl;
    ggm_to_graphlist();
    // std::cout << "Constructing the mDLR objects for " << graphlist.size() << " graphs\n";
    graphlist_to_mDLRlist();
	graph_size = graphlist.size();
	auto master_nodes = master_if_ops.get_ifnodes();
	master_mfreq.resize(master_nodes.size());
	
	for (int i =0; i < master_nodes.size(); i++){
		master_mfreq(i) = nda::dcomplex( 0,(2*master_nodes(i)+1)*M_PI/params.beta);
	}
	
	// std::cout << " \n \n DLR master frequencies are :" << master_nodes <<std::endl;
	reshape_Fk_ggm();
	
}
	
	
void ggm_mDLR::ggm_to_graphlist(){
	// std::cout<< "Loading graphs \n";
	// std::cout << "the  ggm size is :"<< ggm.size();
	for (int i = params.ord_min; i < params.ord_max+1; ++i){
	 for (int j= 0; j< ggm[i].size(); ++j){
		 for (int k=0; k <ggm[i][j].graph_vec.size();++k){
			 // std::cout <<i <<" " << j << " " << k <<std::endl;
			 std::string name = "o"+ std::to_string(i) +"_g"+std::to_string(j) +"_n"+std::to_string(k);
			 AmiGraph::graph_t graph = ggm[i][j].graph_vec[k];
			 graphlist_names.push_back(name);
			 graphlist.push_back(graph);
			 // std::cout<< "Sampling graph " << name << std::endl;
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
		Fk_ggm[i].resize(master_mfreq.size());
		for (int j=0;j< master_mfreq.size();j++){
			Fk_ggm[i][j].resize(mDLR_list[i].MPI_obj.count);			
		}		
	}	
}


void ggm_mDLR::ggm_Fk_solver(){
	for (int i =0; i < graph_size; i++){
	if (rank ==0){
	std::cout << " Computing Frequency kernel for the Graph: " << graphlist_names[i] << std::endl;}
	auto& mDLR = mDLR_list[i];
		for (int j=0; j < master_mfreq.size();j++){
			if (rank==0){
				std::cout<< "Computing -> " << j <<"/" << master_mfreq.size() << " frequency point\n";}	
				auto frequency_kernel= mDLR.evaluate_auxillary_energies(master_mfreq(j)); 
				Fk_ggm[i][j] = std::move(frequency_kernel);
			}
			
		}
	}
	
void ggm_mDLR::intialize_ggm_DLR_W(){
	for (auto &mDLR : mDLR_list){
		mDLR.populate_master_dlrW_from_G0(params.mu);
	    mDLR.transfer_master_DLR_weights_to_dlrR0_elements();
	}	
}
	
	
Bz_container ggm_mDLR::generate_SE(mDLR &_mDLR, std::vector<nda::array<dcomplex,1>> &fk  ){                  
	auto momenta_kernel = _mDLR.compute_momenta_kernel_bz();
    return _mDLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, fk);
}

	
	

	
	