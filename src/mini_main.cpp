#include "dlr_ami.hpp"
using namespace cppdlr;


int main(){
	params_param params;
	std::string loader_file = "params.txt";
	params_loader(loader_file,params);
	
	
	////////////   Load the graphs here //////////////////
	
	AmiGraph::gg_matrix_t ggm;
	std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	g.read_ggmp(params.graph,ggm, params.ord_max);
	std::cout<<"Completed read"<<std::endl;
	std::cout<<std::endl;
	g.ggm_label(ggm,0); 
    AmiGraph::graph_t graph = ggm[2][0].graph_vec[0];	
	

	
	
	
	
	
	
	
	
	
	
	
	
	double NC = 161;
	double beta   = params.beta;
    double eps    = params.eps;
	int kl = params.L;
	double Emax = params.Emax;
	double Uval = params.Uval;
	double lambda = beta*Emax;
	double mu = params.mu;
	double tp = params.tp;
	AmiBase ami;
	//AmiBase::g_prod_t R0=construct_example2();
	
	mDLR multiple_DLR(beta,Uval,eps,Emax,kl,tp,graph);
    
	//////// Computing the frequency kernel ////////////////
	std::cout << "--__--__--__--__--__--__--__--__--__--__--"<< std::endl;
	std::cout <<" Precomputing Computing the frequency kernel \n";
	auto t0 = std::chrono::high_resolution_clock::now();
	
	
	auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	
	
	for (size_t i =0; i<nodes.size();i++){ 
	mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
	}
	
	
	std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	
	for (int i=0; i <nodes.size();i++){
	auto val = mfreq(i);
	std::cout << val <<std::endl;
	auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(val); 
	frequency_kernel_list.push_back(frequency_kernel);
	}
	
	
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	std::cout << " Duration of frequency kernel computation: " <<duration.count() << " ms \n";
	


    std::cout << "--__--__--__--__--__--__--__--__--__--__-- "<<std::endl;
	std::cout << "Preparing for momneta kernel for the first shot" <<std::endl;
	std::cout <<"Populating the 2-dim master_dlrW_in_square grid with weights from G0 with L: " << kl <<std::endl;
	
	multiple_DLR.populate_master_dlrW_from_G0(mu);
	
	std::cout << "Done";
	std::cout << "Transferring DLR weights to gstructs \n";

	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	
	std::cout <<"Done \n\n";
	double U0_density = multiple_DLR.non_interacting_density(mfreq,params.mu);
	std::cout << "Non interacting density n= " << U0_density <<" \n\n";
	params.target_n = U0_density;
    
	int iter = params.iter;
    auto density_lst = nda::array<double,1>(iter-1);
	

    for (int i=1;i<iter;i++){
		std::cout << "--__--__--__--__--__--__--__--__--__--__-- "<<std::endl;
		
		
		std::cout << "computing for momneta kernel for the "<< i <<" shot" <<std::endl;
		
		auto momenta_kernel_bz =multiple_DLR.compute_momenta_kernel_bz();
		auto SE =  multiple_DLR.vdot_freq_momenta_kernel_M(momenta_kernel_bz,frequency_kernel_list);
	
	
		std:: cout << "Target denisty is n= " << params.target_n << std::endl;
		
		 double mu_new = multiple_DLR.adjust_chemical_potential_bisc(params,SE,mfreq,1000);
		
		 auto density = multiple_DLR.compute_density_from_SE(SE,mfreq,mu_new);
		 std::cout<< "Computed density in this iteration n: " << density << std::endl;
		 density_lst(i-1)=density;
		 
		auto GF = multiple_DLR.G_from_DLR_SE_M(SE, mfreq,mu_new);
		
		
		std::string filename_SE = std::format("gf2_data/{}i_shot_SE.txt", i);
		std::string filename_GF = std::format("gf2_data/{}i_shot_GF.txt", i);
		
		std::cout <<   "Writing SE data to " << filename_SE <<std::endl;
		
		multiple_DLR.write_data_momenta(filename_SE,SE,mfreq);
		
		std::cout <<   "Writing GF data to " << filename_GF <<std::endl;
		
		multiple_DLR.write_data_momenta(filename_GF,GF,mfreq);
		
		std::cout << " Repopulating the master DLR with new G" <<std::endl;
	    
		multiple_DLR.repopulate_master_dlrW_from_G(GF);
		
		std::cout << " Tranferring DLR weights to all dlr R0 elements for next iteration" <<std::endl;
		
		multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
		std::cout << "done \n";
       
	}
	
	std::cout<< "Density in each list from first to last \n";
	std::cout << density_lst;

	
}

// int main () {
	
	// params_param params;
	// std::string loader_file = "params.txt";
	// params_loader(loader_file,params);
	
	// double beta   = params.beta;
    // double eps    = params.eps;
	// int kl = params.L;
	// double Emax = params.Emax;
	// double Uval = params.Uval;
	// double lambda = beta*Emax;
	// double mu = params.mu;
	// double tp = params.tp;
 
	// AmiBase ami;
	// AmiBase::g_prod_t R0=construct_example2();
	
	// mDLR multiple_DLR(beta,Uval,eps,Emax,kl,tp,R0);
    
	// Computing the frequency kernel ////////////////
	// std::cout << "--__--__--__--__--__--__--__--__--__--__--"<< std::endl;
	// std::cout <<" Precomputing Computing the frequency kernel \n";
	// auto t0 = std::chrono::high_resolution_clock::now();
	
	
	// auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
	// nda::array<dcomplex,1> mfreq(nodes.size());
	
	// for (size_t i =0; i<nodes.size();i++){ 
	// mfreq(i)=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
	// }
	
	
	
	
	// std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	
	// for (int i=0; i <nodes.size();i++){
	// auto val = mfreq(i);
	// std::cout << val <<std::endl;
	// auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(val); 
	// frequency_kernel_list.push_back(frequency_kernel);
	// }
	
	// nda::array<double,1> energy = {-4,0.1,-1};
	// energy =  energy;
	
	// auto k_kernel = multiple_DLR.evaluate_auxillary_weights(energy);
	
	// for (int i=0; i <nodes.size();i++){
		// auto fk = frequency_kernel_list[i];
		// auto result = nda::dotc(fk,k_kernel);	
		// std::cout << "Mfreq = " << mfreq(i) << "result =" << result << std::endl;
	// }

// }


// int main() {
// params_param params;
// std::string loader_file = "params.txt";
// params_loader(loader_file,params);



	
// AmiGraph g(AmiBase::Sigma, 0);	
// AmiGraph::gg_matrix_t ggm;
// std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;

// g.read_ggmp(params.graph,ggm, params.ord_max);
// std::cout<<"Completed read"<<std::endl;
// std::cout<<std::endl;
// g.ggm_label(ggm,0);  
	
// std::vector<AmiBase::g_prod_t>  R0_Lists; 


// for (int i =  params.ord_min; i <  params.ord_max+1;i++){
	// for (int j =0; j< ggm[i].size();j++){
		 // for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
		// std::cout << "sampling graph with " << "o" << i << "_g" << j << "_n" << k <<std::endl;
		// AmiGraph::edge_vector_t internal_fedges_list;
	    // AmiGraph::graph_t graph = ggm[i][j].graph_vec[k];
		// g.find_internal_fermionic_edges(graph,internal_fedges_list);
		// int size = internal_fedges_list.size();
		// AmiGraph::g_struct_ gprod[size];
			// for (auto fedge : internal_fedges_list){
				// std::cout <<" Eps and Alpha vectors for each are :";
				// print1d(graph[fedge].g_struct_.alpha_);
				// std::cout << " and "; 
				// print1d(graph[fedge].g_struct_.eps_);
				// gprod
				
				// }
		// }
		
	// }

// }





// }