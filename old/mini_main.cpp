


#include "dlr_ami.hpp"
using namespace cppdlr;

int main(int argc, char** argv){
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	

	params_param params;
	std::string loader_file = "../loader/params.txt";
	params_loader(loader_file,params);
	
	
	////////////   Load the graphs here //////////////////
	
	AmiGraph::gg_matrix_t ggm;
	g.read_ggmp(params.graph,ggm, params.ord_max);
	std::cout<<std::endl;
	
	
	
	
	AmiGraph::graph_t graph =ggm[params.ord_max][0].graph_vec[0];
	
	
	
	double NC = 161;
	double beta   = params.beta;
    double eps    = params.eps;
	int kl = params.L;
	double Emax = params.Emax;
	double Uval = params.Uval;
	double lambda = beta*Emax;
	double mu = params.mu;
	double tp = params.tp;
	int iter = params.iter;
	AmiBase ami;
	//AmiBase::g_prod_t R0=construct_example2();
	
	double master_E = 5; double master_eps=1e-7;
	double master_lambda = params.beta*master_E;
    auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
    auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);
	
	mDLR multiple_DLR(beta,Uval,eps,Emax,kl,tp,graph,master_if_ops);
	
	std::cout <<" Testing " <<std::endl;
	ggm_mDLR  mult_mDLR( params, ggm);
   
	//////// Computing the frequency kernel ////////////////
	if (rank ==0){
	std::cout << "--__--__--__--__--__--__--__--__--__--__--"<< std::endl;
	std::cout <<" Precomputing Computing the frequency kernel \n";
	}
	auto t0 = std::chrono::high_resolution_clock::now();
	auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	
	
	for (size_t i =0; i<nodes.size();i++){ 
	mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
	}
	
	
	std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	
	for (int i=0; i <nodes.size();i++){
	auto val = mfreq(i);
	if (rank ==0){
	std::cout << val <<std::endl;
	}
	auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(val); 
	frequency_kernel_list.push_back(frequency_kernel);
	
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	if (rank==0){
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
		std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";
	}
	
	
	multiple_DLR.populate_master_dlrW_from_G0(mu);
	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	
	
	double DENSITY_TOLERANCE = 1e-4;
	int    BISECT_STEPS    = 1000;
    double alpha=0.3;
	int max_iters = iter; 
	auto density_lst = nda::array<double,1>(iter);
    Bz_container SE_old;
	Bz_container SE;
	for (int iter_idx = 1; iter_idx < max_iters; ++iter_idx) {
		// 1) Build the kernel and compute self energy
		auto momenta_kernel = multiple_DLR.compute_momenta_kernel_bz();
		auto SE_new = multiple_DLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, frequency_kernel_list);
		// 2) Mix old SE with new SE 
		
		SE = ( iter_idx < 2) ? SE_new : multiple_DLR.SE_mixer(SE_old,SE_new, alpha);

		// 3) Compute density at current mu, adjust mu if needed
		double density = multiple_DLR.compute_density_from_SE(SE, mfreq, mu);
		double mu_new = (std::abs(params.target_n - density) < DENSITY_TOLERANCE)
						? mu
						: multiple_DLR.adjust_chemical_potential_bisc(params, SE, mfreq, BISECT_STEPS);

		double density_adj = multiple_DLR.compute_density_from_SE(SE, mfreq, mu_new);
		density_lst[iter_idx - 1] = density_adj;

		if (rank == 0) {
			std::cout << "Iteration " << iter_idx
					  << ": density = " << density_adj << std::endl;
		}

		// 4) Compute Green’s function and write out data
		Bz_container GF;
		if (params.DCA ==1){
		  GF = multiple_DLR.G_from_DLR_SE_M_DCA(SE, mfreq, mu_new,NC);
		}
		else {
		  GF = multiple_DLR.G_from_DLR_SE_M(SE, mfreq, mu_new);
		}
		std::string file_SE = std::format("../result/{}i_shot_SE.txt", iter_idx);
		std::string file_GF = std::format("../result/{}i_shot_GF.txt", iter_idx);

		multiple_DLR.write_data_momenta(file_SE, SE, mfreq);
		multiple_DLR.write_data_momenta(file_GF, GF, mfreq);

		// 5) Update master DLR weights for next iteration
		multiple_DLR.repopulate_master_dlrW_from_G(GF);
		multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
		// 6) Old SE back to current SE
		SE_old = SE;
	}

		
		
	

	MPI_Finalize();

	
}



// int main () {
	
	// params_param params;
	// std::string loader_file = "../loader/params.txt";
	// params_loader(loader_file,params);
	
	
	// ////////////   Load the graphs here //////////////////
	
	// AmiGraph::gg_matrix_t ggm;
	// std::cout<<"Attempting to load self-energy graphs from example_graphs"<<std::endl;
	// g.read_ggmp(params.graph,ggm, params.ord_max);
	// std::cout<<"Completed read"<<std::endl;
	// std::cout<<std::endl;
	
	
	
	// for (int i =  params.ord_min; i <  params.ord_max+1;i++){
		// for (int j =0; j< ggm[i].size();j++){
			 // for (int k = 0; k < ggm[i][j].graph_vec.size(); ++k) {
				// std::cout << "Labeling graph with " << "o" << i << "_g" << j << "_n" << k <<std::endl;
				// g.label_systematic(ggm[i][j].graph_vec[k]);
			 // }
		// }
	// }

	// AmiGraph::graph_t graph =ggm[2][0].graph_vec[0];
	
	
	
	// double NC = 161;
	// double beta   = params.beta;
    // double eps    = params.eps;
	// int kl = params.L;
	// double Emax = params.Emax;
	// double Uval = params.Uval;
	// double lambda = beta*Emax;
	// double mu = params.mu;
	// double tp = params.tp;
	// AmiBase ami;
	// //AmiBase::g_prod_t R0=construct_example2();
	
	// mDLR multiple_DLR(beta,Uval,eps,Emax,kl,tp,graph);
   
	
	
    
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
// std::string loader_file = "../loader/params.txt";
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






































// #include <mpi.h>
// #include <iostream>
// #include <vector>
// #include <cppdlr/cppdlr.hpp>

// template<typename T>
 // void print2d( std::vector< std::vector<T>> vec)
// {
    // for ( auto row : vec) {
        // for ( auto elem : row) {
            // std::cout << elem << " ";
        // }
        // std::cout << std::endl;
    // }
	// std::cout << std::endl;
	
// }


// template<typename T> 
 // void print1d( std::vector<T>& vec) {
  // std::cout << "[";
  // for (size_t i = 0; i < vec.size(); ++i) {
    // std::cout << vec[i];
    // if (i != vec.size() - 1) {
      // std::cout << ", ";
    // }
  // }
  // std::cout << "]\n";
// }

// std::vector<std::vector<int>> generate_cartesian_list(){
	// int total_num=3*4*5;
	// int rank, size;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
	// int chunk = total_num/size ;
	// int remainder =  total_num % size;
	// int start = rank*chunk + std::min(rank, remainder);
	// int count = chunk + (rank < remainder ? 1 : 0);
	// int end   = start + count; 
	
	
	// std::vector<int> num_pole_each_dlr = {3,4,5};
    // std::vector<std::vector<int>> cartesian_combo_list;
	// cartesian_combo_list.reserve(total_num);
	
	
	// int i =0;
	// for (int i=start; i < end;i++){
		// std::vector<int> tmp;
		// int previous = 1;
		// for (int j =0;j<  num_pole_each_dlr.size();j++){
			// tmp.push_back( i/previous %( num_pole_each_dlr[j] ));
			// previous = previous*num_pole_each_dlr[j];    
		// }
		// cartesian_combo_list.emplace_back(tmp);	
		
	// }
	// return cartesian_combo_list;
// }

// int main(int argc, char** argv) {
    // MPI_Init(&argc, &argv);
    // int rank, size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	
	
	// int total_num=3*4*5;
	// int chunk = total_num/size ;
	// int remainder =  total_num % size;
	
	// int start = rank*chunk + std::min(rank, remainder);
	// int count = chunk + (rank < remainder ? 1 : 0);
	// int end   = start + count;  
	

	
	// std::cout << "I am rank : " <<  rank  << "with start :" << start  << " end :" << end   <<std::endl;
	// auto combos = generate_cartesian_list();
	// print2d(combos);
	

	// MPI_Finalize();

// 