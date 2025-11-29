


#include "dlr_ami.hpp"
using namespace cppdlr;

int main(int argc, char** argv){
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (rank==0){
	fs::remove_all("data");
    std::string folder_name1 = "data/Summed_data";
	std::string folder_name2 = "data/Individual_data";
	fs::create_directories(folder_name1);
	fs::create_directories(folder_name2);
}





	params_param params;
	std::string loader_file = "params.txt";
	params_loader(loader_file,params);
	
	
	////////////   Load the graphs here //////////////////
	AmiBase::graph_type sigma=AmiBase::Sigma;
	AmiGraph gs(sigma, 0);
	AmiGraph::gg_matrix_t ggm;
	gs.read_ggmp(params.graph,ggm, params.ord_max);
	std::cout<<std::endl;
	double beta   = params.beta;
    double eps    = params.eps;
	int kl = params.L;
	double Emax = params.Emax;
	double Uval = params.Uval;
	double lambda = beta*Emax;
	double mu = params.mu;
	double tp = params.tp;
	int iter = params.iter;
	


	
	
	double master_E = 10; double master_eps=1e-10;
	double master_lambda = params.beta*master_E;
    auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
    auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);
	auto master_bf_ops = imfreq_ops(master_lambda, dlr_rf, Boson);

   	Bz_container scGF;
   	Bz_container scSE;

	ggm_mDLR  mult_mDLR( params,sigma, ggm,master_if_ops);
	mult_mDLR.ggm_Fk_solver();
	mult_mDLR.intialize_ggm_DLR_W();
	mult_mDLR.ScPT_solver(iter,scSE,scGF);
	std::string file_name_AC = "data/Summed_data/Summed_sigma_ac.txt";
	auto p = std::pair<int,int>(0,2);
	mult_mDLR.mDLR_list[0].write_data_momenta_AC_sigma_ij( file_name_AC,scGF,2000,p);
   
    AmiBase::graph_type basetype_2p;
	if (params.compute_2p == 1) {

		if (params.type_2p == 0 ){
		basetype_2p =AmiBase::Pi_phuu;  }

		if (params.type_2p == 1 ){
		basetype_2p =AmiBase::Pi_ppuu;  }

		else {
			std::cout << " Wrong 2p type type with :" <<  params.type_2p << std::endl;

		}


			if (rank==0){

			fs::remove_all("data_ph");
			std::string folder_name3 = "data_ph/Summed_data";
			std::string folder_name4 = "data_ph/Individual_data";	
			fs::create_directories(folder_name3);
			fs::create_directories(folder_name4);	
		}


		if (rank==0){ std::cout << " DOING STUFF for PH \n";}
		AmiGraph::gg_matrix_t ggm_ph;
		AmiGraph gp(basetype_2p, 0);
		gp.read_ggmp(params.graph_2p,ggm_ph, params.ord_max_2p);

		auto nodes =master_bf_ops.get_ifnodes();
		nda::array<dcomplex,1> bfreq_list(nodes.size());
		
		
		for (size_t i =0; i<nodes.size();i++){ 
		bfreq_list(i)=dcomplex(0,(2*nodes[i])*M_PI/beta);
		}
	

		if (rank ==0){
			std::cout << bfreq_list << std::endl;
		}
        Bz_container chi;
		
		ggm_mDLR  mult_ph_mDLR( params,basetype_2p, ggm_ph,master_if_ops);
		mult_ph_mDLR.ggm_Fk_2p_solver(bfreq_list);
		mult_ph_mDLR.intialize_ggm_DLR_W(scGF);
		//mult_ph_mDLR.transfer_ggm_DLR_W(scGF);
		mult_ph_mDLR.single_shot_2p_solver_qext(0,0,bfreq_list, chi);
		std::string file_name_AC = "data_ph/Summed_data/single_shot_AC_ph.txt";
		mult_ph_mDLR.mDLR_list[0].write_data_momenta_AC_chi( file_name_AC,master_bf_ops,chi,2000);








}
    






	MPI_Finalize();

	
}




// #include "dlr_ami.hpp"
// using namespace cppdlr;

// int main(int argc, char** argv){
// 	int size, rank;
// 	MPI_Init(&argc,&argv);
// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &size);
	

// 	params_param params;
// 	std::string loader_file = "../loader/params.txt";
// 	params_loader(loader_file,params);
// 	AmiBase::graph_type sigma=AmiBase::Sigma;
	
// 	////////////   Load the graphs here //////////////////
	
// 	AmiGraph::gg_matrix_t ggm;
// 	AmiGraph gs(sigma, 0);
// 	gs.read_ggmp(params.graph,ggm, params.ord_max);
// 	std::cout<<std::endl;
	
	
	
	
// 	AmiGraph::graph_t graph =ggm[params.ord_max][0].graph_vec[0];
	
	
	
// 	double NC = 161;
// 	double beta   = params.beta;
//     double eps    = params.eps;
// 	int kl = params.L;
// 	double Emax = params.Emax;
// 	double Uval = params.Uval;
// 	double lambda = beta*Emax;
// 	double mu = params.mu;
// 	double tp = params.tp;
// 	int iter = params.iter;
// 	AmiBase ami;
// 	//AmiBase::g_prod_t R0=construct_example2();
	
// 	double master_E = 5; double master_eps=1e-7;
// 	double master_lambda = params.beta*master_E;
//     auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
//     auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);
	
// 	mDLR multiple_DLR(params, sigma, graph,master_if_ops);
	

   
// 	//////// Computing the frequency kernel ////////////////
// 	if (rank ==0){
// 	std::cout << "--__--__--__--__--__--__--__--__--__--__--"<< std::endl;
// 	std::cout <<" Precomputing Computing the frequency kernel \n";
// 	}
// 	auto t0 = std::chrono::high_resolution_clock::now();
// 	auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
// 	nda::array<dcomplex,1> mfreq(nodes.size());
	
	
// 	for (size_t i =0; i<nodes.size();i++){ 
// 	mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
// 	}
	
	
// 	auto  frequency_kernel_list = nda::array<dcomplex, 2>(mfreq.size(), multiple_DLR.MPI_obj.count);
	
// 	for (int i=0; i <nodes.size();i++){
// 	auto val = mfreq(i);
// 	if (rank ==0){
// 	std::cout << val <<std::endl;
// 	}
// 	 frequency_kernel_list(i,nda::range::all) =multiple_DLR.evaluate_auxillary_energies(val); 
	
	
// 	}
// 	auto t1 = std::chrono::high_resolution_clock::now();
// 	if (rank==0){
// 		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
// 		std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";
// 	}
	
	
// 	multiple_DLR.populate_master_dlrW();
// 	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	
	
// 	double DENSITY_TOLERANCE = 1e-4;
// 	int    BISECT_STEPS    = 1000;
//     double alpha=0.3;
// 	int max_iters = iter; 
// 	auto density_lst = nda::array<double,1>(iter);
//     Bz_container SE_old;
// 	Bz_container SE;
// 	for (int iter_idx = 1; iter_idx < max_iters; ++iter_idx) {
// 		// 1) Build the kernel and compute self energy
// 		auto momenta_kernel = multiple_DLR.compute_momenta_kernel_bz();
// 		auto SE_new = multiple_DLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, frequency_kernel_list);
// 		// 2) Mix old SE with new SE 
		
// 		SE = ( iter_idx < 2) ? SE_new : multiple_DLR.SE_mixer(SE_old,SE_new, alpha);

// 		// 3) Compute density at current mu, adjust mu if needed
// 		double density = multiple_DLR.compute_density_from_SE(SE, mfreq, mu);
// 		double mu_new = (std::abs(params.target_n - density) < DENSITY_TOLERANCE)
// 						? mu
// 						: multiple_DLR.adjust_chemical_potential_bisc(params, SE, mfreq, BISECT_STEPS);

// 		double density_adj = multiple_DLR.compute_density_from_SE(SE, mfreq, mu_new);
// 		density_lst[iter_idx - 1] = density_adj;

// 		if (rank == 0) {
// 			std::cout << "Iteration " << iter_idx
// 					  << ": density = " << density_adj << std::endl;
// 		}

// 		// 4) Compute Green’s function and write out data
// 		Bz_container GF;
// 		if (params.DCA ==1){
// 		  GF = multiple_DLR.G_from_DLR_SE_M_DCA(SE, mfreq, mu_new,NC);
// 		}
// 		else {
// 		  GF = multiple_DLR.G_from_DLR_SE_M(SE, mfreq, mu_new);
// 		}
// 		std::string file_SE = std::format("../result/{}i_shot_SE.txt", iter_idx);
// 		std::string file_GF = std::format("../result/{}i_shot_GF.txt", iter_idx);

// 		multiple_DLR.write_data_momenta(file_SE, SE, mfreq);
// 		multiple_DLR.write_data_momenta(file_GF, GF, mfreq);

// 		// 5) Update master DLR weights for next iteration
// 		multiple_DLR.repopulate_master_dlrW_from_G(GF);
// 		multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
// 		// 6) Old SE back to current SE
// 		SE_old = SE;
// 	}

		
		
	

// 	MPI_Finalize();

	
// }



