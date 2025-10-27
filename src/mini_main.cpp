


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
	std::string loader_file = "../loader/params.txt";
	params_loader(loader_file,params);
	
	
	////////////   Load the graphs here //////////////////
	
	AmiGraph::gg_matrix_t ggm;
	g.read_ggmp(params.graph,ggm, params.ord_max);
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
	AmiBase ami;

	
	
	double master_E = 5; double master_eps=1e-7;
	double master_lambda = params.beta*master_E;
    auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
    auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);



	ggm_mDLR  mult_mDLR( params, ggm,master_if_ops);
	mult_mDLR.ggm_Fk_solver();
	mult_mDLR.intialize_ggm_DLR_W();
	mult_mDLR.ScPT_solver(iter);
	// //// temporary lets print SE
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
	
	
// 	////////////   Load the graphs here //////////////////
	
// 	AmiGraph::gg_matrix_t ggm;
// 	g.read_ggmp(params.graph,ggm, params.ord_max);
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
	
// 	mDLR multiple_DLR(beta,Uval,eps,Emax,kl,tp,graph,master_if_ops);
	

   
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
	
	
// 	std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	
// 	for (int i=0; i <nodes.size();i++){
// 	auto val = mfreq(i);
// 	if (rank ==0){
// 	std::cout << val <<std::endl;
// 	}
// 	auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(val); 
// 	frequency_kernel_list.push_back(frequency_kernel);
	
// 	}
// 	auto t1 = std::chrono::high_resolution_clock::now();
// 	if (rank==0){
// 		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
// 		std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";
// 	}
	
	
// 	multiple_DLR.populate_master_dlrW_from_G0(mu);
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
	
	
// 	////////////   Load the graphs here //////////////////
	
// 	AmiGraph::gg_matrix_t ggm;
// 	g.read_ggmp(params.graph,ggm, params.ord_max);
// 	std::cout<<std::endl;
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

	
	
// 	double master_E = 5; double master_eps=1e-7;
// 	double master_lambda = params.beta*master_E;
//     auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
//     auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);

// 	////// testing /////////
// 	AmiGraph::graph_t graph1 =ggm[2][0].graph_vec[0];
// 	AmiGraph::graph_t graph2 =ggm[3][0].graph_vec[0];
// 	AmiGraph::graph_t graph3 =ggm[3][0].graph_vec[1];

// 	mDLR multiple_DLR1(beta,Uval,eps,Emax,kl,tp,graph1,master_if_ops);
// 	mDLR multiple_DLR2(beta,Uval,eps,Emax,kl,tp,graph2,master_if_ops);
// 	mDLR multiple_DLR3(beta,Uval,eps,Emax,kl,tp,graph3,master_if_ops);

// 	auto nodes = master_if_ops.get_ifnodes();
// 	nda::array<dcomplex,1> mfreq(nodes.size());

// 	for (size_t i =0; i<nodes.size();i++){ mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);}
   
// 	std::vector< nda::array<dcomplex,1>> frequency_kernel_list1;
// 	std::vector< nda::array<dcomplex,1>> frequency_kernel_list2;
// 	std::vector< nda::array<dcomplex,1>> frequency_kernel_list3;

// 	for (int i=0; i <nodes.size();i++){ 
// 	auto val = mfreq(i);
// 	if (rank ==0){std::cout << val <<std::endl;}
// 	auto frequency_kernel1=multiple_DLR1.evaluate_auxillary_energies(val); 
// 	auto frequency_kernel2=multiple_DLR2.evaluate_auxillary_energies(val); 
// 	auto frequency_kernel3=multiple_DLR3.evaluate_auxillary_energies(val); 


// 	frequency_kernel_list1.push_back(frequency_kernel1);
// 	frequency_kernel_list2.push_back(frequency_kernel2);
// 	frequency_kernel_list3.push_back(frequency_kernel3);
// 	}

// 	multiple_DLR1.populate_master_dlrW_from_G0(mu);
// 	multiple_DLR1.transfer_master_DLR_weights_to_dlrR0_elements();

// 	multiple_DLR2.populate_master_dlrW_from_G0(mu);
// 	multiple_DLR2.transfer_master_DLR_weights_to_dlrR0_elements();

// 	multiple_DLR3.populate_master_dlrW_from_G0(mu);
// 	multiple_DLR3.transfer_master_DLR_weights_to_dlrR0_elements();


// 	auto momenta_kernel1 = multiple_DLR1.compute_momenta_kernel_bz();
//     auto SE_1 = multiple_DLR1.MPI_vdot_freq_momenta_kernel_M(momenta_kernel1, frequency_kernel_list1);
// 	multiple_DLR1.write_data_momenta("graph1_res.txt", SE_1, mfreq);


//     auto momenta_kernel2 = multiple_DLR2.compute_momenta_kernel_bz();
//     auto SE_2= multiple_DLR2.MPI_vdot_freq_momenta_kernel_M(momenta_kernel2, frequency_kernel_list2);
// 	multiple_DLR2.write_data_momenta("graph2_res.txt", SE_2, mfreq);


// 	auto momenta_kernel3 = multiple_DLR3.compute_momenta_kernel_bz();
//     auto SE_3= multiple_DLR3.MPI_vdot_freq_momenta_kernel_M(momenta_kernel3, frequency_kernel_list3);
// 	multiple_DLR3.write_data_momenta("graph3_res.txt", SE_3, mfreq);








	
	
	
	

// 	// ggm_mDLR  mult_mDLR( params, ggm,master_if_ops);
// 	// mult_mDLR.ggm_Fk_solver();
// 	// mult_mDLR.intialize_ggm_DLR_W();
// 	// std::vector<Bz_container> SE_list(mult_mDLR.graph_size);
// 	// for (int i =0; i< mult_mDLR.graph_size; i++){
	
// 	// 	SE_list[i] = mult_mDLR.generate_SE( mult_mDLR.mDLR_list[i], mult_mDLR.Fk_ggm[i]);	
// 	// 	auto name = mult_mDLR.graphlist_names[i];
// 	// 	std::string file_SE = std::format("../result/SE_{}.txt", name);
// 	// 	//std::string file_SE = std::format("../result/SE_graph_{}.txt", i+1);
// 	// 	mult_mDLR.mDLR_list[i].write_data_momenta(file_SE, SE_list[i], mult_mDLR.master_mfreq);	
// 	// }
	


// 	// std::cout << "It's complete" << std::endl;
// 	// //// temporary lets print SE

// 	MPI_Finalize();

	
// }


    // auto poles = master_if_ops.get_rfnodes()/beta;
    // auto nodes = master_if_ops.get_ifnodes();
	// auto energy = nda::dcomplex(3,0);
    // size_t n = nodes.size();
	// nda::array<dcomplex,1> glist (n);
	// nda::array<dcomplex,1> glist_r (n);
	// nda::array<dcomplex,1> mfreq (n);

	// for (int i =0;i < n;i++){
	// 	auto f= dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
	// 	mfreq(i) = f;
	// 	glist(i) = 1/(f- energy);
	// }

	// //std::cout << "Mfreq values are :" << mfreq <<std::endl << "glist values are :" << glist << std::endl;

	// auto weights= master_if_ops.vals2coefs(beta,glist);

	// for (int i =0;i < n;i++){
	// 	for (int j =0; j< poles.size();j++){
	//      glist_r(i) += weights[j]/(mfreq(i) -poles(j));
	// 	}		
	// }

    // for (int i =0;i < n;i++){
	// std::cout << "mfreq :" << mfreq(i) << "gdata :" << glist(i) << "recovered gdata:" << glist_r(i) <<std::endl;
	// }


