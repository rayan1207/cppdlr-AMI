


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
	AmiGraph::gg_matrix_t ggm_first;
	gs.read_ggmp(params.graph,ggm_first, 2);

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
	


	
	
	double master_E = 5; double master_eps=1e-5;
	double master_lambda = params.beta*master_E;
    auto dlr_rf = build_dlr_rf(master_lambda,master_eps );
    auto master_if_ops = imfreq_ops(master_lambda, dlr_rf, Fermion);
	auto master_bf_ops = imfreq_ops(master_lambda, dlr_rf, Boson);
///////////// guesss //////////////////


//  if (rank==0){ std::cout << "QUICK GF2 Calculations for DCA guess \n";}
//    	Bz_container guess_scGF;
//    	Bz_container guess_scSE;

// 	ggm_mDLR  mult_mDLR_first( params,sigma, ggm_first,master_if_ops);
// 	mult_mDLR_first.ggm_Fk_solver();
// 	mult_mDLR_first.initialize_ggm_DLR_W();
	
// 	mult_mDLR_first.ScPT_solver(iter,guess_scSE,guess_scGF,false);

//  if (rank==0){ std::cout << "\n\n First G obtained succcessfully  \n\n";}


////////// main 1p cal /////////

if (rank ==0) { std::cout << "Beginning higher IPT calculation \n\n";}
	Bz_container scGF;
   	Bz_container scSE;
	ggm_mDLR  mult_mDLR( params,sigma, ggm,master_if_ops);
	mult_mDLR.ggm_Fk_solver();

	if (params.type_1p ==0){
		
		if (rank ==0) { std::cout << "Using Full K resolved IPT \n\n";}
		mult_mDLR.initialize_ggm_DLR_W();
		mult_mDLR.ScPT_solver(iter,scSE,scGF,true);}
	if (params.type_1p ==1){
    
		if (rank ==0) { std::cout << "Using DMFT  IPT \n\n";}
		// auto guess_scGF_local =  average_Bz_container();
		//std::cout << guess_scGF_local;
		mult_mDLR.initialize_local_ggm_DLR_W();
		mult_mDLR.ScPT_DMFT_solver(iter,scSE,scGF,true);

		}

	if (params.type_1p ==2){
		if (rank ==0) { std::cout << "Using SEET  IPT \n\n";}
		mult_mDLR.initialize_ggm_DLR_W();
		// auto guess_scGF_local =  average_Bz_container(guess_scGF);
		mult_mDLR.initialize_local_ggm_DLR_W();
		mult_mDLR.ScPT_LSEET_solver(iter,scSE,scGF,true);}


if (rank==0){ std::cout << " Done\n";}








	std::string file_name_AC = "data/Summed_data/Summed_sigma_ac.txt";
	auto p = std::pair<int,int>(0,0);
	mult_mDLR.mDLR_list[0].write_data_momenta_AC_sigma_ij( file_name_AC,scGF,6000,p);
   
    AmiBase::graph_type basetype_2p;
	if (params.compute_2p == 1) {

		if (params.type_2p == 0 ){
	    std::cout << "computing particle hole channel in 2p sector" << std::endl;
		basetype_2p =AmiBase::Pi_phuu;  }

		if (params.type_2p == 1 ){
		basetype_2p =AmiBase::Pi_ppuu;  }

		// else {
		// 	std::cout << " Wrong 2p type type with :" <<  params.type_2p << std::endl;

		// }


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
		mult_ph_mDLR.initialize_ggm_DLR_W(scGF);
		
		//mult_ph_mDLR.transfer_ggm_DLR_W(scGF);
		mult_ph_mDLR.single_shot_2p_solver_qext(0,0,bfreq_list, chi);
		std::string file_name_AC = "data_ph/Summed_data/single_shot_AC_ph.txt";
		mult_ph_mDLR.mDLR_list[0].write_data_momenta_AC_chi( file_name_AC,master_bf_ops,chi,4000);








}
    






	MPI_Finalize();

	
}




// #include "dlr_ami.hpp"
// using namespace cppdlr;


// struct EnergyItem {
//   double w;   // |x|^2
//   int idx;    // index in original vector
// };

// // Returns kept indices (sorted by descending energy) and kept values
// void l2_energy_truncate(
//     const nda::array<nda::dcomplex,1>& x,
//     double delta,                       // keep 1 - delta energy
//     std::vector<int>& kept_idx,
//     std::vector<nda::dcomplex>& kept_val,
//     double& kept_energy_frac
// ) {
//   const int n = (int)x.size();

//   // Compute weights and total energy
//   std::vector<EnergyItem> items;
//   items.reserve(n);

//   double total = 0.0;
//   for (int i = 0; i < n; ++i) {
//     double w = std::norm(x(i)); // |x|^2, no sqrt
//     total += w;
//     items.push_back({w, i});
//   }

//   if (total == 0.0) {
//     kept_idx.clear();
//     kept_val.clear();
//     kept_energy_frac = 0.0;
//     return;
//   }

//   // Sort by descending energy
//   std::sort(items.begin(), items.end(),
//             [](const EnergyItem& a, const EnergyItem& b) {
//               return a.w > b.w;
//             });

//   // Select smallest prefix reaching target energy
//   const double target = (1.0 - delta) * total;

//   kept_idx.clear();
//   kept_val.clear();
//   kept_idx.reserve(n);
//   kept_val.reserve(n);

//   double acc = 0.0;
//   for (int t = 0; t < n; ++t) {
//     if (items[t].w == 0.0) break; // remaining are all zero
//     acc += items[t].w;
//     kept_idx.push_back(items[t].idx);
//     kept_val.push_back(x(items[t].idx));
//     if (acc >= target) break;
//   }

//   kept_energy_frac = acc / total;
// }

// int main(int argc, char** argv){
// 	int size, rank;
// 	MPI_Init(&argc,&argv);
// 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
// 	MPI_Comm_size(MPI_COMM_WORLD, &size);
	

// 	params_param params;
// 	std::string loader_file = "params.txt";
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
	
	
// 	// for (size_t i =0; i<nodes.size();i++){ 
// 	for (size_t i =0; i<nodes.size(); i++){ 
// 	// mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);
// 	mfreq[i]=dcomplex(0,M_PI/beta);
// 	}
	
	
// 	auto  frequency_kernel_list = nda::array<dcomplex, 2>(mfreq.size(), multiple_DLR.MPI_obj.count);
	
// 	// for (int i=0; i <nodes.size();i++){
// 	// auto val = mfreq(i);
// 	// if (rank ==0){
// 	// std::cout << val <<std::endl;
// 	// }
// 	//  frequency_kernel_list(i,nda::range::all) =multiple_DLR.evaluate_auxillary_energies(val); 
	
	
// 	// }
// 	auto t1 = std::chrono::high_resolution_clock::now();
// 	if (rank==0){
// 		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
// 		std::cout << " Consturctin of frequency kernel took: " <<duration.count() << " ms \n";
// 	}
//     auto something = multiple_DLR.generate_wedge_info();
// 	// std::map<std::vector<int>, int > counts;

// 	// int kN = multiple_DLR.kN;
// 	// for (int i =0;i < kN; i++) {
// 	// 	auto tmp = multiple_DLR.compute_momenta_ind_each_GF(0,M_PI, i );
// 	// 	counts[tmp]++;
// 	// }
// 	// // for (int i =0;i < kN; i++) {
// 	// // 	auto tmp = multiple_DLR.compute_momenta_ind_each_GF(M_PI/2,M_PI/4, i );
// 	// // 	// std::cout << "i= " ; print1d(tmp); std::cout << "\n";
// 	// // 	counts[tmp]++;
// 	// // }

//     // int total_kN = 0;
// 	// int unique_key =0;
// 	// for (const auto&[key,val]: counts){
// 	// 	std::cout << "Vector { ";
//     //     for (int x : key) std::cout << x << " ";
//     //     std::cout << "} appeared " << val << " times.\n";
// 	// 	total_kN += val;
// 	// 	unique_key ++;
// 	// }

// 	// std::cout << " The size of KN:" << kN <<" \n Recovered size " << total_kN << std::endl;
// 	// std::cout << " The size of unique keys: " << unique_key;
	

	

// 	MPI_Finalize();

	
// }




	
	
	// multiple_DLR.populate_master_dlrW();
	// multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	
	
	// double DENSITY_TOLERANCE = 1e-4;
	// int    BISECT_STEPS    = 1000;
    // double alpha=0.3;
	// int max_iters = iter; 
	// auto density_lst = nda::array<double,1>(iter);
    // Bz_container SE_old;
	// Bz_container SE;
	// for (int iter_idx = 1; iter_idx < max_iters; ++iter_idx) {
	// 	// 1) Build the kernel and compute self energy
	// 	auto momenta_kernel = multiple_DLR.compute_momenta_kernel_bz();
	// 	auto SE_new = multiple_DLR.MPI_vdot_freq_momenta_kernel_M(momenta_kernel, frequency_kernel_list);
	// 	// 2) Mix old SE with new SE 
		
	// 	SE = ( iter_idx < 2) ? SE_new : multiple_DLR.SE_mixer(SE_old,SE_new, alpha);

	// 	// 3) Compute density at current mu, adjust mu if needed
	// 	double density = multiple_DLR.compute_density_from_SE(SE, mfreq, mu);
	// 	double mu_new = (std::abs(params.target_n - density) < DENSITY_TOLERANCE)
	// 					? mu
	// 					: multiple_DLR.adjust_chemical_potential_bisc(params, SE, mfreq, BISECT_STEPS);

	// 	double density_adj = multiple_DLR.compute_density_from_SE(SE, mfreq, mu_new);
	// 	density_lst[iter_idx - 1] = density_adj;

	// 	if (rank == 0) {
	// 		std::cout << "Iteration " << iter_idx
	// 				  << ": density = " << density_adj << std::endl;
	// 	}

	// 	// 4) Compute Green’s function and write out data
	// 	Bz_container GF;
	// 	if (params.DCA ==1){
	// 	  GF = multiple_DLR.G_from_DLR_SE_M_DCA(SE, mfreq, mu_new,NC);
	// 	}
	// 	else {
	// 	  GF = multiple_DLR.G_from_DLR_SE_M(SE, mfreq, mu_new);
	// 	}
	// 	std::string file_SE = std::format("../result/{}i_shot_SE.txt", iter_idx);
	// 	std::string file_GF = std::format("../result/{}i_shot_GF.txt", iter_idx);


	// 	multiple_DLR.write_data_momenta(file_SE, SE, mfreq);
	// 	multiple_DLR.write_data_momenta(file_GF, GF, mfreq);

	// 	// 5) Update master DLR weights for next iteration
	// 	multiple_DLR.repopulate_master_dlrW_from_G(GF);
	// 	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
		
	// 	// 6) Old SE back to current SE
	// 	SE_old = SE;
	// }

		
		
	



