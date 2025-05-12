


#include "dlr_ami.hpp"
using namespace cppdlr;



int main(){
	double beta   = 5.0;
    double eps    = 1e-6;
    double e      = 0.5;
	size_t kl = 21;
	double Emax = 6;
	double lambda = beta*Emax;
	
	AmiBase ami;
	AmiBase::g_prod_t R0=construct_example2();
	
	mDLR multiple_DLR(beta,eps,Emax,kl,R0);
	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	
	// multiple_DLR.create_multiple_gstruct();
	// multiple_DLR.generate_cartesian_list();
	// multiple_DLR.generate_auxillary_energy_list();
	//print2d(multiple_DLR.cartesian_combo_list);
	
	
	
	std::vector<std::complex<double>>  mfreq_list;
	std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	auto t0 = std::chrono::high_resolution_clock::now();
	std::cout <<" Computing the frequency kernel \n";
	for (int i=0; i <16;i++){
	auto mfreq = std::complex<double>(0,(2*i+1)*M_PI/beta);
	std::cout << mfreq <<std::endl;
	mfreq_list.push_back(mfreq);
	auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(mfreq);
	frequency_kernel_list.push_back(frequency_kernel);
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	std::cout << " It took time: " <<duration.count() << " ms \n";
	
	
	
	
	
	
	
	
	
	
	
	//nda::array<double,1> energy={-4,0.1,-1};
	//auto weights =multiple_DLR.evaluate_auxillary_weights(energy);
	//std::cout << "result is : " << nda::dotc(frequency_kernel,weights);
	
	
	
	double kx_ext = M_PI; double ky_ext = M_PI;
	
	int internal_dof =  multiple_DLR.N -1;
	std::vector<double> kx_list (multiple_DLR.N);
	std::vector<double> ky_list (multiple_DLR.N);
	kx_list[multiple_DLR.N-1]=kx_ext;
	ky_list[multiple_DLR.N-1]=ky_ext;
	double dk = 2*M_PI/(multiple_DLR.kl-1);
	auto momenta_kernel = nda::zeros<dcomplex>(multiple_DLR.CN);
	
	
	std::cout<<"Computing the momenta kernel \n";
	for (int count=0;count <multiple_DLR.CN;count++){
		auto const &combo =  multiple_DLR.cartesian_combo_list[count];
		for (auto const  &k_combo: multiple_DLR.cartesian_k_combo_list){
         /////////////// put the inside into function
		   for (int i=0; i < internal_dof;i++){
			 kx_list[i] =   multiple_DLR.kvals[ k_combo[static_cast<int>(2*i)]];
			 ky_list[i] =   multiple_DLR.kvals[ k_combo[static_cast<int>(2*i+1)]];
			}
		    auto val = dcomplex(1,0);
			for (int i =0; i< 3;i++){
				auto const &gldr =multiple_DLR.multiple_dlr_structs[i];
				auto const &square = multiple_DLR.multiple_dlr_structs[i].dlrW_in_square;
				auto const &alpha = gldr.ginfo.alpha_;
				double qx = 0;
				double qy = 0;
				for (int j = 0; j < alpha.size(); j++) {
					qx += static_cast<double>(alpha[j]) * kx_list[j];
					qy += static_cast<double>(alpha[j]) * ky_list[j];
					}
				qx = std::fmod(std::abs(qx), 2*M_PI);
				qy = std::fmod(std::abs(qy), 2*M_PI);
				
				int indice1 = static_cast<int> (qx/dk);
				int indice2 = static_cast<int> (qy/dk);
			    // std::cout<< indice1 << indice2;
				auto const &weight_list  = square[indice1][indice2];
				
				val = -1*val* weight_list(combo[i]);
	
			}
			momenta_kernel(count) += val;	
        ////////////////////////////////////////////////////////////////////	
		}
		// std::cout << momenta_kernel_component;
		// momenta_kernel(count) = momenta_kernel_component;
		
	}
	
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::seconds>(t2 - t0);
	std::cout << "Computing the momenta kernel took: " <<duration1.count() << " s \n";
	
	//std::cout << "result is : " << nda::dotc(frequency_kernel_list[0],momenta_kernel)/(-1*std::pow((double) kl,4));

	for (int i =0; i< 16; i++){
		auto result = nda::dotc(frequency_kernel_list[i],momenta_kernel)/(-1*std::pow((double) kl,4));
		std::cout << i+1<<".  mfreq: " << mfreq_list[i] <<" Sigma: " <<  result <<std::endl;
	}


}
	
	
	
	
	
	
	
	
	
	// nda::array<dcomplex ,1> weights;
	
	
	// std::vector<nda::array<dcomplex,1>> G_dlr_list;
	// std::vector<nda::array<dcomplex,1>> G_dlr_w_list;
	
     
	// std::cout<<"1"<< multiple_DLR.cartesian_combo_list.size();
	// for (int i =0; i< multiple_DLR.N;i++){
		// std::cout << " Computing the greens function at energy \n ";
		// auto dlr_R0 =  multiple_DLR.multiple_dlr_structs[i];
		// auto gdlr_R0 = generate_nda_Gdlr_from_energy(dlr_R0, energy[i]);
		// auto weights = dlr_R0.if_ops.vals2coefs(beta, gdlr_R0);
		// std::cout << weights.size()<<std::endl;
		// G_dlr_list.push_back(gdlr_R0);
		// G_dlr_w_list.push_back(weights);
		
		// std::cout<< "For G_dlr with e=  " <<  energy[i] <<"\n";
		// std::cout<< gdlr_R0 <<std::endl;
		// std::cout<< "weights are" <<std::endl;
		// std::cout << weights << std::endl;
	// }
	
	
	// nda::array<dcomplex,1> full_weights (multiple_DLR.CN);
	// for (int i=0; i < multiple_DLR.CN;i++){
		// auto combo = multiple_DLR.cartesian_combo_list[i];
	    // auto result = dcomplex(1,0);
		// for (int j =0; j< multiple_DLR.N ; j++){
			// result = -1*result * G_dlr_w_list[j](combo[j]);
		// }
		// full_weights(i) = result;	
	// }
	
	// auto result = nda::dotc(full_weights,frequency_kernel);
	// std::cout << result;
	
	
	
	
	
	
	
		

// int main(){
    // double beta   = 5.0;
    // double eps    = 1e-10;
    // double e      = 0.5;
	// double Emax = 10.0;
	// double lambda = beta*Emax;
	
	
	
	
	// AmiBase ami;
	// AmiBase::g_prod_t R0=construct_example2();
	
	// dlr_obj dlr_R0; 
	// create_dlr_obj(beta,eps,Emax,dlr_R0,R0[0]);
	
	// std::cout << " Lets print the the imaginary frequecies \n\n";
	// cprint1d(dlr_R0.im_freqs);
	// std::cout << " Lets print the pole locations \n\n";
	// print1d(dlr_R0.pole_locs);
	// nda::array<dcomplex,1> gdlr=  generate_nda_Gdlr_from_energy( dlr_R0,e);
	
	// std::cout << " Constructing Matsubara greens function G_dlr at given energy(e) \n\n";
	// std::cout << gdlr <<std::endl;
	
	
	// std::cout << "Obtain the weights from  G_Dlr \n\n";
	
	// auto weights = dlr_R0.if_ops.vals2coefs(beta, gdlr);
	// std::cout << weights << "\n\n\n";
	
	
	// nda::array<dcomplex,1> gdlr_rec = recover_G_from_poles_n_weights( dlr_R0,weights, dlr_R0.im_freqs); 
	
	// for  (int i =0; i < weights.size(); i++){
	// std::cout << "matusbara freq:" << dlr_R0.im_freqs[i] <<"original: " <<  gdlr[i]<< "Recovered :" << gdlr_rec[i] << "err:" 
		// << gdlr[i]-gdlr_rec[i] <<std::endl;
	// }
	
	// print2d(dlr_R0.evec);
			
// }
// int main(){
    // double beta   = 5.0;
    // double eps    = 1e-10;
    // double e      = 0.5;
	// double Emax = 10.0;
	// double lambda = beta*Emax;
	

    // // 1) Build DLR real‑frequency nodes
    // auto dlr_rf = build_dlr_rf(lambda, eps);
    // imfreq_ops ifops(lambda, dlr_rf, Fermion);

     

	// auto if_nodes = ifops.get_ifnodes();  
	// std::cout << if_nodes;
    // int N        = if_nodes.size();
	
	// std::cout << " Now obtain evaluate the G on the DLR frequencies " <<"\n\n\n";
	
	
	// nda::array<dcomplex,1> G_dlr(if_nodes.size());
	
	// for (int i=0; i < if_nodes.size();i++){
		// auto dlr_n = if_nodes(i);
		// double wn = (2*dlr_n +1)*M_PI/beta;
		// G_dlr(i) = 1/ (dcomplex(0,wn)-e);
	// }
	// std::cout <<"\n";
	
	// std::cout << G_dlr << " \n\n";
	
	
	// std::cout << "Extrating the weight from the G_dlr \n\n";
	// auto weights = ifops.vals2coefs(beta, G_dlr);

	// std::cout << weights <<std::endl;
	// std::cout << "Grabbing the poles from G_dlr \n\n";
	// auto poles = ifops.get_rfnodes()/beta;
	// std::cout << poles <<std::endl;
	
	// auto G_dlr_rec = nda::zeros<dcomplex>( if_nodes.size() );
	
	// for (int i=0; i < if_nodes.size(); i++){
		// auto dlr_n = if_nodes(i);
		// auto wn = dcomplex(0,(2*dlr_n +1)*M_PI/beta);
		
		// for (int j=0; j < poles.size();j++){
			
		// G_dlr_rec(i) = G_dlr_rec(i)+( weights(j)/(wn - poles(j)));
		// }
		
		
		// std::cout << "matusbara freq:" << wn <<"original: " <<  G_dlr(i)<< "Recovered :" << G_dlr_rec(i) << "err:" 
		// << G_dlr(i)-G_dlr_rec(i) <<std::endl;
	// }
	
	// std::cout<<"-----Constructing AMI SPR format -----"<<std::endl;
	
	


// // class instance
	// AmiBase ami;

	// // Problem setup (see ami_example.cpp)
	// AmiBase::g_prod_t R0=construct_example2(); // Sets initial integrand 
	// AmiBase::ami_vars avars=construct_ext_example2(); // Sets 'external' parameter values 

		// //timing info
		// auto t1=std::chrono::high_resolution_clock::now();

	// // Storage objects for S,P,R 
	// AmiBase::S_t S_array;
	// AmiBase::P_t P_array;
	// AmiBase::R_t R_array;

	// // Integration/Evaluation parameters
	// double E_REG=0; // Numerical regulator for small energies.  If inf/nan results try E_REG=1e-8 
	// int N_INT=2;  // Number of Matsubara sums to perform
	// AmiBase::ami_parms test_amiparms(N_INT, E_REG);

	// //Construction Stage
	// ami.construct(test_amiparms, R0, R_array, P_array, S_array);
	// std::complex<double> calc_result=ami.evaluate(test_amiparms,R_array, P_array, S_array,  avars);
	// std::cout<<"Result was "<< calc_result<<std::endl;

// }
