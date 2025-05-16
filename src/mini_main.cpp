


#include "dlr_ami.hpp"
using namespace cppdlr;



int main(){
	params_param params;
	std::string loader_file = "params.txt";
	params_loader(loader_file,params);
	
	double beta   = params.beta;
    double eps    = params.eps;
	int kl = params.L;
	double Emax = params.Emax;
	double Uval = params.Uval;
	double lambda = beta*Emax;
	
	AmiBase ami;
	AmiBase::g_prod_t R0=construct_example2();
	
	mDLR multiple_DLR(beta,Uval,eps,Emax,kl,R0);
    
	//////////////////// Computing the frequency kernel ////////////////
	 std::cout << "--__--__--__--__--__--__--__--__--__--__--"<< std::endl;
	std::cout <<" Precomputing Computing the frequency kernel \n";
	auto t0 = std::chrono::high_resolution_clock::now();
	auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	for (size_t i =0; i<nodes.size();i++){ mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);}
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
	multiple_DLR.populate_master_dlrW_from_G0();
	std::cout << "Done";
	std::cout << "Transferring DLR weights to gstructs \n";
	//multiple_DLR.reshape_dlrW_square_per_kgrid();
	multiple_DLR.transfer_master_DLR_weights_to_dlrR0_elements();
	std::cout <<"Done";
	
    int iter = params.iter;
    for (int i=1;i<iter;i++){
		std::cout << "--__--__--__--__--__--__--__--__--__--__-- "<<std::endl;
		std::cout << "computing for momneta kernel for the "<< i <<" shot" <<std::endl;
		auto momenta_kernel_bz =multiple_DLR.compute_momenta_kernel_bz();
		auto SE =  multiple_DLR.vdot_freq_momenta_kernel_M(momenta_kernel_bz,frequency_kernel_list);
		auto GF = multiple_DLR.G_from_DLR_SE_M(SE, mfreq);
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
		std::cout << "done";
       
	}

	
}