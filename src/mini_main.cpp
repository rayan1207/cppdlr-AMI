


#include "dlr_ami.hpp"
using namespace cppdlr;



int main(){
	double beta   = 5.0;
    double eps    = 1e-6;
    double e      = 0.5;
	size_t kl = 11;
	double Emax = 6;
	double lambda = beta*Emax;
	
	AmiBase ami;
	AmiBase::g_prod_t R0=construct_example2();
	
	mDLR multiple_DLR(beta,eps,Emax,kl,R0);
    
	
	auto nodes = multiple_DLR.master_if_ops.get_ifnodes();
	nda::array<dcomplex,1> mfreq(nodes.size());
	for (size_t i =0; i<nodes.size();i++){ mfreq[i]=dcomplex(0,(2*nodes[i]+1)*M_PI/beta);}
	
	
	

	std::vector< nda::array<dcomplex,1>> frequency_kernel_list;
	auto t0 = std::chrono::high_resolution_clock::now();
	std::cout <<" Computing the frequency kernel \n";
	
	
	
	for (int i=0; i <nodes.size();i++){
	auto val = mfreq(i);
	std::cout << val <<std::endl;
	auto frequency_kernel=multiple_DLR.evaluate_auxillary_energies(val); 
	frequency_kernel_list.push_back(frequency_kernel);
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	std::cout << " Computing frequency kernel took time: " <<duration.count() << " ms \n";
	



	std::cout << "Computing momneta kernel" <<std::endl;
	
	
	nda::array<dcomplex,1> momenta_kernel =multiple_DLR.compute_momenta_kernel_qext(0,M_PI);
	
	
	
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
	std::cout << "Computing the momenta kernel took: " <<duration1.count() << " ms \n";
	
	
	
	for (int i =0; i< nodes.size(); i++){
    auto result = nda::dotc(frequency_kernel_list[i],momenta_kernel)/(-1*std::pow((double) kl,4));
    std::cout << i+1<<".  mfreq: " << mfreq[i] <<" Sigma: " <<  result <<std::endl;	
}


  auto momenta_kernel_bz =multiple_DLR.compute_momenta_kernel_bz();
  auto result =  multiple_DLR.vdot_freq_momenta_kernel_M(momenta_kernel_bz,frequency_kernel_list);

	
}