
#include "dlr_ami.hpp"
using namespace cppdlr;
int main() {


	double beta   = 10;
    double eps    = 1e-7;
	double Emax = 5;
	auto e = dcomplex(2,0);


	double lambda = beta*Emax;
    auto dlr_rf = build_dlr_rf(lambda, eps);

	auto if_ops = imfreq_ops(lambda, dlr_rf, Fermion);
	auto nodes =if_ops.get_ifnodes();
	
	auto all_poles = if_ops.get_rfnodes()/beta;
	
	std::cout << " Poles are : \n" << all_poles << std::endl;
	std::cout << " Nodes are : \n" << nodes << std::endl; 
	
	size_t N = nodes.size();
    auto G = nda::zeros<dcomplex>( N);

	


	for (size_t i = 0; i < N; ++i) {
        dcomplex mf(0.0, (2*nodes(i) + 1) * M_PI / beta);
		G(i) = 1.0 / ( mf - e );
    }
	
	std::cout << " Greens functions are : \n\n";
	std::cout << G <<"\n\n";
	
	auto weights= if_ops.vals2coefs(beta,G);
	
	auto recovered_G = nda::zeros<dcomplex> (N);

	
	for (int i=0; i < N;i++){
		dcomplex iw(0.0, (2*nodes(i) + 1) * M_PI / beta);
		for (int j=0; j < weights.size(); j++){
			recovered_G(i) +=-weights(j)/(iw + all_poles[j]);		
		}			
	}
	std::cout << " Recovered Greens functions are : \n\n";
	std::cout << recovered_G;
	
	
	return 0;
	
	
	
}