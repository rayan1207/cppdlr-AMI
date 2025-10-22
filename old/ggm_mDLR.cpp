#include "dlr_ami.hpp"
using namespace cppdlr;


ggm_mDLR::ggm_mDLR(params_param _params, AmiGraph::gg_matrix_t _ggm): params(_params), ggm(_ggm) {
std::cout << "Initializing ggm_mDLR \n";
ggm_to_graphlist();
}
	
	
	
void ggm_mDLR::ggm_to_graphlist(){
	std::cout<< "Loading graphs \n";
	for (int i = params.ord_min; i < params.ord_max+1; ++i){
	 for (int j= 0; j< ggm[i].size(); ++j){
		 for (int k=0; k <ggm[i][j].graph_vec.size();++k){
			 AmiGraph::graph_t graph = ggm[i][j].graph_vec[k];
			 graph_list.push_back(graph);
			 std::cout<< "Sampling graph o"<<i<<"_g"<<j << "_n"<< k<<".graph"<< std::endl;
			}
		}
	}
	std::cout<< " Loading completed \n";


}	
	
	
	
	
