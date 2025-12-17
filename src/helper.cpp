#include "dlr_ami.hpp"
using namespace cppdlr;


std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r\f\v");
    size_t last = str.find_last_not_of(" \t\n\r\f\v");
    if (first == std::string::npos || last == std::string::npos)
        return "";
    return str.substr(first, (last - first + 1));
}

// Function to parse the parameter name and value from a line
void parseLine(const std::string& line, std::string& paramName, std::string& paramValue) {
    size_t equalPos = line.find('=');
    if (equalPos != std::string::npos) {
        paramName = trim(line.substr(0, equalPos));
        paramValue = trim(line.substr(equalPos + 1));
    }
}

// Function to fill the struct from the provided text file
void params_loader(const std::string& filename, params_param& params) {


    // Open the file and read its contents
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::string paramName, paramValue;
        parseLine(line, paramName, paramValue);
        if (paramName == "Uval")
            params.Uval = std::stod(paramValue);
        else if (paramName == "Emax")
            params.Emax = std::stod(paramValue);
		else if (paramName == "mu")
        params.mu = std::stod(paramValue);
        else if (paramName == "eps")
            params.eps = std::stod(paramValue);
        else if (paramName == "beta")
            params.beta = std::stod(paramValue);
        else if (paramName == "L")
            params.L = std::stoi(paramValue);
		else if (paramName == "iter")
            params.iter = std::stoi(paramValue);
		else if (paramName == "tp")
		params.tp = std::stod(paramValue);
		else if (paramName == "target_n")
			params.target_n = std::stod(paramValue);
		else if (paramName == "graph")
             params.graph = paramValue;	
		else if (paramName == "graph_2p")
             params.graph_2p = paramValue;	
		else if (paramName == "ord_max")
             params.ord_max = std::stoi(paramValue);
        else if (paramName == "ord_min")
             params.ord_min = std::stoi(paramValue);
		else if (paramName == "ord_max_2p")
             params.ord_max_2p = std::stoi(paramValue);
        else if (paramName == "ord_min_2p")
             params.ord_min_2p = std::stoi(paramValue);
		else if (paramName == "mu_L")
             params.mu_L = std::stod(paramValue);
		else if (paramName == "mu_R")
             params.mu_R = std::stod(paramValue);
		else if (paramName == "DCA")
			 params.DCA = std::stoi(paramValue);
		else if (paramName == "patch_N")
			 params.patch_N= std::stoi(paramValue);
		else if (paramName == "compute_2p")
			 params.compute_2p= std::stoi(paramValue);
		else if (paramName == "type_2p")
			 params.type_2p= std::stoi(paramValue);
		else if (paramName == "gfunc")
			 params.gfunc= std::stoi(paramValue);
    }
	
    inputFile.close();
	std::cout << "Loaded params file with values: " << " Beta="  << params.beta << ", Emax="<<  params.Emax 
	<< ", eps = " << params.eps << ", L= " << params.L<< ", mu= " << params.mu << ", tp= "<< params.tp << " and SCS iteration = " << params.iter << "min and max " <<params.ord_min << " " << params.ord_max
	<< "patch_N = " << params.patch_N <<  ". ord_max_2p= " << params.ord_max_2p <<  ", ord_min_2p" << params.ord_min_2p << std::endl<<  ", gloc_sigma =" << params.graph << ", gloc_2p = " << params.graph_2p << std::endl
	<< "compute_2p = " << params.compute_2p << " type_2p = "  << params.type_2p << " gfunc = " << params.gfunc <<  std::endl;
}




void mDLR::write_data_momenta(const std::string& filename,
                            Bz_container& data,
                            nda::array<dcomplex,1>& mfreq  )
{
    // Extract the NDA 1D array for the (i,j) momenta
    size_t k1_l = data.size();
	size_t k2_l = data[0].size();

	if (MPI_obj.rank == 0){
    int size = data[0][0].size();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }



    // Optional: write header
    //file << "# wn    qx     qy     Re       Im\n";
    for (int i =0; i<k1_l;++i){
		for (int j =0; j<k2_l;++j){		
			for (int f = 0; f < size; ++f) {
				double wn_imag = mfreq[f].imag();
				double qx = kvals[i];
				double qy = kvals[j];
				double re = data[i][j](f).real();
				double im = data[i][j](f).imag();

			
				// Write to file
				file << wn_imag << "  "
					 << qx << "  "
					 << qy << "  "
					 << re << "  "
					 << im << "\n";
					}
		}
	}

    file.close();
    std::cout << "Data written to " << filename << std::endl;}
}


	
	




// AmiBase::g_prod_t mDLR::create_R0_from_graph() {
// 	AmiBase::g_prod_t R0;
// 	AmiGraph::edge_vector_t fermionic_edge;
// 	g.find_internal_fermionic_edges(graph,fermionic_edge);
// 	int n = fermionic_edge.size();

// 	std::cout << " Constructing an R0 element with: \n"; 
// 	for (int i =0; i < n;i++){
// 		AmiBase::alpha_t alpha = graph[fermionic_edge[i]].g_struct_.alpha_;
// 		AmiBase::epsilon_t epsilon = graph[fermionic_edge[i]].g_struct_.eps_;
// 		// std::cout << " i = " << i << std::endl;
// 		// print1d(alpha);
// 		// print1d(epsilon);
		
		
// 		AmiBase::g_struct g(epsilon,alpha);
// 		R0.push_back(g);	
// 	}
		
// 	return R0;	
	
// }

std::vector<std::complex<double>> convertToComplex(const std::vector<double> vec) {
    std::vector<std::complex<double>> cplx_vec;
    cplx_vec.reserve(vec.size()); 

    for (auto elem : vec) {
        cplx_vec.emplace_back(elem, 0.0);
    }

    return cplx_vec;
}





double mDLR::hubbard_dispersion(double kx, double ky,double mu){
  double e =-2.0*(std::cos(kx) + std::cos(ky))-4*tp*(std::cos(kx)*std::cos(ky)) -mu;
  return e;
}

double fermi_distribution(double energy, double beta ){
	double arg = (energy*beta);
	return 1.0/(1+std::exp(arg));
}

void mDLR::fill_dlro_pole_info(){
	    total_num =1;
		for (auto dlr_R0 : multiple_dlr_structs){
		total_num =total_num*dlr_R0.pole_num;
		num_pole_each_dlr.push_back(dlr_R0.pole_num);	
	}	
}

void mDLR::fill_dlro_momenta_info(){
	for (int i = 0; i< 2*DOF;i++){
		num_k_each_dlr.push_back(kl);	
	}
	kcombo_element.resize(2*DOF);		
}

nda::array<dcomplex,1> mDLR::LocalG_from_DLR_SE_M(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu) {
	   dcomplex prefactor = dcomplex(kl*kl,0);
	   int mfreq_size = mfreq.size();
	   auto  local_G = nda::zeros<dcomplex>(mfreq_size);
	   for (int m = 0; m < mfreq_size;m++){
		   for (int i = 0; i < kl;i++){
			   for( int j =0; j< kl;j++){
				   double e = hubbard_dispersion(kvals[i],kvals[j],mu);
				   local_G(m) +=  1/(mfreq(m) - e - SE[i][j](m));    
			   }   
		   }    
	   }
	   return local_G/prefactor;	
}


double mDLR::compute_density_from_SE(Bz_container &SE,nda::array<dcomplex,1> &mfreq,double mu){
	auto local_G = LocalG_from_DLR_SE_M(SE,mfreq,mu);
	auto local_weights = master_if_ops.vals2coefs(beta,local_G);
    auto density = nda::dot(local_weights,fd_master_poles).real();
	return 2.0*density;
}
// double mDLR::compute_density_from_SE(Bz_container &SE,
//                                      nda::array<dcomplex,1> &mfreq, // master grid
//                                      double mu)
// {
//     // --- 0. Build local G on the *master* Matsubara grid ---
//     // mfreq is assumed to correspond to master_if_ops.get_ifnodes().
//     auto local_G_master = LocalG_from_DLR_SE_M(SE, mfreq, mu);

//     // --- 1. Get DLR weights for local G in the *master* DLR basis ---
//     // master_if_ops and master_poles are members of mDLR
//     auto w_master = master_if_ops.vals2coefs(beta, local_G_master);
//     // master_poles: nda::array<double,1> of size master_pole_num
//     // corresponds to master_if_ops.get_rfnodes()/beta

//     // --- 2. Define a separate "density DLR" basis (Emax, eps for density only) ---
//     double d_Emax = 10.0;      // tune this
//     double d_eps  = 1e-10;     // tune this

//     double lambda_d = beta * d_Emax;
//     auto dlr_rf_d   = build_dlr_rf(lambda_d, d_eps);
//     cppdlr::imfreq_ops d_if_ops(lambda_d, dlr_rf_d, Fermion);

//     // --- 3. Build the density-DLR Matsubara grid ---
//     auto d_nodes = d_if_ops.get_ifnodes();              // integer nodes
//     int  Nd      = d_nodes.size();

//     nda::array<dcomplex,1> d_mfreq(Nd);
//     for (int i = 0; i < Nd; ++i) {
//         int    n  = d_nodes(i);
//         double wn = (2 * n + 1) * M_PI / beta;         // fermionic Matsubara
//         d_mfreq(i) = dcomplex(0.0, wn);
//     }

//     // --- 4. Reconstruct local G on the *density* grid using master poles+weights ---
//     // G_dens(iω') = sum_j w_master(j) / ( iω' - ε_master(j) )
//     nda::array<dcomplex,1> local_G_dens = nda::zeros<dcomplex>(Nd);
//     for (int i = 0; i < Nd; ++i) {
//         dcomplex iw = d_mfreq(i);
//         dcomplex val(0.0, 0.0);
//         for (int j = 0; j < master_pole_num; ++j) {
//             val += w_master(j) / ( iw - dcomplex(master_poles(j), 0.0) );
//         }
//         local_G_dens(i) = val;
//     }

//     // --- 5. Get DLR weights in the *density* basis ---
//     auto w_dens = d_if_ops.vals2coefs(beta, local_G_dens);

//     // --- 6. Build Fermi factors on the density poles ---
//     auto density_poles = d_if_ops.get_rfnodes() / beta;   // ε_ℓ for density DLR
//     int  Np            = density_poles.size();
//     nda::array<dcomplex,1> fd_density(Np);

//     for (int i = 0; i < Np; ++i) {
//         double e = density_poles(i);
//         fd_density(i) = dcomplex(fermi_distribution(e, beta), 0.0);
//     }

//     // --- 7. Contract weights with Fermi factors: n = 2 * Σ_ℓ w_ℓ f(ε_ℓ) ---
//     double density = nda::dot(w_dens, fd_density).real();
//     return 2.0 * density;
// }

nda::array<dcomplex,1> mDLR::fd_on_master_poles(){
	auto fd_master_poles = nda::array<dcomplex,1>(master_pole_num);
	for (int i =0;i< master_pole_num;i++){
		fd_master_poles(i) = dcomplex(fermi_distribution(master_poles(i), beta),0);
	}
	return fd_master_poles;
}

double mDLR::non_interacting_density(nda::array<dcomplex,1> &mfreq,double mu){
	  dcomplex prefactor = dcomplex(kl*kl,0);
	  int mfreq_size = mfreq.size();
	  auto  local_G = nda::zeros<dcomplex>(mfreq_size);
	   for (int m = 0; m < mfreq_size;m++){
		   for (int i = 0; i < kl;i++){
			   for( int j =0; j< kl;j++){
				   double e = hubbard_dispersion(kvals[i],kvals[j],mu);
				   local_G(m) +=  1/(mfreq(m) - e );    
			   }   
		   }    
	   }
	local_G /= prefactor;
	auto local_weights = master_if_ops.vals2coefs(beta,local_G);
    auto density = nda::dot(local_weights,fd_master_poles).real();
	return 2.0*density;
		
}





double mDLR::adjust_chemical_potential_bisc(params_param &params, Bz_container &SE,nda::array<dcomplex,1> &mfreq, int max){
	
	double a = params.mu_L;
	double b = params.mu_R;
	
	double fa = params.target_n - compute_density_from_SE(SE,mfreq,a);
	double fb = params.target_n - compute_density_from_SE(SE,mfreq,b);
	
	
	if (fa*fb > 0){
		std::cerr<< " Pick the correct lower (a) and upper (b) bound chemical potential \n";
	}
	double c;
	for (int i =0; i < max;i++){
		c = (a+b)/2.0;
		double fc = params.target_n - compute_density_from_SE(SE,mfreq,c);
		if ( std::abs(fa-fb)< 1e-7){
			std::cout << "  Correct chemical potential is found on " << i << "-th iteration "<< std::endl;
			std::cout << " Adjusted mu = " << c << std::endl;
			
			return c;
		}
		
		if (fa*fc < 0){
			b= c;
			fb = fc;	
		}
		
		else {
		 a= c;
		 fa = fc;	
		}	
	}
	
	std::cerr << "Couldnt find chemical potential \n";
	return c;
		
}

Bz_container mDLR::SE_mixer(Bz_container SE_old, Bz_container SE_new, double alpha){
	
	int fk_size = SE_old[0][0].size();
	Bz_container result(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(fk_size)));
	
	for (int i =0;i< kl;i++){
		for (int j =0;j< kl;j++){
			for (int f=0;f< fk_size;f++){	
				result[i][j](f) =  dcomplex(alpha,0)*SE_old[i][j](f) + dcomplex(1-alpha,0)*SE_new[i][j](f);
			
			}
	
		}
		
	}
	return result;
}

Bz_container sum_containers(const std::vector<Bz_container>& all)
{
    if (all.empty()) return {};


    Bz_container result = all.front();

    for (size_t n = 1; n < all.size(); ++n) {
        const auto& cont = all[n];
        for (size_t i = 0; i < cont.size(); ++i) {
            for (size_t j = 0; j < cont[i].size(); ++j) {
                result[i][j] += cont[i][j];  
        	}
   	 	}
	}
	return result;
}


bool SE_converged_abs(const Bz_container &Sigma_old,
                      const Bz_container &Sigma_new,
                      double abs_tol)
{
    double max_abs_diff = 0.0;

    for (std::size_t ix = 0; ix < Sigma_new.size(); ++ix) {
        for (std::size_t iy = 0; iy < Sigma_new[ix].size(); ++iy) {
            const auto &se_old = Sigma_old[ix][iy];
            const auto &se_new = Sigma_new[ix][iy];

            for (std::size_t n = 0; n < se_new.size(); ++n) {
                double diff = std::abs(se_new(n) - se_old(n));
                if (diff > max_abs_diff)
                    max_abs_diff = diff;
            }
        }
    }

    return max_abs_diff < abs_tol;
}
nda::array<dcomplex,1> mDLR::prepare_AC_chi_data(  cppdlr::imfreq_ops &bosonic_if_ops, nda::array<dcomplex,1> &chi,int size){
	double beta = params.beta;
	auto recovered_chi = nda::zeros<dcomplex> (size);
	auto dlr_poles = bosonic_if_ops.get_rfnodes()/beta;
	auto weights= bosonic_if_ops.vals2coefs(beta,chi);
	for (int i=0; i < size;i++){
		auto iw = dcomplex(0,2*i*M_PI/beta);
		for (int j=0; j < weights.size(); j++){
			double  input = beta*dlr_poles[j]/2;
			recovered_chi(i) +=weights(j)*std::tanh(input)/(iw - dlr_poles[j]);		
		}			
	}
	return recovered_chi;
}


nda::array<dcomplex,1> mDLR::prepare_AC_sigma_data(   nda::array<dcomplex,1> &G,int size){
	double beta = params.beta;
	auto recovered_G = nda::zeros<dcomplex> (size);
	auto dlr_poles =master_if_ops.get_rfnodes()/beta;
	auto weights= master_if_ops.vals2coefs(beta,G);
	for (int i=0; i < size;i++){
		auto iw = dcomplex(0,(2*i+1)*M_PI/beta);
		for (int j=0; j < weights.size(); j++){
			recovered_G(i) +=weights(j)/(iw - dlr_poles[j]);		
		}			
	}
	return recovered_G;
}


void mDLR::write_data_momenta_AC_chi( const std::string& filename,cppdlr::imfreq_ops &bosonic_if_ops,
                            Bz_container& bosonic_data,
                             int size)
{
    // Extract the NDA 1D array for the (i,j) momenta
    size_t k1_l = bosonic_data.size();
	size_t k2_l = bosonic_data[0].size();

	if (MPI_obj.rank == 0){


    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }


    file << std::fixed << std::setprecision(17);
    // Optional: write header
    //file << "# wn    qx     qy     Re       Im\n";
    for (int i =0; i<k1_l;++i){
		for (int j =0; j<k2_l;++j){	
			auto ph = 	bosonic_data[i][j];
			auto result =prepare_AC_chi_data(bosonic_if_ops,ph, size);
			for (int f = 0; f < size; ++f) {
				double wn_imag = f;
				double qx = kvals[i];
				double qy = kvals[j];
				double re = result(f).real();
				double im = result(f).imag();

			
				// Write to file
				file << wn_imag << "  "
					 << qx << "  "
					 << qy << "  "
					 << re << "  "
					 << im << "\n";
					}
		}
		
	}

    file.close();
    std::cout << "Data written to " << filename << std::endl;}
}


void mDLR::write_data_momenta_AC_sigma( const std::string& filename,
                            Bz_container& fermionic_data,
                             int size)
{
    // Extract the NDA 1D array for the (i,j) momenta
    size_t k1_l = fermionic_data.size();
	size_t k2_l = fermionic_data[0].size();

	if (MPI_obj.rank == 0){
 

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }


    file << std::fixed << std::setprecision(10);
    // Optional: write header
    //file << "# wn    qx     qy     Re       Im\n";
    for (int i =0; i<k1_l;++i){
		for (int j =0; j<k2_l;++j){	
			auto sigma = 	fermionic_data[i][j];
			auto result =prepare_AC_sigma_data(sigma, size);
			for (int f = 0; f < size; ++f) {
				double wn_imag = (2*f+1)*M_PI/beta;
				double qx = kvals[i];
				double qy = kvals[j];
				double re = result(f).real();
				double im = result(f).imag();

			
				// Write to file
				file << wn_imag << "  "
					 << qx << "  "
					 << qy << "  "
					 << re << "  "
					 << im << "\n";
					}
		}
		
	}

    file.close();
    std::cout << "Data written to " << filename << std::endl;}
}



void mDLR::denoise_FS_points(Bz_container &SE){
	int klx = SE.size();
	int kly = SE[0].size();

	for (int i=0;i < klx;i++){
		for (int j=0;j<kly;j++){
			 auto e =  std::abs(hubbard_dispersion(kvals[i],kvals[j],params.mu));
			if (e < 1e-8){
				 auto &Sigma_k = SE[i][j];
				for (auto &z : Sigma_k ){ z = dcomplex(0, z.imag());
				}
			}
		}
	}
}

void mDLR::symmetrize_fermionic_DLR_array(nda::array<dcomplex, 1> &data) {
    double beta = params.beta;

    auto weights      = master_if_ops.vals2coefs(beta, data);
    auto master_nodes = master_if_ops.get_ifnodes();

    auto abs_nodes = nda::abs(master_nodes);
    auto abs_max   = nda::max_element(abs_nodes);

    int start = -static_cast<int>(abs_max);
    int end   =  static_cast<int>(abs_max);

    auto expanded_nodes = nda::arange(start, end + 1);
    int count = expanded_nodes.size();

    auto expanded_data = nda::zeros<dcomplex>(count);
    auto dlr_poles     = master_if_ops.get_rfnodes() / beta;

    for (int i = 0; i < count; ++i) {
        int n = expanded_nodes(i);
        dcomplex iw(0.0, (2.0 * n + 1.0) * M_PI / beta);

        dcomplex f_pos = 0.0;
        dcomplex f_neg = 0.0;

        for (int j = 0; j < weights.size(); ++j) {
            f_pos += weights(j) / ( iw - dlr_poles(j));
            f_neg += weights(j) / (-iw - dlr_poles(j));
        }

        expanded_data(i) = 0.5 * (f_pos + std::conj(f_neg));
    }

    for (int j = 0; j < master_nodes.size(); ++j) {
        int n = master_nodes(j);
        int i = n - start;
        data(j) = expanded_data(i);
    }
}





void mDLR::symmetrize_fermionic_DLR_Bz(Bz_container &SE){
	int klx = SE.size();
	int kly = SE[0].size();

	for (int i=0;i < klx;i++){
		for (int j=0;j<kly;j++){
			 auto &Sigma_k = SE[i][j];
			 symmetrize_fermionic_DLR_array(Sigma_k);		 
		}
	}
}



void mDLR::write_data_momenta_AC_sigma_ij( const std::string& filename,
                            Bz_container& fermionic_data,
                             int size,std::pair<int, int> ind)
{
    // Extract the NDA 1D array for the (i,j) momenta
    size_t k1_l = fermionic_data.size();
	size_t k2_l = fermionic_data[0].size();

	if (MPI_obj.rank == 0){
 

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

 
    file << std::fixed << std::setprecision(17);
    // Optional: write header
    //file << "# wn    qx     qy     Re       Im\n";
            int i = ind.first; int j= ind.second;
			auto sigma = 	fermionic_data[i][j];
			auto result =prepare_AC_sigma_data(sigma, size);
			for (int f = 0; f < size; ++f) {
				double wn_imag = (2*f+1)*M_PI/beta;
				double qx = kvals[i];
				double qy = kvals[j];
				double re = result(f).real();
				double im = result(f).imag();

			
				// Write to file
				file << wn_imag << "  "
					 << qx << "  "
					 << qy << "  "
					 << re << "  "
					 << im << "\n";
		
		
	}

    file.close();
    std::cout << "Data written to " << filename << std::endl;}
}


// nda::array<dcomplex,1> mDLR::evaluate_auxillary_weights(nda::array<double,1> &energy) {
// 	nda::array<dcomplex ,1> weights;
// 	std::vector<nda::array<dcomplex,1>> G_dlr_list;
// 	std::vector<nda::array<dcomplex,1>> G_dlr_w_list;
// 	for (int i =0; i< N;i++){
// 		auto dlr_R0 =  multiple_dlr_structs[i];
		
// 		auto gdlr_R0 = generate_nda_Gdlr_from_energy(dlr_R0.if_ops, energy[i]);
// 		auto weights = dlr_R0.if_ops.vals2coefs(beta, gdlr_R0);
// 		G_dlr_list.push_back(gdlr_R0);
// 		G_dlr_w_list.push_back(weights);
// 	}
	
	
// 	nda::array<dcomplex,1> full_weights (CN);
// 	for (int i=0; i < CN;i++){
// 		auto const& combo  = cartesian_combo_list[i];
// 	    auto result = dcomplex(1,0);
// 		for (int j =0; j< N ; j++){
// 			result = result * G_dlr_w_list[j](combo[j]);
// 		}
// 		full_weights(i) = result;	
// 	}
	
// 	return full_weights;
// }