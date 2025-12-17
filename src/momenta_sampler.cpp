#include "dlr_ami.hpp"





inline void mDLR::generate_ith_momenta_cartesian_combo(int i, std::vector<int>& result  ) {
    // int DOF = (baseType==AmiBase::Sigma) ? 2*ord : 2*(ord+1);
	int q = i;


    for (int j = 0; j < 2*DOF; ++j) {
        int base     = num_k_each_dlr[j];
        result[j]    = q % base;  
        q           /= base;      
    }
	
}






inline nda::dcomplex mDLR::compute_momenta_one_kCN_kernel(double kx_ext,double ky_ext,std::vector<int> &combo,const std::vector<int> &kcombo)
{
    nda::dcomplex val(1.0, 0.0);
    
    for (int i = 0; i < N; ++i) {

        const auto& info      = multiple_dlr_structs[i].ginfo;
        const int* alpha   = info.alpha_.data();
        int   asz      = info.alpha_.size();
        double qx = alpha[asz - 1] * kx_ext;
        double qy = alpha[asz - 1] * ky_ext;
        for (int j = 0; j < DOF; ++j) {
            double a = static_cast<double>(alpha[j]);
            qx += a * kvals[ kcombo[2*j    ] ];
            qy += a * kvals[ kcombo[2*j + 1] ];
        }


        qx -= std::floor(qx* inv_two_pi) * two_pi;
        qy -= std::floor(qy* inv_two_pi) * two_pi;

        int idx1 = static_cast<int>(std::round(qx * inv_dk));
        int idx2 = static_cast<int>(std::round(qy * inv_dk));
		idx1 = (idx1 % kl + kl) % kl;
		idx2 = (idx2 % kl + kl) % kl;
		
		
        const auto& wlist =
          multiple_dlr_structs[i].dlrW_in_square[idx1][idx2];
        val *=   wlist[ combo[i] ];
    }

    return val;
}


nda::dcomplex mDLR::gfunc(double kx, double ky){
	if (params.gfunc==0){
		return  nda::dcomplex(1,0);
	}
	if (params.gfunc==1){
		return  nda::dcomplex(std::sin(kx),0);
	}
	if (params.gfunc==2){
		return  nda::dcomplex(std::cos(kx)-std::cos(ky),0);
	}
	if (params.gfunc==3){
		return  nda::dcomplex(std::cos(kx)+std::cos(ky),0);
	}
	if (params.gfunc==4){
		return  nda::dcomplex(std::sin(kx)*std::sin(ky),0);
	}
	else {
		return  nda::dcomplex(1,0);
	}
}

nda::dcomplex  mDLR::apply_ggkp(double kx_ext, double ky_ext, const std::vector<int> &kcombo){
	nda::dcomplex val(1.0, 0.0);
	int asz = kkp[0].size();
	for (int i =0; i < 2; i++ ) {
		auto alpha = kkp[i];
	    double qx = alpha[asz - 1] * kx_ext;
        double qy = alpha[asz - 1] * ky_ext;
        for (int j = 0; j < DOF; ++j) {
            double a = static_cast<double>(alpha[j]);
            qx += a * kvals[ kcombo[2*j    ] ];
            qy += a * kvals[ kcombo[2*j + 1] ];
        }

		val *=  gfunc(qx,qy);
	}
	return val;

}




 nda::array<nda::dcomplex,1> mDLR::compute_momenta_kernel_qext(double kx_ext,double ky_ext)
{
    nda::array<nda::dcomplex,1> kernel = nda::zeros<nda::dcomplex>(CN);

    for (int c = 0; c < CN; c++) {
        // const int* combo_ptr = cartesian_combo_list[c].data();
		auto combo = generate_single_CN(c);
        auto sum = nda::dcomplex(0,0);
        for (int k = 0; k < kN; k++) {   
            generate_ith_momenta_cartesian_combo(k,kcombo_element);
            sum += compute_momenta_one_kCN_kernel(
                       kx_ext, ky_ext,
                       combo,
                       kcombo_element);
        }
		
        kernel(c) = sum;
    }

    return kernel;
}

 nda::array<nda::dcomplex,1> mDLR::compute_momenta_kernel_2p_qext(double kx_ext,double ky_ext)
{
    nda::array<nda::dcomplex,1> kernel = nda::zeros<nda::dcomplex>(CN);

    for (int c = 0; c < CN; c++) {
        // const int* combo_ptr = cartesian_combo_list[c].data();
		auto combo = generate_single_CN(c);
        auto sum = nda::dcomplex(0,0);
        for (int k = 0; k < kN; k++) {   
            generate_ith_momenta_cartesian_combo(k,kcombo_element);
            sum += compute_momenta_one_kCN_kernel(
                       kx_ext, ky_ext,
                       combo,
                       kcombo_element)*apply_ggkp(kx_ext,ky_ext,kcombo_element);
        }
		
        kernel(c) = sum;
    }

    return kernel;
}

Bz_container mDLR::compute_momenta_kernel_bz(){
	// int n = (kl+1)/2;
	if (MPI_obj.rank == 0){std::cout<< "Computing momenta kernel \n" ;}

	int n = kl/2 + 1;



	auto reduced_kgrid = nda::array<double,1> (n);
	
	Bz_container M(n,std::vector<nda::array<dcomplex,1>>(n));
	for(size_t i=0; i<n; i++){reduced_kgrid[i] = dk* double(i);}
	int data = n*(n+1)/2;
	int count =1;
	auto t0 = std::chrono::high_resolution_clock::now();
	for (int i =0;i< n; i ++){
		for (int j=0;j<=i;j++){
			double qx = reduced_kgrid(i);
			double qy = reduced_kgrid(j);
			if (MPI_obj.rank == 0){
			std::cout<< "Computing -> " << count <<"/" << data << " data point ";}
			auto momenta_kernel = mDLR::compute_momenta_kernel_qext(qx,qy);
			if (MPI_obj.rank == 0){std::cout<< "(computed) \n" ;}
			M[i][j] = momenta_kernel;
			count++;
		}		
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	if (MPI_obj.rank==0){std::cout << " computation of momenta kernel took: " <<duration.count() << " ms \n";}
	triangle_to_square(M); 
	
	return data_to_full_bz(M);
	 
}


Bz_container mDLR::compute_momenta_kernel_2p_bz(){
	// int n = (kl+1)/2;
	if (MPI_obj.rank == 0){std::cout<< "Computing momenta kernel \n" ;}

	int n = kl/2 + 1;



	auto reduced_kgrid = nda::array<double,1> (n);
	
	Bz_container M(n,std::vector<nda::array<dcomplex,1>>(n));
	for(size_t i=0; i<n; i++){reduced_kgrid[i] = dk* double(i);}
	int data = n*(n+1)/2;
	int count =1;
	auto t0 = std::chrono::high_resolution_clock::now();
	for (int i =0;i< n; i ++){
		for (int j=0;j<=i;j++){
			double qx = reduced_kgrid(i);
			double qy = reduced_kgrid(j);
			if (MPI_obj.rank == 0){
			std::cout<< "Computing -> " << count <<"/" << data << " data point ";}
			auto momenta_kernel = mDLR::compute_momenta_kernel_2p_qext(qx,qy);
			if (MPI_obj.rank == 0){std::cout<< "(computed) \n" ;}
			M[i][j] = momenta_kernel;
			count++;
		}		
	}
	auto t1 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
	if (MPI_obj.rank==0){std::cout << " computation of momenta kernel took: " <<duration.count() << " ms \n";}
	triangle_to_square(M); 
	
	return data_to_full_bz(M);
	 
}





Bz_container mDLR::MPI_vdot_freq_momenta_kernel_M(Bz_container &mk, nda::array<dcomplex,2> &fk){
	// int pp = (baseType==AmiBase::Sigma) ? ord : ord+1;
	int Nf = fk.shape()[0];
	std::vector<dcomplex> local_result;
	std::vector<dcomplex> global_result;
	local_result.resize(kl*kl*Nf);
	global_result.resize(kl*kl*Nf);
	
	for (int i =0;i<kl;i++){
		for (int j=0;j<kl;j++){
			auto const &momenta_kernel = mk[i][j];
			for (int k=0; k<Nf; k++){
				local_result[(i*kl+j)*Nf+k] = prefactor*std::pow(Uval,ord)*nda::dot(fk(k,nda::range::all),momenta_kernel)/(std::pow((double) kl*kl,DOF));
			}			
		}	
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(local_result.data(), global_result.data(), kl*kl*Nf, MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
	
	Bz_container result(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(Nf)));
	
	
	
	//// unflatten the data///
	for (int i =0;i<kl;i++){
		for (int j=0;j<kl;j++){
			auto  &out = result[i][j];
			for (int k =0; k<Nf;k++){
				out(k) = global_result[(i*kl+j)*Nf+k];
			}			
		}	
	}
	
	
	
	return result;	
}

 Bz_container mDLR::MPI_vdot_freq_momenta_kernel(nda::array<dcomplex,1> &momenta_kernel , nda::array<dcomplex,2> &fk){
	// int pp = (baseType==AmiBase::Sigma) ? ord : ord+1;
	int Nf = fk.shape()[0];
	std::vector<dcomplex> local_result;
	std::vector<dcomplex> global_result;
	local_result.resize(Nf);
	global_result.resize(Nf);
	

	for (int k=0; k<Nf; k++){
				local_result[k] = prefactor*std::pow(Uval,ord)*nda::dot(fk(k,nda::range::all),momenta_kernel)/(std::pow((double) kl*kl,DOF));
			}			

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(local_result.data(), global_result.data(), Nf, MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
	
	Bz_container result(1,std::vector<nda::array<dcomplex,1>>(1, nda::array<dcomplex,1>(Nf)));
	
	
	result[0][0] = global_result;
	return result;	
}