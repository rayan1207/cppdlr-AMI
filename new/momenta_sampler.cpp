#include "dlr_ami.hpp"





inline void mDLR::generate_ith_momenta_cartesian_combo(int i, std::vector<int>& result  ) {
    int q = i;
    
    for (int j = 0; j < 2*ord; ++j) {
        int base     = num_k_each_dlr[j];
        result[j]    = q % base;  // one modulo
        q           /= base;      // one division
    }
}






inline nda::dcomplex mDLR::compute_momenta_one_kCN_kernel(double kx_ext,double ky_ext,const int* combo_ptr,const std::vector<int> &kcombo)
{
    nda::dcomplex val(1.0, 0.0);

    for (int i = 0; i < N; ++i) {

        const auto& info      = multiple_dlr_structs[i].ginfo;
        const int* alpha   = info.alpha_.data();
        int   asz      = info.alpha_.size();
        double qx = alpha[asz - 1] * kx_ext;
        double qy = alpha[asz - 1] * ky_ext;
        for (int j = 0; j < ord; ++j) {
            double a = static_cast<double>(alpha[j]);
            qx += a * kvals_ptr[ kcombo[2*j    ] ];
            qy += a * kvals_ptr[ kcombo[2*j + 1] ];
        }


        qx -= std::floor(qx* inv_two_pi) * two_pi;
        qy -= std::floor(qy* inv_two_pi) * two_pi;

        int idx1 = static_cast<int>(std::round(qx * inv_dk));
        int idx2 = static_cast<int>(std::round(qy * inv_dk));
		idx1 = (idx1 % kl + kl) % kl;
		idx2 = (idx2 % kl + kl) % kl;
		
		
        const auto& wlist =
          multiple_dlr_structs[i].dlrW_in_square[idx1][idx2];
        val *=   wlist[ combo_ptr[i] ];
    }

    return val;
} 




 nda::array<nda::dcomplex,1> mDLR::compute_momenta_kernel_qext(double kx_ext,double ky_ext)
{
    nda::array<nda::dcomplex,1> kernel = nda::zeros<nda::dcomplex>(CN);

    for (int c = 0; c < CN; c++) {
        const int* combo_ptr = cartesian_combo_list[c].data();
        auto sum = nda::dcomplex(0,0);
        for (int k = 0; k < kN; k++) {   
            generate_ith_momenta_cartesian_combo(k,kcombo_element);
            sum += compute_momenta_one_kCN_kernel(
                       kx_ext, ky_ext,
                       combo_ptr,
                       kcombo_element);
        }
		
        kernel(c) = sum;
    }

    return kernel;
}

Bz_container mDLR::compute_momenta_kernel_bz(){
	// int n = (kl+1)/2;
	int n = kl/2 + 1;
	auto reduced_kgrid = nda::array<double,1> (n);
	
	Bz_container M(n,std::vector<nda::array<dcomplex,1>>(n));
	for(size_t i=0; i<n; i++){reduced_kgrid[i] = dk* double(i);}
	int data = n*(n+1)/2;
	int count =1;
	for (int i =0;i< n; i ++){
		for (int j=0;j<=i;j++){
			double qx = reduced_kgrid(i);
			double qy = reduced_kgrid(j);
			if (MPI_obj.rank == 0){
			std::cout<< "Computing -> " << count <<"/" << data << " data point \n";}
			auto momenta_kernel = mDLR::compute_momenta_kernel_qext(qx,qy);
			M[i][j] = momenta_kernel;
			count++;
		}		
	}
	triangle_to_square(M);
	
	return data_to_full_bz(M);
	 
}	


Bz_container mDLR::MPI_vdot_freq_momenta_kernel_M(Bz_container mk, std::vector<nda::array<dcomplex,1>> fk){
	
	int Nf = fk.size();
	std::vector<dcomplex> local_result;
	std::vector<dcomplex> global_result;
	local_result.resize(kl*kl*Nf);
	global_result.resize(kl*kl*Nf);
	
	for (int i =0;i<kl;i++){
		for (int j=0;j<kl;j++){
			auto const &momenta_kernel = mk[i][j];
			for (int k=0; k<Nf; k++){
				local_result[(i*kl+j)*Nf+k] = prefactor*std::pow(Uval,ord)*nda::dotc(fk[k],momenta_kernel)/(std::pow((double) kl*kl,ord));
			}			
		}	
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(local_result.data(), global_result.data(), kl*kl*Nf, MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
	
	Bz_container result(kl,std::vector<nda::array<dcomplex,1>>(kl, nda::array<dcomplex,1>(fk.size())));
	
	
	
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