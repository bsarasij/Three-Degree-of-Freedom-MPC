function  [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dm_all,Df0_all,n_w,tvec,anticipation,p,n_meas_dist)

MeasDist_filt_k_true = zeros(anticipation*(p-1)+1,n_meas_dist);

for dist_iter=1:length(alpha_d)
    unfilt_sig=[tvec Dm_all(:,dist_iter)];
    [filt_sig,~] = Type2_Filter_Owais(n_w,alpha_d(dist_iter),unfilt_sig,Df0_all(dist_iter));
    MeasDist_filt_k_true(:,dist_iter)=filt_sig;
end
