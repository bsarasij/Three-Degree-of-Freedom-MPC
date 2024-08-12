function [ref_set,y_init1] = ref_mpc(y_init,y_set,ref_para,p)
ref_set0=y_init;
ny=length(y_set);
disp(y_set)
% ref_set=ones(p*length(alpha_r),1)
for i=1:p
    ref_set0 = diag(ref_para)*ref_set0+(eye(ny)-diag(ref_para))*y_set;
%     disp("ref_set0")
%     disp(ref_set0)
%     disp("ref_set")
    ref_set((i-1)*ny+1:i*ny,1)=ref_set0;
end
y_init1=ref_set(1:ny,:);

