function [filt] = COT_Type1_Filter(alpha_d,Df,Df0)

[nt,nd]=size(Df);
na=length(alpha_d);
Dfil_all = zeros(nt,nd);
unfilt = Df(:,1);
for jj=1:nt  
    Df1=diag(alpha_d)*Df0'+(eye(na)-diag(alpha_d))*Df(jj,1:nd)';
    Df0=Df1';
    Dfil_all(jj,:)=Df1';  
end
filt = Dfil_all; 

 
