function sys_out = ssDel2ABCD(sys)
% discretizes a system with internal delay in a state space model
%
% INPUT:
%   -sys : discrete system, it could be either state space or transfer
%          function
% OUTPUT:
%   -sys_out : state space model with the internal delays discretized
%
% Author : Francesco Campregher
% Maintanance e-mail : campregher99@gmail.com

if nargin < 1
    error("Not enough input argument")
end
if sys.Ts == 0
    error("The input system must be a discrete time system")
end
ss_sys = ss(sys);
if isempty(ss_sys.InternalDelay) && (isempty(ss_sys.InputDelay) || sum(ss_sys.InputDelay) == 0) && (isempty(ss_sys.OutputDelay) || sum(ss_sys.OutputDelay) == 0)
    sys_out = ss_sys;
    return
end

if ~isempty(ss_sys.InternalDelay)
    [H,tau] = getDelayModel(ss_sys);
    
    tau=tau-1;
    
    A=[H.A zeros(size(H.A,1),sum(tau)+length(tau));
       zeros(sum(tau)+length(tau),size(H.A,2)+sum(tau)+length(tau))];
    B=[H.B(1:size(ss_sys.B,1),1:size(ss_sys.B,2));
       zeros(sum(tau)+length(tau),size(ss_sys.B,2))];
    C=[H.C(1:size(ss_sys.C,1),1:size(ss_sys.C,2)) zeros(size(ss_sys.C,1),sum(tau)+length(tau))];
    D=H.d(1:size(ss_sys.D,1),1:size(ss_sys.D,2));
    
    for i = 1:length(tau)
    % Building A matrix expand
        if i == 1
            A(size(H.A,1)+1,1:size(H.A,2))=H.C(size(ss_sys.C,1)+i,:);
            for j = 1:length(tau)
                if j ==1
                    A(size(H.A,1)+1,size(H.A,2)+1:size(H.A,2)+tau(j)+j)=[zeros(1,tau(j)) H.D(size(ss_sys.D,1)+i,size(ss_sys.D,2)+j)];
                else
                    A(size(H.A,1)+1,size(H.A,2)+sum(tau(1:j-1))+j:size(H.A,2)+sum(tau(1:j))+j)=[zeros(1,tau(j)) H.D(size(ss_sys.D,1)+i,size(ss_sys.D,2)+j)];
                end
            end
            A(size(H.A,1)+i+1:size(H.A,1)+i+tau(i),size(H.A,2)+1:size(H.A,2)+tau(i))=eye(tau(i));
        else
            A(size(H.A,1)+sum(tau(1:i-1))+i,1:size(H.A,2))=H.C(size(ss_sys.C,1)+i,:);
            for j = 1:length(tau)
                if j ==1
                    A(size(H.A,1)+sum(tau(1:i-1))+i,size(H.A,2)+1:size(H.A,2)+tau(j)+j)=[zeros(1,tau(j)) H.D(size(ss_sys.D,1)+i,size(ss_sys.D,2)+j)];
                else
                    A(size(H.A,1)+sum(tau(1:i-1))+i,size(H.A,2)+sum(tau(1:j-1))+j:size(H.A,2)+sum(tau(1:j))+j)=[zeros(1,tau(j)) H.D(size(ss_sys.D,1)+i,size(ss_sys.D,2)+j)];
                end
            end
            A(size(H.A,1)+sum(tau(1:i-1))+i+1:size(H.A,1)+sum(tau(1:i-1))+i+tau(i),size(H.A,2)+sum(tau(1:i-1))+i:size(H.A,2)+sum(tau(1:i))+i-1)=eye(tau(i));
        end    
        
    % Building B matrix expand
        if i == 1
            B(size(ss_sys.B,1)+i,:)=H.D(size(ss_sys.D,1)+i,1:size(ss_sys.D,2));
        else
            B(size(ss_sys.B,1)+i+sum(tau(1:i-1)),:)=H.D(size(ss_sys.D,1)+i,1:size(ss_sys.D,2));
        end
        
    % Building C matrix expand
        C(:,size(ss_sys.C,2)+sum(tau(1:i))+i)=H.D(1:size(ss_sys.C,1),size(ss_sys.C,1)+i);
    end
    ss_sys_=ss(A,B,C,D,sys.Ts);
    ss_sys_.InputDelay=ss_sys.InputDelay;
    ss_sys_.OutputDelay=ss_sys.OutputDelay;
    ss_sys=ss_sys_;
end
if sum(ss_sys.InputDelay) ~= 0
    tau=ss_sys.InputDelay;
    

    A=[ss_sys.A zeros(size(ss_sys.A,1),sum(tau));
       zeros(sum(tau),size(ss_sys.A,2)+sum(tau))];
    B=zeros(sum(tau)+size(ss_sys.A,1),size(ss_sys.B,2));
    C=[ss_sys.C zeros(size(ss_sys.C,1),sum(tau))];
    D=ss_sys.D;
    for ind = 1:length(tau)
        if tau(ind) > 0
            if ind == 1
                A(size(ss_sys.A,1)+2:size(ss_sys.A,1)+tau(ind),size(ss_sys.A,1)+1:size(ss_sys.A,1)+tau(ind)-1)=eye(tau(ind)-1);
                B(size(ss_sys.A,1)+1,ind)=1;
            else
                A(size(ss_sys.A,1)+sum(tau(1:ind-1))+2:size(ss_sys.A,1)+sum(tau(1:ind)),size(ss_sys.A,1)+1:size(ss_sys.A,1)+tau(ind)-1)=eye(tau(ind)-1);
                B(size(ss_sys.A,1)+sum(tau(1:ind))+ind,ind)=1;
            end
            A(1:size(ss_sys.A,1),size(ss_sys.A,1)+tau(ind))=ss_sys.B(:,ind);
        else
            B(1:size(ss_sys.A),ind)=ss_sys.B(:,ind);
        end
    end
    ss_sys_=ss(A,B,C,D,sys.Ts);
    ss_sys_.OutputDelay=ss_sys.OutputDelay;
    ss_sys=ss_sys_;
end
if sum(ss_sys.OutputDelay) ~= 0
    tau=ss_sys.OutputDelay;
    A=[ss_sys.A zeros(size(ss_sys.A,1),sum(tau));
       zeros(sum(tau),size(ss_sys.A,2)+sum(tau))];
    B=[ss_sys.B;
       zeros(sum(tau),size(ss_sys.B,2))];
    C=zeros(size(ss_sys.C,1),size(ss_sys.C,2)+sum(tau));
    D=ss_sys.D;
    for ind = 1:length(tau)
        if tau(ind) > 0
            if ind == 1
                A(size(ss_sys.A,1)+1,1:size(ss_sys.A,2))=ss_sys.C(ind,:);
                A(size(ss_sys.A,1)+2:size(ss_sys.A,1)+tau(ind),size(ss_sys.A,2)+1:size(ss_sys.A,2)+tau(ind)-1)=eye(tau(ind)-1);
            else
                A(size(ss_sys.A,1)+sum(tau(1:ind-1))+1,1:size(ss_sys.A,2))=ss_sys.C(ind,:);
                A(size(ss_sys.A,1)+sum(tau(1:ind-1))+2:size(ss_sys.A,1)+sum(tau(1:ind)),size(ss_sys.A,2)+sum(tau(1:ind-1))+1:size(ss_sys.A,2)+sum(tau(1:ind))-1)=eye(tau(ind)-1);
            end
            C(ind,size(ss_sys.A,2)+sum(tau(1:ind)))=1;
        else
            C(ind,1:size(ss_sys.A,2))=ss_sys.C(ind,:);
        end
    end
end
sys_out = ss(A,B,C,D,sys.Ts);
end

