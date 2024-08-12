%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of MATLAB Code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uk] = HMPC_main_bioreactor(data, ts)
persistent time_i Forecast_Signals MeasDist_filt_k_true Dm_past_vec Dist_Forecast_over_p filt_ref_val

disp(time_i)
if isempty(Forecast_Signals)
   out_load= load('Forecast_Signals.mat');
   Forecast_Signals=out_load.Forecast_Signals;
end


%% Settings
%Choose Predictive Model. (MoD-1, ARX- 0)
model_choice=1;
% Choose solver (0 - Quadprog, 1 - CPLEX)
degree=0;

% Disturbance anticipation (0 - OFF, 1 - ON)
anticipation = 0;

% Slack variable and weight (0 - OFF, 1 - ON)
slack_var = 0;
slack_weight= 1000000;

% Simulation Definition
p=50;                               % Prediction horizon
m=20;                               % Control horizon

% Three Degrees of Tuning Parameters
tau_clr=5*ts;
tau_cld=ts*[8;8];
tau_clu=20*ts;
alpha_d=exp(-ts./tau_cld);
alpha_r=exp(-ts./tau_clr);          % 0 is faster than 1 % Setpoint tuning parameter
fa=1-exp(-ts/tau_clu);              % 1 is faster than 0 % Unmeasured disturbance: tuning parameter

% Dimension Declaration
ny=length(alpha_r);
nu=1;                                % number of manipulated inputs
n_meas_dist=length(alpha_d);         % number of measured disturbances
nud=0;                     % Number of discrete inputs


% Filter selection (Type)
unmeas_filter_sel=1;
meas_filter_sel = 2;
n_w = 5;            % Number of filter lags.

% Number of logical variables
ndc=1;

% MLD Matrices (E Matrices)
am=0;         % For making E1, E2, E3, E4, E5 equal to zero or not.
amp_d = 0;    % For making Ed equal to zero or not

% Initializations
Sim_Time=20*60;
y_sp=8*ones(Sim_Time,1); % Define Setpoint
% y_sp(1:400)=7.2;
% y_sp(850:end)=7.5;
y_init=7.2;

% Output range.


% Weights of Process Variables
move_suppression_true=0.001;
out_weight=1;
u_weight=0;

weightu = u_weight;
weightdu =(move_suppression_true);
weighty = (out_weight)';
weightd = zeros(ndc,1);
weightz = zeros(ndc,1);
weight_slack = zeros(length(weighty),1);
weight_slack(length(weighty)) = slack_weight;

% Limits of Process Variables
CO2_u_true=15;
pH_max=10;
pH_min=6;

umin = 0;
umax= inf;%CO2_u_true;
delumin= -inf;
delumax = inf;
ymax = inf;%pH_max;
ymin = -inf;%pH_min;
delta_min = zeros(ndc,1);
delta_max = ones(ndc,1);
z_min = zeros(ndc,1);
z_max = 500*ones(ndc,1);
slack_min = zeros(length(weighty),1);
slack_max = inf*ones(length(weighty),1);
rng(100);

% Iteration counter.
if isempty(time_i)
    time_i = 1;
end

%% Data reading
% Update y array
y_curr = data(1);
Dm_curr = [data(3), data(4)];
u_last=data(2);

%% Disturbance dynamics

% Disturbance model (see Nandola and Rivera IEEE TCST 2011)
if unmeas_filter_sel==1
    Aw_um = zeros(ny);
else
    Aw_um = eye(ny);  % For Type-1 (Single Intergrating) filter, Aw = zeros(ny), For Type-2 (Double Integrating) filter, Aw = eye(ny)
end
Bw_um = eye(ny);
Cw_um = eye(ny);
Unmeas_DistModel = struct('Aw',Aw_um,'Bw',Bw_um, 'Cw', Cw_um);

%% E matrixes
E1=am*degree*ones(10,nu);
E2=am*degree*ones(10,ndc);
E3=am*degree*ones(10,ndc);
E4=am*degree*ones(10,ny);
E5=am*degree*ones(10,1);
Ed=amp_d*ones(size(E1,1),n_meas_dist);

mld_copy.E1=E1;
mld_copy.E2=E2;
mld_copy.E3=E3;
mld_copy.E4=E4;
mld_copy.E5=E5;
mld_copy.Ed=Ed;

%% Optimization formulation parameters as per the dimension of the problem
Weights = struct('wy', weighty, 'wd', weightd, 'wz', weightz, ...
    'wu', weightu, 'wdu', weightdu, 'wslack', weight_slack);
limits = struct('yu', ymax, 'yl', ymin, 'uu', umax, 'ul', umin, ...
    'deluu', delumax, 'delul', delumin, 'delta_lo', delta_min, ...
    'delta_up', delta_max, 'z_up', z_max, 'z_lo', z_min, ...
    'slack_min', slack_min, 'slack_max', slack_max);

%% Measrured Disturbance Signals Generation

setpoints=struct('u',0,'del',zeros(1,1),'z',zeros(1,1)); % Setpoint Generation

%% Disturbance Vector Generation
% Temp_Forecast=Forecast_Signals(:,1);
% Rad_Forecast=Forecast_Signals(:,2);
Dist_Sig_Tstart=0; %Start Time of Controller
Dist_Sig_Tstop=20;
% dist_time_idx=(Dist_Sig_Tstart*60+1:1:(Dist_Sig_Tstop*60+p)+1);
dist_time_vec=(Dist_Sig_Tstart*60*60:ts:(Dist_Sig_Tstop*60+p)*60)';
% Temp_sig_ms=Temp_Forecast(dist_time_idx)-Temp_Forecast(Dist_Sig_Tstart*60+1);
% Rad_sig_ms=Rad_Forecast(dist_time_idx)-Rad_Forecast(Dist_Sig_Tstart*60+1);
% Dm_signal_ant=[dist_time_vec Temp_sig_ms Rad_sig_ms];

%% Measured Disturbance Filter and Anticipation

if isempty(MeasDist_filt_k_true)
     MeasDist_filt_k_true = zeros(anticipation*(p-1)+1,n_meas_dist);
end
% if isempty(Dist_Forecast_over_p)
%     Dist_Forecast_over_p = zeros(p,n_meas_dist);
% end


Dm_init=[17 0];
if anticipation==0
    if time_i<=1
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dm_curr,Dm_init);
        elseif meas_filter_sel==2
             [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dm_curr,Dm_init,n_w,dist_time_vec(time_i,:),anticipation,p,n_meas_dist);
        end
    else
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dm_curr,MeasDist_filt_k_true);
        elseif meas_filter_sel==2
            [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dm_curr,MeasDist_filt_k_true,n_w,dist_time_vec(time_i,:),anticipation,p,n_meas_dist);
        end
    end
    D_forecast_filt=MeasDist_filt_k_true(1,:)';
    D_forecast_unfilt=Dm_curr';
elseif anticipation==1
    Dist_Forecast_over_p=Forecast_Signals(time_i:time_i+p-1,:);
    Dist_Forecast_over_p(1,:)=Dm_curr;
    if time_i<=1
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dist_Forecast_over_p,Dm_init);
        elseif meas_filter_sel==2
             [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dist_Forecast_over_p,Dm_init,n_w,dist_time_vec(time_i:time_i+p-1,:),anticipation,p,n_meas_dist);
        end
    else
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dist_Forecast_over_p,MeasDist_filt_k_true(1,:));
        elseif meas_filter_sel==2
            [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dist_Forecast_over_p,MeasDist_filt_k_true(1,:),n_w,dist_time_vec(time_i:time_i+p-1,:),anticipation,p,n_meas_dist);
        end
    end
    D_f = MeasDist_filt_k_true';
    D_forecast_filt= D_f(:)';
    
    D_unf = Dist_Forecast_over_p';
    D_forecast_unfilt= D_unf(:)';
end

if isempty(Dm_past_vec)
    Dm_past_vec = [Dm_init;Dm_init;Dm_init];
end

Dm_past_vec=circshift(Dm_past_vec,1);
Dm_past_vec(1,:)=MeasDist_filt_k_true(1,:);
Dm_filt_minus1=Dm_past_vec(2,:)';
Dm_filt_minus2=Dm_past_vec(3,:)';


%%
% Dist_Forecast_k_plus_1_to_p=Forecast_Signals(i:i+p-1,:);
% MeasDist_filt_k_plus_1_to_p=zeros(p,n_meas_dist);
% if meas_filter_sel==1
%     [MeasDist_filt_k_plus_1_to_p]=COT_Type1_Filter(alpha_d,Dist_Forecast_k_plus_1_to_p,Dm_curr);
% elseif meas_filter_sel==2
% %      [MeasDist_filt_k_plus_1_to_p]=Type_2_Filter_outer(alpha_d,Dist_Forecast_k_plus_1_to_p,Dm_curr,dist_time_vec(i+1:i+p,:))
%      for dist_iter=1:length(alpha_d)
%         unfilt_sig=[dist_time_vec(i+1:i+p,:) Dist_Forecast_k_plus_1_to_p(:,dist_iter)];
%         [filt_sig,~] = Type2_Filter_Owais(n_w,alpha_d(dist_iter),unfilt_sig,Dm_curr(dist_iter));
%         MeasDist_filt_k_plus_1_to_p(:,dist_iter)=filt_sig;
%      end
% end
% if anticipation==0
%     D_forecast_filt=MeasDist_filt_k_plus_1_to_p(1,:)';
% elseif anticipation==1
%     D_f = MeasDist_filt_k_plus_1_to_p';
%     D_forecast_filt= D_f(:);
% end


% figure;subplot(2,1,1);plot(dist_time_vec(i+1:i+p,:)/60,MeasDist_filt_k_plus_1_to_p(:,1));hold on;plot(dist_time_vec(i+1:i+p,:)/60,Dist_Forecast_k_plus_1_to_p(:,1));subplot(2,1,2);plot(dist_time_vec(i+1:i+p,:)/60,MeasDist_filt_k_plus_1_to_p(:,2));hold on;plot(dist_time_vec(i+1:i+p,:)/60,Dist_Forecast_k_plus_1_to_p(:,2))
% figure;subplot(2,1,1);plot(dist_time_vec(time_i+1:time_i+p,:)/60,MeasDist_filt_k_true(:,1));hold on;plot(dist_time_vec(time_i+1:time_i+p,:)/60,Dist_Forecast_over_p(:,1));subplot(2,1,2);plot(dist_time_vec(time_i+1:time_i+p,:)/60,MeasDist_filt_k_true(:,2));hold on;plot(dist_time_vec(time_i+1:time_i+p,:)/60,Dist_Forecast_over_p(:,2))



%% Type-I Filter for Reference Trajectory
% filt_ref_val=0;
filt_ref=zeros(p,1);
if time_i==1
    filt_ref_val=y_init;
end
% figure;
for p_iter=1:p
    filt_ref_val = diag(alpha_r)*filt_ref_val+(eye(ny)-diag(alpha_r))*y_sp(time_i);
    filt_ref((p_iter-1)*ny+1:p_iter*ny,1)=filt_ref_val;
%     plot(filt_ref);drawnow;
end
filt_ref_val=filt_ref(1);


%%
uk=control_gen(y_curr,u_last,Dm_curr, filt_ref,time_i, Unmeas_DistModel,degree,slack_var,anticipation, Weights,limits,setpoints,p,m,nud,fa,D_forecast_filt,D_forecast_unfilt,Dm_filt_minus1,Dm_filt_minus2,mld_copy,model_choice,ts);

% Iteration counter update.
time_i=time_i+1;

% Flowrates = [uk, 0];
% Controller_name = 'MoD MPC';

end