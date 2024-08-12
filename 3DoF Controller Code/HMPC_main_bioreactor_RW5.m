%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of MATLAB Code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Flowrates,Controller_name] = HMPC_main_bioreactor_RW5(data,ts)

persistent time_i MeasDist_filt_k_true Dm_past_vec ...
Dist_Forecast_over_p filt_ref_val mean_val

addpath(genpath('MoD_OldColors'));

%% Settings
if isempty(mean_val)
    out_mean=load('mean_08_6.mat');
    mean_val=out_mean.mean_08_6;
end


% Choose Predictive Model. (MoD - 1, ARX - 0)
model_choice = 1 ;

% Choose solver (0 - Quadprog, 1 - CPLEX)
degree=0;

% Disturbance anticipation (0 - OFF, 1 - ON)
anticipation = 0;

% Slack variable and weight (0 - OFF, 1 - ON)
slack_var = 0;
slack_weight= 1000000000000;

% Simulation Definition
p=150;                              % Prediction horizon
m=20;                               % Control horizon

% Three Degrees of Tuning Parameter
tau_clr=5*ts;
tau_cld=ts*[1;1];
tau_clu=1*ts;
alpha_d=exp(-ts./tau_cld);
alpha_r=exp(-ts./tau_clr);          % 0 is faster than 1 % Setpoint tuning parameter
fa=1-exp(-ts/tau_clu);              % 1 is faster than 0 % Unmeasured disturbance: tuning parameter

% Dimension Declaration
ny=length(alpha_r);
nu=1;                               % Number of manipulated inputs
n_meas_dist=length(alpha_d);        % Number of measured disturbances
nud=0;                              % Number of discrete inputs

% Filter selection (Type)
unmeas_filter_sel=2;
meas_filter_sel = 2;
n_w = 5;                            % Number of filter lags.

% Number of logical variables
ndc=1;

% MLD Matrices (E Matrices)
am=0;         % For making E1, E2, E3, E4, E5 equal to zero or not.
amp_d = 0;    % For making Ed equal to zero or not

% Initializations
Sim_Time=24*60*60/ts+1;
y_sp=8*ones(Sim_Time,1)-mean_val(1); % Define Setpoint
%y_sp(1:400)=7.2;
%y_sp(850:end)=7.5;

% Weights of Process Variables
move_suppression_true=0.5;
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
pH_max=14;
pH_min=0;

umin = 0-mean_val(2);
umax= inf-mean_val(2);
delumin= -inf;
delumax = inf;
ymax = inf;
ymin = -inf;
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

% Current time and minute.
CT = datetime('now');
CM = minute(CT) + 60*hour(CT);

% Start and ending time (if that condition is set)
ST = 7;
ET = 21;

%% Data reading
% Update y array
y_curr = double(data(3).Value)-mean_val(1);
Dm_curr = double([data(18).Value, data(15).Value])-mean_val(3:4);
u_last = double(data(25).Value)-mean_val(2);

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

E1=am*degree*ones(1,nu);
E2=am*degree*ones(1,ndc);
E3=am*degree*ones(1,ndc);
E4=am*degree*ones(1,ny);
E5=am*degree*ones(1,1);
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
% Dist_Sig_Tstart=0; %Start Time of Controller
% Dist_Sig_Tstop=24;
% dist_time_vec=(Dist_Sig_Tstart*60*60:ts:(Dist_Sig_Tstop*60+p)*60)';
dist_time_vec = 0;

%% Measured Disturbance Filter and Anticipation

if isempty(MeasDist_filt_k_true)
     MeasDist_filt_k_true = zeros(anticipation*(p-1)+1,n_meas_dist);
end

Dm_init = Dm_curr;
if anticipation==0
    if time_i<=1
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dm_curr,Dm_init);
        elseif meas_filter_sel==2
             [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dm_curr,Dm_init,n_w,dist_time_vec,anticipation,p,n_meas_dist);
        end
    else
        if meas_filter_sel==1
            [MeasDist_filt_k_true]=COT_Type1_Filter(alpha_d,Dm_curr,MeasDist_filt_k_true);
        elseif meas_filter_sel==2
            [MeasDist_filt_k_true]=Type_2_Filter_outer(alpha_d,Dm_curr,MeasDist_filt_k_true,n_w,dist_time_vec,anticipation,p,n_meas_dist);
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
    Dm_past_vec = 0*[Dm_init;Dm_init;Dm_init];
end

Dm_past_vec=circshift(Dm_past_vec,1);
Dm_past_vec(1,:)=MeasDist_filt_k_true(1,:);
Dm_filt_minus1=Dm_past_vec(2,:)';
Dm_filt_minus2=Dm_past_vec(3,:)';

%% Type-I Filter for Reference Trajectory

filt_ref=zeros(p,1);
if time_i==1
    filt_ref_val=y_curr;
end

for p_iter=1:p
    filt_ref_val = diag(alpha_r)*filt_ref_val+(eye(ny)-diag(alpha_r))*y_sp(CM+1);
    filt_ref((p_iter-1)*ny+1:p_iter*ny,1)=filt_ref_val;
end
% filt_ref_val=filt_ref(1);

if CM < ST*60 || CM > ET*60
    filt_ref_val = y_curr;
else
    filt_ref_val=filt_ref(1);
end


%%
uk=control_gen_modified(y_curr, u_last, Dm_curr, filt_ref,...
    Unmeas_DistModel, degree, slack_var, anticipation, Weights,...
    limits, setpoints, p, m, nud, fa, D_forecast_filt,...
    D_forecast_unfilt, Dm_filt_minus1, Dm_filt_minus2, mld_copy,...
    model_choice, ts);

% First iteration flag update.
if time_i == 1
    time_i = 2;
end
uk = uk+mean_val(2);

% Only activate in daylight hours
if CM >= ST*60 && CM <= ET*60
    Flowrates = [uk, 0];
else
    Flowrates = [0, 0];
end

Controller_name = 'MoD MPC';

end