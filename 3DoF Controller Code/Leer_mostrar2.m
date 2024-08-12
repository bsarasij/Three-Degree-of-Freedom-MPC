%% Data representation script
% Pablo Otálora Berenguel
%clear
clc
% close all
warning off

%% Data loading
% Select day, month and year.
day_val = 30;
month_val = 05;
year_val = 24;
RW = 5;
file_name = ['Registro_Reactores_' sprintf('%02i', day_val) '_'...
    sprintf('%02i', month_val) '_' sprintf('%02i', year_val) '.csv'];

% Load the dataset with the given options (see function below).
B = leer_registro(file_name);

% Meas_Dist_26 = [B.TempRW5 B.RadGlobal];
% save Meas_Dist_26.mat Meas_Dist_26

% inp_out_db_26_5=[B.pHRW52b(7*60:21*60) B.CaudalCO2RW5(7*60:21*60) B.TempRW5(7*60:21*60) B.RadGlobal(7*60:21*60)];
% inp_out_db_19_5_ms = inp_out_db_19_5-mean(inp_out_db_19_5);
% 
% save inp_out_db_26_5.mat inp_out_db_26_5
% save inp_out_db_14_6_ms.mat inp_out_db_14_6_ms



%% Filter SP

unf_SP = B.pHRW52b;
unf_SP(7*60:min(21*60,length(unf_SP))) = 8;
alpha_r=exp(-1./20);
ny = 1;
filt_ref_val = unf_SP(1);
for p_iter=6:length(B.Tiempo)
    if p_iter < 7*60 || p_iter > 21*60
        filt_ref_val = unf_SP(p_iter);
    else
        filt_ref_val = diag(alpha_r)*filt_ref_val+(eye(ny)-diag(alpha_r))*unf_SP(p_iter-5);
    end
    filt_ref((p_iter-1)*ny+1:p_iter*ny,1)=filt_ref_val;
end

%% Plots
% To plot more than 3 variables at a time, please change the arguments of
% 'tiledlayout' to fit the number of variables.
figure
tiledlayout(6, 1, 'Padding', 'tight', 'TileSpacing', 'tight');

if RW == 5
    % pH
    nexttile
    stairs(B.Tiempo, B.pHRW52b, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(pH2)), pH2, 'LineWidth', 2);
    %stairs(B.Tiempo, unf_SP, '--', 'LineWidth', 2);
    stairs(B.Tiempo, filt_ref(1:length(B.Tiempo)), 'LineWidth', 2);
    grid on
    ylabel('pH [-]');
    ylim([7.4 8.6]);
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);
    legend('pH', 'SP');
    
    % CO2 injection
    nexttile
    stairs(B.Tiempo, B.CaudalCO2RW5, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(CO22)), CO22, 'LineWidth', 2);
    grid on
    ylabel('CO_2 [L/min]');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);
    
    %Temperature
    nexttile
    stairs(B.Tiempo, B.TempRW5, 'LineWidth', 2);
    grid on
    ylabel('T [C]');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);

    % Global irradiance
    nexttile
    stairs(B.Tiempo, B.RadGlobal, 'LineWidth', 2);
    grid on
    ylabel('I_0 [W/m^2]');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);

    %Dillution
    nexttile
    stairs(B.Tiempo, B.SubirNivelRW5, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo, B.BajarNivelRW5, 'LineWidth', 2);
    grid on
    ylabel('Qd [-]');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);

    %Air injection
    nexttile
    stairs(B.Tiempo, B.CaudalAireRW5, 'LineWidth', 2);
    grid on
    ylabel('Qair [L/min]');
    xlabel('Time');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);
else
    % pH
    nexttile
    stairs(B.Tiempo, B.pHRW62b, 'LineWidth', 2);
    hold on;
    yline(8, '--', 'LineWidth', 2);
    grid on
    ylabel('pH [-]');
    ylim([7.6 9]);
    
    % CO2 injection
    nexttile
    stairs(B.Tiempo, B.CaudalCO2RW6, 'LineWidth', 2);
    grid on
    ylabel('CO_2 [L/min]');
    
    %Temperature
    nexttile
    stairs(B.Tiempo, B.TempRW6, 'LineWidth', 2);
    grid on
    ylabel('T [C]');
    xlabel('Time');
end
%%
figure
tiledlayout(4, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile
    stairs(B.Tiempo, B.pHRW52b, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(pH2)), pH2, 'LineWidth', 2);
    %stairs(B.Tiempo, unf_SP, '--', 'LineWidth', 2);
    stairs(B.Tiempo, filt_ref(1:length(B.Tiempo)), 'LineWidth', 2);
    grid on
    ylabel('pH [-]');
    ylim([7.4 8.6]);
    xlim([B.Tiempo(6*60) B.Tiempo(23*60)]);
    legend('pH', 'SP');
    
    % CO2 injection
    nexttile
    stairs(B.Tiempo, B.CaudalCO2RW5, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(CO22)), CO22, 'LineWidth', 2);
    grid on
    ylabel('CO_2 [L/min]');
    xlim([B.Tiempo(6*60) B.Tiempo(23*60)]);
    
    %Temperature
    nexttile
    stairs(B.Tiempo, B.TempRW5, 'LineWidth', 2);
    grid on
    ylabel('T [C]');
    xlim([B.Tiempo(6*60) B.Tiempo(23*60)]);

    % Global irradiance
    nexttile
    stairs(B.Tiempo, B.RadGlobal, 'LineWidth', 2);
    grid on
    ylabel('I_0 [W/m^2]');
    xlim([B.Tiempo(6*60) B.Tiempo(23*60)]);




%%

figure
tiledlayout(2, 1, 'Padding', 'tight', 'TileSpacing', 'tight');

    % pH
    nexttile
    stairs(B.Tiempo, B.pHRW52b, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(pH2)), pH2, 'LineWidth', 2);
    %stairs(B.Tiempo, unf_SP, '--', 'LineWidth', 2);
    stairs(B.Tiempo, filt_ref(1:length(B.Tiempo)), 'LineWidth', 2);
    grid on
    ylabel('pH [-]');
    ylim([7.6 9]);
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);
    legend('pH', 'SP');
    
    % CO2 injection
    nexttile
    stairs(B.Tiempo, B.CaudalCO2RW5, 'LineWidth', 2);
    hold on;
    %stairs(B.Tiempo(1:length(CO22)), CO22, 'LineWidth', 2);
    grid on
    ylabel('CO_2 [L/min]');
    xlim([B.Tiempo(6*60) B.Tiempo(end)]);

%% Other functions

% Function for reading and resampling the data.
function A = leer_registro(nombre, variable_selec)
    % Please change these variables if any other variable is needed. It can
    % also be included as an argument to 'leer_registro'.
    if ~exist('variable_selec', 'var')
        variable_selec = {'Tiempo','pHRW52a','pHRW52b','ODRW52a',...
            'ODRW52b','CaudalCO2RW5','CaudalAireRW5','NivelRW5',...
            'TempRW5','RadGlobal','RadPAR','TempAmbiente',...
            'SubirNivelRW5','BajarNivelRW5',...
            'pHRW62a','pHRW62b','ODRW62a',...
            'ODRW62b','CaudalCO2RW6','CaudalAireRW6','NivelRW6',...
            'TempRW6','SubirNivelRW6','BajarNivelRW6','LitrosCosechadosRW5'};
    end

    % Reading the dataset.
    A = readtable(nombre);

    % Creating the new time column
    A.Tiempo = A.Fecha + A.Hora;

    % Removing the undesired variables and samples without data.
    A = A(:, variable_selec);
    A = A(~any(ismissing(A), 2), :);

    % Retiming the table to one minute sample time.
    A = table2timetable(A);
    A = sortrows(A);
    A = retime(A, 'minutely', 'median');
    A = fillmissing(A, 'linear');
end




