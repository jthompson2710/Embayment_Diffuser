%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 -- James O Thompson (Baylor University)
% Contributors: Chelsea Allison, Ben Andrews, Kenny Befus, Ben Black, 
% Anna Ruefer
% 1D diffusion of H2O and CO2 through an Embayment determined using
% Gradual Declining Fit to Convergence solution
%
%Requires MATLAB 2019a or newer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    runformat = 1; %change to 0 if want to add numbers in the editor below; 1 for interactive dialog boxes
    saveout = 1; %change to 0 is no, 1 is yes
    match = 2; %change to 0 is H2O, 1 is CO2, 2 is both

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if runformat == 0
%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

%%% FILENAME OF EMBAYMENT FILE
    filename = '1D_Example.xlsx';
    outputname = '1D_Example_Test';
    
%%% H2O
    H2O_Concentration_in_Embayment = 3.0;   %Concentration in embayment prior to eruption
    
%%% CO2
    CO2_Concentration_in_Embayment = 500;    %Concentration in embayment prior to eruption

%%% ERUPTION TEMPERATURE 
    TempC = 800;                            %Temperature in Celsius
    TempK = TempC+273;

%%% INITIAL EXSOLVED GAS CONTENT IN WT%
    init_exs = 2.0;                          %Exsolved Gas Content Prior to Eruption (wt%)

%%% Instrument Data ERRORS
    %H2O (percentage)
    H2Oerror = 10.0;
    %CO2 (percentage)
    CO2error = 10.0;

%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

end 
if runformat == 1
%%% FILENAME OF EMBAYMENT FILE
    uiwait(msgbox({'Welcome to this 1D embayment diffusion model';'';...
        'If you have any issues please do not hesitate to contact us';'';...
        'Thank you for using our model and we hope you find it helpful';'';...
        '“All models are wrong, but some are useful” George E. P. Box';''}, 'WELCOME'));
    
    %%% FILENAME OF EMBAYMENT FILE
    [filename,filepath] = uigetfile({'*.xlsx';'*.csv'}, 'Please select Instrument Data Embayment File (.xslx or .csv)');

%%% DIALOG INPUT VALUES
    
    prompt = {'Output Filename:','Enter H_{2}O Concentration in Embayment Prior to Eruption (wt%):',...
        'Enter CO_{2} Concentration in Embayment Prior to Eruption (PPM):',...
        'Enter Eruption Temperature (\circ C):', 'Enter Initial Exsolved Gas Content Prior to Eruption (wt%):',...
        'H_{2}O Instrument Data Error (%):', 'CO_{2} Instrument Data Error (%):'};
    dlgtitle = 'Input Paramenters';
    definput = {'1D_Example_Test','3.0','500.0','800.0','2.0','10.0','10.0'};
    opts.Interpreter = 'tex';
    dims = [1 65];
    input_answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    %%assign user input values to variables
    outputname = input_answer{1};     %output filename
    H2O_Concentration_in_Embayment = str2double(input_answer{2});     %Concentration in embayment prior to eruption
    CO2_Concentration_in_Embayment = str2double(input_answer{3});      %Concentration in embayment prior to eruption
    TempC = str2double(input_answer{4}); %temperature
    TempK = TempC+273; %temperature
    init_exs = str2double(input_answer{5}); %exs content
    H2Oerror = str2double(input_answer{6}); %h20 error 
    CO2error = str2double(input_answer{7}); %co2 error 

end
    %%
% %Load emabyment file and extractaxes, position on mount as columns of data into Matlab
    Validation_Data = readmatrix(filename);
    file_x = Validation_Data(:,1); %second column x position
    vqH2O = (Validation_Data(:,2));
    vqCO2 = (Validation_Data(:,3));

%%% CHECK MODELING POSSIBLE
     if any(vqH2O) == 0 | any(isnan(vqH2O)) == 1 
         uiwait(msgbox({'H_{2}O values are not able to be zero, NaN, or left empty.';'';...
             'There must be H_{2}O data for the model to run and any zero or missing data';...
             'need to be changed to an interpolated value between the neighboring values'...
             ;''}, 'CRITICAL ERROR'));
        return
     end
     
     if any(isnan(vqCO2)) == 1 
        uiwait(msgbox({'CO_{2} values are not able to be NaN or left empty';'';...
            'The model will run without CO_{2} data, however, the column must be filled with zeros';...
            'It is also recommend that any missing CO_{2} data are not filled with zeros.';...
            'It is better if an interpolated value between the neighboring values is used instead'...
            ;''}, 'CRITICAL ERROR'));
        return
     end
     
    %build arrays 
    x_length_actual = max(file_x)-min(file_x);
    x_length = length(file_x)/length(find(file_x==max(file_x)));
    dx = x_length_actual / (x_length-1);
    xq = file_x;
    minvqH2O = min(vqH2O);
    minvqCO2 = min(vqCO2);
    vqH2Oerror = vqH2O * (H2Oerror/100.0);
    vqCO2error = vqCO2 * (CO2error/100.0);
    
    %build the masks
    mask = zeros(x_length,1);
      
    conduit_mask = mask;
    conduit_mask(1,1) = 10;
    emb_mask = mask;
    emb_mask(2:end,1) = 1;
    emb_diff_mask = conduit_mask+emb_mask;
    conduit_mask(1,1) = 1;
%%
    tic
%%% PRELIMINARY SETUP
    
%%% CONSTANTS
    MW_H2O = 18.01528;
    MW_CO2 = 44.0095;

%%% SOLUBILITY, PRESSURE, AND DEPTH SET UP
    F_MDensity = 2380.0; % density in kg/m3
    % based on Liu et al. 2005
    H2Oa1 = 354.94;
    H2Oa2 = 9.623;
    H2Oa3 = -1.5223;
    H2Oa4 = 0.0012439;
    H2Oa5 = -0.0001084;
    H2Oa6 = -0.00001362;
    CO2a1 = 5668.0;
    CO2a2 = 55.99;
    CO2a3 = 0.4133;
    CO2a4 = 0.002041;
    % possible pressures and molar fractions
    H2Od = zeros(1000,5000);
    CO2d = zeros(1000,5000);
    TGPa = ((1:5000)/10000);
    H2OdXf = (1:1000)/1000;
    CO2dXf = flip((1:1000)/1000);
    % calulate master refernce table of solubility at variety of P and X
    for j=1:length(H2OdXf)
        H2Od(j,:) = ( ( (H2Oa1.*((H2OdXf(j).*TGPa.*1e3).^(0.5))) + (H2Oa2.*(H2OdXf(j).*TGPa.*1e3)) + (H2Oa3.*((H2OdXf(j).*TGPa.*1e3).^(1.5))) ) ...
            ./ TempK ) + (H2Oa4.*((H2OdXf(j).*TGPa.*1e3).^(1.5))) + ( (CO2dXf(j).*TGPa.*1e3).*( (H2Oa5.*((H2OdXf(j).*TGPa.*1e3).^(0.5))) ...
            + (H2Oa6.*(H2OdXf(j).*TGPa.*1e3)) ) ) ;
    
        CO2d(j,:) = ( ((CO2dXf(j).*TGPa.*1e3).* ( CO2a1 - (CO2a2.*(H2OdXf(j).*TGPa.*1e3)) ))./TempK ) ...
            + ((CO2dXf(j).*TGPa.*1e3).* ( (CO2a3.*((H2OdXf(j).*TGPa.*1e3).^(0.5))) + ( CO2a4.*((H2OdXf(j).*TGPa.*1e3).^(1.5))) ) );
    end  
    % initial pressure (embayment) based on h2o based on Liu et al. 2005
    [ d_i, x_i ] = min( transpose(abs(H2Od-H2O_Concentration_in_Embayment)) );
    % initial pressure (embayment) based on CO2 based on Liu et al. 2005
    [ cd_i, cx_i ] = min( transpose(abs(CO2d-CO2_Concentration_in_Embayment)) );
    % find where concentrations are the same
    diff_x_i = x_i - cx_i;
    [dMPa, ThousXH2O] = min(abs(diff_x_i));
    pre_vector_index = ThousXH2O;
    % determine X for that concentration
    i_H2OdXf = ThousXH2O / 1000;
    i_CO2dXf = 1-i_H2OdXf;
    % determine pressure and convert from hundreds of MPa to GPa
    GPa_i = x_i(ThousXH2O) / 10000;  
    GPa = GPa_i;
    % calculate inital depth based on pressure
    i_depth = (GPa_i *1e9) / ( F_MDensity * 9.81 );
    
    % Quench pressure (embayment) based on h2o based on Liu et al. 2005
    [ ~, quench_x_i ] = min( transpose(abs(H2Od-minvqH2O)) );
    % Quench pressure (embayment) based on CO2 based on Liu et al. 2005
    [ ~, quench_cx_i ] = min( transpose(abs(CO2d-minvqCO2)) );
    % find where Quench concentrations are the same
    quench_diff_x_i = quench_x_i - quench_cx_i;
    [~, quench_ThousXH2O] = min(abs(quench_diff_x_i));
    % determine Quench pressure and convert from hundreds of MPa to GPa
    GPa_quench = quench_x_i(quench_ThousXH2O) / 10000;  
    % calculate Quench depth based on pressure
    quench_depth = (GPa_quench *1e9) / ( F_MDensity * 9.81 );   
    
    %use fluid composition to determine mass H2O&CO2 in exsolved fraction
    %grams excess vapor w/dissolved volatiles
    i_wtV = init_exs / (1-(0.01*init_exs));
    %mass of H2O & CO2 in exsolved pahse per 100g melt w/dissolved volatiles    
    i_massCO2v = ((i_CO2dXf*MW_CO2) / ((i_H2OdXf*MW_H2O) + (i_CO2dXf*MW_CO2)) )*i_wtV;
    i_massH2Ov = ((i_H2OdXf*MW_H2O) / ((i_H2OdXf*MW_H2O) + (i_CO2dXf*MW_CO2)) )*i_wtV;
    %total H2O & CO2 ("total" includes both dissolved and exsolved)
    i_H2Ot = H2O_Concentration_in_Embayment+i_massH2Ov; % wt%
    i_CO2t = CO2_Concentration_in_Embayment+(10000*i_massCO2v); % ppm
    
    %%
%%% FILL THE MESH AND DOMAIN (1D)
    Endy = x_length;
   %%H2O
    iH2O_C_Con = conduit_mask .* H2O_Concentration_in_Embayment; %Concentration_in_Conduit;
    iH2O_C_Emb = emb_mask .* H2O_Concentration_in_Embayment;
    iH2O_C = iH2O_C_Con + iH2O_C_Emb;
   %%CO2
    iCO2_C_Con = conduit_mask .* CO2_Concentration_in_Embayment; %Concentration_in_Conduit;
    iCO2_C_Emb = emb_mask .* CO2_Concentration_in_Embayment;
    iCO2_C = iCO2_C_Con + iCO2_C_Emb;

%%   
%%% Ascent variation - setup best match
COMP_vqH2O = vqH2O.*emb_mask;
COMP_vqCO2 = vqCO2.*emb_mask;
MC_Ascent_Duration = 0;
Converge20 = -1;

while Converge20 < 0

 %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL   
    time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_C = iH2O_C;
    CO2_C = iCO2_C;
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    r_x = xq+(dx/2);   
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
   
[H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);
 
 COMP_H2O_C = H2O_C.*emb_mask;
 COMP_CO2_C = CO2_C.*emb_mask;

%determine least squares RMS match
 D_H2O_C = sqrt(mean(mean((COMP_H2O_C - COMP_vqH2O).^2))) ./ mean(COMP_vqH2O);
 D_CO2_C = sqrt(mean(mean((COMP_CO2_C - COMP_vqCO2).^2))) ./ mean(COMP_vqCO2);
 
 if MC_Ascent_Duration == 0
     D_H2O_C_All = D_H2O_C;
     D_CO2_C_All = D_CO2_C;
     D_Both_C_All = D_H2O_C+D_CO2_C;
     Converge20 = -1;
 else
     D_H2O_C_All = cat(1, D_H2O_C_All, D_H2O_C);
     D_CO2_C_All = cat(1, D_CO2_C_All, D_CO2_C);
     D_Both_C_All = cat(1, D_Both_C_All, D_H2O_C+D_CO2_C);
     if match == 0
        Converge20 = D_H2O_C_All(end) - D_H2O_C_All(end-1);
     elseif match == 1
        Converge20 = D_CO2_C_All(end) - D_CO2_C_All(end-1);
     else
        Converge20 = D_Both_C_All(end) - D_Both_C_All(end-1);
     end
 end
 
 disp(MC_Ascent_Duration)
 
 MC_Ascent_Duration = MC_Ascent_Duration+20;
end

Ascent_Duration20 = MC_Ascent_Duration-20;
MC_Ascent_Duration = Ascent_Duration20;
Converge5 = -1;

while Converge5 < 0

 %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
    time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_C = iH2O_C;
    CO2_C = iCO2_C;
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    r_x = xq+(dx/2);   
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
   
    
[H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);
 
 COMP_H2O_C = H2O_C.*emb_mask;
 COMP_CO2_C = CO2_C.*emb_mask;

%determine least squares RMS match
 D_H2O_C = sqrt(mean(mean((COMP_H2O_C - COMP_vqH2O).^2))) ./ mean(COMP_vqH2O);
 D_CO2_C = sqrt(mean(mean((COMP_CO2_C - COMP_vqCO2).^2))) ./ mean(COMP_vqCO2);
 
 if MC_Ascent_Duration == Ascent_Duration20
     D_H2O_C_All = D_H2O_C;
     D_CO2_C_All = D_CO2_C;
     D_Both_C_All = D_H2O_C+D_CO2_C;
     Converge5 = -1;
 else
     D_H2O_C_All = cat(1, D_H2O_C_All, D_H2O_C);
     D_CO2_C_All = cat(1, D_CO2_C_All, D_CO2_C);
     D_Both_C_All = cat(1, D_Both_C_All, D_H2O_C+D_CO2_C);
     if match == 0
        Converge5= D_H2O_C_All(end) - D_H2O_C_All(end-1);
     elseif match == 1
        Converge5 = D_CO2_C_All(end) - D_CO2_C_All(end-1);
     else
        Converge5 = D_Both_C_All(end) - D_Both_C_All(end-1);
     end
 end
 
 disp(MC_Ascent_Duration)
 
 MC_Ascent_Duration = MC_Ascent_Duration-5;
  if MC_Ascent_Duration < 0
    MC_Ascent_Duration = -5;
    break
 end
end

Ascent_Duration5 = MC_Ascent_Duration+5;
MC_Ascent_Duration = Ascent_Duration5;
Converge1 = -1;

while Converge1 < 0

 %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
    time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_C = iH2O_C;
    CO2_C = iCO2_C;
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    r_x = xq+(dx/2);   
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
   
    
[H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);
 
 COMP_H2O_C = H2O_C.*emb_mask;
 COMP_CO2_C = CO2_C.*emb_mask;

%determine least squares RMS match
 D_H2O_C = sqrt(mean(mean((COMP_H2O_C - COMP_vqH2O).^2))) ./ mean(COMP_vqH2O);
 D_CO2_C = sqrt(mean(mean((COMP_CO2_C - COMP_vqCO2).^2))) ./ mean(COMP_vqCO2);
 
 if MC_Ascent_Duration == Ascent_Duration5
     D_H2O_C_All = D_H2O_C;
     D_CO2_C_All = D_CO2_C;
     D_Both_C_All = D_H2O_C+D_CO2_C;
     Converge1 = -1;
 else
     D_H2O_C_All = cat(1, D_H2O_C_All, D_H2O_C);
     D_CO2_C_All = cat(1, D_CO2_C_All, D_CO2_C);
     D_Both_C_All = cat(1, D_Both_C_All, D_H2O_C+D_CO2_C);
     if match == 0
        Converge1 = D_H2O_C_All(end) - D_H2O_C_All(end-1);
     elseif match == 1
        Converge1 = D_CO2_C_All(end) - D_CO2_C_All(end-1);
     else
        Converge1 = D_Both_C_All(end) - D_Both_C_All(end-1);
     end
 end
 
 disp(MC_Ascent_Duration)
 
 MC_Ascent_Duration = MC_Ascent_Duration+1;
end

Ascent_Duration1 = MC_Ascent_Duration-1;
MC_Ascent_Duration = Ascent_Duration1;
Converge01 = -1;

while Converge01 < 0

 %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
    time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_C = iH2O_C;
    CO2_C = iCO2_C;
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    r_x = xq+(dx/2);   
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
   
    
[H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);
 
 COMP_H2O_C = H2O_C.*emb_mask;
 COMP_CO2_C = CO2_C.*emb_mask;

%determine least squares RMS match
 D_H2O_C = sqrt(mean(mean((COMP_H2O_C - COMP_vqH2O).^2))) ./ mean(COMP_vqH2O);
 D_CO2_C = sqrt(mean(mean((COMP_CO2_C - COMP_vqCO2).^2))) ./ mean(COMP_vqCO2);
 
 if MC_Ascent_Duration == Ascent_Duration1
     D_H2O_C_All = D_H2O_C;
     D_CO2_C_All = D_CO2_C;
     D_Both_C_All = D_H2O_C+D_CO2_C;
     Converge01 = -1;
 else
     D_H2O_C_All = cat(1, D_H2O_C_All, D_H2O_C);
     D_CO2_C_All = cat(1, D_CO2_C_All, D_CO2_C);
     D_Both_C_All = cat(1, D_Both_C_All, D_H2O_C+D_CO2_C);
     if match == 0
        Converge01 = D_H2O_C_All(end) - D_H2O_C_All(end-1);
     elseif match == 1
        Converge01 = D_CO2_C_All(end) - D_CO2_C_All(end-1);
     else
        Converge01 = D_Both_C_All(end) - D_Both_C_All(end-1);
     end
 end
 
 disp(MC_Ascent_Duration)
 
 MC_Ascent_Duration = MC_Ascent_Duration-0.1;     
 if MC_Ascent_Duration <= 0
     MC_Ascent_Duration = -0.1;
     break
 end
end

Ascent_Duration01 = MC_Ascent_Duration+0.2; %0.2 as answer is previous time step
MC_Ascent_Duration = Ascent_Duration01;

    if MC_Ascent_Duration < 10
        
        if MC_Ascent_Duration <= 0.1
            MC_Ascent_Duration = 0.01;
            Ascent_Duration01 = MC_Ascent_Duration;
        else
            MC_Ascent_Duration = Ascent_Duration01-0.1;
            Ascent_Duration01 = MC_Ascent_Duration;
        end
        Converge001 = -1;

    while Converge001 < 0

     %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
        time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
        dt =0.4*(dx^2/time_step_diffusion);
        Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
        Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
        dt_total = 0;

        % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
        H2O_C = iH2O_C;
        CO2_C = iCO2_C;
        H2O_Cnew = H2O_C;       %later this will change
        CO2_Cnew = CO2_C;       %later this will change
        r_x = xq+(dx/2);   
        gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 


    [H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
        i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
        GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
        init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);

     COMP_H2O_C = H2O_C.*emb_mask;
     COMP_CO2_C = CO2_C.*emb_mask;

    %determine least squares RMS match
     D_H2O_C = sqrt(mean(mean((COMP_H2O_C - COMP_vqH2O).^2))) ./ mean(COMP_vqH2O);
     D_CO2_C = sqrt(mean(mean((COMP_CO2_C - COMP_vqCO2).^2))) ./ mean(COMP_vqCO2);

     if MC_Ascent_Duration == Ascent_Duration01
         D_H2O_C_All = D_H2O_C;
         D_CO2_C_All = D_CO2_C;
         D_Both_C_All = D_H2O_C+D_CO2_C;
         Converge001 = -1;
     else
         D_H2O_C_All = cat(1, D_H2O_C_All, D_H2O_C);
         D_CO2_C_All = cat(1, D_CO2_C_All, D_CO2_C);
         D_Both_C_All = cat(1, D_Both_C_All, D_H2O_C+D_CO2_C);
         if match == 0
            Converge001 = D_H2O_C_All(end) - D_H2O_C_All(end-1);
         elseif match == 1
            Converge001 = D_CO2_C_All(end) - D_CO2_C_All(end-1);
         else
            Converge001 = D_Both_C_All(end) - D_Both_C_All(end-1);
         end
     end
     disp(MC_Ascent_Duration)

     MC_Ascent_Duration = MC_Ascent_Duration+0.01;
    end

    Ascent_Duration001 = MC_Ascent_Duration-0.02; %0.2 as answer is previous time step
    MC_Ascent_Duration = Ascent_Duration001;
    
    end

disp('Found Ascent Rate Now Running Final Model!')
 %%% FINAL RUN FOR MODEL
    time_step_diffusion = H2O_Concentration_in_Embayment* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*MC_Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_C = iH2O_C;
    CO2_C = iCO2_C;
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    r_x = xq+(dx/2);      
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
   
    
[H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf);
 
toc
%%
    close()
    %plot all figures
    if (max(vqCO2) ~= 0) && (max(vqH2O) ~= 0) %both data
    
        figure()
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.6 0.8]);
        subplot(1,2,1)
        errorbar(xq, (vqH2O),vqH2Oerror, 'ko','MarkerSize',4, 'MarkerFaceColor','k','CapSize',0,'color', '[.7 .7 .7]')
        xlabel('Position (\mum)')
        ylabel('H_2O (wt%)')
        axis([0 max(xq) min(min(vqH2O))*0.7 max(max(vqH2O))*1.15])
        hold on
        plot(xq, (H2O_C),'color','[0 0.4470 0.7410]','linewidth',2)
        hold off

        subplot(1,2,2)
        errorbar(xq, (vqCO2),vqCO2error, 'ko','MarkerSize',4, 'MarkerFaceColor','k','CapSize',0,'color', '[.7 .7 .7]')
        xlabel('Position (\mum)')
        ylabel('CO_2 (ppm)')
        axis([0 max(xq) min(min(vqCO2))*0.7 max(max(vqCO2))*1.15])
        hold on
        plot(xq,(CO2_C),'color','[0 0.4470 0.7410]','linewidth',2)
        hold off
        leg1 = legend({'Instrument Data','Model'},'Location','southeast');
        leg1.FontSize = 14;

        sgtitle({append('Initial H_2O (wt%) = ',num2str(H2O_Concentration_in_Embayment), ...
            '    Initial CO_2 (ppm) = ',num2str(CO2_Concentration_in_Embayment),...
            '    Inital Exsolved Phase (wt%) = ',num2str(init_exs)),...
            append('Temperature (^{o}C) = ',num2str(TempC), ...
            '    Time Elapsed (hrs) = ',num2str(MC_Ascent_Duration), ...
            '    Decompression rate (MPa/s) = ','10^{', num2str(log10(1000*(Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e9))),'}',...
            append('Starting Pressure (MPa) = ',num2str(1000*GPa),...
            '    Quench Pressure (MPa) = ',num2str(1000*GPa_quench),...
            '    Quench Depth (m) = ',num2str(quench_depth)), ''})

        drawnow
        savefig(append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_1d_profile.fig'))
        saveas(gcf,append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_1d_profile.jpg'))
        
    elseif (max(vqCO2) == 0) && (max(vqH2O) ~= 0) % only H2O
    
        figure()
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
        errorbar(xq, (vqH2O),vqH2Oerror, 'ko','MarkerSize',4, 'MarkerFaceColor','k')
        xlabel('Position (\mum)')
        ylabel('H_2O (wt%)')
        axis([0 max(xq) min(min(vqH2O))*0.7 max(max(vqH2O))*1.15])
        hold on
        plot(xq, (H2O_C),'color','[0 0.4470 0.7410]','linewidth',2)
        hold off
        leg1 = legend({'Instrument Data','Model'},'Location','southeast');
        leg1.FontSize = 14;
        
        title({'',append('Initial H_2O (wt%) = ',num2str(H2O_Concentration_in_Embayment), ...
            '    Initial CO_2 (ppm) = ',num2str(CO2_Concentration_in_Embayment),...
            '    Inital Exsolved Phase (wt%) = ',num2str(init_exs)),...
            append('Temperature (^{o}C) = ',num2str(TempC), ...
            '    Time Elapsed (hrs) = ',num2str(MC_Ascent_Duration), ...
            '    Decompression rate (MPa/s) = ','10^{', num2str(log10(1000*(Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e9))),'}',...
            append('Starting Pressure (MPa) = ',num2str(1000*GPa_i),...
            '    Quench Pressure (MPa) = ',num2str(1000*GPa_quench),...
            '    Quench Depth (m) = ',num2str(quench_depth)), ''})
    
        drawnow
        savefig(append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_1d_profile.fig'))
        saveas(gcf,append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_1d_profile.jpg'))
    end
    
    if saveout == 1
    %write results out 
    init_exs_string = strrep(num2str(init_exs),'.','_');
    header = {'Position (microns)','H2O Instrument Data (wt%)','CO2 Instrument Data (ppm)','H2O Model (wt%)','CO2 Model (ppm)'};
    oned_data = [xq (vqH2O) (vqCO2) (H2O_C) (CO2_C)];
    oned_header_data = [header; num2cell(oned_data)];
    writecell(oned_header_data,append(outputname,'_M0_',init_exs_string,'_1d_profile.xls'))
    
    header1 = {'Inital H2O (wt%)','Inital CO2 (ppm)','Inital Exsolved Phase [M0] (wt%)','Temperature (Celsius)',...
        'Time Elapsed (hrs)','Decompression rate (MPa/s)', 'Starting Pressure (MPa)','Quench depth (m)','Quench Pressure (MPa)'};
    output_data = [H2O_Concentration_in_Embayment CO2_Concentration_in_Embayment init_exs TempC...
        MC_Ascent_Duration ((Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e6) GPa_i*1000 quench_depth GPa_quench*1000 ];
    output_header_data = [header1; num2cell(output_data)];
    writecell(output_header_data,append(outputname,'_M0_',init_exs_string,'_parameters.xls'))
    end
    
%%
%%%Diffusion Function
function [H2O_C, CO2_C] = EMIDIFF_1D_OPERATOR(H2O_C,CO2_C,H2O_Cnew,CO2_Cnew, ...
    i_H2Ot,i_CO2t,MW_H2O,MW_CO2,H2OdXf,H2Od,CO2d,i_depth,dt_total,Ascent_rate,F_MDensity,...
    GPa,TGPa,Max_Time,dx,dt,xq,x_length,r_x,gamma,TempK,emb_diff_mask,...
    init_exs,pre_vector_index,conduit_mask,emb_mask,i_H2OdXf)

    uH2Od_p = max(H2O_C);
    u_H2OdXf = i_H2OdXf;
 %%% MAIN LOOP - Diffusion of Element only%%
    while dt_total < Max_Time

         %%H2Ot
         %H2Ot concentration matrix
         H2Ot_le2 = double(H2O_C <= 2);
         H2Ot_gt2 = double(H2O_C > 2);
         %H2Ot (<2wt%)
         H2O_D2_le2 = ( H2O_C.* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)./TempK)) ).*emb_diff_mask; %D in um2/s in conduit - 10 fold increase in diffusion
         %Diffusion equation from Ni and Zhang (2008), eq.13 or in abstract, applies to H2Ototal less than 2wt.%
         %H2Ot (>2wt%)
         H2OSX = (H2O_C./MW_H2O)./ ( (H2O_C./MW_H2O) + ((100.0-H2O_C)./32.49) );
         H2O_D2_gt2 = ( H2OSX.* exp( 13.470 + (-49.996.*H2OSX) + (7.0827.*(H2OSX.^0.5)) + (1.8875.*GPa) - ...
               ( ( 9532.3 + (-91933.0.*H2OSX) + (13403.0.*(H2OSX.^0.5)) + (3625.6.*GPa) ) ./ TempK ) ) ) .*emb_diff_mask;%
         %Diffusion equation from Ni and Zhang (2008), eq.12 or in abstract
         H2O_D = (H2O_D2_le2.*H2Ot_le2) + (H2O_D2_gt2.*H2Ot_gt2);
                  
         
         %%CO2
         CO2_D = ((exp(-14.34-((17360-(0.6527.*(GPa.*1e3)))./TempK)+(-0.7172+(1436.8./TempK)).*H2O_C))*(10^12)).*emb_diff_mask;
         %CO2 Diffusion equation from Zhang et al. (2007), eq.28 applies to CO2 in all melts, Cw is water concentration
            

      %Simplifying variables in discretization
        
         up1x = x_length-1; low1x = 2; %starting and ending points in diffusion model x
         mult1x = dt./(xq(low1x:up1x).^gamma.*dx.^2);
         rgammax = r_x(low1x:up1x).^gamma; % same as "(r(i))^gamma"
        %%H2O
         H2O_Dmain = H2O_D(low1x:up1x);  %D(i)
         H2O_Cmain = H2O_C(low1x:up1x); %C(i)
         
        %%CO2
         CO2_Dmain = CO2_D(low1x:up1x);  %D(i)
         CO2_Cmain = CO2_C(low1x:up1x); %C(i)
         
     %Discretized Diffusion equation 
        %%H2O
         H2O_diffopperator = mult1x .* (rgammax.*((H2O_Dmain.*H2O_D(low1x+1:up1x+1))./(H2O_Dmain+H2O_D(low1x+1:up1x+1))).*(H2O_C(low1x+1:up1x+1)-H2O_Cmain)...
             -rgammax.*((H2O_Dmain.*H2O_D(low1x-1:up1x-1))./(H2O_Dmain+H2O_D(low1x-1:up1x-1)).*(H2O_Cmain-H2O_C(low1x-1:up1x-1))));
        
        %%CO2
         CO2_diffopperator = mult1x .* (rgammax.*((CO2_Dmain.*CO2_D(low1x+1:up1x+1))./(CO2_Dmain+CO2_D(low1x+1:up1x+1))).*(CO2_C(low1x+1:up1x+1)-CO2_Cmain)...
             -rgammax.*((CO2_Dmain.*CO2_D(low1x-1:up1x-1))./(CO2_Dmain+CO2_D(low1x-1:up1x-1)).*(CO2_Cmain-CO2_C(low1x-1:up1x-1))));
      
     %Updates Concentration 
        %%H2O
         H2O_Cnew(low1x:up1x) = H2O_Cmain + H2O_diffopperator;
         H2O_Cnew(1) = H2O_Cnew(2);
         H2O_Cnew(end) = H2O_Cnew(end-1);
         H2O_C=H2O_Cnew;
       
        %%CO2
         CO2_Cnew(low1x:up1x) = CO2_Cmain + CO2_diffopperator;
         CO2_Cnew(1) = CO2_Cnew(2);
         CO2_Cnew(end) = CO2_Cnew(end-1);
         CO2_C=CO2_Cnew; 

        %%%update solubility
         %update pressure based on ascent
         GPa = ( (i_depth - (dt_total*Ascent_rate)) * ( F_MDensity * 9.81 ) ) / 1e9 ;
         % update pressure and determine dissolved vector for the given
         % pressure
         [ ~, TGPa_x_i ] = min( (abs(TGPa-GPa)) );
         H2Od_GPa = H2Od(:,TGPa_x_i); % wt%
         CO2d_GPa = CO2d(:,TGPa_x_i); % ppm
         
         % calculaete mass for exsolved vectors
         H2Ov_GPa = i_H2Ot - H2Od_GPa; % wt%
         CO2v_GPa = (i_CO2t - CO2d_GPa); % ppm

         % calculaete moles for exsolved vectors
         moleH2Ov_GPa = (H2Ov_GPa ./ MW_H2O); % mole from wt%
         moleCO2v_GPa = ((CO2v_GPa./10000) ./ MW_CO2); % mole from wt%

         % calculaete Xf from moles for exsolved vectors
         XfH2Ov_GPa =  moleH2Ov_GPa ./ (moleH2Ov_GPa + moleCO2v_GPa);

         % find where Xf vector is the same as Xf position for exsolved
         [~, vector_index] = min( abs( XfH2Ov_GPa - transpose(H2OdXf)) );
         if init_exs < 0.0
             if (vector_index > 998) | (vector_index < 2)  
                 vector_index = pre_vector_index*1.000023; 
                 pre_vector_index = vector_index;
                 vector_index = round(pre_vector_index);         
             else
                pre_vector_index = vector_index;
             end

         else
             if (vector_index > 998) | (vector_index < 2) ;  vector_index = pre_vector_index; end
             pre_vector_index = vector_index;
         end
         u_H2OdXf = vector_index/1000;

         %update concentrations
         %H2O
         H2Od_p = H2Od_GPa(vector_index);
         if H2Od_p > uH2Od_p
             H2Od_p = uH2Od_p;
         end
         uH2Od_p = H2Od_p;
         H2Ov_p = i_H2Ot - H2Od_GPa(vector_index);
         %update diss H2O in conduit based on above
         u_H2O_C_Con = conduit_mask .* H2Od_p;
         u_H2O_C_Emb = emb_mask .* H2O_C;
         H2O_C = u_H2O_C_Con + u_H2O_C_Emb ;        
         
         %CO2
         CO2d_p = CO2d_GPa(vector_index);
         CO2v_p = i_CO2t - CO2d_GPa(vector_index);
         %update diss CO2 in conduit based on above
         u_CO2_C_Con = conduit_mask .* CO2d_p;
         u_CO2_C_Emb = emb_mask .* CO2_C;
         CO2_C = u_CO2_C_Con + u_CO2_C_Emb;       
         
         %update time step for next round based on temperature
         dt = 0.4*(dx^2/(mean(H2O_D(2:end-2))));
         dt_total = dt_total+dt;
         
    end
end
%%%References used in this code
    %Zhang, Y., & Ni, H. (2010). Diffusion of H, C, and O components in silicate melts. Reviews in Mineralogy and Geochemistry, 72(1), 171-225.
    %Ni, H., & Zhang, Y. (2008). H2O diffusion models in rhyolitic melt with new high pressure data. Chemical Geology, 250(1-4), 68-78.
    %Zhang, Y., Xu, Z., Zhu, M., & Wang, H. (2007). Silicate melt properties and volcanic eruptions. Reviews of Geophysics, 45(4).
    %Liu, Y., A. T. Anderson, and C. J. N. Wilson (2007). Melt pockets in phenocrysts and decompression rates of silicic magmas before 
    %       fragmentation, J. Geophys. Res., 112, B06204.