%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 -- James O Thompson (Baylor University)
% Contributors: Chelsea Allison, Ben Andrews, Kenny Befus, Ben Black, 
% Anna Ruefer
% 2D diffusion of H2O and CO2 through an Embayment
%
%NOTE: if you want to redraw the embayment mask delete or rename the
% .mat mask file.
%
%Requires MATLAB 2019a or newer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    runformat = 1; %This is the runnig format. Change to 0 if want to add numbers in the editor below; 1 for interactive dialog boxes
    saveout = 1; %This is the saving option. Change to 0 is no, 1 is yes
    maskmatch = 0; %This sets the data used to create the masks. Change to 0 for H2O or 1 for CO2

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
if runformat == 0
%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

%%% FILENAME OF EMBAYMENT FILE
    filename = '2D_Example.xlsx';
    outputname = '2D_Example_Test_S';
    
%%% H2O
    H2O_Concentration_in_Embayment = 3.0;     %Concentration in embayment prior to eruption in wt%
    
%%% CO2
    CO2_Concentration_in_Embayment = 500.0;      %Concentration in embayment prior to eruption in ppm

%%% FRAGMENTATION DEPTH
    quench_depth = 500; %meters
    
%%% ERUPTION TIMESCALE
    Ascent_Duration = 25;      %Hours to run the model 
  
%%% ERUPTION TEMPERATURE 
    TempC = 800;
    TempK = TempC+273;

%%% INITIAL EXSOLVED GAS CONTENT IN WT%
    init_exs = 2.0;

%%% Instrument Data ERRORS
    %H2O (percentage)
    H2Oerror = 10.0;
    %CO2 (percentage)
    CO2error = 10.0;  
    
%000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

end
if runformat == 1
%%% FILENAME OF EMBAYMENT FILE
    uiwait(msgbox({'Welcome to this 2D embayment diffusion model';'';...
        'If you would like any help please do not hesitate to contact us';'';...
        'Thank you for using our model and we hope you find it helful';'';...
        '“All models are wrong, but some are useful” George E. P. Box';''}, 'WELCOME'));
%%% FILENAME OF EMBAYMENT FILE
    [filename,filepath] = uigetfile({'*.xlsx';'*.csv'}, 'Please select Instrument Data Embayment File (.xslx or .csv)');

%%% DIALOG INPUT VALUES
    
    prompt = {'Output Filename:','Enter H_{2}O Concentration in Embayment Prior to Eruption (wt%):',...
        'Enter CO_{2} Concentration in Embayment Prior to Eruption (PPM):', 'Enter Quench Depth (m):', ...
        'Enter Eruption Temperature (\circ C):', 'Enter Initial Exsolved Gas Content Prior to Eruption (wt%):',...
        'Enter Ascent Time (hrs):', 'H_{2}O Instrument Data Error (%):', 'CO_{2} Instrument Data Error (%):'};
    dlgtitle = 'Input Paramenters';
    definput = {'2D_Example_Test_S','3.0','500.0','500','800.0','2.0','25.0','10.0','10.0'};
    opts.Interpreter = 'tex';
    dims = [1 65];
    input_answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    %%assign user input values to variables
    outputname = input_answer{1};     %output filename
    H2O_Concentration_in_Embayment = str2double(input_answer{2});     %Concentration in embayment prior to eruption
    CO2_Concentration_in_Embayment = str2double(input_answer{3});      %Concentration in embayment prior to eruption
    quench_depth = str2double(input_answer{4}); %meters
    TempC = str2double(input_answer{5}); %temperature
    TempK = TempC+273; %temperature
    init_exs = str2double(input_answer{6}); %exs content
    Ascent_Duration = str2double(input_answer{7});      %Hours to run the model, choosing a close value to cooling time is critical for efficient program run 
    H2Oerror = str2double(input_answer{8}); %h20 error 
    CO2error = str2double(input_answer{9}); %co2 error
end

% %Load 2D emabyment concentration file and extract axes, position on mount as columns of data into Matlab
    Validation_Data = readmatrix(filename);
    file_x = Validation_Data(:,1); %second column x position
    file_y = (Validation_Data(:,2));
    V_H2O = (Validation_Data(:,3));
    V_CO2 = (Validation_Data(:,4));
    x_length_actual = max(file_x)-min(file_x); 
    x_length = length(file_x)/length(find(file_x==max(file_x)));
    dx = x_length_actual / (x_length-1);
    y_length_actual = max(file_y)-min(file_y);
    y_length = length(file_y)/length(find(file_y==max(file_y)));
    dy = y_length_actual / (y_length-1);
    
    %build mesh grid and set up data for modeling (removing erroneous data (e.g., negative contents))
    [xq,yq] = meshgrid((min(file_x):dx:max(file_x)), (min(file_y):dy:max(file_y)));
    vqH2O = flipud(griddata(file_x,file_y,V_H2O,xq,yq));
    vqH2O(vqH2O<0)=1e-4;
    vqCO2 = flipud(griddata(file_x,file_y,V_CO2,xq,yq));
    vqCO2(vqCO2<0)=1e-4;
    [vqH2O_max_col_m, vqH2O_max_col_i] = max(vqH2O);   

%%
%%open data to model and draw mask if required
    if isfile(append(outputname,'_mask.mat'))
        emb_diff_mask = load(append(outputname,'_mask.mat'));
        emb_diff_mask = emb_diff_mask.emb_diff_mask;
        
    else
        uiwait(msgbox({'Please a draw a polygon around the entire embayment including the pipe and conduit';...
            'Right-click when done';'You are able to click off the image'},'ROI Polygon of Embayment'));
        figure()
        if maskmatch == 0
            imshow(vqH2O, 'DisplayRange',[min(vqH2O(:)) max(vqH2O(:))], 'Colormap',parula);
            c = colorbar;
            c.Label.String = 'H_{2}O (wt.%)';
        else
            imshow(vqCO2, 'DisplayRange',[min(vqCO2(:)) max(vqCO2(:))], 'Colormap',parula);
            c = colorbar;
            c.Label.String = 'CO_{2} (ppm)';
        end
        title({'Draw a polygon around the entire embayment including the pipe and conduit', ...
            'You are able to click off the image', 'right-click when done'})
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
        set(gca, 'Position', [0.2 0.2 0.7 0.7]);
        dataarea = drawpolygon('Color','green');
        dataarea_mask = int8(createMask(dataarea));
        datapolyline_mask = dataarea.Position;
        datapolyline_mask = [datapolyline_mask;datapolyline_mask(1,:)];
        close()
        uiwait(msgbox({'Please a draw a polyline through the outer neck of the embayment';...
            'The line must extend across the entire image';'You are able to click off the image';...
            'Right-click when done'},'Outer Neck Location'));
        figure()
        if maskmatch == 0
            imshow(vqH2O, 'DisplayRange',[min(vqH2O(:)) max(vqH2O(:))], 'Colormap',parula);
            c = colorbar;
            c.Label.String = 'H_{2}O (wt.%)';
        else
            imshow(vqCO2, 'DisplayRange',[min(vqCO2(:)) max(vqCO2(:))], 'Colormap',parula);
            c = colorbar;
            c.Label.String = 'CO_{2} (ppm)';  
        end
        title({'Draw a polyline through the outer neck of the embayment', ...
            'The line must extend across the entire image','You are able to click off the image',...
            'right-click when done'})
        hold on
        plot(datapolyline_mask(:,1),datapolyline_mask(:,2), 'r','LineWidth',2)
        hold off
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
        set(gca, 'Position', [0.2 0.2 0.7 0.7]);
        neckline = drawpolyline('Color','green');
        neckline_mask = int8(createMask(neckline));
        close()
        calc_mask = dataarea_mask;
        [neck_idx,neck_idy] = find(neckline_mask==1);
        for n=1:length(neck_idy)
            calc_mask(neck_idx(n):end,neck_idy(n)) = 10;
        end
        emb_diff_mask = calc_mask .* dataarea_mask;
        save(append(outputname,'_mask.mat'),'emb_diff_mask');
    end
    
    %build the masks
    emb_diff_mask = double(emb_diff_mask);

    conduit_mask10 = emb_diff_mask;
    conduit_mask10(conduit_mask10<10) = 0;
    conduit_mask = conduit_mask10;
    conduit_mask(conduit_mask==10) = 1;
    emb_mask = emb_diff_mask;
    emb_mask(emb_mask>1) = 0;
    cry_mask = emb_diff_mask;
    cry_mask(cry_mask>1) = 1;
    cry_mask = double(~cry_mask);

    H2O_emb_diff_mask = emb_diff_mask;
    H2O_emb_diff_mask(H2O_emb_diff_mask==0) = 0.001;
    CO2_emb_diff_mask = emb_diff_mask;
    CO2_emb_diff_mask(CO2_emb_diff_mask==0) = 0.001;
        
    emb_mask_sub = smoothdata(emb_mask,1, 'movmean', 6);
    emb_mask_sub = smoothdata(emb_mask_sub,2, 'movmean', 6);
    emb_mask_sub(emb_mask_sub<1) = 0;
%%
tic
%%% PRELIMINARY SETUP
    
%%% CONSTANTS
    MW_H2O = 18.01528;
    MW_CO2 = 44.0095;

%%% SOLUBILITY, PRESSURE, AND DEPTH SET UP
    F_MDensity = 2380.0; % density in kg/m3
    % based on Lui et al. 2005
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
    % initial pressure (embayment) based on h2o based on Lui et al. 2005
    [ d_i, x_i ] = min( transpose(abs(H2Od-H2O_Concentration_in_Embayment)) );
    % initial pressure (embayment) based on CO2 based on Lui et al. 2005
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
    
    % calculate Quench presssure based on depth
    GPa_quench = (quench_depth * ( F_MDensity * 9.81 )) / 1e9;
    
    %use fluid composition to determine mass H2O&CO2 in exsolved fraction
    %grams excess vapor w/dissolved volatiles
    i_wtV = init_exs / (1-(0.01*init_exs));
    %mass of H2O & CO2 in exsolved pahse per 100g melt w/dissolved volatiles    
    i_massCO2v = ((i_CO2dXf*MW_CO2) / ((i_H2OdXf*MW_H2O) + (i_CO2dXf*MW_CO2)) )*i_wtV;
    i_massH2Ov = ((i_H2OdXf*MW_H2O) / ((i_H2OdXf*MW_H2O) + (i_CO2dXf*MW_CO2)) )*i_wtV;
    %total H2O & CO2 ("total" includes both dissolved and exsolved)
    i_H2Ot = H2O_Concentration_in_Embayment+i_massH2Ov; % wt%
    i_CO2t = CO2_Concentration_in_Embayment+(10000*i_massCO2v); % ppm
    
%%% FILL THE MESH AND DOMAIN (2D)
    Endx = y_length;
    Endy = x_length;
   %%H2O
    H2O_C_Con = conduit_mask .* H2O_Concentration_in_Embayment ; %Concentration_in_Conduit;
    H2O_C_Emb = emb_mask .* H2O_Concentration_in_Embayment ;
    H2O_C_Cry = cry_mask .* (vqH2O*0.9);
    H2O_C = H2O_C_Con + H2O_C_Emb + H2O_C_Cry  ;
   %%CO2
    CO2_C_Con = conduit_mask .* CO2_Concentration_in_Embayment ; %Concentration_in_Conduit;
    CO2_C_Emb = emb_mask .* CO2_Concentration_in_Embayment ;
    CO2_C_Cry = cry_mask .* (vqCO2*0.9);
    CO2_C = CO2_C_Con + CO2_C_Emb + CO2_C_Cry ;
   
 %%% TIME VARIABLES AND MANIPULATIONS FOR MODEL
    time_step_diffusion = H2O_Concentration_in_Embayment * exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)/TempK)); 
    dt =0.4*(dx^2/time_step_diffusion);
    Max_Time = 3600*Ascent_Duration;  % seconds %this is used as a max time, can be raised if needed
    Ascent_rate = (i_depth - quench_depth) / Max_Time; % in meters per second
    dt_total = 0;

    % %MODEL HOW THE DISTURBANCE MIGRATES THRU TIME%%
    H2O_Cnew = H2O_C;       %later this will change
    CO2_Cnew = CO2_C;       %later this will change
    dr = dx;
    r_x = xq+(dx/2);   
    r_y = yq+(dy/2);   
    iter=1 ;
    z=1;
    iter1=1 ;
    z1=1;
    gamma = 2; %Gamma - will be used to change between linear = 0; cylindrical = 1; spherical=2 %treatments; usually gamma will be set to 2 
    uH2Od_p = max(H2O_C);
    u_H2OdXf = i_H2OdXf;
%%
 %%% MAIN LOOP - Diffusion of Element only%%
    figure()
    while dt_total < Max_Time
               
         %%H2Ot
         %H2Ot concentration matrix
         H2Ot_le2 = double(H2O_C <= 2);
         H2Ot_gt2 = double(H2O_C > 2);
         %H2Ot (<2wt%)
         H2O_D2_le2 = ( H2O_C.* exp((9.5279 + 1.8875*GPa -(9698.5+3625.6*GPa)./TempK)) ).*H2O_emb_diff_mask; %D in um2/s in conduit - 10 fold increase in diffusion
         %Diffusion equation from Ni and Zhang (2008), eq.13 or in abstract, applies to H2Ototal less than 2wt.%
         %H2Ot (>2wt%)
         H2OSX = (H2O_C./MW_H2O)./ ( (H2O_C./MW_H2O) + ((100.0-H2O_C)./32.49) );
         H2O_D2_gt2 = ( H2OSX.* exp( 13.470 + (-49.996.*H2OSX) + (7.0827.*(H2OSX.^0.5)) + (1.8875.*GPa) - ...
               ( ( 9532.3 + (-91933.0.*H2OSX) + (13403.0.*(H2OSX.^0.5)) + (3625.6.*GPa) ) ./ TempK ) ) ) .*H2O_emb_diff_mask;%
         %Diffusion equation from Ni and Zhang (2008), eq.12 or in abstract
         H2O_D = (H2O_D2_le2.*H2Ot_le2) + (H2O_D2_gt2.*H2Ot_gt2);

         %%CO2
         CO2_D = ((exp(-14.34-((17360-(0.6527.*(GPa.*1e3)))./TempK)+(-0.7172+(1436.8./TempK)).*H2O_C))*(10^12)).*CO2_emb_diff_mask;
         %CO2 Diffusion equation from Zhang et al. (2007), eq.28 applies to CO2 in all melts, Cw is water concentration            

      %%%Simplifying variables in discretization
         up1x = x_length-1; low1x = 2; %starting and ending points in diffusion model x
         up1y = y_length-1; low1y = 2; %starting and ending points in diffusion model y 
         mult1x = dt./(xq(low1y:up1y,low1x:up1x).^gamma.*dx.^2);
         mult1y = dt./(yq(low1y:up1y,low1x:up1x).^gamma.*dy.^2);   %exactly the same as "dt/(x(i)^gamma*dx^2)" in my main equation
         mult1 = sqrt(mult1x.*mult1y);
         rgammax = r_x(low1y:up1y,low1x:up1x).^gamma; % same as "(r(i))^gamma"
         rgammay = r_y(low1y:up1y,low1x:up1x).^gamma;
         rgamma = sqrt(rgammax.*rgammay); 
        %%H2O
         H2O_Dmain = H2O_D(low1y:up1y,low1x:up1x);  %D(i)
         H2O_Cmain = H2O_C(low1y:up1y,low1x:up1x); %C(i)
         
        %%CO2
         CO2_Dmain = CO2_D(low1y:up1y,low1x:up1x);  %D(i)
         CO2_Cmain = CO2_C(low1y:up1y,low1x:up1x); %C(i)
         
     %Discretized Diffusion equation  - static varibale discretization
        %%H2O
         H2O_diffopperator =  ( (mult1x .* ( (rgammax.*((H2O_Dmain.*H2O_D(low1y:up1y,low1x+1:up1x+1))./(H2O_Dmain+H2O_D(low1y:up1y,low1x+1:up1x+1))).*(H2O_C(low1y:up1y,low1x+1:up1x+1)-H2O_Cmain)...
             -rgammax.*((H2O_Dmain.*H2O_D(low1y:up1y,low1x-1:up1x-1))./(H2O_Dmain+H2O_D(low1y:up1y,low1x-1:up1x-1)).*(H2O_Cmain-H2O_C(low1y:up1y,low1x-1:up1x-1)))) )) + ... 
                                      (mult1y .* ( (rgammay.*((H2O_Dmain.*H2O_D(low1y+1:up1y+1,low1x:up1x))./(H2O_Dmain+H2O_D(low1y+1:up1y+1,low1x:up1x))).*(H2O_C(low1y+1:up1y+1,low1x:up1x)-H2O_Cmain)...
             -rgammay.*((H2O_Dmain.*H2O_D(low1y-1:up1y-1,low1x:up1x))./(H2O_Dmain+H2O_D(low1y-1:up1y-1,low1x:up1x)).*(H2O_Cmain-H2O_C(low1y-1:up1y-1,low1x:up1x)))) )) );
         H2O_diffopperator(isnan(H2O_diffopperator))=0;
        %%CO2
         CO2_diffopperator =  ( (mult1x .* ( (rgammax.*((CO2_Dmain.*CO2_D(low1y:up1y,low1x+1:up1x+1))./(CO2_Dmain+CO2_D(low1y:up1y,low1x+1:up1x+1))).*(CO2_C(low1y:up1y,low1x+1:up1x+1)-CO2_Cmain)...
             -rgammax.*((CO2_Dmain.*CO2_D(low1y:up1y,low1x-1:up1x-1))./(CO2_Dmain+CO2_D(low1y:up1y,low1x-1:up1x-1)).*(CO2_Cmain-CO2_C(low1y:up1y,low1x-1:up1x-1)))) )) + ...
                                      (mult1y .* ( (rgammay.*((CO2_Dmain.*CO2_D(low1y+1:up1y+1,low1x:up1x))./(CO2_Dmain+CO2_D(low1y+1:up1y+1,low1x:up1x))).*(CO2_C(low1y+1:up1y+1,low1x:up1x)-CO2_Cmain)...
             -rgammay.*((CO2_Dmain.*CO2_D(low1y-1:up1y-1,low1x:up1x))./(CO2_Dmain+CO2_D(low1y-1:up1y-1,low1x:up1x)).*(CO2_Cmain-CO2_C(low1y-1:up1y-1,low1x:up1x)))) )) );
         CO2_diffopperator(isnan(CO2_diffopperator))=0;
         
     %Updates Concentration 
        %%H2O
         H2O_Cnew(low1y:up1y,low1x:up1x) = H2O_Cmain + H2O_diffopperator;
         H2O_Cnew(:,1) = H2O_Cnew(:,2);
         H2O_Cnew(:,end) = H2O_Cnew(:,end-1);
         H2O_Cnew(1,:) = H2O_Cnew(2,:);
         H2O_Cnew(end,:) = H2O_Cnew(end-1,:);

         H2O_C=H2O_Cnew;
         H2O_Cm=min(H2O_C(H2O_C>0));
         H2O_C(H2O_C<0) = H2O_Cm;
        
        %%CO2
         CO2_Cnew(low1y:up1y,low1x:up1x) = CO2_Cmain + CO2_diffopperator;
         CO2_Cnew(:,1) = CO2_Cnew(:,2);
         CO2_Cnew(:,end) = CO2_Cnew(:,end-1);
         CO2_Cnew(1,:) = CO2_Cnew(2,:);
         CO2_Cnew(end,:) = CO2_Cnew(end-1,:);

         CO2_C=CO2_Cnew; 
         CO2_Cm=min(CO2_C(CO2_C>0));
         CO2_C(CO2_C<0) = CO2_Cm;

      %%%Plotting Part of Model  
              
                H2O_textpositions = [ 0.05  0.1 0.15    0.2 0.25    0.30    0.3500    0.4    0.45    0.5    0.55    0.6   0.8 0.85    .9 ];
                H2O_ymax = H2O_Concentration_in_Embayment + 1; % ymin = 0;
                H2O_textloc = H2O_ymax.*H2O_textpositions;

                CO2_textpositions = [ 0.05  0.1 0.15    0.2 0.25    0.30    0.3500    0.4    0.45    0.5    0.55    0.6   0.8 0.85    .9 ];
                CO2_ymax = CO2_Concentration_in_Embayment + 10; % ymin = 0;
                CO2_textloc = CO2_ymax.*CO2_textpositions;
                          
                if iter==z  
                
                    
                subplot(2,2,1)
                mesh(xq,yq,flipud(H2O_C),'FaceColor','flat')
                view(2)
                xlim([min(min(xq)),max(max(xq))])
                ylim([min(min(yq)),max(max(yq))])
                title('Modeled')
                c = colorbar;
                c.Label.String = 'H_2O (wt%)';
                caxis([0, max(max(vqH2O))])
                daspect([1 1 1])
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
                axis off
      
                subplot(2,2,2)
                mesh(xq,yq,flipud(CO2_C),'FaceColor','flat')
                view(2)
                xlim([min(min(xq)),max(max(xq))])
                ylim([min(min(yq)),max(max(yq))])
                title('Modeled')
                c = colorbar;
                c.Label.String = 'CO_2 (ppm)';
                caxis([0, max(max(vqCO2))])
                daspect([1 1 1])
                axis off
            
                subplot(2,2,3)
                mesh(xq,yq,flipud(vqH2O),'FaceColor','flat')
                view(2)
                xlim([min(min(xq)),max(max(xq))])
                ylim([min(min(yq)),max(max(yq))])
                title('Instrument Data')
                c = colorbar;
                c.Label.String = 'H_2O (wt%)';
                caxis([0, max(max(vqH2O))])
                daspect([1 1 1])
                axis off
                
                subplot(2,2,4)
                mesh(xq,yq,flipud(vqCO2),'FaceColor','flat')
                view(2)
                xlim([min(min(xq)),max(max(xq))])
                ylim([min(min(yq)),max(max(yq))])
                title('Instrument Data')
                c = colorbar;
                c.Label.String = 'CO_2 (ppm)';
                caxis([0, max(max(vqCO2))])
                daspect([1 1 1])
                axis off
                sgtitle({append('Time Elapsed (hrs) = ',num2str(dt_total/(3600))),append('Decompression rate (MPa/s) = ','10^{', num2str(log10(1000*(Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e9)),'}')})

                drawnow
                             
                z=z+1000; %This controls how often the plot updates, higher numbers is less frames

               end 

             iter=iter+1;

        %%%update solubility
         %update pressure based on ascent
         GPa = ( (i_depth - (dt_total*Ascent_rate)) * ( F_MDensity * 9.81 ) ) / 1e9 ;
         % update pressure and determine dissolved vector for the given
         % pressure
         [ TGPa_d_i, TGPa_x_i ] = min( (abs(TGPa-GPa)) );
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
         XfCO2v_GPa =  moleCO2v_GPa ./ (moleH2Ov_GPa + moleCO2v_GPa);

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
         u_H2O_C_Cry = cry_mask .* H2O_C;
         H2O_C = u_H2O_C_Con + u_H2O_C_Emb + u_H2O_C_Cry ;        
         
         %CO2
         CO2d_p = CO2d_GPa(vector_index);
         CO2v_p = i_CO2t - CO2d_GPa(vector_index);
         %update diss CO2 in conduit based on above
         u_CO2_C_Con = conduit_mask .* CO2d_p;
         u_CO2_C_Emb = emb_mask .* CO2_C;
         u_CO2_C_Cry = cry_mask .* CO2_C;
         CO2_C = u_CO2_C_Con + u_CO2_C_Emb + u_CO2_C_Cry;
         
         %update time step for next round based on temperature
         emb_H2O_D = emb_mask .* H2O_D;
         minH2O_D = mean( mean(emb_H2O_D(emb_H2O_D>0)) );
         dt = 0.4*(dx^2/(minH2O_D));
         dt_total = dt_total+dt;
    end
toc
%%
   %%%draw a new transect to plot
    close()
    uiwait(msgbox({'Please draw polyline for 1D profile through embayment';'Right-click when done'},'Profile'));
    figure()
    if maskmatch == 0
        imshow(vqH2O, 'DisplayRange',[min(vqH2O(:)) max(vqH2O(:))], 'Colormap',parula);
        c = colorbar;
        c.Label.String = 'H_{2}O (wt.%)';
    else
        imshow(vqCO2, 'DisplayRange',[min(vqCO2(:)) max(vqCO2(:))], 'Colormap',parula);
        c = colorbar;
        c.Label.String = 'CO_{2} (ppm)';
    end    
    title({'Draw a polyline for 1D profile through embayment','Right-click when done'})
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.4 0.8]);
    set(gca, 'Position', [0.2 0.2 0.7 0.7]);
    profileline = drawpolyline('Color','green');
    profile_mask = int8(createMask(profileline));
    [profile_idy,profile_idx] = find(profile_mask.');
    close()
    
    %determine coordinates of transect and extract data
    H2O_C_Pro = zeros(1,length(profile_idy));
    vqH2O_Pro = zeros(1,length(profile_idy));
    CO2_C_Pro = zeros(1,length(profile_idy));
    vqCO2_Pro = zeros(1,length(profile_idy));
    for n=1:length(profile_idy)
        H2O_C_Pro(1,n) = H2O_C(profile_idx(n),profile_idy(n));
        vqH2O_Pro(1,n) = vqH2O(profile_idx(n),profile_idy(n));
        CO2_C_Pro(1,n) = CO2_C(profile_idx(n),profile_idy(n));
        vqCO2_Pro(1,n) = vqCO2(profile_idx(n),profile_idy(n));
    end
    
    %plot all figures
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1 0.02 0.4 0.98]);
    
    subplot(3,2,1)
    mesh(xq,yq,flipud(H2O_C),'FaceColor','flat')
    view(2)
    xlim([min(min(xq)),max(max(xq))])
    ylim([min(min(yq)),max(max(yq))])
    title('Modeled')
    c = colorbar;
    c.Label.String = 'H_2O (wt%)';
    caxis([0, max(max(vqH2O))])
    daspect([1 1 1])
    hold on
    line(min(min(xq))+(profile_idy*dx),max(max(yq))-(profile_idx*dx),  max(max(H2O_C))*ones(1,length(profile_idx)),'Color','black' );
    hold off;
    axis off
    % Enlarge figure to full screen.
    
    subplot(3,2,2)
    mesh(xq,yq,flipud(CO2_C),'FaceColor','flat')
    view(2)
    xlim([min(min(xq)),max(max(xq))])
    ylim([min(min(yq)),max(max(yq))])
    title('Modeled')
    c = colorbar;
    c.Label.String = 'CO_2 (ppm)';
    caxis([0, max(max(vqCO2))])
    daspect([1 1 1])
    axis off
    
    subplot(3,2,3)
    mesh(xq,yq,flipud(vqH2O),'FaceColor','flat')
    view(2)
    xlim([min(min(xq)),max(max(xq))])
    ylim([min(min(yq)),max(max(yq))])
    title('Instrument Data')
    c = colorbar;
    c.Label.String = 'H_2O (wt%)';
    caxis([0, max(max(vqH2O))])
    daspect([1 1 1])
    axis off
    
    subplot(3,2,4)
    mesh(xq,yq,flipud(vqCO2),'FaceColor','flat')
    view(2)
    xlim([min(min(xq)),max(max(xq))])
    ylim([min(min(yq)),max(max(yq))])
    title('Instrument Data')
    c = colorbar;
    c.Label.String = 'CO_2 (ppm)';
    caxis([0, max(max(vqCO2))])
    daspect([1 1 1])
    axis off
    
    subplot(3,2,5)
    errorbar((profile_idx-profile_idx(1))*dx, flip(vqH2O_Pro),flip(vqH2O_Pro.*(H2Oerror./100.0)), 'ko','MarkerSize',4, 'MarkerFaceColor','k','CapSize',0,'color', '[.7 .7 .7]')
    xlabel('Position (\mum)')
    ylabel('H_2O (wt%)')
    axis([0 length(profile_idx)*dx 0 max(max(vqH2O_Pro))*1.05])
    title('1D Profile')
    hold on
    plot((profile_idx-profile_idx(1))*dx,flip(H2O_C_Pro),'color','[0 0.4470 0.7410]','linewidth',2)
    hold off
    
    subplot(3,2,6)
    errorbar((profile_idx-profile_idx(1))*dx, flip(vqCO2_Pro),flip(vqCO2_Pro.*(CO2error./100.0)), 'ko','MarkerSize',4, 'MarkerFaceColor','k','CapSize',0,'color', '[.7 .7 .7]')
    xlabel('Position (\mum)')
    ylabel('CO2 (ppm)')
    axis([0 length(profile_idx)*dx 0 max(max(vqCO2_Pro))*1.05])
    title('1D Profile')
    hold on
    plot((profile_idx-profile_idx(1))*dx,flip(CO2_C_Pro),'color','[0 0.4470 0.7410]','linewidth',2)
    hold off
    leg1 = legend({'Instrument Data','Model'},'Location','southeast');
    leg1.FontSize = 10;

    sgt = sgtitle({append('Initial H_2O (wt%) = ',num2str(H2O_Concentration_in_Embayment), ...
        '    Initial CO_2 (ppm) = ',num2str(CO2_Concentration_in_Embayment),...
        '    Inital Exsolved Phase (wt%) = ',num2str(init_exs)),...
        append('Temperature (^{o}C) = ',num2str(TempC), ...
        '    Time Elapsed (hrs) = ',num2str(Ascent_Duration), ...
        '    Decompression rate (MPa/s) = ','10^{', num2str(log10(1000*(Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e9))),'}',...
        append('Starting Pressure (MPa) = ',num2str(1000*GPa_i),...
        '    Quench Pressure (MPa) = ',num2str(1000*GPa_quench),...
        '    Quench Depth (m) = ',num2str(quench_depth)), ''}) ;   
    sgt.FontSize = 10;
    drawnow
    savefig(append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_2d_profile.fig'))
    saveas(gcf,append(outputname,'_M0_',strrep(num2str(init_exs),'.','_'),'_2d_profile.jpg'))   

    if saveout == 1
    %write results out 
    init_exs_string = strrep(num2str(init_exs),'.','_');
    header = {'Position (microns)','H2O Instrument Data (wt%)','CO2 Instrument Data (ppm)','H2O Model (wt%)','CO2 Model (ppm)'};
    oned_data = transpose([transpose((profile_idx-profile_idx(1)).*dx); flip(vqH2O_Pro); flip(vqCO2_Pro); flip(H2O_C_Pro); flip(CO2_C_Pro)]);
    oned_header_data = [header; num2cell(oned_data)];
    writecell(oned_header_data,append(outputname,'_M0_',init_exs_string,'_1d_profile.xls'))

    header1 = {'Inital H2O (wt%)','Inital CO2 (ppm)','Inital Exsolved Phase [M0] (wt%)','Temperature (Celsius)',...
        'Time Elapsed (hrs)','Decompression rate (MPa/s)', 'Starting Pressure (MPa)','Quench depth (m)','Quench Pressure (MPa)'};
    output_data = [H2O_Concentration_in_Embayment CO2_Concentration_in_Embayment init_exs TempC...
        Ascent_Duration ((Ascent_rate * ( F_MDensity * 9.81 ) ) / 1e6) GPa_i*1000 quench_depth GPa_quench*1000 ];
    output_header_data = [header1; num2cell(output_data)];
    writecell(output_header_data,append(outputname,'_M0_',init_exs_string,'_parameters.xls'))
    end

%%
%%%References used in this code
    %Zhang, Y., & Ni, H. (2010). Diffusion of H, C, and O components in silicate melts. Reviews in Mineralogy and Geochemistry, 72(1), 171-225.
    %Ni, H., & Zhang, Y. (2008). H2O diffusion models in rhyolitic melt with new high pressure data. Chemical Geology, 250(1-4), 68-78.
    %Zhang, Y., Xu, Z., Zhu, M., & Wang, H. (2007). Silicate melt properties and volcanic eruptions. Reviews of Geophysics, 45(4).
    %Liu, Y., A. T. Anderson, and C. J. N. Wilson (2007). Melt pockets in phenocrysts and decompression rates of silicic magmas before 
    %       fragmentation, J. Geophys. Res., 112, B06204.
           