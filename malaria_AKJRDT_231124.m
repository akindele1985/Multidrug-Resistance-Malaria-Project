% MIT Software licence
% Copyright 2023 Onifade, Rychtar and Taylor
%
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the “Software”),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom
% the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
% IN THE SOFTWARE.
%
%

function main
% This is MATLAB code for the malaria project
% modified  by Akindele Onifade, Jan Rychtar, Dewey Taylor
% based on the template by Jan Rychtar and  Dewey Taylor
% Email questions to Jan Rychtar rychtarj@vcu.edu

% Version from September 22, 2023
% Written on Matlab 2022a (version 9.12.0.2009381)


%% Technical preliminaries

% Close all figures
close all;

% Clear all variables
clear all;

% Set the random number generator for reproducibility
rng(0)

% Check if the folder for figures exists and create if not
if not(isfolder('Figures'))
    mkdir('Figures')
end

% Set parameters and other variables as global
global r bhm bmh LambdaH LambdaM muH muM beta alphaHinv alphaH rhoH
global sigmaH1inv sigmaH2inv sigmaH1 sigmaH2 deltaH phiH omegaH alphaM alphaMinv x m
% Make the parameter symbolic
syms r bhm bmh LambdaH LambdaM muH muM beta alphaHinv alphaH rhoH
syms sigmaH1inv sigmaH2inv sigmaH1 sigmaH2 deltaH phiH omegaH alphaM alphaMinv x m

% Set auxiliary variables as global and symbolic
global v_EH v_IH v_TH1 v_RH v_TH2 v_DH2 v_EM
syms v_EH v_IH v_TH1 v_RH v_TH2 v_DH2 v_EM



% Set the parameters.
% First column is the name used in matlab,
% second column is the same but as a string
% third is the LaTeX label used in figures
% fourth is the based value
% fifth is the minimal value
% sixth is the max value
Params = {...
    LambdaH	,	'LambdaH'	,	'$\Lambda_h$'	,	22150	,	15000	,	25000	;
    LambdaM	,	'LambdaM'	,	'$\Lambda_m$'	,	1.6*10^10	,	10^10	,	2*10^10	;
    muH	,	'muH'	,	'$\mu_h$'	,	3.56*10^(-5)	,	2*10^(-5)	,	5*10^(-5)	;
    muM	,	'muM'	,	'$\mu_m$'	,	0.05	,	0.04	,	0.1	;
    beta	,	'beta'	,	'$\beta$'	,	0.5	,	0.03	,	1	;
    alphaHinv	,	'alphaHinv'	,	'$\alpha_h^{-1}$'	,	10.5	,	5	,	150	;
    rhoH	,	'rhoH'	,	'$\rho_h$'	,	 0.8	,	0.5	,	1	;
    sigmaH1inv	,	'sigmaH1inv'	,	'$\sigma_{h1}^{-1}$'	,	3	,	1	,	5	;
    sigmaH2inv	,	'sigmaH2inv'	,	'$\sigma_{h2}^{-1}$'	,	7	,	5	,	10	;
    deltaH	,	'deltaH'	,	'$\delta_h$'	,	0.05	,	0.01	,	0.1	;
    phiH ,   'phiH'  , '$\varphi_h$',    0.0055, 0, 0.001 ;
    omegaH,  'omegaH' , '$\omega_h$' , 0.05, 0 , 0.01;
    r, 'r', '$r$', 0.8, 0, 1;
    bmh , 'bmh', '$b_{mh}$', 0.09, 0, 1;
    bhm , 'bhm', '$b_{hm}$', 0.09, 0, 1;
    alphaMinv	,	'alphaMinv'	,	'$\alpha_m^{-1}$'	,	10	,	9	,	20	;
    x, 'x', '$x$, testing prevalence', 0.0, 0, 1;
    };



Comps = {'Sh';
        'Eh1';
        'Ih1';
        'Th1';
        'Rh1';
        'Eh2';
        'Ih2';
        'Th2';
        'Dh2';
        'Rh2';
        'Sm';
        'Em1';
        'Im1';
        'Em2';
        'Im2';};

global Sh Eh1 Ih1 Th1 Rh1 Eh2 Ih2 Th2 Dh2 Rh2 Sm Em1 Im1 Em2 Im2;
syms Sh Eh1 Ih1 Th1 Rh1 Eh2 Ih2 Th2 Dh2 Rh2 Sm Em1 Im1 Em2 Im2;
unknowns = [Sh Ih1 Ih2 Sm Im1 Im2];


% Get the simpler arrays of parameters (individual columns of the Params
paramArray = Params(:,1); %Matlab names
paramNames = Params(:,2); % matlab names as strings
paramLabels = Params(:,3); % Latex labels
paramVals = Params(:,4);  % actual values

% Calculate auxiliary variables
alphaH = 1/alphaHinv;
sigmaH1 = 1/sigmaH1inv;
sigmaH2 = 1/sigmaH2inv;
alphaM = 1/alphaMinv;
m = LambdaM/muM /(LambdaH/muH);


% Calculate variables for the reproduction numbers
v_EH =  muH + alphaH;
v_IH =  rhoH + muH+ deltaH;
v_TH1 =  sigmaH1 + muH + deltaH + phiH;
v_TH2 =  sigmaH1 + muH + deltaH;
v_DH2 =  sigmaH2 + muH + deltaH;
v_EM =  alphaM +  muM;



% Calculate R22
R22 = sqrt(1/muM * bmh*bhm*beta^2*m * alphaH/(alphaH + muH) * alphaM/(alphaM + muM) * ...
    ( 1/v_IH + (1-x)*rhoH/v_IH *1/v_TH2  + (x*rhoH/v_IH + (1-x)*rhoH/v_IH * sigmaH1/v_TH2) * 1/v_DH2 * r));

% Evaluate R22 numerically
R22num = double(subs(R22, paramArray, paramVals ));

R0 = R22;

% Plot dependnce on all parameters
PlotDependence(R0, '$R_0$', 'R0')

% Get R0 at the baseline (x=0)
evaluate(R22)



% THIS PLOTS THE dependence of the equilibria on the x
plotDependence_on_x



return


F =beta*[
    %E1 I1  T1  E2  I2  T2  D2  Em1 Im1  Em2 Im2
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , bmh , 0 , 0 ; % Eh1
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Ih1
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Th1
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , bmh ; % Eh2
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Ih2
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Th2
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Dh2
    0 , bhm*m , r*bhm*m , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; %Em1
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; % Ih2
    0 , 0 , 0 , 0 , bhm*m , bhm*m , r*bhm*m , 0 , 0 , 0 , 0; %EH2
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ; %IH2
    ];

V = [...
    -v_EH , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0;
    alphaH , -v_IH , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
    0 , rhoH , -v_TH1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
    0 , 0 , 0 , -v_EH, 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
    0 , 0 , 0 , alphaH , -v_IH , 0 , 0 , 0 , 0 , 0 , 0 ;
    0 , 0 , 0 , 0 , (1-x)*rhoH , -v_TH2 , 0 , 0 , 0 , 0 , 0 ;
    0 , 0 , phiH , 0 , x*rhoH , sigmaH1 , -v_DH2 , 0 , 0 , 0 , 0 ;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , -v_EM , 0 , 0 , 0 ;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , alphaM , -muM , 0 , 0 ;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -v_EM , 0  ;
    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , alphaM , -muM  ;
    ];


Vinv = inv(V);
%
% for auxcol = 1:size(V,1)
% for auxrow = 1:size(V,1)
%     disp([ string(Vinv(auxcol,auxrow)) '&'])
% end
% disp('\\')
% end
%
%
% return


A = F*inv(V);
%eig(A)

% Calculate values of A after substitution for the parameters
Anum= subs(A, paramArray, paramVals );

% Calculate the eigenvalues (and make them double to see them as decimals)
double(eig(Anum))

R11 = sqrt(1/muM * bmh * bhm * beta^2 * m * alphaH/(alphaH + muH) * alphaM/(alphaM + muM) * ...
    ( 1/v_IH  + rhoH/v_IH * 1/v_TH1* r ));

% Evaluate R11 numerically at parameter values
R11num = double(subs(R11, paramArray, paramVals ));





% % Do the sensitivity
% if R22num>R11num % in this case, R0 = R22
%     sensitivity(R22, '$R_0$', 'R0')
% else % in this case R0 = R11
%     disp('R11 should not be bigger than R22')
%     %disp('R11 is larger')
%     %sensitivity(R11, '$R_0$', 'R0')
% end
% 
% 
% % Do sensitivity based on Marino
% figure('Units','Inches','Position',[0, 0, 6.3, 2.5]);
% t=tiledlayout(1,2)
% t.Padding = 'compact';
% t.TileSpacing = 'compact';
% mytile = 1;
% 
% LHS_PRCC(@getR0, '$R_0$', 'R0', Params, 1000, mytile, mytile+1)




return


%%%%%%%%%

    function output = getR0
        % Returns R0 which in this case is just sqrt of R22

        tempR22 = sqrt(1/muM * bmh*bhm*beta^2*m * alphaH/(alphaH + muH) * alphaM/(alphaM + muM) * ...
            ( 1/v_IH + (1-x)*rhoH/v_IH *1/v_TH2  + (x*rhoH/v_IH + (1-x)*rhoH/v_IH * sigmaH1/v_TH2) * 1/v_DH2 * r));
        output = double(subs(tempR22, paramArray, paramVals ));
    end

    function output = sensitivity(y, y_label, y_name)
        % Performs sensitivity w.r.t. Arriola  and Hyman
        % Takes the respons variable y and gets the index for all params
        % specified in par = paramArray
        % y_lable is how the parameter is called

        % Set what the parameters are
        par = paramArray;

        % For every parameter get the sensitivity index
        for auxpar = 1: length(par)
            % first do it analytically
            sens_index = (par(auxpar)/y) * diff(y,par(auxpar));
            % Then evaluate the symbolic expression at the param values
            num_sens_index(auxpar) = evaluate(sens_index);
        end

        % Order the sensitivity indices
        [ordered_sens_index,idx] = sort(num_sens_index, 'descend');

        % Print out the indices for the table used in LaTeX
        for auxpar = 1: length(par) %for every parameter
            disp([paramLabels{idx(auxpar)} ' & ' num2str(ordered_sens_index(auxpar)) ' \\'])
        end


        % Do the bar graph
        % set the ordered parameters
        temp_parameters = paramLabels(idx);
        % find the indices of significant parameters
        idx_significant = abs(ordered_sens_index)>0.05;
        xPams = categorical(temp_parameters(idx_significant));
        xPams = reordercats(xPams,temp_parameters(idx_significant) );
        figure
        barh(xPams, ordered_sens_index(idx_significant), 'facecolor', 'flat');
        FormatFigure
        SaveFigure([y_name 'bar_plot'])
    end


    function output = PlotDependence(y, y_label, y_name)
        % Takes the response variable y and plots how it depends on the parameters
        % specified in par = paramArray
        % y_lable is how the response variable is called
        % y_name is how it will be saved in the figures

        % Set what the parameters are
        par = paramArray;

        % Make the plots
        for auxpar = 1: length(par) % for every parameter
            %
            % do the plotting only if sensitivity index is large enough
            if 1>0 %abs(num_sens_index(auxpar))>-0.05
                % Get the params other that the focal one
                other_params_idx = ~strcmp(paramNames,paramNames(auxpar));

                % Substitute every other parameter into y
                y_on_par = subs(R0, paramArray(other_params_idx), paramVals(other_params_idx));
                % Plot the figure
                figure
                hold on
                % draw the curve of the variable, get the limits from the
                % params array
                fplot(y_on_par, [Params{auxpar,5} Params{auxpar,6}], 'linewidth',2, 'color',  'k');
                % Draw the horizontal line at y=1
                yline(1, ['k--'], 'linewidth', 1);
                % Specify the x and y axis lables
                xlabel(paramLabels(auxpar))
                ylabel(y_label)
                % Set the y lim to 0 min and let matlab choose the max
                % However, for x, set it to 10
                if strcmp(Params{auxpar,2},'x')
                    ylim([0, 10])
                else
                    ylim([0, inf])
                end

                % Format the figure
                FormatFigure
                % Save the figure
                filename = strjoin({y_name, paramNames{auxpar}}, '_');
                SaveFigure(filename)
            end

        end

    end

    function LHS_PRCC(RespFunction, RespName, FileName, Pars, Samples, tile_for_hist, tile_for_bars)
        % Does PRCC uncertainty and sensitivity analysis
        % Based on Marino, S., Hogue, I. B., Ray, C. J., and Kirschner, D. E. (2008).
        % A methodology for performing global uncertainty and sensitivity analysis in systems
        % biology.Journal of Theoretical Biology, 254(1):178-196
        % Part of the code for PRCC taken from
        % http://malthus.micro.med.umich.edu/lab/usanalysis.html

        % Set the parameters if some are not provided
        if nargin < 5 % if Samples not provided
            % Specify how many times we want to generate the parameter
            % values
            Samples = 1000;
        end

        if nargin < 4  % if the Pars is not specified
            % Do sensititivity on all parameters
            Pars = Params;
        end

        disp('Doing sensitivity analysis based on Marino et al')

        wiggle_room = 0.1; %(how much can the parameter vary)

        % Get the number of parameters, discard the number of columns
        [numPars,~] = size(Pars); %

        % Get the parameter names
        myNames = Pars(:,2);

        % Initiate
        paramsUn = zeros(Samples, numPars); % This stores the values of the parameters
        allResp_Values = zeros(1,Samples); % This stores the values of the response function
        % Generate Latin Hypercube Sampling
        for aux1 = 1:Samples

            for aux = 1:numPars
                % Go through all the parameters
                % Evaluate the command that assigns the parameter with a given distribution with a parameter range
                % We need a column vector, but we will generate a row vector
                % first and transpose it during the assignment so that we can
                % use the eval command

                % Assign parameters at random but within bounds
                eval([myNames{aux} '=' num2str((Pars{aux,4}*2*wiggle_room)*rand+Pars{aux,4}*(1-wiggle_room)) '; '])

                % Evaluate the command that assigns the parameter with a given name the proper base value
                eval(['paramsUn(aux1, aux) = ' myNames{aux} ';'])

                % Calculate auxiliary variables
                %CalculateAuxiliaryVariables
            end
            %temp = paramsUn(aux1, :)
            %pause
            % Evaluate the response variable
            allResp_Values(aux1) = RespFunction();
            %pause
        end

        allResp_Values_truncated = allResp_Values(allResp_Values<20);

        % Plot the histogram to show uncertainty
        %figure
        nexttile(tile_for_hist)
        %allResp_Values
        histogram(allResp_Values_truncated,  'Normalization', 'probability');
        xlabel(RespName, 'Interpreter', 'latex')
        ylabel('Frequency')

        nexttile(tile_for_bars)
        %FormatFigure
        %SaveFigure([FileName '_hist_Marino'])

        % Get tha average value, do not consided NaN (not a number) values
        disp(['The average value of ' RespName ' is ' num2str(nanmean(allResp_Values))]);

        % Plot the bars with PRCC indices
        PRCClocal(paramsUn, allResp_Values', FileName, RespName, Pars)
    end

    function prcc=PRCClocal(values, Y, FileName, YName, Pars)
        % PRCC procedure
        % input:  values: values of the parameters
        %       and Y, response function that we investigate how it depends
        % on the parameters
        % output - the prcc coefficients

        % Based on Marino, S., Hogue, I. B., Ray, C. J., and Kirschner, D. E. (2008).
        % A methodology for performing global uncertainty and sensitivity analysis in systems
        % biology. Journal of Theoretical Biology, 254(1):178-196
        % Part of the code for PRCC taken from
        % http://malthus.micro.med.umich.edu/lab/usanalysis.html

        % Set what we consider important
        importance_threshold = 0.01;

        % Initiate the LHS matrix that will carry all the parameters
        LHSmatrix= values;

        if nargin == 3 % if the last parameter is not specified
            % Do sensititivity on all parameters
            Pars = Params;
        end

        % Get the number of parameters, discard the number of columns
        [numPars,~] = size(Pars); %

        % Get the parameter labels
        myLabels = Pars(:,3);

        % Preallocate PRCC values
        prcc = zeros(1,numPars);

        for i=1:numPars  % Loop for the parameters

            % Remove the ith parameter
            LHStemp=LHSmatrix;
            LHStemp(:,i)=[];
            Ztemp=LHStemp;
            %LHStemp=[];

            % Get the correlation coefficient
            [rho_aux,~]=partialcorr([LHSmatrix(:,i),Y],Ztemp,'type','Spearman');

            % Record the correlation coefficient
            prcc(1,i)=rho_aux(1,2);

        end

        % Sort the correlation coefficients for better plotting
        [sortedPRCC,I] = sort(prcc, 'descend');

        % PLot only the important coefficients
        %figure
        % Determine which ones are important
        importantPRCC = sortedPRCC(abs(sortedPRCC)>importance_threshold);
        if isempty(importantPRCC)
            importantPRCC = sortedPRCC;
        end

        h=barh(importantPRCC);
        ylim([0.5 length(importantPRCC)+0.5]);
        % Ticks at integers and labeled by the parameter names/labels
        ax = gca;
        ax.YTick = 1:length(importantPRCC);
        ax.YTickLabel = myLabels(I(abs(sortedPRCC)>importance_threshold));
        ax.TickLabelInterpreter = 'latex';
        ylabel('Parameters', 'Interpreter', 'latex');
        % Set the x axis between -1 and 1
        xlim([-1 1]);
        xlabel(['PRCC for ' YName],  'Interpreter', 'latex');
        % Make the color gray
        h.FaceColor = [0.5 0.5 0.5];

        %FormatFigure;
        SaveFigure([FileName 'Marino_sens'])
    end

    function output = check_solutions(candidate_solutions)
        % Checks how far are the potential solutions from the solutions
        % We get the values of the right hand sides of the ODEs and see how
        % far are they from 0
        output = max(abs(getRHSs(candidate_solutions)));
    end

    function output = getEquations
        % Returns the equations to solve

        Eh1 = evaluate((rhoH + muH + deltaH)*Ih1/alphaH);
        Th1 = evaluate(rhoH* Ih1/(sigmaH1 +  muH+ deltaH+phiH));
        Rh1 = evaluate(sigmaH1* Th1/(omegaH+muH));
        Eh2 = evaluate((rhoH + muH + deltaH)*Ih2/alphaH);
        Th2 = evaluate((1-x)*rhoH *Ih2/(sigmaH1 + muH+ deltaH));
        Dh2 = evaluate((phiH* Th1 + x*rhoH *Ih2 + sigmaH1* Th2)/ (sigmaH2 +  muH+ deltaH));
        Rh2 = evaluate(sigmaH2* Dh2/ (omegaH+muH));

        Em1 = evaluate(muM* Im1/alphaM);
        Em2 = evaluate(muM* Im2/alphaM);

        Nh = Sh + Eh1 + Ih1 + Th1 + Rh1 + Eh2 + Ih2 + Th2 + Dh2 + Rh2; %sum(comps(1:10)); % all humans
        %Nm = sum(comps(11:15)); % all mosquitoes

        % Set auxiliary variables
        lambdah1 = bmh* beta *Im1;
        lambdah2 = bmh* beta *Im2;
        lambdam1 = bhm*beta*(Ih1+ r* Th1);
        lambdam2 =  bhm*beta*(Ih2+Th2+r* Dh2);

        % Set all the equations
        eq1 = evaluate(0 == (LambdaH + omegaH*(Rh1+Rh2) - muH*Sh)*Nh - (lambdah1+lambdah2)*Sh );
        eq2 = evaluate(0 == lambdah1*Sh - (alphaH + muH)*Eh1*Nh);

        eq6 = evaluate(0 == lambdah2* Sh - (alphaH + muH)*Eh2*Nh);

        eq11 = evaluate(0 == (LambdaM  - muM* Sm)*Nh - (lambdam1+lambdam2)*Sm);
        eq12 = evaluate(0 == lambdam1*Sm - (alphaM +muM)* Em1*Nh);
        eq14 =  evaluate(0 == lambdam2* Sm - (alphaM +muM)* Em2*Nh);

        % Set the output
        output = [eq1, eq2, eq6, eq11, eq12,  eq14];

    end

    function output = numSolve
        % Numerically solves the system of equilibria equations
        % Returns a structured array where the solutions are stored


   

        % Get the system of equations
        eqns = getEquations;

%        range = [0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; ...
%                 0 Inf];

        % Numerically solve the equations and store it in "struct" variable
%         n = 0;
%         positive = false;
%         while ~positive
%            n = n+1  % increase the counter
            %solve the equation

            sols = solve(eqns, unknowns);
            mySh = double(sols.Sh);
            myIh1 = double(sols.Ih1);
            myIh2 = double(sols.Ih2);
            mySm = double(sols.Sm);
            myIm1 = double(sols.Im1);
            myIm2 = double(sols.Im2);


            posSh = (mySh>0)' & (mySh~=evaluate(LambdaH/muH))'; %positive but not DFE
            posIh1 = (myIh1>=0)';
            posIh2 = (myIh2>=0)';
            posSm = (mySm>0)';
            posIm1 = (myIm1>=0)';
            posIm2 = (myIm2>=0)';

            pos = posSh & (posIh1 & posIh2)& posSm &posIm1 & posIm2;  

%             myIh2(pos)
%             
%             %sols = [sh, ih1, ih2, sm, im1, im2];
%             %if all(sols>0)
%                 positive = true;
%            % end
%        end

        %         % Interpret the solutions
        %         numSols = size(sols.(Comps{1}),1); % This is a fix of above suggested by matlab with dynamic fields
        %
        %         disp_Diary(['There are ' num2str(numSols) ' sets of solutions']);
        %
        %         numericalSolutions = zeros(numSols, numComps); % Just initialize the array storing all solutions
        %         % numericalSolutions(i,j) will be the ith set, compartment j
        %
        %         for solSet = 1:numSols
        %             disp_Diary(['In the order ' Comps{:} ', the solution set number ' num2str(solSet) ' is'])
        %
        %             % Actually get the solutions from the structured array
        %             for comp = 1:numComps
        %                 numericalSolutions(solSet, comp) = sols.(Comps{comp})(solSet); % This is a fix suggested by matlab with dynamic fields
        %             end
        %             % Print out the solutions
        %             disp_Diary(num2str(numericalSolutions(solSet,:) ));
        %         end


        output = [mySh(pos) myIh1(pos) myIh2(pos) mySm(pos) myIm1(pos) myIm2(pos)];
    end

    function plotDependence_on_x
        % Store original value of x before we start changing it
        original_x = Params{strcmp(Params(:,2),'x'),4};

        % Set all possible values of x.
        % The dependence is more or less linear, so we do not need many
        % values
        all_x = linspace(0,1,6);
        for auxx = 1:length(all_x)

            % set the value of x to the params array
            Params{strcmp(Params(:,2),'x'),4} = all_x(auxx);
            % Solve the system
            solutions =  numSolve;

            % Unpack solutions
            %ih1(auxx) = solutions(2); % these will actually turn out 0 so we can ignore them
            ih2(auxx) = solutions(3);

        end
        % Restore original value of x
        Params{strcmp(Params(:,2),'x'),4} = original_x ;

        % Calculate th2 and dh2 from ih2
        th2 = evaluate((1-all_x)*rhoH .*ih2/(sigmaH1 + muH+ deltaH));
        dh2 = evaluate(( all_x*rhoH .*ih2 + sigmaH1* th2)/ (sigmaH2 +  muH+ deltaH));

        % Plot figures
        figure
        plot(all_x,ih2, 'k:')
        hold on
        plot(all_x, ih2+th2+dh2, 'r-') % this one gives all 0
        plot(all_x, th2, 'b--') % this one gives all 0
        
        xlabel(Params{strcmp(Params(:,2),'x'),3});
        ylabel('Number of people')

        FormatFigure
        SaveFigure('Equilibrium_on_x')
  
        % Output numerical values for the paper

        disp(['When $x=0$, there are ' num2str(ih2(1)) ' $I_{h2}$,' ...
                num2str(th2(1)) '$T_{h2}$, and ' num2str(dh2(1)) '$D_{h2}$.' ...
               ' When $x=1$, there are ' num2str(ih2(end)) ' $I_{h2}$,' ...
                num2str(th2(end)) '$T_{h2}$, and ' num2str(dh2(end)) '$D_{h2}$'])

    end


    function output = getRHSs(comps)
        % This function takes the compartments
        % input = [Sh, Eh1, Ih1, Th1, Rh1

        % output = returns the right hand sides of the ODEs equations

        % unpack the input
        Sh = comps(1);
        Eh1 = comps(2);
        Ih1 = comps(3);
        Th1 = comps(4);
        Rh1 = comps(5);
        Eh2 = comps(6);
        Ih2 = comps(7);
        Th2 = comps(8);
        Dh2 = comps(9);
        Rh2 = comps(10);
        Sm = comps(11);
        Em1 = comps(12);
        Im1 = comps(13);
        Em2 = comps(14);
        Im2 = comps(15);

        Nh = Sh + Eh1 + Ih1 + Th1 + Rh1 + Eh2 + Ih2 + Th2 + Dh2 + Rh2; %sum(comps(1:10)); % all humans
        %Nm = sum(comps(11:15)); % all mosquitoes

        % Set auxiliary variables
        lambdah1 = bmh* beta *Im1/Nh;
        lambdah2 = bmh* beta *Im2/Nh;
        lambdam1 = bhm*beta*(Ih1+ r* Th1)/Nh;
        lambdam2 =  bhm*beta*(Ih2+Th2+r* Dh2)/Nh;

        % Set all the equations
        eq1 = LambdaH + omegaH*(Rh1+Rh2) - (lambdah1+lambdah2)*Sh - muH*Sh;
        eq2 = lambdah1*Sh - (alphaH + muH)*Eh1;
        eq3 = alphaH* Eh1 - (rhoH + muH + deltaH)*Ih1;
        eq4 = rhoH* Ih1 - (sigmaH1 +  muH+ deltaH+phiH)*Th1;
        eq5 = sigmaH1* Th1   - (omegaH+muH)*Rh1;
        eq6 = lambdah2* Sh - (alphaH + muH)*Eh2;
        eq7 = alphaH * Eh2 - (rhoH + muH + deltaH)*Ih2;
        eq8 =  (1-x)*rhoH *Ih2 - (sigmaH1 + muH+ deltaH)*Th2;
        eq9 =  phiH* Th1 + x*rhoH *Ih2 + sigmaH1* Th2 - (sigmaH2 +  muH+ deltaH)*Dh2;
        eq10 = sigmaH2* Dh2   - (omegaH+muH)*Rh2;
        eq11 = LambdaM - (lambdam1+lambdam2)*Sm - muM* Sm;
        eq12 = lambdam1*Sm - (alphaM +muM)* Em1;
        eq13 = alphaM * Em1 - muM* Im1;
        eq14 =  lambdam2* Sm - (alphaM +muM)* Em2;
        eq15 = alphaM * Em2 - muM * Im2;

        % Set the output
        output = [eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15];

    end

    function output = ODEmodel(t,y)
        % This function provides the ODE equations for the ODE solver
        % It probably does not need t as an input, but we followed matlab
        % example fully
        % the ODEs is dy/dt = right hand sides
        % where y is the whole vector (in our case p, Sh, Eh, Ah, Ch, Mr,
        % Sr, Er, Ir
        % and t is the time

        % Since we already have the right hand sides of the equations
        % written out for checking the equilibrium, this function just
        % calls for that
        output = getRHSs(y)'; % the ode solver requires a column vector
    end

    function plotODE(file_name, tspan, init_conditions)
        % Plots the solutions of the odes and saves the figure to the file
        % file_name (no extension) is a name of the file where the figure will be saved
        % tspan is the interval we will be solving the equations on
        % init_conditions is the initial conditions


        % Numerically solve the ODE
        %[t,y] = ode45(@(t,y) ODEmodel(t,y), tspan, init_conditions);
        [t,y] = ode15s(@(t,y) ODEmodel(t,y), tspan, init_conditions);


        % Unpack the solutions
        %[p, Sh, Eh, Ah, Ch, Mr, Sr, Er, Ir]
        psol = y(:,1);
        Shsol = y(:,2);  % S is the first coordinate of y
        Ehsol = y(:,3);
        Ahsol = y(:,4);
        Chsol = y(:,5);
        Mrsol = y(:,6);
        Srsol = y(:,7);  % S is the first coordinate of y
        Ersol = y(:,8);
        Irsol = y(:,9);

        Nrsol = Srsol + Ersol + Irsol;

        % Change the time to years for plotting
        t = t/365;

        % Plot the solutions, display it as prevalences (i.e., relative to
        % the population size N)
        figure('Units','Inches','Position',[0, 0, 9.4, 2.5]);
        tt=tiledlayout(1,3); % one row, 3 subfigures
        tt.Padding = 'compact'; % as closely packet together as possible
        tt.TileSpacing = 'compact'; % as closely packet together as possible
        % Start the first subfigure
        h=nexttile(tt);
        % Plot multiple graphs into the subfigure
        hold on
        % Set the title
        title('Humans', 'FontSize', 12)
        % plot the graphs
        plot(t,Chsol, 'r-.', 'Linewidth', 2)
        plot(t,Ahsol, 'b--', 'Linewidth', 2)
        plot(t,Shsol, 'k:', 'Linewidth', 2)
        plot(t,Chsol+Ahsol, 'g', 'Linewidth', 2)
        % Add the box around the figure
        box on
        % Set labels
        xlabel('Time (in years)', 'Interpreter','latex', 'FontSize', 12)
        ylabel('Prevalence', 'Interpreter','latex', 'FontSize', 12)

        % New subfigure
        nexttile;
        hold on
        title('Infected kissing bugs', 'FontSize', 12)
        plot(t,Irsol./Nrsol, 'r-', 'Linewidth', 2)
        % Set labels, show legend, format and save
        xlabel('Time (in years)', 'Interpreter','latex', 'FontSize', 12)
        ylabel('Prevalence', 'Interpreter','latex', 'FontSize', 12)
        box on

        % New subfigures
        nexttile
        hold on
        title('Protection level', 'FontSize', 12)
        plot(t,psol, 'k:', 'Linewidth', 2)
        ylim([0,1])
        box on
        xlabel('Time (in years)' , 'Interpreter','latex', 'FontSize', 12)
        ylabel('Level', 'Interpreter','latex', 'FontSize', 12)

        % Save the figure as an eps
        saveas(gcf, file_name, 'eps')

    end

%% Functions independent of the model

    function output = evaluate(X)
        % Takes a symbolic expression X and returns its value when the
        % parameters are plugged in
        output = eval(subs(X,paramArray,paramVals));
    end


    function FormatFigure(h)
        % Formats figures so that all figures are uniformly formatted

        if nargin == 0 % if the function is called without an argument, assume we want to format current figure
            h = gcf;
        end

        % Get axes handle
        haxes=get(h, 'CurrentAxes');

        % Get figure number
        fnbr = get(h,'Number');
        [fig_x,fig_y] = getPosition(fnbr);

        % Set dimensions of the figure
        set(h,'Units','Inches','Position',[fig_x, fig_y, 3, 2.5]);

        % Set default font size
        set(haxes,'FontSize',10)
        % Set interpreter to LaTeX
        set(0,'defaultTextInterpreter','latex')

        % Set the interpreter for latex
        set(haxes, 'TickLabelInterpreter', 'latex');

        % Find all line objects
        hline= findobj(haxes,'Type','line');
        % Set line properties (relatively fat lines)
        set(hline, 'LineWidth'   , 1);

        % Set axes properties
        set(haxes, ...
            'Box'         , 'on'     , ...  %puts a box around the axis
            'TickDir'     , 'in'     , ...
            'TickLength'  , [.02 .01] , ...
            'XMinorTick'  , 'off'      , ...
            'YMinorTick'  , 'off'     , ...
            'LineWidth'   , 1);

        function [xpos,ypos] = getPosition(n)
            % gets the x y position of the figure to plot
            n = mod(n-1,15)+1;
            ypos = floor((n-1)/5)*5;
            xpos = (mod(n-1,5))*3;
        end

    end

    function SaveFigure(filename)
        % Saves a current graphics file to file named filename.pdf

        % Get the handle
        h=gcf;
        set(h,'Units','Inches');
        pos = get(h,'Position');

        % Cut all the white margins
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

        % Actually save
        %print(h,['Figures/' filename, '.pdf'],'-dpdf','-r0')
        saveas(h, ['Figures/' filename],'epsc')

        % This outputs on the screen what figure was saved into what file
        disp_Diary(['Saved a figure ' num2str(get(gcf,'Number')) ' into a file ' ...
            filename '.pdf'])
    end


    function disp_Diary(string)
        % Displays text on the screen and also logs it in the diary
        diary on
        disp(string)
        diary off
    end


end
