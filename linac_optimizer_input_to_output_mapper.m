% This function maps a set of input guesses for the beta and alpha
% functions, runs the solver and then propagates the true beam parameters
% through the lattice that matches the guess to the desired output.  The
% idea is that we want to know if the input maps are one-to-one with the
% output maps that the FEL simplex solver uses.  Currently, the FEL solver
% seems to get confused so we suspect that alpha_i, beta_i are not one to
% one with alpha_f, beta_f.

clear all
close all hidden



% Turn on the plots to see what the optimizer is doing.
variables_struct.plots_on = 0;
% Turn on the step-by-step notification of quad settings.
variables_struct.k_notification = 0;

% The true beam parameters.  These are the beam parameters that you hope
% the optimizer will find.  When running this on the machine, these will be
% the unknown initial beam parameters (IBP) in the machine.
beta_true = 1.35;
alpha_true = 0.60;
gamma_true = (1 + alpha_true^2) / beta_true;
variables_struct.true_values = [beta_true, alpha_true; alpha_true, gamma_true];

% These are the desired values at the output, the match into the FEL.
% Matched Beam Parameters (MBP).
match1 = 2.5;
match2 = -0.5;
match3 = (1 + match2^2) / match1;
variables_struct.matched_values = [match1, match2; match2, match3];

% Seed the k values for the magnets.  This is what the simplex solver uses
% to start the lattice solver.
variables_struct.k1 = 2.0;
variables_struct.k2 = 0.1;
variables_struct.k3 = 1.3;

% Now iterate over input alpha and beta and run the solver to 
Nk = 2^4; % The number of beta values to use.
Nj = 2^4; % The number of alpha values to use.
beta_min = 0.01;
beta_max = 1.5;
alpha_max = 1.0; % Alpha goes from -alpha_max to alpha_max;

beta_list = linspace(beta_min, beta_max, Nk);
alpha_list = linspace(-1.0, 1.0, Nj);

% Build the array to store the output.
beam_param_inputs = zeros(Nk*Nj,2); % 1 is beta, 2 is alpha;
beam_param_results = zeros(Nk*Nj,2); % 1 is beta, 2 is alpha;

NN = 0;

for k = 1 : Nk
    for j = 1 : Nj
        NN = NN + 1;
        
        % The beta, alpha and gamma to feed into the solver.
        test1 = beta_list(k); % A start beta
        test2 = alpha_list(j); % A start alpha
        
        % Calculate gamma.
        test3 = (1 + test2^2)/test1; % Derive the start gamma (not used)
        % This just combines the start points to transfer into the solver.
        variables_struct.input_values = [test1, test2; test2, test3];
        % Now run the solver.
        variables_struct = linac_lattice_solver(variables_struct);
        
        beam_param_inputs(NN,1) = test1;
        beam_param_inputs(NN,2) = test2;
        
        beam_param_results(NN,1) = variables_struct.transported_true_BP(1,1);
        beam_param_results(NN,2) = variables_struct.transported_true_BP(1,2);
        
    end
end

%%

% Now calculate the FEL merti function for the transported true values to
% see what the FET merit function space looks like.

% Unmask the desired variables to make the code below more readable.
bm = variables_struct.matched_values(1,1);
am = variables_struct.matched_values(1,2);

b = beam_param_results(:,1);
a = beam_param_results(:,2);

% Make a width for the merit function
sig_b = 1.5;
sig_a = 1.5;
% Calculate the merit function
merit = exp(-(b-bm).^2/2/sig_b/sig_b).*exp(-(a-am).^2/2/sig_a/sig_a);


%%

figure(13231)
set(gcf, 'Color', 'w')
set(gcf, 'Position', [-1075         261        1017         512])

subplot(1,2,1)
plot(beam_param_inputs(:,1),beam_param_inputs(:,2),'.')
hold on;
plot(variables_struct.true_values(1,1),variables_struct.true_values(1,2), 'rx',...
    'MarkerSize', 20)
hold on;
plot(variables_struct.matched_values(1,1),variables_struct.matched_values(1,2)...
    , 'ko', 'MarkerSize', 20)
xlabel('\beta [m]', 'FontSize', 20)
ylabel('alpha [1]', 'FontSize', 20)
pad = beta_max* 0.05;
xlim([-pad, beta_max + pad])
ylim([-1.1*alpha_max 1.1*alpha_max])
title('Guess Beam Parameters', 'FontSize', 20)

subplot(1,2,2)
plot(beam_param_results(:,1),beam_param_results(:,2),'.')
hold on;
plot(variables_struct.true_values(1,1),variables_struct.true_values(1,2), 'rx',...
    'MarkerSize', 20)
hold on;
plot(variables_struct.matched_values(1,1),variables_struct.matched_values(1,2)...
    , 'ko', 'MarkerSize', 20)
xlabel('\beta [m]', 'FontSize', 20)
ylabel('alpha   [1]', 'FontSize', 20)
xlim([0, 1.1*max(beam_param_results(:,1))])
ylim([1.1*min(beam_param_results(:,2)) 1.1*max(beam_param_results(:,2))])
title('Transported True Beam Parameters', 'FontSize', 20)
legend('Data','True BP','Matched BP')

%%

figure(13232)
set(gcf, 'Color', 'w')
set(gcf, 'Position', [-1684         353         560         420])
plot3(b,a,merit, '.')
xlabel('\beta [m]', 'FontSize', 20)
ylabel('\alpha [m]', 'FontSize', 20)
zlabel('Merit', 'FontSize', 20)
