% This function for trying to optimize a set of quads using matlab
% fminsearch to perform simplex optimation.
%
% The goal of this program is to search the set of initial beam parameters,
% which define the beam parameters before the matching optics which prepare
% the beam for the FEL and use the response of the FEL to determine the
% fitness for the beam for...the FEL.
%
% This whole system, which is going to be called the "optimizer" has two
% main components: a simplex solver and a lattice/optics solver.  The
% simplex solver is the over-arching solver that is used to perform the
% optimization.  The optics solver is the routine that the simplex solver
% uses to simulate the optics between the initial beam and the input into
% the FEL.

clear all
close all hidden

global struct fel_transported_BP
fel_transported_BP.N_FEL = 0;
% This number will be one larger than the N_FEL because of the call to
% linace_lattice_solver at the end of fel_optimize_function.
fel_transported_BP.N_LATTICE = 0;

% Turn on the plots to see what the optimizer is doing.
variables_struct.plots_on = 0;
% Turn on the step-by-step notification of quad settings.
variables_struct.k_notification = 0;


% The true beam parameters.  These are the beam parameters that you hope
% the optimizer will find.  When running this on the machine, these will be
% the unknown initial beam parameters (IBP) in the machine.

beta_true = 1.35;
alpha_true = -0.40;
gamma_true = (1 + alpha_true^2) / beta_true;
variables_struct.true_values = [beta_true, alpha_true; alpha_true, gamma_true];

% These are the desired values at the output, the match into the FEL.
% Matched Beam Parameters (MBP).
match1 = 2.5;
match2 = -0.0;
match3 = (1 + match2^2) / match1;
variables_struct.matched_values = [match1, match2; match2, match3];

% Seed the simplex search with beta and alpha.  Since you don't know what
% the initial beam parameters in the real world you have to start with a
% guess.  There may be ways to use other methods developed by other people
% to seed this better.  A better seed is a faster solve.
test1 = 1.01; % A start beta
test2 = -0.5; % A start alpha
test3 = (1 + test2^2)/test1; % Derive the start gamma (not used)
% This just combines the start points to transfer into the solver.
variables_struct.start_point = [test1, test2];

% Perform the simplex optimization.  The actual simplex solver is contained
% in the fel_merit_function.
variables_struct = fel_optimize_function(variables_struct);

% Print the result.
disp('Solution from Solver:')
variables_struct.T_solved
disp('Desired Solution:')
variables_struct.true_values
disp('FEL Match:')
variables_struct.matched_values
disp('Start Point:')
[test1, test2; test2, test3]


%%
% Now plot the steps the solver goes through.

Nk = fel_transported_BP.N_FEL;

beta_steps_guess = zeros(Nk,1);
alpha_steps_guess = zeros(Nk,1);
beta_steps_trans_true = zeros(Nk,1);
alpha_steps_trans_true = zeros(Nk,1);
MERITS = zeros(Nk,1);
eflag_color = zeros(Nk,3);

for k = 1 : Nk
    
    beta_steps_guess(k,1) = fel_transported_BP.guess_BP{k}(1,1);
    alpha_steps_guess(k,1) = fel_transported_BP.guess_BP{k}(1,2);
    
    beta_steps_trans_true(k,1) = fel_transported_BP.transported_true_BP{k}(1,1);
    alpha_steps_trans_true(k,1) = fel_transported_BP.transported_true_BP{k}(1,2);
    
    MERITS(k,1) = fel_transported_BP.merit{k};
    
    if fel_transported_BP.lattice_solver_eFlag{k} == 1
        eflag_color(k,:) = [0 0 1];
    else
        eflag_color(k,:) = [1 0 0];
    end
    
end

bm = variables_struct.matched_values(1,1);
am = variables_struct.matched_values(1,2);
b = beta_steps_trans_true;
a = alpha_steps_trans_true;

% Make a width for the merit function
sig_b = 1.5;
sig_a = 1.5;
% Calculate the merit function
merit = exp(-(b-bm).^2/2/sig_b/sig_b).*exp(-(a-am).^2/2/sig_a/sig_a);



figure(1323234)
set(gcf, 'Color', 'w')
set(gcf, 'Position', [-1075         261        1017         512])

subplot(1,2,1)
set(gca, 'FontSize', 20)
for k = 1 : Nk
    plot(beta_steps_guess(k), alpha_steps_guess(k), '.',...
        'Color', eflag_color(k,:));
    hold on;
end

plot(variables_struct.true_values(1,1), variables_struct.true_values(1,2), 'rx',...
    'MarkerSize', 20, 'LineWidth', 2)
% hold on;
% plot(variables_struct.matched_values(1,1), variables_struct.matched_values(1,2)...
%     , 'ro', 'MarkerSize', 20, 'LineWidth', 2)
xlabel('\beta [m]', 'FontSize', 20)
ylabel('\alpha [1]', 'FontSize', 20)
title('Guess Beam Parameters', 'FontSize', 20)
legend('Guess BP','True BP')
% Set the colors for the beta/alpha steps based on whether the lattice
% solver worked.
% set(h, {'color'}, num2cell(eflag_color,2))

subplot(1,2,2)
set(gca, 'FontSize', 20)
for k = 1 : Nk
    plot(beta_steps_trans_true(k), alpha_steps_trans_true(k),'.',...
        'Color', eflag_color(k,:))
    hold on;
end
plot(variables_struct.matched_values(1,1), variables_struct.matched_values(1,2)...
    , 'ro', 'MarkerSize', 20, 'LineWidth', 2)
xlabel('\beta [m]', 'FontSize', 20)
ylabel('\alpha [1]', 'FontSize', 20)
title('Transported True Beam Parameters', 'FontSize', 20)
legend('Transported True BP','Matched BP')

%%
% Generate the idea merit function.
NB = 2^6;
B = linspace(0,2*bm, NB);
A = linspace(-2*am-2,2*am+2, NB);
[BB, AA] = meshgrid( B, A);
sig_b = 1.5;
sig_a = 1.5;
merit_true = exp(-(BB-bm).^2/2/sig_b/sig_b).*exp(-(AA-am).^2/2/sig_a/sig_a);


figure(345243)
set(gcf, 'Position', [-1687         315         560         420])
set(gcf, 'Color', 'w')
imagesc(B, A, merit_true)
hold on;
plot(beta_steps_trans_true, alpha_steps_trans_true, 'k.',...
    'MarkerSize',20)
hold on;
plot(beta_steps_trans_true(end), alpha_steps_trans_true(end), 'cx',...
    'MarkerSize',20, 'LineWidth', 2)
xlabel('\beta [m]', 'FontSize', 20)
ylabel('\alpha [1]', 'FontSize', 20)
title('FEL Merit Function and Transported True BP', 'FontSize', 20)
set(gca, 'FontSize', 20)










