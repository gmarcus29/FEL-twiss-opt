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

% The true beam parameters.  These are the beam parameters that you hope
% the optimizer will find.  When running this on the machine, these will be
% the unknown initial beam parameters (IBP) in the machine.
beta_true = 1.35;
alpha_true = 0.58;
gamma_true = (1 + alpha_true^2) / beta_true;
variables_struct.true_values = [beta_true, -alpha_true; -alpha_true, gamma_true];

% These are the desired values at the output, the match into the FEL.
% Matched Beam Parameters (MBP).
variables_struct.matched_values = [4.0, 0; 0, 0.25];

% Seed the simplex search with beta and alpha.  Since you don't know what
% the initial beam parameters in the real world you have to start with a
% guess.  There may be ways to use other methods developed by other people
% to seed this better.  A better seed is a faster solve.
test1 = 2.5; % A start beta
test2 = 1.0; % a start alpha
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
