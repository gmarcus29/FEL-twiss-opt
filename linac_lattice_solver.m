% This function takes in a matched set of beam parameters (MBP) and a guess
% at the true BP called the "input beam parameters" IBP.  It then finds a
% solution for the lattice parameters that transform the IBP into the MPB.

function variables_struct = linac_lattice_solver(variables_struct)

% Setup the fsolver, which is the optics solver. (not simplex)
options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
options.TolFun      = 1e-12;
options.TolX        = 1e-12;
options.MaxIter     = 1e4;
options.MaxFunEvals = 1e4;

% Try to take the IBP and match it to the MBP.  This solves for the k
% values that match the IBP to the MBP.
x0 = [variables_struct.k1, variables_struct.k2];
fun = @(x) linac_lattice_fsolve(x, variables_struct);

% Run the solvers to get the lattice quad focusing values.
% [x,fval,exitflag,output,J] = fsolve(fun, x0, options );
[x,~,eFlag] = fsolve(fun, x0, options );
variables_struct.k1 = x(1);
variables_struct.k2 = x(2);

% Run the true beam parameters through the transport:
M = linac_lattice(variables_struct);
variables_struct.transported_true_BP = M * variables_struct.true_values * M';


%----- Save the solution flag from the lattice solver
% this is to find out whether it is finding a good solution.
global struct fel_transported_BP

% Increment the lattice_solver call count
fel_transported_BP.N_LATTICE = fel_transported_BP.N_LATTICE + 1;

% Save the present guess beam parameters.
fel_transported_BP.lattice_solver_eFlag{fel_transported_BP.N_LATTICE} = eFlag;


%-----



% This is the function which fsolve...solves.
function f = linac_lattice_fsolve(x, variables_struct)

% Unmask the quad values, because I designed the lattic solver to us the
% variable_struct to transfer data around.
variables_struct.k1 = x(1);
variables_struct.k2 = x(2);


X = variables_struct.input_values;
M = linac_lattice(variables_struct);


Mf = M * X * M' ;

f(1) = Mf(1,1) - variables_struct.matched_values(1,1);
f(2) = Mf(1,2) - variables_struct.matched_values(1,2);


% -------------------------------------------------------------------------
% Below this line are all diagnostics, they can be on or off.  Who cares.
% disp([ 'k1: ', num2str(x(1)), ' k2: ', num2str(x(2)) ] )







