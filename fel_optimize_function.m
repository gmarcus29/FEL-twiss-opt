% This is the merit function that the simplex fminsearch function is going
% to try and optimize.

% Going to start with a 3D gaussian of the twiss parameters.

function variables_struct = fel_optimize_function(variables_struct)

% -------------------------------------------------------------------------
% This section sets up the simplex solver.
% The optics/lattice file are called by the solver, so it is in the
% fel_merit member function at the bottom of this file.

% Seed the k values for the magnets.  This is what the simplex solver uses
% to start the lattice solver.
variables_struct.k1 = 2.0;
variables_struct.k2 = 0.1;
variables_struct.k3 = 1.3;

% These are the options for the optimizers - i.e. the simplex solver.
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
options = optimset('Display','none');
options.TolFun = 1e-4;
options.TolX = 1e-4;

% This is the actual simplex solver.  The file which defines the optics is
% in inside the fel_merit function.
fun = @(x) fel_merit(x, variables_struct);
x0 = variables_struct.start_point;
x = fminsearch(fun, x0, options);

% Save the twiss parameters of the result of the simplex search.
beta_solved = x(1);
alpha_solved = x(2);
gamma_solved = (1 + x(2)^2)/x(1);
variables_struct.T_solved =...
    [beta_solved , -alpha_solved ; -alpha_solved , gamma_solved];

% Save the linac lattic quad current solution so you can see what the hell
% it is doing.
gamma = (1 + x(2)^2)/x(1);
variables_struct.input_values = [x(1), -x(2); -x(2), gamma];
variables_struct = linac_lattice_solver(variables_struct);




function merit = fel_merit(x, variables_struct)

gamma = (1 + x(2)^2)/x(1);
variables_struct.input_values = [x(1), -x(2); -x(2), gamma];
variables_struct = linac_lattice_solver(variables_struct);

% Unmask the desired variables to make the code below more readable.
bm = variables_struct.matched_values(1,1);
am = variables_struct.matched_values(1,2);

b = variables_struct.transported_true_BP(1,1);
a = variables_struct.transported_true_BP(1,2);

% Make a width for the merit function
sig_b = 0.5;
sig_a = 0.5;
% Calculate the merit function
merit = -exp(-(b-bm)^2/2/sig_b/sig_b).*exp(-(a-am)^2/2/sig_a/sig_a);

% -------------------------------------------------------------------------
% Code to test/follow the search.  These are diagnostic lines.

disp([ 'k1: ', num2str(variables_struct.k1), ' k2: ', num2str(variables_struct.k2) ] )

bb = linspace(bm-2,bm+2,2^10);
aa = linspace(am-2,am+2,2^10);


[BB,AA] = meshgrid(bb,aa);
LOL = exp(-(BB-bm).^2/2/sig_b/sig_b).*exp(-(AA-am).^2/2/sig_a/sig_a);

figure(79)
set(gcf,'Color', 'w')
set(gcf,'Position', [-1175         340         560         420])
imagesc(bb,aa,LOL)
hold on;
plot(b,a,'co','MarkerSize',20)

figure(80)
set(gcf,'Color', 'w')
set(gcf,'Position', [-1750         340         560         420])
plot(variables_struct.k1,variables_struct.k2,'kd','MarkerSize',20)
xlim([-3 3])
ylim([-3 3])
xlabel('k1')
ylabel('k2')

% pause(0.01)
