% This is the function that "optimize_quads.m" attempts to optimize.

function M = linac_lattice(variables_struct)


Md = @(z) [1, z; 0, 1];
Mq = @(z) [1, 0; -z^2, 1];

% Build a simple lattice.
k1 = variables_struct.k1;
k2 = variables_struct.k2;

M = Md(2.5)...
    * Mq(k2) * Md(0.26) ...
    * Mq(1i*k1) * Md(1.23);



