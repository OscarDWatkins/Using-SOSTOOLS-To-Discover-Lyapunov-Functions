%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITLE: Lyapunov Function Search
%
% The system and substitution are from M. Anghel et al. 2013 - Algorithmic
% Construction of Lyapunov Functions via SOS Mehtods
%
% The process for finding V and expanding the ROA are based of 2024
% paper by Liu et al.- Estimating the region of attraction of wind integrated power
% systems based on improved expanding interior algorithm
%
% DEPENDENCIES:
% Please ensure the following files/libraries are in your MATLAB path
% before running this script:
%
%   * SOSTOOLS-SOSTOOLS400
%   * mosek 
%     In my final report I used SeDumi to be consitent with Anghel, also SeDumi is easier to access
%     However, mosek is much faster
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% GUIDE TO USING SOSTOOLS
% There are two main ways to transform from Trig to Polynomial system - Taylor series or Substitution

% Taylor Series
% Set gn = 0
% Llamda_n can be ignored for the whole script
% Sn are SOS polynomials

% Substitution 
% Choose gn depending on your substitution method
% Llamda_n are just normal polynomials
% Sn are SOS polynomials


% Choosing Degree
% V(x) must be even degree and not have linear terms as V(0) = 0
% sn must be even degree
% In my experience, SOS programs work better with no odd degree terms even if this isn't strictly necessary - this restricts search as well
% Eg use monomials(vars, [2,4]) NOT monomials(vars, 2:4)


%% DEFINING SYSTEM
pvar z1 z2 z3 z4 z5 z6

% Two Machine vs Infinite Bus System from M. Anghel et al.
g_all = [(z3-z2*z3);
    (z1*z3);
    (33.6*z2-67.8*z1-1.89*z3+16.9715*z4+1.9718*z5-1.9718*z1*z4+16.9715*z1*z5-16.9715*z2*z4-1.9718*z2*z5);
    (z6-z5*z6);
    (z4*z6);
    (11.3986*z1+1.2088*z2-98.8604*z4+48.4810*z5-1.2658*z6-1.2088*z1*z4-11.3986*z1*z5+11.3986*z2*z4-1.2088*z2*z5)];

% These Constraints change based on your substitutions
g1 = z1^2 + z2^2 - 2*z2;
g2 = z4^2 + z5^2 - 2*z5;

%% STEP 1
% Set Beta, i = 0 and choose initial p(z)
tol = 1e-5;
vars = [z1; z2; z3; z4; z5; z6];
eps = 1e-5;
beta = 5e-3;
z_sqr = z1^2 + z2^2 + z3^2 + z4^2 + z5^2 + z6^2;
l = eps * z_sqr;
pz = z_sqr;

%% STEP 2
% Solve SOSP0a - Finding Initial V

% Initialise SOS program
prog0a = sosprogram(vars);

% Create Variables for SOS program
[prog0a, V] = sospolyvar(prog0a, monomials(vars, 2:2));

[prog0a, llamda1] = sospolyvar(prog0a, monomials(vars, 0:0));
[prog0a, llamda2] = sospolyvar(prog0a, monomials(vars, 0:1));
[prog0a, llamda1b] = sospolyvar(prog0a, monomials(vars, 0:0));
[prog0a, llamda2b] = sospolyvar(prog0a, monomials(vars, 0:1));
[prog0a, s4] = sossosvar(prog0a, monomials(vars, 0:2));

% Defining Vdot
Vdot = jacobian(V, vars) * g_all;

% Defining Constraints
prog0a = sosineq(prog0a, V - llamda1*g1 - llamda1b*g2 - l);
prog0a = sosineq(prog0a, - s4*(beta - pz) - Vdot - llamda2*g1 - llamda2b*g2 -l);

% Calling Solver
solver_opt.solver = 'mosek';
[prog0a, info] = sossolve(prog0a, solver_opt);

% Evaluate feasibility logic
% - info.pinf == 0     : Not proven primal infeasible
% - info.dinf == 0     : Not proven dual infeasible
% - info.numerr < 2    : No severe numerical errors (0 or 1 is acceptable)
% - info.feasratio > 0.9: Feasibility metric close to 1.0 
feasible = (info.pinf == 0) && (info.dinf == 0) && (info.numerr < 2) && (info.feasratio > 0.9);

if feasible
    % Only extract and clean solution if solver succeeded
    V_sol = sosgetsol(prog0a, V);
    V0a = cleanpoly(V_sol, tol);
    disp("Feasible solution found.");
else
    disp("Problem is infeasible or solver failed.");
    return
end




%% STEP 3
% Solve SOSP0b - Maximising Gamma (gam)

% Bisection Search - If search fails, increase gam-high
gam_low = 0;
gam_high = .8;

while (gam_high - gam_low) > tol
    gam = (gam_high + gam_low)/2;

    % Initialise SOS program
    prog0b = sosprogram(vars);

    % Create Variables for SOS program
    [prog0b, llamda2_prime] = sospolyvar(prog0b, monomials(vars, 0:1));
    [prog0b, llamda2_primeb] = sospolyvar(prog0b, monomials(vars, 0:1));
    [prog0b, s7] = sossosvar(prog0b, monomials(vars, 0:2));
    
    % Defining Vdot
    V0adot = jacobian(V0a, vars) * g_all;

    % Defining Constraints
    prog0b = sosineq(prog0b, - s7*(gam - V0a) - V0adot - llamda2_prime*g1 - llamda2_primeb*g2 - l);

    % Calling Solver
    solver_opt.solver = 'mosek';
    [prog0b,info] = sossolve(prog0b,solver_opt);

    % Evaluate feasibility logic
    % - info.pinf == 0     : Not proven primal infeasible
    % - info.dinf == 0     : Not proven dual infeasible
    % - info.numerr < 2    : No severe numerical errors (0 or 1 is acceptable)
    % - info.feasratio > 0.9: Feasibility metric close to 1.0 (CRITICAL CHECK)
    feasible = (info.pinf == 0) && ...
        (info.dinf == 0) && ...
        (info.numerr < 2) && ...
        (info.feasratio > 0.9);

    if feasible
        gam_low = gam;
    else
        gam_high = gam;
    end
end

% Changes coefficients so gamma = 1    Part of Liu's process
V0b = V0a/gam_low;

displaySymSolution(V0b)

%% STEP 5
% Solve SOSP1a to get s10 and s11
done = false;

V1b(1) = V0b;

% Increase i for better result
for i = 1:50

    V_current = V1b(i);

    % Initialise SOS program
    prog1a = sosprogram(vars);

    % Create Variables for SOS program
    [prog1a, llamda5_prime] = sospolyvar(prog1a, monomials(vars, 0:1));
    [prog1a, llamda5_primeb] = sospolyvar(prog1a, monomials(vars, 0:1));
    [prog1a, s10] = sossosvar(prog1a, monomials(vars, 0:2));

    % Defining Vdot
    V1dot_numeric = jacobian(V_current, vars) * g_all;

    s10_strict = s10;

    % Defining Constraints
    prog1a = sosineq(prog1a, - s10_strict*(1 - V_current) - V1dot_numeric - llamda5_prime*g1 - llamda5_primeb*g2 - l);

    % Calling Solver
    solver_opt.solver = 'mosek';
    [prog1a, info] = sossolve(prog1a, solver_opt);

    % Evaluate feasibility logic
    % - info.pinf == 0     : Not proven primal infeasible
    % - info.dinf == 0     : Not proven dual infeasible
    % - info.numerr < 2    : No severe numerical errors (0 or 1 is acceptable)
    % - info.feasratio > 0.9: Feasibility metric close to 1.0
    feasible = (info.pinf == 0) && (info.dinf == 0) && (info.numerr < 2) && (info.feasratio > 0.9);


    if feasible
        s10_sol = sosgetsol(prog1a, s10_strict);
        s11_sol = 1;

        % Cleaning Solutions
        s10_cln = cleanpoly(s10_sol, tol);
        disp("Feasible")
    else
        disp("Infeasible")
        return 
    end


    %% STEP 6
    % Solve SOSP1b using s10, s11 and V1(i)

    % Initialise SOS program
    prog1b = sosprogram(vars);


    % Create Variables for SOS program
    [prog1b, V1_temp] = sospolyvar(prog1b, monomials(vars, 1:2));
    [prog1b, llamda3] = sospolyvar(prog1b, monomials(vars, 0:0));
    [prog1b, llamda4] = sospolyvar(prog1b, monomials(vars, 0:1));
    [prog1b, llamda5] = sospolyvar(prog1b, monomials(vars, 0:1));
    [prog1b, llamda3b] = sospolyvar(prog1b, monomials(vars, 0:0));
    [prog1b, llamda4b] = sospolyvar(prog1b, monomials(vars, 0:1));
    [prog1b, llamda5b]= sospolyvar(prog1b, monomials(vars, 0:1));
    [prog1b, s9] = sossosvar(prog1b, monomials(vars, 0:0));

    % Defining Vdot
    V1dot_symbolic = jacobian(V1_temp, vars) * g_all;

    % Defining Constraints
    prog1b = sosineq(prog1b, - s9*(1 - V_current) + (1 - V1_temp) - llamda3*g1 - llamda3b*g2);
    prog1b = sosineq(prog1b, V1_temp - llamda4*g1 - llamda4b*g2 - l);
    prog1b = sosineq(prog1b, - s10_cln*(1 - V1_temp) - 1*V1dot_symbolic - llamda5b*g1 - llamda5*g2 - l);

    % Calling Solver
    solver_opt.solver = 'mosek';
    [prog1b, info] = sossolve(prog1b, solver_opt);

    % Evaluate feasibility logic
    % - info.pinf == 0     : Not proven primal infeasible
    % - info.dinf == 0     : Not proven dual infeasible
    % - info.numerr < 2    : No severe numerical errors (0 or 1 is acceptable)
    % - info.feasratio > 0.9: Feasibility metric close to 1.0
    feasible = (info.pinf == 0) && (info.dinf == 0) && (info.numerr < 2) && (info.feasratio > 0.9);

    if feasible
        V1b_sol = sosgetsol(prog1b, V1_temp);
        V1b_cln = cleanpoly(V1b_sol, tol);
        V1b(i+1) = V1b_cln;
    else
        disp("Infeasible")
    return
    end


    if i > 1
        v_diff = V1b(i+1) - V1b(i);
        max_coeff = max(abs(v_diff.coefficient));

        if max_coeff < tol
            disp('Outer loop (Shape iteration) converged! Global success.')
            done = true;
        end
    end

    if done == true
        break;
    end
end

%% =========================================================================
% TRIGONOMETRIC CONVERSION ALGORITHM 
% =========================================================================

% Step 1: Convert the SOSTOOLS polynomial to a multi-line char matrix
char_matrix = char(V1b(i));

% Convert matrix rows to cells and join them into a single 1D continuous string
char_flat = strjoin(cellstr(char_matrix), ' ');

% Step 2: Explicitly define symbolic variables for both coordinate spaces
syms z1 z2 z3 z4 z5 z6
syms x1 x2 x3 x4

V_sym_z = str2sym(char_flat);

% Step 3: Execute the inverse coordinate transformation via symbolic substitution
% This has to change based on substitution method
V_trig_x = subs(V_sym_z, ...
    [z1,         z2,          z3, z4,         z5,          z6], ...
    [sin(x1), 1 - cos(x1),    x2, sin(x3), 1 - cos(x3),    x4]);

% Step 4: FORCE MATLAB to expand combined angles (e.g., cos(x1+x4) -> cos(x1)*cos(x4) - sin(x1)*sin(x4))
% This prevents the Python SMT verifier from choking on multi-variable trig arguments.
V_trig_x = expand(V_trig_x);

disp(' ');
disp('=======================================================');
disp(' Lyapunov Function in Expanded Trig Form (x variables):');
disp(V_trig_x);
disp('=======================================================');

%%
for m = 1:10
    % --- OPTIONAL SANITY CHECK ---
    % Verify that all those constant terms cancel out at the origin
    V_origin_check = subs(V_trig_x, [x1, x2, x3, x4], [0, 0, 0, 0]);
    disp('Value of V(x) at the origin (must be 0):');
    disp(double(V_origin_check));
end

%%

syms x1 x3 x2 x4

%% 1. Define Lyapunov Functions (Symbolic)
% You can compare up to 6 functions 
% This display in the order they are written, so will have to be manually moved around to stop one ROA blocking a smaller one
V1 = V_trig_x;
V2 = 0;
V3 = 0;
V4 = 0;
V5 = 0;
V6 = 0;

%% 2. Configuration & Color Selection
% Order functions from background to foreground. V2 is listed first so it sits in the background.
V_sym_list = {V2, V1, V3, V4, V5, V6};

% Unique RGB colors for each function:
colors = [
    0.500, 0.500, 0.500;  % V2: Medium Grey
    0.000, 0.000, 0.000;  % V1: Black
    0.850, 0.325, 0.098;  % V3: Orange-Red
    0.000, 0.447, 0.741;  % V4: Blue
    0.466, 0.674, 0.188;  % V5: Green
    0.494, 0.184, 0.556   % V6: Purple
    ];

%% 3. Grid Definition & Evaluation
% This gridspace numbers are obviously chosen for this specific system

x1_range = linspace(-2.5, 2.5, 1000);
x3_range = linspace(-2.5, 2.5, 1000);
[X1, X3] = meshgrid(x1_range, x3_range);

%% 4. Combined Plotting
figure; hold on;
set(gca, 'Color', 'w');

% Process and draw each Lyapunov function for V(x) <= 1
for j = 1:length(V_sym_list)
    % Substitute x2 = 0 and x4 = 0
    new_V = subs(V_sym_list{j}, [x2, x4], [0, 0]);

    % Convert symbolic equation to numeric function
    V_num = matlabFunction(new_V, 'Vars', [x1, x3]);

    % Evaluate on the grid
    Vvals = V_num(X1, X3);

    % Expand scalar outputs to full grid size
    if isscalar(Vvals)
        Vvals = repmat(Vvals, size(X1));
    end

    % Render region V <= 1
    draw_filled_blob(x1_range, x3_range, Vvals, 1, colors(j, :));
end

% Axis formatting
axis equal;
axis([-2.5 2.5 -2.5 2.5]);
xlabel('x_1 = \delta_1');
ylabel('x_3 = \delta_2');
box on;
hold off;

%% 5. Helper Functions
function draw_filled_blob(x1_range, x3_range, Zvals, level, faceColor)
C = contourc(x1_range, x3_range, Zvals, [level level]);
idx = 1;
while idx < size(C,2)
    n  = C(2,idx);
    xb = C(1, idx+1 : idx+n);
    yb = C(2, idx+1 : idx+n);
    patch(xb, yb, faceColor, 'EdgeColor', 'k', 'LineWidth', 1.2);
    idx = idx + n + 1;
end
end




