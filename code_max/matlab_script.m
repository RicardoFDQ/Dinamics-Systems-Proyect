%% Parte 1 (preguntas 1, 2, 3, 4 y 11)
clear

R_earth = 6.378e6;

eta = 8e17;
m = 2e3;
r_ref = 8e6;

A = [0 1; 2*eta/(r_ref)^3  0];
B = [0; 1/m];
C = [1 0; 0 1];
D = [0; 0];

x_eq = [r_ref; 0];
u_eq = eta*m/(r_ref)^2;

x0 = [r_ref; 0];

controlability_matrix = ctrb(A, B);
rank_contr_matrix = rank(controlability_matrix);

Q = [1 0; 0 1];
R = 0.0001;
N = [0; 0];

K = lqr(A, B, Q, R, N);

perturbation_standard_deviation = 9.6e7;

% Ganancias backstepping
K1 = 1e-3;
K2 = 4.5;


%% Parte 2 (preguntas 5, 6, 7, 8 y 9)
clear

R_earth = 6.378e6;

eta = 8e17;
m = 2e3;
r_ref = 8e6;
theta_ref = 0;
d_r_ref = 0;
d_theta_ref = 1670 * 1000 / (3600 * R_earth);
u1_eq = m * ((eta/r_ref) - r_ref * d_theta_ref^2);
u2_eq = 2 * m * d_r_ref * d_theta_ref / r_ref;

x_eq = [r_ref; d_r_ref; theta_ref; d_theta_ref];
u_eq = [u1_eq; u2_eq];

A = [0 1 0 0; 
    (2*eta/r_ref^3)+d_theta_ref^2  0  0  2*r_ref*d_theta_ref;
    0 0 0 1;
    2*d_r_ref*d_theta_ref/r_ref^2  -2*d_theta_ref/r_ref  0  -2*d_r_ref/r_ref];
B = [0 0; 1/m 0; 0 0; 0 1/m];
C = eye(4);
D = zeros(4, 2);

controlability_matrix = ctrb(A, B);
rank_contr_matrix = rank(controlability_matrix);
[ctr_decomposition] = ctrbf(A, B, C);


Q = diag([1 1 1e1 1e20]);
R = diag([1e-4 1e6]);


K = lqr(A, B, Q, R);


B_p9 = [0; 0; 0; 1/m];
C_p9 = [1 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 1];

ctr_matrix_p9 = ctrb(A, B_p9);
obs_matrix_p9 = obsv(A, C_p9);

rank_ctr_matrix_p9 = rank(ctr_matrix_p9);
rank_obs_matrix_p9 = rank(obs_matrix_p9);

decomp_obs_p9 = obsvf(A, B_p9, C_p9)