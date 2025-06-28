%% 2
A = [0 1; 2*8e17/8e6^3 0];
B = [0;1/2000];
C = eye(2);

Q = diag([1,1]);
R = 1e-4;

Klqr = lqr(A,B,Q,R,0);


disp(Klqr)

%% 4

Q_k = diag(9.6e7^2*[0,1]);
R_k = diag(1e-5*[1,1]);

[L, ~, ~] = lqe(A, eye(2), C, Q_k, R_k);
disp(L)

%% 5
m = 2e3; %m
mu = 8e17; %m3/s2
r_ref = 8e6; %m
theta_tierra = 1670*1000/3600 / 6.378e6; %m/s
r1 = 0;
%Para la trayectoria definida:

fr = m*(mu/r_ref^2-r_ref*theta_tierra^2);
ftheta = 0;

A = [0 0 1 0;
     0 0 0 1;
     2*mu/r_ref^3+theta_tierra^2, 2*r_ref*theta_tierra 0 0
     0 0 -2*theta_tierra/r_ref 0];

B = [0 0;0 0;1/m 0;0 1/m];

C = eye(4);

Co = ctrb(A,B);
[Ct,Cn] = ctrbf(A,B,C);


Q = diag([1e0;1e1;1;1e20]);
R = diag([1e-4,1e6]);

klqr = lqr(A,B,Q,R);

klqr = klqr.*[1 0 1 1;1 0 1 1]

%% 8
m = 2e3; %m
mu = 8e17; %m3/s2
r_ref = 8e6; %m
theta_tierra = 1670*1000/3600 / 6.378e6; %m/s
theta9 = sqrt(mu/r_ref^3);
%Para la trayectoria definida:

fr = m*(mu/r_ref^2-r_ref*theta_tierra^2);
ftheta = 0;

A = [0 0 1 0;
     0 0 0 1;
     2*mu/r_ref^3+theta_tierra^2, 2*r_ref*theta_tierra 0 0
     0 0 -2*theta_tierra/r_ref 0];

B = [0 0;0 0;1/m 0;0 1/m];

C = eye(4);

Q = diag([1e0;1e1;1e0;1e20]);
R = diag([1e-4,1e6]);

klqr = lqr(A,B,Q,R);

%Perturbaciones y ruidos

pertu_a = 1;
pr = (9.6e7*0.01)^2*pertu_a;
pte = (0.23*0.01)^2*pertu_a;

Q_k = diag([0;0;pr;pte]);

ruido_m = 1;
nr = (8-6.3)*1e5^2*ruido_m;
nte = 1e-5^2*ruido_m;
nrp = 1e4^2*ruido_m;
ntep = 0.727e-5^2*ruido_m;

R_k = diag([nr;nte;nrp;ntep^(1/2)]);


%% 9 Controlabilidad y Observabilidad

A = [0 0 1 0;
     0 0 0 1;
     2*mu/r_ref^3+theta_tierra^2, 2*r_ref*theta_tierra 0 0
     0 0 -2*theta_tierra/r_ref 0];

Bn = [0;0;0;1/m];
Cn = diag([1 0 0 1]);

ctr9 = ctrb(A,Bn);
obs9 = obsv(A,Cn);

rank(ctr9)
rank(obs9)