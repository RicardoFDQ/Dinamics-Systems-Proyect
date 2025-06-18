%% 2
A = [0 1; 2*8e17/8e6^3 0];
B = [0;1/200];
C = eye(2);

Q = diag([1,0.1]);
R = 1;

Klqr = lqr(A,B,Q,R,0);

K = place(A,B,[-1,-2]);

disp(K)