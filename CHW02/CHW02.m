G1 = gsp_graph([0 1 1 0 0 0 0 1; 1 0 1 1 1 0 0 1; 1 1 0 1 0 0 0 0; 0 1 1 0 1 1 0 1; 0 1 0 1 0 1 1 1; 0 0 0 1 1 0 1 0; 0 0 0 0 1 1 0 0; 1 1 0 1 1 0 0 0],[0 0; 1 1; 0.5 3; 4 3;4 1; 5 2.5; 4.5 0.7; 3 0]);
%gsp_plot_graph(G1);
G1 = gsp_compute_fourier_basis(G1);
U = G1.U;
x = 2*U(:,1) + U(:,2);
y = awgn(x,10);
figure(1);
gsp_plot_signal(G1,x);
figure(2);
gsp_plot_signal(G1,y);

V = G1.W * G1.W;
sigma = eig(V,"vector");
Wnorm = (1/sqrt(sigma(8)))*G1.W;
[U1, D1] = eig(Wnorm);
xhat1 = U1^(-1)*x;
xhat2 = (G1.U)^(-1)*x;
lambdas = D1*ones(8,1);
figure(3);
gsp_plot_signal_spectral(G1,xhat2);
figure(4);
stem(lambdas, xhat1);

Z1 = xhat2(1)* G1.U(:, 1) + xhat2(2,1)*G1.U(:, 2);
figure(5);
gsp_plot_signal(G1,Z1);


Z3= xhat1(7)* G1.U(:, 7) + xhat1(8)* G1.U(:, 8);
figure(6);
gsp_plot_signal(G1, Z3);

yhat1 = U1^(-1)*y;
yhat2 = (G1.U)^(-1)*y;


Z2 = yhat2(1)* G1.U(:, 1) + yhat2(2)*G1.U(:, 2);
figure(7);
gsp_plot_signal(G1, Z2)

Z4 = yhat1(7)*G1.U(:, 7) + yhat1(8)* G1.U(:,8);
figure(8);
gsp_plot_signal(G1, Z4);

h1 = [1; 1; 0; 0; 0; 0; 0; 0];
h2 = [0; 0; 0; 0; 0; 0; 1; 1];

snr(x,Z2)
snr(x,Z4)

S1 = ((G1.U)^(-1))*diag(h1)*(G1.U);
if S1*G1.L == G1.L*S1
    disp('is LSI.')
else
    disp('not LSI.')
end

S2 = (U1^(-1)) * diag(h2) * U1;
if S2*Wnorm == Wnorm*S2
    disp('is LSI.')
else
    disp('not LSI.')
end


