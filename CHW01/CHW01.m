G1 =  gsp_graph([0 1.6 2.4; 1.6 0 0.8; 2.4 0.8 0], [1 1; 2 2; 3 1]);
G2 = gsp_graph([0 0.7 1.1 2.3; 0.7 0 0 0; 1.1 0 0 0; 2.3 0 0 0], [0 0; 0 1; -1 -1; 1 -1]);
figure(1);
gsp_plot_graph(G1);
figure(2);
gsp_plot_graph(G2);

Gs = gsp_graph_product(G1, G2);
Gs.coords = [0 0; 0.1 1; 0 2; -0.1 3; 0 4; 0.1 5; 1 0; 1.1 1; 1 2; 0.9 3; 1 4; 1.1 5];
figure(3);
gsp_plot_graph(Gs);

disp(Gs.A);
disp(Gs.L);

param.rule = 'kronecker';
Gt = gsp_graph_product(G1, G2, param);
Gt.coords = [0 0; 0.1 1; 0 2; -0.1 3; 0 4; 0.1 5; 1 0; 1.1 1; 1 2; 0.9 3; 1 4; 1.1 5];
figure(4);
gsp_plot_graph(Gt);

disp(Gt.A);
disp(Gt.L);

MyG = Gt;
signal = 20*rand(12, 1)-10;
figure(5);
gsp_plot_signal(MyG, signal);


MyG = gsp_compute_fourier_basis(MyG);
disp('Eigenvectors are:');
disp(MyG.U);
disp('Eigenvalues are:');
disp(MyG.e);

figure(6);
gsp_plot_signal_spectral(MyG,signal);

figure(7);
gsp_plot_signal(MyG, MyG.U(:,2));
figure(8);
gsp_plot_signal(MyG, MyG.U(:,3));
figure(9);
gsp_plot_signal(MyG, MyG.U(:,11));
figure(10);
gsp_plot_signal(MyG, MyG.U(:,12));

GL = gsp_logo;
signal1 = zeros(1130,1);
for i = 1:1130
    if GL.coords(i , 1) < 205
        signal1(i) = -1;
    else
        if GL.coords(i , 1) > 390
            signal1(i) = -0.5;
        else 
            signal1(i) = 1;
        end
    end
end

figure(11);
gsp_plot_signal(GL, signal1);

GL = gsp_compute_fourier_basis(GL);
disp('Eigenvectors are:');
disp(GL.U);
disp('Eigenvalues are:');
disp(GL.e);

B = zeros(1130, 2);
for i = 1:1130
    B(i, :) = [GL.U(i, 2), GL.U(i,3)];
end

figure(12);
gsp_plot_graph(gsp_graph(zeros(1130),B));


C = zeros(1130, 3);
for i = 1:1130
    C(i, :) = [GL.U(i, 2), GL.U(i,3), GL.U(i,4)];
end

figure(13);
gsp_plot_graph(gsp_graph(zeros(1130),C));