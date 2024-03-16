clear;
clc;

function x = MyVertexSelection(G)
G = gsp_compute_fourier_basis(G);
G.U
u_max = G.U(:, n);
V_1 = [0];
for i=1: n
    if u_max(i) >= 0
        V_1(end+1) = i;
    end
end
V_1(V_1 == 0) = [];
x = V1;
end


function x = MySKReduction(G,V1)
n1 = (size(G.W));
n = n1(1);
Q = 4*n*log10(n);
V1c = [];
V = [];
for i=1:n
V(i) = i;
end
V = V';
NewV1=zeros(n, 1);
for i=1: length(V1)
NewV1(V1(i)) = V1(i);
end

V1c = V - NewV1;
V1c(V1c == 0)= [];

newL = G.L(V1, V1) - (G.L(V1, V1c)*(G.L(V1c, V1c)^(-1))*(G.L(V1c, V1)));
newW = -newL + diag(diag(newL));
newNewW = zeros(length(V1),length(V1));
Z = [];
for i = 1:length(V1)
    Z(i) = i;
end
Z = Z';
resistanceDistance = zeros(length(V1),length(V1));
for i=1:length(V1)
    for j=1:length(V1)
       resistanceDistance(i,j) = (kroneckerDelta(Z,sym(i))-kroneckerDelta(Z,sym(j)))'*(pinv(newL))*(kroneckerDelta(Z,sym(i))-kroneckerDelta(Z,sym(j))); 
    end
end

P = zeros(length(V1),length(V1));
for i=1:length(V1)
    for j=1:length(V1)
        P(i,j) = (resistanceDistance(i,j)*newW(i,j))/sum(sum(resistanceDistance.*newW));
    end
end
for t=1:Q
    Pvector = reshape(P,[1, length(V1)^2]);
    y = randsample(length(V_1)^2,1,true,Pvector);
    j = rem(y, length(V1));
    i = ((y-j)/length(V1))+1;
    if j == 0
        j = length(V1);
    end
    if i == j
    newNewW(i,j) = 0;
    else
        newNewW(i,j) = newNewW(i,j) + newW(i,j)/(Q*P(i,j));
        newNewW(j,i) = newNewW(i,j);
    end
end
x = gsp_graph(newNewW);
end


function y = MyHfilter(x, G)
G = gsp_compute_fourier_basis(G); 
y = G.U * ((diag((G.e)*2+1))^(-1)) * (G.U)^(-1) * x;
end

function y = MyDS(x, V1)
y = [];
for i= 1 : length(V1)
y(i) = x(V1(i));
end
end

function f = MyInterpolate(downX, V1, G, epsilon, flag)
G = gsp_compute_fourier_basis(G);
n1 = size(G);
n = n1(1);
V = [];
for i=1:n
V(i) = i;
end
V = V';
NewV1=zeros(n, 1);
for i=1: length(V1)
NewV1(V1(i)) = V1(i);
end

V1c = V - NewV1;
V1c(V1c == 0)= [];


phi=zeros(n, length(V1));

if flag == 0
for j = 1 : length(V1)
sum = zeros(n,1);
    for i= 1 : n 
sum = sum + ((1/(G.e(i) + epsilon))* conj(G.U(j, i)) * G.U(:, i));
    end
    phi(:, j) = sum;
end
Lbar = G.L+ epsilon*eye(length(V1));
alfa = ((Lbar(V1, V1)) - (Lbar(V1, V1c)* ((Lbar(V1c, V1c))^(-1))) * Lbar(V1c, V1)) * downX;

f = phi * alfa;
end
if flag == 1
f = G.U * pinv(G.U(V1, :))* downX;
end

end

function [L_1, m_0, y_0, x_1] = MyAnalysis(G, x, epsilon, flag)
m_0 = MyVertexSelection(G);
L_1 = MySKReduction(G, m_0).L;
z = MyHfilter(x, G);
x_1 = MyDS(z,m_0);
y_0 = (MyInterpolate(x_1, m_0, G, epsilon,flag)) - x;
end

function [cL, cy, cm, cx] = MyPyramidAnalysis(G, x, N, epsilon, flag)
cL = cell(N);
cy = cell(N);
cm = cell(N);
cx = cell(N);
[cL(1), cm(1), cy(1), cx(1)]= MyAnalysis(G, x, epsilon,flag);
for i = 1 : N-1
   [cL(i+1), cy(i+1), cm(i+1), cx(i+1)] = MyAnalysis(G(i), x(i), epsilon, flag);
end
end

function x_0 = MySynthesis(x_1, G, y, m, epsilon)
n1 = size(G);
n = n1(1);
V1c = [];
V = [];
for i=1:n
V(i) = i;
end
V = V';
NewV1=zeros(n, 1);
for i=1: length(V1)
NewV1(V1(i)) = V1(i);
end

V1c = V - NewV1;
V1c(V1c == 0)= [];
phi=zeros(n, length(m));
for j = 1 : length(m)
sum = zeros(n,1);
    for i= 1 : n 
sum = sum + ((1/(G.e(i) + epsilon))* conj(G.U(j, i)) * G.U(:, i));
    end
    phi(:, j) = sum;
end
x_0 = [phi * phi(:, m)^(-1), eye(length(m))] *[x_1; y];

end

function x = MyPyramidSynthesis(N, cG, cy, cm, xN, epsilon)
x = xN;
for i = 1:N
x = MySynthesis(x, cG{i}, cy{i}, cm{i}, epsilon);
end
end
