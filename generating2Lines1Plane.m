D = 3; %Dimension of ambient space
n = 3; %Number of subspaces
d1 = 1; d2 = 1; d3 = 2;%d1 and d2: dimension of subspace 1 and 2
N1 = 50; N2 = 50; N3 = 50; %N1 and N2: number of points in subspace 1 and 2
X1 = randn(D,d1) * randn(d1,N1); %Generating N1 points in a d1 dim. subspace
X2 = randn(D,d2) * randn(d2,N2); %Generating N2 points in a d2 dim. subspace
X3 = randn(D,d3) * randn(d3,N3);

X = [X1 X2 X3];

X1 = X1';
X2 = X2';
X3 = X3';

scatter3(X1(:,1), X1(:,2), X1(:,3));
hold on
scatter3(X2(:,1), X2(:,2), X2(:,3));
hold on
scatter3(X3(:,1), X3(:,2), X3(:,3));

