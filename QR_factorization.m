close all
clear
clc
% 数据初始化。
A = [1 2 3 4; 1 4 5 6; 1 5 6 7; 1 8 9 10; 1 11 12 13];
B = [1 2 3 1 4; 1 4 1 6 6; 1 5 6 7 5; 1 8 9 10 1; 1 6 2 6 3];
% [q r] = qrs(B)
% [q r] = qrg(B)
[q,r,p] = qrd(B)
% 利用命令行输入语句：
% A = input('请输入原始矩阵：')

% 基于Gran-Schmidt正交化方法，适用于满秩的方阵（非奇异矩阵）
function [Q, R] = qrs(A)
[m, n] = size(A);
p = zeros(m, n);
Q = zeros(m, n);
R = zeros(m, n);
for i = 1:n
    if i == 1
        p(:, i) = A(:, i);
        Q(:, i) = p(:, i)/norm(p(:, i), 2);
        R(1, 1) = norm(p(:, i), 2);
    else
        temp = zeros(n, 1);
        for j = 1:i-1
            R(j, i) = (A(:, i)' * p(:, j)) / (norm(p(:, j), 2) * norm(p(:, j), 2));
            temp = temp +  R(j, i) * p(:, j);
            R(j, i) = R(j, i) .* R(j, j);
        end
        p(:, i) = A(:, i) - temp;
        Q(:, i) = p(:, i) / norm(p(:, i), 2);
        R(i, i) = norm(p(:, i), 2);
    end
end
end
% 基于householder矩阵的列主元的QR分解算法。
function [Q,R,p] = qrd(A)
[m,n] = size(A);
Q = eye(m);
p = eye(n);
for j=1:n
    c(j) = A(1:m,j)' * A(1:m,j);
end
[cr,r]=max(c);
for k=1:n-1
    H_t=eye(m);
    if(cr<=0),break;end
    c([k r])=c([r k]);
    % p记录了原始矩阵的列交换信息。
    p(:,[k r])=p(:,[r k]);
    % 交换列使首列范数最大。
    A(1:m,[k r])=A(1:m,[r k]);
    H=hst(A(k:m,k));
    A(k:m,k:n)=H * A(k:m,k:n);
    H_t(k:end,k:end)=H;
    Q = H_t * Q;
    for j=k+1:n
        c(j)=c(j)-A(k,j)^2;
    end
    [cr,r]=max(c(k+1:n));
    r=r+k;
end
Q=Q';
R=A;
end
% 对x的householder矩阵构造算法。
function [H]=hst(x)
xmod=sqrt(x' * x);
alpha=-sign(x(1))*xmod;
x(1)=x(1)+alpha;
u=x/sqrt(sum(x' * x));
H=eye(length(x))-2*(u * u');
end
% 基于Givens变化的QR分解算法（只能用于方阵，且没有列主元）。
function [Q,R]=qrg(A)
[N,M]=size(A);
R=zeros(N);
T=eye(N);
B=A;
for j=1:N-1
    Tj=eye(N+1-j);T_t=eye(N);
    B=B(min(j,2):end,min(j,2):end);
    b=B(:,1);
    for i=2:N+1-j
        temp=eye(N+1-j);
        cs=sqrt(b(1)^2+b(i)^2);
        c=b(1)/cs;
        s=b(i)/cs;
        temp(1,1)=c;
        temp(i,i)=c;
        temp(i,1)=-s;
        temp(1,i)=s;
        b = temp * b;
        Tj = temp * Tj;
    end
    B = Tj * B;
    R(j,j:end)=B(1,:);
    T_t(j:end,j:end)=Tj;
    T = T_t * T;
end
R(j+1,j+1:end)=B(end,end);
Q=T';
end
