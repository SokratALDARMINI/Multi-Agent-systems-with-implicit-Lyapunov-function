%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB program is attached with Simulink model. it is built depending on the articles:
%              1- "Yu, Zhiyong, Shuzhen Yu, and Haijun Jiang.
%                 "Consensus of multi-agent systems with finite-time and fixed-time observation."
%                  Information Sciences 512 (2020): 909-928."
%              2- "Lopez-Ramirez, Francisco, et al. "Finite-time and fixed-time observer design:
%                  Implicit Lyapunov function approach." Automatica 87 (2018): 52-60."
%              3- "Zimenko, Konstantin, et al. "A note on delay robustness for homogeneous systems
%                  with negative degree." Automatica 79 (2017): 178-184."               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of connections between the agents
N=5; %number of agents
G1=[1 1 2 3 4];
G2=[2 5 3 4 5];
G=graph(G1,G2);
% Agent model
A=[0 1;0 0];
C=[1 1];
B=[0;1];
rank([C;C*A;C*A*A]) % Check observability
%initial values
x1_0=[1;2];
x2_0=[-1;2];
x3_0=[3;-2];
x4_0=[-2;3];
x5_0=[6;0];
on_off=0;%parameter to add noise or delete it
s=sqrt(0.1);%noise standard deviation
Td=0%delay value
%% Plot of Graph that connects between the agents
figure();
plot(G,'LineWidth',2);
title('Connections between the agents');
set(gca,'FontSize',20)
%% transformation to get chain of integrator form
TITA=[1 1;0 -1];
Ac=TITA*A/TITA;
Cc=C/TITA;
c0=Cc(1);
F=Ac(:,1);
C_tilda=[1 0];
A_tilda=[0 Ac(1,2);0 0];
%% Solution of LMI to for finite time stability of the observer
%Parameters
beta=0.1;
mui=0.4;
row=0.6;
tou=1.2;

Hv=[1 0;0 1-mui/(1+mui)];
E=@(x)(([x^(1-mui/(1+mui)) 0;0 x^(1-2*mui/(1+mui))]-x*eye(2)));

clear LMIs
P=sdpvar(2,2);
S=sdpvar(2,2);
Y=sdpvar(2,1);
alpha_2=sdpvar(1,1);% it is named as X in the report
LMIs=[P>=eye(2)*0.1];
LMIs=LMIs+[S>=eye(2)*0.1];
LMIs=LMIs+[alpha_2>=0.1];
LMIs=LMIs+[[P*A_tilda+Y*C_tilda+A_tilda'*P+C_tilda'*Y'+row*P+beta*(Hv*P+P*Hv),P;P,-row*S]<=0];
LMIs=LMIs+[(P-C_tilda'*alpha_2*C_tilda)>=0];
LMIs=LMIs+[(Hv*P+P*Hv)>=eye(2)];
LMIs=LMIs+[[tou*alpha_2,Y';Y,P]>=0];
for v=0:0.01:1
LMIs=LMIs+[(E(v)*S*E(v)-P/tou)<=0];
end
solvesdp(LMIs);
P=double(P)
S=double(S)
Y=double(Y)
alpha_2=double(alpha_2)
alpha=sqrt(alpha_2)
H=P\Y
%Check that ((E(v)*S*E(v)-P/tou)<=0) for more accurate grid
o=zeros(100001,1);
i=0;
for v=0:1e-5:1
    i=i+1;
    o(i)=max(eig((E(v)*S*E(v)-P/tou)));
end
figure();
plot([0:1e-5:1],o,'LineWidth',2);
grid on;
title('$$\max (eig\{ \Xi (v)S\Xi (v) - P/\tau \} )$$ (sampling step for $$v$$ is $$10^{-5}$$)', 'Interpreter', 'LaTeX','FontSize',20);
set(gca,'FontSize',20)
%% Solution of LMI to for consensus convergance
%Parameters
ep=0.1;
L=full(laplacian(G));%Laplacian matrix of the graph
Lambda2=[0 1 0 0 0]*eig(L);%second eigenvalue of the laplacian matrix
X=sdpvar(2,2);
clear LMIs;
LMIs=[X>=0];
LMIs=LMIs+[[2*Lambda2*B*B'-A*X-X*A',sqrt(ep)*X;sqrt(ep)*X,eye(2)]>=0];
solvesdp(LMIs);
P2=inv(double(X));
K=B'*P2; %Feedback gain
%% Matrices and parameters used to define the triggering function
ep0=0.1;
sigma=0.98;
eta=0.1;
L_kron_eye2=kron(L,eye(2));
M=(5*eye(5)-ones(5,5))/5;
M_kron_eye2=kron(M,eye(2));
eye2_kron_K=kron(eye(5),K);
L_kron_K=kron(L,K);
c1=sigma*ep/(2*norm(P2*B*B'*P2))
c2=ep0/(2*norm(P2*B*B'*P2))
%% Maximum time for finite time observation
%e0=[1;2];%for first agent
%e0=[-1;2];%for second agent
%e0=[3;-2];%for third agent
%e0=[-2;3];%for fourth agent
e0=[6;0];%for fifth agent
f=@(x)(e0'*[x^(-1) 0;0 x^(-1+mui/(1+mui))]*P*[x^(-1) 0;0 x^(-1+mui/(1+mui))]*e0-1);
Vmin=0.000005;
a=0.001;
b=20;
for i=1:1:1000
if(f(b)>0)
    a=b;b=2*b;
elseif (f(a)<0)
    b=a;a=max(a/2,Vmin);
else
    c=(a+b)/2;
    if(f(c)<0)
        b=c;
    else
        a=max(c,Vmin);
    end
end
if (abs(b-a)<0.000001)
    break;
end
end
V=b;
Te0=V^(mui/(1+mui))/(mui/(1+mui)*beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% the next parts to be implemented after running the model in Simulink
%% Observation errors
t=out.e1(:,1);
ex1=zeros(length(t),5);
ex2=zeros(length(t),5);
ex1(:,1)=out.e1(:,2);
ex1(:,2)=out.e2(:,2);
ex1(:,3)=out.e3(:,2);
ex1(:,4)=out.e4(:,2);
ex1(:,5)=out.e5(:,2);

ex2(:,1)=out.e1(:,3);
ex2(:,2)=out.e2(:,3);
ex2(:,3)=out.e3(:,3);
ex2(:,4)=out.e4(:,3);
ex2(:,5)=out.e5(:,3);
figure();
plot(t,ex1,'LineWidth',2);
title('Observation error for state x_1');
legend('e_{11}','e_{21}','e_{31}','e_{41}','e_{51}');
axis([0 30 -6 3])
set(gca,'FontSize',20)
ex1(log(abs(ex1))<-20.77)=exp(-20.77);
figure();
semilogy(t,abs(ex1),'LineWidth',2);
title('Observation error for state x_1 (logarithmic scale)');
legend('e_{12}','e_{22}','e_{32}','e_{42}','e_{52}');
set(gca,'FontSize',20)

figure();
plot(t,ex2,'LineWidth',2);
title('Observation error for state x_2');
legend('e_{11}','e_{21}','e_{31}','e_{41}','e_{51}');
axis([0 30 -6 3])
set(gca,'FontSize',20)
ex2(log(abs(ex2))<-20.77)=exp(-20.77);
figure();
semilogy(t,abs(ex2),'LineWidth',2);
title('Observation error for state x_2 (logarithmic scale)');
legend('e_{12}','e_{22}','e_{32}','e_{42}','e_{52}');
set(gca,'FontSize',20)
%% Triggering instants
trig=out.trig(:,2);
o=zeros(length(t),1);
for i=1:1:length(t)
    if(abs(trig(i))>1e2)
        o(i)=1;
    end
end
figure();
plot(t,o)
axis([0 30 -0.5 1.5])
grid on;
title(['Triggering instants: ',int2str(sum(o)),' instants'],'FontSize',20);
set(gca,'FontSize',20)
%% Consensus of agents
t=out.x(:,1);
x1=out.x(:,2:2:11);
x2=out.x(:,3:2:11);
figure();
plot(t,x1,'LineWidth',2);
title('State x_1');
legend('x_{11}','x_{21}','x_{31}','x_{41}','x_{51}');
axis([0 30 -10 50])
set(gca,'FontSize',20)

figure();
plot(t,x2,'LineWidth',2);
title('State x_2');
legend('x_{12}','x_{22}','x_{32}','x_{42}','x_{52}');
set(gca,'FontSize',20)

%% Phase plot for delay case:
t=out.e1(:,1);
ex1=zeros(length(t),5);
ex2=zeros(length(t),5);
ex1(:,1)=out.e1(:,2);
ex1(:,2)=out.e2(:,2);
ex1(:,3)=out.e3(:,2);
ex1(:,4)=out.e4(:,2);
ex1(:,5)=out.e5(:,2);

ex2(:,1)=out.e1(:,3);
ex2(:,2)=out.e2(:,3);
ex2(:,3)=out.e3(:,3);
ex2(:,4)=out.e4(:,3);
ex2(:,5)=out.e5(:,3);
figure();
plot(ex1(:,1),ex2(:,1),'LineWidth',2);
grid on;
hold on;
plot(ex1(:,2),ex2(:,2),'LineWidth',2);
plot(ex1(:,3),ex2(:,3),'LineWidth',2);
plot(ex1(:,4),ex2(:,4),'LineWidth',2);
plot(ex1(:,5),ex2(:,5),'LineWidth',2);
title('Observation errors in phase portrait');
legend('e_{1}','e_{2}','e_{3}','e_{4}','e_{5}');
xlabel('e_{i1}');
ylabel('e_{i2}');
set(gca,'FontSize',20)
