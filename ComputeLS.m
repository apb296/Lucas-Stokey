%% LS - BGP 
% This program computes the LS allocation for separable BGP case : u(c,1-n)=psi log (c) + (1-psi) log (1-n)
clc
clear all
close all
%% PARAMETERS
% Set the technology and preference parameters.

psi=.69;
beta=.9;
sSize=2;
g=[.1 .2];
pi=[.5 .5 ; .5 .5];

Para.psi=psi; % psi is the relative weight on consmumption in the utility function
Para.beta=beta; % time discount factor
Para.sSize=sSize; % This is the dimension of the markov shock
Para.g=g; % The vector g is the value of the expenditure shock in each state s
Para.pi=pi;  % The probability transition matrix for the markov shock

% This matrix is used to solve for x
A=eye(sSize);
for s=1:sSize
    for sprime=1:sSize
        if ~(s==sprime)
        A(s,sprime)=-beta*pi(s,sprime);
        else
                   A(s,sprime)=1-beta*pi(s,sprime);
        end 
    end
end
Para.A=A;
options=optimset('Display','off');

%% Compute the FB

c_FB=psi*(1-g);
n_FB=c_FB+g;
%%
% Solve the x associated with the FB allocation using the recursive
% implemntability consition
%%
LHS=psi-(1-psi).*(n_FB./(1-n_FB));
x_FB=A\LHS';
%%
% Now use the time 0 budget constraint to get the level of assets required
% to support the FB
%%
b_FB=x_FB'.*c_FB./psi;

%% Set the grid for initial b_ 
bMin=min(b_FB)*1.1;
bMax=-bMin;
bGridSize=25;
b_Grid=linspace(bMin,bMax,bGridSize);

%% Solve the allocation for a range of b_,s0=1
% Initial Guess : LSAllocation = [n0 n(s) phi]
s0=1;
LSAllocation0=[n_FB(s0) n_FB 0];

for bind=1:bGridSize
    b_=b_Grid(bind);
LSAllocation(bind,:)=fsolve(@(z) ResFOC(z,b_,s0,Para),LSAllocation0,options);
LSAllocation0=LSAllocation(bind,:);
end

%% Plots
% The red line depicts the minumum level of assets necessary to support the
% FB. The dotted black line refers to s=1 or the low goverment expenditure
% shock

%% TIME 0 labor tax
tax=@(n,s) 1-((1-psi)/(psi)).*(n-g(s))./(1-n);
figure()    
plot(b_Grid,tax(LSAllocation(:,1),s0),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$\tau(s0)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 0 labor tax')

%% TIME 1 onwards labor tax
figure()
plot(b_Grid,tax(LSAllocation(:,2),1),':k','LineWidth',2)
hold on
plot(b_Grid,tax(LSAllocation(:,3),2),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$\tau(s)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
legend('g_l','g_h')
title('Time 1 - labor tax')

%% Multiplier - Phi
figure()
plot(b_Grid,LSAllocation(:,4),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$\phi$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Implementability Multiplier')

%% TIME 0 Consumption and Labor policies
figure()    
plot(b_Grid,LSAllocation(:,1)-g(s0),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$c(s0)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 0 consumption')

figure()
plot(b_Grid,LSAllocation(:,1),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$n(s0)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 0 labor supply')


figure()
uc =@(n,s) psi./(n-g(s));
plot(b_Grid,uc(LSAllocation(:,1),1),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$u_c(s0)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 0 marginal utility')


%% TIME 1 Consumption and Labor policies
figure()    
plot(b_Grid,LSAllocation(:,2)-g(1),':k','LineWidth',2)
hold on
plot(b_Grid,LSAllocation(:,3)-g(2),'k','LineWidth',2)
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$c(s)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
legend('g_l','g_h')

title('Time 1 consumption')

figure()    
plot(b_Grid,LSAllocation(:,2),':k','LineWidth',2)
hold on
plot(b_Grid,LSAllocation(:,3),'k','LineWidth',2)

xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$n(s)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 1 labor supply')


figure()    
uc =@(n,s) psi./(n-g(s));
plot(b_Grid,uc(LSAllocation(:,2),1),':k','LineWidth',2)
hold on
plot(b_Grid,uc(LSAllocation(:,3),2),'k','LineWidth',2)

xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$u_c(s)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 1 marginal utility')
%% Time 1 Assets     
for bind=1:bGridSize
    S=sSize;
% Retrive the solution
n0=LSAllocation(bind,1);
n=LSAllocation(bind,2:2+S-1);
phi=LSAllocation(bind,end);
c0=n0-g(s0);
l0=1-n0;
uc0=psi/c0;
ul0=(1-psi)/l0;
ucc0=-psi/c0^2;
ull0=-(1-psi)/l0^2;
c=n-g;
l=1-n;
uc=psi./c;
ul=(1-psi)./l;
ucc=-psi./(c.^2);
ull=-(1-psi)./(l.^2);
% compute x from the time -1 implemntability
LHS=uc.*c-ul.*n;
x=A\LHS';
b(bind,:)=x'./uc;
ArrowSec0(bind,:)=beta.*pi(s0,:).*uc/uc0;
ArrowSec_l(bind,:)=beta.*pi(s0,:).*uc/uc(1);
ArrowSec_2(bind,:)=beta.*pi(s0,:).*uc/uc(2);
Q0(bind,:)=sum(ArrowSec0(bind,:));
Q1(bind,:)=[sum(ArrowSec_l(bind,:)) sum(ArrowSec_2(bind,:)) ];

end
figure()
plot(b_Grid,b(:,1),':k','LineWidth',2)
hold on
plot(b_Grid,b(:,2),'k','LineWidth',2)

xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$b(s)$','Interpreter','Latex','FontSize',14)
vline([min(b_FB)],':r')
title('Time 1 assets')

%% Bond Prices
    figure()
    plot(b_Grid,Q0,'b','LineWidth',2)
    hold on
    plot(b_Grid,Q1(:,1),':k','LineWidth',2)
    hold on
    plot(b_Grid,Q1(:,2),'k','LineWidth',2)
    legend('$Q_0(s0)$','$Q_1(s_l)$','$Q_1(s_h)$')
    h = legend;
set(h, 'interpreter', 'latex')

    title('One period bonds')
    xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
    ylabel('$Q(s)$','Interpreter','Latex','FontSize',14)
%%
% The above figure plots bond prices as a function of initial assets. The
% blue line is the bond price in time 0 where s0=s_l. The black lines
% (dotted= s_l) are one period maturity bonds in time 1 conditional on the state s 

%% MRS - rho(s_1|s0)
figure()
plot(b_Grid,ArrowSec0(:,1)./pi(s0,1),':k','LineWidth',2)
hold on
plot(b_Grid,ArrowSec0(:,2)/pi(s0,2),'k','LineWidth',2)
legend('$\rho_0(s_l| s0)$','$\rho(s_h|s0)$')
h = legend;
set(h, 'interpreter', 'latex')
title('Pricing Kernel - time 0')
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$\rho(s_1|s_0)$','Interpreter','Latex','FontSize',14)
%%
% The above figure plots the pricing kernel at time 0 as a function of
% initial assets.  The dotted (solid) line plots MRS[sl|s0] (MRS[sh|s0) 

%% Pricing Kernel - rho(s'|s)
figure()
plot(b_Grid,ArrowSec_l(:,2)./pi(1,2),':k','LineWidth',2)
hold on
plot(b_Grid,ArrowSec_l(:,1)./pi(1,1),'b','LineWidth',2)
hold on
plot(b_Grid,ArrowSec_2(:,1)./pi(2,1),'k','LineWidth',2)
legend('$\rho_1(s_h| s_l)$','$\rho(s|s)$','$\rho(s_l|s_h)$')
h = legend;
set(h, 'interpreter', 'latex')
title('Pricing Kernel - time 1')
xlabel('$b_{-1}$','Interpreter','Latex','FontSize',14)
ylabel('$\rho(s_2|s_1)$','Interpreter','Latex','FontSize',14)

%%
% The above figure plots the pricing kernel at time 1 as a function of
% initial assets.  The dotted (solid) line plots MRS[sh|sl] (MRS[sl|sh). Note that if the state
% remains same we have MRS=beta


