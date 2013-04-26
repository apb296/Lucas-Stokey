function [ res ] = ResFOC(LSAllocation,b_,s0,Para)
% This function computes the residual for a given LS allocation :
% n0,n(s),phi

% Retrive the paramters from the para struct.
psi=Para.psi;
beta=Para.beta;
g=Para.g;
pi=Para.pi;
S=Para.sSize;
A=Para.A;

% Retrive the guessed solution
n0=LSAllocation(1);
n=LSAllocation(2:2+S-1);
phi=LSAllocation(end);

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


if and(((n0>g(s0)) && (n0<1)),(n>g) & (n<1))
   
% Time 0 FOC
res(1)=(uc0-ul0)-phi*( uc0+ c0*ucc0-ul0-n0*(-ull0)-ucc0*b_ );
% Time 1 FOC
res(2:2+S-1)=(uc-ul)-phi*(uc+(n-g).*ucc-ul-n.*(-ull));
% Time 0 budget
res(2+S)=uc0*c0+ beta*pi(s0,:)*x-ul0*n0-uc0*b_;
    else
        res= [ abs(n0-g(s0)) abs(n-g) phi]*100+10;
end

end

