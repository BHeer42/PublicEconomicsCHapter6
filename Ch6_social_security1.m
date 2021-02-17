function Ch6_social_security1();
% 
%	Chapter 6.3.2 PAYG System
%   Author: Burkhard Heer
%   Last Change: June 5, 2019


clc;
close all;

%
% Step 1: Parameterization 
%

ls=0.3;			% labor supply
alpha=0.36;		% production elasticity of capital
beta0=0.40;		% discount factor 
n=0.1;			% population growth rate
nt=20;          % number of transition periods
% disutility parameters: see Example 2
nu0=257.15;
nu1=3.33;

%
% Step 2:   Computaion of the steady state for the
%           case 1: no pension
%
tau=0;

par.ls = ls;
par.alpha = alpha;
par.beta0 = beta0;
par.n = n;
par.nu0 = nu0
par.nu1 = nu1;
par.tau = tau;
par.nt = nt;

% compute a grid for the solution for various k
% eq. (6.19)
kmin=0.001;
kmax=5.0;
nk=1000;
k=linspace(kmin,kmax,nk);
yk=zeros(nk,1);

for i=1:1:nk
    par.i = i;
	yk(i)=ksteady(k(i),par);
	i
    yk(i)
end

clc;
yk = abs(yk);
% minimum of yk as starting value
[M, I] = min(yk)

kinitial = k(I);
kss = fsolve(@(x)ksteady(x,par),kinitial);

disp('ksteady(kss,par) ');
ksteady(kss,par)
kss
yss=kss^(alpha)*ls^(1-alpha)

c1ss=1/(1+beta0)*(wage(kss,ls,par)*ls+(n-interest(kss,ls,par))/(1+interest(kss,ls,par))*tau*wage(kss,ls,par)*ls);
c2ss=beta0*c1ss*(1+interest(kss,ls,par));
utilss=log(c1ss)+beta0*log(c2ss)-nu0*ls^(1+nu1)/(1+nu1)


%
% Step 3:   Computaion of the steady state for the
%           case 2: PAYG pension 
% 
tau=0.3;
par.tau = tau;

for i=1:1:nk
    par.i = i;
	yk(i)=ksteady(k(i),par);
	i
    yk(i)
end

clc;
yk = abs(yk);
% minimum of yk as starting value
[M, I] = min(yk)

kinitial = k(I);
kssd = fsolve(@(x)ksteady(x,par),kinitial);

pause;

disp('ksteady(kssd,par) ');
ksteady(kssd,par)
kssd
yssd=kssd^(alpha)*ls^(1-alpha)

c1ssd=1/(1+beta0)*(wage(kssd,ls,par)*ls+(n-interest(kssd,ls,par))/(1+interest(kssd,ls,par))*tau*wage(kssd,ls,par)*ls);
c2ssd=beta0*c1ssd*(1+interest(kssd,ls,par));
utilssd=log(c1ssd)+beta0*log(c2ssd)-nu0*ls^(1+nu1)/(1+nu1)
cec=exp( (utilssd-utilss)/(1+beta0) )-1
pause;

%
% Step 4:
% computation of the dynamics
%
kt=zeros(nt+1,1);
utilt=kt;
periods=linspace(0,nt,21);
kt(1)=kssd;
utilt(1)=utilssd;
cect=zeros(nt+1,1);
cect(1)=exp( (utilt(1)-utilss)/(1+beta0) )-1;


tau0=0.3;
par.tau0 = tau0;    % social security rate in period t
tau1=0.3;           % social security rate in period t+1
par.tau1 = tau1;
for i=1:1:nt
    if i==2 
        tau1=0; % first generation still has to pay contributions for pension of the old
                % but will not receive a pension itself
        par.tau1 = tau1;
    elseif i==3
        tau0=0;     % generation born in period t=2 (corresponding to i=3) is the first generation
                    % that does not pay social security contributions
        par.tau0=tau0;
    end	 
	k0=kt(i);
    par.k0 = k0;
    k1 = fsolve(@(x)kdyn(x,par),k0); 
	d0 = tau0*wage(k0,ls,par)*ls;   % social security contribution in period t=i-1
    d1 = tau1*wage(k1,ls,par)*ls;   % social security contribution in period t+1 = i
    
    c1 = 1/(1+beta0)* (wage(k0,ls,par)*ls -d0 + d1*(1+n)/(1+interest(k1,ls,par)));
	c2=beta0*c1*(1+interest(k1,ls,par));
	utilt(i+1)=log(c1)+beta0*log(c2)-nu0*ls^(1+nu1)/(1+nu1);
	cect(i+1)=exp( (utilt(i+1)-utilss)/(1+beta0) )-1;
    kt(i+1)=k1;
end

clc;


figure
plot(periods,kt);
xlabel('Period t');
ylabel('Capital K_t');
pause;

figure
plot(periods,cect*100);
ylabel('Welfare \Delta_t');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% interest rate 
function [y] = interest(k,l,par)
    alpha = par.alpha;
	y=alpha*k^(alpha-1)*l^(1-alpha);
end

%% wage 
function [w] = wage(k,l,par)
    alpha = par.alpha;
	w=(1-alpha)*k^alpha*l^(-alpha);
end


% procedure that computes k_t+1 given k_t*/
% input: k1 - next period capital stock 
% output: the value of the equilibrium condition
function [y] =  kdyn(x,par)
    beta0 = par.beta0;
	n = par.n;
    ls = par.ls;
    tau0 = par.tau0;
    tau1 = par.tau1;
    k1 = x;
    k0 = par.k0;
	d0 = tau0*wage(k0,ls,par)*ls;
    d1 = tau1*wage(k1,ls,par)*ls;
    
	y=k1*(1+n)- (1-tau0)*wage(k0,ls,par)*ls*beta0/(1+beta0)+ 1/(1+beta0)*(1+n)/(1+interest(k1,ls,par)) * d1;
end

%% ksteady: eq. (6.19) 
function [y] = ksteady(k0,par)
    
    tau = par.tau;
    ls = par.ls;
    i = par.i;
    n = par.n;
    beta0 = par.beta0;
    
	d=tau*wage(k0,ls,par)*ls;
	y = k0*(1+n)- wage(k0,ls,par)*ls*beta0/(1+beta0)+1/(1+beta0)*(1+beta0+beta0*interest(k0,ls,par)+n)/(1+interest(k0,ls,par)) * d;	 
end

