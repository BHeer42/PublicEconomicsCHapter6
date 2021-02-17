function Ch6_social_security2()
% 
%	Chapter 6.3.3 PAYG System
%   Social security and elastic labor
%   Author: Burkhard Heer
%   Last Change: June 18, 2019


clc;
close all;


% Parameterization 
%

lss=0.3;		% labor supply steady state
alpha=0.36;		% production elasticity of capital
beta0=0.40;		% discount factor 
n=0.1;			% population growth rate
nt=5;          % number of transition periods
% disutility parameters: see Example 2
nu1=3.33;
% case 1: no pension
tau=0;      % labor income tax in steady state
tau0=0;     % labor income tax during transition in period t
tau1=0;     % labor income tax during transition in period t+1

par.lss = lss;
par.alpha = alpha;
par.beta0 = beta0;
par.n = n;
par.nu1 = nu1;
par.tau = tau;
par.tau0 = tau0;
par.tau1 = tau1;

par.nt = nt;
par.display = 0;    % if 1, the solution from findkt() is graphed

% computation of the steady state and calibration of nu0
kss=( beta0/(1+beta0) * (1-alpha)/ (1+n)  )^(1/(1-alpha))*lss;
wss=wage(kss,lss,par);
rss=interest(kss,lss,par);
c1ss=1/(1+beta0)*wss*lss;
nu0=wss/ ( lss^nu1 * c1ss); 
par.nu0 = nu0;
xinitial = [kss, lss]
disp('ksteady(kss|lss)=0? ');
ksteady(xinitial,par) 
pause;


kss;
lss;
tau;
par.kss = kss;
par.lss = lss;
c1ss;
yss=kss^(alpha)*lss^(1-alpha)
c2ss=beta0*(1+rss)*c1ss
utilss=log(c1ss)+beta0*log(c2ss)-nu0*lss^(1+nu1)/(1+nu1)
pause;




par.k1 = kss;
par.l1 = lss;
disp('kdyn in steady state =0?');
kdyn(xinitial,par)
pause;

% steady state case 2: pension

par.tau = 0.3;
par.tau0 = 0.3;
par.tau1 = 0.3;
tau = par.tau;

% steady state with pensions
xss1 = fsolve(@(x)ksteady(x,par),xinitial);
disp('ksteady(x)=0 for \tau = 30\%?');
ksteady(xss1,par);
kssp=xss1(1)
lssp=xss1(2)
par.k1=kssp; 
par.l1=lssp; 
par.kssp = kssp;
par.lssp = lssp;
wssp=wage(kssp,lssp,par);
rssp=interest(kssp,lssp,par);
c1ssp=1/(1+beta0)*((1-tau)*wssp*lssp+(1+n)/(1+rssp)*tau*wssp*lssp)
tau

yssp = kssp^(alpha)*lssp^(1-alpha);
c2ssp=beta0*(1+rssp)*c1ssp;
utilssp=log(c1ssp)+beta0*log(c2ssp)-nu0*lssp^(1+nu1)/(1+nu1);

cec=exp( (utilssp-utilss)/(1+beta0) )-1;
cec
pause;

par.tau0 = 0.3;
par.tau1 = 0.3;
disp('kdyn in steady state =0 for tau=30\%?');
kdyn(xss1,par)
pause;


% computation of the dynamics
kt=zeros(nt+1,1);
labort=kt;
periods=linspace(0,nt,nt+1);

kinit=0.99*kss;
kt(nt+1)=kss;
labort(nt+1)=lss;
disp('computation of the solution for the pension case');
findkt(kinit,par)   
pause;

ksolution = fsolve(@(x)findkt(x,par),kinit);
ksolution
pause;
disp('value of findkt: ');
findkt(ksolution,par)
pause;

% display solution
par.display = 1;
findkt(ksolution,par);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = findkt(x,par)
% input: k_T
% output: k_0

    nt = par.nt;
% value of k in the last period; afterwards k is equal to kss
    kt = zeros(nt+1,1);
    labort = zeros(nt+1,1);
	kt(nt+1)=x;
    par1 = par;
	
% we need to find labort(nt+1)	
% the optimal labor supply follows from the first-order condition of labor 

    par1.tau0 = 0;
    par1.tau1 = 0;
    kss = par.kss;
    lss = par.lss;
	w1=wage(kss,lss,par1);
	r1=interest(kss,lss,par1);
	d1= par.tau1*w1*lss;
	
	par1.k0 = x;
	par1.r1 = r1;
    par1.d1 = d1;
    
    l0 = fsolve(@(x)findl(x,par1),lss);
    
	labort(nt+1)=l0;

    x1 = [kss, lss];
    
    for i=nt:-1:1
        i
        if i==3
            par.tau1=0;
            par.tau0=0;
        elseif i==2     % i=2: period t=1; young generation has to pay contributions, but does not receive pension in old age	
            par.tau0=0.3;
        elseif i==1;
            par.tau1 = 0.3;
        end	
        
		       
		k1=kt(i+1); 
        l1=labort(i+1);
        par1.k1 = k1;
        par1.l1 = l1;
        
        
        x0 = fsolve(@(x)kdyn(x,par1),x1);
		x1 = x0;
		k0=x0(1); 
        l0=x0(2);
		kt(i)=k0;
		labort(i)=l0;
    end
    
    % display solution
    if par.display==1
        periods=linspace(0,nt,nt+1);
        figure
        plot(periods,kt);
        xlabel('Period');
        ylabel('Capital stock k_t');
        x=horzcat(periods',kt);
        x
    end

    
	y=kt(1)-par1.kssp;
end
    
%% first-order condition labor
function [y] = findl(l,par)

    tau0 = par.tau0;
    
    nu0 = par.nu0;
    nu1 = par.nu1;
    beta0 = par.beta0;
    n = par.n;
    r1 = par.r1;
    d1 = par.d1;
    l0 = l;
    k0 = par.k0;
	w0=wage(k0,l0,par);
	y=(1-tau0)*w0-nu0*l0^nu1/(1+beta0)*( (1-tau0)*w0*l0 + (1+n)/(1+r1)*d1 );
    
end


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
    nu0 = par.nu0;
    nu1 = par.nu1;
    tau0 = par.tau0;
    tau1 = par.tau1;
    k1 = par.k1;
    l1 = par.l1;
   
    
	k0 = x(1);
	l0 = x(2);
	y=zeros(2,1);
	r1=interest(k1,l1,par);
	w0=wage(k0,l0,par);
	w1=wage(k1,l1,par);
	d1=tau1*w1*l1;
    
		
	y(1)=k1*(1+n)-beta0/(1+beta0)*(1-tau0)*w0*l0+1/(1+beta0)*(1+n)/(1+r1)*d1;
	y(2)=(1-tau0)*w0-nu0*l0^nu1/(1+beta0)*( (1-tau0)*w0*l0 + (1+n)/(1+r1)*d1 );
   
end

%% ksteady: steady state
function [y] = ksteady(x,par)
    
    tau = par.tau;
    n = par.n;
    beta0 = par.beta0;
    nu0 = par.nu0;
    nu1 = par.nu1;
    
    
	kp=x(1);
	lp=x(2);
	r1=interest(kp,lp,par);
	w0=wage(kp,lp,par);
	w1=wage(kp,lp,par);
 	d1=lp*tau*w1;
	y=zeros(2,1);
	y(1)=kp*(1+n)-beta0/(1+beta0)*(1-tau)*w0*lp+1/(1+beta0)*(1+n)/(1+r1)*d1;
	y(2)=(1-tau)*w0-nu0*lp^nu1/(1+beta0)*( (1-tau)*w0*lp + (1+n)/(1+r1)*d1 );
    
end

