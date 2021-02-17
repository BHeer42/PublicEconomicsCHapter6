@ ----------------------------  Ch6_social_security2.g -------------------------------

computes the Example 2 in Chapter 2: Social security and elastic labor

author: Burkhard Heer

this version: June 18, 2019s


------------------------------------------------------------------------------ @

new; clear all; cls;
library pgraph, user;
MachEps=machEpsilon;
#include toolbox.src;



alpha=0.36;		// production elasticity of capital
beta0=0.40;		// discount factor 
nu1=3.33;		// inverse of Frisch elasticity
n=0.1;			// population growth rate 

// case 1: no pension
tau=0;

// compute a grid for the solution for various k
lss=0.3;

kss=( beta0/(1+beta0) * (1-alpha)/ (1+n)  )^(1/(1-alpha))*lss;
wss=wage(kss,lss);
rss=rate(kss,lss);
c1ss=1/(1+beta0)*wss*lss;
nu0=wss/ ( lss^nu1 * c1ss); 

k1=kss; l1=lss; k0=0; l0=0; w0=0; r1=rate(kss,lss); w1=wage(kss,lss); d1=0; // initialization because they get values assigned in procedures
"ksteady(kss|lss)=0?";
ksteady(kss|lss); wait;

"kss: " kss;
"lss: " lss;
"tau: " tau;
"c1: " c1ss;
"y: " kss^(alpha)*lss^(1-alpha);
c2ss=beta0*(1+rss)*c1ss;
utilss=ln(c1ss)+beta0*ln(c2ss)-nu0*lss^(1+nu1)/(1+nu1);
wait;

tau0=0;
tau1=0;
"kdyn in steady state =0?";
kdyn(kss|lss);
wait;



// case j=2: pension

tau=0.3;
// steady state with pensions
x0=kss|lss;
{x,crit}=FixVMN1(x0,&ksteady);
"ksteady(x)=0?";
ksteady(x);
"steady state values k,l: " x';
wait;
kssp=x[1];
lssp=x[2];
k1=kssp; l1=lssp; 
"kssp: " kssp;
"lssp: " lssp;
wssp=wage(kssp,lssp);
rssp=rate(kssp,lssp);
c1ssp=1/(1+beta0)*((1-tau)*wssp*lssp+(1+n)/(1+rssp)*tau*wssp*lssp);
"tau: " tau;
"c1p: " c1ssp;
"y: " kssp^(alpha)*lssp^(1-alpha);
c2ssp=beta0*(1+rssp)*c1ssp;
utilssp=ln(c1ssp)+beta0*ln(c2ssp)-nu0*lssp^(1+nu1)/(1+nu1);

cec=exp( (utilssp-utilss)/(1+beta0) )-1;
"cec: " cec;
wait;

tau0=0.3;
tau1=0.3;
"kdyn in steady state =0?";
kdyn(kssp|lssp);
wait;


// computation of the dynamics
nt=5;		// number of transition periods
kt=zeros(nt+1,1);
labort=kt;
periods=seqa(0,1,nt+1);

// initialization
tau0=0;		// tau in period t 
tau1=0;		// tau in period t+1  

j=1;
kinit=0.99*kss;
kt[nt+1]=kss;
labort[nt+1]=lss;
"computation of the solution for the pension case";
findkt(kinit);
wait;
{ksolution,crit}=FixVMN1(kinit,&findkt);
"solution k: " ksolution; wait;
"value of findkt: ";
findkt(ksolution);
wait;

struct plotControl myPlot;
myplot=plotGetDefaults("xy");

plotSetLegend(&myplot,"off");
plotSetGrid(&myplot,"off");
plotsetxlabel(&myplot,"period","verdana",20,"black");
plotsetylabel(&myplot,"capital stock","verdana",20,"black");
plotXY(myplot,periods,kt);


// procedure to find labor supply l0 for given k0,k1,l1

proc findl(l0);
	local y;
	w0=wage(k0,l0);
	y=(1-tau)*w0-nu0*l0^nu1/(1+beta0)*( (1-tau)*w0*l0 + (1+n)/(1+r1)*d1 );
	retp(y);
endp;


/* procedure that computes k_t+1 given k_t*/
// input: k0 - next period capital stock 
// output: the value of the equilibrium condition

proc kdyn(x0);
	local crit,y;
	k0=x0[1];
	l0=x0[2];
	y=zeros(2,1);
	r1=rate(k1,l1);
	w0=wage(k0,l0);
	w1=wage(k1,l1);
	
	d1=tau1*w1*l1;
		
	y[1]=k1*(1+n)-beta0/(1+beta0)*(1-tau0)*w0*l0+1/(1+beta0)*(1+n)/(1+r1)*d1;
	y[2]=(1-tau0)*w0-nu0*l0^nu1/(1+beta0)*( (1-tau0)*w0*l0 + (1+n)/(1+r1)*d1 );
	retp(y);
endp;



// input: k_T
// output: k_0
proc findkt(x);
	local y,crit,i,x0;

// value of k in the last period; afterwards k is equal to kss
	kt[nt+1]=x;
	
// we need to find labort[nt+1,j]	
// the optimal labor supply follows from the first-order condition of labor 

	j=1;		// no pensions in new steady state
	tau0=0;
	tau1=0;
	w1=wage(kss,lss);
	r1=rate(kss,lss);
	d1=0;
	
	k0=x;
	
	{l0,crit}=FixVMN1(lss,&findl);
	labort[nt+1]=l0;


	i=nt+1;
	do until i==1;
		i=i-1; i;
		
		// i=1: Period 0 with steady state, tau=30% for young and old
		// i=2: Period 1, young generation still has to finance old agents (tau0=30% so that d_1=tau_1 w_1 l), but will not receive
		//		a pension in old age in period 2 (i=3) with tau1=0% so that d_2=0
		if i==2; tau1=0; endif;	// first generation still has to pay contributions for pension of the old
		if i==3; tau0=0; endif;	// first generation which does not have to pay contributions
		k0=kt[i,1];
		
		k1=kt[i+1]; l1=labort[i+1];
		{x0,crit}=FixVMN1(k1|l1,&kdyn);
		if crit[1]/=0; "no convergence of FixVMN1"; wait; endif;
		k0=x0[1]; l0=x0[2];
		kt[i]=k0;
		labort[i]=l0;
		"i~k0~l0~k1~l1: " i~k0~l0~k1~l1; 
	endo;
	retp(kt[1]-kssp);	// how close it the value kt[1] to its initial value kssp?
endp;

// computes the wage
proc wage(k,l);
	retp( (1-alpha)* k^(alpha) * l^(-alpha));
endp;

// computes the interest rate
proc rate(k,l);
	retp( alpha* k^(alpha-1) * l^(1-alpha) );
endp;

// computes steady state for pension case
// kp -- steady state value of capital stock with pensions
// lp -- steady state value of labor with pensions
proc ksteady(x);
	local kp,lp,y,d1;
	kp=x[1];
	lp=x[2];
	r1=rate(kp,lp);
	w0=wage(kp,lp);
	w1=wage(kp,lp);
	d1=lp*tau*w1;
	y=zeros(2,1);
	y[1]=kp*(1+n)-beta0/(1+beta0)*(1-tau)*w0*lp+1/(1+beta0)*(1+n)/(1+r1)*d1;
	y[2]=(1-tau)*w0-nu0*lp^nu1/(1+beta0)*( (1-tau)*w0*lp + (1+n)/(1+r1)*d1 );
	retp(y);
endp;

proc(0)=set_xy();

    graphset;
    _pmcolor = 0|0|0|0|0|0|0|0|15; 
    _pframe = { 1, 1};
    _pcolor = { 12, 2, 1, 0}; 
    _pltype = { 6, 3, 5, 6 };
    _pcross=0;
    //_pgrid = { 2, 2}; 
    _ptitlht = 0.25;
    _plwidth = 10;
    _paxht=0.25;
    _pnumht=0.20;
    _pdate=0;
    fonts("simplex complex microb simgrma");
    pause(0.5);

retp();
endp;
