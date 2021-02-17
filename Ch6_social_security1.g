@ ----------------------------  Ch6_social_security1.g -------------------------------

computes the Example 1 in Chapter 6: Social security and inelastic labor

also: transition from steady state with tau=0.3 to tau=0

labor supply = 0.30

author: Burkhard Heer

this version: May 28, 2019

------------------------------------------------------------------------------ @

new; clear all; cls;
library pgraph, user;
MachEps=machEpsilon;
#include toolbox.src;


_plotstyle=1;       // 1 -- newer versions of Gauss


//
//	Step 1: Parameterization
//

ls=0.3;			// labor supply
alpha=0.36;		// production elasticity of capital
beta0=0.40;		// discount factor 
n=0.1;			// population growth rate 

// disutility parameters: see Example 2
nu0=257.15;
nu1=3.33;

// case 1: no pension
tau=0;

//
// Step 2: Computation of initial value for non-linear eqs system
//


// compute a grid for the solution for various k
kmin=0.001;
kmax=5.0;
nk=1000;
k=seqa(kmin,(kmax-kmin)/(nk-1),nk);
yk=zeros(nk,1);
i=0;
do until i==nk;
	i=i+1;
	yk[i]=ksteady(k[i]);
	i~yk[i];
endo;

// Initial guess
i0=minindc(abs(yk)); // locates the index of the absolute minimum in yk
kinit=k[i0];	
"Kinit: " kinit;
"i0 : " i0;
wait;

//
// Step 3: Computation of the final steady state with tau=0%
//

{kss,crit}=FixVMN1(kinit,&ksteady);
if crit[1]/=0; "no convergence of FixVMN1"; wait; endif;

"retcode FixVMN1: " crit;
wait;
"ksteady(kss)";
ksteady(kss);
"kss: " kss;
"y: " kss^(alpha)*ls^(1-alpha);
c1ss=1/(1+beta0)*(wage(kss,ls)*ls+(n-rate(kss,ls))/(1+rate(kss,ls))*tau*wage(kss,ls)*ls);
c2ss=beta0*c1ss*(1+rate(kss,ls));
utilss=ln(c1ss)+beta0*ln(c2ss)-nu0*ls^(1+nu1)/(1+nu1);
wait;

//
// Step 4: computation of steady state with pensions
//

tau=0.3;
// compute a grid for the solution for various k
kmin=0.001;
kmax=5.0;
nk=1000;
k=seqa(kmin,(kmax-kmin)/(nk-1),nk);
yk=zeros(nk,1);
i=0;
do until i==nk;
	i=i+1;
	yk[i]=ksteady(k[i]);
	i~yk[i];
endo;

// Initial guess
i0=minindc(abs(yk)); // locates the index of the absolute minimum in yk
kinit=k[i0];	

{kssd,crit}=FixVMN1(kinit,&ksteady);
if crit[1]/=0; "no convergence of FixVMN1"; wait; endif;

"retcode FixVMN1: " crit;
wait;
"ksteady(kssd)";
ksteady(kssd);
"kssd: " kssd;
"y: " kssd^(alpha)*ls^(1-alpha);
c1ssd=1/(1+beta0)*(wage(kssd,ls)*ls+(n-rate(kssd,ls))/(1+rate(kssd,ls))*tau*wage(kssd,ls)*ls);
c2ssd=beta0*c1ssd*(1+rate(kssd,ls));
utilssd=ln(c1ssd)+beta0*ln(c2ssd)-nu0*ls^(1+nu1)/(1+nu1);
wait;
cec=exp( (utilssd-utilss)/(1+beta0) )-1;
"cec: " cec;

//
// Step 4: computation of the dynamics
//

nt=20;		// number of transition periods
kt=zeros(nt+1,1);
utilt=zeros(nt,1);
periods=seqa(0,1,21);
periods1=seqa(0,1,20);
kt[1,1]=kssd;
utilt[1,1]=utilssd;
cect=zeros(nt,1);
cect[1]=exp( (utilt[1]-utilss)/(1+beta0) )-1;
i=0;
tau0=0.3;	// tax rate at young age
tau1=0.3;	// tax rate at old age
do until i==nt;
	i=i+1;
	// i=1: Period 0 with steady state, tau=30% for young and old
	// i=2: Period 1, young generation still has to finance old agents (tau0=30% so that d_1=tau_1 w_1 l), but will not receive
	//		a pension in old age in period 2 (i=3) with tau1=0% so that d_2=0
	if i==2; tau1=0; endif;	// first generation still has to pay contributions for pension of the old
	if i==3; tau0=0; endif;	// first generation which does not have to pay contributions
	k0=kt[i,1];
	{k1,crit}=FixVMN1(k0,&kdyn);
        
	d0=tau0*wage(k0,ls)*ls;
	d1=tau1*wage(k1,ls)*ls;
    c1 = 1/(1+beta0)* (wage(k0,ls)*ls -d0 + d1*(1+n)/(1+rate(k1,ls)));
	c2=beta0*c1*(1+rate(k1,ls));
	utilt[i,1]=ln(c1)+beta0*ln(c2)-nu0*ls^(1+nu1)/(1+nu1);
	cect[i,1]=exp( (utilt[i]-utilss)/(1+beta0) )-1;
	if crit[1]/=0; "no convergence of FixVMN1"; wait; endif;
    kt[i+1,1]=k1;
endo;


if _plotstyle==0;
	set_xy;
	k_font= "Capital K]t[";
	delta_font= "Welfare \204 \68 \201]t[";
	xlabel("Period t");
	ylabel(k_font);
	XY(periods,kt);

	wait;

	ylabel(delta_font);
	XY(periods1,cect*100);
	
else;
	graphset;
    struct PlotControl myPlot;
    myplot=plotGetDefaults("xy");	
	plotSetLegend(&myPlot,"off");
	plotSetGrid(&myPlot,"off");
	plotSetLineThickness(&myPlot,5);
    plotSetXLabel(&myplot,"Period t", "verdana", 20, "black");
    plotSetYLabel(&myplot,"Capital Stock");
    plotXY(myplot,periods,kt);
     wait;
	 
    plotSetYLabel(&myplot,"Welfare");
    plotXY(myplot,periods1,cect*100);	
endif;


@ --------------------------   Procedures ---------------------------------- @


/* procedure that computes the steady state*/
// input: k0 - steady state capital stock 
// output: the value of the equilibrium condition
proc ksteady(k0);
	local y,d;
	d=tau*wage(k0,ls)*ls;
	y=k0*(1+n)-wage(k0,ls)*ls*beta0/(1+beta0)+1/(1+beta0)*(1+beta0+beta0*rate(k0,ls)+n)/(1+rate(k0,ls)) * d;
	retp(y);
endp;


/* procedure that computes k_t+1 given k_t*/
// input: k1 - next period capital stock 
// output: the value of the equilibrium condition
proc kdyn(k1);
	local y, d0, d1;
//	d0=tau0*wage(k0,ls)*ls;
	d1=tau1*wage(k1,ls)*ls;
	y=k1*(1+n)-(1-tau0)*wage(k0,ls)*ls*beta0/(1+beta0)+ 1/(1+beta0)*(1+n)/(1+rate(k1,ls)) * d1;
	retp(y);
endp;


// computes the wage
proc wage(k,l);
	retp( (1-alpha)* k^(alpha) * l^(-alpha) );
endp;

// computes the interest rate
proc rate(k,l);
	retp( alpha* k^(alpha-1)*l^(1-alpha) );
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
