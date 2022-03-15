function [pF,sOutF] = psus(func,d,t_star,n,p,outdist,inpdist)
%%% Inputs:
%%%     func    - single input function handle to the probabilistic code
%%%     d       - dimensionality of the input space (scalar)
%%%     t_star  - critical threshold (scalar)
%%%     n       - target number of samples per P-SuS level (scalar)
%%%     p       - target level probability (scalar)
%%%     outdist - name of the output distribution (char)
%%%     inpdist - name and parameter tuple for the input distribution 
%%%               (1-by-2 cell). If left not specified, a uniform distro on
%%%               [0,1] is used.
%%%
%%% Outputs:
%%% 	pF	  - probability of failure structure, with deterministic pF and
%%% 			information about the pF distribution - 
%%%				mean and variance under independence and perfect dependence
%%%				between levels.
%%% 	sOutF - structure containing full information about the psus run
%%%
%%% Only a normal proposal distribution is used for now.
%%% The function requires Statistics and Machine Learning Toolbox in MATLAB
%%% to be installed.

%% PREPARE PRELIMINARIES
[getMV, excdFun] = checkOut(outdist); %Check that output distribution is
                                      %available and get moment transform
                                      %and membership functions
zeroProb = false;
logcAcc = @(pTarg) rand(length(pTarg),1) < pTarg; %Acceptance function

dFunc = @(x,s) normrnd(x,s); %Construct proposal distribution

p0N = p*n;

%% PREPARE INPUTS
if ~exist('inpdist','var'), inpdist = {'unif',[0,1]}; end
dist = constructD(inpdist);

%% SAMPLE DATA
x = dist.random(n,d);

%% OBTAIN RESPONSE - Assume 2 parameter location scale for now
[p1,p2] = func(x); %Get output distribution parameters
D = getMV(p1,p2);
y = D(:,1);
u = D(:,2);

%% RANK RESPONSES
[xS,pS,yS,uS] = psort(outdist,[p1,p2],y,u,x,50);

%% Output
inpPar = struct('func',func,'outdist',outdist,'dim',d,'t_star',t_star,...
    'p_0',p,'N',n);
sOut = struct('x',[],'y',[],'u',[],'pars',[],'indF',[],'indFi',[],'t_i',[],...
	'p_star',[],'p_ij',[],'N_i',[],'N_C',[],'p_Ci',[],'N_F',[],'p_Fi',[]);

%% Set loop
L = 1; %Conditional Level
nGen = n; %Samples at uncond level
nPt = n; %To correctly compute pF if no conditional levels are needed

%% Run loop
while true %nF < n*p
    %% Record failure
    pExcdF = excdFun(pS,t_star); %Probability of exceeding threshold
    indF = logcAcc(pExcdF);
    nF = sum(indF);
    
    %% Compute moments of counting distro
    mn(L,1) = sum(pExcdF); %Mean - can be used in both Poisson and Gaussian approx.
    vr(L,1) = sum( pExcdF.*(1-pExcdF) ); %Variance - for Gaussian approx.;
    
    %% Compute scaling constant - N_F
    CF(L,1) = min( mn(L)/3/sqrt(vr(L)), (nGen(L)-mn(L))/3/sqrt(vr(L)) );
    
    %% Fill in new data
    sOut = fillOut(sOut,L,xS,pS,yS,uS,indF,[],[],pExcdF,[],...
                   nGen(L),[],nF,[],[],[],mn(L),vr(L),CF(L));
                                       
    if nF > n*p, break; end
	
	%% CALCULATE LEVEL
	level = yS(p0N);
    
	%% Next level probabilities
	pInFi = excdFun(pS,level); %Probability of exceedance
	
	%% Counting distribution moments
	mn(L,1) = sum(pInFi); %Mean - can be used in both Poisson and Gaussian approx.
	vr(L,1) = sum( pInFi.*(1-pInFi) ); %Variance - for Gaussian approx.;
    
    %% Compute scaling constant - N_C
    C(L,1) = min( mn(L)/3/sqrt(vr(L)), (nGen(L)-mn(L))/3/sqrt(vr(L)) );
    
	%% Choose seeds
	indFi = logcAcc(pInFi);
	indFi = find(indFi,floor(mn(L)),'first');
	
	nPt(L,1) = length(indFi); 
	
	seeds = xS(indFi,:);
	
	%% Fill in update
	sOut = fillOut(sOut,L,[],[],[],[],[],indFi,level,[],pInFi,[],nPt(L),...
                          [],mn(L),vr(L),C(L),[],[],[]);

	%% Use MMA to populate the conditional level
    [condSamp, p1, p2] = pmma(dFunc,excdFun,dist,func,d,seeds,...
                              pS(indFi,1),pS(indFi,2),level,n,nPt(L));
    
    p1 = p1(:);
    p2 = p2(:);
    
    L = L+1;
	
    %% Restructuring 'seeds' for sorting
    rows = numel(condSamp)/d; %Find number of rows for reshaping
    condSamp = reshape(permute(condSamp,[1,3,2]),rows,d); %Reshape samples appropriately
    
	nGen(L,1) = size(condSamp,1);
    
    D = getMV(p1,p2);
    y = D(:,1);
    u = D(:,2);

    [xS,pS,yS,uS] = psort(outdist,[p1,p2],y,u,condSamp,50);
    
    %% Timing
    if L == 17
        warning('Probability of failure is zero to machine precision\n%s','Exiting...');
        zeroProb = true;
        break;
    end
end

%% Calculate probability of failure
L = L-1; %Adjust to correct number of levels

if ~zeroProb
    % Independence among Bernoulli's
    pF.pF = prod( nPt(1:L)./nGen(1:L) ) * nF/nGen(L+1);
    pF.mean = prod(mn./nGen);
    pF.var = varindepprod(mn,vr)/prod(nGen.^2);
    
    % Maximal allowable dependence
    C = [sOut(:).p_Ci];
    C = [[C(:).C]';sOut(end).p_Fi.C];
    pF.Cvar = varindepprod(mn,C.^2.*vr)/prod(nGen.^2);
    
else
    pF.pF = 0;
    pF.mean = 0;
    pF.var = inf;
end

sOutF.Results = sOut;
sOutF.pF = pF;
sOutF.Inputs = inpPar;
end

function out = fillOut(out,L,x,pars,y,u,indF,indFi,lev,pExF,pExFi,nGen,nC,nF,...
                       mn,vr,C,mnF,vrF,CF)
  
if ~isempty(x), out(L).x = x; end
if ~isempty(pars), out(L).pars = pars; end
if ~isempty(y), out(L).y = y; end
if ~isempty(u), out(L).u = u; end

if ~isempty(indF), out(L).indF = indF; end
if ~isempty(indFi), out(L).indFi = indFi; end

if ~isempty(lev), out(L).t_i = lev; end
if ~isempty(pExF), out(L).p_star = pExF; end
if ~isempty(pExFi), out(L).p_ij = pExFi; end

if ~isempty(nGen), out(L).N_i = nGen; end
if ~isempty(nC), out(L).N_C = nC; end
if ~isempty(nF), out(L).N_F = nF; end

if ~isempty(mn), out(L).p_Ci.mn = mn; end
if ~isempty(vr), out(L).p_Ci.vr = vr; end
if ~isempty(C), out(L).p_Ci.C = C; end

if ~isempty(mnF), out(L).p_Fi.mn = mnF; end
if ~isempty(vrF), out(L).p_Fi.vr = vrF; end
if ~isempty(CF), out(L).p_Fi.C = CF; end

end

function [getMV, excdFun] = checkOut(outdist)
outDists = {'Normal'; 'Uniform'; 'Logistic'; 'Laplace'; 'Student'};

name = outDists{strncmpi(outdist,outDists,4)};

switch name
    case 'Normal'
        getMV = @(m,s) [m,s.^2];
        excdFun = @(p,lev) normcdf( (p(:,1)-lev) ./ p(:,2) );
    case 'Uniform'
        getMV = @(a,b) [(a+b)/2, (b-a).^2/12];
        excdFun = @(p,lev) 1-unifcdf(lev,p(:,1),p(:,2));
    case 'Logistic'
        getMV = @(m,s) [m,(s*pi).^2/3];
        excdFun = @(p,lev) 1-cdf('logistic',lev,p(:,1),p(:,2));
    case 'Laplace'
        getMV = @(m,s) [m,2*s.^2];
        excdFun = @(p,lev) 1-lapcdf(lev,p(:,1),p(:,2));
    case 'Student'
        getMV = @(m,nu) [m,nu./(nu-2)];
        excdFun = @(p,lev) 1-tcdf(lev-p(:,1),p(:,2));
    otherwise
        error('Unrecognized output distribution.')
end

end

function pd = constructD(distInfo)
continuousDist = {'Beta'; 'BirnbaumSaunders'; 'Burr'; 'Exponential';...
'ExtremeValue'; 'Gamma'; 'GeneralizedExtremeValue'; 'GeneralizedPareto';...
'HalfNormal'; 'InverseGaussian'; 'Logistic'; 'Loglogistic'; 'Lognormal';...
'Nakagami'; 'Normal'; 'Rayleigh'; 'Rician'; 'Stable'; 'Triangular';...
'Uniform'; 'Weibull'};

name = continuousDist{strncmpi(distInfo{1},continuousDist,4)};

pars = distInfo{2};
switch length(pars)
    case 1
        pd = makedist(name, pars);
    case 2
        pd = makedist(name, pars(1), pars(2));
    case 3
        pd = makedist(name, pars(1), pars(2), pars(3));
    case 4
        pd = makedist(name, pars(1), pars(2), pars(3), pars(4));
    otherwise
        error('No supported distribution has %d parameters.');
end

end

function p = lapcdf(lev,m,s)

if lev <= m
    p = 0.5*exp( (lev-m)./s );
else
    p = 1 - 0.5*exp( (m-lev)./s );
end

end