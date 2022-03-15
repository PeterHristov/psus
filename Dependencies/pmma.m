function [seeds,par1,par2,pA] = pmma(propFunc,excdFunc,dist,func,d,seeds,par1,par2,level,nL,nC)
%%% Modified Metropolis algorithm for distributional data.
%%% Inputs:
%%% 	propFunc - proposal pdf
%%% 	excdFun  - exceedence function
%%% 	dist 	 - distribution for model inputs
%%% 	func     - performance function
%%% 	d        - input dimensionality
%%%		seeds    - seeds for level
%%% 	par1	 - first parameter of the distribution
%%% 	par2	 - second parameter of the distribution
%%% 	y        - mean of response
%%% 	u        - variance of response
%%% 	level    - moments of interim threshold variable
%%% 	nL       - minimum number of samples per level (keep seeds)
%%% 	nC       - number of seeds
%%%
%%% Outputs:
%%% 	seeds - seeds for the next level of ps
%%% 	par1  - first parameter of the distribution of the seeds
%%% 	par2  - second parameter of the distribution of the seeds
%%% 	pA	  - probability of acceptance
%%%
%%%	The function requires Statistics and Machine Learning Toolbox in MATLAB
%%% to be installed.


nS = ceil( (nL-nC)/nC ); %Number of states per chain when keeping seeds

s = repmat(std(seeds),size(seeds,1),1); %Proposal distribution std
propFunc = @(x) propFunc(x,s);

pA = zeros(nS,d);%Initializing array for the probability of acceptance

for k = 1:nS
    %% Random walk
    urand = rand(nC,d); %Uniform numbers to check acceptance
    
    pstar = propFunc(seeds(:,:,k)); %Step for random walk
    
    r = dist.pdf(pstar)./dist.pdf(seeds(:,:,k)); %Uniform pdfs 
    accept = urand < r; %Acceptance criterion
    pA(k,:) = mean(accept); %Probability of acceptance
    zeta = seeds(:,:,k); %Copy to zeta
    zeta(accept) = pstar(accept); %Replace the old samples with the accepted ones where appropriate
    
    %% Domain acceptance
    par1(:,k+1) = par1(:,k);
    par2(:,k+1) = par2(:,k);

    [par1(any(accept,2),k+1), par2(any(accept,2),k+1)] = func( zeta(any(accept,2),:) ); %Evaluate the objective function at new points
    
    seeds(:,:,k+1) = seeds(:,:,k); %Populate the next state of the chain with the previous one by default
    
    pInFi = excdFunc([par1(:,k+1),par2(:,k+1)],level); %Prob. y is in Fi
    
    urandF = rand(nC,1); %Uniform RV's to check acceptance in F
    inFi = urandF < pInFi;

	%% Selection
    seeds(inFi,:,k+1) = zeta(inFi,:); %Replace seeds that have passed with zeta...
    par1(~inFi, k+1) = par1(~inFi, k); %... replace parameter values too...
    par2(~inFi, k+1) = par2(~inFi, k); %...
end
end
