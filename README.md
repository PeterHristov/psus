# psus
This repository contains code for the probabilistic subset simulation algorithm described in "Subset simulation for probabilistic computer models" (to appear).

All codes are written in MATLAB and as such require at least MATLAB 2016b to run. You will also need to have the Statistics and Machine Learning Toolbox by MathWorks installed.
I have provided the main code, psus.m, dependencies in the `Dependencies` folder and some `Examples`.  

## Getting started
### Instalation
There's no installation required to use psus, provided you have a working copy of MATLAB 2016b or above. As usual, you can clone the code by typing

    git clone https://github.com/PeterHristov/psus.git

in your favourite git terminal.

### Use
You can check out the demos in the `Examples` folder to get started. If you save your scripts in the `psus` folder you will need to add the `Dependencies` folder to your MATLAB search path using

	addpath('./Dependencies')
	
otherwise you will also need to add the psus folder to the search path

	addpath(genpath('<path to the cloning folder>psus'))

Let's take a look at the `psus_demo.m` file, whcih is meant to be the demo driver

	%% Set up
	d = 2;
	t_star = 230;
	N = 1000;
	p0 = 0.1;
	mode = 'neg';

	outD = 'norm';
	inpD = {'unif',[0,1]};

	%% Run 1 iter
	func = @(x)branin_uncert(x,1,mode); %Largest level of uncertainty
	[~,sOut{1}] = psus(func,d,t_star,N,p0,outD,inpD);
		
	%% Run multiple iter's
	iter = [1,5,10,15,25,50,75,100,250:250:1000]';

	for i = 2:length(iter)
		func = @(x)branin_uncert(x,iter(i),mode); %Shrinking level of uncert
		[~,sOut{i,1}] = psus(func,d,t_star,N,p0,outD,inpD);
		disp(iter(i))
	end

There are several things happening in a few sections. Let's go through them one by one.

	%% Set up
	d = 2;
	t_star = 230;
	N = 1000;
	p0 = 0.1;
	mode = 'neg';
	
	outD = 'norm';
	inpD = {'unif',[0,1]};
	
The first section is responsible for setting the main parameters of the test problem and the algorithm. I have tried to use a notation consistent with the paper. In its initial version `psus_demo.m` analyses the well-known branin fucntion, which has been *distributionalised* (as described in the paper). This fucntion is two-dimensional, so `d = 2`. The `mode` and `t_star` parameters have been chosen to illustrate the ideas in the paper. The idea behind `mode = 'neg'` is that one of the two failure domains will disappear at low fidelities (set via the second parameter in the model). The `inpD` and `outD` variables hold information about the input and output distribution of the model, respectively, in this case a standard uniform and a normal. As intended, the no parameters are passed for the output distribution, since these vary with the model input. Refer to the paper for more information on the admissible distrbutions (or try to run with whatever distributions and the code should produce an informative error if you don't get it right). The other parameters, `N` and `p0` set the behaviour of the PSuS. Again, check the paper to understand how they affect the results.

	%% Run 1 iter
	func = @(x)branin_uncert(x,1,mode); %Largest level of uncertainty
	[~,sOut] = psus(func,d,t_star,N,p0,outD,inpD);
	
Here I am running PSuS with the function at its lowest fidelity (`iter = 1`). I am choosing to suppress the probability of failure output (the first output), because this is contained in the `sOut` struct. The first output of `psus.m` summarises contains the probability of failure and running the code with this output only offers a little speed up as no strctures are written. **A note on efficiency!** This is a research piece of code and while I have made a sincere effort to avoid unnecessary overheads and bottlenecks, I have not spent time on explicitly optimising the code.

Since `psus.m` only accepts functions with a single input, you must cast your test function in that form using a function handle. Of course you can write your own functions to test the algorithm (this would be very much appreciated), but you have to have in mind the following if you decide to use the code as is

 - Test functions can only have a single input parameter (see above).
 - Test functions should have two outputs, not an 2-by-1 or 1-by-2 array of outputs. These are the two parameters of the `outD` distribution.
 - The functions must be vectorised, as the code does not provide a one-at-a-time running feature.
 
 
	%% Run multiple iter's
	iter = [1,5,10,15,25,50,75,100,250:250:1000]';

	for i = 2:length(iter)
		func = @(x)branin_uncert(x,iter(i),mode); %Shrinking level of uncert
		[~,sOut(i,1)] = psus(func,d,t_star,N,p0,outD,inpD);
		disp(iter(i))
	end
	
This block runs PSuS for progressively finer quality approximations (`iter` from `2` to `1000`) to demonstrate the fact that PSuS generalises SuS.

### Contributions
The main contribution you can make to this effort is to use it with probabilistic computer models and get back to me at p[dotty]hristov[at]liv[dotter]ac[dotted]uk
#### License
The code is distributed under GPL v3. Please check the license file.
