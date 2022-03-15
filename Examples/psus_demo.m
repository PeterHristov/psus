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
