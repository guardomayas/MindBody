%%%%% ITClh %%%%%%%

% Modified script that imports eprime extracted data and fits logistic
% regression model assuming hyperbolic discounting (Mazur model) from Joe
% Kable with revisions. Silvia 08.15.14

function info = ITClh(fileName)
% First load ITC data to workspace
addpath('/Users/silvia/OneDrive - Universidad del rosario/MATLAB/bData');
% addpath('D:\silvia.lopez\Documents\MATLAB\bData');
%addpath(['/Users/silvia/Dropbox (Personal)/MATLAB/ITC'...
    %'/Analysis/individual/cbdm_itc']);
    load(fileName,'data','sub','ses')
    
rawdata = data;    
% Take out non answered trials
rawdata(isnan(rawdata(:,5)),:) = [];

v1 = rawdata(:,1);
v2 = rawdata(:,3);
d1 = rawdata(:,2);
d2 = rawdata(:,4);
choice = rawdata(:,5);

percent_imp= 1-sum(choice)/length(choice);

fprintf('Subject ID# = %g\n',sub);
fprintf('Session # = %g\n',ses);

info = fit_discount_model(rawdata,choice,v1,d1,v2,d2,sub,ses);

fprintf('Discount Rate = %g\n',info.b(2));
fprintf('Noise = %g\n',info.b(1));
fprintf('Percent Impulsive Choice = %g\n',percent_imp);

if info.b(2) < 0.0001 || info.b(2) > 0.2512
   fprintf('WARNING: Calculated k is outside the range and is not trustworthy\n');
end

thisdir = pwd;
% % make a data directory if necessary
% if ~isdir(fullfile(thisdir,'modLH'))
%     disp('Making data directory');
%     mkdir('modLH');
% end
% 
% % save output file as structure containing all fields of 'info'
% outfile = fullfile('modLH',sprintf('DRLH_S%d_%d.mat',subject,session));
% save(outfile,'-struct','info');
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'estLH'))
    disp('Making data directory');
    mkdir('estLH');
end

% save output file as structure containing all fields of 'info'
outfile = fullfile('estLH',sprintf('LH_S%d_%d.mat',sub,ses));
save(outfile,'-struct','info');
end

%----- FIT DISCOUNTING MODEL - LOGISTIC FIT OF HYPERBOLIC DISCOUNTING
%     info = fit_discount_model(choice,v1,d1,v2,d2)
%
%     Fits a probabilistic discounting model to binary choices by maximum likelihood.
%
%     Calls local functions: log-likelihood, choice probability and
%     discount
%
%     INPUTS
%     choice    - Dependent variable. The data should be *ungrouped*,
%                   such that CHOICE is a column of 0s and 1s, where 1 indicates 
%                   a choice of OPTION 2.
%     v1        - values of option 1 (ie, sooner option)
%     d1        - delays of option 1
%     v2        - values of option 2 (ie, later option)
%     d2        - delays of option 2
%
%     OUTPUTS
%     info      - data structure with following fields:
%                     .nobs      - number of observations
%                     .nb        - number of parameters
%                     .optimizer - function minimizer used
%                     .exitflag  - see FMINSEARCH
%                     .a         - fitted parameters; note that for all the
%                                  available models, the first element of A
%                                  is a noise term for the choice
%                                  function, the remaining elements are
%                                  parameters for the selected discount
%                                  functions. eg., for dfn='exp', A(2) is
%                                  the time constant of the hyperbolic
%                                  decay.
%                     .LL        - log-likelihood evaluated at maximum
%                     .LL0       - restricted (minimal model) log-likelihood
%                     .AIC       - Akaike's Information Criterion 
%                     .BIC       - Schwartz's Bayesian Information Criterion 
%                     .r2        - pseudo r-squared
%                   This is a struct array if multiple discount functions are fit.
%     p         - Estimated choice probabilities evaluated at the values & 
%                   delays specified by the inputs v1, v2, p1, p2. This is
%                   a cell array if multiple models are fit.
%
%     REVISION HISTORY:
%     brian  03.10.06 written
%     brian  03.14.06 added fallback to FMINSEARCH, multiple fit capability
%     kenway 11.29.06 added CI evaluation for FMINSEARCH (which returns
%     Hessian)
%     joe kable 03.12.09 modified to work with revised version of
%     choice_prob and discount
%     khoi 06.26.09 simplified 
%     silvia 03.24.16 switched to fmincon to constrain parameter values

function info = fit_discount_model(datamat,choice,v1,d1,v2,d2,subject,session)

nobs = length(choice);
alpha = 0.05;

%Initialize parameters
b0 = [0.005 0.01]';

% Fit model, using FMINCON
% if exist('fminunc','file')
%   try
     optimizer = 'fmincon';
     OPTIONS = optimoptions(@fmincon,'MaxIter',1000);
     lowB = [0 0.001];
     upB = [8 6.4];
%      OPTIONS = optimoptions(@fmincon,'Display','off','Algorithm','sqp','MaxIter',1000);
     [b,negLL,exitflag,convg,~,~,H] = fmincon(@local_negLL,b0,[],[],[],[],lowB,upB,[],OPTIONS,choice,v1,d1,v2,d2);
%      if exitflag ~= 1 % trap occasional linesearch failures
%         optimizer = 'fminsearch';
%         fprintf('FMINUNC failed to converge, switching to FMINSEARCH\n');
%      end         
%   catch
%      optimizer = 'fminsearch';
%      fprintf('Problem using FMINUNC, switching to FMINSEARCH\n');
%   end
% else
%   optimizer = 'fminsearch';
% end
% 
% if strcmp(optimizer,'fminsearch')
%   optimizer = 'fminsearch';
%   OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,...
%      'DiffMinChange',1e-4,'Maxiter',1000,'MaxFunEvals',2000);
%   [b,negLL,exitflag,convg] = fminsearch(@local_negLL,b0,OPTIONS,choice,v1,d1,v2,d2);
% end

% if exitflag ~= 1
%   fprintf('Optimization FAILED, #iterations = %g\n',convg.iterations);
% else
%   fprintf('Optimization CONVERGED, #iterations = %g\n',convg.iterations);
% end

fprintf('#iterations = %g\n',convg.iterations);
     
% Choice probabilities (for OPTION 2)
p = choice_prob(v1,d1,v2,d2,b);
% Unrestricted log-likelihood
LL = -negLL;
% Restricted log-likelihood
LL0 = sum((choice==1).*log(0.5) + (1 - (choice==1)).*log(0.5));

% Confidence interval, requires Hessian from FMINCON
try
   invH = inv(-H);
   se = sqrt(diag(-invH));
catch
end

info.subject = subject;
info.session = session;
info.nobs = nobs;
info.nb = length(b);
info.optimizer = optimizer;
info.exitflag = exitflag;
info.b = b;

try
   info.se = se;
   info.ci = [b-se*norminv(1-alpha/2) b+se*norminv(1-alpha/2)]; % Wald confidence
   info.tstat = b./se;
catch
end

info.p = p;
info.LL = LL;
info.LL0 = LL0;
info.AIC = -2*LL + 2*length(b);
info.BIC = -2*LL + length(b)*log(nobs);
info.r2 = 1 - LL/LL0;
info.percentimp = 1-sum(choice)/length(choice); 
info.percentdelayed = sum(choice)/length(choice); 
info.rawdata = datamat; %data table
% filesep '/Users/silvia/Dropbox/MATLAB/ITC/Analysis/individual'

end

%----- LOG-LIKELIHOOD FUNCTION
function sumerr = local_negLL(beta,choice,v1,d1,v2,d2)

p = choice_prob(v1,d1,v2,d2,beta);

% Trap log(0)
ind = p == 1;
p(ind) = 0.9999;
ind = p == 0;
p(ind) = 0.0001;
% Log-likelihood
err = (choice==1).*log(p) + (1 - (choice==1)).*log(1-p);
% Sum of -log-likelihood
sumerr = -sum(err);
end


%----- DISCOUNT FUNCTION - HYPERBOLIC
%     y = discount(v,d,kappa)
%
%     INPUTS
%     v     - values
%     d     - delays
%     kappa  - discount rate
%
%     OUTPUTS
%     y     - discounted values
%
%     REVISION HISTORY:
%     brian lau 03.14.06 written
%     khoi 06.26.09 simplified 

function y = discount(v,d,kappa)

y = v./(1+kappa*d);

end

%----- CHOICE PROBABILITY FUNCTION - LOGIT
%     p = choice_prob(v1,d1,v2,d2,beta);
%
%     INPUTS
%     v1    - values of option 1 (ie, sooner option)
%     d1    - delays of option 1
%     v2    - values of option 2 (ie, later option)
%     d2    - delays of option 2
%     beta  - parameters, noise term (1) and discount rate (2)
%
%     OUTPUTS
%     p     - choice probabilities for the **OPTION 2**
%
%     REVISION HISTORY:
%     brian lau 03.14.06 written
%     khoi 06.26.09 simplified 

function p = choice_prob(v1,d1,v2,d2,beta)

u1 = discount(v1,d1,beta(2));
u2 = discount(v2,d2,beta(2));

% logit, smaller beta = larger error
p = 1 ./ (1 + exp(beta(1)*(u1-u2)));
end

