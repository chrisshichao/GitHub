%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              DSGE Macroeconomics Modeling 
%
%              Chao Shi                
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Using the Tauchen-Hussey method, find an equalibrium numerically,
%       that is, find decision rules, and value functions. Compute and draw
%       stationary distribution of wealth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  parameters set up for Tauchen-Hussey
%  variables with u(y) for permanet and with v(y) for transitory
%

xbar_u    = 2.00;        % unconditional mean for permanet part
xbar_v    = 2.00;        % unconditional mean for transitory part
rho_u     = 0.00;        % rho for permanet part
rho_v     = 0.97;        % rho for transitory part
sigma_z_u = 0.05;        % variance for permanet part
sigma_z_v = 0.02;        % variance for transitory part
n_u       = 5;           % 5-state Markov chain for permanet part
n_v       = 2;           % 2-state Markov chain for transitory part

%
%  applying Tauchen-Hussey to find probability matrix and income
%  for both parts
%

%  permanet part Tauchen-Hussey
[s_u,p_u]=tauch_hussAR1(xbar_u,rho_u,sigma_z_u,n_u);
%  transitory part Tauchen-Hussey
[s_v,p_v]=tauch_hussAR1(xbar_v,rho_v,sigma_z_v,n_v);

%
%  generating cross sectional probability matrix by kronecker product
%  also generating cross sectional income by sum of sum
%

y=zeros(10,1);          % set up income vector with all zeros
count=1;                % fill up income vecotr by sum of sum
for i = 1 : 5
    for j = 1 : 2
        y(count) = 0.97+s_u(i)+s_v(j);
        count=count+1;
    end
end 
y=y';
p=kron(p_u,p_v);        % probability matrix by kronecker product

%
%  parameters set up for economic environment
%

nm        = 10;              % number of total Markov chain
sigma     = 1.50;            % risk aversion              
beta      = 0.98;            % subjective discount factor 
prob      = p   ;            % probability matrix
alpha     = 0.33;            % capital's share of income
A_initial = 10.0;            % initial value for aggregate asset
g         = 0.20;            % relaxation parameter
for i = 1:10                 % different income fraction 
    theta(1,i)=i/10;
end;

%
%  form asset grid
%   
max_a     = 15;                      % maximum value of asset grid   
inc_a     = 0.15;                    % size of asset grid increments
nap      = round(max_a/inc_a);      % number of grid points
kgrid     = [ 0.15:inc_a:max_a ]';   % form asset grid

%
%  iteration to find converge point for agregate asset
%

liter     = 1;               % iteration times
maxiter   = 100;             % upper bond for iteration
tol       = 0.001;           % converge tolerance
differ    = 10;              % difference between a' and a         
A         = A_initial;       % value for inital asset
tme       = cputime;         % record cpu time

disp('ITERATING ON ASSET PATH');
disp('');
disp('    Iteration Seconds   Differ    MeanA    OldA');

while  (differ > tol) && (liter <= maxiter);
   
   %
   %  calculate rental rate of asset
   %
   %%%r=.05;
   r = 0.05*(alpha) * A^(alpha-1);
   % 
   %  revise utility function such that negetive or zero utility values
   %  will be replaced by a large negetive number so that those values
   %  won't be chosen as max value      
   %
   %  calculte different utility scheme
   %
   
   util = repmat(-10000,[nap,nap,nm]);
   for ii=1:nm;
       for i=1:nap;
           ap=(i-1)*inc_a;
           for j=1:nap;
               app = (j-1)*inc_a;
               cons = y(1,ii)*theta(1,ii) + (r + 1)*ap - app;
               if cons > 0;
                   util(j,i,ii)=(cons)^(1-sigma)/(1-sigma);
               end;
           end;
       end;
   end;
   
   %
   %  set up value function and decision rules
   %
   
   v       = zeros(nap,nm);           % value function with all zeros
   decis   = zeros(nap,nm);           % decision rule with all zeros
   [cs,rs] = size(util);              % obtain dimension of utiliy matrix
   test    = 10;                      % test statistic for converge
   
   %
   %  iterate on Bellman equation and get the decision 
   %  rules and the value function at the optimum         
   %
   
   while test ~= 0;
       for ii=1:nm;
           for i=1:cs;
               sum_r=zeros(nap,1);
               for j=1:nm;
                   sum_r=sum_r+prob(ii,j)*v(:,j);
               end;
               u(:,i,ii)=util(:,i,ii)+beta*sum_r;
           end; 
           [tv_temp, tdecis_temp]=max(u(:,:,ii));
           tdecis(:,ii) = tdecis_temp';
           tv(:,ii) = tv_temp';
       end;
       test=max(any(tdecis-decis));
       v=tv;
       decis=tdecis;      
   end;
   decis   = 0.5*(decis-1)*inc_a;         % update decision rule
   
   %
   %   form transition matrix as in class mentioned
   %   trans will be a 1000x1000 matrix in this case
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 
   %  
   
   g_new   = zeros(cs,cs,nm);         % zeros marix with dimension 100x100
   for ii=1:nm;
       for i=1:cs
           g_new(i,tdecis(i,ii),ii)=1;
       end       
   end;

   for ii=1:nm;
       for i=1:nm;
           first_index_start=(ii-1)*cs+1;
           second_index_start=(i-1)*cs+1;
           first_index_end=ii*cs;
           second_index_end=i*cs;
           trans(first_index_start:first_index_end,second_index_start:second_index_end)=prob(ii,i)*g_new(:,:,ii);
       end;
   end;
   
   %
   %  generate probst matrix and update it through
   %  transition matrix
   %
   
   trans   = trans';                         % transpose transition matrix for calculation
   probst  = (1/(nm*nap))*ones(nm*nap,1);    % create probst matrix with dim 100x1
   test    = 1;                              % test statistic for converge
   while test > 10^(-8);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1;
   end;
   
   %
   %  calculate new aggregate asset and meanA
   %
   
   vect_decsic     = decis(:);                % vectorize the decision rule for plot
   meanA           = probst'*vect_decsic;     % calculate meanA
   
   %
   %  calculate measure over (A,y) pairs
   %
   
   distri_shell    = zeros(cs,nm);            % create the shell for probst
   distri_shell(:) = probst;                  % filling the shell
   
   %
   %  calculate stationary distribution of A
   %
   
   [v1,d1]         = eig(prob');              % eigenvalue of prob matrix      
   [dmax,imax]     = max(diag(d1));           % find max value on diagonal
   probst1         = v1(:,imax);              % find corresponding column
   probst1_div     = sum(probst1);            % sum of all colmuns          
   probst1         = probst1/probst1_div;     % probst1 weighted
   proba           = sum(distri_shell');      % stationary distribution of asset'
   proba           = proba';
   
   %
   %  form difference of meanA and A_old to update A
   %
   
   A_old           = A;                       % set old asset value      
   A_new           = g*meanA + (1-g)*A_old;   % calculate new asset value
   differ          = abs((A_old-meanA)/A_old);% difference used to test convergence
   A               = A_new;                   % update asset value
   
   disp([ liter cputime-tme differ meanA A_old ]);
   liter           = liter+1;
end;
tme2               = cputime;                 % record cpu time after iteration

disp('----------------------------------------------------------------------');
fprintf('The Aiyagari method takes %u iterations to converge.\n',liter);
fprintf('Total time elapsed is %.4f seconds with converge meanA = %.4f.\n',tme2-tme,meanA);
disp('----------------------------------------------------------------------');
disp('');
disp('final period value functions are:');
disp([v(100,:)]);
disp('final period decision rules are:');
disp([decis(100,:)]);

figure(1)
plot(kgrid',v(:,1),kgrid',v(:,2),kgrid',v(:,3),kgrid',v(:,4),kgrid',v(:,5),...
     kgrid',v(:,6),kgrid',v(:,7),kgrid',v(:,8),kgrid',v(:,9),kgrid',v(:,10))
title('VALUE FUNCTION PATH');
xlabel('Period');
ylabel('Value Function');

figure(2)
plot(kgrid',decis(:,1),kgrid',decis(:,2),kgrid',decis(:,3),kgrid',decis(:,4),...
     kgrid',decis(:,5),kgrid',decis(:,6),kgrid',decis(:,7),kgrid',decis(:,8),...
     kgrid',decis(:,9),kgrid',decis(:,10));
title('DECISION RULE PATH');
xlabel('Period');
ylabel('Decision Rule');

figure(3)
plot(kgrid,proba);
title('DISTRIBUTION OF ASSET');
xlabel('Period');
ylabel('Asset');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Imagine that the shocks increase in two speperate cases:
%       one with an increase in variance of permanet shocks by 0.02,
%       and the other with an increase in variance of tranitory shocks by
%       0.02. Notice that in two cases, the magnitude of an increase in
%       variance of shocks is the same
%        i) Compute the distribution of wealth in two cases
%       ii) Compute the decision rule for consumption in two cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  All other steps are indentical to part (a)
%  except now case1 for sigma_z_u=0.07 and case2 for sigma_z_v=0.04
%
%  parameters set up for Tauchen-Hussey
%  variables with u(y) for permanet and with v(y) for transitory
%

xbar_u    = 2.00;        % unconditional mean for permanet part
xbar_v    = 2.00;        % unconditional mean for transitory part
rho_u     = 0.00;        % rho for permanet part
rho_v     = 0.97;        % rho for transitory part
sigma_z_u = 0.07;        % variance for permanet part
sigma_z_v = 0.02;        % variance for transitory part
n_u       = 5;           % 5-state Markov chain for permanet part
n_v       = 2;           % 2-state Markov chain for transitory part

%
%  applying Tauchen-Hussey to find probability matrix and income
%  for both parts
%

%  permanet part Tauchen-Hussey
[s_u,p_u]=tauch_hussAR1(xbar_u,rho_u,sigma_z_u,n_u);
%  transitory part Tauchen-Hussey
[s_v,p_v]=tauch_hussAR1(xbar_v,rho_v,sigma_z_v,n_v);

%
%  generating cross sectional probability matrix by kronecker product
%  also generating cross sectional income by sum of sum
%

y=zeros(10,1);          % set up income vector with all zeros
count=1;                % fill up income vecotr by sum of sum
for i = 1 : 5
    for j = 1 : 2
        y(count) = 0.97+s_u(i)+s_v(j);
        count=count+1;
    end
end 
y=y';
p=kron(p_u,p_v);        % probability matrix by kronecker product

%
%  parameters set up for economic environment
%

nm        = 10;              % number of total Markov chain
sigma     = 1.50;            % risk aversion              
beta      = 0.98;            % subjective discount factor 
prob      = p   ;            % probability matrix
alpha     = 0.33;            % capital's share of income
A_initial = 10.0;            % initial value for aggregate asset
g         = 0.20;            % relaxation parameter
for i = 1:10                 % different income fraction 
    theta(1,i)=i/10;
end;

%
%  form asset grid
%   
max_a     = 15;                      % maximum value of asset grid   
inc_a     = 0.15;                    % size of asset grid increments
nap      = round(max_a/inc_a);      % number of grid points
kgrid     = [ 0.15:inc_a:max_a ]';   % form asset grid

%
%  iteration to find converge point for agregate asset
%

liter     = 1;               % iteration times
maxiter   = 50;              % upper bond for iteration
tol       = 0.01;           % converge tolerance
differ    = 10;              % difference between a' and a         
A         = A_initial;       % value for inital asset
tme       = cputime;         % record cpu time

disp('ITERATING ON ASSET PATH');
disp('');
disp('    Iteration Seconds   Differ    MeanA    OldA');

while  (differ > tol) && (liter <= maxiter);
   
   %
   %  calculate rental rate of asset
   %
   %%%r=.05;
   r = 0.03*(alpha) * A^(alpha-1);
   % 
   %  revise utility function such that negetive or zero utility values
   %  will be replaced by a large negetive number so that those values
   %  won't be chosen as max value      
   %
   %  calculte different utility scheme
   %
   
   util = repmat(-10000,[nap,nap,nm]);
   for ii=1:nm;
       for i=1:nap;
           ap=(i-1)*inc_a;
           for j=1:nap;
               app = (j-1)*inc_a;
               cons = y(1,ii)*theta(1,ii) + (r + 1)*ap - app;
               if cons > 0;
                   util(j,i,ii)=(cons)^(1-sigma)/(1-sigma);
               end;
           end;
       end;
   end;
   
   %
   %  set up value function and decision rules
   %
   
   v       = zeros(nap,nm);           % value function with all zeros
   decis   = zeros(nap,nm);           % decision rule with all zeros
   [cs,rs] = size(util);              % obtain dimension of utiliy matrix
   test    = 10;                      % test statistic for converge
   
   %
   %  iterate on Bellman equation and get the decision 
   %  rules and the value function at the optimum         
   %
   
   while test ~= 0;
       for ii=1:nm;
           for i=1:cs;
               sum_r=zeros(nap,1);
               for j=1:nm;
                   sum_r=sum_r+prob(ii,j)*v(:,j);
               end;
               u(:,i,ii)=util(:,i,ii)+beta*sum_r;
           end; 
           [tv_temp, tdecis_temp]=max(u(:,:,ii));
           tdecis(:,ii) = tdecis_temp';
           tv(:,ii) = tv_temp';
       end;
       test=max(any(tdecis-decis));
       v=tv;
       decis=tdecis;      
   end;
   decis   = (decis-1)*inc_a;         % update decision rule
   
   %
   %   form transition matrix as in class mentioned
   %   trans will be a 1000x1000 matrix in this case
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 
   %  
   
   g_new   = zeros(cs,cs,nm);         % zeros marix with dimension 100x100
   for ii=1:nm;
       for i=1:cs
           g_new(i,tdecis(i,ii),ii)=1;
       end       
   end;

   for ii=1:nm;
       for i=1:nm;
           first_index_start=(ii-1)*cs+1;
           second_index_start=(i-1)*cs+1;
           first_index_end=ii*cs;
           second_index_end=i*cs;
           trans(first_index_start:first_index_end,second_index_start:second_index_end)=prob(ii,i)*g_new(:,:,ii);
       end;
   end;
   
   %
   %  generate probst matrix and update it through
   %  transition matrix
   %
   
   trans   = trans';                         % transpose transition matrix for calculation
   probst  = (1/(nm*nap))*ones(nm*nap,1);    % create probst matrix with dim 100x1
   test    = 1;                              % test statistic for converge
   while test > 10^(-8);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1; 
   end;
   
   %
   %  calculate new aggregate asset and meanA
   %
   
   vect_decsic     = decis(:);                % vectorize the decision rule for plot
   meanA           = probst'*vect_decsic;     % calculate meanA
   
   %
   %  calculate measure over (A,y) pairs
   %
   
   distri_shell    = zeros(cs,nm);            % create the shell for probst
   distri_shell(:) = probst;                  % filling the shell
   
   %
   %  calculate stationary distribution of A
   %
   
   [v1,d1]         = eig(prob');              % eigenvalue of prob matrix      
   [dmax,imax]     = max(diag(d1));           % find max value on diagonal
   probst1         = v1(:,imax);              % find corresponding column
   probst1_div     = sum(probst1);            % sum of all colmuns          
   probst1         = probst1/probst1_div;     % probst1 weighted
   proba           = sum(distri_shell');      % stationary distribution of asset'
   proba           = proba';
   
   %
   %  form difference of meanA and A_old to update A
   %
   
   A_old           = A;                       % set old asset value      
   A_new           = g*meanA + (1-g)*A_old;   % calculate new asset value
   differ          = abs((A_old-meanA)/A_old);% difference used to test convergence
   A               = A_new;                   % update asset value
   
   disp([ liter cputime-tme differ meanA A_old ]);
   liter           = liter+1;
end;
tme2               = cputime;                 % record cpu time after iteration



figure(1)
plot(kgrid,proba);
title('COMPARISION FOR DISTRIBUTION OF ASSET');
xlabel('Period');
ylabel('Asset');
hold on;


xbar_u    = 2.00;        % unconditional mean for permanet part
xbar_v    = 2.00;        % unconditional mean for transitory part
rho_u     = 0.00;        % rho for permanet part
rho_v     = 0.97;        % rho for transitory part
sigma_z_u = 0.05;        % variance for permanet part
sigma_z_v = 0.04;        % variance for transitory part
n_u       = 5;           % 5-state Markov chain for permanet part
n_v       = 2;           % 2-state Markov chain for transitory part

%
%  applying Tauchen-Hussey to find probability matrix and income
%  for both parts
%

%  permanet part Tauchen-Hussey
[s_u,p_u]=tauch_hussAR1(xbar_u,rho_u,sigma_z_u,n_u);
%  transitory part Tauchen-Hussey
[s_v,p_v]=tauch_hussAR1(xbar_v,rho_v,sigma_z_v,n_v);

%
%  generating cross sectional probability matrix by kronecker product
%  also generating cross sectional income by sum of sum
%

y=zeros(10,1);          % set up income vector with all zeros
count=1;                % fill up income vecotr by sum of sum
for i = 1 : 5
    for j = 1 : 2
        y(count) = 0.97+s_u(i)+s_v(j);
        count=count+1;
    end
end 
y=y';
p=kron(p_u,p_v);        % probability matrix by kronecker product

%
%  parameters set up for economic environment
%

nm        = 10;              % number of total Markov chain
sigma     = 1.50;            % risk aversion              
beta      = 0.98;            % subjective discount factor 
prob      = p   ;            % probability matrix
alpha     = 0.33;            % capital's share of income
A_initial = 10.0;            % initial value for aggregate asset
g         = 0.20;            % relaxation parameter
for i = 1:10                 % different income fraction 
    theta(1,i)=i/10;
end;

%
%  form asset grid
%   
max_a     = 15;                      % maximum value of asset grid   
inc_a     = 0.15;                    % size of asset grid increments
nap      = round(max_a/inc_a);      % number of grid points
kgrid     = [ 0.15:inc_a:max_a ]';   % form asset grid

%
%  iteration to find converge point for agregate asset
%

liter     = 1;               % iteration times
maxiter   = 50;              % upper bond for iteration
tol       = 0.01;           % converge tolerance
differ    = 10;              % difference between a' and a         
A         = A_initial;       % value for inital asset
tme       = cputime;         % record cpu time

disp('ITERATING ON ASSET PATH');
disp('');
disp('    Iteration Seconds   Differ    MeanA    OldA');

while  (differ > tol) && (liter <= maxiter);
   
   %
   %  calculate rental rate of asset
   %
   %%%r=.05;
   r = 0.05*(alpha) * A^(alpha-1);
   % 
   %  revise utility function such that negetive or zero utility values
   %  will be replaced by a large negetive number so that those values
   %  won't be chosen as max value      
   %
   %  calculte different utility scheme
   %
   
   util = repmat(-10000,[nap,nap,nm]);
   for ii=1:nm;
       for i=1:nap;
           ap=(i-1)*inc_a;
           for j=1:nap;
               app = (j-1)*inc_a;
               cons = y(1,ii)*theta(1,ii) + (r + 1)*ap - app;
               if cons > 0;
                   util(j,i,ii)=(cons)^(1-sigma)/(1-sigma);
               end;
           end;
       end;
   end;
   
   %
   %  set up value function and decision rules
   %
   
   v       = zeros(nap,nm);           % value function with all zeros
   decis   = zeros(nap,nm);           % decision rule with all zeros
   [cs,rs] = size(util);              % obtain dimension of utiliy matrix
   test    = 10;                      % test statistic for converge
   
   %
   %  iterate on Bellman equation and get the decision 
   %  rules and the value function at the optimum         
   %
   
   while test ~= 0;
       for ii=1:nm;
           for i=1:cs;
               sum_r=zeros(nap,1);
               for j=1:nm;
                   sum_r=sum_r+prob(ii,j)*v(:,j);
               end;
               u(:,i,ii)=util(:,i,ii)+beta*sum_r;
           end; 
           [tv_temp, tdecis_temp]=max(u(:,:,ii));
           tdecis(:,ii) = tdecis_temp';
           tv(:,ii) = tv_temp';
       end;
       test=max(any(tdecis-decis));
       v=tv;
       decis=tdecis;      
   end;
   decis   = (decis-1)*inc_a;         % update decision rule
   
   %
   %   form transition matrix as in class mentioned
   %   trans will be a 1000x1000 matrix in this case
   %   The eigenvector associated with the unit eigenvalue
   %   of trans' is  the stationary distribution. 
   %  
   
   g_new   = zeros(cs,cs,nm);         % zeros marix with dimension 100x100
   for ii=1:nm;
       for i=1:cs
           g_new(i,tdecis(i,ii),ii)=1;
       end       
   end;

   for ii=1:nm;
       for i=1:nm;
           first_index_start=(ii-1)*cs+1;
           second_index_start=(i-1)*cs+1;
           first_index_end=ii*cs;
           second_index_end=i*cs;
           trans(first_index_start:first_index_end,second_index_start:second_index_end)=prob(ii,i)*g_new(:,:,ii);
       end;
   end;
   
   %
   %  generate probst matrix and update it through
   %  transition matrix
   %
   
   trans   = trans';                         % transpose transition matrix for calculation
   probst  = (1/(nm*nap))*ones(nm*nap,1);    % create probst matrix with dim 100x1
   test    = 1;                              % test statistic for converge
   while test > 10^(-8);
      probst1 = trans*probst;
      test = max(abs(probst1-probst));
      probst = probst1;
   end;
   
   %
   %  calculate new aggregate asset and meanA
   %
   
   vect_decsic     = decis(:);                % vectorize the decision rule for plot
   meanA           = probst'*vect_decsic;     % calculate meanA
   
   %
   %  calculate measure over (A,y) pairs
   %
   
   distri_shell    = zeros(cs,nm);            % create the shell for probst
   distri_shell(:) = probst;                  % filling the shell
   
   %
   %  calculate stationary distribution of A
   %
   
   [v1,d1]         = eig(prob');              % eigenvalue of prob matrix      
   [dmax,imax]     = max(diag(d1));           % find max value on diagonal
   probst1         = v1(:,imax);              % find corresponding column
   probst1_div     = sum(probst1);            % sum of all colmuns          
   probst1         = probst1/probst1_div;     % probst1 weighted
   proba           = sum(distri_shell');      % stationary distribution of asset'
   proba           = proba';
   
   %
   %  form difference of meanA and A_old to update A
   %
   
   A_old           = A;                       % set old asset value      
   A_new           = g*meanA + (1-g)*A_old;   % calculate new asset value
   differ          = abs((A_old-meanA)/A_old);% difference used to test convergence
   A               = A_new;                   % update asset value
   
   disp([ liter cputime-tme differ meanA A_old ]);
   liter           = liter+1;
end;
tme2               = cputime;                 % record cpu time after iteration

figure(1)
plot(kgrid,proba,'r');
legend('permanent','transitory','Location','northeast');
hold off;
