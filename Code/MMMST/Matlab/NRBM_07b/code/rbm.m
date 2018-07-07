function [best_vishid,best_hidbiases,best_visbiases] = rbm(X,numhid,varargin)
%function [vishid,hidbiases,visbiases] = rbm(X,numhid,varargin)
% 
% RBM for one layer
% if input is real valued data (with variance = 1) then set 'gaussian' to 1
%
% options : maxepoch(50) batchsize(100) verbosity(1) gaussian(0)
%

%
% Copyright (C) 2012, by Trinh-Minh-Tri Do minhtrido@gmail.com
%
%   This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

maxepoch	= get_varargin(varargin,'maxepoch',50);
batchsize	= get_varargin(varargin,'batchsize',100);
verbosity	= get_varargin(varargin,'verbosity',1);
gaussian	= get_varargin(varargin,'gaussian',1);
permutation	= get_varargin(varargin,'permutation',1);
lrate		= get_varargin(varargin,'lrate',1e-1);
rbm_init	= get_varargin(varargin,'rbm_init','randn');
D		= get_varargin(varargin,'D',3);
momentum	= get_varargin(varargin,'momentum',0);
samplingV	= get_varargin(varargin,'samplingV',1);
gsd		= get_varargin(varargin,'gsd',1);
supportX	= max(max(abs(X)));
XHoldOut	= get_varargin(varargin,'XHoldOut',[]);
weightcost	= get_varargin(varargin,'weightcost',0);

XHoldOut	= XHoldOut'; 

starting_time = cputime;
T = 1;
A = 1;
iter0 = 1e+7*100 / batchsize; % iter0 = 1M is good for gaussian visible unit
NNotBetter = 3;
if gaussian
	epsilonw      = lrate / size(X,1);   % Learning rate for weights 
	epsilonvb     = lrate / size(X,1);   % Learning rate for biases of visible units 
	epsilonhb     = lrate / size(X,1);   % Learning rate for biases of hidden units 
	if weightcost==0
		weightcost  = 2e-5;
	end
else
	epsilonw      = lrate;   % Learning rate for weights 
	epsilonvb     = lrate;   % Learning rate for biases of visible units 
	epsilonhb     = lrate;   % Learning rate for biases of hidden units 
	if weightcost==0
		weightcost  = 0.0002;   
	end
end
restart   = 1;
if momentum
	initialmomentum  = 0.5;
	finalmomentum    = 0.9;
else
	initialmomentum  = 0.0;
	finalmomentum    = 0.0;
end

[numdims n]=size(X);

numbatches = -floor(-n/batchsize);
if permutation
	r  = randperm(ceil(n));
else
	r = 1:n;
end
be = round((0:(numbatches-1))*(n/numbatches))+1;
en = [be(2:end)-1,n];

bestError = inf;
iterNotBetter = 0;

if restart ==1,
  restart=0;
  epoch=1;

% Initializing symmetric weights and biases.
  vishid     = 0.1*randn(numdims, numhid);
  hidbiases  = zeros(1,numhid);

  visbiases  = zeros(1,numdims);

  B=1;
  if ~gaussian
	scalar = 1.0;
  else
	scalar = 1 / gsd;
  end
  if gaussian
    if strcmp(rbm_init,'autoencoder')
      [vishid,hidbiases,hidvis,visbiases] = autoencoderG(X,numhid);
    elseif strcmp(rbm_init,'interval')
      B = numhid / numdims; % B bits per input unit
      for i=1:numdims
	for b=1:B
          j = (i-1)*B + b;
          vishid(i,j)  = 1;
          hidbiases(j) = -B/2 + b + 0.5;
        end
      end
      visbiases(i) = -B/2;
    elseif strcmp(rbm_init,'subset')
      ns= numdims / D; % number of segment
      B = numhid / ns; % number of hidden unit per segment
      for i=1:ns
	bvis = (i-1)*D+1; evis = i*D;
        bhid = (i-1)*B+1; ehid = i*B;
	dispV(1,verbosity,sprintf('Pre-init Gaussian RBM step %d of %d : visible %d to %d , hidden %d to %d',i,ns,bvis,evis,bhid,ehid));
        [vishid(bvis:evis,bhid:ehid) hidbiases(bhid:ehid) visbiases(bvis:evis)] = rbm(X(bvis:evis,:),B,...
		'maxepoch',maxepoch,...
		'batchsize',batchsize,...
		'verbosity',verbosity,...
		'gaussian',gaussian,...
		'permutation',permutation,...
		'lrate',lrate,...
		'rbm_init','randn');
	vishid(bvis:evis,bhid:ehid) = vishid(bvis:evis,bhid:ehid) / scalar;
      end
    end
  end
  if scalar~=1
    Xzoom = X*scalar; % Xzoom now has support B
    if ~isempty(XHoldOut)
	XHoldOut = XHoldOut*scalar;
    end
  else
    Xzoom = X;
  end

%  poshidprobs = zeros(batchsize,numhid);
%  neghidprobs = zeros(batchsize,numhid);
%  posprods    = zeros(batchsize,numhid);
%  negprods    = zeros(batchsize,numhid);
  vishidinc  = zeros(numdims,numhid);
  hidbiasinc = zeros(1,numhid);
  visbiasinc = zeros(1,numdims);
end

best_vishid = vishid;
best_hidbiases = hidbiases;
best_visbiases = visbiases;

dispV(1,verbosity,sprintf('RBM (gaussian %d momentum %d gsd %g scalar %g lrate %g) %d x %d (XHoldOut [%d x %d])',gaussian,momentum,gsd,scalar,lrate,size(X,1),numhid,size(XHoldOut))); 

iter = 1;

for epoch = epoch:maxepoch,
 errsum=0;errsum2=0;
 for batch = 1:numbatches
 iter = iter+1;
%%%%%%%%% START POSITIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  data = Xzoom(:,r(be(batch):en(batch)))';
  numcases = size(data,1);
  poshidprobs = 1./(1 + exp(-data*vishid- repmat(hidbiases,numcases,1)));    %[H2 x v] * [v * n] = [H2 x n]
  posprods    = data' * poshidprobs; % [v * n ][n * H2]
  poshidact   = sum(poshidprobs);
  posvisact = sum(data);

%%%%%%%%% END OF POSITIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1
  if gaussian
	recdata = poshidprobs*vishid' + repmat(visbiases,numcases,1);
  else
	recdata = 1./ (1 + exp(-poshidprobs*vishid' - repmat(visbiases,numcases,1)));
  end
  err= sum(sum( (data-recdata).^2 )) / (scalar^2);
else
  err = 0;
end

%%%%%%%%% START NEGATIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  neghidprobs = poshidprobs;
  for t=1:T
    poshidstates =  neghidprobs > rand(numcases,numhid);
    if gaussian
      negdata = poshidstates*vishid' + repmat(visbiases,numcases,1);
      if samplingV
	      negdata = negdata + A*randn(numcases,numdims);
      end
    else
      negdata = 1./ (1 + exp(-poshidstates*vishid' - repmat(visbiases,numcases,1)));
      if samplingV
	      negdata = negdata > rand(numcases,numdims);
      end
    end
    neghidprobs = 1 ./ (1 + exp(-negdata*vishid - repmat(hidbiases,numcases,1)));
  end
  err2= sum(sum( (data-negdata).^2 )) / (scalar^2);

  negprods  = negdata'*neghidprobs;
  neghidact = sum(neghidprobs);
  negvisact = sum(negdata); 

%%%%%%%%% END OF NEGATIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  errsum = err + errsum;
  errsum2 = err2 + errsum2;
   if iter>1000,
     momentum=finalmomentum;
   else
     momentum=initialmomentum;
   end;

%%%%%%%%% UPDATE WEIGHTS AND BIASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    stepsize = iter0 / (iter+iter0);
    vishidinc = momentum*vishidinc + stepsize*epsilonw*( (posprods-negprods)/numcases - weightcost*vishid);
    visbiasinc = momentum*visbiasinc + stepsize*(epsilonvb/numcases)*(posvisact-negvisact);
    hidbiasinc = momentum*hidbiasinc + stepsize*(epsilonhb/numcases)*(poshidact-neghidact);

    vishid = vishid + vishidinc;
    visbiases = visbiases + visbiasinc;
    hidbiases = hidbiases + hidbiasinc;

%%%%%%%%%%%%%%%% END OF UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  end

  if isempty(XHoldOut)
	rmseHoldOut = 0;
  else
	ntmp = size(XHoldOut,1);
 	tmp = 1./(1 + exp(-XHoldOut*vishid- repmat(hidbiases,ntmp,1)));    %[H2 x v] * [v * n] = [H2 x n]
  	if gaussian
        	recHoldOut = tmp*vishid' + repmat(visbiases,ntmp,1);
  	else
        	recHoldOut = 1./ (1 + exp(-tmp*vishid' - repmat(visbiases,ntmp,1)));
	end
  	rmseHoldOut = sqrt(sum(sum( (XHoldOut-recHoldOut).^2 )) / ((scalar^2)*ntmp*numdims));
  end

  dispV(1,verbosity,sprintf('epoch %4i (T=%d) %6.2fs - RMSE pixel %.4f - HoldOut RMSE pixel %.4f ', epoch, T,cputime-starting_time,sqrt(errsum/(n*numdims)),rmseHoldOut)); 
  if ~isempty(XHoldOut)
	err = rmseHoldOut;
  else
	err = errsum2;
  end
  if err<bestError
     bestError = err;
     iterNotBetter = 0;
     best_vishid = vishid;
     best_hidbiases = hidbiases;
     best_visbiases = visbiases;
  end
  if err>bestError*0.9999
     iterNotBetter = iterNotBetter +1;
     if (iterNotBetter>=NNotBetter)&(~isempty(XHoldOut))
        break;
%        T = T+1;
     end
  end
end;

best_vishid = best_vishid*scalar;
