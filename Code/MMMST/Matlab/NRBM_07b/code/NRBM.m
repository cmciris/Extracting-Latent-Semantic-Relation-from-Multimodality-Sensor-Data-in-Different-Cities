function [wbest stats] = NRBM(fhandle,auxdata,lambda,w0,varargin)
%function [wbest stats] = NRBM(fhandle,auxdata,lambda,w0,varargin)
%
% Input :
%         fhandle : handle of R/dR computing function (value and gradient)
%         auxdata : data used by R (auxdata is constant in the definition of R)
%          lambda : regularization parameter (>0)
%              w0 : [dim x 1] initial solution
%
% Optional options : set them in varargin e.g. wbest = NRBM(fhandle,auxdata,lambda,w0,'maxiter','100','epsilon',1e-3);
%            wreg : regularization point                                  (default=   w0)
%             reg : regularization weight                                 (default= ones(dim,1))
%
%         freport : handle of report function                             (default=  [])
%         maxiter : maximum number of iteration                           (default=  500)
%           maxCP : maximum number of cutting plane in approx. problem    (default=  200) ! NRBM working memory = maxCP x dim
%         epsilon : relative tolerance (wrt. value of f)                  (default= 0.01)
%       Rpositive : tell if R(w) is always positive                       (default=    1)
%         Rconvex : tell if R(w) is convex                                (default=    0)
%    computeGapQP : verify the solution of approximation                  (default=    1)
%       verbosity : level of verbosity from 0 to 2                        (default=    1)
%              LS : activate line search                                  (default=    1)
%     CPMinApprox : build cutting plane at minimier of approx. problem    (default=    0)
%
% Output :
%           wbest : best observed solution
%           stats : struct containing
%                       fbest : best observed value of f
%                      numite : number of iterations
%                    numfeval : number of function evaluations
%                  timeReport : time for reporting
%                   timeFGrad : time for function evaluation / computing gradient
%                    timeNRBM : time of nrbm
%                   timeTotal : total time
%
% This function optimize the following unconstrained problem :
%
%     min 0.5*lambda*((w-wreg).*reg)*((w-wreg).*reg)' + R(w)
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
%

% Last revised Nov 23, 2012


	dim = numel(w0);
	
	wreg      = get_varargin(varargin,    'wreg',          w0);
	reg       = get_varargin(varargin,     'reg', ones(dim,1));
	freport   = get_varargin(varargin, 'freport',          []);

% we introduce the new variable 
%     wnew = (w-wreg).* reg 
% <=> w    =  wnew ./ reg + wreg
% then R(w) = R( wnew ./ reg + wreg)
% defining a new function Rn(wnew) = R( wnew ./ reg + wreg)
% the minimization problem is equivalent to
%
%     min 0.5*lambda*wnew*wnew' + Rn(wnew)                         
%
% and the gradient is computed by: 
%   d Rn(wnew) / d wnew = d R( wnew ./ reg + w0) / d wnew = (d R(w) / d w) * ( d w / d wnew ) = (d R(w) / d wnew) ./ reg
%

	w0 = w0(:);wreg=wreg(:);reg=reg(:);
	wn0 = (w0 - wreg).* reg;
	
	[wnbest stats] = NRBM_kernel(@fnew,@freportnew,{auxdata fhandle freport wreg reg},lambda,wn0,varargin{:});
	wbest = wnbest ./ reg + wreg;

end

%==========================================================================
function [Remp gradwn] = fnew(wnew,auxdatanew)
	[auxdata fhandle freport wreg reg] = deal(auxdatanew{:});
	w = wnew ./ reg + wreg;
	[Remp gradw] = feval(fhandle,w,auxdata);
	gradwn = gradw(:) ./ reg;
end

%==========================================================================
function str = freportnew(wnew,auxdatanew)
	[auxdata f freport wreg reg] = deal(auxdatanew{:});
	if isempty(freport)
		str = '';
	else
		w = wnew ./ reg + wreg;
		str = feval(freport,w,auxdata);
	end
end

%=====================================================================================
function [wbest stats] = NRBM_kernel(fhandle,freport,auxdata,lambda,w0,varargin);
%
% Basic non-convex regularized bundle method for solving the unconstrained problem
% min_w 0.5*lambda*w*w' + R(w)
%

	%=============================================
	% get hyper parameters
	%=============================================
	maxiter      = get_varargin(varargin,      'maxiter',  500);
	maxCP        = get_varargin(varargin,        'maxCP',  200);
	epsilon      = get_varargin(varargin,      'epsilon', 0.01);
	Rpositive    = get_varargin(varargin,    'Rpositive',    1);
	Rconvex      = get_varargin(varargin,      'Rconvex',    0);
	computeGapQP = get_varargin(varargin, 'computeGapQP',    1);
	verbosity    = get_varargin(varargin,    'verbosity',    1);
	LS           = get_varargin(varargin,           'LS',    1);
	CPMinApprox  = get_varargin(varargin,  'CPMinApprox',    0);

	%=============================================
	% initialization
	%=============================================
	dimW     = numel(w0);
	dispV(2,verbosity,sprintf('NRBM:: Start! dimW=%d',dimW));
	time0	= cputime;
	timeNRBM	= 0.0;
	timeFGrad	= 0.0;
	timeReport	= 0.0;

	gVec     = zeros(dimW,maxCP);		% set of subgradients
	bias     = zeros(maxCP,1);
	Q        = zeros(maxCP,maxCP);		% precomputed gVec*gVec'
	inactive = ones(maxCP,1) * maxCP;	% number of consecutive iterations in which the cutting plane is inactive (alpha=0);
	cp_ite   = ones(maxCP,1); 			% the iteration where cp is built
	cp_ite(maxCP) = maxiter+1;			% the last slot is reserved for aggregation cutting plane
	alpha    = zeros(maxCP,1); 			% lagrangian multipliers

	s= 0;
	distCum	= zeros(maxCP,1);		% Cumulate distance 
	sumdist	= 0;
	w		= w0;			        % current solution
	wbest    = zeros(dimW,1);fbest = inf; Rbest = inf;tbest=1;jbest=1;dual = -inf;
	strbest = '';

	% workspace for cutting plane to be added
	wolfe = struct('a1', 0.5000,'a0', 0.0100,'c1', 1.0000e-04,'c2', 0.9000,'maxiter', 5,'amax', 1.1000);
	newW     = zeros(dimW,wolfe.maxiter+1); 
	newGrad  = zeros(dimW,wolfe.maxiter+1);
	newF     = zeros(1,wolfe.maxiter+1,1);
	
	% first function evaluation and fisrt cutting plane
	timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
	[newF(1),newGrad(:,1)] = FGgrad_objective(w0,  {auxdata fhandle lambda});
	timeFGrad	= timeFGrad + cputime - time0;time0=cputime;
	fstart		= newF(1);
	newW(:,1)	= w0;
	nbNew		= 1; % one cutting plane to be add in next iteration
	astar		= wolfe.a0;
	numfeval	= 1;
	gap			= inf;
	
	% ----------------------
	% main loop
	% ----------------------
	for t=1:maxiter
		% find memory slots of new cutting planes
		listCPold = find(inactive<maxCP);
		[v,idx] = sort(inactive*maxiter*10 - cp_ite,'descend'); % sort old CPs by number of inactive then by iteration number
		listCPnew = idx(1:nbNew);

		listCPold = setdiff(listCPold,listCPnew); % performance is ok for hundreds of elements, but this line must be optimized for large maxCP
		inactive(listCPnew) = 0;
		cp_ite(listCPnew) = t;
		listCP = find(inactive<maxCP);

		% -------------------------------------
		% precompute Q(.,.) for new cutting planes
		% -------------------------------------
		gVec(:,listCPnew) = newGrad(:,1:nbNew) - lambda*newW(:,1:nbNew);
		Q(listCPnew,:) = gVec(:,listCPnew)' * gVec ;
		Q(:,listCPnew) = Q(listCPnew,:)';
	
		% -----------------------------------------------------------------------------------------------
		% precompute Q(.,.) for aggregation cutting plane, this code could be optimized by working only on Q
		% -----------------------------------------------------------------------------------------------
		Q(maxCP,:) = gVec(:,maxCP)' * gVec;
		Q(:,maxCP) = Q(maxCP,:)';

		% -------------------------------------
		% adding each cutting plane to bundle
		% -------------------------------------
		wbestold = wbest;
		for k=1:nbNew
			reg  = 0.5 * lambda * newW(:,k)' * newW(:,k);
			Remp = newF(k) - reg;
			j = listCPnew(k);
			bias(j) = Remp - gVec(:,j)' * newW(:,k);
%			wVec(:,j) = newW(:,k);
			fcurrent = newF(k);
			dispV(2,verbosity,sprintf('t=%d k=%d j=%d fcurrent=%.3e reg=%.3e Remp=%.3e',t,k,j,fcurrent,reg,Remp));
			distCum(j) = norm(wbest - newW(:,k))^2;
			if(fbest > fcurrent)
				fbest	= fcurrent;
				Rbest   = Remp;
				wbest	= newW(:,k);
				tbest	= t;
				jbest	= j;
				dist = (wbest - wbestold)'*(wbest - wbestold);
				sumdist = sumdist+dist;
  				dispV(2,verbosity,sprintf('norm_wbest=%g dist=%g sumdist=%g',norm(wbest),dist,sumdist));
				distCum = distCum + dist;
				distCum(j) = 0;
			end

			if (~Rconvex)
				% solving conflict
				if (jbest==j) % descent step
					list = [listCPold;listCPnew(1:k-1)];
					for i=list'
						score = 0.5*lambda*wbest'*wbest + gVec(:,i)'*wbest + bias(i);
						gamma = max(0, score - fbest + 0.5 * lambda * distCum(i));
						bias(i) = bias(i)-gamma;
					end
				else % null step
					% estimate g_t at w_tbest
					dist = distCum(j);
					score = 0.5*lambda*dist + gVec(:,j)'*wbest + bias(j);
					if (score > Rbest) % conflict
						% trying to solve conflict by descent g_t so that g_t(w_t) = fbest
						U = Rbest - 0.5*lambda*dist - gVec(:,j)'*wbest;
						L = fbest - reg - gVec(:,j)'*newW(:,k);
						dispV(2,verbosity,sprintf('NULL_STEP_CONFLICT Rbest=%g score=%g L=%g U=%g dist=%g',Rbest,score,L,U,dist));
						if (L<=U)
							dispV(2,verbosity,'NULL_STEP_CONFLICT LEVEL_1');
							bias(j) = L; 
						else
							dispV(2,verbosity,'NULL_STEP CONFLICT LEVEL_2');
							gVec(:,j) = - lambda * wbest;
							bias(j) = fbest - (reg + gVec(:,j)'*newW(:,k));
							Q(j,:) = gVec(:,j)' * gVec;
							Q(:,j) = Q(j,:)';
						end
						score = 0.5*lambda*dist + gVec(:,j)'*wbest + bias(j);
						dispV(2,verbosity,sprintf('new_score=%g',score));
					end
				end
			end
		end

		% ------------------------------------------
		% Solving QP program
		% ------------------------------------------
		t1=cputime;
		[alpha(listCP),dual] = minimize_QP(lambda,Q(listCP,listCP),bias(listCP),Rpositive,epsilon);

		% ------------------------------------------
		% get QP program solution
		% update w and counting inactive
		% ------------------------------------------
		listA = find(alpha(listCP)>0);
		listI = find(alpha(listCP)==0);
		w = wsum_row(gVec,-alpha(listCP(listA))/lambda,listCP(listA));
		inactive(listCP(listA)) = 0;
		inactive(listCP(listI)) = inactive(listCP(listI)) + 1;
		inactive(jbest)=0; % make sure that the best point is always in the set
		
		% -----------------------------------
		% gradient aggregation
		% -----------------------------------
		gVec(:,maxCP) = - lambda * w;
		bias(maxCP)	= dual + 0.5 * lambda * w'*w;
		distCum(maxCP) = alpha(listCP)' * distCum(listCP);
		inactive(maxCP)=0; % make sure that aggregation cp is always active

		% ----------------------------------------------
		% estimate the gap of approximated dual problem
		% ----------------------------------------------
		if computeGapQP
			score = (w' * gVec)' + bias;
			if (Rpositive)
				primal = 0.5*lambda*w'*w + max( [score(listCP);0]);
			else
				primal = 0.5*lambda*w'*w + max(score(listCP));
			end
			gapQP = primal - dual;
			dispV(2,verbosity,sprintf('quadratic_programming: primal=%.3e dual=%.3e gap=%.2f%%',primal,dual,gapQP*100.0/(abs(primal))));
		end
		dispV(2,verbosity,sprintf('Time_QP_and_update=%.2f seconds',cputime-t1));
		gap   = fbest-dual;
		timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
		if((tbest==t)&&(~isempty(freport)))
			strbest = feval(freport, wbest, auxdata);
		end
		timeReport	= timeReport + cputime - time0;time0=cputime;

		%-------------------------------------------------------
		% output
		%-------------------------------------------------------
		dispV(1,verbosity,sprintf('t=%3d nfeval=%3d f=%.3e f*=%.3e R*=%.3e gap=%5.2f%% - (%s) elapsed time %.0f + %.0f + %.0f = %.0fs',t,numfeval,fcurrent,fbest,Rbest,gap*100/abs(fbest),strbest,timeNRBM,timeFGrad,timeReport,timeNRBM+timeFGrad+timeReport));
		if ((gap/abs(fbest) < epsilon)||(gap<1e-6)||(t>=maxiter))
			break;
		end
		if (~LS)
			nbNew = 1;
			timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
			[newF(1) newGrad(:,1)] = FGgrad_objective(w,  {auxdata fhandle lambda});
			timeFGrad	= timeFGrad + cputime - time0;time0=cputime;
			newW(:,1)   = w;
			numfeval    = numfeval + 1;
			continue;
		end		

		%-------------------------------------------------------
		% doing line search from wbest to w
		%-------------------------------------------------------
		search_direction = w-wbest; norm_dir = norm(search_direction);
		if (CPMinApprox) || (t==1)
			astart = 1.0;
		else
			astart = min(astar / norm_dir,1.0); % astar = astar^0.97;
			if astart==0	astart = 1.0;end
		end
		timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
		[astar wLineSearch fLineSearch gLineSearch w1 f1 g1 nfeval] = myLineSearchWolfe(wbest, fbest, gVec(:,jbest)+lambda*wbest, search_direction, astart, wolfe.amax, wolfe.c1, wolfe.c2, wolfe.maxiter, @FGgrad_objective, {auxdata fhandle lambda});
		timeFGrad	= timeFGrad + cputime - time0;time0=cputime;
		numfeval = numfeval + nfeval;
		if f1~=fLineSearch
			nbNew = 2;
			newF(1:2)      = [f1,fLineSearch];
			newW(:,1:2)    = [w1,wLineSearch];
			newGrad(:,1:2) = [g1,gLineSearch];
		else
			nbNew = 1;
			newF(1)      = fLineSearch;
			newW(:,1)    = wLineSearch;
			newGrad(:,1) = gLineSearch;
		end
		if (fbest<=fLineSearch) && (astart ~= 1)
			numfeval = numfeval + 1;
			nbNew = nbNew + 1;
			timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
			[newF(nbNew),newGrad(:,nbNew)] = FGgrad_objective(w,  {auxdata fhandle lambda});
			timeFGrad	= timeFGrad + cputime - time0;time0=cputime;
			newW(:,nbNew) = w;
		end
		astar = astar*norm_dir; % true step length
		dispV(2,verbosity,sprintf('step_length=%.3e',astar));
	end
	timeNRBM	= timeNRBM + cputime - time0;time0=cputime;
	dispV(1,verbosity,sprintf('DONE_NRBM numfeval=%d timeNRBM=%.2fs timeFGrad=%.2fs timeReport=%.2fs timeTotal=%.2fs',numfeval,timeNRBM,timeFGrad,timeReport,timeNRBM+timeFGrad+timeReport));
	stats.fbest = fbest;
	stats.numite = t;
	stats.numfeval   = numfeval;
	stats.timeNRBM	 = timeNRBM;
	stats.timeFGrad	 = timeFGrad;
	stats.timeReport = timeReport;
	stats.timeTotal  = timeNRBM + timeFGrad + timeReport;
end	

%====================================================================
function [fval grad tt] = FGgrad_objective(w, auxdata_R_lambda)
	[auxdata R lambda] = deal(auxdata_R_lambda{:});
	[Rval grad] = feval(R,w,auxdata);
	grad = grad(:); % make sure that we have column vector
	reg = 0.5 * lambda * (w'*w);
	fval = Rval + reg;
	grad = grad + lambda * w;
	tt = [];
end

%====================================================================
function [alpha,dual] = minimize_QP(lambda,Q,B,Rpositive,EPS)
	T = size(Q,1);
	if (T+Rpositive)==1
		alpha = 1;
		dual = -0.5 * Q / lambda + B;
		return;
	end
	if (Rpositive)
		B = [B;0];
		Q = [[Q,zeros(T,1)];zeros(1,T+1)];
	end
	SCALE = abs(max(max(abs(Q)))) / (1000.0 * lambda);
	Q = Q / SCALE;
	B = B / SCALE;
	listSolver = {'libqp' 'imdm' 'kowalczyk' 'keerthi'};% there are other solvers (such as 'kozinec' 'keerthi' 'mdm') but these four are OK
	for k=6:10
		for j=1:numel(listSolver)
			t0= cputime;
			if strcmpi(listSolver{j},'libqp')
				opt = struct('MaxIter',10^k,'TolRel',1e-2*EPS);
				[alpha,stat] = libqp_splx(Q/lambda,-B,1,ones(1,numel(B)),0,[],opt);
				dual = stat.QP;
				niter = stat.nIter;
			else
				[alpha,dual,stat] = gmnp(Q/lambda,-B',{'solver',listSolver{j},'tmax',10^k,'tolrel',1e-2*EPS});
				niter = stat.t;
			end
			% ======================================================
			% verify the solution of approximated dual problem
			% because gmnp returns invalid solution sometime
			% ======================================================
			eps = 1e-4;
			if (min(alpha)<-eps)||(abs(sum(alpha)-1)>eps)
				stat.exitflag = 0;
			end;
			if stat.exitflag>0
				break;
			end
		end
		if stat.exitflag~=0
			break;
		end
	end
	if stat.exitflag==0
		disp('WARNING Solving QP of approx. problem failled, do not reach enough accuracy !!!');
	end
	B = B * SCALE;
	dual = -dual*SCALE;
	if (Rpositive)
		alpha = alpha(1:T);
	end
end