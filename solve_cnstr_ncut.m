function [clusters, cncut] = solve_cnstr_ncut(W, fold, gdeg,cdeg, gamma, maxiterations, verbosity)

 
	[ix,jx,wval]=find(W);

    Wtriu = triu(W);
    
    rvalold=zeros(length(ix),1); 
    pars.MAXITER=100;
    pars.epsilon = 1E-6; 

	counter=0;	
	FctValOld=inf;
	FctValOuter=[]; 
    fnew=fold; 

    n = size(W,1);

    dual_Obj = inf;
    MAX_ITERS_FLAG = true;
    
    Obj = -1;
    
    if gamma*k == 0
        disp('starting the unconstrained NCut minimization');
    else
        disp('starting the volume constrained NCut minimization');
    end
        
    while MAX_ITERS_FLAG %(counter<maxiterations)
	  
    		sval = wval.*abs(fnew(ix)-fnew(jx));

        Pfnew = fnew - (fnew'*gdeg/sum(gdeg));
        FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(gdeg'*abs(Pfnew));
	
	  
        if(FctVal>=FctValOld)
            
            if(verbosity)
               disp(['Functional has not decreased. Old: ',num2str(FctValOld,'%1.16f'),' - New: ',num2str(FctVal,'%1.16f'),'. Increasing number of inner iterations to ', num2str(pars.MAXITER*2)]); 
            end

            pars.MAXITER=pars.MAXITER*2;
            if(pars.MAXITER>=maxiterations)     
                break;
            end

            fold=foldback;
            FctOld=FctValOld;
            FctVal=FctValOld;
            
        else
            
			fold = fnew;
			FctValOuter=[FctValOuter,FctVal];
            			
            foldback=fold;
        	FctOld=FctValOld;            
            FctValOld = FctVal;
            
        end
        	  
        % Compute the subgradient of the denominator.
        Deg = spdiags(deg,0,size(fold,1),size(fold,1));
        Pfold = fold - (fold'*deg/sum(deg));
        vec = Deg*sign(Pfold) - deg*sum( Deg*sign( Pfold ))/sum(deg);        

        % Compute the subgradient of R_2(f), without factor 1/2.
        r2 = 2*gamma*sparse(qix,1,qval.*sign(fold(qix) - fold(qjx)), size(W,1),1);

        % FISTA call here.
        deg_unit_W = sparse(ix, 1, 1);
        max_deg = max(deg_unit_W);

        if gamma*volQ == 0
            [fnew,rvalold,dual_Obj,niter]=mex_solve_inner_problem(Wtriu,FctVal*full(vec)/2, rvalold,pars.MAXITER,pars.epsilon,4,sqrt(max_deg));
        else
           [fnew,rvalold,v,dual_Obj,niter]=mex_solve_cnstr_inner_problem ...
               (Wtriu,FctVal*full(vec)/2,rvalold,pars.MAXITER,pars.epsilon,4,full(r2)/2,v,gamma*volQ/2, sqrt(max_deg), degree);
        end           

       if size(unique(fnew))==1
           break;
       end
       dual_Obj = -dual_Obj;

        fnew = fnew/norm(fnew);

%            Obj = sum(abs(D*fnew))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;
        Obj = sum(wval.*abs(fnew(ix)-fnew(jx)))+gamma*volQ*(max(fnew)-min(fnew))/2 - fnew'*r2/2 - FctVal*fnew'*vec/2;

        if isnan(dual_Obj) || isnan(Obj) || counter == 100
            disp('stop');
        end

        g = fnew;        
        
        if verbosity && (counter<100 || rem(counter,10) == 0)
            disp(['iter: ', num2str(counter), 'FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, g)), 'FctVal after opt_thres: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' iters: ', num2str(niter)]);
        end
          
        totIters = counter;
        if abs(dual_Obj) < 1e-8
                MAX_ITERS_FLAG = false;
        end       
        
        if counter > maxiterations
            MAX_ITERS_FLAG = false;
        end

		counter=counter+1;
    end

    disp(['iter: ', num2str(counter), 'FctVal: ', num2str(fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fnew)), ' primal Obj: ', num2str(Obj), ' dual Obj: ', num2str(dual_Obj), ' duality gap is: ', num2str(Obj - dual_Obj)]);
	% compute most recent Functional value
	sval = wval.*abs(fnew(ix)-fnew(jx));
    Pfnew = fnew - (fnew'*deg/sum(deg));
	FctVal = (sum(sval) - gamma* sum(qval.*abs(fnew(qix)-fnew(qjx))) + gamma*volQ * (max(fnew) - min(fnew)))/(deg'*abs(Pfnew));
    
    if(FctVal<FctValOld && length(unique(fnew)) > 1)
        fold = fnew;
        FctValOuter=[FctValOuter,FctVal];
    end

    fprintf('\n Final functional value: %f\n', fctval_cnstr_one_spec_Q(W, deg, Q, gamma, fold));
end

