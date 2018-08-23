function [lnprior]=logprior(N,val,priorRange,pos,delta,modelProblem,typeSampling,key,keyCount,influentialReactions,meanVector,valVector)

% signalling pathway example

if modelProblem==8
    
    lnprior=0.0;
    
    for n=1:N
        if val(n)==delta
            lnprior=lnprior+log(0.5);
        else
            lnprior=lnprior+log(0.5)+lognormal(val(n),0,1/0.1);
        end
        
    end
    
elseif modelProblem==9
    
    lnprior=0.0;
    meanVector=[0.0; 0.0; 0.0; 0.0; 0.0; 0.0];    
    
    
    for n=1:N
        if val(n)==delta
            lnprior=lnprior+log(0.5);
        else
            lnprior=lnprior+log(0.5)+lognormal(val(n),meanVector(n),1/0.1);
        end
        
    end
    
elseif modelProblem==11
    
        
    lnprior=0.0;
    
    for n=1:N
        if val(n)==delta
            lnprior=lnprior+log(0.5);
        else
            lnprior=lnprior+log(0.5)+lognormal(val(n),meanVector(n),valVector(n));
        end
        
    end

    %if val(4)~=delta && val(6)~=delta && val(7)~=delta && val(9)~=delta && val(12)~=delta && val(14)~=delta
    %   lnprior=lnprior+log(100);   
    %end        
    
else
    
    lnprior=0.0;
    for n=1:N
        if val(n)==delta
            lnprior=lnprior+log(0.5);
        else
            lnprior=lnprior+log(0.5)-log(priorRange(2,n)-priorRange(1,n));
        end
    end
end

if lnprior==-Inf
    lnprior=-1E300;
end

if lnprior==Inf
    lnprior=1E300;
end





