function [acceptance]=calcAcceptanceRatio(logN1,logD1)

if logD1==-inf
    ratio=1.0;
else
    ratio=exp(logN1-logD1);
end

acceptance=min([1 ratio]);
     
