function [logPDF] = logMultiVariateNormal(Ysamples,centeringLocation,eigenVectorMatrix,eigenValueMatrix,numberOfVariables)

 % eigenVector is the eigen matrix of the Sigma
 % eigenValueMatrix is the eigen value matrix of the Hessian
 % matrix
 
 % covariance is inverse of the Hessian


 logPDF=-numberOfVariables/2*log(2*22/7)+1/2*log(det(eigenValueMatrix))-1/2*(Ysamples-centeringLocation)*eigenVectorMatrix*eigenValueMatrix*eigenVectorMatrix'*(Ysamples-centeringLocation)';

end

