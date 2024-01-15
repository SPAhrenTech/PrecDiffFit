function y=getCovarianceInverseM(M)
%Returns the inverse of the covariance matrix
   MtM=M'*M;   
   y=inv(MtM);
end

