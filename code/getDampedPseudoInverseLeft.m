function y=getDampedPseudoInverseLeft(M,nu)
%Returns the damped left pseudoinverse
   MtM=M'*M;   
   wMtM=MtM+nu*diag(diag(MtM));
   iwMtM=inv(wMtM);
   Mt=M';
   y=iwMtM*Mt;
end

