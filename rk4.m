function rk4out = rk4(F,initcond,Jm,k,tN,t0)

t=t0:k:tN;

auxmat(:,1)=initcond;

for i = 1 : (tN-t0)/k
    
    a = k*feval(F, auxmat(:,i),     Jm, k*(i-1)     );
    b = k*feval(F, auxmat(:,i)+a/2, Jm, k*(i-1)+k/2 );
    c = k*feval(F, auxmat(:,i)+b/2, Jm, k*(i-1)+k/2 );
    d = k*feval(F, auxmat(:,i)+c,   Jm, k*(i-1)+k   );
    
    auxmat(:,i+1) = auxmat(:,i) + 1/6*( a + 2*b + 2*c + d );
end

rk4out=[t;auxmat(:,:)];
% First row is t. Second row is r. Third row is rdot.

end