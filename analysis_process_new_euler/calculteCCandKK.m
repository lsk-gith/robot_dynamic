function [CC]=calculteCCandKK(dh,Pc,m,Ic,g,R,P,G)
ddq=[0,0,0,0,0,0];
dq=[0,0,0,0,0,0];
for i=1:6
    dq(i)=1;
    [F,N] = outsideEquation(dh,dq,ddq,Pc,m,Ic,g,R,P);
    [CG]= insideEquation(dh,F,N,Pc,R,P);
    CG = CG';
    CG = CG - G;
    for j=1:6
       CC(j,i)=CG(j);
       dq(i)=0;
    end
end
end