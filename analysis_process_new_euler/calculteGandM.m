function [G,M]=calculteGandM(dh,Pc,m,Ic,g,R,P)
dq=[0,0,0,0,0,0]';
ddq=[0,0,0,0,0,0]';
[F,N] = outsideEquation(dh,dq,ddq,Pc,m,Ic,g,R,P);
[G]= insideEquation(dh,F,N,Pc,R,P);
G = G';
for i=1:6
    ddq(i)=1;
    [F,N] = outsideEquation(dh,dq,ddq,Pc,m,Ic,g,R,P);
    [MG]= insideEquation(dh,F,N,Pc,R,P);
    MG = MG';
    MG = MG - G;
    for j=1:6
        M(j,i)=MG(j);
        ddq(i)=0;
    end
end
end