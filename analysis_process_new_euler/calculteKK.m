function [KK]=calculteKK(dh,G,Pc,m,Ic,g,R,P,CC)
% global nn ll;
ddq=[0,0,0,0,0,0];
dq = [0,0,0,0,0,0];
if isa(m(1),'double')
    KK = [];
end
i = 1;
for nn = 1:5
    for ll=(nn+1):6
        dq(nn)=1;dq(ll)=1;
        [F,N] = outsideEquation(dh,dq,ddq,Pc,m,Ic,g,R,P);
        [CK]= insideEquation(dh,F,N,Pc,R,P);
        CK=CK';
        CK=CK-G-CC(:,nn)-CC(:,ll);
        KK(:,i)= CK;
        dq(nn)=0;dq(ll)=0;
        i = i + 1;
    end
end
end