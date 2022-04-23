function [R,P] = compute_frame_transform(dh,thetalist)
%%dh为改进坐标系下n乘以4的矩阵n为自由度个数,thetalist为theta的数值
[n,~] = size(dh);
alpha = dh(:,1);
a = dh(:,2);
d = dh(:,3);
theta = dh(:,4);
theta = theta + thetalist;

for i=1:n
    if alpha(i)==0
        Salpha(i)=0;
        Calpha(i)=1;
    elseif abs(alpha(i))==pi/2
        Calpha(i)=0;
        if alpha(i)==pi/2
            Salpha(i)=1;
        elseif alpha(i)==-pi/2
            Salpha(i)=-1;
        end
    elseif abs(alpha(i))==pi
        Salpha(i)=0;
        Calpha(i)=-1;
    else
        Salpha(i)=sin(alpha(i));
        Calpha(i)=cos(alpha(i));
    end
end


for i=1:n
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*Calpha(i)  cos(theta(i))*Calpha(i)  -Salpha(i)   -d(i)*Salpha(i)
              sin(theta(i))*Salpha(i)  cos(theta(i))*Salpha(i)  Calpha(i)    d(i)*Calpha(i)
              0                              0                              0                  1];
              R(:,:,i)=(T(1:3,1:3,i));
              P(:,i)=(T(1:3,4,i));
end
end
