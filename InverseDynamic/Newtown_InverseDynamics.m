function [taulist]= Newtown_InverseDynamics(q, dq, ddq, g,dh, m, Pc, Ic, f_tip)
%@input:q:joint angle for every link,  q(6X1).
%@input:dq:joint angle velocity for every link,  dq(6X1).
%@input:ddq:joint angle acculate for every link,  ddq(6X1).
%@input:g:9.81.
%@input:dh stanard for modified DH, dh_list = [alpha; a; d; theta];],  DH(4X6).
%@input:m:the mess of link,  m(1X6).
%@input:Pc sandard the mess center of link,example:Pc(1,1) standsrd for the center of x axis on first link, Pc(2,1)standsrd for the center of y axis on first link,Pc(3,1)standsrd for the center of z axis on first link.  Pc(3X6).
%@input:Ic sandard the inertia tensor of link, example:Ic(:,:,1)standard for the first link inertia tensor, and Ic(:,:,1) is a symmetry matrix, Ic(3X3X6).
%@input:f_tip:external force and torque,f_tip(1,:) stanard for the force, f_tip(2,:) stanard for the torque, f_tip(2X3)

%@output:taulist : the every link need torque,  taulist(6X1)

%R：3x3,旋转矩阵，P：3x6,后一连杆坐标系在前一连杆坐标系中的位置
%w：3x6,连杆角速度，dw：3x6,连杆角加速度，dv：3x6,连杆线加速度，dvc：3x6,连杆质心线加速度
%Ic:3x3x6,等同于inertia_tensor_list
%Pc:3x6, mass_center_list的转置
%F:3x6,各轴所受合力，N:3x6,各轴所受合力矩
%f:3x6,前一轴给后一轴的力，n:3x6,前一轴给后一轴的力矩
alpha = dh(1,:);
a = dh(2,:);
d = dh(3,:);
theta = dh(4,:);
Z=[0;0;1];
%转换矩阵建立

theta = theta + q';
T=zeros(4,4,6);R=zeros(3,3,6);P=zeros(3,6);
for i=1:6
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*cos(alpha(i))  cos(theta(i))*cos(alpha(i))  -sin(alpha(i))   -d(i)*sin(alpha(i))
              sin(theta(i))*sin(alpha(i))  cos(theta(i))*sin(alpha(i))  cos(alpha(i))    d(i)*cos(alpha(i))
              0                              0                              0                  1];
    R(:,:,i)=T(1:3,1:3,i);
    P(:,i)=T(1:3,4,i);
end

%运动学正向递推
w0 = zeros(3,1); dw0 = zeros(3,1);
dv0 = [0;0;g];
w = zeros(3,6); dw = zeros(3,6);
dv = zeros(3,6); dvc = zeros(3,6);
F = zeros(3,6); N = zeros(3,6);

%i = 0
w(:,1) = R(:,:,1)' * w0 + dq(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dq(1) * Z) + ddq(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
dvc(:,1) = cross(dw(:,1), Pc(:,1))+cross(w(:,1), cross(w(:,1), Pc(:,1))) + dv(:,1);
for i = 1:5
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dq(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dq(i+1) * Z)+ ddq(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
   dvc(:,i+1) = cross(dw(:,i+1), Pc(:,i+1)) + cross(w(:,i+1), cross(w(:,i+1), Pc(:,i+1))) + dv(:,i+1);
end

for i = 1:6
   F(:,i)=m(i)*dvc(:,i) ;
   N(:,i)=Ic(:,:,i) * dw(:,i) + cross(w(:,i), Ic(:,:,i) * w(:,i));
end

%动力学逆向递推
%先计算杆6的力和力矩
taulist = zeros(6,1);
f=zeros(3,6); n=zeros(3,6);

f(:,6) = F(:,6) + f_tip(1,:)';
n(:,6) = N(:,6) + f_tip(2,:)' + cross(Pc(:,6), F(:,6));
taulist(6) = n(:,6)' * Z;
%再计算杆5到1的力和力矩
for i=5:-1:1
   f(:,i) = R(:,:,i+1) * f(:,i+1) + F(:,i);
   n(:,i) = N(:,i) + R(:,:,i+1) * n(:,i+1) + cross(Pc(:,i), F(:,i))...
            + cross(P(:,i+1), R(:,:,i+1) * f(:,i+1));
   taulist(i) = n(:,i)' * Z;
   
end     

end