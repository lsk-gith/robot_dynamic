function [Jee,M,C,G] = Lagrange_InverseDynamics(dh,rob,Pc,Ic,g,m,q,dq)
%@input: dh stanard for modified DH, dh_list = [alpha; a; d; theta]; DH(4X6).
%@input: int rob the 0 standard for rotation, the 1 standard for translation. but also length(rob) equal size(dh,2);  rob(1X6).
%@input:Pc sandard the mess center of link, example: Pc(1,1) standsrd for the center of x axis on first link, Pc(2,1)standsrd for the center of y axis on first link,Pc(3,1)standsrd for the center of z axis on first link  Pc(3X6).
%@input:Ic sandard the inertia tensor of link,example:Ic(:,:,1)standard for the first link inertia tensor, and Ic(:,:,1) is a symmetry matrix   Ic(3X3X6).
%@input:g:9.81.
%@input:m:the mess of link,  m(1X6).
%@input:q:joint angle for every link,  q(1X6).
%@input:dq:joint angle velocity for every link,  dq(1X6).

%@output:Pee: the end position for robot, respect standard for the position of x y z axis.   Pee(3X1).
%@output:Ree: the end rotation for robot, it's is a rotation traslation matrix. Ree(3X3).
%@output:Jee: jacobians for robot contain linear and angle part, Jee(6X6).
%@output:M:M(q),  M(6X6).
%@output:C:C(q,dq),   C(6X6).
%@output:G:G(q),   G(1X6).
%% parameters
alpha = dh(1,:);
a = dh(2,:);
d = dh(3,:);
theta = dh(4,:);
N = size(dh,2); 
for i=1:N
   pc{i}=Pc(:,i);
end
for i=1:N
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

%% R
for i=1:N
    if rob(i)==0    % Rotation matrix of revolute joints
        R{i}=[cos(q(i)) -sin(q(i)) 0;
             sin(q(i))*Calpha(i) cos(q(i))*Calpha(i) -Salpha(i);
             sin(q(i))*Salpha(i) cos(q(i))*Salpha(i)  Calpha(i)];
    elseif rob(i)==1    % Rotation matrix of prismatic joints
        R{i}=[cos(theta(i)) -sin(theta(i)) 0;
             sin(theta(i))*Calpha(i) cos(theta(i))*Calpha(i) -Salpha(i);
             sin(theta(i))*Salpha(i) cos(theta(i))*Salpha(i)  Calpha(i)];
    end
end

%% P
for i=1:N
    if rob(i)==0    % Position vector of revolute joints
        p{i}=[a(i);-Salpha(i)*d(i);Calpha(i)*d(i)];
    elseif rob(i)==1    % Position vector of prismatic joints
        p{i}=[a(i);-Salpha(i)*q(i);Calpha(i)*q(i)];
    end
end


for i=1:N
    pic{i}=p{i}+R{i}*pc{i};     
end

poc{1}=pic{1};
for i=2:N
    nic=pic{i};      
    for j=i-1:-1:1
        nic=p{j}+R{j}*nic;
    end
    poc{i}=nic;
end
Pee = poc{N};

%% R0{i} stand for 0--> ith framework rotation transition matrix
R0{1}=R{1};
for i=2:N
    R0{i}=simplify(R0{i-1}*R{i});
end
Ree = R0{N};

%% linear part of jacobians
for i=1:N
    Jv{i}=simplify(jacobian(poc{i},q));
end

%% angular part of jacobians
if rob(1)==0
    Jo{1}=simplify([R0{1}(:,3),zeros(3,N-1)]);
elseif rob(1)==1
    Jo{1}=zeros(3,N);
end

for i=2:N
    if rob(i)==1
        Jo{i}=zeros(3,N);
    elseif rob(i)==0
        Jo{i}=Jo{i-1};
        Jo{i}(:,i)=R0{i}(:,3);
    end
end
Jee = [Jv{N};Jo{N}];

%% M
for i = 1:size(rob,2)
    In{i} = Ic(:,:,i);
end

for i=1:N
    Mass{i}=simplify(m(i)*Jv{i}.'*Jv{i}+Jo{i}.'*R0{i}*In{i}*R0{i}.'*Jo{i});
end

M=0;
for i=1:N
    M=simplify(Mass{i}+M);
end

%% C
for k=1:N
   for s=1:N
       c(1)=.5*((diff(M(k,s),q(1))+diff(M(k,1),q(s))-diff(M(1,s),q(k)))*dq(1));
      for i=2:N
       c(i)=.5*((diff(M(k,s),q(i))+diff(M(k,i),q(s))-diff(M(i,s),q(k)))*dq(i))+c(i-1);
      end
      C(k,s)=simplify(c(N));
   end
end

%% G
P(1)=m(1)*[0,0,g]*poc{1};
for i=2:N
    P(i)=P(i-1)+m(i)*[0,0,g]*poc{i};
end
P=simplify(P(N));
for i=1:N
    G(i,:)=simplify(diff(P,q(i)));
end

%% dynamic equations
%  M(q)*ddq + C(q,dq)dq + G(q) = u
end