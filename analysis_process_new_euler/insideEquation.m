function [tau] =insideEquation(dh,F,N,Pc,R,P)
number = size(dh,1);
Z = [0,0,1]';
R67 = eye(3,3);
f77 = zeros(3,1);%外力
n77 = zeros(3,1);%外力矩
P67 = zeros(3,1);
f(:,number) = R67 * f77 + F(:,number);
n(:,number) = N(:,number) + R67 * n77 + cross(Pc(:,6), F(:,6)) + cross(P67, R67 * f77);
tau(number) = n(:,number)' * Z;

for i=number-1:-1:1
   f(:,i) = R(:,:,i+1) * f(:,i+1) + F(:,i);
   n(:,i) = N(:,i) + R(:,:,i+1) * n(:,i+1) + cross(Pc(:,i), F(:,i))...
            + cross(P(:,i+1), R(:,:,i+1) * f(:,i+1));
   tau(i) = n(:,i)' * Z;
end

end