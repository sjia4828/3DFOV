function [newx,newy,newz,C]=Ellipsoidvisualisation(Q,T_order,T_nor,X,Y,Z)



[R,D]=eig(Q);

% form the ellipsoid centorid (0,0,0)
% scaling para S
% S=1.5*T_order;
S=0.5;
[x,y,z] = ellipsoid(0,0,0,(1-D(1,1))*S,(1-D(2,2))*S,(1-D(3,3))*S,50);

% rotate data with orientation matrix R and Spatial position X,Y,Z
a = kron(R(:,1),x); 
b = kron(R(:,2),y); 
c = kron(R(:,3),z);

data = a+b+c; 
n = size(data,2);

newx = data(1:n,:)+X; 
newy = data(n+1:2*n,:)+Y;
newz = data(2*n+1:end,:)+Z;



C1=colormap('cool');
[m,n]=size(C1);
C_index=fix(m*double(T_nor))+1;
if C_index >= m
   C_index=m;
end
C_scale=C1(C_index,:);
C=cat(3,ones(size(x)).*C_scale(1),ones(size(x)).*C_scale(2),ones(size(x)).*C_scale(3));

end