function [dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz]=gradient_tps3_3D(x,y,z,X,Y,Z,U,V,W)%coef_matu,coef_matv,coef_matw)
% This function is called by gradient_rbf_3D.m
% x,y,z: Current point where gradient is being evaluated (all are single
%        values)
% X,Y,Z,U,V,W: Mesh around the point where the gradient is being evaluated.
%              All meshes should be column vectors

%%% Compute the coefficient matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the polynomial (P) and RBF-phi (B) matrices
P = zeros(length(X),4);
B = zeros(length(X),length(Y));
% Form the polynomial (P) matrix
for i=1:length(X)
    P(i,:)=[1,X(i),Y(i),Z(i)];%,X(i)^2,X(i)*Y(i),Y(i)*Y(i)];
end
% Form the phi matrix
for i=1:length(X)
    for j=1:length(X)% index for columns
       dummy=((X(i)-X(j))^2 +(Y(i)-Y(j))^2 + (Z(i)-Z(j))^2);
        B(i,j)=dummy*log(sqrt(dummy+eps));
    end
end

G=[B,P;
transpose(P),zeros(size(P,2),size(P,2))];
% Get the coefficient matrix
coef_matu=G\[U;zeros(size(P,2),1)];
coef_matv=G\[V;zeros(size(P,2),1)];
coef_matw=G\[W;zeros(size(P,2),1)];

%%% Compute the gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the interpolated values
dudx_interp=0;dvdx_interp=0;dwdx_interp=0;dudy_interp=0;dvdy_interp=0;dwdy_interp=0;
dudz_interp=0;dvdz_interp=0;dwdz_interp=0;

% For each mesh point, compute the interpolated gradient
for k = 1:length(X)
    dummy = ((x-X(k))^2 +(y-Y(k))^2+(z-Z(k))^2);
    dummy1 = (1+2*log(sqrt(dummy+eps)));
    dudx_interp=(x-X(k))*dummy1*coef_matu(k)+dudx_interp;
    dudy_interp=(y-Y(k))*dummy1*coef_matu(k)+dudy_interp;
    dudz_interp=(z-Z(k))*dummy1*coef_matu(k)+dudz_interp;
    dvdx_interp=(x-X(k))*dummy1*coef_matv(k)+dvdx_interp;
    dvdy_interp=(y-Y(k))*dummy1*coef_matv(k)+dvdy_interp;
    dvdz_interp=(z-Z(k))*dummy1*coef_matv(k)+dvdz_interp;
    dwdx_interp=(x-X(k))*dummy1*coef_matw(k)+dwdx_interp;
    dwdy_interp=(y-Y(k))*dummy1*coef_matw(k)+dwdy_interp;
    dwdz_interp=(z-Z(k))*dummy1*coef_matw(k)+dwdz_interp;
end
% Compute the final gradient value based on the interpolated values and the
% coefficient values
Pddx=[0,1,0,0];Pddy=[0,0,1,0];Pddz=[0,0,0,1];%,x^2,x*y,y^2];%,x^2,x*y,y^2];
dudx=dudx_interp+Pddx*coef_matu(end-size(Pddx,2)+1:end);
dvdx=dvdx_interp+Pddx*coef_matv(end-size(Pddx,2)+1:end);  
dwdx=dwdx_interp+Pddx*coef_matw(end-size(Pddx,2)+1:end); 
dudy=dudy_interp+Pddy*coef_matu(end-size(Pddy,2)+1:end);  
dvdy=dvdy_interp+Pddy*coef_matv(end-size(Pddy,2)+1:end);  
dwdy=dwdy_interp+Pddy*coef_matw(end-size(Pddy,2)+1:end);
dudz=dudz_interp+Pddz*coef_matu(end-size(Pddy,2)+1:end);  
dvdz=dvdz_interp+Pddz*coef_matv(end-size(Pddy,2)+1:end);  
dwdz=dwdz_interp+Pddz*coef_matw(end-size(Pddy,2)+1:end);

end


            