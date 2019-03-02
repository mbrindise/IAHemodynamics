function [V1,V2,V3,D,A] = KLDs3v(U1,U2,U3)
% Perform Karhunen-Loeve / Proper-Orthogonal decomposition on the two 
% dimensional vector signal U_i integrated on an associated GRIDTYPE.
% The gridtype for 3D is limited to 'even' only! For 2D, other grid types
% are available. These grid types can be ported into the 3D, but this has
% not been completed

% U to be decomposed needs to be reordered as a K by NT array, where 
% If U is originally 2D, reorder U as a single strip first.  
% 
% V1: Eigenvectors corresponding to U velocity
% V2: Eigenvectors corresponding to V velocity
%  D: Diagonal matrix of eigenvalues
%  A: Time Coefficients

[K1,K2,K3,NT] = size(U1);
if size(U2)~=[K1,K2,K3,NT]
    error('U1 and U2 must be the same sizes')
end

% Only an EVEN grid type may be used for the 3D code
w_yj = ones(K1,K2,K3);
c1_j     = ones(K1,1);
c1_j(1)  = 0.5;
c1_j(K1) = 0.5;
c2_j     = ones(K2,1);
c2_j(1)  = 0.5;
c2_j(K2) = 0.5;
c3_k     = ones(1,1,K3);
c3_k(1,1,1)   = 0.5;
c3_k(1,1,K3)   = 0.5;
c_jp = c1_j*c2_j.';
c_j = c_jp.*c3_k;
clear c1_j c2_j c3_k;

%unfold 2D fields into 1D strips
U1    = reshape(permute(U1,[2,1,3,4]),K1*K2*K3,NT,1);
U2    = reshape(permute(U2,[2,1,3,4]),K1*K2*K3,NT,1);
U3    = reshape(permute(U3,[2,1,3,4]),K1*K2*K3,NT,1);
w_yj = reshape(permute(w_yj,[2,1,3]),K1*K2*K3,1);
c_j  = reshape(permute( c_j,[2,1,3]),K1*K2*K3,1);

%concatenate fields
[V,D,A] = KLDs([U1;U2;U3],'custom',[w_yj;w_yj;w_yj],[c_j;c_j;c_j]);
D = D*NT; % Scale the eigenvalues properly

%need to reorder 1-D eigenvectors into 3D eigenmodes
V1 = V(      1:  K1*K2*K3,:); % Obtain U values back into own matrix
V2 = V(K1*K2*K3+1:2*K1*K2*K3,:); % Obtain V values back into own matrix
V3 = V(2*K1*K2*K3+1:3*K1*K2*K3,:); % Obtain W values back into own matrix
V1 = permute(reshape(V1,K2,K1,K3,NT),[2,1,3,4]); % Return U matrix to 4D (3D-Space, 1D-Time)
V2 = permute(reshape(V2,K2,K1,K3,NT),[2,1,3,4]); % Return V matrix to 4D (3D-Space, 1D-Time)
V3 = permute(reshape(V3,K2,K1,K3,NT),[2,1,3,4]); % Return V matrix to 4D (3D-Space, 1D-Time)

