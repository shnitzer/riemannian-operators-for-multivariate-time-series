function [ S ] = FixedGeodes_eff( A,B,t,Dim )
% Input:
%   A,B - PSD matrices of rank Dim
%   0<t<1 - Desired point along the geodesic
%   Dim - rank of A or B.
% Output:
%   S - the point t along the geodesic streching from A to B.

[U1,S1] = eigs(A,Dim);
[U2,S2]  = eigs(B,Dim);

VA = U1(:,1:Dim);
VB = U2(:,1:Dim);

[OA,SAB,OB] = svd(VA.'*VB);
UA          = VA * OA;
UB          = VB * OB;
theta       = acos(diag(SAB));
tmp         = UB * pinv(diag(sin(theta)));
X           = ( speye(size(A)) * tmp - UA * (UA.' * tmp) );
U           = UA * diag(cos(theta*t)) + X * diag(sin(theta*t));

RB2         = OB.' * S2(1:Dim,1:Dim)       * OB; % UB.' * B * UB;
RA          = OA.' * sqrt(S1(1:Dim,1:Dim)) * OA;
RAm1        = OA.' / sqrt(S1(1:Dim,1:Dim)) * OA;
R2          = RA * expm( t*logm(RAm1 * RB2 * RAm1) ) * RA;
S           = U * R2 * U.';

end

