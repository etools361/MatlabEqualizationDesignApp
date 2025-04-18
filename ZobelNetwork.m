%--------------------------------------------------------------------------
% Edited by bbl
% Date: 2023-07-14(yyyy-mm-dd)
% T线圈maximally flat amplitude (MFA) or maximally flat envelope delay (MFED) 
%--------------------------------------------------------------------------
% syms w K w0 real
% syms w real
syms K L1 C1 real
syms s
% delta = 6;
% w0 = 2*pi*0.1e9;
% K = 10^(delta/20);
R1 = 50*(K-1)/(K+1);
R2 = 100*K/(K^2-1);
% s  = -1i*w;
L1 = 50^2*C1;
% C1 = 1/(s^2*L1);
% syms s C1 R1 R2 L1 real

A=1/(s*C1);
B=R1;
CC=R1;
D=R2+s*L1;
E=50;

% KK = B*CC*A + B*D*A + B*E*A + D*CC*A - E*CC*A - E*E*A - E*E*B - E*E*CC
% simplify(KK)

D0 = B*CC*A + B*D*A + B*E*A + B*E*CC + D*CC*A + D*E*A + D*E*B + D*E*CC;
D12   = D*A + D*B + D*CC + B*CC;
H = D12/D0*E
H2 = simplify(H)
% M = 50*C1*(K-1);
% H2 = (s + 1/M)/(s + K/M);
% a = 1/(50*C1*(K-1))
% b = K*a


%----------------------------negative slope -----------------------
syms K L1 C1 real
syms s
R1 = 50*(K-1)/(K+1);
R2 = 100*K/(K^2-1);
L1 = 50^2*C1;

A=s*L1;
B=R1;
CC=R1;
D=R2+1/(s*C1);
E=50;

% KK = B*CC*A + B*D*A + B*E*A + D*CC*A - E*CC*A - E*E*A - E*E*B - E*E*CC
% simplify(KK)

D0 = B*CC*A + B*D*A + B*E*A + B*E*CC + D*CC*A + D*E*A + D*E*B + D*E*CC;
D12   = D*A + D*B + D*CC + B*CC;
H = D12/D0*E
H2 = simplify(H)

M = (K-1)/(50*C1)
H2 = (1/K)*(s + M)/(s + M/K)
% a = M
% b = M/K --> K = a/b
% a = b*K;
% a = w0*(K-1);
% b = w0*(K-1)/K
% w0 = 1/(50*C1);






