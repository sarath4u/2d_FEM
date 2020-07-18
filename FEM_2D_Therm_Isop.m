%====================================================================================%
% A FINITE ELEMENT SOLUTION ON A STRUCTURED MESH FOR 2D STEADY STATE HEAT CONDUCTION %
%====================================================================================%

Lx = 10;                                                  % Length in the x direction
Ly = 10;                                                  % Length in the y direction
nx = 51;                                                    % no of divisions in x
ny = 51;                                                    % no of divisions in y

% at this point the divisions are uniform

xcord = linspace(0,Lx,nx);                                 % x co-ordinates
ycord = linspace(0,Ly,ny);                                 % y co-ordinates
nodalcords = combvec(xcord,ycord)';                        % nodal coordinates
n_elmnts = (nx-1)*(ny-1);                                  % no of elements
elements = zeros(n_elmnts, 4);                             % elements
force = 1;

% =================================================================================== %

% 7___8___9
% |(3)|(4)|
% 4___5___6
% |(1)|(2)|
% 1___2___3  *Elemenets and Nodes*

%%% To be adapted to read from the geometry file

for j = 1 : 1 : n_elmnts
    q = fix((j - 1) /(nx - 1));
    elements(j, 1) = j + q;
    elements(j, 2) = j + q + 1; 
    elements(j, 3) = j + q + nx + 1;
    elements(j, 4) = j + q + nx;
end

% Nodal Vector
nodevec = (1 : nx * ny)';
% ==================================================================================== %

%%% Boundaries

%  ___B3___
% |        |
% B4       B2    
% |        | 
% |___B1___|

% Change here for different BCs

B1 = 1:nx;
B2 = nx : nx : ny * nx;
B3 = ny * nx - nx + 1 : ny * nx;
B4 = 1: nx : ny*nx;

% 1. Dirichlet BCs

bndrynds = unique([B1 B2 B3 B4])';      % nodes on Dirichlet boundary
%bndrynds = unique([B1 B2])';           % nodes on Dirichlet boundary

cnt = 1;
nonbndrynds = zeros(size(nodevec, 1) - size(bndrynds, 1), 1);  % nodes NOT on the Dirichlet boundary
for nn = 1 : nx*ny
    res = ismember(nn, bndrynds);
    if res == 0
        nonbndrynds(cnt) = nn;
        cnt = cnt + 1;
    end
end

% 2. Nuemann BCs

trn = 1;                                            % traction/flux

%Nubound = unique([B3 B4])';                        % Nodes on the boundary
Nubound = unique([])';

Nuele = ones(size(elements));
eleTracMat = zeros(size(elements));

% ======================================================================================== %

elmatpos = (combvec (1 : 4, 1 : 4))';

k = diag(ones(1,2));                        % element stiffness matrix
wg = 1;                                     % Weights
gp = [-0.5774   -0.5774;
       0.5774   -0.5774; 
       0.5774    0.5774;                
      -0.5774    0.5774];                   % Gauss Points
K = zeros(nx * ny);                         % Stiffness Matrix for assembly
F = zeros(nx * ny, 1);                      % Force vector for assembly
Traction = zeros(nx * ny, 1);

% ========================================================================================= %

for i = 1 : 1 : n_elmnts
    el = elements(i, :);

% Isoparametric Shape Functions

%  (-1,1)4______3(1,1)
%        |      |
%        |      |    
%        |      | 
% (-1,-1)1______2(1,-1)

% N1 = (1+eps)(1-hota)/4 %
% N2 = (1+eps)(1+hota)/4 %
% N3 = (1-eps)(1+hota)/4 %
% N4 = (1-eps)(1-hota)/4 %

% Je = [ (x13 + x24)/4 + (x21 + x43)hota/4   (y13 + y24)/4 + (y21 + y43)hota/4
%        (x21 + x34)/4 + (x21 + x43)eps /4   (y21 + y34)/4 + (y21 + y43)eps /4 ]
% Je is not constant 

    x1 = nodalcords(el(1),1);
    y1 = nodalcords(el(1),2);
    x2 = nodalcords(el(2),1);
    y2 = nodalcords(el(2),2);
    x3 = nodalcords(el(3),1);
    y3 = nodalcords(el(3),2);
    x4 = nodalcords(el(4),1);
    y4 = nodalcords(el(4),2);
    
    x13 = x1 - x3;
    x24 = x2 - x4;
    x21 = x2 - x1;
    x43 = x4 - x3;
    
    y13 = y1 - y3;
    y24 = y2 - y4;
    y21 = y2 - y1;
    y43 = y4 - y3;
    
    for mm = 1 : 4
        res = ismember(elements(i, mm), Nubound);
        if res == 0
            Nuele(i, mm) = 0;  
        end
    end
    % ================================================================================ %
    %%% Applying traction
    
    if Nuele(i, 1) && Nuele(i, 2)
        eleTracMat(i, :) = eleTracMat(i, :) + trn * [1 0 0 1] * x13 / 2;
    end
    if Nuele(i, 2) && Nuele(i, 3)
        eleTracMat(i, :) = eleTracMat(i, :) + trn * [1 1 0 0] * y21 / 2;
    end
    if Nuele(i, 3) && Nuele(i, 4)
        eleTracMat(i, :) = eleTracMat(i, :) + trn * [0 1 1 0] * (x24 - x43) / 2;
    end
    if Nuele(i, 1) && Nuele(i, 4)
        eleTracMat(i, :) = eleTracMat(i, :) + trn * [0 0 1 1] * (y21 - y13) / 2;
    end
    eleTrac = eleTracMat(i, :)';
    
    %%% Traction vector
    
    for jj = 1 : size(elements, 1)
        Traction(elements(jj, :), 1) = Traction(elements(jj, :), 1) + eleTrac(:, 1);
    end
    % =============================================================================== %
    
    for pt = 1 : 1 : size(gp)
        N = [(1 + gp(pt, 1)) * (1 - gp(pt, 2)) / 4 ...
             (1 + gp(pt, 1)) * (1 + gp(pt, 2)) / 4 ...
             (1 - gp(pt, 1)) * (1 + gp(pt, 2)) / 4 ...
             (1 - gp(pt, 1)) * (1 - gp(pt, 2)) / 4 ];
        
        Je = [  (x13 + x24)/4 + (x21 + x43) * gp(pt, 2)/4 ...
                (y13 + y24)/4 + (y21 + y43) * gp(pt, 2)/4; ... 
                (x21 - x43)/4 + (x21 + x43) * gp(pt, 1)/4 ...
                (y21 - y43)/4 + (y21 + y43) * gp(pt, 1)/4 ];
            
        Jin = inv(Je);

% Be = J^(-1)/4 *  [ 1-hota 1+hota -1-hota -1+hota
%                   -1-eps  1+eps   1-eps  -1+eps ]    not a constant

        Be = (Jin / 4) * [ 1 - gp(pt, 2) 1 + gp(pt, 2) -1 - gp(pt, 2) -1 + gp(pt, 2); ...
                          -1 - gp(pt, 1) 1 + gp(pt, 1)  1 - gp(pt, 1) -1 + gp(pt, 1)];

% *******   ********   ******** %
%  Stiffness Matrix generation  %
% *******   ********   ******** %                  
              
%      __ 
%      \ 
% Ke = /_ wg.Be'.k.Be.|Je|                      element stiffness matrix
       
        Ke = wg * Be' * k * Be * det(Je);
%      __ 
%      \ 
% fe = /_ wg.N.|Je|                             element force vector
        
        fe = wg * N' * det(Je);
    end
    
    matpos = (combvec(el,el))';
    
%   Assembly of the stiffness matrix
%  
    for j = 1 : 1 : size(matpos,1)
        K(matpos(j, 1), matpos(j,2)) = K(matpos(j, 1), matpos(j,2)) ...
            + Ke(elmatpos(j, 1),elmatpos(j,2));
    end
    
% *******   ********   ******** %
%    Force vector generation    %
% *******   ********   ******** %

    for ii = 1 : size(elements, 1)
        F(elements(ii, :), 1) = F(elements(ii, :), 1) + fe(:, 1);
    end
    
    
end

%%% Applying Dirichlet Boundary Condition

% Reducing the Stiffness matrix and the Force vector

for jj = size(bndrynds, 1): -1 :1
   K(bndrynds(jj), :) = [];
   K(:, bndrynds(jj)) = [];
   F(bndrynds(jj)) = [];
   Traction(bndrynds(jj)) = [];
end

%%% nodal solution

F = F + Traction;
u = K\F;

% Filling the dirichlet condition in the solution vector

temp = zeros(nx*ny,1);
count = 1;
for kk = 1 : size(temp, 1)
    res = ismember(kk, nonbndrynds);
    if res == 1
        temp(kk, 1) = u(count, 1);
        count = count + 1;
    else
        temp(kk, 1) = 0.0;                  %% Dirichlet condition
    end
end

% Re-assembling in the Global solution Matrix

T = zeros(nx,ny);
row = 1;
col = 1;
count = 1;
while row <= nx
    T(row, col) = temp(count,1);
    if col == ny
        row = row + 1;
        col = 0;
    end
    col = col + 1;
    count = count + 1;
end

% Plotting the result

surf(xcord,ycord,T)
colormap jet

%========================================================================================%
