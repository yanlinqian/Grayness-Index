function H = hist2w(X,W,y,z,C)
%% Weighted 2D histogram
% function H = hist2w(X,W,y,z,C)
%
% Inputs: X - Nx2 data
%         W - Nx1 weights
%         y - vector of M(1) bin centers along 1st dim of X
%         z - vector of M(2) bin centers along 2nd dim of X
%         C - M(1)xM(2)x3 RGB coloring or Kx3 colormap (default: gray)
% Output: H - M(1)xM(2) bar heights
% 
% Grid vectors (y, z) must have linear spacing
% If no output arguments specified, plots 2D colored histogram.
%
% Johannes Traa - March 2014

%% check inputs
if nargin < 5 || isempty(C); C = linspace(0.5,1,200)'*ones(1,3); end

%% prepare
M = [length(y) length(z)]; % # bins
w = [y(2)-y(1) z(2)-z(1)]; % bar widths

%% calculate weighted histogram
H = zeros(M(1),M(2),1); % bar heights

% boundary indicators
yu = (X(:,1) <= y(2)      - w(1)/2); % under y range
ya = (X(:,1) >= y(M(1)-1) + w(1)/2); % above y range
zu = (X(:,2) <= z(2)      - w(2)/2); % under z range
za = (X(:,2) >= z(M(2)-1) + w(2)/2); % above z range

% corners
H(1,1)       = sum(W(yu & zu)); % corner (-,-)
H(M(1),1)    = sum(W(ya & zu)); % corner (+,-)
H(1,M(2))    = sum(W(yu & za)); % corner (-,+)
H(M(1),M(2)) = sum(W(ya & za)); % corner (+,+)

% sides
for i=2:M(2)-1 % left, right
  sh = ((X(:,2) >= z(i) - w(2)/2) & (X(:,2) <= z(i) + w(2)/2)); % horizontal slice
  H(1,i)    = sum(W(sh & yu)); % left
  H(M(1),i) = sum(W(sh & ya)); % right
end
for j=2:M(1)-1 % bottom, top
  sv = ((X(:,1) >= y(j) - w(1)/2) & (X(:,1) <= y(j) + w(1)/2)); % vertical slice
  H(j,1)    = sum(W(sv & zu)); % bottom
  H(j,M(2)) = sum(W(sv & za)); % top
end

% interior
for i=2:M(2)-1
  for j=2:M(1)-1
    sh = ((X(:,2) >= z(i) - w(2)/2) & (X(:,2) <= z(i) + w(2)/2)); % horizontal slice
    sv = ((X(:,1) >= y(j) - w(1)/2) & (X(:,1) <= y(j) + w(1)/2)); % vertical slice
    H(j,i) = sum(W(sh & sv));
  end
end

%% plot
if nargout == 0
  figure
  bar3c(H,y,z,[],C);
  axis([min(y) max(y) min(z) max(z) 0 1.1*max(H(:))])
  clear H B
  return
end



end








function H = bar3c(Z,x,y,w,C)
%% 3D bar plot from matrix (colored by height)
% function H = bar3c(Z,x,y,w,C)
%
% Input: Z - NxM matrix
%        x - 1xN grid (default: 1:N)
%        y - 1xM grid (default: 1:M)
%        w - 1x2 bar half-widths (default: full half-width)
%        C - Kx3 colormap or 
%            NxMx3 bar colors (default: colormap(gray))
% Output: H - matrix of surf handles
%
% bar3c creates a 3D colored bar plot by plotting each bar as an individual
% surface object with its own color.
%
% Example (Gaussian bar plot):
%   M = 50; % grid resolution
%   x = linspace(-3,3,M); % x-grid
%   y = linspace(-3,3,M); % y-grid
%   C = [0.3 -0.2; -0.2 0.6]; % covariance
%   [X,Y] = meshgrid(x,y);
%   XY = [X(:) Y(:)]; % grid pairs
%   Z = 1/sqrt(det(2*pi*C))*exp(-1/2*sum((XY/C).*XY,2));
%   Z = reshape(Z,[M,M]);
%   figure
%   bar3c(Z,x,y,[],jet(50))
%
% Johannes Traa - July 2013

[N,M] = size(Z);

%% check inputs
if nargin < 2 || isempty(x); x = 1:N; end
if nargin < 3 || isempty(y); y = 1:M; end
if nargin < 4 || isempty(w); w = [x(2)-x(1) y(2)-y(1)]/2; end
if nargin < 5 || isempty(C); C = colormap(gray); end

%% x and y info for surfs
Hx = cell(1,N); % x values of surf
Hy = cell(1,M); % y values of surf
for i=1:N
  X = x(i) + [-w(1) -w(1) w(1) w(1)];
  Hx{i} = repmat(X,[4,1]);
end
for j=1:M
  Y = y(j) + [-w(2) -w(2) w(2) w(2)]';
  Hy{j} = repmat(Y,[1,4]);
end

%% bar coloring
if ismatrix(C)
  r = ceil((Z-min(Z(:)))/(max(Z(:))-min(Z(:)))*size(C,1));
  r(r(:)==0) = 1;
end

%% make bars
hold on

H = zeros(N,M); % surf handles
Hz = zeros(4,4); %+min(Z(:)); % z values of surf
for i=1:N
  for j=1:M
    Hz(2:3,2:3) = ones(2,2)*Z(i,j);
    
    % bar color
    if ismatrix(C); C_ij = C(r(i,j),:);
    else            C_ij = reshape(C(i,j,:),[1,3]);
    end
    
    % plot bar
    H(i,j) = surf(Hx{i},Hy{j},Hz,'FaceColor',C_ij,'EdgeColor','k');
  end
end

view(-50,40)
axis vis3d tight
grid on
drawnow

%% output
if nargout < 1; clear H; end


end