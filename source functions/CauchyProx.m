function u = CauchyProx(x, gamma, mu)  %,a,b,c,d
% Cauchy proximal operator
%
% u = arg min_u (||u - x||_2)^2/(2*mu) - log(gamma/(gamma^2 + u^2)) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
%      x      - Input for the proximal operator
%      gamma  - Cauchy scale parameter
%      mu     - FB algorithm step size. It can be chosen in relation with
%               Lipschitz constant. 
% OUTPUT :
%      u - the result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  a    = 1;
%  b    = -x;
%  c    = gamma^2+2*mu;
%  d    = -gamma^2*x;
%  p    = c/a - b.*b/a/a/3. 
%  q    = (-2.*b.*b.*b/a/a/a + 9.*b.*c/a/a - 27.*d/a) / 27. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% 
% Copyright (C) Oktay Karakus,PhD
%               and
%               Alin Achim, PhD
% University of Bristol, UK
% o.karakus@bristol.ac.uk
% alin.achim@bristol.ac.uk
% April 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REFERENCE
%
% [1] O Karakus, P Mayo, and A Achim. "Convergence Guarantees for 
%     Non-Convex Optimisation with Cauchy-Based Penalties"
%       IEEE Transactions on Signal Processing, 2020.
%
% [2] T Wan, N Canagarajah, and A Achim. "Segmentation of noisy colour 
%     images using Cauchy distribution in the complex wavelet domain." 
%     IET Image Processing 5.2 (2011): 159-170.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = gamma.^2 + 2*mu - x.^2/3;

q = 2*x.^3/27 + gamma.^2.*x - (x/3)*(gamma^2 + 2*mu);

DD = p.^3/27 + q.^2/4;

s = (abs(q/2 + sqrt(DD)).^(1./3.)).*sign(q/2 + sqrt(DD));

t = (abs(q/2 - sqrt(DD)).^(1./3.)).*sign(q/2 - sqrt(DD));

u = x/3 + s + t;