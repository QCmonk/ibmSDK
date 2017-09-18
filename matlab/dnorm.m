function result = dnorm(N,A)
beep off
%DNORM Diamond norm
%
%   dnorm(J) returns the diamond norm (completely-bounded induced 1-norm) of a
%            Choi matrix J corresponding to the difference of two CP maps.
%
%   The normalization of J must be such that partial prace of the first
%   half of J always yields the identity for trace preserving
%   superoperators.
% 
%  Copyright 2014 Raytheon BBN Technologies
%  
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%  
%      http://www.apache.org/licenses/LICENSE-2.0
%  
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%
J = N-A;
dimsq = size(J,1);
dim = round(sqrt(dimsq));
Jh = (J'+J)/2;

cvx_begin sdp
  cvx_precision best

  variable W(dimsq,dimsq) hermitian;
  variable rho(dim,dim) hermitian;

  maximize( trace(Jh'*W) );
  subject to
    W >= 0;
    rho >= 0;
    W - kron(eye(dim),rho) <= 0;
    trace(rho) == 1;

cvx_end
result = cvx_optval;
end

