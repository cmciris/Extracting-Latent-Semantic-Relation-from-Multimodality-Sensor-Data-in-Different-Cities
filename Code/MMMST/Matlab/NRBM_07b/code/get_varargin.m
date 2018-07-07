function val = get_varargin(varargin,label,default_val)
% function val = get_varargin(varargin,label,default_val)
% returns the value of a parameter from varargin (cell array) 
% Example : 
%           lambda = get_varargin(varargin,'lambda',0.01)
%           va = {'alpha',0.1,'beta',0.2,'delta',0.3,'gamma',0.4}
%           delta = get_varargin(va,'delta',0.01)
%           -> delta = 0.3

%
% Copyright (C) 2012, by Trinh-Minh-Tri Do minhtrido@gmail.com
%
%   This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  if nargin<3
    val = [];
  else
    val = default_val;
  end

  for i=1:floor(numel(varargin)/2)
    str = varargin{(i-1)*2+1};
    if isstr(str)&&strcmpi(str,label)
      val = varargin{i*2};
    end
  end

