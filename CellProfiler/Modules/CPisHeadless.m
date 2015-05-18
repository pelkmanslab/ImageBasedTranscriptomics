function isHeadless = CPisHeadless()
%CPISHEADLESS Return true if no MATLAB desktop is present.
%
% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Authors:
%   yauhen.yakimovich@uzh.ch
%
%   Use this function to avoid memory matlab plotting in "headless" mode,
%   e.g. in a command line or when running on HPC-cluster.

%
% $Revision: 1 $

isHeadless = ~usejava('desktop');

end

