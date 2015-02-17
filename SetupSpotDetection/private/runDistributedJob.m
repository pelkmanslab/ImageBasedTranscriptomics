function [strBrutusSubmission] = runDistributedJob(strBrutus,strFunction,varargin)
% Function to run distributed jobs; Will recognize, whether Brutus is used
% or not and allows to execute the function outside of brutus. Thus the
% same code can be used (which is useful for debugging - especially if a
% small minority of jobs crashes for unexpected reasons);

% creates brutus commands for job submission from matlab
% strBrutus     Brutus specific string.: eg.:
%               'bsub -W 8:00 -R "rusage[mem=2000]"'
% strFunction   Name of function that should be executed
% varargin      various input elements for strFunction
%               (arguments can be string or number)
%
% also see: MakeMatlabJobSubmisison for an earlier version, which did not
% discriminate between brutus and local computers.


if checkIfIAmWorkingOnACluster == true;   % when working on brutus, use brutus specific settings
    
    NEWLINE = sprintf('\n');
    
    strBrutusSubmission = [strBrutus '<<EOS;' NEWLINE...
        'bmatlab <<M_PROG;' NEWLINE...
        strFunction '('];
    
    for j=1:length(varargin)
        if isnumeric(varargin{j})
            strToAdd = num2str(varargin{j});
        elseif ischar(varargin{j})
            strToAdd = varargin{j};
        else
            error(['Do not recognize format of ' varargin{j}])
        end
        if j==1
            strBrutusSubmission = [strBrutusSubmission '''' strToAdd]; %#ok<AGROW>
        else
            strBrutusSubmission = [strBrutusSubmission ''', ''' strToAdd]; %#ok<AGROW>
        end
    end
    
    strBrutusSubmission = [strBrutusSubmission ''')' NEWLINE...
        'M_PROG' NEWLINE...
        'EOS'];
    
    % execute
    system(strBrutusSubmission);
    
else  % run locally
    
    strLocalSubmission = [strFunction '('];
    
    for j=1:length(varargin)
        if isnumeric(varargin{j})
            strToAdd = num2str(varargin{j});
        elseif ischar(varargin{j})
            strToAdd = varargin{j};
        else
            error(['Do not recognize format of ' varargin{j}])
        end
        if j==1
            strLocalSubmission = [strLocalSubmission '''' strToAdd]; %#ok<AGROW>
        else
            strLocalSubmission = [strLocalSubmission ''', ''' strToAdd]; %#ok<AGROW>
        end
    end
    
    strLocalSubmission = [strLocalSubmission ''')'];
    
    % execute
    eval(strLocalSubmission);
end


end