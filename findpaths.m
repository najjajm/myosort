%% FINDPATHS
% Finds paths between connected nodes in a graph
%
% SYNTAX
%   paths = findpaths(A)
%
% REQUIRED INPUTS
%   A (square array) - adjacency matrix 
%
% OPTIONAL INPUTS: none
%
% VARIABLE INPUTS: none
%
% OUTPUTS
%   paths (cell array) - list of connected nodes (separated by cells)
%
% EXAMPLE(S)
%
%   % generate a list of connected nodes
%   nNodes = 10;
%   x = nchoosek(1:nNodes,2);
%   y = rand(size(x,1),1); 
%   y(y>.9)=1; 
%   y(y~=1)=0;
%
%   % make adjacency matrix
%   A = false(nNodes);
%   A(sub2ind(size(A), x(y==1,1), x(y==1,2))) = true;
%
%   % find paths
%   paths = findpaths(A);
%
%
% IMPLEMENTATION
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% SEE ALSO:

% Authors: Najja Marshall
% Emails: njm2149@columbia.edu
% Dated: August 2017

function paths = findpaths(A, varargin)
%% Parse inputs

% initialize input parser
P = inputParser;
P.FunctionName = 'FINDPATHS';

% add required, optional, and parameter-value pair arguments
addRequired(P, 'A', @(x) size(x,1)==size(x,2));

% clear workspace (parser object retains the data while staying small)
parse(P, A, varargin{:});
clear ans varargin


%% Main loop

paths = {[]};
rowNo = 1;

while nnz(A) > 0
    
    % find connections from node in current row
    edges = find(A(rowNo,:));
    
    if ~isempty(edges)
        
        % add unique nodes to path list
        paths{end} = [rowNo, edges(~ismember(edges,paths{end}))];
        
        % remove connections from adjacency matrix
        A(rowNo,edges) = 0;
        
        % check connections to and from each node in the path
        nodeIdx = 2;
        while nodeIdx <= length(paths{end})
            
            % current node
            node = paths{end}(nodeIdx);
            
            % add connections from current node
            edges = find(A(node,:));
            paths{end} = [paths{end}, setdiff(edges,paths{end})];
            
            % remove from adjacency matrix
            A(node,edges) = 0;
            
            % add connections to current node
            edges = find(A(:,node));
            paths{end} = [paths{end}, setdiff(edges',paths{end})];
            
            % remove from adjacency matrix
            A(edges,node) = 0;
            
            nodeIdx = nodeIdx+1;
        end
        
        % add cell to path list if any edges remain
        if nnz(A) > 0
            paths = [paths; {[]}];
        end
    end
    
    % increment row to check in adjacency matrix
    rowNo = rowNo+1;
end

% sort list of nodes
paths = cellfun(@sort,paths,'uni',false);