%% FINDPATHS find paths
% Identifies all sets of connected nodes in an undirected graph.
%
% SYNTAX
%   [grp,solo] = findpaths(A)
%
% REQUIRED INPUTS
%   A (square array) - adjacency matrix for an undirected graph. Where
%       a_{ij} is true or 1 indicates a connection between nodes i and j.
%
% OPTIONAL INPUTS: none
%
% VARIABLE INPUTS: none
%
% OUTPUTS
%   grp (cell array) - list of connected nodes
%   solo (cell array) - list of unconnected nodes
%
% EXAMPLE
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
%   % find connected nodes
%   grp = findpaths(A);
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

function [grp,solo] = findpaths(A, varargin)
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

grp = {[]};
rowNo = 1;

while nnz(A) > 0
    
    % find connections from node in current row
    edges = find(A(rowNo,:));
    
    if ~isempty(edges)
        
        % add unique nodes to path list
        grp{end} = [rowNo, edges(~ismember(edges,grp{end}))];
        
        % remove connections from adjacency matrix
        A(rowNo,edges) = 0;
        
        % check connections to and from each node in the path
        nodeIdx = 2;
        while nodeIdx <= length(grp{end})
            
            % current node
            node = grp{end}(nodeIdx);
            
            % add connections from current node
            edges = find(A(node,:));
            grp{end} = [grp{end}, setdiff(edges,grp{end})];
            
            % remove from adjacency matrix
            A(node,edges) = 0;
            
            % add connections to current node
            edges = find(A(:,node));
            grp{end} = [grp{end}, setdiff(edges',grp{end})];
            
            % remove from adjacency matrix
            A(edges,node) = 0;
            
            nodeIdx = nodeIdx+1;
        end
        
        % add cell to path list if any edges remain
        if nnz(A) > 0
            grp = [grp; {[]}];
        end
    end
    
    % increment row to check in adjacency matrix
    rowNo = rowNo+1;
end

% sort list of nodes
grp = cellfun(@sort,grp,'uni',false);

% unsorted
solo = num2cell(setdiff(1:size(A,1),cell2mat(grp')))';