%% PLOTMARKERS - makes plotting markers
% This simple script generates a large set of ordered plot markers

function M = plotmarkers()

M = struct('line',{[]},'mark',{[]});

markers = {[],'.','*','o','+','x','s','d','^','v','>','<','p','h'};
lineStyles = {'-','--',':','-.'};
colors = {'k','b','r','g','c','m'};

M.line = cell(length(lineStyles) * length(colors) * length(markers),1);

idx = 1;
for i = 1:length(markers)
    for j = 1:length(lineStyles)
        for k = 1:length(colors)
            M.line{idx} = [colors{k}, lineStyles{j}, markers{i}];
            idx = idx + 1;
        end
    end
end

markers = markers(2:end);

M.mark = reshape(repmat(markers, length(colors), 1),...
    length(colors)*length(markers), 1);

M.mark = cellfun(@(x,y) strcat(x,y), repmat(colors', length(markers), 1),...
    M.mark, 'uniformoutput', false);