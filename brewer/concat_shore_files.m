% fnames = {'shore_1_3','shore_4_8','shore_9_13'};

n_files = 13; % numel(fnames);

for ifile = 8:n_files % 2 % 
    
    tic
    display(['loading file ' num2str(ifile)])
%     shore_i = load(['F:\Sensorgnome\Sjoerd_data_2017\' fnames{ifile}]);
    shore_i = load(['E:\Sensorgnome\Sjoerd_data_2017\shore_' ...
        num2str(ifile)]);
    toc
    
%     names = fieldnames(shore_i.recs);
    names = fieldnames(shore_i);
    n_vars = numel(names);
    
    if ifile == 8
        shore_all = cell(1,n_vars);
    end

    gg = 0;
%     for irow = 1:numel(shore_i.recs)

        for jj = 3:n_vars
%            shore_all{jj} = [shore_all{jj}' eval(['shore_i.recs(' ...
%                num2str(irow) ').' names{jj}])']';
            var_i = eval(['shore_i.' names{jj}]);
            % remove pesky NA values
            if strcmp(class(var_i),'cell')
                var_i(strcmp(var_i,'NA')) = 'NaN';
            end
           shore_all{jj} = horzcat(shore_all{jj}',var_i')';
        end

%     end

end

save('shore_2013_2016','shore_all')
cols = 1:83; % [3,4,5,7,8,36,37,39,40,64,72];
for ii = 1:numel(cols)
display([num2str(cols(ii)) ' ' names{cols(ii)} ' ' ...
    num2str(size(shore_all{cols(ii)},1)) ' ' class(shore_all{cols(ii)})])
end