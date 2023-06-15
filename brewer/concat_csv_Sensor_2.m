% fnames = {'shore_1_3','shore_4_8','shore_9_13'};

n_files = 2; % 13; % numel(fnames);

for ifile = 1:n_files % 2 % 
    
    tic
    display(['loading file ' num2str(ifile)])
%     shore_i = load(['F:\Sensorgnome\Sjoerd_data_2017\' fnames{ifile}]);
    shore_i = load(['E:\Sensorgnome\Sjoerd_data_2017\shore_' ...
        num2str(ifile)]);
    toc
    
%     names = fieldnames(shore_i.recs);
    names = fieldnames(shore_i);
    n_vars = numel(names);
    n_recs = 0;
    n_rows = size(eval(['shore_i.' names{1}]),1);
    
    if ifile == 1

        shore_all = struct(shore_i); % cell(1,n_vars);
        n_recs = n_rows;
        
    else
         rows = (n_recs+1):(n_recs+n_rows);
        for ii=1:numel(names)
              d=[];
           
              eval(['d=shore_i.' names{ii} ';'])
%               end
              eval([sprintf('shore_all.%s',names{ii}) '=d;'])
        end
      n_recs = rows(end);
    end


end

save('shore_2013_2016','shore_all')
cols = 1:83; % [3,4,5,7,8,36,37,39,40,64,72];
% for ii = 1:numel(cols)
% display([num2str(cols(ii)) ' ' names{cols(ii)} ' ' ...
%     num2str(size(shore_all{cols(ii)},1)) ' ' class(shore_all{cols(ii)})])
% end