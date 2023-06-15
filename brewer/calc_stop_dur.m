                       % first determine stop durn previous
                       first_beg = find(stop_nr_beg ...
                           (ii,js,k_y,:) == ...
                           current_st_nr,1,'first');
                       last_beg = find(stop_nr_beg ...
                           (ii,js,k_y,:) == ...
                           current_st_nr,1,'last');
                       % if previous rec was there already
                       % new stop?
                       %(need to check for new flyby within day)
%                        if ~isempty(first_beg) 
%                             prev_rec = first_beg < jd && first_beg > 1;
%                        else
%                             prev_rec = 0;
%                        end
                       
                       if ~isempty(first_beg) 
                           
                           prev_rec = first_beg < jd && first_beg > 1;
                           if prev_rec == 1
                                prev_new = stop_nr_new(ii,js,k_y, ...   
                                first_beg -1) == stop_nr_beg(ii,js,k_y,jd);
                                t_first = lst_t_day_new ...
                                (ii,js,k_y,first_beg -1);
                           else
                               prev_rec = 0;
                                t_first = fst_t_day_stop ...
                                (ii,js,k_y,first_beg);
                           end
                           
                          sd = lst_t_day_stop(ii,js,k_y, ...
                              last_beg) - t_first;
                          stop_durs(ii,js,k_y,current_st_nr) = sd;
                          dtcts = squeeze(dtct_stop_beg(ii,js,k_y, ...
                              first_beg:last_beg));
                          propn_dtct(ii,js,k_y,current_st_nr) = ...
                              nansum(dtcts+prev_rec)/ ...
                              ceil(last_beg - first_beg +1 +prev_rec);
                          % find max gaps
%                          if sum(dtcts == 0)>0
%                               keyboard
%                           end
%                            dtcts(isnan(dtcts)) = 0;
                          dtct_str= bwconncomp(dtcts==0);
                          str_size = ...
                              cellfun('prodofsize',dtct_str.PixelIdxList);
                          non_dtcts = zeros(size(dtcts));
                          for i_ob = 1:dtct_str.NumObjects
                            non_dtcts(dtct_str.PixelIdxList{i_ob}) = ...
                                str_size(i_ob);
                          end
                          max_non_dtct(ii,js,k_y,current_st_nr) = ...
                              max(non_dtcts);
                     
                       else % only same day last stop ended
                           
%                            first_new = find(stop_nr_new ...
%                            (ii,js,k_y,:) == ...
%                            current_st_nr,1,'first');
%                             prev_rec = 0;
                            
                             last_new = find(stop_nr_new ...
                           (ii,js,k_y,:) == current_st_nr,1,'last');
                           stop_durs(ii,js,k_y,current_st_nr) = ...
                           lst_t_day_new(ii,js,k_y,last_new) - ...
                           fst_t_day_new(ii,js,k_y,last_new);
                              propn_dtct(ii,js,k_y,current_st_nr) = 1;
                               max_non_dtct(ii,js,k_y,current_st_nr) = 0;
                              
                       end