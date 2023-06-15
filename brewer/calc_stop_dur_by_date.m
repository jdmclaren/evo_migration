           % first determine stop durn previous
           first_beg = find(stop_nr_beg ...
               (ii,k_y,:) == ...
               current_st_nr,1,'first');
           last_beg = find(stop_nr_beg ...
               (ii,k_y,:) == ...
               current_st_nr,1,'last');
           
           
           % if previous rec was there already
           % (previously) a new stop on date with earlier stop recs?
           % (need to check e.g. for new flyby within day)
           if ~isempty(first_beg) 

               prev_rec = first_beg < jd && first_beg > 1;
               if prev_rec == 1
                   last_new = find(stop_nr_new ...
                   (ii,k_y,:) == ...
                   current_st_nr,1,'last');
                  prev_new = ~isempty(last_new);
%                     prev_new = stop_nr_new(ii,k_y, ...   
%                     first_beg -1) == stop_nr_beg(ii,k_y,jd);
                    if prev_new == 1
                        t_first = lst_t_day_new ...
                        (ii,k_y,last_new);
                    else
                        t_first = fst_t_day_stop ...
                        (ii,k_y,first_beg);
                    end
               else
                   prev_rec = 0;
                    t_first = fst_t_day_stop ...
                    (ii,k_y,first_beg);
               end
               
               % now calc stop d (know there was pev rec thus also last
               % one)
              stop_dur = lst_t_day_stop(ii,k_y, ...
                  last_beg) - t_first;
              stop_durs(ii,k_y,current_st_nr) = stop_dur;
              % calculate detct through stopover
              dtcts = squeeze(dtct_stop_beg(ii,k_y, ...
                  first_beg:last_beg));
              propn_dtct(ii,k_y,current_st_nr) = ...
                  (nansum(dtcts)+prev_rec)/ ...
                  ceil(last_beg - first_beg +1 +prev_rec);
              % find max gaps (use matlab pixel trick);
              dtct_str= bwconncomp(dtcts==0);
              str_size = ...
                  cellfun('prodofsize',dtct_str.PixelIdxList);
              non_dtcts = zeros(size(dtcts));
              for i_ob = 1:dtct_str.NumObjects
                non_dtcts(dtct_str.PixelIdxList{i_ob}) = ...
                    str_size(i_ob);
              end
              max_non_dtct(ii,k_y,current_st_nr) = ...
                  max(non_dtcts);

           else % No prev begin-date recs for this stopover 
               % Could be one of two things
                % 1) all recs on same date last stop ended
                % 2) all recs on same day after season changed
                % i.e. current_st_nr is updated but prev stopover durn not
                % yet processed

                % if first date only former
                % If not check season of 
                % previous stopover event 

                if jd > 1 

                    % check I didn't mess up change of seasons
                    % -stopover should be processed before current_st_nr
                    % is updated
                    if current_st_nr > 1
                        
                        if j_seasn(current_st_nr) ~= ...
                          j_seasn(current_st_nr-1)
                            
                          % should not come here! Not keeping track of 
                          % last record for date in earlier season 
                          % (could be many missing in between)
                          display(['error: change of season in' ...
                          calc_stop_dur subr.'])
                          keyboard
                      
                        end
%                       % change season, find so only recs on current date jd
%                       stop_durs(ii,k_y,current_st_nr-1) = ...
%                        lst_t_day_stop(ii,k_y,jd) - ...
%                        fst_t_day_stop(ii,k_y,jd);
%                         % in all cases with only one day of detections,
%                         % there are no dates missing detections
%                        propn_dtct(ii,k_y,current_st_nr-1) = 1;
%                        max_non_dtct(ii,k_y,current_st_nr-1) = 0;
                      
                    end
                      
                         last_new = find(stop_nr_new ...
                       (ii,k_y,:) == current_st_nr,1,'last');
                       stop_durs(ii,k_y,current_st_nr) = ...
                       lst_t_day_new(ii,k_y,last_new) - ...
                       fst_t_day_new(ii,k_y,last_new);
                   
                       % since only one day of detections,
                        % there are no dates missing detections
                       propn_dtct(ii,k_y,current_st_nr) = 1;
                       max_non_dtct(ii,k_y,current_st_nr) = 0;
                       
%   dt                  end
                    
                end

           end
           
           if  propn_dtct(ii,k_y,current_st_nr) == 1 & ...
                  stop_durs(ii,k_y,current_st_nr) > 70
              
              keyboard
              
           end