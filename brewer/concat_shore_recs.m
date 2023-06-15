recs = [];
for i = 9:13 % 1:3 % 
    
    recs_i = load(['shore_' num2str(i) '.mat']);
    recs = [recs recs_i];
    
end

save('shore_9_13','recs','-v7.3')