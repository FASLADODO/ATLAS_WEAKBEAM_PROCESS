function [ temp_filter_input2 ] = CIFilter(TBLHC2)
index_sig = abs(TBLHC2(:,end))>0;
temp_filter_input2 = TBLHC2;
filter_range =100;
foot_step = 50;
distance = max(temp_filter_input2(:,1));
seg_distance = 0:foot_step:distance+foot_step;
seg_distance = [seg_distance,max(seg_distance)+filter_range];
record = temp_filter_input2(:,1).*0;
for i = 1:length(seg_distance)
    index = (temp_filter_input2(:,1) >= seg_distance(i) - filter_range/2).*(temp_filter_input2(:,1) < seg_distance(i)+filter_range/2).*index_sig > 0;
    if ~isempty(index)
        seg_mean_evelation = mean(temp_filter_input2(index,2));
        seg_std_evelation = std(temp_filter_input2(index,2));
        record = record + ((abs(temp_filter_input2(:,2) - seg_mean_evelation) > 2*seg_std_evelation).*index);
        temp_filter_input2(record>0,end) = 0;
        index_sig = abs(temp_filter_input2(:,end))>0;        
    end
end
index_rec = record.*(abs(TBLHC2(:,end))>0)>0;
temp_filter_input2(index_rec,end) = 0;
index_final = find(abs(temp_filter_input2(:,end))>0);
temp_filter_input2 = BPC( index_final, temp_filter_input2);
end

