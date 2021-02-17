function [ temp_filter_input2 ] = BPC( index_final, temp_filter_input2)
i = 1;
while i < length(index_final)-1 
   for j = i+1 : length(index_final) 
        if temp_filter_input2(index_final(i),end)~=temp_filter_input2(index_final(j),end)
            slope_r = atand((temp_filter_input2(index_final(j),1)-temp_filter_input2(index_final(j-1),1))/(temp_filter_input2(index_final(j),2)-temp_filter_input2(index_final(j-1),2)));
            if slope_r < 0
                slope = - (90 - abs(slope_r));
            else
                slope = 90 - slope_r;
            end
            pos = [(temp_filter_input2(index_final(j-1),1)+temp_filter_input2(index_final(j),1))/2,(temp_filter_input2(index_final(j-1),2)+temp_filter_input2(index_final(j),2))/2];
            half_length = (sqrt(power(temp_filter_input2(index_final(j),1)-temp_filter_input2(index_final(j-1),1),2)+power(temp_filter_input2(index_final(j),2)-temp_filter_input2(index_final(j-1),2),2)))/2;
            A = half_length;
            if abs(slope) < 5
                B = A/6;
            else
                B = A/2;
            end
            theta = slope;
            TransX = temp_filter_input2(:,1) - pos(1);
            TransY = temp_filter_input2(:,2) - pos(2);
            TransX2 = TransX * cosd(theta) + TransY * sind(theta);
            TransY2 = -TransX * sind(theta) + TransY * cosd(theta);
            dis = power(TransX2,2)/power(A,2) + power(TransY2,2)/power(B,2);
            Neighbors=find(dis<1);
            temp_filter_input2(Neighbors,end)=-1;
            break;
        end
   end
   i = j;
end
end

