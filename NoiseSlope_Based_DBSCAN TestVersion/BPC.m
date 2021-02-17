function [ temp_filter_input2 ] = BPC( index_final, temp_filter_input2)
i = 1;
%% Test
figure()
index_final2 = find(temp_filter_input2(:,end)>0);
scatter(temp_filter_input2(:,1),temp_filter_input2(:,2),'.');
hold on
scatter(temp_filter_input2(index_final2,1),temp_filter_input2(index_final2,2),'.');

while i < length(index_final)-1 
   for j = i+1 : length(index_final) 
        if temp_filter_input2(index_final(i),end)~=temp_filter_input2(index_final(j),end)
            slope_r = atand((temp_filter_input2(index_final(j),1)-temp_filter_input2(index_final(j-1),1))/(temp_filter_input2(index_final(j),2)-temp_filter_input2(index_final(j-1),2)));
            if slope_r < 0
                slope = - (90 - abs(slope_r));
            else
                slope = 90 - slope_r;
            end
            x1 = temp_filter_input2(index_final(j-1),1);
            y1 = temp_filter_input2(index_final(j-1),2);
            x2 = temp_filter_input2(index_final(j),1);
            y2 = temp_filter_input2(index_final(j),2);
            slope_e = atand((y2-y1)/(x2-x1));
%             [ellipse] = DrawEllipse( x1,y1,x2,y2,slope_e);
%             plot(ellipse(:,1),ellipse(:,2),'r');
%             plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'r-');
            pos = [(x1+x2)/2,(y1+y2)/2];
            half_length = sqrt((x2-x1)^2+(y2-y1)^2)/2;
%             if half_length<=50
            A = half_length;
            B=A/4;
            [ellipse] = DrawEllipse( x1,y1,x2,y2,slope_e,A,B);
            plot(ellipse(:,1),ellipse(:,2),'r');
%                 if abs(slope) < 5
%                     B = A/6;
%                 else
%                     B = A/2;
%                 end
                theta = slope;                
                TransX = temp_filter_input2(:,1) - pos(1);
                TransY = temp_filter_input2(:,2) - pos(2);
                TransX2 = TransX * cosd(theta) + TransY * sind(theta);
                TransY2 = -TransX * sind(theta) + TransY * cosd(theta);
%                 S1 = abs(TransX2)<=A;
%                 S2 = abs(TransY2)<=B;
%                 Neighbors=S1.*S2>0;
                dis = power(TransX2,2)/power(A,2) + power(TransY2,2)/power(B,2);
                Neighbors = dis<1;
%                 scatter(temp_filter_input2(Neighbors,1),temp_filter_input2(Neighbors,2),'g.');              
                if B < 5
                    temp_filter_input2(Neighbors,end)=-1;
                else
                    temp_filter_input2= SSDF( temp_filter_input2,Neighbors,slope_e,slope);
                end
                break;
%             end
        end
   end
   i = j;
end
end

