function [Cleaned_Regions]=Visual_Data_Partition(dqF)
    %Input1: Data Partition in Raw in vector -->R nx1
    
    %Output1: Renumbered data partitions

    %Originally Written by Soner Günaydin, January 2022
    %Simplified and changed by Yuri Shardt, March 2022

    %Define the parameter
    Insufficient_Data=[];%Plotting the Insufficient Data and Output1
    Segmentation_Filter=[];%Clean the data and get the absolut data
    color=[];%for generating randomly colors for ALL partition data
    Partition_Change=[];%count the changes
    model_counter=[];%Counts the partition data
    New_Part_Nr=[];%show all Partition_Number
    %%
    %Default setting if nothing is selected in the input
     
        Parameter_Input=0;
    
        Trash=1;
   
    %%
    %this loop split the insufficient Data and the Segmented Data and
    %converts
    %negative Values into absolute Values
    a=1;
    i=1;
    while i<length(dqF)     
        if dqF(i) ~= dqF(i+1)
            Insufficient_Data(a,:)=.1;
        else
            Segmentation_Filter(a,:)=abs(dqF(i));
            Insufficient_Data(a,:)=0;
        end
            i=i+1;
            a=a+1;
    end
    %%
    %%Plotting Settings start
    
    figure('Name','Partition_Data');%Name the Figure
    grid on;%Plotting with Grid
    hold on;%that all plots contain in one graph
    %Link='.-';%Define the Connection Between the values
    
    %%Plotting settings End
    %%
    %This loop finds all changes to a variable, so the 
    %color assignment of the partition number should follow;
    
    i=1;%return the value for the next loop
    k=1;%return the value for the next loop

    for i=1:length(Segmentation_Filter)
        if i==length(Segmentation_Filter)
                Partition_Change(k,:)=Segmentation_Filter(i,:);
                break;
        end         
        if Segmentation_Filter(i,:) ~= Segmentation_Filter(i+1,:) && Segmentation_Filter(i,:)~=0
           Partition_Change(k,:)=Segmentation_Filter(i,:);
           k=k+1;
        end
    end
    %%
    %this loop generate randomly colors with the size in Partition_Change and
    %take a sum how much Partition Number exists
    
    i=1;%return the value for the next loop
    while i<=length(Partition_Change)
        color(i,:)=[rand(1) rand(1) rand(1)];
        model_counter(i,:)=i;
        i=i+1;
    end
    %%
    %this loop assigns the Partition_Number with its color and new
    %Partition Number 
    
    Col=0;%default Value for the Partition and Colors
    i=1;%return the value for the next loop
    j=1;%return the value for the next loop
    counter=1;
       
%     while i<=length(Partition_Change)   
%         j=1;
%         Col=Partition_Change(i,:);
%         while j<=length(Partition_Change) 
%             if Col==Partition_Change(j,:) 
%                color(i,:)=color(j,:);
%                New_Part_Nr(i,:)=counter;%model_counter(j,:);
%                counter=counter+1;
%             end
%             j=j+1;
%         end
%         i=i+1;
%     end
    change=0;
    for i=1:length(Segmentation_Filter)-1
        if Segmentation_Filter(i+1)-Segmentation_Filter(i)>0 && change==0
            New_Part_Nr(i,:)=0;
            color(i,:)=color(end,:);
            change=1;
        end
        if Segmentation_Filter(i+1)-Segmentation_Filter(i)==0 && change==1
            New_Part_Nr(i,:)=counter;
            color(i,:)=color(counter,:);
        end
        if Segmentation_Filter(i+1)-Segmentation_Filter(i)==0 && change==0
            New_Part_Nr(i,:)=0;
            color(i,:)=color(end,:);
        end
        if Segmentation_Filter(i+1)-Segmentation_Filter(i)<0 && change==1
            New_Part_Nr(i,:)=counter;
            color(i,:)=color(counter,:);
            counter=counter+1;
            change=0;
        end
    end
    i=length(Segmentation_Filter);
        if Segmentation_Filter(i-1)-Segmentation_Filter(i)==0 && change==1
            New_Part_Nr(i,:)=counter;
            color(i,:)=color(counter,:);
        elseif Segmentation_Filter(i-1)-Segmentation_Filter(i)==0 && change==0
            New_Part_Nr(i,:)=0;
            color(i,:)=color(end,:);
        elseif Segmentation_Filter(i)-Segmentation_Filter(i-1)<0 && change==1
            New_Part_Nr(i,:)=counter;
            color(i,:)=color(counter,:);
            counter=counter+1;
            change=0;
        else
            New_Part_Nr(i,:)=counter;
            color(i,:)=color(counter,:);
        end
        New_Part_Nr=[zeros(1,size(New_Part_Nr,2));New_Part_Nr];
    
    
    %%
%     %this Loop Plots the Regions with their partition number and the color
%     i=1;%return the value for the next loop
%     j=1;%return the value for the next loop
%     a=1;%moving the color to the next number
%     k=1;%Definded Variable for the filtered Regions
%     
%     while i<=length(Segmentation_Filter)       
%         i=j;%Jump to the next Partitionnumber       
%         while Segmentation_Filter(j,:)==Segmentation_Filter(j+1,:) 
%            Partition_Data(k,:)=Segmentation_Filter(j,:);%Filtered Region
%            j=j+1;  
%            k=k+1;
%            if j==size(Segmentation_Filter,1)%secure, that the loop are defined on the last Number
%                break%Break the loop if the size of Segmentation_Filter is reached
%            end
%         end        
%         if Segmentation_Filter(j,:)~=0 && length(i:j)>=Parameter_Input
%             New_Region(i:j,:)=New_Part_Nr(a,:);
%             Temp_Region=New_Region(i:j,:);%Set the defined region
%             plot((i:j),Temp_Region,Link,'Color',color(a,:)) %plot the area
%             Temp_Region=[];%reset the Region 
%             a=a+1;%moving the color forward
%         end
%         j=j+1;       
%         if  j>=size(Segmentation_Filter,1)%secure, that the loop are defined on the last Number
%             break%Break the loop if the size of Segmentation_Filter is reached
%         end
%     end
%     %%
%     %this loop generates the legend for the Partition Number
%     i=1;%return the value for the next loop
%     j=1;%return the value for the next loop
%     
%     while j<=(a)-1
%         leg(j)={['Partition # ',num2str(j)]};
%         j=j+1;
%     end
%     
%     if Trash==1
%         leg(j)={'Insufficient Data'};%Legend for the Insufficient Data
%     end
%     %%
%     %this loop plots the Insufficent Data
%     while i<=length(Insufficient_Data)
%         if  Trash~=1
%             break;
%         end
%         if Insufficient_Data(i)~=0 
%             plot(i,Insufficient_Data(i)*0,'*-r');
%         end
%             i=i+1;
%     end
%     %%
%     %Plotting Settings begin
%     %This Settings define the legends 
%     lgd=legend(leg);
%     lgd.Box="off";
%     lgd.NumColumns = 100;
%     lgd.AutoUpdate='off';
%     lgd.Orientation = 'vertical';
%     lgd.Location='southoutside';
    plot(New_Part_Nr,'o');
    %This settings define the title xlabel and ylabel
    title('Data Partition');
    xlabel('Data Point');
    ylabel("Partition Number");
    %Plotting Settings end
    hold off;
    Cleaned_Regions=New_Part_Nr;
    %%
    %for the Output1, filtered Data without insufficient Data and small
    %Partition Number

%     a=1;%return the value for the next loop
%     j=1;%return the value for the next loop
%     for i=1:length(New_Region)
%          
%         if New_Region(i,:)~=0 && length(i:j)<=Parameter_Input
%             Regions_Without_Ins(a,:)=New_Region(i,:);
%             a=a+1;
%         else
%             j=i;
%         end      
%     end
end