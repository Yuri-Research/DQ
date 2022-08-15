function [Cleaned_Regions]=simplifyDataPartitions(dqF)
    %Input1: Raw Data Partition as a vector -->R^{nx1}
    
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
    %This loop splits the insufficient data and the segmented Data and
    %converts negative values into absolute values
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
    change=0;
    counter=1;
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
    Cleaned_Regions=New_Part_Nr;
 end