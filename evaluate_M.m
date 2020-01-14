function MR =  evaluate_M(du_List,NodeNum)
    tql = tabulate(du_List);
    Nr = tql(:,2);
    summer=0;
    for i =1:length(Nr)
        summer=summer+Nr(i)*(Nr(i)-1);
    end
    MR = (1-summer/(NodeNum*(NodeNum-1)))^2;
end