function result = CalPhaseScore(x0, intlist)
% x0: [int12 int13]
    intlist(:,2) = intlist(:,2).*x0(1);
    intlist(:,3) = intlist(:,3).*x0(2);
    temp = CalPhaseListLite(intlist);
    md = temp(:,4);
    md=md(md>0.3 & md<2);
    result = std(md);
end