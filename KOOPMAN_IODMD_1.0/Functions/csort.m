function [A,I]=csort(AA,i)

[AAA,I]=sort(AA(:,i));

A=AA(I,:);

