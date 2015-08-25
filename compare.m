function compare(k);
% s: sight 
% k: 0 turn off, 1 turn on
%clear all
if k==1, status = dos('./gipl'); end
m1=load('dump/result.txt');
m2=load('mesres.txt');


calc=m1(:,5:end);
mes=m2(:,2:end);

Depths=[0.001 0.072 0.125 0.2 0.277 0.354 0.424 0.506 0.583 0.741 0.885 1.1]; 
num=[3,5,8,11];  % indexes

n=length(calc(1:end,1));
%    1  2  3  4  5  6  7  8  9 10 11 12
%    J  F  M  A  M  J  J  A  S  O  N  D
mon=[31 28 31 30 31 30 31 31 30 31 30 31];
beg=6; k=1; 
for i=1+beg:12+beg
    if i==12, curr=i; else curr=rem(i,12);  end
    if k==1, month(k)=mon(curr); 
    else month(k)=month(k-1)+mon(curr); 
    end
    k=k+1; 
end

yy=fix(n/365);
for i=1:yy
for j=1:12
  k=j+(i-1)*12;
  xtick(k)=month(j)+365*(i-1);
end
xxtick(i)=(365-182)+365*(i-1);
end

nx=length(xtick);
xmonth={ 'A' 'S' 'O' 'N' 'D' 'J' 'F' 'M' 'A' 'M' 'J' 'J'};
xyear={ '2008' '2009' '2010' '2011' '2012'};
figure(1); 
for i=1:4
    subplot(4,1,i); 
  depth=num(i);
	plot(1:n,calc(:,depth),'b-','LineWidth',2); hold on % calculated
	plot(1:n,mes(1:n,depth),'r-','LineWidth',2); % mesaured
	err(i)=mean(abs(calc(:,depth)-mes(1:n,depth)));
	set(gca,'XTick',xtick); set(gca,'XTickLabel',xmonth)
	title(['Depth=',strcat(num2str(Depths(depth)),' m'), ',   MAE=',num2str(err(i))]); 
	ylabel 'Temperature (^oC)'; hold off; grid on;
	axis([1 n min(mes(1:n,depth))-3 max(mes(1:n,depth))+3])
end
xlabel 'Time (day)'
err(5)=sum(err(:))/4; disp([err]);
legend('Simulated','Observed',1);
