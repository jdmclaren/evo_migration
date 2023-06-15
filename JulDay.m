function jd=JulDay(date);
% julday function       convert date to Julian Day

% input  : - matrix in which the first column is the Day in the month
%            the second column is the month, the third is the year
%            and the fourth column is fraction of day (U.T.)
% output : - matrix of Julian Days.
%    By  Eran O. Ofek           January 1994
%--------------------------------------------------------------------
if nargin>1,
   error('1 arg only');
end
y = date(:,3);
m = date(:,2);
d = date(:,1);
f = date(:,4);
l = length(y);
for i=1:l,
   B = 0;
   if m(i)<3,
      y(i) = y(i) - 1;
      m(i) = m(i) + 12;
   end
   if ((y(i)>1582) | ((y(i)==1582) & (m(i)>10)) | ((y(i)==1582) & (m(i)==10) & (d(i)>=15)))   
      A = floor(y(i)./100);
      B = 2 - A + floor(A./4);
   end
   jd(i) = floor(365.25.*(y(i) + 4716)) + floor(30.6001.*(m(i) + 1)) + d(i) + B - 1524.5 + f(i);
end
