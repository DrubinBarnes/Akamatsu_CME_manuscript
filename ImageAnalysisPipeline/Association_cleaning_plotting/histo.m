function h = histo(h,dat)

if (h.count==0)     
  h.numbins = ceil((h.range(2)-h.range(1))/h.binwidth);
  h.hist = zeros(1,h.numbins);
  h.vals = 1:h.numbins;
  h.vals = h.range(1)+(h.vals-0.5)*h.binwidth;
end

if (dat>h.range(1) && dat<h.range(2))
  bin = ceil((dat-h.range(1))/h.binwidth);
  h.hist(bin) = h.hist(bin)+1;
  h.count = h.count+1;
end