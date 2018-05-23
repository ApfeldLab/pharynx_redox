%  set up a mesh for smoothing weather data

%  load the lat/long data, lon first in degrees and minutes, and lat next

load LatLongsNum.txt

%  convert to decimal values and reverse latitude

Lat = LatLongsNum(:,3) + LatLongsNum(:,4)./60;
Lng = LatLongsNum(:,1) + LatLongsNum(:,2)./60;
Lat = max(Lat) - Lat;

%  load the daily temperature .mat file

load daily

%  display place names

disp([num2str(1:35)', place])

%  plot the points with place names

figure(1)
plot(Lat,Lng,'o')
for i=1:35
    text(Lat(i),Lng(i),place(i,:))
end

%  Delaunay mesh for all 35 points

Tri = delaunay(Lat,Lng);

triplot(Tri, Lat, Lng);

%  drop five stations that are too close to other stations

%  ptdrop Halifax, Baggotville, Montreal, Toronto, Vancouver

ptdrop = [2, 9, 12, 14, 26];
nptdrop = length(ptdrop);

%  construct the reduced set of points and their indices

LatLngptdrop = zeros(35-nptdrop,2);
LatLngptindx = zeros(35-nptdrop,1);
m = 0;
for i=1:35
    if ~any(ptdrop == i)
        m = m + 1;
        LatLngptdrop(m,1) = Lat(i);
        LatLngptdrop(m,2) = Lng(i);
        LatLngptindx(m)   = i;
    end
end

%  replot the mesh with numbered points

Tri = delaunay(LatLngptdrop);

figure(2)
triplot(Tri, LatLngptdrop(:,1), LatLngptdrop(:,2));
for i=LatLngptindx
    text(Lat(i),Lng(i),num2str(i));
end

%  now drop triangles that cover the exterior of Canada

tridrop = [6, 9, 10, 11, 22, 30, 38];

ntridrop = length(tridrop);

ntri = size(Tri,1);

Tridrop = zeros(ntri - ntridrop,3);
m = 0;
for i=1:ntri
    if ~any(tridrop == i)
        m = m + 1;
        Tridrop(m,:) = Tri(i,:);
    end
end

%  plot the final mesh with place names

figure(2)
triplot(Tridrop, LatLngptdrop(:,1), LatLngptdrop(:,2));
for i=LatLngptindx
    text(Lat(i),Lng(i),place(i,:));
end

%  set up the spatial smoothing algorithm

p = LatLngptdrop;
e = [];
t = Tridrop;

np = size(p,1);
nt = size(t,1);

%  add required paths

addpath('../..')

%  set up the FEM basis object and plot it

order = 1;
order = 2;

basisobj = create_FEM_basis(p, e, t, order);

%  set up a dummy FEM functional data object

precfd = fd(zeros(getnbasis(basisobj),1),basisobj);
tempfd = precfd;

%  set up the monthly precipitaton data

indrng = [  1, 31;
           32, 59;
           60, 90;
           91,120;
          121,151;
          152,181;
          182,212;
          213,243;
          244,273;
          274,304;
          305,334;
          335,365];
      
monthtemp = zeros(35,12);
monthprec = zeros(35,12);
for i=1:12
    monthtemp(:,i) = mean(tempav(indrng(i,1):indrng(i,2),:))';
    monthprec(:,i) = mean(precav(indrng(i,1):indrng(i,2),:))';
end
      
%  use log10 precipitation as data to be smoothed
%  if precipitation used, surface is dominated by Prince Rupert

precdata = log10(monthprec);
tempdata = monthtemp;

%  select the index sequence

% monthindex = 1:12;   %  calendar year
% monthindex = [7:12, 1:6];  % mid-summer to early-summer
monthindex = [5:12, 1:4];  % growing season to early spring

%  smooth the logdata and plot the surface for each month

precwrd = 1;

for imonth=monthindex
    
    if precwrd
        precdatai = [(1:np)', precdata(LatLngptindx,imonth)];
    else
        tempdatai = [(1:np)', tempdata(LatLngptindx,imonth)];
    end
    
    if order == 2
        lambda = 0;
        if precwrd
            [precfdi, laplacefdi] = smooth_FEM_fd(precdatai,precfd,lambda);
        else
            [tempfdi, laplacefdi] = smooth_FEM_fd(tempdatai,tempfd,lambda);
        end
    else
        lambda = 0;
        if precwrd
            precfdi = smooth_FEM_fd(precdatai,precfd,lambda);
        else
            tempfdi = smooth_FEM_fd(tempdatai,tempfd,lambda);
        end
    end
    
    %  plot the smoothed logdata
    
    figure(imonth)
    subplot(1,1,1)
    if precwrd
        plot(precfdi, 0, [], [], 51)
    else
        plot(tempfdi, 0, [], [], 51)
    end
    view(-10, 70)
    if precwrd
        title(['\fontsize{13} log10 precipitation for month ',num2str(imonth)])
    else
        title(['\fontsize{13} temperature for month ',num2str(imonth)])
    end
    colorbar
    
%     if order == 2
%         figure(imonth+12)
%         subplot(1,1,1)
%         plot(laplacefdi)
%         title(['\fontsize{13} laplacian of log10 precipitation for month ',num2str(imonth)])
%     end
%     colorbar
    
    pause
    
end


