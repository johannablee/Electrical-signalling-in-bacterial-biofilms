% Polar coordinate data analysis of biofilms by Johanna Blee
%8/12/2017
%This code takes a video containing a time lapse of a biofilm with a 
%radially propagating signal and converts the signal into polar coordinates
%and radial averages. The moments of the radial average signal are found, 
%and used to calculate mean, std, kurtosis and skewness. 
%The average radial signal is also fitted with Gaussians. (Requires wmean.m)


%Create VideoReader
biofilmVideo = VideoReader('timelapsevideo.avi'); %Insert video file name here 
ii = 1;
px=70/220; %pixels calibration


%read in video of timelapse
while hasFrame(biofilmVideo)
   img = readFrame(biofilmVideo);
   images{ii}=img;  %ii=number of frames in video
   red(:,:,ii)=img(:,:,1); %creates an array of red parts of image
   green(:,:,ii)=img(:,:,2); %creates an array of green parts of image
   blue(:,:,ii)=img(:,:,3); %creates an array of blue parts of image 
   gray{ii}=0.2989*red(:,:,ii)+0.5870*green(:,:,ii)+0.1140*blue(:,:,ii); %creates a grayscale image
   Gray(:,:,ii)=0.2989*red(:,:,ii)+0.5870*green(:,:,ii)+0.1140*blue(:,:,ii);
   Avin1(:,ii)= mean(Gray(:,:,ii));%creates vector for the average intensity at one distance through the biofilm
   Avin2(ii)= mean(Avin1(:,ii));  %creates a mean intensity over all pixels in time
   ii = ii+1;
end

 %creates matrices of average intensity at each rho for each image
[ab ba ca]=size(Gray(:,:,:));


%converting to polar coordinates 
for jk=1:ca 
    for k=1:ab
        for p=1:ba
            %thresholding to remove cell trap and background-just get left with cells
            if max(Gray(k,p,:))<7
                rem(k,p,jk)=NaN;
            else
                rem(k,p,jk)=Gray(k,p,jk);
            end 
        end 
    end     
    M=double(rem(:,:,jk));
    X0=334; Y0=239; %central biofilm coordinates
    [Y X z]=find(M); 
    X=X-X0; Y=Y-Y0;
    theta = atan2(Y,X); %creates theta values
    rho = sqrt(X.^2+Y.^2); %creates rho values

    % Determine the minimum and the maximum x and y values:
    rmin = min(rho); tmin = min(theta);
    rmax = max(rho); tmax = max(theta);

    % Define the resolution of the grid:
    rres=300; % # of grid points for R coordinate. (change to needed binning)
    tres=300; % # of grid points for theta coordinate (change to needed binning)

    F = TriScatteredInterp(rho,theta,z,'natural');

    %Evaluate the interpolant at the locations (rhoi, thetai).
    %The corresponding value at these locations is Zinterp:

    [rhoi,thetai] = meshgrid(linspace(rmin,rmax,rres),linspace(tmin,tmax,tres));
    Zinterp = F(rhoi,thetai); %this is matrices of intesnity values of image in rho and theta instead of x and y
    values(:,:,jk)=Zinterp;
    [ik ij]=size(Zinterp);

    for iji=1:ij
        J=Zinterp(:,iji);
        Rhoav(iji)=nanmean(J); %creates an average over all theta for each value of rho.
    end 
    intrho(:,jk)=Rhoav;

end 
 
[aa, bb]=size(intrho);
T=10*linspace(1,bb,bb);%creates vector for time. Change according to time step size of timelapse video
X=px*((sqrt((X0*X0)+(Y0*Y0)))/aa)*linspace(1,aa,aa); %creates vector for the distance of each point
%Finding the signals mean, std, kurtosis,skewness for each theta and rho
%and then finding averages
for iji=1:ij
    for gg=1:ik
        xxy1(:)=values(gg,iji,:);
        radialmean1(gg,iji)=wmean(T,xxy1);%signals mean at each theta and rho
        radialstd1(gg,iji)=sqrt(nanvar(T,xxy1));%signal std at each theta and rho
        ww1=nansum(xxy1);
        xx1=sum(T);
        wwxx1=(xxy1.*T);
        radialskewness1(gg,iji)=((sum(xxy1.*((T-radialmean1(gg,iji))/radialstd1(gg,iji)).^3))/(ww1));%signal skewness at each theta and rho
        radialk1(gg,iji)=((sum(xxy1.*((T-radialmean1(gg,iji))/radialstd1(gg,iji)).^4))/(ww1));%signal kurtosis at each theta and rho
        radialA1(gg,iji)=trapz(T,xxy1);%signal area at each theta and rho
        errskewness1(gg,iji)=sqrt(6*ww1*(ww1-1)/((ww1-2)*(ww1+1)*(ww1+3))); %standard error on the skewness
        errkurtosis1(gg,iji)=4*(ww1^2-1)*errskewness1(gg,iji) / ((ww1-3)*(ww1+5));%standard error on the kurtosis
        for tt=1:length(T)-1
            allerrorAsqrd1(tt)=(((T(tt+1)-T(tt))^2/12)*((xxy1(tt+1)-xxy1(tt))/(T(tt+1)-T(tt))))^2;
        end
        errA11(gg,iji)=sqrt(sum(allerrorAsqrd1));%standard error on the signal area
        
    
    end
    avskewness(iji)=nanmean(radialskewness1(:,iji));% average radial skewness
    stdskewness(iji)=nanstd(radialskewness1(:,iji));% standard deviation on radial skewness
    avkurtosis(iji)=nanmean(radialk1(:,iji));% average radial kurtosis
    stdkurtosis(iji)=nanstd(radialk1(:,iji));% standard deviation on radial kurtosis
    avA(iji) = nanmean(radialA1(:,iji));%average radial signal area
    stdA(iji) = nanstd(radialA1(:,iji));%standard deviation on radial signal area
    avstd(iji)=nanmean(radialstd1(:,iji)); % average standard deviation of radial average
    stdstd(iji)=nanstd(radialstd1(:,iji));% standard deviation of the signals standard deviation of radial average
    avmean(iji)=nanmean(radialmean1(:,iji));% Radial average of signals mean time
    stdmean(iji)=nanstd(radialmean1(:,iji));%standard deviation of radial average of signals mean time of signals mean time
    covmean(iji)=stdmean(iji)/avmean(iji);%coefficient of variation of the signals radial mean 
end 



intrho=intrho-(min(min(intrho))); % background subtraction
intrho=intrho/max(max(intrho));%signal normalisation
[aa bb]=size(intrho);
startcell=15/(px*((sqrt((X0*X0)+(Y0*Y0)))/aa)); % start of flow trap in pixels
endcell=30/(px*((sqrt((X0*X0)+(Y0*Y0)))/aa)); % end of flow trap in pixels
biofilmedge=84/(px*((sqrt((X0*X0)+(Y0*Y0)))/aa)); % biofilm edge
for jj=1:aa  %this is to remove flow trap to insure not affecting results and biofilm so only get biofilm cells
    if jj > startcell && jj < endcell
        intrho1(jj,:)= NaN(1,bb);
        avmean1(jj)=NaN;
        stdmean1(jj)=NaN;
        covmean1(jj)=NaN;
    elseif jj > biofilmedge
        intrho1(jj,:)= NaN(1,bb);
        avmean1(jj)=NaN;
        stdmean1(jj)=NaN;
        covmean1(jj)=NaN;
    else
        intrho1(jj,:)=intrho(jj,:);
        avmean1(jj)=avmean(jj);
        stdmean1(jj)=avstd(jj);
        covmean1(jj)=covmean(jj);
    end 
end 



%fits gaussians to signals of all radial averages and also calculates
%signal mean, std, kurtsosis, skewness and signal area of all radial
%avergaes.
for lll=1:aa
    xxy=intrho1(lll,:);
    try
       %attempt gaussian fitting also possible to use nlinfit if
       %curvefittingtoolbox not available
       f = fit(T.',xxy.','gauss1'); %f(x) =  a1*exp(-((x-b1)/c1)^2)
       g=coeffvalues(f); % gives coeffiencients of fits
       h=confint(f);%give confidence intervals of fits
       a1(lll)=g(1);%amplitude
       b1(lll)=g(2);%signal mean time
       c1(lll)=g(3);%signal width
    catch
       %if fails to fit gaussians set fit parameters to NaN
       f=NaN;
       h=NaN;
       a1(lll)=NaN;
       b1(lll)=NaN;
       c1(lll)=NaN;
    end
    
    %same as before to find mean,std,skewness,kurtosis and signal
    %area
    radialmean(lll)=wmean(T,xxy);
    try
       radialstd(lll)=sqrt(var(T,xxy));
    catch
       radialstd(lll)=NaN;
    end
    ww=sum(xxy);
    xx=sum(T);
    wwxx=(xxy.*T);
    radialskewness(lll)=((sum(xxy.*((T-radialmean(lll))/radialstd(lll)).^3))/(ww));
    radialk(lll)=((sum(xxy.*((T-radialmean(lll))/radialstd(lll)).^4))/(ww));
    radialA(lll)=trapz(T,xxy);
    errskewness(lll)=sqrt(6*ww*(ww-1)/((ww-2)*(ww+1)*(ww+3)));
    errkurtosis(lll)=4*(ww^2-1)*errskewness(lll) / ((ww-3)*(ww+5));
    for tt=1:length(T)-1
        allerrorAsqrd(tt)=(((T(tt+1)-T(tt))^2/12)*((xxy(tt+1)-xxy(tt))/(T(tt+1)-T(tt))))^2;
    end
    errA1(lll)=sqrt(sum(allerrorAsqrd));
end

%contour plot of signal
figure;
set(gcf,'color','w');
contourf(T,X,intrho1);
xlabel('Time (min)');
ylabel('Radial Distance (\mum)');
colorbar;
ylim([0 85]);
ylabel(colorbar,'Tht (a.u)');
xlim([110 550]);



%quick plots to check variables but for figures and fitting I export
%variables to origin
figure
set(gcf,'color','w');
plot(radialmean,X)

figure
set(gcf,'color','w');
plot(X,radialA);





