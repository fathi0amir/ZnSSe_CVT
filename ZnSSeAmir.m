%ZnSe 
%Calculating the ZnSSe parial pressures and optimal growth temperature
%clear
format long
c=[0.5 1 1.5 2 2.5 5 10 15];%The concentration of Iodine (Transport Agent)
x=0.5;%Determining the composition
%Here we calculate the pressure for different Iodine concentrations
for j=1:1:900%This loop takes care of Temperature vaiation
 t(j)=600+j;
   k2(j)=10^(7.64-5849./t(j)-4154/(t(j).^2)-0.83*log10(t(j))-1.5e-4*t(j));%Equilibrium Constant for ZnSe+I2->ZnI2+1/2Se2
   k3(j)=10.^(4.34-7879./t(j)+4264/(t(j).^2)+0.33*log10(t(j))+2.0e-5*t(j));%Equilibrium Constant for I2->2I
   k1(j)=10.^(8.8-7539./t(j)+8745/(t(j).^2)-1.19*log10(t(j))+2.18e-6*t(j));%Equilibrium Constant for ZnS+I2->ZnI2+1/2S2
   kmix(j)=k1(j).^x.*k2(j).^(1-x);%The total equlibrium constant
   for i=1:1:8%This loop takes care of different concentrations
            dyy=0.0;%A dummy variable to use in Newton Approximation
            dy=0.0001;%Another dummy variable to use in Newton Approximation
            for num=1:100%This loop mainly concerns with precision matter
               if abs(dy-dyy)<=1e-10%The accuracy limit
                break
                else 
                 dyy=dy;
                     f=2*2^(1/3)*x^(-x/3)*(1-x)^(-(1-x)/3)*(kmix(j)/k3(j))^(2/3)*dyy^(4/3)+2*dyy^2/k3(j)+dyy-6.4663042e-4*c(i)*t(j);%Difining a fuction for Iodine pressure
                     fp=4/3*2*2^(1/3)*x^(-x/3)*(1-x)^(-(1-x)/3)*(kmix(j)/k3(j))^(2/3)*dyy^(1/3)+4*dyy/k3(j)+1;%The Derivative of the above function
                     dy=dyy-f/fp;
               end
            end
%Now we calculate the pressures
       pi1(j,i)=dy;
       pi2(j,i)=(pi1(j,i)^2)/k3(j);
       pzni2(j,i)=((kmix(j)/k3(j))*(2/x)^(x/2)*(2/(1-x))^((1-x)/2))^(2/3)*pi1(j,i).^(4/3);
       pse2(j,i)=((1-x)/2)*pzni2(j,i);
       ps2(j,i)=(x/2)*pzni2(j,i);
       pi0(j,i)=2*pi2(j,i)+pi1(j,i)+2*pzni2(j,i);
       alpha(j,i)=pzni2(j,i)/pi0(j,i);%Defining Alpha
       pt(j,i)=pzni2(j,i)+pi2(j,i)+pi1(j,i)+pse2(j,i);
   end
end 
%And now we find the Delta Alpha as function of Delta T
for ii=1:1:8%Loop for different concentrations
            for jj=1:1:893%And this is for Temperature
              dalpha(jj,ii)=alpha(jj+7,ii)-alpha(jj,ii);
              dt(jj)=-273+(t(jj));
            end
 end