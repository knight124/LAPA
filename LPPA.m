
clear
clc
close all

nx = 50;
ny = 50;
num_itteration = 10000;
line
dx =0.0001;% m
dy = dx;
V= zeros (nx,ny);
bc = [1,0,0,nan];% v (left top right bottom) NAN =free

Ex = zeros(nx,ny);
Ey = zeros(nx,ny);
% sets boundrt conditions 
V(1,:)=(~isnan(bc(2)))*bc(2);
V(end,:)=(~isnan(bc(4)))*bc(4);
V(:,1)=(~isnan(bc(1)))*bc(1);
V(:,end)=(~isnan(bc(3)))*bc(3);
V(isnan(V))=0;

imf = 1% use image filter  method 

if ~imf
    for l = 1:num_itteration
        %stores the previous run
        pv = V;
        % finds the finite difrence 
        for m = 1:ny
            for n=1:nx
                if (n == 1&&m ==1)
                    if (isnan(bc(1)))
                        V(m,n)=  (V(m+1,n)+V(m,n+1))/2;
                        
                    else
                        V(m,n)=bc(1);
                    end
                elseif (n == 1&&m ==ny)
                    if (isnan(bc(1)))
                        V(m,n)=   (V(m-1,n)+V(m,n+1))/2;
                    else
                        V(m,n)=bc(1);
                    end
                elseif (n == nx&&m ==1)
                    if (isnan(bc(3)))
                        V(m,n)=   (V(m+1,n)+V(m,n-1))/2;
                    else
                        V(m,n)=bc(3);
                    end
                elseif (n == nx&&m ==ny)
                    if (isnan(bc(3)))
                        V(m,n)=   (V(m-1,n)+V(m,n-1))/2;
                    else
                        V(m,n)=bc(3);
                    end
                elseif (n ==1)
                    if (isnan(bc(1)))
                        V(m,n)=   (V(m-1,n)+V(m+1,n)+V(m,n+1))/3;
                        
                    else
                        V(m,n)=bc(1);
                    end
                elseif (n ==nx)
                    if (isnan(bc(3)))
                        V(m,n)=   (V(m-1,n)+V(m+1,n)+V(m,n-1))/3;
                        
                    else
                        V(m,n)=bc(3);
                    end
                elseif (m ==1)
                    if (isnan(bc(2)))
                        V(m,n)=   (V(m,n-1)+V(m+1,n)+V(m,n+1))/3;
                        
                    else
                        V(m,n)=bc(2);
                    end
                elseif (m ==nx)
                    if (isnan(bc(4)))
                        V(m,n)=  (V(m,n-1)+V(m-1,n)+V(m,n+1))/3;
                        
                    else
                        V(m,n)=bc(4);
                    end
                else
                    V(m,n)=   (V(m+1,n)+V(m,n+1)+V(m-1,n)+V(m,n-1))/4;
                    
                end
                
            end
            
        end
        %calculates the efeild 
        for p=1:nx
            if p==1
                Ex(:,p) = -(-V(:,p)+V(:,p+1))/dx;
            elseif p==nx
                Ex(:,p) = -(-V(:,p-1)+V(:,p))/dx;
            else
                Ex(:,p) = -((-1/2)*V(:,p-1)+(1/2)*V(:,p+1))/dx;
            end
        end
        for p=1:ny
            if p==1
                Ey(p,:) = -(-V(p,:)+V(p+1,:))/dy;
            elseif p==nx
                Ey(p,:) = -(-V(p-1,:)+V(p,:))/dy;
            else
                Ey(p,:) = -((-1/2)*V(p-1,:)+(1/2)*V(p+1,:))/dy;
            end
        end
        %plots the results 
        if (mod(l,100)==1)
            xlabel('x position')
            ylabel('y position')
            zlabel('Voltage (v)')
            figure(2)
            quiver(1:1:nx,1:1:ny,Ex,Ey)
            xlabel('x position')
            ylabel('y position')
            pause
        end
        % exits if small enough diffrence
        if mean(mean(V-pv))<=1e-8
            break
        end
        
    end
    
else
    for l= 1:num_itteration
        pv = V;
        % finds the next finite diffrence itreation
        V = imboxfilt(V);
        temp = V;
        % resets boundry conditions 
        V(1,:)=(~isnan(bc(2)))*bc(2);
        V(end,:)=(~isnan(bc(4)))*bc(4);
        V(:,1)=(~isnan(bc(1)))*bc(1);
        V(:,end)=(~isnan(bc(3)))*bc(3);
        V(isnan(V))=temp(isnan(V));
        % finds the efeild
        for p=1:nx
            if p==1
                Ex(:,p) = -(-V(:,p)+V(:,p+1))/dx;
            elseif p==nx
                Ex(:,p) = -(-V(:,p-1)+V(:,p))/dx;
            else
                Ex(:,p) = -((-1/2)*V(:,p-1)+(1/2)*V(:,p+1))/dx;
            end
        end
        for p=1:ny
            if p==1
                Ey(p,:) = -(-V(p,:)+V(p+1,:))/dy;
            elseif p==nx
                Ey(p,:) = -(-V(p-1,:)+V(p,:))/dy;
            else
                Ey(p,:) = -((-1/2)*V(p-1,:)+(1/2)*V(p+1,:))/dy;
            end
        end
        % plots results 
        if (mod(l,100)==1)
            figure(1)
            surf(V)
            xlabel('x position')
            ylabel('y position')
            zlabel('Voltage (v)')
            figure(2)
            quiver(1:1:nx,1:1:ny,Ex,Ey)
            xlabel('x position')
            ylabel('y position')
            pause
        end
        % exits if small enough diffrence
        if mean(mean(V-pv))<=1e-8
            break
        end
    end
end
figure(1)
surf(V)
xlabel('x position')
ylabel('y position')
zlabel('Voltage (v)')
figure(2)
quiver(1:1:nx,1:1:ny,Ex,Ey)
xlabel('x position')
ylabel('y position')