clc
clear;
%Define the input parameter
Re = 100; 
dt = 0.001;
nx = 7; 
ny = 7; 
length_x = 6;
length_y = 6;
hx = length_x/(nx-1);
hy = length_y/(ny-1);
x = linspace(0,length_x*(1+1/(nx-1)),nx+1);
y = linspace(0,length_y*(1+1/(ny-1)),ny+1);

[Xx,Xy]=meshgrid(x(1:end-1),y);
[Yx,Yy]=meshgrid(y,x(1:end-1));
[Px,Py] = meshgrid(x,y);

%%
%Initialize the variables
%Final collocated mesh variable
u_final(nx,ny) = 0; 
v_final(nx,ny) = 0;
p_final(nx,ny) = 1;
u_final(1,:)=1;
%Staggered variables
%Staggered variables used for iteration
u(nx+1,ny) = 0;
v(nx,ny+1) = 0;
p(nx+1,ny+1) = 1;
u(1,:) = 2;
%Staggered variables used for update value
u_new(nx+1,ny) = 0;
v_new(nx,ny+1) = 0;
p_new(nx+1,ny+1) = 1;
u_new(1,:) = 2;
%%
%Solve the governing equation
error = 1;
iteration = 0;
error_req = 1e-7; %final error require for stable solution
figure(1) %for error monitoring

while error > error_req
%x-momentum eq. -interior
   for i = 2 : nx
      for j = 2 : ny-1
           advection_x = ((0.5*(u(i,j)+u(i,j+1)))^2-(0.5*(u(i,j-1)+u(i,j)))^2)/(Px(i,j+1)-Px(i,j));
           advection_y = (0.25*(u(i-1,j)+u(i,j))*(v(i-1,j)+v(i-1,j+1)) - 0.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))) / (Yy(i,j) - Yy(i-1,j));
           diffusion =  (1/Re) * ...
           (((u(i, j+1) - u(i, j)) / (Xx(i, j+1) - Xx(i, j)) - (u(i, j) - u(i, j-1)) / (Xx(i, j) - Xx(i, j-1))) / (Px(i, j+1) - Px(i, j)) + ...
           ((u(i+1, j) - u(i, j)) / (Xy(i+1, j) - Xy(i, j)) - (u(i, j) - u(i-1, j)) / (Xy(i, j) - Xy(i-1, j))) / (0.5 * (Xy(i+1, j) - Xy(i-1, j))));
           u_new(i,j)= u(i,j)+dt*(diffusion - (advection_x + advection_y));
      end
   end
   
%x-momentum eq. boundary
u_new(1,:) = 2 - u_new(2,:);    %uN --- %u_new(2,:) vua solve o buoc tren
u_new(nx+1,:) = - u_new(nx,:);    %uS --- %u_new(n,:) vua solve o buoc tren
u_new(2:nx,1) = 0;           %uW
u_new(2:nx,ny) = 0;       %uE


%y-momentum eq. -interior
for i = 2 : nx-1 % V internal velocity co 3 hang` va 4 cot
    for j = 2 : ny
        advection_y = ((0.5*(v(i-1,j)+v(i,j)))^2-((0.5*(v(i,j)+v(i+1,j)))^2))/(Py(i+1,j)-Py(i,j)); 
        advection_x = ((0.25*(v(i,j+1)+v(i,j))*(u(i+1,j)+u(i,j)))-(0.25*(v(i,j)+v(i,j-1))*(u(i+1,j-1)+u(i,j-1))))/(Xx(i,j)-Xx(i,j-1));
        diffusion = (1/Re) * ...
    ( ( (v(i+1,j) - v(i,j)) / (Yy(i+1,j) - Yy(i,j)) - (v(i,j) - v(i-1,j)) / (Yy(i,j) - Yy(i-1,j)) ) / (Py(i,j) - Py(i-1,j)) + ...
    ( ( (v(i,j+1) - v(i,j)) / (Yx(i,j+1) - Yx(i,j)) - (v(i,j) - v(i,j-1)) / (Yx(i,j) - Yx(i,j-1)) ) / (0.5 * (Yx(i,j+1) - Yx(i,j))) ));
        v_new(i,j)= v(i,j)+ dt*(diffusion - (advection_x + advection_y));
    end
end

%y-momentum eq. boundary
v_new(1,2:ny) = 0;           %vN
v_new(nx,2:ny) = 0;         %vS
v_new(:,1) = - v_new(:,2);      %vW --- %v_new(:,2) vua solve o buoc tren
v_new(:,ny+1) = - v_new(:,ny);    %vE --- %v_new(:,n) vua solve o buoc tren



%Continuity eq. - interior
for i = 2:nx
    for j = 2:ny
        
    M = (-1/(Px(i,j+1) - Px(i,j))) + (-1/(Px(i,j) - Px(i,j-1))) + ...
    (Xx(i,j) - Xx(i,j-1)) / (Yy(i,j) - Yy(i-1,j)) * ...
    ((-1/(Py(i+1,j) - Py(i,j))) + (-1/(Py(i,j) - Py(i-1,j))));
   
    N = (-p(i,j+1)/(Px(i,j+1) - Px(i,j))) + (-p(i,j-1)/(Px(i,j) - Px(i,j-1))) + ...
    (Xx(i,j) - Xx(i,j-1)) / (Yy(i,j) - Yy(i-1,j)) * ...
    ((-p(i+1,j)/(Py(i+1,j) - Py(i,j))) + (-p(i-1,j)/(Py(i,j) - Py(i-1,j))));
    
    D = (1/dt) * ((u_new(i,j) - u_new(i,j-1)) / (Xx(i,j) - Xx(i,j-1)) + (v_new(i,j) - v_new(i-1,j)) / (Yy(i,j) - Yy(i-1,j)));
        
    p_new(i,j) = (N-(Xx(i,j)-Xx(i,i-1))*D)/M;
    end
end


%Continuity eq. - Boundary
p_new(1,:) = p_new(2,:);            %pN
p_new(nx+1,:) = p_new(nx,:);        %pS
p_new(:,1) = p_new(:,2);            %pW
p_new(:,nx+1) = p_new(:,ny);        %pE

%Update Velocity
% error_u = 0;
% sigma_u = 0;
for i = 1 : nx
    for j = 1: ny -1
        u(i,j) = u_new(i,j)-dt*(p_new(i,j+1)-p_new(i,j))/(Px(i,j+1)-Px(i,j));
%         error_u = error_u + (u(i,j)-u_new(i,j));
%         sigma_u = sigma_u + u(i,j);
    end
end    


% error_v = 0;
% sigma_v = 0;
for i = 1 : nx -1
    for j = 1: ny
        v(i,j) = v_new(i,j)-dt*(p_new(i+1,j)-p_new(i,j))/(Py(i+1,j)-Py(i,j));
%         error_v = error_v + (v(i,j)-v_new(i,j));
%         sigma_v = sigma_v + v(i,j);
    end
end
% error = 0;
% error = error + error_u/sigma_u + error_v/sigma_v;

error = 0; 
for i = 2 : nx-1       %internal U and V
    for j = 2 : ny-1
        error = error + abs((u_new(i,j)-u_new(i,j-1))/hx + (v_new(i-1,j)-v_new(i,j))/hy);
    end
end

%Error monitoring after every few timstep
if(rem(iteration,2))==0
    figure(1);
    semilogy(iteration, error, '-ko');
    hold on
    xlabel('Iteration');
    ylabel('Residual Error');
end

iteration = iteration + 1;
end
%%
%Create the contour
x_dom = ((1:nx)-1).*hx;
y_dom = length_y-((1:ny)-1).*hy;
[X,Y]= meshgrid(x_dom,y_dom);
figure(21);
vel_final = sqrt(u_final.^2+v_final.^2);
contourf(X',Y',vel_final,21,'LineStyle','none');
colorbar
colormap('jet')
xlabel('x')
ylabel('y')
figure(21);
hold on
quiver(X',Y',u_final,v_final,10,'k-')
     xlim([0,length_y]);
     ylim([0,length_x]);
