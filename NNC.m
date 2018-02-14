

MESHDEF=500;
%% Part 1A

b=1;
A=[5 b;b 5];
fx=zeros(201)
for x1=-100:100
    for x2=-100:100
        x=[x1;x2];
        fx(x1+101,x2+101)=x'*A*x;
        
    end
    
end
figure
surf(-100:100,-100:100,fx)
title('Surface Plot b=1')
figure
contour(-100:100,-100:100,fx)
title('Contour Plot b=1')

b=5;
A=[5 b;b 5];
fx=zeros(201)
for x1=-100:100
    for x2=-100:100
        x=[x1;x2];
        fx(x1+101,x2+101)=x'*A*x;
        
    end
    
end
figure
surf(-100:100,-100:100,fx)
title('Surface Plot b=5')
figure
contour(-100:100,-100:100,fx)
title('Contour Plot b=5')

b=10
A=[5 b;b 5];
fx=zeros(201)
for x1=-100:100
    for x2=-100:100
        x=[x1;x2];
        fx(x1+101,x2+101)=x'*A*x;
        
    end
    
end
figure
surf(-100:100,-100:100,fx)
title('Surface Plot b=10')
figure
contour(-100:100,-100:100,fx)
title('Contour Plot b=10')

%% Part 1B

b=1;
A=[5 b;b 5];
disp('eigen values of A where b=1')
eig(A)
disp('positive definitive?0?')
[~,p] = chol(A)

b=5;
A=[5 b;b 5];
disp('eigen values of A where b=5')
eig(A)
disp('positive definitive?0?')

[~,p] = chol(A)

b=10;
A=[5 b;b 5];
disp('eigen values of A where b=10')
eig(A)
disp('positive definitive?0?')
[~,p] = chol(A)

%% Part 1C

disp('if one of its eigenvalues of a mat is nonpositive, it must be non-positive definite')
disp('all meet upper-left submatrices Ak have strictly positive determinant')
disp('b=10 a eigenvalue of -5, therefore its not positive')
disp('b=1 yes, b=5 semi b=10 no')
%% Part 1D
disp('PD-bowl, SPD- Onedim, NonPD-saddle')
%% Part 1E

disp('b=1,positive definitive')
disp('f(X)=5x1^2+2x1x2+5x2^2')
disp('f(lX)=l^25x1^2+2l^2x1x2+l^2*5x2^2=l^2f(X)')
disp('b=1, nonhomogenious')

disp('b=5,non positive definitive')
disp('f(X)=5x1^2+10x1x2+5x2^2')
disp('b=10,non positive definitive')
disp('f(X)=5x1^2+20x1x2+5x2^2')


disp('No norms in any case')

%% Part 2A
mu = [0 0];
sigma = eye(2);
R = chol(sigma);
m=repmat(mu,500,1);
data = m+ randn(500,2)*R;

scatter(data(:,1),data(:,2));
pbaspect([1 1 1])

hold on
M=sigma^-1;

% calculate Mah
for i=1:500
    d(i)=sqrt(data(i,:)*M*data(i,:)');

end
x=data(:,1);
y=data(:,2);
z=d';


%plot3(x,y,z,'.','markersize',12)

[xi,yi] = meshgrid(min(x):(max(x)-min(x))/MESHDEF:max(x),min(y):(max(y)-min(y))/MESHDEF:max(y));
zi = griddata(x,y,z,xi,yi);
%surf(xi,yi,zi);
[c,h] = contour(xi,yi,zi);
clabel(c,h);
%% Part 2B

T=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)]*[3 0;0 1]
new_data=data*T;
sigma=cov(new_data);
figure
scatter(new_data(:,1),new_data(:,2));
daspect([1 1 1])
hold on
M=sigma^-1;

% calculate Mah
for i=1:500
    d(i)=sqrt(new_data(i,:)*M*new_data(i,:)');

end
x=new_data(:,1);
y=new_data(:,2);
z=d';


%plot3(x,y,z,'.','markersize',12)

[xi,yi] = meshgrid(min(x):(max(x)-min(x))/MESHDEF:max(x),min(y):(max(y)-min(y))/MESHDEF:max(y));
zi = griddata(x,y,z,xi,yi);
%surf(xi,yi,zi);
[c,h] = contour(xi,yi,zi);
clabel(c,h);

%% Part 2C

A=[5 0;0 2];
disp('A will become the inverse of Cov matrix')
sigma=A^-1;
disp('var of var(x1_new)=0.2, compared to var(x1)=1')
disp('var of var(x2_new)=0.5, compared to var(x2)=1')
disp('according to identity Var(cX)=c^2Var(X)')
cx1=sqrt(.2/1)
cx2=sqrt(.5/1)
T=[cx1 0;0 cx2]

super_new_data=data*T;

figure
scatter(super_new_data(:,1),super_new_data(:,2));
daspect([1 1 1])
hold on
M=sigma^-1;

% calculate Mah
for i=1:500
    d(i)=sqrt(super_new_data(i,:)*A*super_new_data(i,:)');

end
x=super_new_data(:,1);
y=super_new_data(:,2);
z=d';


%plot3(x,y,z,'.','markersize',12)

[xi,yi] = meshgrid(min(x):(max(x)-min(x))/MESHDEF:max(x),min(y):(max(y)-min(y))/MESHDEF:max(y));
zi = griddata(x,y,z,xi,yi);
%surf(xi,yi,zi);
[c,h] = contour(xi,yi,zi);
clabel(c,h);

%% Computer Assignment

load('data.mat')
load('label.mat')
% Below quoted was a test algo, it didn't work that well, better one is
% after
% % need to create 9 classes
% % get all true 1s, build a really noice 1 by avreaging all 1s
% %getting indecies of 1s
% C=zeros(28,28,10);
% for iC=0:9
%     
%     numz=find(labelTrain==iC);
%     num_mod=zeros(28);
%     for i=1:length(numz)
%         num_mod=num_mod+imageTrain(:,:,numz(i));
%     end
%     num_mod=num_mod/length(numz);
%     plotMNIST(num_mod);
%     C(:,:,iC+1)=num_mod;
% end
%
% now we have the models to compare to
% gonna classify now

% TESTID=1
% 
% for TESTID=1:500
% 
%     UUT=imageTest(:,:,TESTID);
%     for iC=0:9
%         dC(iC+1)=sqrt(sum(sum(abs(UUT-C(:,:,iC+1)).^2)));
%     end
%     [M,I] = min(dC);
%     answer(TESTID)=I-1;
% end
% answer=answer';
% quoting it out, the above method yielded 22% err rate, which is bad

%
%% Better Algo
for TESTID=1:500


    UUT=imageTest(:,:,TESTID);

    for iT=1:5000
        dT(iT)=sqrt(sum(sum(abs(UUT-imageTrain(:,:,iT)).^2)));
    end
    [M,I] = min(dT);
    answer(TESTID,1)=M;
    answer(TESTID,2)=I;
    answer(TESTID,3)=labelTrain(I);
end


%%
diff=answer(:,3)-labelTest;
perf=[labelTest diff];


error_ind=find(perf(:,2)~=0);
err_cat=labelTest(error_ind);
    
n=hist(err_cat);

countpercat=hist(labelTest);

err_rate=n./countpercat*100;
plot(0:9,err_rate)
ylabel('error rate')
xlabel('catagory')
%% Part 2
disp('total error rate is')
err_total=sum(n)/sum(countpercat)*100
disp('percent')
%% Part 3
% first error
badboi=error_ind(1)

figure('rend','painters','pos',[10 10 900 500])

subplot(1,2,1)
imshow(imageTest(:,:,badboi))
title('Test Image')
xlabel(sprintf('the pair matched with d=%f',answer(badboi,1)))
subplot(1,2,2)
imshow(imageTrain(:,:,answer(badboi,2)))
title('Matched Image')
xlabel('the 4 matches the shape of the 9 minus the gap. could be consusing for human even')

%% second error
badboi=error_ind(2)

figure('rend','painters','pos',[10 10 900 500])

subplot(1,2,1)
imshow(imageTest(:,:,badboi))
title('Test Image')
xlabel(sprintf('the pair matched with d=%f',answer(badboi,1)))
subplot(1,2,2)
imshow(imageTrain(:,:,answer(badboi,2)))
title('Matched Image')
xlabel('Another 49 mess up, same reason, the general shape is very close')

%% third error
badboi=error_ind(3)

figure('rend','painters','pos',[10 10 900 500])

subplot(1,2,1)
imshow(imageTest(:,:,badboi))
title('Test Image')
xlabel(sprintf('the pair matched with d=%f',answer(badboi,1)))
subplot(1,2,2)
imshow(imageTrain(:,:,answer(badboi,2)))
title('Matched Image')
xlabel('yet anoher one, this one is very bad to human even, just bad writing')

%% forth error
badboi=error_ind(4)

figure('rend','painters','pos',[10 10 900 500])

subplot(1,2,1)
imshow(imageTest(:,:,badboi))
title('Test Image')
xlabel(sprintf('the pair matched with d=%f',answer(badboi,1)))
subplot(1,2,2)
imshow(imageTrain(:,:,answer(badboi,2)))
title('Matched Image')
xlabel('OMG 49 is the death of this method, same reason')

%% fifth error
badboi=error_ind(5)

figure('rend','painters','pos',[10 10 900 500])

subplot(1,2,1)
imshow(imageTest(:,:,badboi))
title('Test Image')
xlabel(sprintf('the pair matched with d=%f',answer(badboi,1)))
subplot(1,2,2)
imshow(imageTrain(:,:,answer(badboi,2)))
title('Matched Image')
xlabel('the first stroke of the 2 almost matched 1 perfectly, even the error from 2nd stroke is not enough to trip')

