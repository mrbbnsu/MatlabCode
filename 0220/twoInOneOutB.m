%% 二输入 L，P - F
%% 完整算法
%% 1130 增加插值算法 -- 插值不如ANN
%% 1201 新理解
%% 1208 两个数据混合
clear
clc
clf
istrain =0;
isinter = 0;
layer = [15 10];
fileID = fopen('0215_105308.txt');
% data=fread(fo,'int');
data=fread(fileID,'int');
data=reshape(data',18,[])';

mode = data(:,13);
l = data(:,1);
f = 5e4-data(:,11);
f = smooth(f,15);
p = data(:,9); % set value
P = data(:,10);% Actual value
adcpressure = data(:,12);
hyslmax= data(:,4);
hyslmin = data(:,5);
hyspmax = data(:,14);
hyspmin = data(:,15);
dir = data(:,6);
index = (l>=hyslmax)&(mode~=8);
l(index)=hyslmax(index);
index = (l<=hyslmin)&(mode~=8);
l(index)=hyslmin(index);
index = (p==hyspmin-50) & mode ==41;
hyspmin(index) =p(index);
%% 上升曲线
index = find((mode==4 & dir==0) | (mode ==5 &dir==0));
plot(p)
hold on
plot(l)
% plot3(p,l,f,'o')
plot(index,l(index),'o')
hold on
plot(index,p(index),'o')

% plot(f(index))
% hold off
% plot3(p(index),l(index),f(index),'o')

if(istrain)
    % [fit1, gof] = fit( [l(index) p(index)],f(index), 'linearinterp', 'Normalize', 'on' );
    [~,~,~,fit1,~]=MLAnalysis([l(index) P(index)],f(index),layer);
    save fit1 fit1
else
    load fit1
end

% plot(fit1)
% error
%% 定气压 下降线
index =( find(mode==4));
% plot3(p,l,f,'o')
clf
plot(l(index))
hold on
plot(p(index))

index = find(mode==4 & dir==1 & p==hyspmax);
plot(index,l(index),'o')
hold on
plot(index,p(index),'o')
F2FixP=[hyslmax(index) l(index) hyspmax(index) f(index)];
if(isinter)
    %     F2P = @(x)scatteredInterpolant(F2FixP(:,1:3),F2FixP(:,end),x,'nearest');
    F2P = scatteredInterpolant(F2FixP(:,1:3),F2FixP(:,end));
else
    if(istrain)
        % [fit1, gof] = fit( [l(index) p(index)],f(index), 'linearinterp', 'Normalize', 'on' );
        [~,~,~,F2P,~]=MLAnalysis(F2FixP(:,1:3),F2FixP(:,end),layer);
        save F2P F2P
    else
        load F2P
    end
end
% FF=F2P(F2FixP(:,1),F2FixP(:,2),F2FixP(:,3) );
% error
% clf
% plot(FF)
% hold on
% plot(f(index),'o')
% plot(f(index))
%% 定长度 下降线
index =( find(mode==5));
% plot3(p,l,f,'o')
clf
plot(index,l(index)/10)
hold on
plot(index,p(index))
plot(index,P(index))
index = find(mode==5 & dir==1 & l ==hyslmax & p<hyspmax);
plot(index,l(index)/10,'o')
hold on
plot(index,p(index),'o')
plot(index,P(index),'o')
F2FixL=[hyspmax(index) P(index) hyslmax(index) f(index)];
if(isinter)
    F2L= scatteredInterpolant(F2FixL(:,1:3),F2FixL(:,end));
    %     F2L = @(x)griddatan(F2FixL(:,1:3),F2FixL(:,end),x,'nearest');
else
    %  [ANNY3,LinearY,r3,F2L,B3]=MLAnalysis(F2FixL(:,1:3),F2FixL(:,end),[10 6]);
    if(istrain)
        % [fit1, gof] = fit( [l(index) p(index)],f(index), 'linearinterp', 'Normalize', 'on' );
        [~,~,~,F2L,~]=MLAnalysis(F2FixL(:,1:3),F2FixL(:,end),layer);
        save F2L F2L
    else
        load F2L
    end
end
%  F2L = @(palpha,pbeta,length)griddatan(F2FixL(:,1:3),F2FixL(:,end),[palpha,pbeta,length],'nearest');
% FF=F2P(F2FixP(:,1),F2FixP(:,2),F2FixP(:,3) );
% clf
% plot(FF)
% hold on
% plot(f(index),'o')
% plot(f(index))




if(isinter)
    fun2u = @(lalpha,lbeta,pressure) F2P( [lalpha,lbeta,pressure]);
    fun2v = @(palpha,pbeta,length) F2L([palpha,pbeta,length]);
    fun3u = @(lalpha,lbeta,pressure) F3P([lalpha,lbeta,pressure]) ;%- F2P([max(l),lbeta,pressure]')
    fun3v = @(palpha,pbeta,length) F3L([palpha,pbeta,length]) ;%- F2L([max(p),pbeta,length]')
else
    fun2u = @(lalpha,lbeta,pressure) F2P( [lalpha,lbeta,pressure]');
    fun2v = @(palpha,pbeta,length) F2L([palpha,pbeta,length]');
    fun3u = @(lalpha,lbeta,pressure) F3P([lalpha,lbeta,pressure]') ;%- F2P([max(l),lbeta,pressure]')
    fun3v = @(palpha,pbeta,length) F3L([palpha,pbeta,length]') ;%- F2L([max(p),pbeta,length]')

    fun2u = @(lalpha,lbeta,pressure) calANN(F2P,[lalpha,lbeta,pressure]');
    fun2v = @(palpha,pbeta,length) calANN(F2L,[palpha,pbeta,length]');
%     fun3u = @(lalpha,lbeta,pressure) calANN(F3P,[lalpha,lbeta,pressure]') ;%- F2P([max(l),lbeta,pressure]')
%     fun3v = @(palpha,pbeta,length) calANN(F3L,[palpha,pbeta,length]') ;%- F2L([max(p),pbeta,length]')
end
fun1 = @(l,pressure)calANN(fit1,[l,pressure]');
%% 验证数据
% data=importdata("data13.txt");
% mode = data(:,2);
% l = data(:,1);
% f = 5e4-data(:,13);
% p = data(:,11); % set value
% P = data(:,12);% Actual value
% adcpressure = data(:,12);
% hyslmax= data(:,4);
% hyslmin = data(:,5);
% hyspmax = data(:,14);
% hyspmin = data(:,15);
% dir = data(:,6);
index =( mode==8);
l = l(index);
p = p(index);
P = P(index);
f= f(index);
% l(index)=hyslmax(index);
% index = (l<=hyslmin)&(mode~=8);
% l(index)=hyslmin(index);
% index = (p==hyspmin-50) & mode ==41;
% hyspmin(index) =p(index);
index = find(mode==8);
clf
% plot(l)
% hold on
% l=smooth(l);

l(l<0)=0;
% u = l(index);
% v = p(index);
% actualf=f(index);
u=smooth(l);
v=p;
actualf =f;
subplot 311
u=u-0;
plot([u v])
hold on

% 标记转折点
[~,uMax]=findpeaks(u);
[~,vMax]=findpeaks(v);
[~,umin]=findpeaks(-u);
[~,vmin]=findpeaks(-v);
plot(uMax,u(uMax),'o',vMax,v(vMax),'o')
plot(umin,u(umin),'o',vmin,v(vmin),'o')
vm=1;
vM=1;
uM=1;
um=1;
v = smooth(p);
fcal=[];
uvL=[];
udir=1;
vdir=1;
DF(1)=0;
DFV=[];
DFU=[];
USum = 0;
VSum = 0;
orderflag = 0;
tic
for i=1:2000
    unM = length(uM);
    unm = length(um);
    vnM = length(vM);
    vnm = length(vm);
    uvMn(i,:)=[unM unm vnM vnm];
    nDF = size(uvL,1);

    if(i>2)
        u1 = u(i-1);
        v1 = v(i-1);
        u2 = u(i);
        v2 = v(i);
    end
    if(udir==1)
        uM(end)=i;
    else
        um(end)=i;
        uM(end)=i;
    end
    if(vdir==1)
        vM(end)=i;
    else
        vm(end)=i;
        vM(end)=i;
    end
    DIR(i,:)=[udir vdir udir*10+vdir];
    if(size(uvL,1)>1)
        udirlast = uvL(end,1);
        vdirlast = uvL(end,2);
    end
    DF = fun1(u(i),v(i));

    if(size(DFU)>0)
        % 定v变u v定的是上一个转折点的值
        % V还没转折过,lastvert 是当前点
        ulastvert = uvL(end,3);
        if(size(DFV,1)==0)
            vlastvert=i;
        else
            vlastvert = ulastvert;
        end
        %         if(udir==-1) %% 下降
        %
        %         else   %% 上升
        %             DFU(end) = -fun3u(u(um(end)),u(ulastvert),v(i))+fun3u(u(um(end)),u2,v(i)); %% (u430,v)--(u,v) 定v u增加
        %         end
        for uk = length(uM)-1: length(uM)-1
            DFU(uk) = fun2u(u(uM(uk)),u(um(uk+1)),v(i))-fun2u(u(uM(uk+1)),u(um(uk+1)),v(i));%% 定v u减小
        end

    end
    if(size(DFV)>0)

        lastvert = uvL(end,3);
        if(size(DFU,1)==0)
            lastvert=i;
        end

        %         if(vdir==-1)
        % %             DFV(end) = -fun2v(v(vM(end)),v(lastvert),u2)+fun2v(v(vM(end)),v2,u2);%% 定v u减小
        %
        %
        %         else
        %             DFV(end) = -fun3v(v(vm(end)),v(lastvert),u2)+fun3v(v(vm(end)),v2,u2); %% (u430,v)--(u,v) 定v u增加
        %         end
        %         DFV(end) = fun2v(v(vM(1)),v(vm(2)),u(i))-fun2v(v(vM(2)),v(vm(2)),u(i));%% 定v u减小
        for vk = length(vM)-1: length(vM)-1
            DFV(vk) = fun2v(v(vM(vk)),v(vm(vk+1)),u(i))-fun2v(v(vM(vk+1)),v(vm(vk+1)),u(i));%% 定v u减小
        end


    end

    fcal(i) = sum(DF)+sum(DFV)+sum(DFU);
    USum(i) = sum(DFU);
    VSum(i) = sum(DFV);

    if(ismember(i,uMax))
        uvL=[uvL;udir vdir i 1];
        um(end+1)=0;
        uM(end+1)=0;
        udir =-1; % 下降
        DFU(end+1)=0;
        if(size(DFV,1)>0)
            %             DFV(end+1)=0;
        end
    elseif(ismember(i,umin))
        uvL=[uvL;udir vdir i 1];
        %         uM(end+1)=0;
        udir =1; % 上升
        %         DFU(end+1)=0;
        if(size(DFV,1)>0)
            %             DFV(end+1)=0;
        end
    end
    if(ismember(i,vMax))
        uvL=[uvL;udir vdir i 2];
        vm(end+1)=0;
        vM(end+1)=0;
        vdir =-1; % 下降
        DFV(end+1)=0;
        if(size(DFU,1)>0)
            %             DFU(end+1)=0;
        end
    elseif(ismember(i,vmin))
        uvL=[uvL;udir vdir i 2];
        %         vM(end+1)=0;
        vdir =1; % 上升
        %
        %         DFV(end+1)=0;
        if(size(DFU,1)>0)
            %             DFU(end+1)=0;
        end
    end



end
toc
subplot 312



fcal=smooth(fcal);
plot(fcal,'k*')
hold on
plot(actualf)

actualf=actualf(1:length(fcal));
rms(fcal-actualf)
corrcoef([fcal actualf])
difffcal = smooth(diff(fcal),9);
fcalFilter = cumsum(difffcal)+fcal(1);
% plot(fcalFilter,'bs')
subplot 313

difff=diff(f);
vindex = [uMax; vMax ;umin ;vmin]+1;
plot(difffcal)
hold on
plot(difff)
plot(vindex,difff(vindex),'ko')


figure


plot(fcal,'k','LineWidth',2)
hold on
plot(actualf,'LineWidth',2)
grid on
legend('Predicted value','Actual value')
title('Force')

set(gcf,'position',[100,100,1200,500]);


% plot(umin,u(umin),'o',vmin,v(vmin),'o')
error
ffit = fit1(p(index),l(index));

error
clf
ffit = fillmissing(ffit,'movmedian',100);
plot(ffit-ffit(1)+1800)
hold on
plot(f(index))
% xlim([7000 8000])
%