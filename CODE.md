clear 
clc
format long
 
patch 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%信号的长度
n   =  100;    
%感知矩阵A的行数m >= k*log(n/k)
m   =  64;                  
%信号的稀疏度 s 从 5, 10, 15, 20, 25
s    =  12;                      %数据中非零元素的个数;

% 控制矩阵的一致性(Coherence)的参数F = 5, 10, 15, 20
F       =   10;
% 感知矩阵的构造
A           =   (1/sqrt(m)) .* cos((2*pi/ F) .* rand(m,n));

% 构造原始的稀疏信号
xtrue             =  zeros(n,1);         % 原始信号；
Indtrue         =  randperm(n);     %随机生成脚标数；
%产生出只有s个非零元素的原始信号；
xtrue(Indtrue(1:s))  =   randn(s,1) ;   
% 构造观测数据
b               =   A*xtrue;

% 最小二范数解及其二范数的数值
Aplus            =     pinv(A);
PA                 =     eye(n) - Aplus*A;
xplus             =     Aplus*b;
L2xplus         =     norm(xplus,2)^2;

% 初始化迭代
wj_XNext                           =         zeros(n,1);
wj_YNull                            =         zeros(n,1);
wj_rho_Y                            =         100;    % 这个参数的初始化是需要给出理论上的数值范围
IterMax                              =         10*n;
wj_IterK                                   =         1;
% Lagrangian 向量
wj_LambdaXPrev               =         zeros(m,1);
wj_LambdaYPrev               =         zeros(m,1);
wj_LambdaPlusPrev           =         zeros(n,1);
% 终止判断条件
wj_TerminalCond                   =          1.0;
wj_RelErr                                =          1e-8;

% 迭代过程产生的相对误差
wj_RelErr_X2True                                          =   zeros(IterMax,1);
wj_RelErr_Xk2Xkp1                                       =   zeros(IterMax,1);
wj_ObjVal4L1oL2_Xk                                    =   zeros(IterMax,1);
wj_ObjVal4L1_Xk                                          =   zeros(IterMax,1);
wj_ObjVal4L1oL2_Yk                                    =   zeros(IterMax,1);
wj_ObjVal4L2_YkXk                                      =   zeros(IterMax,1);
wj_ObjVal4L2_Yk                                          =   zeros(IterMax,1);
% 求解命令fsolve的时间
% wj_fsolve_TOC       =   0.0;
wj_cvx_TOC               = 0.0;

% 主要的循环过
tic
while ( wj_IterK <= IterMax && wj_TerminalCond > wj_RelErr)
    %% 第一步：更新x^{k+1}
    % 保存上一步的迭代解 x^{k} 和x^{k+1}
    wj_XPrev        =   wj_XNext;
    % X-子问题的参数wj_rho_X 是Y-子问题的参数wj_rho_Y
    wj_rho_X         =   wj_rho_Y;
    % 迭代产生的wj_u^{k}
    wj_UIterK        =   wj_YNull + xplus - (1/wj_rho_X).*wj_LambdaPlusPrev;
    % 求解x^{k+1}
    %{
    6 solvers initialized (* = default):
    Gurobi     9.00       {cvx}/gurobi/a64
 *  Mosek      9.1.9      {cvx}/mosek/a64
    Mosek_2    9.1.9      /home/askuyue/wMWork/MATLAB-Add-Ons/mosek/7/toolbox/r2013a
    Mosek_3    9.1.9      /home/askuyue/wMWork/MATLAB-Add-Ons/mosek/7/toolbox/r2013a
    SDPT3      4.0        {cvx}/sdpt3
    SeDuMi     1.3.4      {cvx}/sedumi
    %}
%        cvx_solver  gurobi  
%         cvx_solver  mosek
%        cvx_save_prefs
%      cvx_solver sedumi 
        tic
        cvx_begin  %quiet
                variable xcvx(n);
                variable xcvxup(n);
                minimize( norm(xcvx,1)/sqrt( norm(wj_YNull,2)^2 + L2xplus) + (wj_rho_X/2)*sum( (xcvx - wj_UIterK).^2) );
                subject to 
                     A*xcvx  == b; 
%                  minimize( sum(xcvxup)  + (wj_rho_X * sqrt( norm(wj_YNull,2)^2 + L2xplus)/2)*sum( (xcvx - wj_UIterK).^2) );
%                           xcvx      >= -xcvxup;
%                           xcvx      <= xcvxup;
%                           xcvxup  >= zeros(n,1);
        cvx_end
        wj_cvx_TOC              =           toc + wj_cvx_TOC;
        wj_XNext                  =           xcvx;
%% solver_sBP in TFOCS
%             wj_mu                       =            wj_rho_X * sqrt( norm(wj_YNull,2)^2 + L2xplus);
%             wj_XNext                  =             solver_sBP(A,b,wj_mu,wj_UIterK,zeros(m,1));
 
    %% 第二步：更新y^{k+1}
    % 更新wj_rho_Y  
    wj_rho_Y    =   30*norm(wj_XNext,1)/norm(xplus,2)^3;
    % 迭代产生的wj_v^{k}
    wj_VIterK   =   wj_XNext - xplus + (1/wj_rho_Y).*wj_LambdaPlusPrev;
    % 求解y^{k+1} 近似的线性方程组求解
     wj_YNull   =   (1 - norm(wj_XNext,1)/( wj_rho_Y* norm(wj_XNext,2)^3 ))^(-1) .* (PA * wj_VIterK);         
%      %% levenberg-marquardt method
%      wj_Fk                  =        @(y) (wj_rho_Y - norm(wj_XNext,1)/ norm(wj_XNext,2)^3 ).* y - wj_rho_Y .* (PA * wj_VIterK);
%      wj_YNull                  =        zeros(n,1);
%      tic
%      wj_options         =         optimoptions('fsolve','Display','iter-detailed','algorithm','levenberg-marquardt');
%      wj_YNull             =         fsolve(wj_Fk,wj_YNull,wj_options);
%      wj_fsolve_TOC    =         toc + wj_fsolve_TOC;
     
      %% 更新Lagrange乘子Lambda
       wj_LambdaPlusPrev           =           wj_LambdaPlusPrev +  wj_rho_Y.*(wj_XNext - wj_YNull - xplus);
       
      %% 计算迭代过程中的相对误差
      % 第k步的解和真实解之间的相对误差
      wj_RelErr_X2True(wj_IterK)                                    =               norm(wj_XNext - xtrue,2)/norm(xtrue,2);
      % 第k步的解和第k+1之间的相对误差
      wj_RelErr_Xk2Xkp1(wj_IterK)                                  =              norm(wj_XNext - wj_XPrev,2) / norm(wj_XPrev,2);  
      % 计算第k步的L1/L2的关于X的比值
      wj_ObjVal4L1oL2_Xk(wj_IterK)                               =              norm(wj_XNext,1)/norm(wj_XNext,2);
      % 计算第k步的L1范数的关于X的数值
      wj_ObjVal4L1_Xk(wj_IterK)                                     =              norm(wj_XNext,1);
      % 计算第k步的L1/L2的关于变量Y的比值
      wj_ObjVal4L1oL2_Yk(wj_IterK)                               =               norm(wj_YNull,1)/norm(wj_YNull,2);
      % 计算第k步的Xk和Yk之间的相对误差
      wj_ObjVal4L2_YkXk(wj_IterK)                                 =               norm(wj_XNext - wj_YNull,2);
      % 计算第k步的Yk的二范数
      wj_ObjVal4L2_Yk(wj_IterK)                                     =               norm(wj_YNull,2);
       %% 更新迭代参数和迭代次数
       wj_TerminalCond                    =           wj_RelErr_Xk2Xkp1(wj_IterK);
       wj_IterK                                   =           wj_IterK + 1;
               
end
wj_ADMMdp_TOC   =   toc;
wj_IterK       =   wj_IterK -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Scale-invariant approach for sparse signal recovery
% Rahimi, Wang, Dong, Lou : rwdl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 初始化迭代
rwdl_XNext                            =         zeros(n,1);
rwdl_YNext                            =         zeros(n,1);
rwdl_ZNext                            =         zeros(n,1);
rwdl_rho_1                             =         wj_rho_Y;    % 这个参数的初始化是需要给出理论上的数值范围
rwdl_rho_2                             =         wj_rho_Y;    % 这个参数的初始化是需要给出理论上的数值范围
rwdl_IterK                               =         1;
% Lagrangian 向量
rwdl_LambdaVPrev                =         zeros(n,1);
rwdl_LambdaWPrev               =         zeros(n,1);
% 终止判断条件
rwdl_TerminalCond                   =          1.0; 
rwdl_RelErr                                =           wj_RelErr;

% 迭代过程产生的相对误差
rwdl_RelErr_X2True                                          =    zeros(IterMax,1);
rwdl_RelErr_Xk2Xkp1                                       =    zeros(IterMax,1);
rwdl_RelErr_Xk2Yk                                           =    zeros(IterMax,1);
rwdl_RelErr_Xk2Zk                                           =    zeros(IterMax,1);
rwdl_ObjVal4L1oL2_Xk                                    =    zeros(IterMax,1);
rwdl_ObjVal4L1_Xk                                          =    zeros(IterMax,1);
% 求解命令fsolve的时间
% wj_fsolve_TOC       =   0.0;
rwdl_cvx_TOC               = 0.0;

% 主要的循环过
tic
while ( rwdl_IterK <= IterMax && rwdl_TerminalCond > rwdl_RelErr)
        %% 第一步：更新x^{k+1}
        % 保持上一步的结果
        rwdl_XPrev                     =         rwdl_XNext;
        % 计算中间变量fk
        rwdl_fxk                          =         (rwdl_rho_1/(rwdl_rho_1 + rwdl_rho_2)) .* (rwdl_YNext - rwdl_LambdaVPrev ./ rwdl_rho_1) +...
                                                           (rwdl_rho_2/(rwdl_rho_1 + rwdl_rho_2)) .* (rwdl_ZNext - rwdl_LambdaWPrev ./ rwdl_rho_2);
        % 计算xkp1 利用PA=eye(n) - A'* (A*A')^(-1)*A;
        rwdl_XNext                     =        PA * rwdl_fxk +  xplus;
        %% 第二步：更新y^{k+1}
        % 计算中间变量
        rwdl_ck                             =      norm(rwdl_ZNext,1);
        rwdl_dk                             =      rwdl_XNext + rwdl_LambdaVPrev ./ rwdl_rho_1;
        rwdl_etak                          =      norm(rwdl_dk,2);
        rwdl_Dk                            =      rwdl_ck / (rwdl_rho_1 * rwdl_etak^3);
        rwdl_Ck                             =      ( (27*rwdl_Dk + 2 + sqrt((27*rwdl_Dk+2)^2 - 4))/2 )^(1/3);
        rwdl_taok                          =      1/3 + (rwdl_Ck + 1/rwdl_Ck)/3;
        if rwdl_etak == 0
            rwdl_YNext                   =      randn(n,1);
            rwdl_YNext                   =      rwdl_YNext .* ( (rwdl_ck/rwdl_rho_1)^(1/3) / norm(rwdl_YNext,2));
        else
            rwdl_YNext                   =      rwdl_taok .* rwdl_dk;
        end
        %% 第三步：更新z^{k+1}
        rwdl_MidZ                       =       rwdl_XNext + rwdl_LambdaWPrev ./ rwdl_rho_2;
%         rwdl_ZNext                      =       sign(rwdl_MidZ) .* max(abs(rwdl_MidZ) - 1/(rwdl_rho_2*norm(rwdl_YNext,2)),0);
        tic
        cvx_begin  %quiet
                variable zrwdlcvx(n);
                minimize( norm(zrwdlcvx,1)/norm(rwdl_YNext,2)+ (rwdl_rho_2/2)*sum( (zrwdlcvx - rwdl_MidZ).^2) );
                subject to 
                     A*zrwdlcvx  == b; 
        cvx_end
        rwdl_cvx_TOC                        =         rwdl_cvx_TOC + toc;
        rwdl_ZNext                            =         zrwdlcvx;
        %% 第四步：更新LambdaV W
        rwdl_LambdaVPrev                =         rwdl_LambdaVPrev + rwdl_rho_1 .* (rwdl_XNext - rwdl_YNext) ;
        rwdl_LambdaWPrev               =         rwdl_LambdaWPrev + rwdl_rho_2 .* (rwdl_XNext - rwdl_ZNext);
       %% 计算迭代过程中的算法行为
       % 迭代过程产生的相对误差
        rwdl_RelErr_X2True(rwdl_IterK)                                          =    norm(rwdl_XNext - xtrue,2)/norm(xtrue,2);
        rwdl_RelErr_Xk2Xkp1(rwdl_IterK)                                       =    norm(rwdl_XNext - rwdl_XPrev,2)/norm(rwdl_XPrev,2);
        rwdl_RelErr_Xk2Yk(rwdl_IterK)                                           =    norm(rwdl_XNext - rwdl_YNext,2);
        rwdl_RelErr_Xk2Zk(rwdl_IterK)                                           =    norm(rwdl_XNext - rwdl_ZNext,2);
        rwdl_ObjVal4L1oL2_Xk(rwdl_IterK)                                    =    norm(rwdl_XNext,1)/norm(rwdl_XNext,2);
        rwdl_ObjVal4L1_Xk(rwdl_IterK)                                          =    norm(rwdl_XNext,1);
       %% 更新迭代参数和迭代次数
       rwdl_TerminalCond                    =           norm(rwdl_XNext - rwdl_XPrev,2) / norm(rwdl_XPrev,2);
       rwdl_IterK                                   =           rwdl_IterK + 1;
end
rwdl_ADMMrw_TOC         =       toc;

%{ 
%% Plots for ADMMdp
figure(777)
plot(1:n,xtrue,'rd',1:n,wj_XNext,'g-+','LineWidth',2)
legend('$x^{t}$','$x^{k+1}$','Interpreter','latex')
xlabel('Indices')

figure(1)
subplot(2,2,1)
semilogy(1:wj_IterK,wj_RelErr_X2True(1:wj_IterK),'r-.','LineWidth',2)
legend('ADMMdp - $\frac{||x^k-x^{t}||_2}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(a) Relative errors')

subplot(2,2,2)
semilogy(1:wj_IterK,wj_RelErr_Xk2Xkp1(1:wj_IterK),'m--','LineWidth',2)
legend('ADMMdp - $\frac{||x^{k+1}-x^{k}||_2}{||x^{k}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(b) Relative errors')

subplot(2,2,3)
plot(1:wj_IterK,wj_ObjVal4L1oL2_Xk(1:wj_IterK),'b:','LineWidth',2)
hold on
plot(1:wj_IterK,(norm(xtrue,1)/norm(xtrue,2)).*ones(wj_IterK,1),'c:','LineWidth',2)
legend('ADMMdp - $\frac{||x^{k}||_1}{||x^{k}||_2}$','$\frac{||x^{t}||_1}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(c) Objective functions of L1/L2')

subplot(2,2,4)
plot(1:wj_IterK,wj_ObjVal4L1_Xk(1:wj_IterK),'k-','LineWidth',2)
hold on 
plot(1:wj_IterK,norm(xtrue,1).*ones(wj_IterK,1),'g-','LineWidth',2)
legend('ADMMdp - $||x^{k}||_1$','$||x^{t}||_1$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(d) Objective functions of L1')

figure(2)
subplot(1,2,1)
plot(1:wj_IterK,wj_ObjVal4L1oL2_Yk(1:wj_IterK),'r-.','LineWidth',2)
legend('ADMMdp - $\frac{||y^{k}||_1}{||y^{k}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(e) Objective functions of L1/L2')

subplot(1,2,2)
plot(1:wj_IterK,wj_ObjVal4L2_YkXk(1:wj_IterK),'m--','LineWidth',2)
legend('ADMMdp - $||x^{k}-y^{k}||_2$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(f) Absolute errors') 

figure(3)
plot(1:wj_IterK,wj_ObjVal4L2_Yk(1:wj_IterK),'rd-','LineWidth',2)
legend('ADMMdp - $||y^{k}||_2$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(g) Objective functions of L2')
hold on
plot(1:wj_IterK,wj_ObjVal4L1_Xk(1:wj_IterK) ./wj_ObjVal4L1oL2_Xk(1:wj_IterK),'gs-','LineWidth',2)
legend('ADMMdp - $||x^{k}||_2$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(h) Objective functions of L1') 



%% Plots for mADMMrw
rwdl_IterK              =       rwdl_IterK - 1;
figure(4)
subplot(2,2,1)
semilogy(1:rwdl_IterK,rwdl_RelErr_X2True(1:rwdl_IterK),'r-.','LineWidth',2)
legend('mADMMrw - $\frac{||x^k-x^{t}||_2}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(i) Relative errors')

subplot(2,2,2)
semilogy(1:rwdl_IterK,rwdl_RelErr_Xk2Xkp1(1:rwdl_IterK),'m--','LineWidth',2)
legend('mADMMrw - $\frac{||x^{k+1}-x^{k}||_2}{||x^{k}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(j) Relative errors')

subplot(2,2,3)
plot(1:rwdl_IterK,rwdl_ObjVal4L1oL2_Xk(1:rwdl_IterK),'b:','LineWidth',2)
hold on
plot(1:wj_IterK,(norm(xtrue,1)/norm(xtrue,2)).*ones(wj_IterK,1),'c:','LineWidth',2)
legend('mADMMrw - $\frac{||x^{k}||_1}{||x^{k}||_2}$','$\frac{||x^{t}||_1}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(k) Objective functions of L1/L2')

subplot(2,2,4)
plot(1:rwdl_IterK,rwdl_ObjVal4L1_Xk(1:rwdl_IterK),'k-','LineWidth',2)
hold on 
plot(1:rwdl_IterK,norm(xtrue,1).*ones(rwdl_IterK,1),'g-','LineWidth',2)
legend('mADMMrw - $||x^{k}||_1$','$||x^{t}||_1$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(l) Objective functions of L1')


figure(5)
plot(1:rwdl_IterK,rwdl_RelErr_Xk2Yk(1:rwdl_IterK),'k-','LineWidth',2)
hold on 
plot(1:rwdl_IterK,rwdl_RelErr_Xk2Zk(1:rwdl_IterK),'g-','LineWidth',2)
legend('mADMMrw - $||x^{k} - y^{k}||_1$','$||x^{k} - z^{k}||_1$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(m) Relative errors')

figure(999)
plot(1:n,xtrue,'rd',1:n,rwdl_XNext,'g-+','LineWidth',2)
legend('mADMMrw - $x^{t}$','$x^{k+1}$','Interpreter','latex')
xlabel('Indices')
%}
%% Plots in Paper for ADMMdp and mADMMrw
clc
figure(11111)
subplot(2,2,1)
semilogy(1:wj_IterK,wj_RelErr_X2True(1:wj_IterK),'g-.','LineWidth',2)
hold on
semilogy(1:rwdl_IterK,rwdl_RelErr_X2True(1:rwdl_IterK),'b--','LineWidth',2)
legend('ADMMdp - $\frac{||x^k-x^{t}||_2}{||x^{t}||_2}$','mADMMrw - $\frac{||x^k-x^{t}||_2}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(a) Relative errors')

subplot(2,2,2)
semilogy(1:wj_IterK,wj_RelErr_Xk2Xkp1(1:wj_IterK),'g-.','LineWidth',2)
hold on
semilogy(1:rwdl_IterK,rwdl_RelErr_Xk2Xkp1(1:rwdl_IterK),'b--','LineWidth',2)
legend('ADMMdp - $\frac{||x^{k+1}-x^{k}||_2}{||x^{k}||_2}$','mADMMrw - $\frac{||x^{k+1}-x^{k}||_2}{||x^{k}||_2}$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(b) Relative errors')

subplot(2,2,3)
plot(1:wj_IterK,wj_ObjVal4L1oL2_Xk(1:wj_IterK),'g-.','LineWidth',2)
hold on
plot(1:rwdl_IterK,rwdl_ObjVal4L1oL2_Xk(1:rwdl_IterK),'b--','LineWidth',2)
hold on
plot(1:rwdl_IterK,(norm(xtrue,1)/norm(xtrue,2)).*ones(rwdl_IterK,1),'r-','LineWidth',2)
legend('Interpreter','latex')
legend('ADMMdp - $\frac{||x^{k}||_1}{||x^{k}||_2}$','mADMMrw - $\frac{||x^{k}||_1}{||x^{k}||_2}$',...
            '$\frac{||x^{t}||_1}{||x^{t}||_2}$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(c) Objective functions of L1/L2')

subplot(2,2,4)
plot(1:wj_IterK,wj_ObjVal4L1_Xk(1:wj_IterK),'g-.','LineWidth',2)
hold on 
plot(1:rwdl_IterK,rwdl_ObjVal4L1_Xk(1:rwdl_IterK),'b--','LineWidth',2)
hold on 
plot(1:rwdl_IterK,norm(xtrue,1).*ones(rwdl_IterK,1),'r-','LineWidth',2)
legend('ADMMdp - $||x^{k}||_1$','mADMMrw - $||x^{k}||_1$','$||x^{t}||_1$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(d) Objective functions of L1')


%% AE Plots in Paper for ADMMdp and mADMMrw
figure(22222)
subplot(2,2,1)
plot(1:wj_IterK,wj_ObjVal4L2_YkXk(1:wj_IterK),'g-.','LineWidth',2)
legend('ADMMdp - $||x^{k}-y^{k}||_2$','Interpreter','latex')
xlim([1 wj_IterK])
xlabel('Iteration number')
title('(e) Absolute errors') 

subplot(2,2,2)
semilogy(1:rwdl_IterK,rwdl_RelErr_Xk2Yk(1:rwdl_IterK),'r-','LineWidth',2)
hold on 
semilogy(1:rwdl_IterK,rwdl_RelErr_Xk2Zk(1:rwdl_IterK),'b-.','LineWidth',1)
hold on
semilogy(1:rwdl_IterK,rwdl_RelErr_Yk2Zk(1:rwdl_IterK),'k:','LineWidth',1)
legend('mADMMrw - $||x^{k} - y^{k}||_2$','mADMMrw - $||x^{k} - z^{k}||_2$','mADMMrw - $||y^{k} - z^{k}||_2$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(f) Absolute errors')
 
subplot(2,2,3)
semilogy(1:30,wj_ObjVal4L2_YkXk(1:30),'g-.','LineWidth',2)
legend('ADMMdp - $||x^{k}-y^{k}||_2$','Interpreter','latex')
xlim([1 30])
xlabel('Iteration number')
title('(g) Absolute errors') 

subplot(2,2,4)
semilogy(1:rwdl_IterK,rwdl_RelErr_Xk2Yk(1:rwdl_IterK),'r-','LineWidth',2)
hold on
semilogy(1:rwdl_IterK,rwdl_RelErr_Yk2Zk(1:rwdl_IterK),'k:','LineWidth',2)
legend('mADMMrw - $||x^{k} - z^{k}||_2$','mADMMrw - $||y^{k} - z^{k}||_2$','Interpreter','latex')
xlim([1 rwdl_IterK])
xlabel('Iteration number')
title('(h) Absolute errors')


norm(rwdl_XNext - xtrue,2)
norm(wj_XNext - xtrue,2)


 
