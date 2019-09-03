function plot_maker

Llx = 50*pi;
K = 256;
ep = .1;
dt = 5e-3;
om = -1;
sig = 1e-5;
Nens = 512;
k0 = 1;
[Om,cg,ad,anl] = param_maker(k0,om,sig);
mwid = sqrt(2*anl/ad);
disp('Minimal Stable Envelope Width')
disp(mwid)
tf = mwid;

nowdths = 4;
widths = linspace(.5*mwid,1.5*mwid,nowdths);
nls_avgs = zeros(nowdths,1);
stds_nls = zeros(nowdths,1);
act_avgs = zeros(nowdths,1);
stds_act = zeros(nowdths,1);
int_avgs = zeros(nowdths,1);
stds_int = zeros(nowdths,1);
for jj=1:nowdths
    [nls_avg,std_nls,act_avg,std_act,int_avg,std_int] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,widths(jj),k0,Nens);
    nls_avgs(jj) = nls_avg;
    stds_nls(jj) = std_nls;
    act_avgs(jj) = act_avg;
    stds_act(jj) = std_act;
    int_avgs(jj) = int_avg;
    stds_int(jj) = std_int;
end

%widths = 1.05*mwid;
%disp(widths)
%[nls_avg,std_nls,act_avg,std_act,int_avg,std_int] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,widths,k0,Nens);

figure(1)
%tstr = sprintf('$\omega=$%f, $\epsilon=$%f', om, ep);
%title(tstr)

subplot(2,1,1);
hold on
plot(widths,int_avgs,'k-',widths,nls_avgs,'k--',widths,act_avgs,'k:','LineWidth',2)
ax = gca;
line([mwid mwid],get(ax,'YLim'),'Color',[0 0 0],'LineWidth',2)
hold off
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlabel('$\sigma$','Interpreter','LaTeX','FontSize',30)
ylabel('Mean','Interpreter','LaTeX','FontSize',30)

subplot(2,1,2);
hold on
plot(widths,stds_int,'k-',widths,stds_nls,'k--',widths,stds_act,'k:','LineWidth',2)
ax = gca;
line([mwid mwid],get(ax,'YLim'),'Color',[0 0 0],'LineWidth',2)
hold off
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlabel('$\sigma$','Interpreter','LaTeX','FontSize',30)
ylabel('Standard Deviation','Interpreter','LaTeX','FontSize',30)

legend({'$Initial$','$NLS$','$Dysthe$'},'Interpreter','LaTeX','FontSize',30)