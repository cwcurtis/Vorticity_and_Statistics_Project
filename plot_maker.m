function plot_maker

Llx = 5*pi;
K = 512;
dk = pi/Llx;
ep = .05;
dt = 2e-3;
om = 1.12;
sig = 0;
Nens = 512;
k0 = 1;
[Om,cg,ad,anl] = param_maker(k0,om,sig);
mwid = sqrt(anl/ad);
disp('Timescale is')
disp(1/anl)
disp('Minimal Stable Envelope Width')
disp(mwid)
tf = ceil(max(1.5/abs(anl),1));

%{
nowdths = 16;
widths = linspace(1.05*mwid,2*sqrt(2)*mwid,nowdths);
nls_avgs = zeros(nowdths,1);
stds_nls = zeros(nowdths,1);
ktss_nls = zeros(nowdths,1);
bfi_nls = zeros(nowdths,1);
dysthe_avgs = zeros(nowdths,1);
stds_dysthe = zeros(nowdths,1);
ktss_dysthe = zeros(nowdths,1);
bfi_dysthe = zeros(nowdths,1);
int_avgs = zeros(nowdths,1);
stds_int = zeros(nowdths,1);
ktss_int = zeros(nowdths,1);
for jj=1:nowdths
    [nls_avg,nls_std,nls_kts,nls_bfi,dysthe_avg,dysthe_std,dysthe_kts,dysthe_bfi,int_avg,int_std,int_kts] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,widths(jj),k0,Nens);
    nls_avgs(jj) = nls_avg;
    stds_nls(jj) = nls_std;
    ktss_nls(jj) = nls_kts;
    bfi_nls(jj) = mwid/widths(jj);
    dysthe_avgs(jj) = dysthe_avg;
    stds_dysthe(jj) = dysthe_std;
    ktss_dysthe(jj) = dysthe_kts;
    bfi_dysthe(jj) = mwid/widths(jj);
    int_avgs(jj) = int_avg;
    stds_int(jj) = int_std;
    ktss_int(jj) = int_kts;
end

[pnls,Snls] = polyfit(bfi_nls,ktss_nls,2);
[pdysthe,Sdysthe] = polyfit(bfi_dysthe,ktss_dysthe,2);
fnls = pnls(1)*bfi_nls.^2 + pnls(2)*bfi_nls + pnls(3);
fdysthe = pdysthe(1)*bfi_dysthe.^2 + pdysthe(2)*bfi_dysthe + pdysthe(3);

disp('Leading p coefficients are:')
disp([pnls(1) pdysthe(1)])

plot(bfi_nls,ktss_nls,'k*',bfi_nls,fnls,'k-.',bfi_dysthe,ktss_dysthe,'ko',bfi_dysthe,fdysthe,'k:','LineWidth',2)
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlabel('BFI','Interpreter','LaTeX','FontSize',30)
ylabel('Kurtosis','Interpreter','LaTeX','FontSize',30)
legend({'$NLS$','$NLS~Fit$','$Dysthe$','$Dysthe~Fit$'},'Interpreter','LaTeX','FontSize',30)

%}

widths = 1.05*mwid;
%widths = 2*sqrt(2)*mwid;
disp(widths)
[~,~,~,~,~,~,~,~,~,~,~] = nls_Dysthe_comparison(Llx,K,ep,tf,dt,om,sig,widths,k0,Nens);

%{
figure(1)
%tstr = sprintf('$\omega=$%f, $\epsilon=$%f', om, ep);
%title(tstr)

subplot(3,1,1);
hold on
plot(widths,int_avgs,'k-',widths,nls_avgs,'k--',widths,dysthe_avgs,'k:','LineWidth',2)
ax = gca;
line([mwid mwid],get(ax,'YLim'),'Color',[0 0 0],'LineWidth',2)
hold off
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
set(gca,'xtick',[])
xlim([min(widths) max(widths)])
ylabel('Mean','Interpreter','LaTeX','FontSize',30)

subplot(3,1,2);
hold on
plot(widths,stds_int,'k-',widths,stds_nls,'k--',widths,stds_dysthe,'k:','LineWidth',2)
ax = gca;
line([mwid mwid],get(ax,'YLim'),'Color',[0 0 0],'LineWidth',2)
hold off
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
set(gca,'xtick',[])
xlim([min(widths) max(widths)])
ylabel('Std. Dev.','Interpreter','LaTeX','FontSize',30)

subplot(3,1,3);
hold on
plot(widths,ktss_int,'k-',widths,ktss_nls,'k--',widths,ktss_dysthe,'k:','LineWidth',2)
ax = gca;
line([mwid mwid],get(ax,'YLim'),'Color',[0 0 0],'LineWidth',2)
hold off
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlim([min(widths) max(widths)])
xlabel('$\sigma$','Interpreter','LaTeX','FontSize',30)
ylabel('Kurtosis','Interpreter','LaTeX','FontSize',30)

legend({'$Initial$','$NLS$','$Dysthe$'},'Interpreter','LaTeX','FontSize',30)
%}