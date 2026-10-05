
%% 

% Note that here use the left hemisphere only
nregs=308;
nregs_lh=152;

%load data
gene_regional_expression = readmatrix('expression.csv');
X=gene_regional_expression(1:nregs_lh,:); % Predictors
Y_data=readmatrix('stat.csv'); % Response variable
Y=Y_data(1:nregs_lh);

% z-score:
X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
dim=15;
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); % no need to do this but it keeps outputs tidy

%%% plot correlation of PLS component 1 with t-statistic (from Cobre as an example):
figure
plot(XS(:,1),Y,'r.')
[R0,p0]=corrcoef(XS(:,1),Y) 
xlabel('XS scores for PLS component 1','FontSize',14);
ylabel('t-statistic- lh','FontSize',14);
grid on

% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:

rep=5000;
for dim=1:4
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);    
    for j=1:rep        
        %j        
        order=randperm(size(Y,1));        
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=cumsum(100*PCTVAR(2,1:dim));        
        Rsq(j) = temp(dim);    
    end
dim
Rsq_cum(dim)=Rsquared
psq_cum(dim)=length(find(Rsq>=Rsquared))/rep
end
figure
plot(1:dim, psq_cum,'ok','MarkerSize',4,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on

dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
%% 



% Bootstrap to get the gene list:

genes = readcell("gene_name.csv");  % this needs to be imported first
geneindex=1:15632;

%number of bootstrap iterations:
bootnum=5000;

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);

%print out results
csvwrite('PLS1_ROIscores.csv',XS(:,1));
csvwrite('PLS2_ROIscores.csv',XS(:,2));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

numGenes1 = length(PLS1w);
pvals1 = zeros(numGenes1, 1);
for k = 1:numGenes1
    [h,pz1,ci1,zval]=ztest(PLS1w(k),0,PLS1sw(k));
    pvals1(k)=pz1;
end

numGenes2 = length(PLS2w);
pvals2 = zeros(numGenes2, 1);
for k = 1:numGenes2
    [h,pz2,ci2,zval]=ztest(PLS2w(k),0,PLS2sw(k));
    pvals2(k)=pz2;
end
% pvals1 = 2 * (1 - normcdf(abs(temp1)));
% pvals2 = 2 * (1 - normcdf(abs(temp2)));
pvals1_fdr = mafdr(pvals1, 'BHFDR', true);
pvals2_fdr = mafdr(pvals2, 'BHFDR', true);

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
pvals1_fdr = pvals1_fdr(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
pvals2_fdr = pvals2_fdr(ind2);

%print out results
% later use first column of these csv files for pasting into GOrilla (for
% bootstrapped ordered list of genes) 
fid1 = fopen('PLS1_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid1, '%s, %d, %f, %f\n', PLS1{i}, geneindex1(i), Z1(i), pvals1_fdr(i));
end
fclose(fid1)

fid2 = fopen('PLS2_geneWeights.csv','w')
for i=1:length(genes)
  fprintf(fid2, '%s, %d, %f, %f\n', PLS2{i}, geneindex2(i), Z2(i), pvals2_fdr(i));
end
fclose(fid2)
