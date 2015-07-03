clear all
close all

load BI_para_steric_sim
load steric_sim_trajcode

%I need the full trajectory


jmraw=load('X:\computing\fortran\Yang\BIforBDtirm\newb1_with_hmcfixedpartialw\elecexp\exptirm-2s_1mMnaoh-001.stats.txt');

S_Full=traj_noisedcode;

% first read initialize PI, Aij, Emit
%P_noise=exp(-fel);
%P_noise=P_noise/sum(P_noise);
%P_noise=P_noise';
%P=P';
%Pi=P_noise;
Pi=P;
A=prog;
%A=Ptotal;
[PP DD]=eigs(A,1);
A_T=A';
Emit=pn;
Emit_T=Emit';
nSeg=100;
seqLen=floor(length(S_Full)/nSeg);
EMSteps=1;
noise_strength=sigma;
defaultSigma=noise_strength/(x(2)-x(1));
%intermediate parameter
% pi_up: the nominator for estimating pi
% A_up: the nominator for estimating A (the jump probability matrix)
% mu_up: the nominator for estimating mu
% mu_down: the denominator for estimating mu

for iter=1:EMSteps
    disp(iter)
    %first do the step using the alpha-beta algorithm
    likehood=0.0;
    Pi_up=zeros(nbin,1);
    A_up=zeros(nbin,nbin);
    mu_up=zeros(nbin,1);
    sigma_up=zeros(nbin,1);
    sigma_down=zeros(nbin,1);
    count=zeros(nbin,1);
    for i=1:nSeg
        disp(i)
        startp=1+(i-1)*seqLen;
        endp=i*seqLen;
        S=S_Full(startp:endp);
        for j=1:length(S)
            if j == 1
            alpha(:,j)= Pi.*Emit_T(:,S(j));
            c(j)=sum(alpha(:,j));
            alpha(:,j)=alpha(:,j)/c(j);
            else
            alpha(:,j)=(A*alpha(:,j-1)).*Emit_T(:,S(j));
            c(j)=sum(alpha(:,j));
            alpha(:,j)=alpha(:,j)/c(j);           
            end
        end
        for j=length(S):-1:1
            if j==length(S)
                beta(:,j)=ones(nbin,1);
            else
                beta(:,j)=A_T*(beta(:,j+1).*Emit_T(:,S(j+1)));
                beta(:,j)=beta(:,j)/c(j+1);
            end
            
        end
        gamma=alpha.*beta;
        
        for j=1:length(S)-1
  eta(:,:,j)=1/c(j+1)*(repmat((Emit_T(:,S(j+1)).*beta(:,j+1))',nbin,1).*A_T.*repmat(alpha(:,j),1,nbin));
        end
        
        Pi_up=Pi_up+gamma(:,1);
 %       Pi_down=Pi_down+sum(gamma(:,1));      
        A_up=A_up+sum(eta(:,:,:),3);
        
        mu_up=gamma(:,:)*S;
        mu_down=sum(gamma,2);
        mu=(mu*(i-1)+mu_up./mu_down)/i;
        mu(isnan(mu))=find(isnan(mu));
        
        count=count+sum(gamma,2);
        
        for ii=1:nbin
        
        sigma_up(ii)=sigma_up(ii)+sum(gamma(ii,:)'.*(S-mu(ii)).*(S-mu(ii)));
        sigma_down(ii)=sigma_down(ii)+mu_down(ii);
        end
        % now we accumulate term from short sequence in order to do the M
        % step
        
          tempsum=0.0; 
        for j=1:length(S)
        temp1=gamma(:,j)'*log(Emit_T(:,S(j)));
        temp1(isnan(temp1))=0.0;
        tempsum=tempsum+temp1;
        end
        temp2=sum(eta(:,:,:),3).*log(A);
        %here I set 0*inf=nan =0
        temp2(isnan(temp2))=0.0;
        
        temp3=gamma(:,1).*log(Pi);
        temp3(isnan(temp3))=0.0;
        likehood=likehood+sum(temp3)+sum(sum(temp2))+tempsum;
        
    end
     likehoodset(iter)=likehood;
    % now we are in the M step
      Pi=Pi_up/sum(Pi_up);
      A=A_up;
      for j=1:nbin
          A(:,j)=A(:,j)/sum(A(:,j));
      end
      A_T=A';
      sigma=sigma_up./sigma_down;
      sigma=sqrt(sigma);
      % now construct a emission matrix
      [V, lamb] =eigs(A,1);
      PcalSet(:,iter)=V;
      lambSet(iter)=lamb;
for i=1:nbin
    if(count(i)<100)
        mu(i)=i;
        sigma(i)=defaultSigma;
    end
    
    for j=1:nbin
        pn(j,i)=exp(-(j-mu(i))^2/2/sigma(i)^2);
    end
end
for i=1:nbin
    pn(:,i)=pn(:,i)/sum(pn(:,i));
end
    Emit=pn;
    Emit_T=Emit';

end






figure(1)
plot(x-5,w_real)
%V=V/sum(V);
%fel_est=-log(V);

p_equ=count/sum(count);
fel_est=-log(p_equ);
fel_est=fel_est-min(fel_est);

fel_est=fel_est-min(fel_est);
hold on
plot(x,fel_est,'r');
hold on
PP=PP/sum(PP);
fel_noise=-log(PP);
fel_noise=fel_noise-min(fel_noise);

%plot(x,fel_noise);
plot(x,felnoised)
xlabel('h/nm')
ylabel('U/kT')
legend('theory','estimated','fel noised')
figure(2)
plot(-likehoodset);

figure(3)
plot(x,sigma);
polish;
save hmm_steric_sim;
