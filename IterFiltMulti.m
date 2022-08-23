function [IMF] = IterFiltMulti(sig,options)
% Multivariate Iterative Filtering (MIF)

% sig - Multichannel or multivariate signal contain Nc channel 
%       and L sample to each channel
%       sig --- Nc x L
% IMF - is a cell array, each cell contain the multivariate intrinsic mode function (MIMFs)
%       corresponding to all particular channel. For example cell (1,1) contaion first MIMF corresponding
%       to all channel.

% Please cite the following paper if are using this code or
% part of the code.
%
% [1] Das, Kritiprasanna, and Ram Bilas Pachori. "Schizophrenia 
% detection technique using multivariate iterative filtering and
% multichannel EEG signals." Biomedical Signal Processing and 
% Control 67 (2021): 102525.
% [2] Cicone, Antonio, Jingfang Liu, and Haomin Zhou. "Adaptive 
% local iterative filtering for signal decomposition and 
% instantaneous frequency analysis." Applied and Computational
% Harmonic Analysis 41.2 (2016): 384-411.
% (Ripped from IF_v8_3.m by A Cicone)
%
%
% For any queries or help plese feel free to write a mail to 
% kpdas95@gmail.com. I will be hapy to help.

%%
L=size(sig,2); % Length of the signal
Nc=size(sig,1); % Number of the channel
E_sig=sum(sig.^2,2)'; % Signal Energy
E_imf=zeros(1,Nc); % Energy of the MIMFs
SC_E=zeros(1,Nc); % Stopping Criteria depend on energy
SC_Th=.97; % Threshold to stop iteration
IMF={};

logM=zeros(1,options.IF.NIMFs);
f=sig;
N = length(f);
if size(f,1)>size(f,2)
    f = f.';
end


%IMF =zeros(options.IF.NIMFs,N);

nameFile=sprintf('%1.0d',sum(round(clock*1000)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Iterative Filtering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
load('prefixed_double_filter','MM');


%Create a signal without zero regions and compute the number of extrema
f=sig;
for i=1:Nc
    f_pp=[];
    f_pp=f(i,:);
    f_pp(abs(f_pp)<=10^-18)=[];
    if(SC_E(i)<SC_Th)
        maxmins_pp=Maxmins_v3_3(f_pp,options.IF.extensionType);
        diffMaxmins_pp=diff(maxmins_pp);
        N_pp_t(i)=length(f_pp);
        k_pp_t(i) = length(maxmins_pp);
    else
        N_pp_t(i)=length(f_pp);
        k_pp_t(i) = 0;
    end
end

[k_pp,k_ppI]=max(k_pp_t); %Finding maximum no of extream exist through all channels
N_pp=N_pp_t(k_ppI); % chosing the length of the signal corresponds to the channel having maximum no of extremas.

countIMFs=0;
% M=[20,66,200];
M=[];
%Stops depend on extrema point, Sum of energy cover by IMFs
while countIMFs < options.IF.NIMFs && k_pp>=options.IF.ExtPoints && toc<options.maxTime && (sum(SC_E<SC_Th)>=1 )
    countIMFs=countIMFs+1;
    
    SD=1;
    SD_tmp=ones(1,Nc);
    h=f; % Choose proper channel having maximum no of extrema
    
    
    if isempty(M) || length(M)<countIMFs
        
        if isa(options.IF.alpha,'char')
            if strcmp(options.IF.alpha,'ave') % Using an average mask length
                m = 2*round(N_pp/k_pp*options.IF.Xi);
            elseif strcmp(options.IF.alpha,'Almost_min') % Using an almost min mask length
                if 2*round(options.IF.Xi*prctile(diffMaxmins_pp,30))<2*round(N_pp/k_pp*options.IF.Xi)
                    m = 2*round(options.IF.Xi*prctile(diffMaxmins_pp,30));
                else
                    m = 2*round(N_pp/k_pp*options.IF.Xi);
                end
                if countIMFs>1
                    if m<=logM(countIMFs-1)
                        m=ceil(logM(countIMFs-1)*1.1);
                    end
                end
                
            else
                disp(' Value of alpha not recognized')
                return
            end
        else % using a fixed value alpha
            m = 2*round(options.IF.Xi*(max(diffMaxmins_pp)*options.IF.alpha+min(diffMaxmins_pp)*(1-options.IF.alpha)));
        end
    else
        m=M(countIMFs);
    end
    
    inStepN=0;

    logM(countIMFs)=m;
    a = get_mask_v1(MM,m);
    %rec_a(countIMFs,1:length(a))=a;
    
    ExtendSig=1==0;
    if N < length(a) % we need to extend the signal
        ExtendSig=1==1;
        Nxs=ceil(length(a)/N);
        N_old=N;
        if rem(Nxs,2)==0
            Nxs=Nxs+1;
        end
        h_n=[];
        for ii=1:Nxs
            h_n=[h_n h];
        end
        h=h_n;
        N=Nxs*N;
    end
    
    Nza=N-length(a);
    if rem(Nza,2)==0
        a = [zeros(1,Nza/2) a zeros(1,Nza/2)];
        ifftA=real(fft([a((length(a)-1)/2+1:end) a(1:(length(a)-1)/2)]));
        % figure,plot(circshift(a,(length(a)-1)/2+1)-ifft(real(fft(circshift(a,(length(a)-1)/2+1)))),'r')
    else
        a = [zeros(1,(Nza-1)/2) a zeros(1,(Nza-1)/2+1)];
        %csA=circshift(a,(length(a))/2+1);
        ifftA=real(fft([a((length(a))/2:end) a(1:(length(a))/2-1)]));
        % figure,plot(circshift(a,(length(a))/2+1)-ifft(real(fft(circshift(a,(length(a))/2+1)))),'r')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Change
    
    
    if options.plots>0 %&& rem(inStepN,5)==0
        if gcf > 30
            close all
        end
        figN=figure;
        set(figN,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    % Stops depend on standard deviation
    while sum(SD_tmp>options.IF.delta)>1 && inStepN < options.IF.MaxInner
        inStepN=inStepN+1;
        fftH=fft(h.').' ;
        h_ave=ifft((ifftA.*fftH).').';
        
        %%%%%%%%%%%%%%%% Updating stopping criterium %%%%%%%%%%%%%%%%%
        
        %SD=norm(h_ave(k_ppI,:))^2/norm(h(k_ppI,:))^2;
        
        SD_tmp=(sum(h_ave.*h_ave,2)./sum(h.*h,2))';
        
        
        %%%%%%%%%%%%%%%%%% generating f_n %%%%%%%%%%%%%%%%%%
        
        h=h-h_ave;      

    end
    
    if ExtendSig % we reduce the signal
        N=N_old;
        h=h(:,N*(Nxs-1)/2+1:N*((Nxs-1)/2+1));
    end
    if inStepN >= options.IF.MaxInner
        disp('Max # of inner steps reached')
        %return
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    E_imf=E_imf+sum(h.*h,2)';
    SC_E=E_imf./E_sig;
    
 
    IMF{1,countIMFs}= h;
    f=f-h;
    for i=1:Nc
        f_pp=[];
        f_pp=f(i,:);
        f_pp(abs(f_pp)<=10^-18)=[];
        if(SC_E(i)<SC_Th)
            maxmins_pp=Maxmins_v3_3(f_pp,options.IF.extensionType);
            if isempty(maxmins_pp)
                break
            end
            diffMaxmins_pp=diff(maxmins_pp);
            N_pp_t(i)=length(f_pp);
            k_pp_t(i) = length(maxmins_pp);
        else
            N_pp_t(i)=length(f_pp);
            k_pp_t(i) = 0;
        end
    end
    
    [k_pp,k_ppI]=max(k_pp_t); %Finding maximum no of extream exist through all achannels
    N_pp=N_pp_t(k_ppI); % chosing the length of the signal corresponds to the chaneel having maximum no of extremas.
    
    
end %End of while, outer loop
%disp('Out of while loop');

IMF{1,countIMFs+1} =  f;
logM = logM(1:countIMFs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_imf=E_imf+(sum(f.*f,2))';
SC_E=E_imf./E_sig;
%ttt=cputime-tt;
toc;

%Plotting IMF
% for i=1:size(IMF,2)
%     figure(i);
%     plot([IMF{1, i}']);
% end


%Auxiliar functions

%disp(logM)


end
function a=get_mask_v1(y,k)
%
% Rescale the mask y so that its length becomes 2*k+1.
% k could be an integer or not an integer.
% y is the area under the curve for each bar

n=length(y);
m=(n-1)/2;

if k<=m % The prefixed filter contains enough points
    
    if mod(k,1)==0     % if the mask_length is an integer
        
        a=zeros(1,2*k+1);
        
        for i=1:2*k+1
            s=(i-1)*(2*m+1)/(2*k+1)+1;
            t=i*(2*m+1)/(2*k+1);
            
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            %t2=ceil(t)-t;
            
            if floor(t)<1
                disp('Ops')
            end
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        
    else   % if the mask length is not an integer
        new_k=floor(k);
        extra = k-new_k;
        c=(2*m+1)/(2*new_k+1+2*extra);
        
        a=zeros(1,2*new_k+3);
        
        t=extra*c+1;
        t1=t-floor(t);
        %t2=ceil(t)-t;
        if k<0
            disp('Ops')
            a=[];
            return
        end
        a(1)=sum(y(1:floor(t)))+t1*y(floor(t));
        
        for i=2:2*new_k+2
            s=extra*c+(i-2)*c+1;
            t=extra*c+(i-1)*c;
            %s1=s-floor(s);
            s2=ceil(s)-s;
            
            t1=t-floor(t);
            
            
            a(i)=sum(y(ceil(s):floor(t)))+s2*y(ceil(s))+t1*y(floor(t));
        end
        t2=ceil(t)-t;
        
        a(2*new_k+3)=sum(y(ceil(t):n))+t2*y(ceil(t));
    end
else % We need a filter with more points than MM, we use interpolation
    dx=0.01;
    % we assume that MM has a dx = 0.01, if m = 6200 it correspond to a
    % filter of length 62*2 in the physical space
    f=y/dx; % function we need to interpolate
    dy=m*dx/k;
    b=interp1(0:m,f(m+1:2*m+1),0:m/k:m);
    if size(b,1)>size(b,2)
        b=b.';
    end
    if size(b,1)>1
        fprintf('\n\nError!')
        disp('The provided mask is not a vector!!')
        a=[];
        return
    end
    a=[fliplr(b(2:end)) b]*dy;
    if abs(norm(a,1)-1)>10^-14
        fprintf('\n\n Warning!\n\n')
        fprintf(' Area under the mask equals %2.20f\n',norm(a,1))
        fprintf(' it should be equal to 1\n We rescale it using its norm 1\n\n')
        a=a/norm(a,1);
    end
end

end


function varargout = Maxmins_v3_3(f,extensionType)
% Based on version 3
% Minor revisions: 1) added for constant extention the checking for Mins and
%                     Maxs emptiness
%                  2) completed the code for the periodical case

% Based on Version 2.
% Modified the way zero-derivative regions are handled.
%
% Identify the maxima and minima of a signal f

tol=10^-15;

if nargin == 1, extensionType = 'p'; end
N = length(f);
Maxs = zeros(1,N);
Mins = zeros(1,N);
df = diff(f);


h = 1;
%cIn=0;
if strcmp(extensionType,'p')% && df(1) == 0% && df(end) == 0 && f(N)-f(1)==0
    while h<N && abs(df(h)) <= tol % Condion for no maxima minima
        %        cIn=cIn+1;
        h=h+1;
    end
    %     if df(h) < 0
    %         Initial_df=-1;
    %     else
    %         Initial_df=+1;
    %     end
    if h==N
        if nargout<=1
            varargout{1}=[];  % No maxima minima found
        elseif nargout==2
            varargout{1}=[];
            varargout{2}=[];
        end
        return
    end
end

cmaxs=0;
cmins=0;

if strcmp(extensionType,'c') && abs(df(1)) <= tol
    while abs(df(h)) <= tol
        h=h+1;
    end
    if df(h) < -tol
        cmaxs=cmaxs+1;
        Maxs(cmaxs)=h;
    elseif df(h) > +tol
        cmins=cmins+1;
        Mins(cmins)=h;
    end
end

c = 0;

N_old=N;
if strcmp(extensionType,'p')
    df=diff([f f(2:h+1)]);
    N=N+h;
end

last_df=[];
for i=h:N-2
    if   df(i)*df(i+1) <= tol && df(i)*df(i+1) >= -tol
        if df(i) < -tol
            last_df=-1;
            posc = i;
        elseif df(i) > tol
            last_df=+1;
            posc = i;
        end
        c = c + 1;
        if df(i+1) < -tol
            if last_df==+1
                cmaxs=cmaxs+1;
                Maxs(cmaxs)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        if df(i+1) > tol
            if last_df==-1
                cmins=cmins+1;
                Mins(cmins)=mod(posc+floor((c-1)/2)+1,N_old);
            end
            c=0;
        end
        
    end
    if   df(i)*df(i+1) < -tol
        if df(i) < -tol && df(i+1) > tol % Condition for minima
            cmins=cmins+1;
            Mins(cmins)=mod(i+1,N_old);
            if Mins(cmins)==0
                Mins(cmins)=1;
            end
            last_df=-1;
        elseif df(i) > tol && df(i+1) < -tol % Condition for maxima
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=mod(i+1,N_old);
            if Maxs(cmaxs)==0
                Maxs(cmaxs)=1;
            end
            last_df=+1;
        end
    end
end
if c > 0
    %     if strcmp(extensionType,'p')
    %         % we deal with the boundary
    %         df_0=f(N)-f(1);
    %         if df_0==0
    %             if Initial_df < 0
    %                 if last_df==+1
    %                     cmaxs=cmaxs+1;
    %                     Maxs(cmaxs)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             elseif Initial_df > 0
    %                 if last_df==-1
    %                     cmins=cmins+1;
    %                     Mins(cmins)=mod(posc+floor((c+cIn-1)/2)+1,N);
    %                 end
    %             end
    %         else
    %             disp('Code missing!')
    %         end
    %     end
    if strcmp(extensionType,'c')
        if last_df > 0
            cmaxs=cmaxs+1;
            Maxs(cmaxs)=posc+1;
        else
            cmins=cmins+1;
            Mins(cmins)=posc+1;
        end
    end
    %%For debugging purpose
    if cmins==0 || ~islogical(cmins)
        disp('cmins is not logical');
        %cmins =1;
    end
    %%%%%%%%%%%%
    if ~cmins==0
        if Mins(cmins)==0
            Mins(cmins)=N;
        end
    end
    if ~cmaxs==0
        if Maxs(cmaxs)==0
            Maxs(cmaxs)=N;
        end
    end
end

Maxs=Maxs(1:cmaxs);
Mins=Mins(1:cmins);
maxmins=sort([Maxs Mins]);

% if strcmp(extensionType,'p') % we deal with a periodic signal
%     disp('Code to be completed')
%     if isempty(maxmins)
%         maxmins = 1;
%     else
%         if maxmins(1)~=1 && maxmins(end)~=N
%             if (f(maxmins(end)) > f(maxmins(end)+1) && f(maxmins(1)) > f(maxmins(1)-1)) || (f(maxmins(end)) < f(maxmins(end)+1) && f(maxmins(1)) < f(maxmins(1)-1))
%                 maxmins=[1 maxmins];
%             end
%         end
%     end
% else
if strcmp(extensionType,'c')
    if not(isempty(maxmins)) && not(isempty(Mins)) && not(isempty(Maxs))
        if maxmins(1) ~= 1 && maxmins(end) ~= N && df(1)~=0 && df(end)~=0
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
            maxmins = [1, maxmins, N];
        elseif maxmins(1) ~= 1 && df(1)~=0
            maxmins = [1, maxmins];
            if Maxs(1) < Mins(1)
                Mins=[1 Mins];
            else
                Maxs=[1 Maxs];
            end
        elseif  maxmins(end) ~= N && df(end)~=0
            maxmins = [maxmins, N];
            if Maxs(end) < Mins(end)
                Maxs=[Maxs N];
            else
                Mins=[Mins N];
            end
        end
    end
elseif strcmp(extensionType,'r')
    disp('Code to be completed')
    if isempty(maxmins)
        maxmins = [1, N];
    else
        if maxmins(1) ~= f(1) && maxmins(end) ~= f(end)
            maxmins = [1, maxmins, N];
        elseif f(maxmins(1)) ~= f(1)
            maxmins = [1, maxmins];
        elseif  f(maxmins(end)) ~= f(end)
            maxmins = [maxmins, N];
        end
    end
end

if nargout<=1
    varargout{1}=maxmins;
elseif nargout==2
    varargout{1}=Maxs;
    varargout{2}=Mins;
end

end
