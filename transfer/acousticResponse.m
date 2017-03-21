function acousticResponse(N,GR_DWELL,gr,sty,dirA,first,timeStart,figNr, dirsSimu)
Ndir=size(dirsSimu,1)
%INPUT:     N       :   number of sampling points
%           GR_DWELL:   dwell time of the system (temporal distance of
%           samples)
%           sty     :   style used for legends
%           dirA    :   Array used for the scan directions
%           (Read/Phase/Slice)
%           first   :   first example -> used for plotting
%           timeStart:  Start time of sequence (relevant as slice gradients
%           typically before 0)
%           figNr   :   used to plot into right figure
%           
%           TO BE SUPPLIED IN ADDITION: Scanner individual transfer
%           function (mp_spl_trf_xyz.txt)

 fs=(1/GR_DWELL) ; %sampling frequency in 1 / s 
      
    [B,A]=adsgn(fs);
    [H,W]=aspec(B,A,fs);
    
    load('mp_spl_trf_x.dat');  % supply frequency response function here
    load('mp_spl_trf_y.dat');   
    load('mp_spl_trf_z.dat');

    irfh(dirA(1),:,:)= mp_spl_trf_x;
    irfh(dirA(2),:,:)= mp_spl_trf_y;
    irfh(dirA(3),:,:)= mp_spl_trf_z;

    fs_i = irfh(1,2,1)-irfh(1,1,1);
    irf(:,:,1)=irfh(:,:,1);
    irf(:,:,2)=irfh(:,:,2)+i*irfh(:,:,3);
    
    bin_vals = [0:floor((N/2-1))];
    fax_Hz = (bin_vals./(N/fs));%/ 100 from s to s aka from Hz to Hz
   
    fs_s = fax_Hz(2)-fax_Hz(1);
    tt=[0:GR_DWELL:GR_DWELL*(N-1)]';
    tt=tt+timeStart;

    figure(figNr+1)

    for diri=1:Ndir
        dir=dirsSimu(diri);
        ft_gr(dir,:)=fft(gr(:,dir),N)./N; %complex
        up_ft_gr(dir,:)=interp1((fax_Hz),((ft_gr(dir,1:length(fax_Hz)))),[0:2:irf(1,end,1)],'linear');
        spl(dir,:)=(irf(dir,:,2).*(up_ft_gr(dir,:)./GR_DWELL));
    end


    for diri=1:Ndir
        dir=dirsSimu(diri);

        subplot(2,Ndir,diri) 
        plot(tt*1000,gr(:,dir),'color',sty{dir},'linewidth',2);
        xlabel('time in ms')
        ylabel('gradient strenght in T/m');
        title('gradient waveforms');
        hold on;
        
        subplot(2,Ndir,diri+Ndir)
        plot(fax_Hz, (abs((ft_gr(dir,1:length(fax_Hz))))),'color',sty{dir},'linewidth',2);
      %  axis([0, 8000, 0, max(max(abs(ft_gr(dir,1:length(fax_Hz)))),0.0001)]);

        xlabel('Frequency (Hz)')
        ylabel('Magnitude');
        title('Mag. spectrum (Hz)');
        hold on;
    end

    figure(figNr+2)
    legendx={'IRF_y','A-weighting','g_y'};
    legendy={'IRF_x','A-weighting','g_x'};
    legendz={'IRF_z','A-weighting','g_z'};
    
    for diri=1:Ndir
        dir=dirsSimu(diri);
        subplot(2,Ndir,diri) 
        axis([0, 8000, 0, max(max(abs(irf(dir,:,2))))]);     
        mimi=min(min(min(abs(ft_gr(dir,1:length(fax_Hz)))),min(min(abs(irf(dir,:,2))))));
        mama=max(max(max(abs(ft_gr(dir,1:length(fax_Hz)))),max(max(abs(irf(dir,:,2))))));  
        
        if first==0
            plot(irf(dir,:,1),abs(irf(dir,:,2)),'color',[120/256,120/256,120/256],'linewidth',1);
            %hold on;
            %plot(W/2/pi, 20*log10(abs(H)).*3,'color',[130/256,130/256,80/256],'linewidth',2);
        end
        hold on;
       plot(fax_Hz, abs(ft_gr(dir,1:length(fax_Hz)))./(max(max(abs(ft_gr(dir,1:length(fax_Hz))))))*10,'color',sty{dir},'linewidth',2);

        
        xlabel('Frequency (Hz)')
        ylabel('Magnitude');
   
        xlabel('Frequency (Hz)')
        ylabel('Magnitude');
        title('Transfer functions and gradient spectra');
        hold on
        
        subplot(2,Ndir,diri+Ndir)
        plot(irf(1,:,1),(abs(spl(dir,:))),'color',sty{dir},'linewidth',2);
        hold on;
        
        xlabel('Frequency (Hz)')
        ylabel('Magnitude');
        title('IRF *g');
    end
end