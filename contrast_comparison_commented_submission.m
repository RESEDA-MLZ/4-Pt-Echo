%the paper has been submitted to Measurement Science and Technology
%and is currently published on arxiv.org https://arxiv.org/abs/2103.06623

clear all                                                                   % clears all variables

tic                                                                         % starts timer
I_mean = 1;                                                                 % makes sure that I0 * sin(x) + I_mean > 0
trys = 500;                                                                 % sets the number of repetitions of the experiment 
NoEvents = round(logspace(1,7,180));                                        % list of how many events are tested
Contrasts = 0.05:0.05:0.95;                                                 % list of all evaluted contrasts (must be between of 0 and 1)
phi_0 = 60;                                                                 % sets the initial phase in degrees (°)

events = length(NoEvents);                                                  % defines empty variables for later use
phi_rec          = zeros(events,trys,length(Contrasts));
phi_fit_4        = zeros(events,trys,length(Contrasts));
phi_fit_16       = zeros(events,trys,length(Contrasts));
contrast_rec     = zeros(events,trys,length(Contrasts));
contrast_rec_alt = zeros(events,trys,length(Contrasts));
contrast_fit_4   = zeros(events,trys,length(Contrasts));
contrast_fit_16  = zeros(events,trys,length(Contrasts));



for m = 1:length(Contrasts)                                                 %loop for contrasts defined in line 8
    
    I_0 = I_mean*Contrasts(m);                                              % calculates the amplitude
    
    prob_I_4 = zeros(4,1);                                                  % defines empty variables for later use
    prob_I_16 = zeros(16,1);
    
    for n = 2:17
        prob_I_16(n) = I_mean/16 - I_0/(2*pi) * (cosd(22.5*(n)+phi_0) - cosd(22.5*(n-1)+phi_0));    %calculates the probabity for a single event to be in a certain time bin for 16 time bins (Eq 4 in the paper)
    end    
    
    
    for n = 2:5
        prob_I_4(n) = I_mean/4 - I_0/(2*pi) * (cosd(90*(n)+phi_0) - cosd(90*(n-1)+phi_0));  % calculates the probabity for a single event to be in a certain time bin for 4 time bins
    end    
    
    cumprob_4 = cumsum(prob_I_4);                                          % cumulative sum for a single event for 16 time bins
    cumprob_16 = cumsum(prob_I_16);                                        % cumulative sum for a single event for 4 time bins 
    
       
    for l = 1:trys                                                         % loop for list of applied runs
        
        for k = 1:events                                                   % loop for applied events (see list in line 7)

            rng('shuffle')                                                     % optionally shuffles the random generator
            emnit = I_mean.*rand(NoEvents(k),1);                               % normalize the probability

            if k == 1                                                          % genertes empty variables in the first iteration

                I_4 = zeros(4,1);
                I_16 = zeros(16,1);

            end

            I_4_tmp  = histc(emnit,cumprob_4);                                 % histogramms the random generated values into a temporal variable
            I_16_tmp = histc(emnit,cumprob_16);

            I_16 = I_16 + I_16_tmp(1:16);                                      % sums the histogramms for each loop run
            I_4 = I_4 + I_4_tmp(1:4);

            [FittedCurve_4, contrast_fit_4(k,l,m), phi_fit_4(k,l,m), Average] = Fit_Sinus(1:size(prob_I_4)-1,I_4');     %evaluates contrast phi by fitting a sine to the histogrammed events
            [FittedCurve_16, contrast_fit_16(k,l,m), phi_fit_16(k,l,m), Average] = Fit_Sinus(1:size(prob_I_16)-1,I_16');

            phi_fit_4(k,l,m)  = phi_fit_4(k,l,m);                              % writes the deduced parapmeters for each run into this variable
            phi_fit_16(k,l,m) = phi_fit_16(k,l,m);
            phi_fit_4(k,l,m)  = phi_fit_4(k,l,m)*180/pi;
            phi_fit_16(k,l,m)  = phi_fit_16(k,l,m)*180/pi;

            if contrast_fit_4(k,l,m) <= 0                                      % catches and corrects cases with negativ contrasts for 4 time bins

                contrast_fit_4(k,l,m) = abs(contrast_fit_4(k,l,m));
                phi_fit_4(k,l,m) = phi_fit_4(k,l,m) - 180;

            end

            if contrast_fit_16(k,l,m) <= 0                                     % catches and corrects cases with negativ contrasts for 16 time bins 

                contrast_fit_16(k,l,m) = abs(contrast_fit_16(k,l,m));
                phi_fit_16(k,l,m) = phi_fit_16(k,l,m) - 180;

            end

            I1 = I_4(1);                                                       % rewrites the histgrammed events for convenience into new variables (names equal to eq 5a and 5d in the paper)
            I2 = I_4(2);
            I3 = I_4(3);
            I4 = I_4(4);

            a = I1 + I4;                                                       % equations 6a to 6d in the paper
            b = I2 + I3;
            c = I1 + I2;
            d = I3 + I4;

            phi_rec(k,l,m) = 180/pi*atan2((a-b),(c-d));                        % equation 7 in the paper 

            contrast_rec(k,l,m) = (c-d)./(c+d).*pi/(2.*cosd(phi_rec(k,l,m)));  % equations 9a and 9b
            contrast_rec_alt(k,l,m) = (a-b)./(a+b).*pi/(2.*sind(phi_rec(k,l,m)));        


            %==============displays reslts during the simulation =============
            %{
            if and(l == 1, k == 1)

                f1 = figure;

            end

            figure(f1)
            plot(RoughX,I,'o')
            hold on
            plot(FineX, sum(I)/16.*(1 + contrast_rec(k,l) .* sind(FineXDeg + phi_rec(k,l))),'b');
            plot(FineX, sum(I)/16.*(1 + contrast_rec_alt(k,l) .* sind(FineXDeg + phi_rec(k,l))),':r');
            plot(FineX, sum(I)/16.*(1 + contrast_fit(k,l) .* sind(FineXDeg + phi_fit(k,l))),'k');        
            title(sprintf('run %d out of %d, \\phi_0 = %d°, c = %0.2g, No. of Events: %d out of %d',l,trys,phi_0,I_0/I_mean,NoEvents(k),NoEvents(end)))
            hold off

            %}
        end
        
        if l == 1
            
            f2 = figure;
        
        end
        
        
        phi_fit_4(phi_fit_4 > 360) = phi_fit_4(phi_fit_4 > 360) - 360;      % shifts the phase case sensitive by 2*pi if phi > 2pi; equation 10 in the paper
        phi_fit_16(phi_fit_16 > 360) = phi_fit_16(phi_fit_16 > 360) - 360;
        
        phi_rec(phi_rec <= -180) = phi_rec(phi_rec <= -180) + 360;          % shifts the phase case sensitive by 2*pi if phi < -pi; equation 10 in the paper
        phi_fit_4(phi_fit_4 <= -180) = phi_fit_4(phi_fit_4 <= -180) + 360;
        phi_fit_16(phi_fit_16 <= -180) = phi_fit_16(phi_fit_16 <= -180) + 360;
        
        figure(f2);
        set(f2,'Name',sprintf('run %d out of %d, \\phi_0 = %d°, c = %0.2g',l,trys,phi_0,I_0/I_mean))
        subplot(2,2,1);
        semilogx(NoEvents,phi_rec(:,l,m) - 360/4,'b', NoEvents,phi_fit_4(:,l,m) - 360/8,'r',NoEvents,phi_fit_16(:,l,m) - 360/32,'g')
        %semilogx(NoEvents,phi_rec(:,l,m) - 360/16,'b', NoEvents,phi_fit_4(:,l,m) - 360/4,'r',NoEvents,phi_fit_16(:,l,m) - 360/16,'g')    
        ylabel(sprintf('\\phi (°)'))
        axis([0 NoEvents(end) phi_0-5 phi_0+5])
        
        
        subplot(2,2,3); 													 % plots intermediate results during simulation
        hold off
        loglog(NoEvents,2*std(phi_rec(:,1:l,m),0,2),'b')
        hold on
        loglog(NoEvents,2*std(phi_fit_4(:,1:l,m),0,2),'r')
        loglog(NoEvents,2*std(phi_fit_16(:,1:l,m),0,2),'g')
        xlabel(sprintf('#Events'))
        ylabel(sprintf('2\\sigma_{\\phi} (°)'))
        
        
        subplot(2,2,2); 													 % plots intermediate results during simulation
        semilogx(NoEvents,contrast_rec(:,l,m),'b',NoEvents,contrast_fit_4(:,l,m),'r',NoEvents,contrast_fit_16(:,l,m),'g')    
        ylabel(sprintf('C'))
        modify_errorbars(0)
        axis([0 NoEvents(end) I_0/I_mean-0.2 I_0/I_mean+0.2])
        
        subplot(2,2,4); 													 % plots intermediate results during simulation
        hold off
        loglog(NoEvents,2*std(contrast_rec(:,1:l,m),0,2),'b')
        hold on
        loglog(NoEvents,2*std(contrast_fit_4(:,1:l,m),0,2),'r')
        loglog(NoEvents,2*std(contrast_fit_16(:,1:l,m),0,2),'g')
        %loglog(NoEvents,2*std(contrast_rec_alt(:,1:l,m),0,2),'k')
        xlabel(sprintf('#Events'))
        ylabel(sprintf('2\\sigma_{C}'))
        axis([0 NoEvents(end) 2E-4 2])
    
    end
    
end

if cosd(mean(mean(mean(phi_rec)))) >= sind(mean(mean(mean(phi_rec))))       % case sensitive handling of C_rec (eq 10)
    
    contrast_rec_merged = contrast_rec;
    contrast_rec_merged(phi_rec == -90) = contrast_rec_alt(phi_rec == -90); % catches and corrects cases with negativ contrasts
    contrast_rec_merged(phi_rec == 90) = contrast_rec_alt(phi_rec == 90);
    
else
    
   contrast_rec_merged = contrast_rec_alt;                                  % catches and corrects cases with negativ contrasts
   contrast_rec_merged(phi_rec == 0) = contrast_rec(phi_rec == 0);
   contrast_rec_merged(phi_rec == 180) = contrast_rec(phi_rec == 180);
    
end  


%calculates the standard deviations for varying events(line8), varying contrasts (line 9, fixed phase (line10) and writes them into variables

Z_rec_merged = reshape(2*std(contrast_rec_merged(:,1:l,:),0,2),length(NoEvents),length(Contrasts));      
Z_fit_4 = reshape(2*std(contrast_fit_4(:,1:l,:),0,2),length(NoEvents),length(Contrasts));
Z_fit_16 = reshape(2*std(contrast_fit_16(:,1:l,:),0,2),length(NoEvents),length(Contrasts));
Z_phi_rec = reshape(2*std(phi_rec(:,1:l,:),0,2),length(NoEvents),length(Contrasts));
Z_phi_fit_4 = reshape(2*std(phi_fit_4(:,1:l,:),0,2),length(NoEvents),length(Contrasts));
Z_phi_fit_16 = reshape(2*std(phi_fit_16(:,1:l,:),0,2),length(NoEvents),length(Contrasts));

[X,Y] = meshgrid(Contrasts,NoEvents);

toc