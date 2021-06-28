function [f]= MPC_progress(part,subpart,f,si,r)

    %% 0. Initialising
    if part==0
        if subpart==1    
            f = waitbar(0,' '); 
            titleHandle = get(findobj(f,'Type','axes'),'Title');
            set(titleHandle,'FontSize',15)
            set(findall(f),'Units', 'normalized');
            set(f,'Position', [0.2 0.5 0.6 0.1]);
            waitbar(0,f,'(0/6) Initialising Program'); 
        elseif subpart==2    
            waitbar(1/7*0.5,f,'(0/6) Initialising Program: (1/2) Defining Directories'); 
        elseif subpart==3    
            waitbar(1/7,f,'(0/6) Initialising Program: (2/2) Defining Simulation Parameters');
        end
    
    %% 1. Assessing data
    elseif part==1
        if subpart==1    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(1/6) Assessing data: (1/4) Evaluating Identification and Validation Data'); 
        elseif subpart==2    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(1/6) Assessing data: (2/4) Loading Data'); 
        elseif subpart==3    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(1/6) Generating Animations: (3/4) Loading Data');
        elseif subpart==4    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(1/6) Generating Animations: (4/4) Loading Data');
        end
    
    %% 2. Dynamic Mode Decomposition
    elseif part==2
        if subpart==1    
            waitbar(1/7*part+1/7*1/5*subpart,f,'(2/6) Dynamic Mode Decomposition: (1/5) Defining Core P arameters'); 
        elseif subpart==2    
            waitbar(1/7*part+1/7*1/5*subpart,f,'(2/6) Dynamic Mode Decomposition: (2/5) Reading & Processing Identification Data'); 
        elseif subpart==3    
            waitbar(1/7*part+1/7*1/5*subpart,f,'(2/6) Dynamic Mode Decomposition: (3/5) Reading & Processing Validation Data');
        elseif subpart==4    
            waitbar(1/7*part+1/7*1/5*subpart,f,'(2/6) Dynamic Mode Decomposition: (4/5) Remove Base Flow');
        elseif subpart==5   
            if isempty(si)
                waitbar(1/7*part+1/7*1/5*subpart,f,'(2/6) Dynamic Mode Decomposition: (5/5) Performing Singular Value Decomposition');
            else   
                waitbar(1/7*part+1/7*1/5*subpart,f,['(2/6) Dynamic Mode Decomposition: (5/5) Generating & Saving Models (' ,num2str(si) ,' / ',num2str(r), ')']);
            end
        end
    
     %% 3. Validation
    elseif part==3
        if subpart==1    
            waitbar(1/7*part+1/7*1/2*subpart,f,['(3/6) Models Validaiton: (1/2) Validating Models (' ,num2str(si) ,' / ',num2str(r), ')']); 
        elseif subpart==2    
            waitbar(1/7*part+1/7*1/2*subpart,f,'(3/6) Model Validaiton: (2/2) Saving Model Fit for Identification and Validation '); 
        end
    
     %% 4. Dynamical Analysis
    elseif part==4
        if subpart==1    
            waitbar(1/7*part+1/7*1/2*subpart,f,'(4/6) Dynamical Analysis: (1/3) Evaluating DMD Models Modal Dynamical Properties (frequency, energy, damping)'); 
        elseif subpart==2    
            waitbar(1/7*part+1/7*1/2*subpart,f,'(4/6) Dynamical Analysis: (2/3) Evaluating Proper Orthogonal Decomposition (POD) Modes '); 
        elseif subpart==3
             waitbar(1/7*part+1/7*1/2*subpart,f,['(4/6) Dynamical Analysis: (3/3) Evaluating Three Dimensional Proper Orthogonal Decomposition (POD) Modes  (' ,num2str(si) ,' / ',num2str(r), ')']); 
        end
     
     %% 5. Comparison
    elseif part==5
        if subpart==1    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(5/6) Flow Field Reconstruction: Reconsruction of Flow based on DMD Dynamical Properities '); 
        elseif subpart==2    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(5/6) Flow Field Reconstruction: Comparing DMD Flow Reconstruction and SOWFA Flow '); 
        elseif subpart==3    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(5/6) Flow Field Reconstruction: Evaluating Model Error');
        elseif subpart==4    
            waitbar(1/7*part+1/7*1/4*subpart,f,'(5/6) Flow Field Reconstruction: Saving Meaningfull Results to Directory');
        end
    end
end