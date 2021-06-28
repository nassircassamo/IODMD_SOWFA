function [ ] = printresults(f,P,LambdaDiag,damping,method,dirdmd_val)

nomedoficheiro='/results_dmd';
method=num2str(method);
extensao=strcat(method, '.txt');
nomedoficheiro=strcat(nomedoficheiro, extensao);
permissao= 'wt'; 
fid=fopen([dirdmd_val nomedoficheiro], permissao); 

   if fid == -1
        disp('Error openning file'); 
   else
   fprintf(fid,['Mode number  |  Pole location ', char(955) ,...
       ' [Real] and [Imaginary]  |  Frequency ' ,char(937) ,' [Hz]  |  '...
       'Frequency ' ,char(937) ,' [St: fD/U]  |   '...
       'Damping ' ,char(958) ,' [ ]   |   '...
       'Energy || ' ,char(958) ,' ||   |   ' ...
       'Energy fraction [ ] ']);

    for si= 1:length(f)
        fprintf(fid,'\n %10.0f %15.4f %18.4f  %20.6f %20.4f %20.4f %20.2f %20.2f', char(si),...
            real(LambdaDiag(si)),imag(LambdaDiag(si)), f(si), f(si)*178/9, damping(si), P(si),P(si)/sum(P)*100 );
    end
          
    result = fclose(fid); %fecho do ficheiro
    if result ~=0
        frintf('Erroor closing file');  
    end
        fprintf('\n The File %14s has been created and saved with success \n', nomedoficheiro);
   end
end

