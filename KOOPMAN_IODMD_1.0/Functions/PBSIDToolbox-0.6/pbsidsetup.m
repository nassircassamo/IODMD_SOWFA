% Add current directory to the current search path
addpath(pwd);
disp('PBSIDToolbox succesfully added to current search path.')

% Add current directory to the startup path
status = savepath;
if status == 1
 warning('pbsidsetup:AddToPathdefFailed',['Could not add PBSIDToolbox to the MATLAB search path defaults. Try to add the toolbox path (',pwd,') manually to pathdef.m.'])
else
 disp('PBSIDToolbox succesfully added to MATLAB search path defaults.')
end

compilation_ok = 1;
% compile functions
try
    disp('Compiling lpvitr...')
    mex('-O','@idafflpv/private/lpvitr.c');
    system('move lpvitr.* @idafflpv/private/');
catch exception
    compilation_ok = 0;
    disp('Could not compile @idafflpv/private/lpvitr.c');
    disp('You can still use the slower ".m" version.');
    disp(['Compilation error: ',exception.message]);
end

try
    disp('Compiling sfun_rttime...')
    mex('-O','simulink/sfun_rttime.c');
    system('move sfun_rttime.* simulink/');
catch exception
    compilation_ok = 0;
    disp('Could not compile simulink/sfun_rttime.c');
    disp(['Compilation error: ',exception.message]);
end

% compile functions SPGL1
try
    disp('Compiling SPGL1 oneProjector to enable BPDN regularization...')
    mex private/spgl1oneProjector/oneProjectorMex.c private/spgl1oneProjector/oneProjectorCore.c private/spgl1oneProjector/heap.c -output private/oneProjectorMex -DNDEBUG
catch exception
    compilation_ok = 0;
    disp('Could not compile oneProjector.');
    disp('You can still use the slower ".m" version.');
    disp(['Compilation error: ',exception.message]);
end
if compilation_ok
    disp('Compilation of mex files was succesful.')
    disp('PBSID setup completed');
else
    disp('PBSID setup completed, but without compilation of MEX-files, resulting in slow algorithms.');
    disp('If the compiler was missing, please install one (e.g. the free Visual Studio C++ for Windows), and rerun the PBSID setup.');
end
