% Add current map to the path
addpath(pwd);
disp('LTIToolbox succesfully added to the path.')

% Compilation of mex-files
MATLAB_PATH = matlabroot;
COMPILE_OPTIONS = '';
v = ver('matlab');
matver = sscanf(v.Version, '%d.%d.%d')';
COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -DMATLAB_VERSION=0x' sprintf('%02d%02d', matver(1), matver(2)) ];
MATLAB_VERSION = matver(2);

if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer) ...
        || strcmpi('MACI', computer) || strcmpi('MAC', computer) ...
        || strcmpi('MACI64', computer)
    % GNU/Linux (x86-32 or x86-64) or MacOS (Intel or PPC)
    LAPACK_PATH = ' -lmwlapack';
    if MATLAB_VERSION < 5;
        BLAS_PATH = '';
    else
        BLAS_PATH = ' -lmwblas';
    end
    if strcmpi('GLNX86', computer) || strcmpi('GLNXA64', computer)
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs',' -DNON_UNIX_STDIO'];
    else
        COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DSkip_f2c_Undefs'];
    end
elseif strcmpi('PCWIN', computer) || strcmpi('PCWIN64', computer)
    % Windows (x86-32 or x86-64)
    if strcmpi('PCWIN', computer)
        if MATLAB_VERSION < 6;
            MANUFACTURER = 'lcc';
        else
            cc = mex.getCompilerConfigurations('Any','Selected');
            MANUFACTURER = cc.Manufacturer;
        end
        switch lower(MANUFACTURER)
            case {'lcc'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'lcc', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DLCCWIN32'];
            case {'microsoft'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'microsoft', 'libmwlapack.lib');
            case {'sybase'}
                BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwblas.lib');
                LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win32', 'watcom', 'libmwlapack.lib');
                COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DWATCOMWIN32'];
            otherwise
                disp('Try "mex -setup", because BLAS/LAPACK library is not available!')
        end
    else
        BLAS_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
        LAPACK_PATH = fullfile(MATLAB_PATH, 'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
    end
    if MATLAB_VERSION < 5;
        BLAS_PATH = ''; % On <= 7.4, BLAS in included in LAPACK
    else
        BLAS_PATH = [' "' BLAS_PATH '"'];
    end
    LAPACK_PATH = [' "' LAPACK_PATH '"'];
    COMPILE_OPTIONS = [COMPILE_OPTIONS,' -DMSDOS',' -DUSE_CLOCK',' -DNO_ONEXIT'];
else
    error('Unsupported platform')
end

% Large array dims for 64 bits platforms appeared in Matlab 7.3
if (strcmpi('GLNXA64', computer) || strcmpi('PCWIN64', computer) ...
        || strcmpi('MACI64', computer)) ...
        && ~(MATLAB_VERSION < 3)
    COMPILE_OPTIONS = [ COMPILE_OPTIONS ' -largeArrayDims' ];
end

% Comment next line to suppress optimization
COMPILE_OPTIONS = [ ' -O' COMPILE_OPTIONS ];

% Comment next line to suppress compilation debugging info
%COMPILE_OPTIONS = [ ' -v' COMPILE_OPTIONS ];

disp('Compiling ltiitr...')
eval(['mex ',COMPILE_OPTIONS,' src/ltiitr.c']);
disp('Compiling fcordom...')
eval(['mex ', COMPILE_OPTIONS, ' src/fcordom.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling fdordom...')
eval(['mex ', COMPILE_OPTIONS, ' src/fdordom.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling simlns...')
eval(['mex ', COMPILE_OPTIONS, ' src/simlns.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling dac2bd...')
eval(['mex ', COMPILE_OPTIONS, ' src/dac2bdc.c',' src/ib01cd.c',...
    ' src/ib01qd.c', ' src/ib01rd.c', ' src/tb01wd.c', ' src/mb01td.c',...
    ' src/mb01sd.c', ' src/ma02ad.c', ' src/mb02ud.c', ' src/mb03ud.c',...
    ' src/mb04od.c', ' src/mb04oy.c', ' src/select.c', ' src/pow_dd.c',...
    ' src/pow_ii.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling dinit...')
eval(['mex ', COMPILE_OPTIONS, ' src/dinit.c', ' src/ib01cd.c', ' src/ib01qd.c',...
    ' src/ib01rd.c', ' src/tb01wd.c', ' src/mb01td.c', ' src/mb01sd.c',...
    ' src/ma02ad.c', ' src/mb02ud.c', ' src/mb03ud.c',' src/mb04od.c',...
    ' src/mb04oy.c', ' src/select.c', ' src/pow_dd.c', ' src/pow_ii.c',...
    BLAS_PATH, LAPACK_PATH]);
disp('Compiling dmodpo...')
eval(['mex ', COMPILE_OPTIONS, ' src/dmodpo.c', ' src/ib01bd.c', ' src/ib01pd.c',' src/ib01px.c',' src/ib01py.c',...
    ' src/sb02nd.c', ' src/sb02rd.c', ' src/sb02mt.c', ' src/sb02mr.c', ' src/sb02ms.c', ' src/sb02mv.c',...
    ' src/sb02mw.c', ' src/sb02qd.c', ' src/sb02ru.c', ' src/sb02sd.c', ' src/sb03qx.c', ' src/sb03qy.c',...
    ' src/sb03mx.c', ' src/sb03my.c', ' src/sb03sx.c', ' src/sb03sy.c', ' src/sb03mv.c', ' src/sb03mw.c',...
    ' src/sb04px.c', ' src/ma02ad.c', ' src/ma02ed.c', ' src/mb01sd.c', ' src/mb01ru.c', ' src/mb01ud.c',...
    ' src/mb01rx.c', ' src/mb01ry.c', ' src/mb01vd.c', ' src/mb02pd.c', ' src/mb02qy.c', ' src/mb02ud.c',...
    ' src/mb03od.c', ' src/mb03ud.c', ' src/mb04kd.c', ' src/mb04od.c','  src/mb04oy.c', ' src/pow_dd.c',...
    ' src/select.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling dordpo...')
eval(['mex ', COMPILE_OPTIONS, ' src/dordpo.c', ' src/ib01nd.c',' src/ib01md.c', ' src/ib01my.c', ' src/ma02ad.c',...
    ' src/ma02ed.c', ' src/ma02fd.c', ' src/mb03od.c', ' src/mb03ud.c', ' src/mb04id.c' ,' src/mb04iy.c',...
    ' src/mb04oy.c', ' src/mb04od.c', ' src/pow_dd.c', ' src/d_sign.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling fastr...')
eval(['mex ', COMPILE_OPTIONS, ' src/fastr.c', ' src/ib01md.c', ' src/ib01my.c', ' src/ma02ed.c', ' src/ma02fd.c',...
    ' src/mb04id.c', ' src/mb04od.c', ' src/mb04oy.c', ' src/d_sign.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling fastsvd...')
eval(['mex ', COMPILE_OPTIONS, ' src/fastsvd.c', ' src/mb03ud.c', BLAS_PATH, LAPACK_PATH]);
disp('Compiling ltifrf...')
if MATLAB_VERSION <= 5;
    eval(['mex ', COMPILE_OPTIONS, ' src/ltifrf.c', ' src/tb05ad.c', ' src/mb02rz.c',...
        ' src/mb02sz.c', ' src/mb02tz.c', ' src/err.c', ' src/endfile.c', ' src/close.c',...
        ' src/cabs.c', ' src/sig_die.c', ' src/d_imag.c', ' src/d_cnjg.c',...
        ' src/dcabs1.c', ' src/z_abs.c', ' src/z_div.c', BLAS_PATH, LAPACK_PATH]);
else
    eval(['mex ', COMPILE_OPTIONS, ' src/ltifrf.c', ' src/tb05ad.c', ' src/mb02rz.c',...
        ' src/mb02sz.c', ' src/mb02tz.c', ' src/err.c', ' src/endfile.c', ' src/close.c',...
        ' src/cabs.c', ' src/sig_die.c', ' src/d_imag.c', ' src/d_cnjg.c',...
        ' src/z_abs.c', ' src/z_div.c', BLAS_PATH, LAPACK_PATH]);
end
disp('Compilation of mex files was succesfull.')
disp('Run ltidemo to check correctness of compilation.')