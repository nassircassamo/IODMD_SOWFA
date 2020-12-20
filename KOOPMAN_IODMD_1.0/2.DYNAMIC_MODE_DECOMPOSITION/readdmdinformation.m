function [rotSpeed, nacelleYaw, time1,rotorAzimuth,pitch,powerGenerator]=readdmdinformation(dirName)

for j=1:1:length(dirName)
    cd(dirName{j})  
    cd ..
    
    [nTurbine,time1,dt,nVal,powerRotor] = readTurbineOutputGlobal(dirName{j},'rotorPower');
    [time2,pitch] = readPitchData(strcat(dirName{j},'/1000/bladePitch')); %pitch=pitch{:};
    [nTurbine,time3,dt,nVal,powerGenerator] = readTurbineOutputGlobal(dirName{j},'generatorPower');
    [nTurbine,time4,dt,nVal,thrust] = readTurbineOutputGlobal(dirName{j},'rotorAxialForce');
    [nTurbine,time4,dt,nVal,thrustv] = readTurbineOutputGlobal(dirName{j},'rotorVerticalForce');
    [nTurbine,time4,dt,nVal,thrusth] = readTurbineOutputGlobal(dirName{j},'rotorHorizontalForce');
    [nTurbine,time5,dt,nVal,torqueGen] = readTurbineOutputGlobal(dirName{j},'generatorTorque');
    [nTurbine,time6,dt,nVal,torqueRotor] = readTurbineOutputGlobal(dirName{j},'rotorTorque');
    [nTurbine,time7,dt,nVal,rotSpeed] = readTurbineOutputGlobal(dirName{j},'rotorSpeed');
    [nTurbine,time7,dt,nVal,rotSpeedF] = readTurbineOutputGlobal(dirName{j},'rotorSpeedFiltered');
    [nTurbine,time8,dt,nVal,nacYaw] = readTurbineOutputGlobal(dirName{j},'nacelleYaw');
    [nTurbine,time8,dt,nVal,rotorAzimuth] = readTurbineOutputGlobal(dirName{j},'rotorAzimuth');
    %[nTurbine2,bla,time9,dt9,nVal9,bladeV] = readTurbineOutputSectional(dirName{j},'bladePointVmag');
    [nTurbine,timeyaw,dt,nVal,nacelleYaw] = readTurbineOutputGlobal(dirName{j}, 'nacelleYaw');
    %% Read all data concerning the validation from second experiment
    %Data from experiment has exact same structure

%     [nTurbine,time1_val,dt,nVal,powerRotor_val] = readTurbineOutputGlobal(dirName_val{j},'rotorPower');
%     [time2_val,pitch_val] = readPitchData(strcat(dirName_val{j},'/1000/bladePitch')); %pitch=pitch{:};
%     [nTurbine,time3_val,dt,nVal,powerGenerator_val] = readTurbineOutputGlobal(dirName_val{j},'generatorPower');
%     [nTurbine,time4_val,dt,nVal,thrust_val] = readTurbineOutputGlobal(dirName_val{j},'rotorAxialForce');
%     [nTurbine,time4_val,dt,nVal,thrustv_val] = readTurbineOutputGlobal(dirName_val{j},'rotorVerticalForce');
%     [nTurbine,time4_val,dt,nVal,thrusth_val] = readTurbineOutputGlobal(dirName_val{j},'rotorHorizontalForce');
%     [nTurbine,time5_val,dt,nVal,torqueGen_val] = readTurbineOutputGlobal(dirName_val{j},'generatorTorque');
%     [nTurbine,time6_val,dt,nVal,torqueRotor_val] = readTurbineOutputGlobal(dirName_val{j},'generatorTorque');
%     [nTurbine,time7_val,dt,nVal,rotSpeed_val] = readTurbineOutputGlobal(dirName_val{j},'rotorSpeed');
%     [nTurbine,time7_val,dt,nVal,rotSpeedF_val] = readTurbineOutputGlobal(dirName_val{j},'rotorSpeedFiltered');
%     [nTurbine,time8_val,dt,nVal,nacYaw_val] = readTurbineOutputGlobal(dirName_val{j},'nacelleYaw');
%     [nTurbine,time8_val,dt,nVal,rotorAzimuth_val] = readTurbineOutputGlobal(dirName_val{j},'rotorAzimuth');
%     [nTurbine2,bla,time9_val,dt9,nVal9,bladeV_val] = readTurbineOutputSectional(dirName_val{j},'bladePointVmag');
end
