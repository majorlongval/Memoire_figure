function DS = experiment_loader(t,f)
%EXPERIMENT_LOADER Takes the data in input as a table and gives back a struct of the
% the none numerical data must be removed
% all the timestamps are put to zero with reference to the first value
DS.time    = t.accelerometerTimestamp_sinceReboots-...
             t.accelerometerTimestamp_sinceReboots(1); 
DS.fs = 1/(DS.time(2)-DS.time(1));

DS.accel_x = 9.81*(t.accelerometerAccelerationXG-...
                   t.accelerometerAccelerationXG(1000));
DS.accel_xf = lowpass(DS.accel_x,f,DS.fs);

DS.accel_y = 9.81*(t.accelerometerAccelerationYG-...
                   t.accelerometerAccelerationYG(1000));
DS.accel_yf = lowpass(DS.accel_y,f,DS.fs);

DS.accel_z = 9.81*(t.accelerometerAccelerationZG-...
                   t.accelerometerAccelerationZG(1000));
DS.accel_zf = lowpass(DS.accel_z,f,DS.fs);

DS.gyro_x  = (t.gyroRotationXrads);
DS.gyro_xf = lowpass(DS.gyro_x,f,DS.fs);

DS.gyro_y  = (t.gyroRotationYrads-t.gyroRotationYrads(1000));
DS.gyro_yf = lowpass(DS.gyro_y,f,DS.fs);

DS.gyro_z  = (t.gyroRotationZrads-t.gyroRotationZrads(1000));
DS.gyro_zf = lowpass(DS.gyro_z,f,DS.fs);

DS.roll    = (t.motionRollrad-t.motionRollrad(1000))*180/pi;
DS.rollf = lowpass(DS.roll,f,DS.fs);

DS.pitch   = (t.motionPitchrad-t.motionPitchrad(1000))*180/pi;
DS.pitchf = lowpass(DS.pitch,f,DS.fs);

DS.yaw     = (t.motionYawrad-t.motionYawrad(1000))*180/pi;
DS.yawf = lowpass(DS.yaw,f,DS.fs);

end

