function sensorData = sensorModel(s,pos,m)
% to genreated simulated sensor data based on the source term, the sensor
% position and sensor characteristics.

conc = plumeModel(s,pos); %根据羽流模型计算污染物分布

% add noise 
datasize = size(conc);
error = m.sig_pct * conc .* randn(datasize); % add noise or fluctuations 
% 传感器的误差设置。
sensorData = conc + error;

% 传感器未响应的状态
% not detection if below the threshold
sensorData(sensorData<m.thresh) = 0;

% not detectin due to the mis-detection rate
sensorData(rand(datasize)<(1-m.Pd)) = 0;

end

