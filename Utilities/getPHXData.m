function phx_data = getPHXData()

dummy = load('phx_tlm_v9.mat');
phx = dummy.phx_tlm;

time = phx.nonchan.imu.sclk;
t0 = 896225523.82790;           % Time of entry (according to Blanchard paper)
time = time - t0;

accel_body = phx.nonchan.imu.cruise_frame.filtered.accel_cg;
ang_rate_body = phx.nonchan.imu.cruise_frame.filtered.rate;

mass = phx.nonchan.imu.cruise_frame.mass_prop.mass;
cg = phx.nonchan.imu.cruise_frame.mass_prop.cg_c;
imu = cg;
imu(:,1) = 1.284866;
imu(:,2) = -0.5195119;
imu(:,3) = -0.228428;

radartime = phx.nonchan.radar.proc.time;
radartime = radartime - t0;
alt = phx.nonchan.radar.proc.altitude.m;


phx_data.IMU.time = time;
phx_data.IMU.accel = accel_body;
phx_data.IMU.gyro = ang_rate_body;
phx_data.mass = mass;
phx_data.cg_loc = cg;
phx_data.IMU.imu_loc = imu;
phx_data.radar.time = radartime(1:end-6);
phx_data.radar.alt = alt(1:end-6);

return