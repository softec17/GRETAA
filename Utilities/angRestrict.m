function ang_out = angRestrict(ang_in,ang_low,ang_hi)
% Assuming that angles are in radians

ang_out = unwrap(ang_in);   % Restrict
ang_out = mod(ang_out,ang_hi);

return