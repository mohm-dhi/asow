function u_wind_z = u_wind_correction(u_wind, z, ref)

for zi=1:length(z)

    if z(zi)>= -1*ref
        u_wind_z(zi,:) = u_wind .* ((ref+z(zi))/ref);
    else
        u_wind_z(zi,:) = u_wind * 0;
    end

end

end