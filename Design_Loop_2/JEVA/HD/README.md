### Data_Prep_For_JEVA.m

The input file used in this script is the dfs0 file which contains Water Elevation, U Velocity, V Velocity, Current Speed and Current Direction columns for different stations in the following format:

Water Elevation (station 1), Water Elevation (station 2), Water Elevation (station 3), ..., Water Elevation (station N), U Velocity (station 1), U Velocity (station 2), U Velocity (station 3), ..., U Velocity (station N), V Velocity (station 1), V Velocity (station 2), V Velocity (station 3), ..., V Velocity (station N), Current Speed (station 1), Current Speed (station 2), Current Speed (station 3), ..., Current Speed (station N), Current Direction (station 1), Current Direction (station 2), Current Direction (station 3), ..., Current Direction (station N)

In order to have consistant m_structure objects with the ones required for JEVA script, the following lines should be replaced at the end of m_tide file:

	% components as m_structures
    if isWL
        S.WL    = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,1)],'WLtot',1,X.bins);
        S.WLTde = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,2)],'WLtide',1,X.bins);
        S.WLRes = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,3)],'WLres',1,X.bins);
    else
        S.CS    = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,1)],'CStot',1,X.bins);
        S.CSTde = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,2)],'CStide',1,X.bins);
        S.CSRes = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,3)],'CSres',1,X.bins);
        S.CD    = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,4)],'CDtot',1,Y.bins);
        S.CDTde = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,5)],'CDtide',1,Y.bins);
        S.CDRes = m_structure(X.name,X.xyz,X.ttt,X.legend,[X.time,data_out(:,6)],'CDres',1,Y.bins);
    end
	
	
The change is made at the item of the parameters. The most important change is related to residual water level which is changed from "WL_Residual" to "WLres" because accoring to m_item, the unit of "WL_Residual" is "m" but the unit of "WLres" is "m+vertical reference" (e.g. "mMSL").


