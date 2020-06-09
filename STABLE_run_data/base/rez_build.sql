SELECT
REZ_Build.RezID AS Region,
WindLimitOnshore,
WindLimitOffshore,
SolarLimit,
TotalBuildCap,
ExpansionLimit,
TransExpansionCost*1000 AS REZoneExpansionCost -- conversion from $/kW to $/MW
-- ConnectCost  -- given in $/kW
FROM REZ_Build
