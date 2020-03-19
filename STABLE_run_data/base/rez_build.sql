SELECT
REZ_Build.RezID AS Region,
WindLimitOnshore,
WindLimitOffshore,
SolarLimit,
TotalBuildCap,
ExpansionLimit,
TransExpansionCostMW AS TransExpansionCost
-- ConnectCostMW AS ConnectCost
FROM REZ_Build
