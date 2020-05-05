---- AEMO_Registration ------
-----------------------------

DROP TABLE IF EXISTS AEMO_Reg_GSL;
CREATE TABLE AEMO_Reg_GSL (
    Participant              TEXT,
    StationName              TEXT,
    Region                   TEXT,
    DispatchType             TEXT,
    Category                 TEXT,
    Classification           TEXT,
    FuelSourcePrimary        TEXT,
    FuelTypeMajor            TEXT,
    FuelTypeMinor            TEXT,
    FuelSourceDescriptor     TEXT,
    TechnologyTypePrimary    TEXT,
    TechnologyTypeSecondary  TEXT,
    TechnologyTypeDescriptor TEXT,
    UnitCount                INTEGER,
    Aggregation              CHAR,
    DUID                     TEXT    UNIQUE,
    RegCapMW                 REAL,
    MaxCapMW                 REAL,
    MaxROCperMin             INTEGER,
    CHECK (DispatchType IN ('Generator', 'Load') ),
    CHECK (Category IN ('Market', 'Non-Market') ),
    CHECK (Classification IN ('Scheduled', 'Non-Scheduled', 'Semi-Scheduled') ),
    CHECK (TechnologyTypePrimary IN ('Load', 'Combustion', 'Renewable', 'Storage', 'Transmission') ) 
);

DROP TABLE IF EXISTS Technology_Types;
CREATE TABLE Technology_Types AS
 SELECT * FROM STABLE_data.Technology_Types
;

DROP TABLE IF EXISTS Tech_TxZone_Cap;
CREATE TABLE Tech_TxZone_Cap AS
--select sum(Capacity) from (
SELECT 
 Technology_Types.TechTypeID, Technology_Types.TechTypeDescriptor, Technology_Types.FuelType,
 AEMO_Reg_GSL_Coords.Region AS DemandRegion, AEMO_Reg_GSL_Coords.ZoneID AS RegionID,
--RegCapMW
sum(RegCapMW) AS Capacity
FROM
 Technology_Types JOIN AEMO_Reg_GSL_Coords
  ON 
  AEMO_Reg_GSL_Coords.TechnologyTypeSecondary = Technology_Types.TechTypeDescriptor
--  AEMO_Reg_GSL_Coords.TechnologyTypeDescriptor = STABLE_data.Technology_Types.TechTypeDescriptor
AND
 AEMO_Reg_GSL_Coords.FuelTypeMinor = Technology_Types.FuelType
WHERE
AEMO_Reg_GSL_Coords.DispatchType='Generator'
AND
AEMO_Reg_GSL_Coords.TechnologyTypePrimary != 'Renewable'
AND
AEMO_Reg_GSL_Coords.ZoneID IS NOT NULL
GROUP BY TechTypeID, Region, ZoneID
--AND AEMO_Reg_GSL.Category='Market'
--GROUP BY Technology_Types.TechTypeDescriptor, Technology_Types.FuelType, AEMO_Reg_GSL.Region
--ORDER BY FuelType
--ORDER BY Region, Capacity DESC
--)
;

DROP TABLE IF EXISTS Tech_REZone_Cap;
CREATE TABLE Tech_REZone_Cap AS
--SELECT sum(Capacity) FROM (
SELECT 
Technology_Types.TechTypeID, Technology_Types.TechTypeDescriptor, Technology_Types.FuelType,
Region AS DemandRegion, REZid AS RegionID,
--DUID
sum(RegCapMW) AS Capacity
FROM
Technology_Types JOIN (
SELECT
DUID, TechnologyTypeDescriptor, TechnologyTypeSecondary, FuelTypeMinor, DispatchType, TechnologyTypePrimary, Region, Renew_REZ.REZid, RegCapMW
FROM
AEMO_Reg_GSL JOIN Renew_REZ USING (DUID)
) AS GSLtemp
ON  Technology_Types.TechTypeDescriptor = GSLtemp.TechnologyTypeSecondary
--ON  Technology_Types.TechTypeDescriptor = TechnologyTypeDescriptor
AND Technology_Types.FuelType = GSLtemp.FuelTypeMinor
GROUP BY TechTypeID, REZid
--)
;

DROP VIEW IF EXISTS Tech_All_Cap;
CREATE VIEW IF NOT EXISTS Tech_All_Cap AS
SELECT * FROM Tech_REZone_Cap
UNION
SELECT * FROM Tech_TxZone_Cap
;

-- Inspection:

SELECT DISTINCT TechTypeID FROM Tech_TxZone_Cap
;


SELECT sum(Capacity) FROM Tech_All_Cap
;


-----------------------------
----- ISP_2019_Input --------

DROP TABLE IF EXISTS REZ_Build;
CREATE TABLE REZ_Build (
    RezID                TEXT,
    WindLimitOnshore     REAL,
    WindLimitOffshore    REAL,
    SolarLimit           REAL,
    TotalBuildCap        REAL,
    ExpansionLimit       REAL,
    TransExpansionCost   REAL,
    ConnectCost          REAL
);

DROP TABLE IF EXISTS Tech_Parameters;
CREATE TABLE Tech_Parameters (
    StID              TEXT    PRIMARY KEY,
    Station           TEXT,
    TechType          TEXT,
    Status            TEXT,
    FuelType          TEXT,
    Efficiency        REAL,
    AuxiliaryLoad     REAL,
    FixedOPEX         REAL,
    VariableOPEX      REAL,
    CapCostPower      REAL,
    CapCostEnergy     REAL,
    EconomicLife      INTEGER,
    TechnicalLife     INTEGER,
    FullOutage        REAL,
    AvgPlanMainten    INTEGER,
    ConnectCost       REAL,
    FeasibilityPeriod REAL,
    ConstructionTime  REAL,
    TotalLeadTime     REAL,
    LoadFollowCost    REAL,
    StartLevel        REAL,
    MaxEtoP
);

DROP VIEW IF EXISTS STABLE_Format;
--CREATE VIEW STABLE_Format AS
DROP TABLE IF EXISTS STABLE_Format;
CREATE TABLE STABLE_Format AS
SELECT 
StID,
Efficiency,
round(FixedOPEX*1000,0) AS FixedCost, -- conversion from $/kW/yr to $/MW/yr
VariableOPEX AS VariableCost,
CapCostPower*1000 AS OvernightCostPower, -- conversion from $/kW to $/MW
round((CapCostPower/MaxEtoP)*1000,2) AS OvernightCostEnergy, -- conversion from $/kW to $/MW
0 as CurtailmentCost,
LoadFollowCost AS LoadIncreaseCost,
LoadFollowCost AS LoadDecreaseCost,
TechnicalLife AS InvestLifetime,
EconomicLife AS InvestRecovery,
30000 AS MaxCapacity,
NULL AS MaxEnergy,
StartLevel,
MaxEtoP AS MaxEnergyToPowerRatio
FROM Tech_Parameters
;

SELECT * FROM STABLE_Format;

SELECT DISTINCT REZid from Traces WHERE WindTrDUID = '' AND  SolarTrDUID = '';

SELECT DISTINCT REZid, WindTrDUID FROM Traces WHERE WindTrDUID != '' ;

SELECT DISTINCT REZid, SolarTrDUID FROM Traces WHERE SolarTrDUID != ''  GROUP BY REZid ;


--- Construct Fuel Price Data by Region

DROP TABLE IF EXISTS Fuel_BioDist2030;
CREATE TABLE Fuel_BioDist2030 AS
SELECT
TechType,
FuelType,
RegionID,
FYE2030*3.6 AS FuelCost -- Note: Conversion from $/GJ to $/MWh, and hard-coded year of 2030 etc.
FROM
(
SELECT
*
FROM Fuel_Prices
WHERE
Region='NEM'
AND CostScenario='Neutral'
)
JOIN 
(
SELECT RegionID, SuperRegion
FROM
STABLE_data.Regions
)
ON (SuperRegion = Region)
;

DROP TABLE IF EXISTS Map_ISP_Names;
CREATE TABLE Map_ISP_Names AS
SELECT StID, ISPname, ISP_Names.StationID, ZoneID, Map_Generators.FuelSourceDescriptor
FROM 
ISP_Names
JOIN
AEMO_Registration.Map_Generators
ON
ISP_Names.StationID = Map_Generators.STATIONID
WHERE
STID LIKE 'E%'
--AND
--Map_Generators.FuelSourcePrimary='Fossil'
GROUP BY StID
ORDER BY StID
;

SELECT
ISPname, FuelType, TechType, ZoneID, FuelCostGJ 
FROM (
SELECT 
FuelType,
TechType,
Status,
Region,
Generator,
FYE2030 AS FuelCostGJ -- Note: hard-coded year of 2030 etc.
FROM
Fuel_Prices
WHERE
CostScenario='Neutral'
AND Status='Existing'
) JOIN Map_ISP_Names ON ISPname = Generator 
;

DROP TABLE IF EXISTS Fuel_CoalGas2030;
CREATE TABLE Fuel_CoalGas2030 AS
SELECT 
TechType, FuelType, ZoneID AS RegionID,
round(median(FuelCostGJ*3.6),2) AS FuelCost -- Note: Conversion from $/GJ to $/MWh
FROM
(
SELECT 
FuelType,
TechType,
Status,
Region,
Generator,
FYE2030 AS FuelCostGJ -- Note: hard-coded year of 2030 etc.
FROM
Fuel_Prices
WHERE
CostScenario='Neutral'
AND Status='Existing'
)
JOIN
Map_ISP_Names ON ISPname = Generator 
GROUP BY FuelType, TechType, ZoneID
;

CREATE TABLE Fuel_RecipGas2030 (
    TechType TEXT,
    FuelType TEXT,
    RegionID,
    FuelCost
);

DROP VIEW IF EXISTS Fuel_ExistAll2030;
CREATE VIEW Fuel_ExistAll2030 AS
SELECT * FROM
(
SELECT * FROM Fuel_BioDist2030
 UNION 
SELECT * FROM Fuel_CoalGas2030
)
;

DROP TABLE Fuel_NewGas2030;
CREATE TABLE Fuel_NewGas2030 AS
SELECT
TechType, FuelType, RegionID,
FuelCostGJ*3.6 AS FuelCost -- Note: Conversion from $/GJ to $/MWh
FROM
(
SELECT
FuelType,
TechType,
Status,
Region,
Generator,
FYE2030 AS FuelCostGJ
FROM Fuel_Prices
WHERE
Status='NewEntrant'
AND CostScenario='Neutral'
AND FuelType LIKE '%Gas'
)
JOIN
(
SELECT RegionID, DemandRegion
FROM
STABLE_data.Regions
WHERE
RegionType='TxZone'
)
ON (DemandRegion = Region)
;

CREATE VIEW Fuel_NewAll2030 AS
SELECT * FROM
(
SELECT * FROM Fuel_NewGas2030
  UNION
SELECT
TechType,
FuelType,
Region AS RegionID,
FYE2030*3.6 AS FuelCost
FROM Fuel_Prices
WHERE
Status='NewEntrant'
AND CostScenario='Neutral'
AND FuelType LIKE '%Coal'
)
;



DROP VIEW IF EXISTS Fuel_Map;
CREATE VIEW Fuel_Map AS
SELECT * FROM
(
SELECT
StID, Tech_Parameters.TechType, Status, Tech_Parameters.FuelType, RegionID, FuelCost
FROM
Tech_Parameters JOIN Fuel_ExistAll2030 
 ON (
 Tech_Parameters.TechType = Fuel_ExistAll2030.TechType 
 AND
 lower(Tech_Parameters.FuelType) = lower(Fuel_ExistAll2030.FuelType)
 )
WHERE
Status='GenericExisting'
UNION
SELECT
StID, Tech_Parameters.TechType, Status, Tech_Parameters.FuelType, RegionID, FuelCost
FROM
Tech_Parameters JOIN Fuel_NewAll2030 
 ON (
 Tech_Parameters.TechType = Fuel_NewAll2030.TechType 
 AND
 lower(Tech_Parameters.FuelType) = lower(Fuel_NewAll2030.FuelType)
 )
WHERE
Status='NewEntrant'
)
;


-------- STABLE_data --------
-----------------------------


DROP TABLE IF EXISTS Technology_Types;
CREATE TABLE Technology_Types (
    TechTypeID         TEXT    PRIMARY KEY,
    TechTypeName       TEXT,
    TechTypeDescriptor TEXT,
    FuelType           TEXT   , -- REFERENCES Fuel_Sources (FuelTypeMinor),
    TechTypeCategory   TEXT    REFERENCES Technology_Categories (TechCategory),
    Existing           BOOLEAN,
    NewEntrant         BOOLEAN,
    IsRenewable        BOOLEAN,
    IsDispatchable     BOOLEAN,
    IsStorage          BOOLEAN,
    IsHybrid           BOOLEAN
);

DROP TABLE IF EXISTS Tech_Map;
CREATE TABLE Tech_Map (
    TechID     TEXT PRIMARY KEY,
    RegionType TEXT,
    TechTypeID TEXT,
    Status     TEXT,
    ISP_StID   TEXT
);

CREATE TABLE REZ_Trace_Map (
    RezID     TEXT PRIMARY KEY,
    Wind_Exi  TEXT,
    Solar_Exi TEXT,
    Wind_New  TEXT,
    Solar_New TEXT
);


SELECT
RegionID,
TechID as TechTypeID,
Efficiency,
FixedCost,
VariableCost,
OvernightCostPower,
OvernightCostEnergy,
CurtailmentCost,
LoadIncreaseCost,
LoadDecreaseCost,
InvestLifetime,
InvestRecovery,
MaxCapacity,
MaxEnergy,
StartLevel,
MaxEnergyToPowerRatio
FROM
Regions
JOIN
(
SELECT
TechID,
RegionType,
--ISP_StID,
STABLE_Format.*
FROM
Tech_Map JOIN ISP_2019_Input.STABLE_Format
ON ( ISP_StID = STABLE_Format.StID )
)
USING (RegionType)
;

SELECT
TechID,
Technology_Types.*
FROM 
Technology_Types
JOIN
Tech_Map
USING (TechTypeID)
;

DROP TABLE IF EXISTS Tech_All_Cap;
CREATE TABLE IF NOT EXISTS Tech_All_Cap AS
SELECT * FROM AEMO_Registration.Tech_REZone_Cap
UNION
SELECT * FROM AEMO_Registration.Tech_TxZone_Cap
;

SELECT sum(Capacity) FROM Tech_All_Cap;

DROP TABLE IF EXISTS Technologies;
CREATE TABLE Technologies AS
SELECT * FROM ISP_2019_Input.STABLE_Format
;
--SELECT * FROM ISP_2019_Input.STABLE_Format;


DROP TABLE IF EXISTS REZ_Build;
CREATE TABLE REZ_Build AS
SELECT
  RezID,
  WindLimit,
  SolarLimit,
  TotalBuildCap,
  ExpansionLimit,
  TransExpansionCost*1000 AS TransExpansionCostMW, -- conversion from $/kW to $/MW
  ConnectCost*1000 AS ConnectCostMW  -- conversion from $/kW to $/MW
FROM ISP_2019_Input.REZ_Build;

SELECT DISTINCT TechTypeID FROM Tech_Map
WHERE RegionType='TxZone' AND Status LIKE '%Existing'
;

SELECT DISTINCT TechTypeID FROM Tech_Map
WHERE RegionType='REZone' AND Status LIKE '%Existing'
;

DROP VIEW Existing_Cap;
--SELECT sum(Capacity) FROM (
CREATE VIEW Existing_Cap AS
SELECT TechID, RegionID, Capacity
FROM
(
SELECT * FROM
Tech_Map JOIN Tech_All_Cap USING (TechTypeID) 
WHERE RegionType='TxZone' AND Status LIKE '%Existing'
UNION
SELECT *
FROM
Tech_Map JOIN Tech_All_Cap USING (TechTypeID) 
WHERE RegionType='REZone' AND Status LIKE '%Existing'
)
--)
;

SELECT sum(Capacity) FROM Existing_Cap
;


SELECT TechID AS TechTypeID, RegionID, ISP_StID, Capacity as CapExist
--, STABLE_Format.*
FROM
(Existing_Cap JOIN Tech_Map USING (TechID))
--JOIN ISP_2019_Input.STABLE_Format ON (ISP_StID = ISP_2019_Input.STABLE_Format.StID )
;

------ STABLE_run_data ------
-----------------------------

-- Construct the run data for a STABLE instance

DROP TABLE IF EXISTS STABLE_run_data.Regions;
CREATE TABLE Regions AS
 SELECT * FROM STABLE_data.Regions
;

DROP TABLE IF EXISTS STABLE_run_data.Technology_Types;
CREATE TABLE Technology_Types AS
 SELECT * FROM STABLE_data.Technology_Types
;

DROP TABLE IF EXISTS REZ_Trace_Map;
CREATE TABLE REZ_Trace_Map AS
 SELECT * FROM STABLE_data.REZ_Trace_Map
;

DROP TABLE IF EXISTS Technologies;
CREATE TABLE Technologies AS
 SELECT * FROM ISP_2019_Input.STABLE_Format
;

DROP TABLE IF EXISTS Fuel_Map;
CREATE TABLE Fuel_Map AS
 SELECT * FROM ISP_2019_Input.Fuel_Map
;


DROP TABLE IF EXISTS Tech_Map;
CREATE TABLE Tech_Map AS
 SELECT * FROM STABLE_data.Tech_Map
;

DROP TABLE IF EXISTS Existing_Cap;
CREATE TABLE Existing_Cap AS
 SELECT * FROM STABLE_data.Existing_Cap
;

SELECT sum(Capacity) FROM Existing_Cap
;


DROP VIEW IF EXISTS Fuel_Costs_All;
CREATE VIEW Fuel_Costs_All AS
SELECT TechID, TechTypeID, RegionID, Tech_Map.Status, 
NULL AS FuelType,
0 AS FuelCost
FROM
Regions JOIN Tech_Map USING (RegionType)
JOIN Technology_Types USING (TechTypeID)
WHERE TechTypeCategory != 'Combustion'
UNION
SELECT TechID, TechTypeID, RegionID, Tech_Map.Status, FuelType,
CASE WHEN TechTypeID LIKE '%Waste' THEN 0 ELSE FuelCost END AS FuelCost
FROM
Tech_Map JOIN Fuel_Map ON (Tech_Map.ISP_StID = Fuel_Map.StID)
;

SELECT TechID, TechTypeID, RegionID, Tech_Map.Status, 
NULL AS FuelType,
0 AS FuelCost
FROM
Regions JOIN Tech_Map USING (RegionType)
JOIN Technology_Types USING (TechTypeID)
WHERE TechTypeCategory != 'Combustion'
;

SELECT TechID, TechTypeID, RegionID, Tech_Map.Status, FuelType, FuelCost
FROM
Tech_Map JOIN Fuel_Map ON (Tech_Map.ISP_StID = Fuel_Map.StID)
;

DROP VIEW IF EXISTS Tech_New_Parameters;
CREATE VIEW Tech_New_Parameters AS
SELECT Tech_Map.TechID,
Tech_Map.Status,
Technology_Types.TechTypeCategory,
Technology_Types.IsRenewable,
Technology_Types.IsDispatchable,
Technology_Types.IsStorage,
Technology_Types.IsHybrid,
Fuel_Costs_All.RegionID,
NULL as CapExist,
FuelCost,
Technologies.*
FROM
Tech_Map JOIN Technologies ON (Tech_Map.ISP_StID = Technologies.StID )
LEFT JOIN Fuel_Costs_All ON (Tech_Map.TechID = Fuel_Costs_All.TechID)
JOIN Technology_Types ON (Tech_Map.TechTypeID = Technology_Types.TechTypeID)
WHERE Tech_Map.Status='NewEntrant'
;

--SELECT sum(CapExist) FROM (
DROP VIEW IF EXISTS Tech_Existing_Parameters;
CREATE VIEW Tech_Existing_Parameters AS
SELECT Tech_Map.TechID, 
Tech_Map.Status,
Technology_Types.TechTypeCategory,
Technology_Types.IsRenewable,
Technology_Types.IsDispatchable,
Technology_Types.IsStorage,
Technology_Types.IsHybrid,
Existing_Cap.RegionID,
Existing_Cap.Capacity as CapExist,
FuelCost,
Technologies.*
FROM
(Existing_Cap JOIN Tech_Map USING (TechID))
JOIN Technologies ON (Tech_Map.ISP_StID = Technologies.StID )
LEFT JOIN Fuel_Costs_All ON (Tech_Map.TechID = Fuel_Costs_All.TechID AND Existing_Cap.RegionID = Fuel_Costs_All.RegionID)
JOIN Technology_Types ON (Tech_Map.TechTypeID = Technology_Types.TechTypeID)
--)
;

DROP TABLE IF EXISTS Tech_All_Parameters;
CREATE TABLE Tech_All_Parameters (
    TechID                TEXT,
    Status                TEXT,
    TechTypeCategory      TEXT,
    IsRenewable           NUM,
    IsDispatchable        NUM,
    IsStorage             NUM,
    IsHybrid              NUM,
    RegionID              TEXT,
    CapacityExist         REAL,
    FuelCost              REAL,
    StID                  TEXT,
    Efficiency            REAL,
    FixedCost             REAL,
    VariableCost          REAL,
    OvernightCostPower    REAL,
    OvernightCostEnergy   REAL,
    CurtailmentCost       REAL,
    LoadIncreaseCost      REAL,
    LoadDecreaseCost      REAL,
    InvestLifetime        INT,
    InvestRecovery        INT,
    MaxCapacity           REAL,
    MaxEnergy             REAL,
    StartLevel            REAL,
    MaxEnergyToPowerRatio REAL
);

INSERT INTO Tech_All_Parameters
SELECT * FROM Tech_Existing_Parameters
UNION
SELECT * FROM Tech_New_Parameters
;

SELECT
RegionID as Nodes,
DemandRegion,
IncludeLevel as IncludeFlag
FROM Regions
WHERE DemandRegion IS NOT NULL
AND IncludeLevel=1
;

DROP TABLE IF EXISTS Arcs;
CREATE TABLE Arcs (
  ZoneConnectorName TEXT,
  FromNode TEXT,
  ToNode   TEXT,
  ForwardLimit  REAL,
  ReverseLimit  REAL
)
;

DROP TABLE IF EXISTS REZ_Build;
CREATE TABLE REZ_Build AS
SELECT * FROM STABLE_data.REZ_Build
;

DROP TABLE IF EXISTS TxZ_Connect_Cost;
CREATE TABLE TxZ_Connect_Cost (
    TechTypeID TEXT,
    TechID TEXT,
    DemandRegion TEXT,
    ISP_StID TEXT,
    ConnectCost REAL
)
;

DROP TABLE IF EXISTS TxZ_Connect_Cost;
--CREATE TABLE TxZ_Connect_Cost AS
SELECT TechID, DemandRegion,
ConnectCost*1000 as ConnectCost -- Conversion to $/MW
FROM ISP_2019_Input.TxZ_Connect_Cost
;



DROP TABLE IF EXISTS Tech_Scenario;
CREATE TABLE Tech_Scenario AS
SELECT * FROM STABLE_data.Tech_Scenario
;

DROP TABLE IF EXISTS Demand_Scenario;
CREATE TABLE Demand_Scenario AS
SELECT * FROM STABLE_data.Demand_Scenario
;

DROP TABLE IF EXISTS Tech_Inertia_Value;
CREATE TABLE Tech_Inertia_Value AS
SELECT * FROM STABLE_data.Tech_Inertia_Value;

DROP TABLE IF EXISTS Tech_Inertia_Value;
CREATE TABLE Tech_Inertia_Value AS
SELECT * FROM STABLE_data.Tech_Inertia_Value;

DROP TABLE IF EXISTS Tech_Inertia_Storage;
CREATE TABLE Tech_Inertia_Storage AS
SELECT * FROM STABLE_data.Tech_Inertia_Require;

--DROP TABLE IF EXISTS Carbon_Price;
DROP TABLE IF EXISTS Carbon_Param;
CREATE TABLE Carbon_Param AS
SELECT * FROM STABLE_data.Carbon_Param
;


CREATE TABLE Carbon_Content AS
SELECT DISTINCT 
FuelType,
0 AS CarbonContent
FROM Technology_Types;

DROP TABLE IF EXISTS H2_Tech;
CREATE TABLE H2_Tech AS
SELECT * FROM STABLE_data.H2_Tech;

DROP TABLE IF EXISTS RE_Targets;
CREATE TABLE RE_Targets AS
SELECT * FROM STABLE_data.RE_Targets;


-- In STABLE_data DB:
CREATE TABLE H2_Tech (
    TechID                TEXT,
    Status                TEXT,
    TechTypeDescriptor    TEXT,
    FuelType              TEXT,
    TechTypeCategory      TEXT,
    IsRenewable           NUM,
    IsDispatchable        NUM,
    IsStorage             NUM,
    IsHybrid              NUM,
    RegionID              TEXT,
    CapacityExist         REAL,
    FuelCost              REAL,
    StID                  TEXT,
    Efficiency            REAL,
    FixedCost             REAL,
    VariableCost          REAL,
    OvernightCostPower    REAL,
    OvernightCostEnergy   REAL,
    CurtailmentCost       REAL,
    LoadIncreaseCost      REAL,
    LoadDecreaseCost      REAL,
    InvestLifetime        INT,
    InvestRecovery        INT,
    MaxCapacity           REAL,
    MaxEnergy             REAL,
    StartLevel            REAL,
    MaxEnergyToPowerRatio REAL
);

DROP TABLE IF EXISTS RE_Targets;
CREATE TABLE RE_Targets (
    Region TEXT,
    ScenarioName TEXT,
    ScenarioYear INT,
    MinRET REAL
);

DROP TABLE IF EXISTS VPP_Battery_Uptake;
CREATE TABLE VPP_Battery_Uptake (
    Year INT,
    Region TEXT,
    ValueType TEXT,
    Unit TEXT,
    VPP_Scenario INT,
    Data REAL
)
;

DROP TABLE IF EXISTS Map_Output;
CREATE TABLE Map_Output (
TechID TEXT,
BlackCoal INTEGER,
BrownCoal INTEGER,
GasOther INTEGER,
CCGT INTEGER,
OCGT INTEGER,
Distillate INTEGER,
Biomass INTEGER,
Hydro INTEGER,
Wind INTEGER,
SolarPV INTEGER,
SolarThermal INTEGER,
PumpedHydro INTEGER,
Battery INTEGER,
BatteryVPP INTEGER
)
;


-- In STABLE_run_data DB:
DROP TABLE IF EXISTS H2_Tech;
CREATE TABLE H2_Tech AS
SELECT * FROM STABLE_data.H2_Tech;

DROP TABLE IF EXISTS RE_Targets;
CREATE TABLE RE_Targets AS
SELECT * FROM STABLE_data.RE_Targets;

DROP TABLE IF EXISTS VPP_Battery_Uptake;
CREATE TABLE VPP_Battery_Uptake AS
SELECT * FROM STABLE_data.VPP_Battery_Uptake;

DROP TABLE IF EXISTS Map_Output;
CREATE TABLE Map_Output AS
SELECT * FROM STABLE_data.Map_Output;
