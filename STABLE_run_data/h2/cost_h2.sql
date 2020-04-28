SELECT
StID AS Technologies,
Tech_Costs.Scenario,
Tech_Costs.FYE2020*1000 AS FYE2020, -- conversion from $/kW to $/MW
Tech_Costs.FYE2025*1000 AS FYE2025, -- conversion from $/kW to $/MW
Tech_Costs.FYE2030*1000 AS FYE2030, -- conversion from $/kW to $/MW
Tech_Costs.FYE2035*1000 AS FYE2035, -- conversion from $/kW to $/MW
Tech_Costs.FYE2040*1000 AS FYE2040, -- conversion from $/kW to $/MW
Tech_Costs.FYE2045*1000 AS FYE2045, -- conversion from $/kW to $/MW
Tech_Costs.FYE2050*1000 AS FYE2050  -- conversion from $/kW to $/MW
FROM
Tech_Costs
WHERE StID IN ('N_RecipH2', 'N_Electrolyser')
