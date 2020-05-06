# %% Write the formatted output to file:

import CSV
import XLSX

prepare_df_xlsx(df) = ( collect(DataFrames.eachcol(df)), DataFrames.names(df) )

scenario_timestamp = scen_settings[:scen]
# scenario_timestamp = "Raw-Output-Template-Testing"

# cf. results_summary.jl
# We annotate with the label "Direct" since we are not summarising the results into technology groupings:
xlsx_output_file = joinpath(resultsdir,"STABLE_summary-Direct-$(scenario_timestamp).xlsx")

CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Demand.csv"),dfDict["load"])
CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Gen-TxZ.csv"),resSplit[:TxZ_GEN])
CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Gen-REZ.csv"),resSplit[:REZ_GEN_CU])

CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Storage.csv"),resSplit[:STORAGE])
CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Flow.csv"),resSplit[:FLOW])


# Write capacity data and create place-holders for larger datasets:
XLSX.openxlsx(xlsx_output_file, mode="w") do xf
    XLSX.rename!(xf[1],"CAPACITY_GEN")
        XLSX.writetable!(xf["CAPACITY_GEN"], prepare_df_xlsx(resSplit[:CAPACITY_GEN])...)
    XLSX.addsheet!(xf,"CAPACITY_STO")
        XLSX.writetable!(xf["CAPACITY_STO"], prepare_df_xlsx(resSplit[:CAPACITY_STO])...)
    XLSX.addsheet!(xf,"CAPACITY_REZ_EXP")
        XLSX.writetable!(xf["CAPACITY_REZ_EXP"], prepare_df_xlsx(resSplit[:CAPACITY_REZ_EXP])...)
    XLSX.addsheet!(xf,"CAPACITY_SYNC")
        XLSX.writetable!(xf["CAPACITY_SYNC"], prepare_df_xlsx(resSplit[:N_SYNC])...)
    XLSX.addsheet!(xf,"CAPACITY_H2")
        XLSX.writetable!(xf["CAPACITY_H2"], prepare_df_xlsx(resSplit[:CAPACITY_H2])...)
end

# Use inter-region flow DataFrame and filter for each region to obtain _net_ flow _from_ region:
interflow_dict = Dict{String, DataFrame}()
for DR in dtr.sets[:DemandRegions]
    interflow_dict[DR] = @linq resSplit[:INTERFLOW] |>
              where(:FromRegion .== DR) |>
              by(:Hours, Level = sum(:FLOW))
end

XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
    for DR in dtr.sets[:DemandRegions]
        XLSX.addsheet!(xf,"INTERFLOW_"*DR)
           XLSX.writetable!(xf["INTERFLOW_"*DR], prepare_df_xlsx(interflow_dict[DR])...)
    end
end

#=
XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
    XLSX.addsheet!(xf,"OPDEMAND")
    XLSX.addsheet!(xf,"GEN_TxZ")
        # XLSX.writetable!(xf["GEN_TxZ"], prepare_df_xlsx(resSplit[:TxZ_GEN])...)
    XLSX.addsheet!(xf,"GEN_REZ")
        # XLSX.writetable!(xf["GEN_REZ"], prepare_df_xlsx(resSplit[:REZ_GEN])...)
    XLSX.addsheet!(xf,"STORAGE")
        # XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
    XLSX.addsheet!(xf,"FLOW")
        # XLSX.writetable!(xf["FLOW"], prepare_df_xlsx(resSplit[:FLOW])...)
end
=#

# CSV.write(joinpath(resultsdir,"$(scenario_timestamp)-Dispatch.csv"),resSplit[:DISPATCH])
# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["DISPATCH"], prepare_df_xlsx(resSplit[:DISPATCH])...)
# end

# XLSX.openxlsx(xlsx_output_file, mode="rw") do xf
#     XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
# end

# XLSX.writetable(xlsx_output_file,
#     CAPACITY_GEN=prepare_df_xlsx(df_N_TECH_Include),
#     CAPACITY_STO=prepare_df_xlsx(resSplit[:N_STO_all])
#     )
