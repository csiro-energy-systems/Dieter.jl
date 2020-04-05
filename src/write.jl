# %% Write the formatted output to file:

import CSV
import XLSX

prepare_df_xlsx(df) = ( collect(DataFrames.eachcol(df)), DataFrames.names(df) )


if Base.Sys.isapple()
      resultsdir = joinpath(ENV["HOME"],"Documents/Projects/ESM/","results_STABLE")
elseif Base.Sys.iswindows()
      resultsdir = joinpath("F:\\STABLE\\","results_STABLE")
end

xlsx_output_file = joinpath(resultsdir,"STABLE_summary-$(run_timestamp).xlsx")

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Demand.csv"),dfDict["load"])
CSV.write(joinpath(resultsdir,"$(run_timestamp)-Gen-TxZ.csv"),resSplit[:TxZ_GEN])
CSV.write(joinpath(resultsdir,"$(run_timestamp)-Gen-REZ.csv"),resSplit[:REZ_GEN_CU])

CSV.write(joinpath(resultsdir,"$(run_timestamp)-Storage.csv"),resSplit[:STORAGE])
CSV.write(joinpath(resultsdir,"$(run_timestamp)-Flow.csv"),resSplit[:FLOW])


# Write capacity data and create place-holders for larger datasets:
XLSX.openxlsx(xlsx_output_file, mode="w") do xf
    XLSX.rename!(xf[1],"OPDEMAND")
    XLSX.addsheet!(xf,"CAPACITY_GEN")
        XLSX.writetable!(xf["CAPACITY_GEN"], prepare_df_xlsx(resSplit[:CAPACITY_GEN])...)
    XLSX.addsheet!(xf,"CAPACITY_STO")
        XLSX.writetable!(xf["CAPACITY_STO"], prepare_df_xlsx(resSplit[:CAPACITY_STO])...)
    XLSX.addsheet!(xf,"CAPACITY_REZ_EXP")
        XLSX.writetable!(xf["CAPACITY_REZ_EXP"], prepare_df_xlsx(resSplit[:CAPACITY_REZ_EXP])...)
    XLSX.addsheet!(xf,"CAPACITY_SYNC")
        XLSX.writetable!(xf["CAPACITY_SYNC"], prepare_df_xlsx(resSplit[:N_SYNC])...)
    XLSX.addsheet!(xf,"GEN_TxZ")
        # XLSX.writetable!(xf["GEN_TxZ"], prepare_df_xlsx(resSplit[:TxZ_GEN])...)
    XLSX.addsheet!(xf,"GEN_REZ")
        # XLSX.writetable!(xf["GEN_REZ"], prepare_df_xlsx(resSplit[:REZ_GEN])...)
    XLSX.addsheet!(xf,"STORAGE")
        # XLSX.writetable!(xf["STORAGE"], prepare_df_xlsx(resSplit[:STORAGE])...)
    XLSX.addsheet!(xf,"FLOW")
        # XLSX.writetable!(xf["FLOW"], prepare_df_xlsx(resSplit[:FLOW])...)
    for DR in DemandRegions
        XLSX.addsheet!(xf,"INTERFLOW_"*DR)
           XLSX.writetable!(xf["INTERFLOW_"*DR], prepare_df_xlsx(df_regflow[DR])...)
    end
end


# CSV.write(joinpath(resultsdir,"$(run_timestamp)-Dispatch.csv"),resSplit[:DISPATCH])
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
