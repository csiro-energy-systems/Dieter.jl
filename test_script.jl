rdir = "./results_of_test_case/"

df = joinpath(rdir,"N.feather") |> Feather.read

techs=["Solar PV", "Wind Onshore", "Wind Offshore", "Biomasse"]

sector=:h2

# capacity = get_result(rdir, :N, filter_by=create_filter_dict(sector))
capacity = joinpath(rdir,"N.feather") |> Feather.read


renewable = filter(row-> row[:fuel] in techs, capacity)
sort!(renewable, [:min_res, sector])

renewable = by(renewable, [:min_res, sector, :fuel], :Value => sum)


color_dict = Dict()
default_colour = ColorSchemes.leonardo;
# ColorScheme([Colors.RGB(0.0, 0.0, 0.0), Colors.RGB(1.0, 1.0, 1.0)],
      # "custom", "twotone, black and white");
color_dict["Wind Offshore"] = [default_colour[10]];
color_dict["Wind Onshore"] = [default_colour[20]];
color_dict["Solar PV"] = [default_colour[30]];

plt = plot(legend=:topleft)
fuels = unique(renewable[!,:fuel])
len_sector = unique(capacity[!,sector]) |> length

# for f in fuels
f = "Solar PV"
f = "Wind Onshore"
f = "Wind Offshore"

df = renewable[renewable.fuel .== f, :]
c = c_gradient(color_dict[f], len_sector) |> permutedims
x = df[!,:min_res]
y =  df[!,:Value_sum]
g = df[!,sector]
l = vcat(f, fill("", len_sector-1)) |> permutedims
t = "Generation Capacity"
if sum(df[!,:Value_sum]) > 0
    # plot!(x, y; group=g, c=c, label=l, width=2, title=t, grid=false, kwargs...)
    plot!(x, y; group=g, c=c, label=l, width=10, title=t, grid=false, marker = (:hexagon, 20, 0.6, :green, stroke(3, 0.2, :black, :dot)))
end
# end
