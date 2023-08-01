using JuliaPhyx
jp = JuliaPhyx

res = jp.PhyOutput()

@time channel = jp.getchan(res, 385, 0, "max")