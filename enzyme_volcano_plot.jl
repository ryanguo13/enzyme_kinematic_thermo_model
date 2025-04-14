include("thermo_kinetics.jl")

# 默认参数设置
params = ThermoKineticsParams(
    ΔGT = -40.0,   # kJ/mol
    ΔG1 = -15.0,   # kJ/mol (确保ΔG2合理)
    R = 8.314,      # J/mol·K
    T = 298.0,      # K
    α1 = 0.5,
    α2 = 0.6,
    k10 = 1.0,     # 更合理的基准速率(1/s)
    k20 = 0.1,     # 更合理的基准速率(1/s)
    ET = 1.0        # mM
)

# 底物浓度列表 - 调整为对数刻度[0.1, 1, 10, 100]以匹配参考图像
substrate_concentrations = [0.1, 1.0, 10.0, 100.0]

# 绘制火山图
p = plot_activity_volcano(params, substrate_concentrations)

# 保存图像
savefig(p, "result/enzyme_activity_volcano.svg")
println("火山图已保存为 result/enzyme_activity_volcano.svg")
