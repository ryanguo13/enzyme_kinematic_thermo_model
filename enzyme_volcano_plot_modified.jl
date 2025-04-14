include("thermo_kinetics_modified.jl")
include("compare_models.jl")

# 默认参数设置 - 与Python代码保持一致
params = ThermoKineticsParams(
    ΔGT = -40.0,   # kJ/mol
    ΔG1 = -15.0,   # kJ/mol
    R = 8.314,     # J/mol·K
    T = 300.0,     # K，与Python代码一致
    α1 = 0.5,      # 与Python代码一致
    α2 = 0.5,      # 与Python代码一致
    k10 = 1.0,     # 与Python代码一致
    k20 = 1.0,     # 与Python代码一致
    ET = 0.01,     # mM，与Python代码一致
    model_type = "three_step"  # 默认使用三步骤模型
)

# 底物浓度列表 - 调整为对数刻度[0.1, 1, 10, 100]以匹配参考图像
substrate_concentrations = [0.1, 1.0, 10.0, 100.0]

# 绘制三步骤模型火山图
p_three = plot_activity_volcano(params, substrate_concentrations)

# 保存三步骤模型火山图
savefig(p_three, "result/enzyme_activity_volcano_three_step.png")
println("三步骤模型火山图已保存为 result/enzyme_activity_volcano_three_step.png")

# 更新参数为四步骤模型
params_four = ThermoKineticsParams(
    ΔGT = -40.0,   # kJ/mol
    ΔG1 = -15.0,   # kJ/mol
    R = 8.314,     # J/mol·K
    T = 300.0,     # K
    α1 = 0.5,      
    α2 = 0.5,      
    k10 = 1.0,     
    k20 = 1.0,     
    ET = 0.01,     # mM
    model_type = "four_step"  # 使用四步骤模型
)

# 绘制四步骤模型火山图
p_four = plot_activity_volcano(params_four, substrate_concentrations)

# 保存四步骤模型火山图
savefig(p_four, "result/enzyme_activity_volcano_four_step.png")
println("四步骤模型火山图已保存为 result/enzyme_activity_volcano_four_step.png")

# 绘制三步骤和四步骤模型的对比图
p_comparison = compare_volcano_plots(substrate_concentrations, params)

# 保存对比图
savefig(p_comparison, "result/enzyme_activity_volcano_comparison.png")
println("三步骤和四步骤模型的对比图已保存为 result/enzyme_activity_volcano_comparison.png")

# 输出参数设置，便于比对
println("\n参数设置:")
println("ΔGT = ", params.ΔGT, " kJ/mol")
println("ΔG1 = ", params.ΔG1, " kJ/mol")
println("R = ", params.R, " J/mol·K")
println("T = ", params.T, " K")
println("RT = ", params.R / 1000 * params.T, " kJ/mol")
println("α1 = ", params.α1)
println("α2 = ", params.α2)
println("k10 = ", params.k10, " 1/s")
println("k20 = ", params.k20, " 1/s")
println("ET = ", params.ET, " mM")

# 计算结果
results = calculate_kinetics(params)
println("\n计算结果:")
println("Km = ", results.Km, " μM")
println("v_max = ", results.v_max, " μM/s")