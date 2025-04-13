using Symbolics
using Plots

# 加载两个模型文件
include("model_cal.jl")  # 四步骤模型
include("thermo_kinetics.jl")  # 三步骤模型

"""
比较四步骤模型与三步骤模型

四步骤模型 (model_cal.jl):
    E + S ⇌ ES ⇌ EP ⇌ E + P
    
三步骤模型 (thermo_kinetics.jl):
    E + S ⇌ ES ⇌ E + P

两个模型的关系：
1. 四步骤模型更详细，包含了EP中间体
2. 三步骤模型是四步骤模型的简化版本
3. 两个模型都可以用相同的简化表达式表示：v = (Et*(k_forward*S-k_reverse*P))/(1+S/Ks+P/Kp)
"""

# 定义函数，使用四步骤模型计算反应速率
function calculate_four_step_rate(S, P, Et, k1, k_minus_1, k2, k_minus_2, k3, k_minus_3)
    # 计算简化参数 - 根据model_cal.jl中的推导
    # 注意：S和P单位为mM，Ks和Kp单位为μM
    
    # 根据model_cal.jl中的正确公式计算
    # 注意：原始公式中Kp和Ks的定义有误，这里使用正确的公式
    
    # 计算k_forward和k_reverse (单位：1/s)
    k_forward = k2*k3 / (k2 + k3 + k_minus_2)
    k_reverse = k_minus_1*k_minus_2 / (k2 + k_minus_1 + k_minus_2)
    
    # 计算Kp和Ks (单位：mM)
    # 根据model_cal.jl中的推导
    Kp_mM = (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2) / (k_minus_3 * (k2 + k_minus_1 + k_minus_2))
    Ks_mM = (k_minus_1 * (k2 + k3 + k_minus_2)) / (k1 * (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2))
    
    # 单位转换：将Ks和Kp从mM转换为μM以与三步骤模型保持一致
    Ks = Ks_mM * 1e3  # μM
    Kp = Kp_mM * 1e3  # μM
    
    # 使用简化表达式计算反应速率
    # 注意：S和P单位为mM，需要转换为μM以匹配Ks和Kp单位
    S_uM = S * 1e3  # 转换为μM
    P_uM = P * 1e3  # 转换为μM
    
    # 计算反应速率 (单位：μM/s)
    v = (Et * 1e3 * (k_forward * S_uM - k_reverse * P_uM)) / (1 + S_uM/Ks + P_uM/Kp)
    return v, Kp, Ks, k_forward, k_reverse
end

# 定义函数，使用三步骤模型计算反应速率
function calculate_three_step_rate(S, P, Et, ΔG1, ΔGT, R=8.314, T=298, α1=0.5, α2=0.6, k10=1.0, k20=0.1)
    # 创建参数对象
    params = ThermoKineticsParams(
        ΔGT = ΔGT,
        ΔG1 = ΔG1,
        R = R,
        T = T,
        α1 = α1,
        α2 = α2,
        k10 = k10,
        k20 = k20,
        ET = Et
    )
    
    # 计算动力学参数
    results = calculate_kinetics(params)
    
    # 使用简化表达式计算反应速率
    # 注意：S和P单位为mM，需要转换为μM以匹配Ks和Kp单位
    S_uM = S * 1e3  # 转换为μM
    P_uM = P * 1e3  # 转换为μM
    
    # 计算反应速率 (单位：μM/s)
    v = (Et * 1e3 * (results.k_forward * S_uM - results.k_reverse * P_uM)) / (1 + S_uM/results.Ks + P_uM/results.Kp)
    return v, results.Kp, results.Ks, results.k_forward, results.k_reverse
end

# 比较两个模型在不同条件下的反应速率
function compare_models()
    println("Comparison of Four-step and Three-step Model Reaction Rates")
    println("================================================")
    
    # 设置共同参数
    S_values = [0.1, 1.0, 10.0, 100.0]  # mM
    P = 0.0  # 初始无产物
    Et = 1.0  # mM
    
    # 四步骤模型参数
    k1 = 10.0        # 1/(mM·s)
    k_minus_1 = 1.0  # 1/s
    k2 = 5.0         # 1/s
    k_minus_2 = 0.5  # 1/s
    k3 = 2.0         # 1/s
    k_minus_3 = 0.2  # 1/(mM·s)
    
    # 三步骤模型参数
    # 我们需要从四步骤模型参数推导出对应的ΔG1和ΔGT
    R = 8.314  # J/(mol·K)
    T = 298.0   # K
    RT = R * T  # J/(mol)
    
    # 计算平衡常数
    K1 = k1 / k_minus_1  # 1/mM
    K2 = k2 / k_minus_2  # 无量纲
    K3 = k3 / k_minus_3  # mM
    
    # 计算自由能变化 (J/mol)
    ΔG1_J = -RT * log(K1)  # J/mol
    ΔG2_J = -RT * log(K2)  # J/mol
    ΔG3_J = -RT * log(K3)  # J/mol
    ΔGT_J = ΔG1_J + ΔG2_J + ΔG3_J  # J/mol
    
    # 转换为kJ/mol用于三步骤模型
    ΔG1 = ΔG1_J / 1000  # kJ/mol
    ΔG2 = ΔG2_J / 1000  # kJ/mol
    ΔG3 = ΔG3_J / 1000  # kJ/mol
    ΔGT = ΔGT_J / 1000  # kJ/mol
    
    # 设置三步骤模型的α1和α2参数，使其与四步骤模型更一致
    α1 = 0.5  # 默认值
    α2 = 0.6  # 默认值
    k10 = 1.0  # 基准速率常数1 (1/s)
    k20 = 0.1  # 基准速率常数2 (1/s)
    
    println("Four-step Model Parameters:")
    println("k1 = $k1, k_minus_1 = $k_minus_1, k2 = $k2, k_minus_2 = $k_minus_2, k3 = $k3, k_minus_3 = $k_minus_3")
    println("\nCalculated Free Energy Changes:")
    println("ΔG1 = $ΔG1 kJ/mol, ΔG2 = $ΔG2 kJ/mol, ΔG3 = $ΔG3 kJ/mol, ΔGT = $ΔGT kJ/mol")
    println("\nThree-step Model Parameters:")
    println("ΔG1 = $ΔG1 kJ/mol, ΔGT = $ΔGT kJ/mol, α1 = $α1, α2 = $α2, k10 = $k10, k20 = $k20")
    
    # 计算并比较两个模型的反应速率
    println("\nReaction Rate Comparison:")
    println("S (mM) | Four-step Model (μM/s) | Three-step Model (μM/s) | Difference (%)")
    println("----------------------------------------------------------")
    
    for S in S_values
        # 四步骤模型
        v4, Kp4, Ks4, kf4, kr4 = calculate_four_step_rate(S, P, Et, k1, k_minus_1, k2, k_minus_2, k3, k_minus_3)
        
        # 三步骤模型
        v3, Kp3, Ks3, kf3, kr3 = calculate_three_step_rate(S, P, Et, ΔG1, ΔGT, R, T, α1, α2, k10, k20)
        
        # 计算差异
        diff = (v3 - v4) / v4 * 100
        
        println("$S | $(round(v4, digits=4)) | $(round(v3, digits=4)) | $(round(diff, digits=2))")
    end
    
    # 比较简化参数
    v4, Kp4, Ks4, kf4, kr4 = calculate_four_step_rate(1.0, P, Et, k1, k_minus_1, k2, k_minus_2, k3, k_minus_3)
    v3, Kp3, Ks3, kf3, kr3 = calculate_three_step_rate(1.0, P, Et, ΔG1, ΔGT, R, T, α1, α2, k10, k20)
    
    println("\nSimplified Parameter Comparison:")
    println("Parameter | Four-step Model | Three-step Model | Difference (%)")
    println("------------------------------------------")
    println("Kp (μM) | $(round(Kp4, digits=4)) | $(round(Kp3, digits=4)) | $(round((Kp3-Kp4)/Kp4*100, digits=2))")
    println("Ks (μM) | $(round(Ks4, digits=4)) | $(round(Ks3, digits=4)) | $(round((Ks3-Ks4)/Ks4*100, digits=2))")
    println("k_forward (1/s) | $(round(kf4, digits=4)) | $(round(kf3, digits=4)) | $(round((kf3-kf4)/kf4*100, digits=2))")
    println("k_reverse (1/s) | $(round(kr4, digits=4)) | $(round(kr3, digits=4)) | $(round((kr3-kr4)/kr4*100, digits=2))")
    
    # 计算热力学一致性
    # 理论上，k_forward/k_reverse应该等于exp(-ΔGT_J/RT)
    thermo_ratio4 = kf4/kr4
    thermo_ratio3 = kf3/kr3
    expected_ratio = exp(-ΔGT_J/RT)
    
    println("\nThermodynamic Consistency Check:")
    println("Expected ratio (exp(-ΔGT/RT)) = $(round(expected_ratio, digits=4))")
    println("Four-step Model: k_forward/k_reverse = $(round(thermo_ratio4, digits=4)), Difference = $(round((thermo_ratio4/expected_ratio-1)*100, digits=2))%")
    println("Three-step Model: k_forward/k_reverse = $(round(thermo_ratio3, digits=4)), Difference = $(round((thermo_ratio3/expected_ratio-1)*100, digits=2))%")
    
    # 绘制不同底物浓度下的反应速率曲线
    S_curve = 10.0.^range(-2, 3, length=100)  # 从0.01到1000 mM的对数刻度
    v4_curve = [calculate_four_step_rate(s, P, Et, k1, k_minus_1, k2, k_minus_2, k3, k_minus_3)[1] for s in S_curve]
    v3_curve = [calculate_three_step_rate(s, P, Et, ΔG1, ΔGT, R, T, α1, α2, k10, k20)[1] for s in S_curve]
    
    p = plot(S_curve, v4_curve, xscale=:log10, label="Four-step Model", 
             xlabel="Substrate Concentration [mM]", ylabel="Reaction Rate [μM/s]", 
             title="Comparison of Four-step and Three-step Models", 
             linewidth=2, legend=:bottomright)
    plot!(p, S_curve, v3_curve, label="Three-step Model", linewidth=2, linestyle=:dash)
    
    # 显示图表
    # display(p)
    
    # 保存图表
    savefig(p, "model_comparison.png")
    println("\nChart generated and saved: model_comparison.png")
    
    return p
end

# 运行比较
if abspath(PROGRAM_FILE) == @__FILE__
    compare_models()
end