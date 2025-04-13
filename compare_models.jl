# 比较三步骤和四步骤模型的火山图

include("thermo_kinetics_modified.jl")

"""
比较三步骤和四步骤模型的火山图

参数:
- S_range: 底物浓度范围 (mM)
- params: 基本参数设置

返回:
- 包含两种模型火山图的对比图
"""
function compare_volcano_plots(S_range::Vector{Float64}, params::ThermoKineticsParams)
    # 生成ΔG1和ΔGT网格
    ΔG1_range = LinRange(-20.0, 20.0, 100)  # 调整范围与参考图像一致
    ΔGT_range = LinRange(-40.0, 40.0, 100)  # 调整范围与参考图像一致
    
    # 准备存储结果
    v_matrix_three = zeros(length(ΔG1_range), length(ΔGT_range), length(S_range))
    v_matrix_four = zeros(length(ΔG1_range), length(ΔGT_range), length(S_range))
    km_lines_three = []
    km_lines_four = []
    
    # 计算三步骤模型
    three_step_params = ThermoKineticsParams(
        ΔGT = params.ΔGT, ΔG1 = params.ΔG1,
        R = params.R, T = params.T,
        α1 = params.α1, α2 = params.α2,
        k10 = params.k10, k20 = params.k20,
        ET = params.ET,
        model_type = "three_step"
    )
    
    # 计算四步骤模型
    four_step_params = ThermoKineticsParams(
        ΔGT = params.ΔGT, ΔG1 = params.ΔG1,
        R = params.R, T = params.T,
        α1 = params.α1, α2 = params.α2,
        k10 = params.k10, k20 = params.k20,
        ET = params.ET,
        model_type = "four_step"
    )
    
    # 计算三步骤模型的活性
    for (i, ΔG1) in enumerate(ΔG1_range)
        for (j, ΔGT) in enumerate(ΔGT_range)
            # 更新参数
            new_params = ThermoKineticsParams(
                ΔGT = ΔGT, ΔG1 = ΔG1,
                R = params.R, T = params.T,
                α1 = params.α1, α2 = params.α2,
                k10 = params.k10, k20 = params.k20,
                ET = params.ET,
                model_type = "three_step"
            )
            
            # 计算结果
            res = calculate_kinetics(new_params)
            
            # 计算各底物浓度下的活性
            for (k, S) in enumerate(S_range)
                # 计算酶活性（单位：μM/s）
                RT = params.R / 1000 * params.T  # kJ/mol
                g1 = exp(ΔG1/RT)
                gT = exp(ΔGT/RT)
                K = params.k20/params.k10/gT^params.α2 * g1^(params.α1+params.α2-1)
                Km_val = g1*(1+K)
                Km = Km_val * 1e3  # 转换为μM
                
                S_uM = S * 1000.0  # mM到μM
                v = params.k20 * g1^params.α2 / gT^params.α2 * S_uM * params.ET * 1000.0 / (S_uM + Km)  # ET从mM转为μM
                
                # 存储对数值，处理零值和负值
                v_matrix_three[i,j,k] = v > 0 ? log10(v) : -9  # 设置最小值为-9（对应10^-9 μM/s）
            end
        end
    end
    
    # 计算四步骤模型的活性
    for (i, ΔG1) in enumerate(ΔG1_range)
        for (j, ΔGT) in enumerate(ΔGT_range)
            # 更新参数
            new_params = ThermoKineticsParams(
                ΔGT = ΔGT, ΔG1 = ΔG1,
                R = params.R, T = params.T,
                α1 = params.α1, α2 = params.α2,
                k10 = params.k10, k20 = params.k20,
                ET = params.ET,
                model_type = "four_step"
            )
            
            # 计算结果
            res = calculate_kinetics(new_params)
            
            # 计算各底物浓度下的活性
            for (k, S) in enumerate(S_range)
                # 计算酶活性（单位：μM/s）
                S_uM = S * 1000.0  # mM到μM
                Ks = res.Ks
                k_forward = res.k_forward
                
                v = params.ET * 1000.0 * k_forward * S_uM / (1 + S_uM/Ks)  # ET从mM转为μM
                
                # 存储对数值，处理零值和负值
                v_matrix_four[i,j,k] = v > 0 ? log10(v) : -9  # 设置最小值为-9（对应10^-9 μM/s）
            end
        end
    end
    
    # 计算三步骤模型的Km = [S]线
    for (k, S) in enumerate(S_range)
        km_line_points = []
        for ΔG1 in ΔG1_range
            RT = params.R / 1000 * params.T  # kJ/mol
            S_uM = S * 1000.0  # 从mM转换为μM
            
            g1 = exp(ΔG1/RT)
            term1 = S_uM/g1 - 1
            
            if term1 > 0
                term2 = term1 / (params.k20/params.k10) * g1^(1-params.α1-params.α2)
                ΔGT_val = RT * log(term2)
                
                # 只在合理范围内添加点
                if ΔGT_val >= -40.0 && ΔGT_val <= 40.0
                    push!(km_line_points, (ΔG1, ΔGT_val))
                end
            end
        end
        push!(km_lines_three, km_line_points)
    end
    
    # 计算四步骤模型的Km = [S]线
    for (k, S) in enumerate(S_range)
        km_line_points = []
        for ΔG1 in ΔG1_range
            # 创建临时参数对象计算Ks
            for test_ΔGT in LinRange(-40.0, 40.0, 50)
                temp_params = ThermoKineticsParams(
                    ΔGT = test_ΔGT,
                    ΔG1 = ΔG1,
                    R = params.R,
                    T = params.T,
                    α1 = params.α1,
                    α2 = params.α2,
                    k10 = params.k10,
                    k20 = params.k20,
                    ET = params.ET,
                    model_type = "four_step"
                )
                
                res = calculate_kinetics(temp_params)
                Ks = res.Ks / 1000.0  # 转换为mM以与S比较
                
                # 如果Ks接近S，则添加点
                if abs(log10(Ks) - log10(S)) < 0.1
                    push!(km_line_points, (ΔG1, test_ΔGT))
                    break
                end
            end
        end
        push!(km_lines_four, km_line_points)
    end
    
    # 绘制比较图
    # 设置颜色范围
    clim = (-9, 2)  # 最小值（对应10^-9 μM/s）到最大值（对应10^2 = 100 μM/s）
    
    # 创建2x4网格图（左侧三步骤模型，右侧四步骤模型）
    plots = []
    labels = ["a", "b", "c", "d", "e", "f", "g", "h"]  # 子图标签
    
    # 绘制三步骤模型子图
    for (k, S) in enumerate(S_range)
        # 三步骤模型
        p_three = contourf(ΔG1_range, ΔGT_range, v_matrix_three[:,:,k]',
                 xlabel="ΔG₁ [kJ/mol]", ylabel="ΔG_T [kJ/mol]",
                 title="三步骤模型 [S]=$(S) mM",
                 color=:jet, levels=50,
                 linewidth=0.5, contour_labels=false,
                 clims=clim)
        
        # 添加子图标签
        annotate!(p_three, -18, 35, text(labels[k], :left, 14, :bold))
        
        # 绘制Km = [S]线
        if !isempty(km_lines_three[k])
            x_vals = [point[1] for point in km_lines_three[k]]
            y_vals = [point[2] for point in km_lines_three[k]]
            plot!(p_three, x_vals, y_vals, line=(:black, :dash, 2), label="")
        end
        
        push!(plots, p_three)
        
        # 四步骤模型
        p_four = contourf(ΔG1_range, ΔGT_range, v_matrix_four[:,:,k]',
                xlabel="ΔG₁ [kJ/mol]", ylabel="ΔG_T [kJ/mol]",
                title="四步骤模型 [S]=$(S) mM",
                color=:jet, levels=50,
                linewidth=0.5, contour_labels=false,
                clims=clim)
        
        # 添加子图标签
        annotate!(p_four, -18, 35, text(labels[k+4], :left, 14, :bold))
        
        # 绘制Km = [S]线
        if !isempty(km_lines_four[k])
            x_vals = [point[1] for point in km_lines_four[k]]
            y_vals = [point[2] for point in km_lines_four[k]]
            plot!(p_four, x_vals, y_vals, line=(:black, :dash, 2), label="")
        end
        
        push!(plots, p_four)
    end
    
    # 组合子图为2x4布局
    final_plot = plot(plots..., layout=(2,4), size=(1600, 800), margin=10mm, dpi=300)
    
    # 添加共享的颜色条，设置标题和范围
    plot!(final_plot, colorbar=true, colorbar_title="log v [μM/s]", clims=clim,
          right_margin=15mm)  # 增加右侧边距以容纳颜色条
    
    return final_plot
end

# 如果直接运行此文件，则生成比较图
if abspath(PROGRAM_FILE) == @__FILE__
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
        ET = 0.01      # mM，与Python代码一致
    )
    
    # 底物浓度列表 - 调整为对数刻度[0.1, 1, 10, 100]以匹配参考图像
    substrate_concentrations = [0.1, 1.0, 10.0, 100.0]
    
    # 绘制比较图
    p = compare_volcano_plots(substrate_concentrations, params)
    
    # 保存图像
    savefig(p, "enzyme_activity_volcano_comparison.png")
    println("三步骤和四步骤模型的火山图比较已保存为 enzyme_activity_volcano_comparison.png")
end