using Symbolics  # 符号计算支持
using Plots
using Measures  # 用于设置图表边距

struct ThermoKineticsParams
    ΔGT::Float64       # 总反应自由能变化 (kJ/mol)
    ΔG1::Float64       # 酶-底物复合物自由能变化 (kJ/mol)
    R::Float64         # 气体常数 (J/mol·K)，默认8.314
    T::Float64         # 温度 (K)，默认298
    α1::Float64        # BEP关系敏感系数1
    α2::Float64        # BEP关系敏感系数2
    k10::Float64       # 基准速率常数1
    k20::Float64       # 基准速率常数2
    ET::Float64        # 总酶浓度 (mM)
    model_type::String # 模型类型："three_step"或"four_step"
    
    # 关键字参数构造函数
    function ThermoKineticsParams(;
        ΔGT, ΔG1,
        R=8.314, T=300,  # 修改为300K，与Python代码一致
        α1=0.5, α2=0.5,  # 修改α2为0.5，与Python代码一致
        k10=1.0, k20=1.0, # 修改k20为1.0，与Python代码一致
        ET=0.01,         # 修改ET为0.01，与Python代码一致
        model_type="three_step") # 默认为三步骤模型
        new(ΔGT, ΔG1, R, T, α1, α2, k10, k20, ET, model_type)
    end
end

# 定义包含所有计算结果的结构体
# 使用isdefined检查避免重复定义
if !@isdefined(ThermoKineticsResults)
    struct ThermoKineticsResults
        ΔG2::Float64       # 产物形成自由能变化 (kJ/mol)
        ΔG3::Float64       # 仅四步骤模型：EP到E+P的自由能变化 (kJ/mol)
        Ea1::Float64       # 正反应活化能 (J/mol)
        Ea1r::Float64      # 逆反应活化能 (J/mol)
        Ea2::Float64       # 产物形成活化能 (J/mol)
        k1::Float64        # 正反应速率常数 (1/s)
        k1r::Float64       # 逆反应速率常数 (1/s)
        k2::Float64        # 产物形成速率常数 (1/ms)
        k3::Float64        # 仅四步骤模型：EP到E+P的速率常数 (1/s)
        k3r::Float64       # 仅四步骤模型：E+P到EP的速率常数 (1/(mM·s))
        Km::Float64        # 米氏常数 (μM)
        v_max::Float64     # 最大酶活性 (μM/s)
        # 简化模型参数
        Kp::Float64        # 产物解离常数 (μM)
        Ks::Float64        # 底物解离常数 (μM)
        k_forward::Float64 # 正向催化速率常数 (1/s)
        k_reverse::Float64 # 逆向催化速率常数 (1/s)
        model_type::String # 使用的模型类型
    end
end

function calculate_kinetics(params::ThermoKineticsParams)
    # 修改RT计算，与Python代码保持一致
    RT = params.R / 1000 * params.T  # 修改为kJ/mol，与Python代码一致
    
    # 热力学参数验证和转换
    # 移除硬编码验证，改为输出警告
    if params.ΔG1 >= 0
        @warn "ΔG1 = $(params.ΔG1) kJ/mol 为正值（吸能过程），计算结果可能不符合实际情况"
    end
    if params.ΔGT >= 0
        @warn "ΔGT = $(params.ΔGT) kJ/mol 为正值（吸能过程），计算结果可能不符合实际情况"
    end
    
    # 使用kJ/mol单位，与Python代码一致
    ΔG1 = params.ΔG1
    ΔGT = params.ΔGT
    
    # 计算g1和gT，与Python代码保持一致
    g1 = exp(ΔG1 / RT)
    gT = exp(ΔGT / RT)
    
    if params.model_type == "three_step"
        # 三步骤模型计算
        ΔG2 = ΔGT - ΔG1
        ΔG3 = 0.0  # 三步骤模型中不存在ΔG3
        
        # 验证自由能变化关系，但不输出警告以减少输出
        # 如果ΔG2为正值，反应可能不利于产物形成，但这是热力学允许的
        
        # 计算各步骤的平衡常数，确保热力学一致性
        K1 = exp(-ΔG1 / RT)  # 1/mM
        K2 = exp(-ΔG2 / RT)  # 无量纲
        
        # 计算三步骤模型的速率常数，使用与四步骤模型一致的公式
        # 第一步：E + S ⇌ ES
        k1 = params.k10 * exp(params.α1 * ΔG1 / RT)  # 正向 (1/(mM·s))
        k_minus_1 = params.k10 * exp((params.α1 - 1) * ΔG1 / RT)  # 逆向 (1/s)
        
        # 第二步：ES ⇌ E + P
        k2 = params.k20 * exp(params.α2 * ΔG2 / RT)  # 正向 (1/s)
        k_minus_2 = params.k20 * exp((params.α2 - 1) * ΔG2 / RT)  # 逆向 (1/(mM·s))
        
        # 在三步骤模型中，我们需要模拟四步骤模型的行为
        # 假设k3很大，k_minus_3很小，使EP中间体迅速转化为E+P
        k3 = 1000.0  # 一个较大的值
        k_minus_3 = k3 / K2  # 确保热力学一致性：k3/k_minus_3 = K2
        
        # 计算米氏常数Km（单位μM）
        # 使用热力学一致的公式
        Km_val = k_minus_1 / k1  # 单位mM
        Km = Km_val * 1e3  # 转换为μM
        
        # 计算最大反应速率v_max（单位μM/s）
        v_max = k2 * params.ET * 1e3  # ET单位为mM，转换为μM
        
        # 计算Kp（产物解离常数，单位μM）
        Kp_val = k_minus_2 / k2 * 1e3  # 转换为μM
        
        # 计算Ks（底物解离常数，单位μM）
        Ks_val = Km_val * 1e3  # 转换为μM
        
        # 计算k_forward（正向催化速率常数，单位1/s）
        k_forward_val = k2
        
        # 计算k_reverse（逆向催化速率常数，单位1/s）
        k_reverse_val = k_minus_1 * k_minus_2 / k1 / k2
        
    else  # "four_step" 模型
        # 四步骤模型计算
        # 在四步骤模型中，我们需要将总自由能变化ΔGT分配给三个步骤
        # 假设ΔG1已经给定，我们需要分配ΔG2和ΔG3
        # 为简单起见，我们假设ΔG2 = ΔG3 = (ΔGT - ΔG1) / 2
        
        ΔG2 = (ΔGT - ΔG1) / 2
        ΔG3 = (ΔGT - ΔG1) / 2
        
        # 计算各步骤的平衡常数
        K1 = exp(-ΔG1 / RT)  # 1/mM
        K2 = exp(-ΔG2 / RT)  # 无量纲
        K3 = exp(-ΔG3 / RT)  # mM
        
        # 计算四步骤模型的速率常数
        # 第一步：E + S ⇌ ES
        k1 = params.k10 * exp(params.α1 * ΔG1 / RT)  # 正向 (1/(mM·s))
        k_minus_1 = params.k10 * exp((params.α1 - 1) * ΔG1 / RT)  # 逆向 (1/s)
        
        # 第二步：ES ⇌ EP
        k2 = params.k20 * exp(params.α2 * ΔG2 / RT)  # 正向 (1/s)
        k_minus_2 = params.k20 * exp((params.α2 - 1) * ΔG2 / RT)  # 逆向 (1/s)
        
        # 第三步：EP ⇌ E + P
        # 假设第三步的α值也为0.5
        α3 = 0.5
        k3 = params.k20 * exp(α3 * ΔG3 / RT)  # 正向 (1/s)
        k_minus_3 = params.k20 * exp((α3 - 1) * ΔG3 / RT)  # 逆向 (1/(mM·s))
        
        # 计算四步骤模型的简化参数
        # 根据model_cal.jl中的推导
        Kp_val = (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2) / (k_minus_3 * (k2 + k_minus_1 + k_minus_2)) * 1e3  # 转换为μM
        Ks_val = (k_minus_1 * (k2 + k3 + k_minus_2)) / (k1 * (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2)) * 1e3  # 转换为μM
        k_forward_val = k2*k3 / (k2 + k3 + k_minus_2)
        k_reverse_val = k_minus_1*k_minus_2 / (k2 + k_minus_1 + k_minus_2)
        
        # 计算米氏常数Km（单位μM）
        # 在四步骤模型中，Km = Ks*(1 + [P]/Kp)
        # 假设[P] = 0，则Km = Ks
        Km = Ks_val
        
        # 计算最大反应速率v_max（单位μM/s）
        # v_max = k_forward * ET
        v_max = k_forward_val * params.ET * 1e3  # ET单位为mM，转换为μM
    end
    
    # 验证热力学一致性：k_forward/k_reverse = exp(-ΔGT/RT)
    # 理论上，k_forward/k_reverse应该等于exp(-ΔGT/RT)
    thermo_ratio = k_forward_val/k_reverse_val
    expected_ratio = exp(-ΔGT/RT)
    
    # 只有在差异非常大时才输出警告（超过5%），减少警告输出
    if abs(thermo_ratio/expected_ratio - 1.0) > 0.05
        @warn "热力学一致性检查：k_forward/k_reverse = $(thermo_ratio)，期望值 = $(expected_ratio)，差异 = $(100*(thermo_ratio/expected_ratio-1))%"
    end
    
    # 计算活化能 (J/mol)
    Ea1 = -params.R * params.T * log(k1 / params.k10)  # 第一步正向反应活化能
    Ea1r = -params.R * params.T * log(k_minus_1 / params.k10)  # 第一步逆向反应活化能
    Ea2 = -params.R * params.T * log(k2 / params.k20)  # 第二步反应活化能
    
    return ThermoKineticsResults(ΔG2, ΔG3, Ea1, Ea1r, Ea2, k1, k_minus_1, k2, k3, k_minus_3, Km, v_max, 
                                Kp_val, Ks_val, k_forward_val, k_reverse_val, params.model_type)
end

# 火山图可视化函数
function plot_activity_volcano(params::ThermoKineticsParams, S_range::Vector{Float64})
    
    # 生成ΔG1和ΔGT网格
    ΔG1_range = LinRange(-20.0, 20.0, 100)  # 调整范围与参考图像一致
    ΔGT_range = LinRange(-40.0, 40.0, 100)  # 调整范围与参考图像一致
    
    # 准备存储结果
    v_matrix = zeros(length(ΔG1_range), length(ΔGT_range), length(S_range))
    # 准备存储Km = [S]线的数据
    km_lines = []
    
    # 计算各点活性
    for (i, ΔG1) in enumerate(ΔG1_range)
        for (j, ΔGT) in enumerate(ΔGT_range)
            # 确保ΔG2在合理范围内
            ΔG2 = ΔGT - ΔG1
            
            # 更新参数
            new_params = ThermoKineticsParams(
                ΔGT = ΔGT, ΔG1 = ΔG1,
                R = params.R, T = params.T,
                α1 = params.α1, α2 = params.α2,
                k10 = params.k10, k20 = params.k20,
                ET = params.ET,
                model_type = params.model_type
            )
            
            # 计算结果
            res = calculate_kinetics(new_params)
            
            # 计算各底物浓度下的活性
            for (k, S) in enumerate(S_range)
                # 计算酶活性（单位：μM/s）
                S_uM = S * 1000.0  # mM到μM
                
                if params.model_type == "three_step"
                    # 使用热力学一致的表达式（三步骤模型）
                    RT = params.R / 1000 * params.T  # kJ/mol
                    
                    # 计算各步骤的平衡常数
                    K1 = exp(-ΔG1 / RT)  # 1/mM
                    K2 = exp(-(ΔGT - ΔG1) / RT)  # 无量纲
                    
                    # 计算速率常数
                    k1 = params.k10 * exp(params.α1 * ΔG1 / RT)  # 正向 (1/(mM·s))
                    k_minus_1 = params.k10 * exp((params.α1 - 1) * ΔG1 / RT)  # 逆向 (1/s)
                    k2 = params.k20 * exp(params.α2 * (ΔGT - ΔG1) / RT)  # 正向 (1/s)
                    
                    # 计算米氏常数
                    Km = k_minus_1 / k1 * 1e3  # 转换为μM
                    
                    # 计算酶活性
                    v = k2 * S_uM * params.ET * 1000.0 / (S_uM + Km)  # ET从mM转为μM
                else
                    # 四步骤模型计算
                    # 使用简化表达式：v = (ET * (k_forward * S - k_reverse * P)) / (1 + S/Ks + P/Kp)
                    # 假设P = 0
                    Ks = res.Ks
                    k_forward = res.k_forward
                    
                    v = params.ET * 1000.0 * k_forward * S_uM / (1 + S_uM/Ks)  # ET从mM转为μM
                end
                
                # 存储对数值，处理零值和负值
                v_matrix[i,j,k] = v > 0 ? log10(v) : -9  # 设置最小值为-9（对应10^-9 μM/s）
            end
        end
    end
    
    # 计算Km = [S]线
    # 对每个底物浓度，计算满足Km = [S]的ΔG1和ΔGT关系
    for (k, S) in enumerate(S_range)
        km_line_points = []
        for ΔG1 in ΔG1_range
            # 根据热力学关系计算当Km = S时对应的ΔGT
            RT = params.R / 1000 * params.T  # kJ/mol
            
            # 注意：S单位为mM，需要转换为μM以匹配Km单位
            S_uM = S * 1000.0  # 从mM转换为μM
            
            if params.model_type == "three_step"
                # 三步骤模型的Km = [S]线
                # 求解方程：S*1000 = g1*(1+K)，其中K = k20/k10/gT^α2 * g1^(α1+α2-1)
                # 变形为：ΔGT = RT * ln((S*1000/g1 - 1)/(k20/k10) * g1^(1-α1-α2))
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
            else  # 四步骤模型
                # 四步骤模型的Km = [S]线
                # 在四步骤模型中，Km = Ks*(1 + [P]/Kp)
                # 假设[P] = 0，则Km = Ks
                # 当Ks = S时，我们需要求解ΔGT
                
                # 创建临时参数对象计算Ks
                temp_params = ThermoKineticsParams(
                    ΔGT = 0.0,  # 临时值，会在循环中更新
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
                
                # 尝试不同的ΔGT值，找到使Ks ≈ S的值
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
        end
        push!(km_lines, km_line_points)
    end
    
    # 绘制火山图
    plots = []
    labels = ["a", "b", "c", "d"]  # 子图标签
    
    # 设置颜色范围为[-9, 2]，与参考图像一致
    clim = (-9, 2)
    
    for (k, S) in enumerate(S_range)
        # 使用jet配色方案，与论文中的图像一致
        p = contourf(ΔG1_range, ΔGT_range, v_matrix[:,:,k]',
                     xlabel="ΔG₁ [kJ/mol]", ylabel="ΔG_T [kJ/mol]",
                     title="[S] = $(S) mM",  # 添加单位并使格式更清晰
                     color=:jet, levels=50,  # 增加等高线数量以获得更平滑的过渡
                     linewidth=0.5, contour_labels=false,
                     clims=clim)  # 设置颜色范围
        
        # 添加子图标签
        annotate!(p, -18, 35, text(labels[k], :left, 14, :bold))
        
        # 绘制Km = [S]线
        if !isempty(km_lines[k])
            x_vals = [point[1] for point in km_lines[k]]
            y_vals = [point[2] for point in km_lines[k]]
            # 只在第一个子图显示标签，其他子图不显示标签
            if k == 1
                plot!(p, x_vals, y_vals, line=(:black, :dash, 2), label="Km = [S]")
            else
                plot!(p, x_vals, y_vals, line=(:black, :dash, 2), label="")
            end
        end
        
        push!(plots, p)
    end
    
    # 组合子图为2x2布局
    final_plot = plot(plots..., layout=(2,2), size=(800, 800), margin=10mm)
    
    # 添加共享的颜色条，设置标题和范围
    plot!(final_plot, colorbar=true, colorbar_title="log v [μM/s]", clims=clim,
          right_margin=15mm)  # 增加右侧边距以容纳颜色条
    
    return final_plot
end

# 示例用法
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
    
    results = calculate_kinetics(params)
    
    # 输出结果验证
    println("\n计算结果验证:")
    println("ΔG2 = ", results.ΔG2, " kJ/mol")
    println("Km = ", results.Km, " μM") 
    println("v_max = ", results.v_max, " μM/s")
    
    # 输出简化模型参数
    println("\n简化模型参数:")
    println("Kp = ", results.Kp, " μM")
    println("Ks = ", results.Ks, " μM")
    println("k_forward = ", results.k_forward, " 1/s")
    println("k_reverse = ", results.k_reverse, " 1/s")
    
    # 计算并输出简化表达式
    println("\n简化表达式: v = (Et*(k_forward*S-k_reverse*P))/(1+S/Ks+P/Kp)")
    
    # 计算不同底物浓度下的反应速率
    S_values = [0.1, 1.0, 10.0, 100.0]  # mM
    P = 0.0  # 初始无产物
    println("\n不同底物浓度下的反应速率:")
    for S in S_values
        # 使用与Python代码一致的表达式计算v
        RT = params.R / 1000 * params.T  # kJ/mol
        g1 = exp(params.ΔG1/RT)
        gT = exp(params.ΔGT/RT)
        K = params.k20/params.k10/gT^params.α2 * g1^(params.α1+params.α2-1)
        Km_val = g1*(1+K)
        Km = Km_val * 1e3  # 转换为μM
        
        S_uM = S * 1000.0  # mM到μM
        v = params.k20 * g1^params.α2 / gT^params.α2 * S_uM * params.ET * 1000.0 / (S_uM + Km)  # ET从mM转为μM
        
        println("[S] = ", S, " mM, v = ", v, " μM/s, log(v) = ", log10(v))
    end
end