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
    
    # 关键字参数构造函数
    function ThermoKineticsParams(;
        ΔGT, ΔG1,
        R=8.314, T=298,
        α1, α2, k10, k20, ET)
        new(ΔGT, ΔG1, R, T, α1, α2, k10, k20, ET)
    end
end

# 定义包含所有计算结果的结构体
# 使用isdefined检查避免重复定义
if !@isdefined(ThermoKineticsResults)
    struct ThermoKineticsResults
        ΔG2::Float64       # 产物形成自由能变化 (kJ/mol)
        Ea1::Float64       # 正反应活化能 (J/mol)
        Ea1r::Float64      # 逆反应活化能 (J/mol)
        Ea2::Float64       # 产物形成活化能 (J/mol)
        k1::Float64        # 正反应速率常数 (1/s)
        k1r::Float64       # 逆反应速率常数 (1/s)
        k2::Float64        # 产物形成速率常数 (1/ms)
        Km::Float64        # 米氏常数 (μM)
        v_max::Float64     # 最大酶活性 (μM/s)
        # 简化模型参数
        Kp::Float64        # 产物解离常数 (μM)
        Ks::Float64        # 底物解离常数 (μM)
        k_forward::Float64 # 正向催化速率常数 (1/s)
        k_reverse::Float64 # 逆向催化速率常数 (1/s)
    end
end

function calculate_kinetics(params::ThermoKineticsParams)
    # 单位转换（kJ → J）并验证输入
    ΔGT = params.ΔGT * 1000
    ΔG1 = params.ΔG1 * 1000
    
    # 热力学参数验证和转换
    # 移除硬编码验证，改为输出警告
    if params.ΔG1 >= 0
        @warn "ΔG1 = $(params.ΔG1) kJ/mol 为正值（吸能过程），计算结果可能不符合实际情况"
    end
    if params.ΔGT >= 0
        @warn "ΔGT = $(params.ΔGT) kJ/mol 为正值（吸能过程），计算结果可能不符合实际情况"
    end
    
    # 单位转换（kJ → J）
    ΔG1 = params.ΔG1 * 1000
    ΔGT = params.ΔGT * 1000
    ΔG2 = ΔGT - ΔG1
    
    # 验证自由能变化关系
    # 使用更宽松的范围，避免过多警告
    if !(ΔG2 < 0)
        @warn "ΔG2 = $(ΔG2/1000) kJ/mol 为正值，反应可能不利于产物形成"
    end
    
    # BEP关系计算活化能（修正活化能基准值计算）
    # 假设指前因子A1=1e13 s⁻¹, A2=1e12 s⁻¹ (典型酶反应值)
    A1 = 1e13
    A2 = 1e12
    Ea10 = -log(params.k10/A1) * params.R * params.T  # 正确阿伦尼乌斯方程推导
    Ea20 = -log(params.k20/A2) * params.R * params.T
    Ea1 = Ea10 + params.α1 * ΔG1
    Ea1r = Ea10 + params.α1 * (ΔG1 + params.R * params.T)
    Ea2 = Ea20 + params.α2 * (ΔG2 + params.R * params.T)
    
    # 阿伦尼乌斯方程计算速率常数
    g1 = exp(ΔG1 / (params.R * params.T))
    gT = exp(ΔGT / (params.R * params.T))
    k1 = params.k10 * g1^(-params.α1)
    k1r = params.k10 * g1^(1 - params.α1)  # α1r = 1 - α1
    k2 = params.k20 * gT^params.α2
    
    # 保持mM单位体系计算
    # 计算平衡常数（无量纲）
    K = (params.k20 * gT) / (params.k10 * g1^params.α1)
    
    # 基准转换常数
    C0 = 1.0  # mM基准
    
    # 单位转换基准
    mmol_to_mol = 1000.0  # 1M = 1000mM
    
    # 最终修正的计算模型
    
    # 1. 速率常数单位转换（基准k10/k20单位为1/s）
    # 保持原始单位，不进行转换
    k10 = params.k10  # 1/s
    k20 = params.k20  # 1/s
    
    # 2. 计算米氏常数Km（单位μM）
    # 使用热力学关系计算
    RT = params.R * params.T  # J/mol
    Km_val = exp(ΔG1/(RT)) * (1 + (k20/k10)*exp(ΔGT/(RT)))
    Km = Km_val * 1e3  # 转换为μM
    
    # 3. 计算k2（单位1/s）
    k2_val = k20 * exp(-params.α2*ΔG2/(RT))
    
    # 4. 计算最大反应速率（单位μM/s）
    vmax_val = k2_val * params.ET * 1e3  # ET单位为mM，转换为μM
    v_max = vmax_val  # μM/s
    
    # 计算简化模型参数（基于model_cal.jl的推导）
    # 将三步骤模型（E+S ⇌ ES ⇌ E+P）与四步骤模型（E+S ⇌ ES ⇌ EP ⇌ E+P）对应
    # 根据热力学关系，我们可以将速率常数与自由能变化关联
    
    # 从自由能变化计算速率常数
    # 使用阿伦尼乌斯方程和BEP关系
    RT = params.R * params.T  # J/mol
    
    # 计算三步骤模型的速率常数
    # 第一步：E + S ⇌ ES
    k1 = params.k10 * exp(-params.α1 * ΔG1/RT)  # 正向 (1/(mM·s))
    k_minus_1 = params.k10 * exp((1-params.α1) * ΔG1/RT)  # 逆向 (1/s)
    
    # 第二步：ES ⇌ E + P
    k2 = params.k20 * exp(-params.α2 * ΔG2/RT)  # 正向 (1/s)
    k_minus_2 = params.k20 * exp((1-params.α2) * ΔG2/RT)  # 逆向 (1/(mM·s))
    
    # 计算简化参数，使其与四步骤模型对应
    # 在四步骤模型中：
    # Kp = (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2) / (k_minus_3 * (k2 + k_minus_1 + k_minus_2))
    # Ks = (k_minus_1 * (k2 + k3 + k_minus_2)) / (k1 * (k2*k3 + k3*k_minus_1 + k_minus_1*k_minus_2))
    # k_forward = k2*k3 / (k2 + k3 + k_minus_2)
    # k_reverse = k_minus_1*k_minus_2 / (k2 + k_minus_1 + k_minus_2)
    
    # 在三步骤模型中，我们需要模拟四步骤模型的行为
    # 假设k3很大，k_minus_3很小，使EP中间体迅速转化为E+P
    # 设定一个大的k3值和小的k_minus_3值来模拟四步骤模型
    k3_eff = 1000.0  # 一个较大的值
    k_minus_3_eff = 0.001  # 一个较小的值
    
    # 计算Kp（产物解离常数，单位μM），使用与四步骤模型相同的公式
    Kp_val = (k2*k3_eff + k3_eff*k_minus_1 + k_minus_1*k_minus_2) / (k_minus_3_eff * (k2 + k_minus_1 + k_minus_2)) * 1e3  # 转换为μM
    
    # 计算Ks（底物解离常数，单位μM），使用与四步骤模型相同的公式
    Ks_val = (k_minus_1 * (k2 + k3_eff + k_minus_2)) / (k1 * (k2*k3_eff + k3_eff*k_minus_1 + k_minus_1*k_minus_2)) * 1e3  # 转换为μM
    
    # 计算k_forward（正向催化速率常数，单位1/s），使用与四步骤模型相同的公式
    k_forward_val = k2*k3_eff / (k2 + k3_eff + k_minus_2)
    
    # 计算k_reverse（逆向催化速率常数，单位1/s），使用与四步骤模型相同的公式
    k_reverse_val = k_minus_1*k_minus_2 / (k2 + k_minus_1 + k_minus_2)
    
    # 验证热力学一致性：k_forward/k_reverse = exp(-ΔGT/RT)
    # 理论上，k_forward/k_reverse应该等于exp(-ΔGT/RT)
    thermo_ratio = k_forward_val/k_reverse_val
    expected_ratio = exp(-ΔGT/RT)
    
    # 如果比值相差超过1%，输出警告
    if abs(thermo_ratio/expected_ratio - 1.0) > 0.01
        @warn "热力学一致性检查：k_forward/k_reverse = $(thermo_ratio)，期望值 = $(expected_ratio)，差异 = $(100*(thermo_ratio/expected_ratio-1))%"
    end
    
    return ThermoKineticsResults(ΔG2/1000, Ea1, Ea1r, Ea2, k1, k1r, k2, Km, v_max, 
                                Kp_val, Ks_val, k_forward_val, k_reverse_val)
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
                ET = params.ET
            )
            
            # 计算结果
            res = calculate_kinetics(new_params)
            
            # 计算各底物浓度下的活性
            for (k, S) in enumerate(S_range)
                # 计算酶活性（单位：μM/s）
                # 使用简化表达式计算：v = (Et*(k_forward*S-k_reverse*P))/(1+S/Ks+P/Kp)
                # 在这个火山图中，我们假设P=0（初始反应条件）
                P = 0.0  # 产物浓度为零
                
                # 单位转换：S和P从mM转换为μM以匹配Ks和Kp
                S_uM = S * 1000.0  # mM到μM
                P_uM = P * 1000.0  # mM到μM
                
                # 使用简化表达式计算酶活性
                # 注意：S和P单位为mM，而Ks和Kp单位为μM
                # 因此需要将S和P转换为μM，或将Ks和Kp转换为mM
                # 这里我们选择将S和P转换为μM
                numerator = new_params.ET * 1000.0 * (res.k_forward * S_uM - res.k_reverse * P_uM)  # ET从mM转为μM
                denominator = 1.0 + S_uM/res.Ks + P_uM/res.Kp
                v = numerator / denominator  # 结果单位为μM/s
                
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
            # 使用简化的关系：ln(Km) = ΔG1/(RT) + ln(1 + k2/k1*exp(ΔGT/(RT)))
            # 当Km = S时，求解ΔGT
            RT = params.R * params.T  # J/mol
            k_ratio = params.k20 / params.k10
            
            # 求解方程：S*1000 = exp(ΔG1/RT) * (1 + k_ratio*exp(ΔGT/RT))
            # 变形为：ΔGT = RT * ln((S*1000/exp(ΔG1/RT) - 1)/k_ratio)
            # 注意：S单位为mM，需要转换为μM以匹配Km单位
            S_uM = S * 1000.0  # 从mM转换为μM
            term1 = S_uM / exp(ΔG1*1000/RT)  # 确保使用正确的RT单位（J/mol而非kJ/mol）
            if term1 > 1  # 确保对数参数为正
                ΔGT_val = RT * log((term1 - 1)/k_ratio)
                if ΔGT_val >= -40.0 && ΔGT_val <= 40.0  # 确保在图表范围内
                    push!(km_line_points, (ΔG1, ΔGT_val))
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
        # 使用viridis配色方案，与参考图像类似
        p = contourf(ΔG1_range, ΔGT_range, v_matrix[:,:,k]',
                     xlabel="ΔG₁ [kJ/mol]", ylabel="ΔG_T [kJ/mol]",
                     title="[S] = $(S) mM",  # 添加单位并使格式更清晰
                     color=:viridis, levels=30,  # 增加等高线数量以获得更平滑的过渡
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
        v = (params.ET * (results.k_forward * S - results.k_reverse * P)) / 
            (1.0 + S/results.Ks + P/results.Kp)
        println("[S] = ", S, " mM, v = ", v, " μM/s, log(v) = ", log10(v))
    end
end
