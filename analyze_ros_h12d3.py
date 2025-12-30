"""分析 H12D3 的 ROS 動態"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model_ros import ALL_PARAMS, UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, SIMULATION, get_env_for_treatment

# 分析當前參數下 H12D3 的 ROS 動態
p = UVAParams()
env = get_env_for_treatment('H12D3')

fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
C_buf_init = Xd_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6
ROS_init_mM = 0.01
initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, ROS_init_mM]

transplant_day = SIMULATION['transplant_offset']
simulation_days = SIMULATION['days']
t_start = transplant_day * 86400
t_end = (transplant_day + simulation_days) * 86400

# 詳細模擬
sol = solve_ivp(
    uva_sun_derivatives_ros,
    (t_start, t_end),
    initial_state,
    args=(p, env),
    method='RK45',
    max_step=60,
    t_eval=np.linspace(t_start, t_end, 1000)
)

if sol.success:
    t_days = sol.t / 86400
    ROS = sol.y[5, :]
    Stress = sol.y[4, :]

    print("=" * 80)
    print("H12D3 ROS 動態分析")
    print("=" * 80)
    print(f"\n當前 ROS 參數:")
    print(f"  ros_prod_coeff = {p.ros_prod_coeff:.0e}")
    print(f"  ros_vmax_apx = {p.ros_vmax_apx}")
    print(f"  ros_vmax_cat = {p.ros_vmax_cat}")
    print(f"  K_ros_apx = {p.K_ros_apx} mM")
    print(f"  K_ros_cat = {p.K_ros_cat} mM")
    print(f"  ros_damage_gain = {p.ros_damage_gain}")
    print(f"  K_ros_damage = {p.K_ros_damage} mM")

    print(f"\nROS 動態:")
    print(f"  初始 ROS: {ROS[0]:.4f} mM")
    print(f"  最大 ROS: {np.max(ROS):.4f} mM")
    print(f"  最終 ROS: {ROS[-1]:.4f} mM")
    print(f"  最終 Stress: {Stress[-1]:.2f}")

    # 計算 ROS 生成和清除速率 (以 12h 照射期間)
    I_UVA = 22.0  # W/m²
    ros_prod = p.ros_prod_coeff * I_UVA
    ros_level = 1.0  # 假設 1 mM
    ros_apx = p.ros_vmax_apx * ros_level / (p.K_ros_apx + ros_level)
    ros_cat = p.ros_vmax_cat * ros_level / (p.K_ros_cat + ros_level)
    ros_decay = p.ros_passive_decay * ros_level

    print(f"\n速率分析 (假設 ROS = 1 mM):")
    print(f"  生成速率 (22 W/m²): {ros_prod:.6f} mM/s = {ros_prod*3600:.4f} mM/h")
    print(f"  APX 清除 (@1mM): {ros_apx:.6f} mM/s = {ros_apx*3600:.4f} mM/h")
    print(f"  CAT 清除 (@1mM): {ros_cat:.6f} mM/s = {ros_cat*3600:.4f} mM/h")
    print(f"  被動衰減 (@1mM): {ros_decay:.6f} mM/s = {ros_decay*3600:.4f} mM/h")
    print(f"  總清除: {(ros_apx + ros_cat + ros_decay)*3600:.4f} mM/h")
    print(f"  淨變化: {(ros_prod - ros_apx - ros_cat - ros_decay)*3600:.4f} mM/h")

    # 計算穩態 ROS
    # 在穩態: ros_prod = ros_apx + ros_cat + ros_decay
    # 這是一個複雜方程，簡化分析
    print(f"\n問題診斷:")
    if ros_prod < ros_apx + ros_cat + ros_decay:
        print(f"  ❌ 生成速率太低，ROS 無法累積")
        print(f"  建議: 提高 ros_prod_coeff 或降低清除速率")
    else:
        print(f"  ✓ 生成速率高於清除，ROS 可累積")

    # 計算需要的 ros_prod_coeff
    target_steady_ros = 5.0  # 目標穩態 ROS
    # 簡化假設：主要靠 APX 清除 (因為 K_ros_apx << K_ros_cat)
    needed_clearance = p.ros_vmax_apx * target_steady_ros / (p.K_ros_apx + target_steady_ros)
    needed_prod = needed_clearance / I_UVA
    print(f"\n建議 ros_prod_coeff (目標穩態 ROS = {target_steady_ros} mM): {needed_prod:.2e}")

print("\n" + "=" * 80)
